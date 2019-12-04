# Functions needed to process process and run gibbs sampler

# Get expected value of prior probability on full mixing proportions
# Based on empirical concordance with various paths in red_class
get_prior_prop <- function(red_class, fits, d, n, dist_tol = 0, MAP = TRUE) {
    mus <- map(fits, "mu")
    
    if (MAP) {
        labels <- purrr::map(fits, "cluster") %>% simplify2array    
    } else {
        labels <- purrr::map(fits, "post_prob") %>%
            lapply(function(X) apply(X, 1, function(Y) base::sample(1:ncol(X), size = 1, prob = Y))) %>%
            simplify2array
    }
    
    prior_count <-
        cget_prior_count(red_class, mus, labels, d, n, dist_tol)
    (prior_count / n) / (sum(prior_count / n))
}

# This function will filter out classes that are really similar to each other
# since they probably are not that informative and this task of computing
# prior proportions (alpha) is super burdensome when there are too many classes
reduce_by_hamming <- function(red_class, hamming_tol = 1, force_canonical = TRUE) {
    M <- nrow(red_class)
    fullidx <- seq(M)
    out <- rep(NA, M)
    
    count <- 1
    while(length(fullidx) != 0) {
        candidate_indices <- creduce_by_hamming(as.matrix(red_class), fullidx-1, hamming_tol, M)
        candidate_indices <- fullidx[which(candidate_indices == 1)]
        
        if(length(candidate_indices) == 1) {
            keepidx <- candidate_indices
        } else {
            keepidx <- sample(candidate_indices, size = 1)    
        }
        
        out[count] <- keepidx
        fullidx <- fullidx[!(fullidx %in% candidate_indices)]
        count <- count + 1
    }
    
    outidx <- out[!is.na(out)]
    
    # If we have canonical behavior, keep it
    if(force_canonical) {
        allneg <- which(sapply(1:M, function(X) all(red_class[X,] == -1)))
        allone <- which(sapply(1:M, function(X) all(red_class[X,] == 0)))
        allpos <- which(sapply(1:M, function(X) all(red_class[X,] == 1)))
        outidx <- unique(c(outidx, allneg, allpos, allone))
    }
    
    red_class[outidx,]
}

# Using cluster labels to determine the means, covariances
get_hyperparams <- function(fits, d, red_class, dat, clusters, var_quantile = 0.75) {
    mu0_temp_pos <- mu0_temp_neg <- rep(NA, d)
    Psi0_temp_pos <- Psi0_temp_neg <- matrix(rep(NA, d ^ 2), nrow = d, ncol = d)
    spl <- names(fits) %>%
        str_split(pattern = "_") %>%
        purrr::map(as.numeric) %>%
        abind::abind(along = 2) %>%
        t
    mus <- purrr::map(fits, "mu") %>% purrr::map(abs)
    sigmas <- purrr::map(fits, "sigma")
    rhos <- purrr::map(fits, "rho") %>% purrr::map(abs)
    combos <- purrr::map(fits, "combos")
    
    # For each dimension
    for (i in seq(d)) {
        idx <- t(apply(spl, 1, function(X)
            X == i))
        muvec_pos <- sigmavec_pos <-
            muvec_neg <- sigmavec_neg <- rep(0, sum(idx))
        count <- 1
        for(j in 1:nrow(idx)) {
            if(any(idx[j,])) {
                pidx <- which(idx[j,])
                subcl <- combos[[j]][,pidx]
                
                if(any(subcl == 1)) {
                    count2 <- count3 <- 0
                    posidx <- which(subcl == 1)
                    for(k in posidx) {
                        if(sum(clusters[[j]] == k) > 0) {
                            muvec_pos[count] <- muvec_pos[count] + mean(dat[clusters[[j]] == k, i])
                            count2 <- count2 + 1
                        }
                        if(sum(clusters[[j]] == k) > 1) {
                            if(var(dat[clusters[[j]] == k, i]) > sigmavec_pos[count]) {
                                sigmavec_pos[count] <- var(dat[clusters[[j]] == k, i])  
                            }
                            # sigmavec_pos[count] <- sigmavec_pos[count] + var(dat[clusters[[j]] == k, i])
                            count3 <- count3 + 1
                        }
                    }
                    muvec_pos[count] <- muvec_pos[count] / count2
                    #sigmavec_pos[count] <- sigmavec_pos[count] / count3
                }
                
                if(any(subcl == -1)) {
                    count2 <- count3 <- 0
                    negidx <- which(subcl == -1)
                    for(k in negidx) {
                        if(sum(clusters[[j]] == k) > 0) {
                            muvec_neg[count] <- muvec_neg[count] + mean(dat[clusters[[j]] == k, i])
                            count2 <- count2 + 1
                        }
                        if(sum(clusters[[j]] == k) > 1) {
                            #sigmavec_neg[count] <- sigmavec_neg[count] + var(dat[clusters[[j]] == k, i])
                            if(var(dat[clusters[[j]] == k, i]) > sigmavec_neg[count]) {
                                sigmavec_neg[count] <- var(dat[clusters[[j]] == k, i])  
                            }
                            count3 <- count3 + 1
                        }
                    }
                    muvec_neg[count] <- muvec_neg[count] / count2
                    #sigmavec_neg[count] <- sigmavec_neg[count] / count3
                }
                count <- count + 1
            }
        }
        mu0_temp_pos[i] <- mean(muvec_pos[muvec_pos != 0])
        mu0_temp_neg[i] <- mean(muvec_neg[muvec_neg != 0]) # will be NaN if subset has length 0
        Psi0_temp_pos[i,i] <- quantile(sigmavec_pos[sigmavec_pos != 0], var_quantile)
        Psi0_temp_neg[i,i] <- quantile(sigmavec_neg[sigmavec_neg != 0], var_quantile)
    }
    
    if(any(is.na(diag(Psi0_temp_pos)))) {
        mn <- mean(diag(Psi0_temp_pos), na.rm = TRUE)
        diag(Psi0_temp_pos)[is.na(diag(Psi0_temp_pos))] <- mn
    }
    
    if(any(is.na(diag(Psi0_temp_neg)))) {
        mn <- mean(diag(Psi0_temp_neg), na.rm = TRUE)
        diag(Psi0_temp_neg)[is.na(diag(Psi0_temp_neg))] <- mn
    }
    
    for (i in seq(nrow(spl))) {
        subcl <- combos[[i]]
        posassocidx <- apply(subcl, 1, function(X) all(X == 1) | all(X == -1))
        negassocidx <- apply(subcl, 1, function(X) (X[1] == 1 & X[2] == -1) | (X[1] == -1 & X[2] == 1))
        if(sum(posassocidx) > 0) {
            for(j in which(posassocidx)) {
                Psi0_temp_pos[spl[i, 1], spl[i, 2]] <-
                    Psi0_temp_pos[spl[i, 2], spl[i, 1]] <- cov(dat[clusters[[i]] == j, spl[i,]])[1,2]
            }
        }
        if(sum(negassocidx) > 0) {
            for(j in which(negassocidx)) {
                Psi0_temp_neg[spl[i, 1], spl[i, 2]] <-
                    Psi0_temp_neg[spl[i, 2], spl[i, 1]] <- cov(dat[clusters[[i]] == j, spl[i,]])[1,2]
            }
        }
    }
    
    mu0 <- matrix(rep(0, length(red_class)),
                  nrow = nrow(red_class),
                  ncol = ncol(red_class))
    Psi0 <- base::array(rep(0, length(Psi0_temp_pos) * nrow(red_class)),
                        dim = c(d, d, nrow(red_class)))
    
    # Update mu0 based on association patterns
    for (i in seq(ncol(red_class))) {
        mu0[red_class[, i] == 1, i] <- mu0_temp_pos[i]
        mu0[red_class[, i] == -1, i] <- mu0_temp_neg[i]
    }
    
    # Update Psi0 based on association patterns
    for (i in seq(nrow(red_class))) {
        diag(Psi0[, , i])[red_class[i,] == 1] <- diag(Psi0_temp_pos)[red_class[i,] == 1]
        diag(Psi0[, , i])[red_class[i,] == -1] <- diag(Psi0_temp_neg)[red_class[i,] == -1]
        diag(Psi0[, , i])[red_class[i,] == 0] <- rep(1, sum(red_class[i,] == 0))
        
        # Don't update rho if there arent multiple associations
        assoc_idx <- red_class[i, ] != 0
        if (sum(assoc_idx) > 1) {
            wai <- which(assoc_idx)
            
            # for every pairwise associated combination
            cmb <- combn(wai, 2)
            for (j in seq(ncol(cmb))) {
                if(prod(red_class[i, cmb[,j]]) == 1) {
                    Psi0[cmb[1, j], cmb[2, j], i] <-
                        Psi0[cmb[2, j], cmb[1, j], i] <- Psi0_temp_pos[cmb[1,j], cmb[2,j]]
                } else {
                    Psi0[cmb[1, j], cmb[2, j], i] <-
                        Psi0[cmb[2, j], cmb[1, j], i] <- Psi0_temp_neg[cmb[1,j], cmb[2,j]]
                }
            }
        }
    }
    list("mu0" = mu0, "Psi0" = Psi0)
}

# Count how many observations belong to each class
get_kappa <- function(z, nclass) {
    kappa <- rep(0, nclass)
    runlens <- rle(sort(z))
    kappa[runlens$values] <- runlens$lengths
    kappa
}