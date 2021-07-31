merge_classes_multichain <- function(n_groups, chain_list, burnin_list) {
    if(typeof(burnin_list) != "list") {
        burnin_list <- map(1:length(chain_list), ~burnin_list)
    }
    newclustlabs <- cumsum(map_int(chain_list, ~ ncol(.x$prop_chain)))

    mu <- map2(chain_list, burnin_list, ~
                   map_dfc(.x$mu_chains, ~ colMeans(.x[-.y, ]), .y = .y) %>%
                   as.matrix %>%
                   unname) %>%
        do.call(what = cbind, args = .)
    prop <- map2(chain_list, burnin_list, ~ colMeans(.x$prop_chain[-.y, ])) %>%
        unlist %>%
        unname
    sig_ests <-
        map2(chain_list, burnin_list, ~
                 lapply(.x$Sigma_chains, function(X)
                     apply(X[, , -.y], c(1, 2), mean)) %>%
                 simplify2array) %>%
        abind::abind(along = 3)

    z_chain <- map2(chain_list, burnin_list, ~ .x$z_chain[, -.y])
    for(i in 2:length(newclustlabs)) {
        zold <- 1
        for(z in (newclustlabs[i-1]+1):newclustlabs[i]) {
            z_chain[[i]][z_chain[[i]] == zold] <- z
            zold <- zold + 1
        }
    }
    rles <- map(z_chain, ~ apply(.x, 1, function(X) rle(sort(X))))
    z <- map(rles, ~ map_int(.x, ~ .x$values[which.max(.x$lengths)])) %>%
        unlist

    nm <- ncol(mu)
    dm <- nrow(mu)

    cluster_dist <- matrix(0, nm, nm)
    combos <- combn(unique(z), 2)
    for (i in 1:ncol(combos)) {
        cluster_dist[combos[1, i], combos[2, i]] <-
            cluster_dist[combos[2, i], combos[1, i]] <-
            get_KL_distance(
                mu1 = mu[, combos[1, i]],
                mu2 =  mu[, combos[2, i]],
                Sigma1 = sig_ests[, , combos[1, i]],
                Sigma2 = sig_ests[, , combos[2, i]]
            )
    }
    rownames(cluster_dist) <- colnames(cluster_dist) <- 1:nm

    rmidx <- rowSums(cluster_dist) == 0

    if(sum(!rmidx) < n_groups) {
      n_groups <- sum(!rmidx)
      warning(paste0("n_groups should be less than or equal to the number of non-empty clusters.\nYou have ",
                  sum(!rmidx), " non-empty clusters, so merging with that number instead."))
    }

    cluster_dist <- cluster_dist[!rmidx,!rmidx]
    cl <- hclust(as.dist(cluster_dist), method = "complete")

    merge_idx <- cutree(cl, k = n_groups)
    merge_prop <- rep(0, length(unique(merge_idx)))
    merge_mu <- matrix(NA, length(unique(merge_idx)), dm)
    
    merge_sigma <- array(NA, dim = c(dm, dm, length(unique(merge_idx))))
    
    
    for (i in sort(unique(merge_idx))) {
        subidx <- as.numeric(names(merge_idx[merge_idx == i]))
        sub_prop <- prop[subidx]
        sub_mu <- mu[, subidx]
        
        sub_sigma <- sig_ests[,,subidx]
        
        merge_prop[i] <- sum(sub_prop)

        if (length(sub_prop) == 1) {
            merge_mu[i, ] <- sub_mu
            
            merge_sigma[,,i] <- sub_sigma
        } else {
            merge_mu[i, ] <- sub_mu %*% (sub_prop / merge_prop[i])
            
            for(j in seq_along(sub_prop)) {
                sub_sigma[,,j] <- sub_sigma[,,j] * (sub_prop[j] / merge_prop[i])
            }
            merge_sigma[,,i] <- apply(sub_sigma, c(1,2), sum)
        }
    }
    merge_prop <- merge_prop / sum(merge_prop)

    outz <- z
    for (i in sort(unique(merge_idx))) {
        outz[z %in% as.numeric(names(merge_idx[merge_idx == i]))] <- i
    }

    list(
        "merged_z" = outz,
        "merged_mu" = merge_mu,
        "merged_sigma" = merge_sigma,
        "merged_prop" = merge_prop,
        "clustering" = cl
    )
}


# Clustering the columns
compute_distances_between_conditions_multichain <- function(chain_list, burnin_list) {
    if(typeof(burnin_list) != "list") {
        burnin_list <- map(1:length(chain_list), ~burnin_list)
    }

    n_vec <- map_int(chain_list, ~ nrow(.x$z_chain))
    sample_prop_vec <- n_vec / sum(n_vec)

    weighted_corrs <- purrr::pmap(.l = list(chain_list, burnin_list, sample_prop_vec), ~ {
        # Get parameter estimates after burn-in
        prop <- colMeans(..1$prop_chain[-(..2), ]) * ..3

        sig_ests <- map(..1$Sigma_chains, function(xx) apply(xx[, , -(..2)], c(1, 2), mean)) %>%
            simplify2array

        # get correlations from covariances
        correlations <- array(apply(sig_ests, 3, cov2cor), dim(sig_ests))

        # weight correlations by mixing weights
        imap(prop, ~ .x * correlations[, , .y]) %>%
            Reduce(`+`, x = .)
    }) %>%
        Reduce(`+`, x = .)

    # Convert to distance
    return(sqrt(1-(weighted_corrs) ^ 2))
}


# Clustering the rows
compute_distances_between_clusters_multichain <- function(chain_list, burnin_list) {
    if(typeof(burnin_list) != "list") {
        burnin_list <- map(1:length(chain_list), ~ burnin_list)
    }

    newclustlabs <- cumsum(map_int(chain_list, ~ ncol(.x$prop_chain)))

    # Get parameter estimates after burn-in
    mu <- map2(chain_list, burnin_list, ~
                   map_dfc(.x$mu_chains, ~ colMeans(.x[-.y, ]), .y = .y) %>%
                   as.matrix %>%
                   unname) %>%
        do.call(what = cbind, args = .)
    prop <- map2(chain_list, burnin_list, ~ colMeans(.x$prop_chain[-.y, ])) %>%
        unlist %>%
        unname
    z_chain <- map2(chain_list, burnin_list, ~ .x$z_chain[, -.y])

    for(i in 2:length(newclustlabs)) {
        zold <- 1
        for(z in (newclustlabs[i-1]+1):newclustlabs[i]) {
            z_chain[[i]][z_chain[[i]] == zold] <- z
            zold <- zold + 1
        }
    }
    rles <- map(z_chain, ~ apply(.x, 1, function(X) rle(sort(X))))
    z <- map(rles, ~ map_int(.x, ~ .x$values[which.max(.x$lengths)])) %>%
        unlist

    sig_ests <-
        map2(chain_list, burnin_list, ~
                 lapply(.x$Sigma_chains, function(X)
                     apply(X[, , -.y], c(1, 2), mean)) %>%
                 simplify2array) %>%
        abind::abind(along = 3)

    nm <- ncol(mu)
    dm <- nrow(mu)


    cluster_dist <- matrix(0, nm, nm)
    combos <- combn(unique(z), 2)
    for (i in 1:ncol(combos)) {
        cluster_dist[combos[1, i], combos[2, i]] <-
            cluster_dist[combos[2, i], combos[1, i]] <-
            get_KL_distance(
                mu1 = mu[, combos[1, i]],
                mu2 =  mu[, combos[2, i]],
                Sigma1 = sig_ests[, , combos[1, i]],
                Sigma2 = sig_ests[, , combos[2, i]]
            )
    }
    rmidx <- rowSums(cluster_dist) == 0
    cluster_dist <- cluster_dist[!rmidx, !rmidx]

    labels <- apply(sign(mu[,!rmidx]), 2, function(X)
        paste0("(", paste0(X, collapse = ","), ")"))

    rownames(cluster_dist) <- colnames(cluster_dist) <- labels

    return(cluster_dist)
}

# Get row reordering (for bi-clustering heatmaps)
get_row_reordering_multichain <- function(row_clustering, chain_list, burnin_list, dat_list) {
   if(typeof(burnin_list) != "list") {
       burnin_list <- map(1:length(chain_list), ~burnin_list)
   }

  dat <- as.matrix(dplyr::bind_rows(dat_list))

  # Get MAP class labels based on posterior samples
  # rename cluster labels after aggregatng
  z_list <- map2(chain_list, burnin_list, ~ get_MAP_z(.x, .y))
  newclustlabs <- cumsum(map_int(chain_list, ~ ncol(.x$prop_chain)))
  for(i in 2:length(newclustlabs)) {
      # zold_count <- 1
      # zold <- sort(row_clustering$order)[zold_count]
      zold <- 1
      for(z in (newclustlabs[i-1]+1):newclustlabs[i]) {
          z_list[[i]][z_list[[i]] == zold] <- z
          zold <- zold + 1
      }
  }



  # cluster number
  nm <- length(row_clustering$order)

  # Get a row label based on row clustering
  z <- unlist(z_list)
  rm(z_list)

  count <- 1
  for(znew in sort(unique(z))) {
      z[z == znew] <- count
      count <- count + 1
  }

  newz <- z
  count <- 1
  for (i in seq(nm)) {
    idx <- z == row_clustering$order[i]
    if (sum(idx) > 0) {
      if (sum(idx) > 1) {
        reord <-
          hclust(as.dist(sqrt(1 - cor(
            t(dat[idx,]), method = "pearson"
          ) ^ 2)))$order
        newz[idx] <- (count:(count + sum(idx) - 1))[reord]
      } else {
        newz[idx] <- count:(count + sum(idx) - 1)
      }
      count <- count + sum(idx)
    }
  }
  newz
}


# test_consistency_multichain <- function(chain,
#                              burnin,
#                              u = NULL,
#                              b = 0.5,
#                              with_zero = FALSE,
#                              agnostic_to_sign = FALSE) {
#   z_chain <- chain$z_chain[,-burnin]
#   labs <- map(chain$mu_chains, ~ sign(.x[1,])) %>%
#     do.call(`rbind`, .)
#   n <- nrow(z_chain)
#   D <- ncol(labs)
#
#   if(with_zero & !is.null(u)) {
#     warning("Testing consistency and including the null class, but a threshold u was specified. Ignoring u and setting u = D.")
#   }
#
#   if(!is.null(u)) {
#     if(u > D) {
#       stop("u cannot be greater than the dimension of the data.")
#     }
#   }
#
#
#   if(agnostic_to_sign) {
#     labs <- abs(labs)
#   }
#
#   # Testing for consistency across all dimensions
#   if(with_zero) {
#     u <- D
#
#     # Find labels where all entries are equal
#     consistent_lab_idx <- apply(labs, 1, function(X) diff(range(X)) == 0)
#
#     # If there are no consistent labels, there are no consistent observations
#     if(sum(consistent_lab_idx) == 0) {
#       return(rep(FALSE, n))
#     }
#
#     consistent_lab_idx <- which(consistent_lab_idx)
#   } else { # Testing for replicability of a signal
#
#     consistent_pos_lab_idx <- apply(labs, 1, function(X) sum(X == 1) >= u)
#     consistent_neg_lab_idx <- apply(labs, 1, function(X) sum(X == -1) >= u)
#
#     # If there are no consistent labels, there are no consistent observations
#     if(sum(consistent_pos_lab_idx) == 0 & sum(consistent_neg_lab_idx) == 0) {
#       return(rep(FALSE, n))
#     }
#
#     consistent_lab_idx <-
#       sort(unique(c(
#         which(consistent_neg_lab_idx),
#         which(consistent_pos_lab_idx)
#       )))
#   }
#
#   # Find probabilities that each observation is assigned a consistent label
#   consistent_prob <-
#     colMeans(apply(z_chain, 1, function(X)
#       X %in% consistent_lab_idx))
#
#   if(length(b) == 1) {
#     return(consistent_prob > b)
#   } else {
#     replicable <- map(b, ~ consistent_prob > .x) %>%
#       set_names(paste(b))
#     return(replicable)
#   }
# }
