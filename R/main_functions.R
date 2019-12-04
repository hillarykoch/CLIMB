# generate data for general GMM
rGMM <- function(n, prop, mu, sigma) {
    # n: sample size
    # prop: cluster mixing proportions
    # k: number of clusters
    # d: dimension of data (number of replicates)
    # mu: list of mean vectors of clusters
    # sigma: list of variance covariance matrices of clusters

    prop <- prop / sum(prop)
    k <- length(prop)
    d <- ifelse(length(prop) == 1, nrow(sigma), nrow(sigma[[1]]))
    rles <- sample(seq_along(prop), n, replace = TRUE, prob = prop) %>%
        sort %>%
        rle
    num <- rep(0, k)
    num[rles$values] <- rles$lengths

    tag <- rep(1:k, times = num)
    y <- lapply(seq_along(prop),
                function(X)
                    if(num[X] > 0) {
                        rmvnorm(num[X], mean = mu[[X]], sigma = sigma[[X]])    
                    }) %>%
        abind(along = 1) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X)
            paste0("y.", X)))

    list(
        "data" = y,
        "cluster" = tag
    )
}

# generate data for constrained version of GMCM
rconstr_GMCM <- function(n, prop, mu, sigma, rho, d, combos = NULL) {
    # n: sample size
    # prop: mixing proportion of each cluster
    # d: dimension of data (number of replicates)
    # mu: mean of "reproducible" components
    # sigma: variance of "reproducible" components
    # rho: correlation between replicates of same component
    # d: number of replicates
    # k: number of clusters

    if (is.null(combos)) {
        # Generate all combinations of replication, given in any order
        combos <- expand.grid(rep(list(-1:1), d)) %>% as.matrix
    }

    k <- nrow(combos)
    if (k != length(prop)) {
        stop("length(prop) must be equal to total number of clusters.")
    }

    sigma <- c(sigma, 1, sigma)
    mu <- c(-mu, 0, mu)
    prop <- prop / sum(prop)

    if (sum(prop == 0) > 0) {
        keepidx <- prop != 0
        combos <- combos[keepidx, ]
        prop <- prop[keepidx]
        k <- sum(keepidx)
    }

    rles <- sample(seq_along(prop), n, replace = TRUE, prob = prop) %>%
        sort %>%
        rle
    num <- rep(0, k)
    num[rles$values] <- rles$lengths
    tag <- rep.int(1:k, times = num)

    y <- lapply(seq_along(prop),
                function(X)
                    if(num[X] > 0) {
                        rmvnorm(
                            num[X],
                            mean = mu[combos[X, ] + 2],
                            sigma = get_constr_sigma(diag(sigma[combos[X, ] +
                                                                    2]),
                                                     rho,
                                                     combos[X, ])
                    )}) %>%
        abind(along = 1) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X)
            paste0("y.", X)))

    list(
        "data" = y,
        "cluster" = tag
    )
}


r_assoc <- function(d, num_patterns) {
    combos <- matrix(rep(NA, d * num_patterns),
                     nrow = num_patterns,
                     ncol = d)

    # Always have the all 0 pattern
    combos[1, ] <- rep(0, d)

    count <- 2
    while (any(is.na(combos[num_patterns, ]))) {
        assoc <- sample(-1:1, replace = TRUE, size = d)
        if (!any(apply(combos, 1, function(X)
            all(X == assoc)), na.rm = TRUE)) {
            combos[count,] <- assoc
            count <- count + 1
        }
    }
  combos
}