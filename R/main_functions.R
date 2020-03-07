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

# generate data for second constrained version of GMCM
rconstr0_GMM <-
    function(n, prop, mu, sigma, rho, d, combos = NULL) {
        # n: sample size
        # prop: mixing proportion of each cluster
        # d: dimension of data (number of replicates)
        # mu: mean of "reproducible" components
        #        a 2-vector for negative, positive association
        # sigma: variance of "reproducible" components,
        #        a 2-vector for negative, positive association
        # rho: correlation between replicates of same component
        #        a 3-vector for negative, cross, positive association
        # d: number of replicates
        # k: number of clusters
        
        # Generate all combinations of replication, given in any order
        if (is.null(combos)) {
            combos <- expand.grid(rep(list(-1:1), d)) %>% as.matrix
        }
        
        k <- nrow(combos)
        if (k != length(prop)) {
            stop("length(prop) must be equal to total number of clusters.")
        }
        
        if (length(mu) != 2) {
            stop("length(mu) must equal 2.")
        }
        if (mu[1] >= 0 | mu[2] <= 0) {
            stop("mu[1] must be < 0 and mu[2] must be >= 0.")
        }
        if (length(sigma) != 2) {
            stop("length(sigma) must equal 2.")
        }
        if (any(sigma <= 0)) {
            stop("elements of sigma must be positive.")
        }
        if (length(rho) != 3) {
            stop("length(rho) must equal 3.")
        }
        if (rho[1] < 0 | rho[3] < 0 | rho[2] > 0) {
            stop("rho[1] and rho[3] must be >= 0, rho[2] must be <= 0.")
        }
        
        sigma <- c(sigma[1], 1, sigma[2])
        mu <- c(mu[1], 0, mu[2])
        prop <- prop / sum(prop)
        
        mu_combos <- replace(combos, combos == -1, mu[1]) %>%
            replace(combos == 1, mu[3]) %>%
            as.matrix
        sig_combos <- replace(combos, combos == -1, sigma[1]) %>%
            replace(combos == 1, sigma[3]) %>%
            replace(combos == 0, 1) %>%
            as.matrix
        rho_combos <- rep(0, k) %>%
            replace(apply(combos, 1, sum) == -2, rho[1]) %>%
            replace(apply(combos, 1, sum) == 0 &
                        !apply(combos, 1, function(X)
                            any(X == 0)),
                    rho[2]) %>%
            replace(apply(combos, 1, sum) == 2, rho[3])
        
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
                                mean = mu_combos[X, ],
                                sigma =
                                    get_constr0_sigma(sig_combos[X,],
                                                      combos[X,],
                                                      rho)
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
