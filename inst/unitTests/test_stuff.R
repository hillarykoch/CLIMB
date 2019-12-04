test_mahalanobis <- function() {
    n <- 1500
    prop <- c(1 / 3, 1 / 3, 1 / 3)
    mu <- list(c(1, 1), c(-1, -1), c(1 - sqrt(2), 3 + sqrt(2)))
    sigma <- list(matrix(c(1, .75, .75, 1), nrow = 2, byrow = T),
                  matrix(c(2, -.8, -.8, 1), nrow = 2, byrow = T),
                  diag(1:2))

    sim <- rGMCM(n, prop, mu, sigma)$data
    RUnit::checkEqualsNumeric(mahalanobis(sim, mu[[1]], sigma[[1]]),
                              Mahalanobis(as.matrix(sim), mu[[1]], sigma[[1]]))
}

test_dmvnorm <- function() {
    n <- 1500
    prop <- c(1 / 3, 1 / 3, 1 / 3)
    mu <- list(c(1, 1), c(-1, -1), c(1 - sqrt(2), 3 + sqrt(2)))
    sigma <- list(matrix(c(1, .75, .75, 1), nrow = 2, byrow = T),
                  matrix(c(2, -.8, -.8, 1), nrow = 2, byrow = T),
                  diag(1:2))

    sim <- rGMCM(n, prop, mu, sigma)$data

    RUnit::checkEqualsNumeric(dmvnorm(sim, mu[[1]], sigma[[1]]),
                              cdmvnorm(as.matrix(sim), mu[[1]], sigma[[1]]))
}

# R code for E step
R_Estep_pGMM <- function(dat, prop, mu, sigma) {
    num <-
        lapply(seq(mu), function(X)
            dmvnorm(dat, mu[[X]], sigma[[X]]) * prop[X]) %>%
        abind(along = 2)
    denom <- rowSums(num)
    h_est <- num / denom
}

test_get_prior_count <- function(fits, red_class, d, n) {
    # Compute prior with R code
    spl <- strsplit(names(fits), "_") %>% map(as.numeric)
    prop_path <- rep(0, nrow(red_class))
    combos <- map(fits, "mu")
    combos <-
        lapply(combos, function(X)
            replace(X, X > 0, 1) %>% replace(X < 0,-1))

    for (i in seq(nrow(red_class))) {
        tags <- list()
        for (j in seq_along(spl)) {
            idx <-
                apply(combos[[j]], 1, function(X)
                    all(X == red_class[i, spl[[j]]]))
            tags[[j]] <- fits[[j]]$cluster == which(idx)
        }
        prop_path[i] <-
            mean(abind(tags, along = 2) %>% apply(1, all))
    }
    prop_path <- prop_path / sum(prop_path)

    # Compute with C++ function
    cprior_prop <- get_prior_prop(red_class, fits, d, n)

    RUnit::checkEqualsNumeric(prop_path, cprior_prop)
}

test_get_true_assoc_idx <- function(red_class, true_assoc) {
    true_assoc <- as.matrix(true_assoc)
    trueidx <- rep(NA, nrow(true_assoc))
    for (i in seq(nrow(true_assoc))) {
        tst <- rep(NA, nrow(red_class))
        for (j in 1:nrow(red_class)) {
            tst[j] <- all(true_assoc[i, ] == red_class[j, ])
        }
        trueidx[i] <- which(tst)
    }

    Ridx <- sort(trueidx)
    Cidx <- sort(cget_true_assoc_idx(red_class, true_assoc))

    RUnit::checkEqualsNumeric(Ridx, Cidx)
}
