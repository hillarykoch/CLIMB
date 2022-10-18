get_pairwise_fits <-
    function(z,
             nlambda = 10,
             parallel = TRUE,
             ncores = 10,
             bound = 0.0,
             flex_mu = FALSE,
             ...) {
        combos <- combn(ncol(z), 2)
        n <- nrow(z)
        lambda <- sqrt(log(n)) * 10 ^ seq(-1, 0.5, length.out = nlambda)
        z <- as.data.frame(z)

        if (parallel) {
            if(ncores > ncol(combos)) {
                warning("You are requesting more cores than this analysis requires.")
            }
            cluster <- makeCluster(ncores)
            registerDoParallel(cluster)
            clusterEvalQ(cluster, library(CLIMB))
            
            if(!flex_mu) {
                fits <-
                    foreach(j = 1:ncol(combos), .packages = "foreach") %dopar% {
                        tryCatch(expr =
                        fconstr_pGMM(
                            x = z[, combos[, j]],
                            lambda = lambda,
                            tol = 1e-04,
                            itermax = 200,
                            penaltyType = "SCAD",
                            bound = bound
                        ), error = function(e) NA)
                    }    
            } else {
                fits <-
                    foreach(j = 1:ncol(combos), .packages = "foreach") %dopar% {
                        tryCatch(expr =
                        fconstr0_pGMM(
                            x = z[, combos[, j]],
                            lambda = lambda,
                            tol = 1e-04,
                            itermax = 200,
                            penaltyType = "SCAD",
                            bound = bound
                        ), error = function(e) NA)
                    }    
            }
            
        } else {
            if(!flex_mu) {
                fits <- lapply(1:ncol(combos), function(j)
                    tryCatch(expr =
                    fconstr_pGMM(
                        x = z[, combos[, j]],
                        lambda = lambda,
                        tol = 1e-04,
                        itermax = 200,
                        penaltyType = "SCAD",
                        bound = bound
                    ), error = function(e) NA))    
            } else {
                fits <- lapply(1:ncol(combos), function(j)
                    fconstr0_pGMM(
                        x = z[, combos[, j]],
                        lambda = lambda,
                        tol = 1e-04,
                        itermax = 200,
                        penaltyType = "SCAD",
                        bound = bound
                    ))
            }
        }
        names(fits) <-
            apply(combos, 2, function(X)
                paste(X, collapse =	"_"))

        fits
    }

get_prior_weights <- function(reduced_classes, fits, parallel = FALSE, ncores = 20, delta = NULL) {
    n <- length(fits[[1]]$cluster)
    d <- ncol(reduced_classes)
    if(is.null(delta)) {
        delta <- 0:choose(d,2)
    }
    reduced_classes <- as.matrix(reduced_classes)

    if(parallel) {
        cluster <- makeCluster(ncores)
        registerDoParallel(cluster)
        clusterEvalQ(cluster, library(CLIMB))

        p <- foreach(delta = delta, .packages = "foreach") %dopar% {
            get_prior_prop(reduced_classes, fits, d, n, dist_tol = delta, MAP = FALSE)
        }
    } else {
        p <- lapply(delta, function(X) get_prior_prop(reduced_classes, fits, d, n, dist_tol = X, MAP = FALSE))
    }
    p
}

get_hyperparameters <- function(z, fits, reduced_classes, prior_weights) {
    n <- length(fits[[1]]$cluster)
    D <- as.numeric(strsplit(tail(names(fits),1), "_")[[1]][2])

    labels <- purrr::map(fits, "post_prob") %>%
        lapply(function(X) apply(X, 1, function(Y) base::sample(1:ncol(X), size = 1, prob = Y)))

    rc <- dplyr::filter(reduced_classes, prior_weights * n > D)

    hyp <- get_hyperparams(fits, D, rc, z, clusters = labels, var_quantile = 0.75)
    for(i in 1:dim(hyp$Psi0)[3]) {
        idx <- diag(hyp$Psi0[,,i]) == 1
        subPsi <- hyp$Psi0[!idx, !idx, i]
        if(length(subPsi) > 1) {
            while(!(LaplacesDemon::is.positive.definite(subPsi))) {
                subPsi <- subPsi + diag(nrow(subPsi))
            }
        }
        hyp$Psi0[!idx, !idx, i] <- subPsi
    }

    alpha <- prior_weights[prior_weights * n > D] / sum(prior_weights[prior_weights * n > D])
    kappa0 <- round(n * alpha)

    list("Psi0" = hyp$Psi0, "mu0" = hyp$mu0, "alpha" = alpha, "kappa0" = kappa0)
}
