#----------------------------------------------
# Gelman-rubin convergence diagnostics for parallel mcmcs
# Calculate PSRFs for all sampled params (Geman-Rubin convergence diagnostics)
calc_PSRF <- function(ests_list, vars_list, M, N) {
    grandmean <- Reduce(`+`, ests_list)/M
    B <- purrr::map(ests_list, ~ (.x - grandmean)^2) %>%
        Reduce(`+`, .) %>%
        `*` ((M+1) / (M * (M-1)))
    W <- purrr::map2(vars_list, N, ~ (.y - 1) / .y * .x) %>%
        Reduce(`+`, .) %>%
        `/` (M)
    Vhat <- W + B
    PSRF <- (Vhat / W) %>%
        replace(is.nan(.), NA)
    
    return(PSRF)
}

calculate_gelmanRubin <- function(chain_list, burnin_list) {
    M <- length(chain_list)
    if(!is.list (burnin_list)) {
        burnin_list <- list(burnin_list)
    }
    
    if(length(burnin_list) == 1) {
        burnin_list <- rep(burnin_list, M)
    }
    
    chain_list <- purrr::map2(chain_list, burnin_list, ~ {
        .x$mu_chains <-
            purrr::map(.x$mu_chains, function(.chain, .burnin)
                .chain[-.burnin, ], .burnin = .y)
        .x$Sigma_chains <-
            purrr::map(.x$Sigma_chains, function(.chain, .burnin)
                .chain[, , -.burnin], .burnin = .y)
        .x$prop_chain <- .x$prop_chain[-.y,]
        .x$z_chain <- .x$z_chain[, -.y]
        return(.x)
    })
    
    # number of post burn-in iterations per chain
    N <- purrr::map_int(chain_list, ~nrow(.x$prop_chain))
    
    # Get parameter estimates after burn-in for each chain
    mu_ests_list <- purrr::map(chain_list, "mu_chains") %>%
        purrr::map(~ purrr::map(.x, colMeans) %>%
                       do.call(what = rbind, args = .))
    Sigma_ests_list <- purrr::map(chain_list, "Sigma_chains") %>%
        purrr::map(~ purrr::map(.x, ~apply(.x, c(1,2), mean)) %>%
                       simplify2array())
    prop_ests_list <- purrr::map(chain_list, "prop_chain") %>%
        purrr::map(colMeans)
    
    # Get parameter estimate variance after burn-in for each chain
    mu_vars_list <- purrr::map(chain_list, "mu_chains") %>%
        purrr::map(~ purrr::map(.x, ~apply(.x, 2, var)) %>%
                       do.call(what = rbind, args = .))
    Sigma_vars_list <- purrr::map(chain_list, "Sigma_chains") %>%
        purrr::map(~ purrr::map(.x, ~apply(.x, c(1,2), var)) %>%
                       simplify2array())
    prop_vars_list <- purrr::map(chain_list, "prop_chain") %>%
        purrr::map(~apply(.x, 2, var))
    
    #---------------------------------------------------------------------------
    # Calculate PSRFs for all sampled params
    # Returns NAs for dimensions associated with 0 label, since these do not get sampled during MCMC
    
    
    PSRF_mu <- calc_PSRF(mu_ests_list, mu_vars_list, M, N)
    PSRF_Sigma <- calc_PSRF(Sigma_ests_list, Sigma_vars_list, M, N)
    PSRF_prop <- calc_PSRF(prop_ests_list, prop_vars_list, M, N)
    
    list(
        PSRF_mu = PSRF_mu,
        PSRF_Sigma = PSRF_Sigma,
        PSRF_prop = PSRF_prop
    )
}