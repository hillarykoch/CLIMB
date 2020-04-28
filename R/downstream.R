# Computes pairwise KL distance between 2 MVNs
get_KL_distance <- function(mu1, mu2, Sigma1, Sigma2) {
    k <- nrow(Sigma1)
    KL1 <- 0.5 * (sum(diag(solve(Sigma1) %*% Sigma2)) - k + log(det(Sigma1) / det(Sigma2)))
    KL2 <- 0.5 * (sum(diag(solve(Sigma2) %*% Sigma1)) - k + log(det(Sigma2) / det(Sigma1)))

    0.5 * (KL1 + KL2)
}

merge_classes <- function(n_groups, chain, burnin) {
    mu <- map(chain$mu_chains, ~ colMeans(.x[-burnin, ])) %>%
        bind_cols() %>%
        as.matrix %>%
        unname
    prop <- colMeans(chain$prop_chain[-burnin, ])
    z_chain <- chain$z_chain[, -burnin]
    sig_ests <-
        map(chain$Sigma_chains, ~ apply(.x[, , -burnin], c(1, 2), mean)) %>% simplify2array
    nm <- ncol(mu)
    dm <- nrow(mu)

    rles <- apply(z_chain, 1, function(X) rle(sort(X)))
    z <- map_int(rles, ~ .x$values[which.max(.x$lengths)])

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
                  sum(!rmidx), " empty clusters, so merging with that number instead."))
    }

    cluster_dist <- cluster_dist[!rmidx,!rmidx]
    cl <- hclust(as.dist(cluster_dist), method = "complete")

    merge_idx <- cutree(cl, k = n_groups)
    merge_prop <- rep(0, length(unique(merge_idx)))
    merge_mu <- matrix(NA, length(unique(merge_idx)), dm)
    for (i in sort(unique(merge_idx))) {
        subidx <- as.numeric(names(merge_idx[merge_idx == i]))
        sub_prop <- prop[subidx]
        sub_mu <- mu[, subidx]
        merge_prop[i] <- sum(sub_prop)

        if (length(sub_prop) == 1) {
            merge_mu[i, ] <- sub_mu
        } else {
            merge_mu[i, ] <- sub_mu %*% (sub_prop / merge_prop[i])
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
        "merged_prop" = merge_prop,
        "clustering" = cl
    )
}


# Clustering the columns
compute_distances_between_conditions <- function(chain, burnin) {
  # Get parameter estimates after burn-in
  prop <- colMeans(chain$prop_chain[-burnin, ])
  sig_ests <-
    map(chain$Sigma_chains, ~ apply(.x[, , -burnin], c(1, 2), mean)) %>%
    simplify2array

  # get correlations from covariances
  correlations <- array(apply(sig_ests, 3, cov2cor), dim(sig_ests))

  # weight correlations by mixing weights
  weighted_corrs <- imap(prop, ~ .x * correlations[, , .y]) %>%
    Reduce(`+`, x = .)

  # Convert to distance
  sqrt(1-(weighted_corrs) ^ 2)
}


# Clustering the rows
compute_distances_between_clusters <- function(chain, burnin) {
  # Get parameter estimates after burn-in
  mu <- map(chain$mu_chains, ~ colMeans(.x[-burnin, ])) %>%
    bind_cols() %>%
    as.matrix %>%
    unname
  prop <- colMeans(chain$prop_chain[-burnin, ])
  z_chain <- chain$z_chain[, -burnin]
  sig_ests <-
    map(chain$Sigma_chains, ~ apply(.x[, , -burnin], c(1, 2), mean)) %>%
    simplify2array
  nm <- ncol(mu)
  dm <- nrow(mu)

  rles <- apply(z_chain, 1, function(X) rle(sort(X)))
  z <- map_int(rles, ~ .x$values[which.max(.x$lengths)])


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


  cluster_dist[!rmidx, !rmidx]
}
