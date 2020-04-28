# Get most likely class label based on posterior samples
get_MAP_z <- function(chain, burnin) {
  rles <- apply(chain$z_chain[,-burnin], 1, function(X) rle(sort(X)))
  z <- map_int(rles, ~ .x$values[which.max(.x$lengths)])
}

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
    sig_ests <-
        map(chain$Sigma_chains, ~ apply(.x[, , -burnin], c(1, 2), mean)) %>% simplify2array
    z <- get_MAP_z(chain, burnin)
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
  rmidx <- rowSums(cluster_dist) == 0
  cluster_dist <- cluster_dist[!rmidx, !rmidx]

  labels <- apply(sign(mu[,!rmidx]), 2, function(X)
    paste0("(", paste0(X, collapse = ","), ")"))

  rownames(cluster_dist) <- colnames(cluster_dist) <- labels

  cluster_dist
}

# Get row reordering (for bi-clustering heatmaps)
get_row_reordering <- function(row_clustering, chain, burnin) {
  # Get MAP class labels based on posterior samples
  z <- get_MAP_z(chain, burnin)

  # cluster number
  nm <- length(row_clustering$order)

  # Get a row label based on row clustering
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

test_consistency <- function(chain,
                             burnin,
                             u,
                             b = 0.5,
                             with_zero = FALSE,
                             agnostic_to_sign = FALSE) {
  z_chain <- chain$z_chain[,-burnin]
  labs <- map(chain$mu_chains, ~ sign(.x[1,])) %>%
    do.call(`rbind`, .)
  n <- nrow(z_chain)
  D <- ncol(labs)

  if(with_zero & !is.null(u)) {
    warning("Testing consistency and including the null class, but a threshold u was specified. Ignoring u and setting u = D.")
  }

  if(u > D) {
    stop("u cannot be greater than the dimension of the data.")
  }

  if(agnostic_to_sign) {
    labs <- abs(labs)
  }

  # Testing for consistency across all dimensions
  if(with_zero) {
    # Find labels where all entries are equal
    consistent_lab_idx <- apply(labs, 1, function(X) diff(range(X)) == 0)

    # If there are no consistent labels, there are no consistent observations
    if(sum(consistent_lab_idx) == 0) {
      return(rep(FALSE, n))
    }

    consistent_lab_idx <- which(consistent_lab_idx)
  } else { # Testing for replicability of a signal

    consistent_pos_lab_idx <- apply(labs, 1, function(X) sum(X == 1) >= u)
    consistent_neg_lab_idx <- apply(labs, 1, function(X) sum(X == -1) >= u)

    # If there are no consistent labels, there are no consistent observations
    if(sum(consistent_pos_lab_idx) == 0 & sum(consistent_neg_lab_idx) == 0) {
      return(rep(FALSE, n))
    }

    consistent_lab_idx <-
      sort(unique(c(
        which(consistent_neg_lab_idx),
        which(consistent_pos_lab_idx)
      )))
  }

  # Find probabilities that each observation is assigned a consistent label
  consistent_prob <-
    colMeans(apply(z_chain, 1, function(X)
      X %in% consistent_lab_idx))

  consistent_prob > b
}
