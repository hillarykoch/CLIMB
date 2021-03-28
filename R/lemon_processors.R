# Functions needed to process fconstr_pGMCM to then dump into an LGF file

# Get the 2-mers from the fits
get_h <- function(fits) {
    mus <- purrr::map(fits, "mu") %>%
        lapply(function(X)
            replace(X, X > 0, 1)) %>%
        lapply(function(X)
            replace(X, X < 0,-1))
    mus
}

# Find out which associations exist in each dimension
# (may be all, but if filterable, can make problem easier)
filter_h <- function(h, d) {
  spl <- strsplit(names(h), split = "_") %>%
    purrr::map(as.numeric)

  filt <- list()
  for (i in seq(d)) {
    lapply(seq_along(spl), function(X) {
      idx <- spl[[X]] == i
      if (any(idx)) {
        h[[X]][, idx]
      }
    }) %>%
      unlist() %>%
      unique() %>%
      sort() -> filt[[i]]
  }
  names(filt) <- seq(d)
  filt
}

# Extract which results correspond to consecutive dimensions (e.g. 1_2, 2_3, but not 1_3)
get_consecutive <- function(h, non_consec = FALSE) {
  nms <- names(h)
  consec <- strsplit(names(h), split = "_") %>%
    purrr::map(as.numeric) %>%
    map(diff) %>%
    `==`(1)
  if (!non_consec) {
    h[consec]
  } else {
    h[!consec]
  }
}

# Make LGF file to be read by LEMON DigraphReader
write_LGF <- function(h, d, path) {
    cat("Writing LGF file...")
    
    # Make node section of LGF file (first pass)
    filt_h <- filter_h(h, d)
    l <- length(unlist(filt_h))
    dims <- purrr::map(filt_h, length)
    node <-
        data.frame(
            "label" = c(0, seq_along(unlist(filt_h)), l + 1),
            "dim" = c(0, lapply(seq(dims), function(X)
                rep(X, times = dims[[X]])) %>% unlist, d + 1),
            "assoc" = c(0, unlist(filt_h), 0)
        )
    
    # Make arc section of LGF file (first pass)
    consec <- get_consecutive(h)
    
    src <- trg <- list()
    for (i in seq_along(consec)) {
        src[[i]] <-
            sapply(consec[[i]][, 1], function(X)
                node$label[node$assoc == X & node$dim == i])
        trg[[i]] <-
            sapply(consec[[i]][, 2], function(X)
                node$label[node$assoc == X & node$dim == i + 1])
    }
    
    sources <- dplyr::filter(node, dim == 1)$label
    targets <- dplyr::filter(node, dim == d)$label
    arcs <- data.frame(
        "source" = c(rep(0, length(sources)),
                     abind::abind(src),
                     targets),
        "target" = c(sources,
                     abind::abind(trg),
                     rep(l + 1, length(targets)))
    )
    
    # Recursively eliminate nodes that don't have sources however many arcs back
    keepidx <- c(0, which(seq(l) %in% arcs$target & seq(l) %in% arcs$source), l+1)
    src_bind <- abind::abind(src)
    trg_bind <- abind::abind(trg)
    
    while(length(keepidx) != nrow(node)) {
        # Filter out nodes that we don't need to keep
        node <- node[node$label %in% keepidx, ]
        
        # make the arcs table based on filtered nodes
        sources <- node$label[node$dim == 1][
            node[node$dim == 1,]$label %in% arcs$source &
                node[node$dim == 1,]$label %in% arcs$target]
        targets <- node$label[node$dim == d][
            node[node$dim == d,]$label %in% arcs$source &
                node[node$dim == d,]$label %in% arcs$target]
        
        arckeepidx <- which(src_bind %in% keepidx & trg_bind %in% keepidx)
        src_bind <- src_bind[arckeepidx]
        trg_bind <- trg_bind[arckeepidx]
        
        arcs <- data.frame(
            "source" = c(
                rep(0, length(sources)),
                src_bind,
                targets
            ),
            "target" = c(
                sources,
                trg_bind,
                rep(l + 1, length(targets))
            )
        )
        
        keepidx <- c(0, which(seq(l) %in% arcs$target & seq(l) %in% arcs$source), l+1)
    }
    
    # Write the LGF file
    readr::write_tsv(data.frame("@nodes"),
                     file = path,
                     col_names = FALSE)
    readr::write_tsv(node,
                     file = path,
                     col_names = TRUE,
                     append = TRUE)
    cat("\n", file = path, append = TRUE)
    
    readr::write_tsv(
        data.frame("@arcs"),
        file = path,
        col_names = FALSE,
        append = TRUE
    )
    cat("\t\t -\n", file = path, append = TRUE)
    readr::write_tsv(
        arcs,
        file = path,
        col_names = FALSE,
        append = TRUE,
        na = ""
    )
    cat("\n", file = path, append = TRUE)
    readr::write_tsv(
        data.frame("@attributes"),
        file = path,
        col_names = FALSE,
        append = TRUE
    )
    attrib <- data.frame("type" = c("source", "target"),
                         "label" = c(0, length(unlist(filt_h)) + 1))
    readr::write_tsv(attrib,
                     file = path,
                     col_names = FALSE,
                     append = TRUE)
    cat("done!\n")
}

split_LGF <- function(path) {
    cat("Splitting LGF file into 2...")
    lgf_temp <- read_tsv(path, skip = 1)
    lgf_nodes <- dplyr::slice(lgf_temp, 1:(which(lgf_temp$label == "@arcs") - 1))
    lgf_arcs <- dplyr::slice(lgf_temp, (which(lgf_temp$label == "@arcs") + 2):(which(lgf_temp$label == "@attributes") - 1)) %>%
        setNames(c("start", "end", "dash")) %>%
        dplyr::select(1:2)
    lgf_attributes <- dplyr::slice(lgf_temp, (n()-1):n())
    
    # up to this layer will be included in the first run
    split_layer <- ceiling((max(lgf_nodes$dim)-1) / 2)
    
    
    # two node sets
    lgf_nodes1 <- dplyr::filter(lgf_nodes, dim <= split_layer) %>%
        bind_rows(tibble(label = as.character(max(as.numeric(.$label)) + 1), dim = split_layer + 1, assoc = "0"))
    lgf_nodes2 <- dplyr::filter(lgf_nodes, dim > split_layer) %>%
        bind_rows(tibble(
            label = as.character(min(as.numeric(.$label))-1),
            dim = split_layer,
            assoc = "0")) %>%
        arrange(label) %>%
        mutate(dim = dim - split_layer)
    
    # two arc sets
    lgf_arcs1 <- dplyr::filter(lgf_arcs, end < max(as.numeric(lgf_nodes1$label))) %>%
        bind_rows(
            tibble(
                start = dplyr::filter(lgf_nodes1, dim == split_layer)$label,
                end = max(as.numeric(lgf_nodes1$label))
            ))
    
    lgf_arcs2 <- bind_rows(
        tibble(start = as.numeric(lgf_nodes2$label)[2] - 1,
               end = as.numeric(dplyr::filter(lgf_nodes2, dim == 1)$label)),
        dplyr::filter(lgf_arcs, as.numeric(start) >= as.numeric(lgf_nodes2$label)[2]) %>%
            mutate(across(where(is.character), .fn = as.numeric))
    )
    
    
    
    # Write the LGF file
    if(str_detect(path, "\\.[[:alpha:]]*")) {
        path1 <- str_replace(path, "\\.", "1\\.")
        path2 <- str_replace(path, "\\.", "2\\.")
    } else {
        path1 <- paste0(path, "1")
        path2 <- paste0(path, "2")
    }
    
    readr::write_tsv(data.frame("@nodes"),
                     file = path1,
                     col_names = FALSE)
    readr::write_tsv(data.frame("@nodes"),
                     file = path2,
                     col_names = FALSE)
    readr::write_tsv(lgf_nodes1,
                     file = path1,
                     col_names = TRUE,
                     append = TRUE)
    readr::write_tsv(lgf_nodes2,
                     file = path2,
                     col_names = TRUE,
                     append = TRUE)
    cat("\n", file = path1, append = TRUE)
    cat("\n", file = path2, append = TRUE)
    
    readr::write_tsv(
        data.frame("@arcs"),
        file = path1,
        col_names = FALSE,
        append = TRUE
    )
    readr::write_tsv(
        data.frame("@arcs"),
        file = path2,
        col_names = FALSE,
        append = TRUE
    )
    
    cat("\t\t -\n", file = path1, append = TRUE)
    cat("\t\t -\n", file = path2, append = TRUE)
    
    readr::write_tsv(
        lgf_arcs1,
        file = path1,
        col_names = FALSE,
        append = TRUE,
        na = ""
    )
    readr::write_tsv(
        lgf_arcs2,
        file = path2,
        col_names = FALSE,
        append = TRUE,
        na = ""
    )
    
    cat("\n", file = path1, append = TRUE)
    cat("\n", file = path2, append = TRUE)
    
    readr::write_tsv(
        data.frame("@attributes"),
        file = path1,
        col_names = FALSE,
        append = TRUE
    )
    
    readr::write_tsv(
        data.frame("@attributes"),
        file = path2,
        col_names = FALSE,
        append = TRUE
    )
    
    attrib1 <- data.frame("type" = c("source", "target"),
                          "label" = c(0, as.numeric(tail(lgf_nodes1$label, 1))))
    attrib2 <- data.frame("type" = c("source", "target"),
                          "label" = c(as.numeric(lgf_nodes2$label[1]), as.numeric(tail(lgf_nodes2$label, 1))))
    
    readr::write_tsv(attrib1,
                     file = path1,
                     col_names = FALSE,
                     append = TRUE)
    readr::write_tsv(attrib2,
                     file = path2,
                     col_names = FALSE,
                     append = TRUE)
    cat("done!\n")
    
    split_layer
}

# get paths using LEMON, then convert to latent association vectors
get_paths <- function(filepath) {
    cat("Finding latent classes...")
    path_build <- cgetPaths(filepath = filepath)
    cat("done!\n")
    path_build
}

# Prune paths that are discordant with "non-consecutive" pairwise estimates
prune_paths <- function(h, assoc_mx, split_layer = 0) {
    nonconsec <- get_consecutive(h, non_consec = TRUE)
    labs <- names(nonconsec)
    keepers <-
        matrix(0, nrow = nrow(assoc_mx), ncol = length(nonconsec))

    for (i in seq_along(nonconsec)) {
        pair <- strsplit(labs[i], split = "_") %>%
            `[[` (1) %>%
            as.numeric %>%
            `-` (split_layer)
        keepers[, i] <-
            crowMatch(assoc_mx[, pair], nonconsec[[i]])[, 1]
    }

    # If a row is ever not a keeper (contains at least one 0), remove it
    prunes <- apply(keepers, MARGIN = 1, function(X)
        any(X == 0))
    assoc_mx[!prunes, ]
}

combine_prunes <- function(prune1, prune2, h, split_layer) {
    join_pair <- h[paste0(split_layer, "_", split_layer + 1)][[1]]
    
    map_dfr(1:nrow(join_pair), ~ {
        l_candidates <- prune1[prune1[,split_layer] == join_pair[.x,1],]
        r_candidates <- prune2[prune2[,1] == join_pair[.x,2],]
        suppressMessages(tidyr::expand_grid(
            as_tibble(l_candidates),
            as_tibble(r_candidates),
            .name_repair = "unique"
        ))
    }) %>%
        as.matrix
} 

# Put everything together in one function here, get_reduced_classes
get_reduced_classes <- function(fits, d, filepath = "lgf.txt", split_in_two = TRUE) {
    h <- get_h(fits)
    filt <- filter_h(h, d)
    write_LGF(h, d, filepath)
    
    if(!split_in_two) {
        paths <- get_paths(filepath)
        assoc <- cassociate(paths, filepath, length(unlist(filt)))
        prune_paths(h, assoc) # can probably replace this with C code, cprune_paths
    } else {
        split_layer <- suppressWarnings(split_LGF(filepath))
        
        # Write the split LGF files
        if(str_detect(filepath, "\\.[[:alpha:]]*")) {
            filepath1 <- str_replace(filepath, "\\.", "1\\.")
            filepath2 <- str_replace(filepath, "\\.", "2\\.")
        } else {
            filepath1 <- paste0(filepath, "1")
            filepath2 <- paste0(filepath, "2")
        }
        
        paths1 <- get_paths(filepath1)
        paths2 <- get_paths(filepath2)
        
        assoc1 <- cassociate(paths1, filepath1, length(unlist(filt[1:split_layer])))
        assoc2 <- cassociate(paths2, filepath2, length(unlist(filt[(split_layer):length(filt)])))
        
        rmidx1 <- sapply(as.character((split_layer+1):d), function(X) stringr::str_detect(names(h), pattern = paste0("^", X, "_|_", X, "$"))) %>%
            apply(1, any)
        rmidx2 <- sapply(as.character(seq(split_layer)), function(X) stringr::str_detect(names(h), pattern = paste0("^", X, "_|_", X, "$"))) %>%
            apply(1, any)
        prune1 <- prune_paths(h[!rmidx1], assoc1, split_layer = 0)
        prune2 <- prune_paths(h[!rmidx2], assoc2, split_layer = split_layer)
        
        prune_final <- prune_paths(h,
                                   combine_prunes(prune1, prune2, h, split_layer = split_layer),
                                   split_layer = 0)
        colnames(prune_final) <- NULL
        prune_final
    }
}

# Among all reduced classes, get the indices which correspond to the truth
# For simulated data only
get_true_assoc_idx <- function(red_class, true_assoc) {
  true_assoc <- as.matrix(true_assoc)
  trueidx <- cget_true_assoc_idx(red_class, true_assoc)
  keepidx <- trueidx != 0
  sort(trueidx[keepidx])
}
