#!/usr/bin/env Rscript
# Author: Matias Becker Burgos (matias.beckerburgos@unil.ch)
# Description: This script calculates the codiversification between two phylogenetic trees heavily inspired by the Sprockett (Moeller) method.
# Used Claude to help me with the implementation.
#
# Usage:
# Rscript codiv.R [options]
#
# Required arguments:
#   --host_tree: Path to the host phylogenetic tree in Newick or Nexus format.
#   --symbiont_tree: Path to the symbiont phylogenetic tree in Newick or Nexus format.
#   --associations: Path to a tab-delimited file containing host-symbiont associations TSV (columns: host, symbiont).
#   --output: path for output TSV results file.
#
# Optional arguments:
#   --min_hosts: Min. number of distinct host per subtree (default: 3).
#   --min_symbionts: Min. number of distinct symbiont per subtree (default: 3).
#   --max_symbionts: Max symbiont tips per subtree (default: 500). Caps large root subtrees.
#   --global_perms:  Permutations for the global FDR test (default: 999).
#   --fdr_alpha:     FDR significance threshold for summary (default: 0.05).
#   --skip_fdr:      Flag. Skip permutation FDR calculation entirely.
#   --fdr_min_symbionts: Min N_Symbionts for a node to contribute to the global null max (default: 10).
#   --fdr_null_output: Path to save raw null distribution TSV (default: <output>_null_dist.tsv).
#   --span_fraction: Max. fraction of total tree span for subtree (default: 0.1).
#   --permutations: number of permutations (dfault: 99).
#   --seed: random seed for reproducibility (default: 42).
#   --cores: Number of CPU cores to use (default: 4).
#   --methods: comma-separated hommola/paco/parafit (default: hommola,paco,parafit).
#   --focus_host: Optional comma-separated list of host names (default: NULL).
#   --no_subtree_feat: If set, do not calculate subtree features (default: FALSE).
#.  --no_continue: IF set, do not continue from previous runs (default: FALSE).
#   --quiet: If set, suppress progress messages (default: FALSE).
#   --tree_format: Format of input trees (newick or nexus, default: auto).
#   --null_traversals: Number of full traversals on permuted host trees for empirical null model (default: 0, disabled).
#   --null_alpha: Significance threshold used to count 'hits' in each null traversal (default: 0.05).
#   --null_r_threshold: Minimum Hommola r (or analogous) to count a node as a hit (default: 0, any sig. node counts).
#   --skip_null: Flag. Skip null model traversal even if null_traversals > 0.
#   --occupancy: Flag. Run node-wise occupancy scan (parsimony-based phylogenetic signal test).
#   --occ_permutations: Permutations for occupancy scan (default: 999).
#   --occ_min_tips: Min symbiont tips in subtree for occupancy test (default: 5).
#   --occ_output: Path for occupancy scan output TSV (default: <output>_occupancy.tsv).
#   --rarefy_depth: Symbionts to sample per host for rarefied tests (default: 0 = auto,
#                   uses the smallest per-host count within each subtree).
#
# Example usage:
# Rscript codiv.R \
#   --host_tree host_tree.nwk \
#   --symbiont_tree symbiont_tree.nwk \
#   --associations associations.tsv \
#   --output results.tsv \
#   --cores 4 \
#   --methods hommola,paco,parafit \
#   --permutations 99 \

# Load required libraries
suppressPackageStartupMessages({
    library(ape)
    library(castor)
    ## adephylo removed — using castor::get_all_pairwise_distances (C++, much faster)
    library(TreeDist)
    library(parallel)
    library(tibble)
    library(tidyr)
    library(dplyr)
    library(vegan)
})

# Function to parse command-line arguments
parse_args <- function() {
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
        cat(readLines(con = textConnection(
            grep("^#", readLines(sys.frame(0)$ofile), value = TRUE)
            )), sep = "\n")
        quit(status = 0)
    }
    
    # default parameters
    defaults <- list(
        min_hosts = 3,
        min_symbionts = 3,
        max_symbionts = 500,
        span_fraction = 0.1,
        permutations = 99,
        seed = 42,
        cores = 4,
        methods = "hommola,paco,parafit",
        focus_host = NULL,
        no_subtree_feat = FALSE,
        no_continue = FALSE,
        quiet = FALSE,
        tree_format = "auto",
        global_perms      = 999L,
        fdr_alpha         = 0.05,
        fdr_min_symbionts = 10L,
        fdr_null_output   = NULL,
        skip_fdr          = FALSE,
        null_traversals  = 0L,
        null_alpha       = 0.05,
        null_r_threshold = 0,
        skip_null        = FALSE,
        occupancy        = FALSE,
        occ_permutations = 999L,
        occ_min_tips     = 5L,
        occ_output       = NULL,
        rarefy_depth     = 0L   ## 0 = auto (min per-host count within subtree)
    )
    # parse arguments
    i <- 1
    while (i <= length(args)) {
        key <- sub("^--", "", args[i])
        # flags without values
        if (key %in% c("no_subtree_feat", "no_continue", "quiet", "skip_fdr", "skip_null", "occupancy")) {
            defaults[[key]] <- TRUE
            i <- i + 1L
            next        
    } 
    if (i + 1L > length(args)) stop(sprintf("Flag --%s requires a value", key))
        value <- args[i + 1L]
        if (key %in% c("min_hosts", "min_symbionts", "permutations", "seed", "cores", "global_perms", "null_traversals", "occ_permutations", "occ_min_tips", "fdr_min_symbionts", "rarefy_depth")) {
            defaults[[key]] <- as.integer(value)
        } else if (key %in% c("span_fraction", "fdr_alpha", "null_alpha", "null_r_threshold")) {
            defaults[[key]] <- as.numeric(value)
        } else {
            defaults[[key]] <- value
        }
        i <- i + 2L
    }

    required <- c("host_tree", "symbiont_tree", "associations", "output")
    missing <- required[vapply(required, function(x) is.null(defaults[[x]]), logical(1))]
    if (length(missing)) stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))

    defaults
}

# Helper Functions:

## Tree reading helper funtion
read_tree <- function(path, format = "auto") {
    if (!file.exists(path)) stop(sprintf("Tree file not found: %s", path))
    if (format == "auto") {
        first <- readLines(path, n = 3L)
        format <- if (any(grepl("^(#NEXUS|#nexus)", first))) "nexus" else "newick"
    }
    if (format == "nexus") {
        tree <- read.nexus(path)
    } else {
        tree <- read.tree(path)
    }
    numeric_tips <- tree$tip.label[grepl("^[0-9]+$", tree$tip.label)]
    if (length(numeric_tips) > 0) {
        message(sprintf("Dropping %d numeric tip labels (bootstrap values misread as tips): %s ...",
                        length(numeric_tips), paste(head(numeric_tips, 3), collapse = ", ")))
        tree <- ape::drop.tip(tree, numeric_tips)
    }
    tree
}

## Input validation helper function
check_inputs <- function(host_tree, symbiont_tree, associations, min_hosts, min_symbionts, span_fraction, permutations) {
    if (!inherits(host_tree, "phylo")) stop("Host tree is not a valid phylo object.")
    if (!inherits(symbiont_tree, "phylo")) stop("Symbiont tree is not a valid phylo object.")
    if (!is.data.frame(associations) || !all(c("Host", "Symbiont") %in% colnames(associations))) {
        stop("Associations file must be a tab-delimited file with columns 'Host' and 'Symbiont'.")
    }
    if (min_hosts < 2) stop("Minimum number of hosts per subtree must be at least 2.")
    if (min_symbionts < 2) stop("Minimum number of symbionts per subtree must be at least 2.")
    if (span_fraction <= 0 || span_fraction > 1) stop("Span fraction must be between 0 and 1.")
    if (permutations < 1) stop("Number of permutations must be at least 1.")
}

## Host-symbiont association matrix construction helper function
host_symbiont_links <- function(df) {
    hosts <- sort(unique(df$Host))
    symbionts <- sort(unique(df$Symbiont))
    mat <- matrix(0L, nrow = length(hosts), ncol = length(symbionts), dimnames = list(hosts, symbionts))
    for (i in seq_len(nrow(df))) mat[df$Host[i], df$Symbiont[i]] <- 1L
    mat
}

## Collapse tree helper function
collapse_monophyletic <- function(tree, Host_to_symbiont_df) {
    tips_to_drop <- c()
    subtree_list <- lapply(seq_len(ape::Nnode(tree)), function(node) 
        castor::get_subtree_at_node(tree, node)$subtree)
    names(subtree_list) <- tree$node.label

    for (i in seq_along(subtree_list)) {
        i_tree <- subtree_list[[i]]
        i_hosts <- Host_to_symbiont_df$Host[
            match(i_tree$tip.label, Host_to_symbiont_df$Symbiont)
        ]
        if (length(unique(i_hosts)) == 1L) {
            dists <- castor::get_all_distances_to_root(i_tree)[seq_len(length(i_tree$tip.label))]
            tips_to_drop <- unique(c(tips_to_drop, i_tree$tip.label[dists != max(dists)]))
        }
    }
    ape::drop.tip(tree, tips_to_drop)
}

## Tree shape statistics helper function
safe_tree_stats <- function(tree) {
    tree <- ape::collapse.singles(tree)
    if (!ape::is.binary(tree)) tree <- ape::multi2di(tree)
    if (!ape::is.rooted(tree)) tree <- ape::root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
    ## castor::tree_imbalance requires >= 3 tips
    if (length(tree$tip.label) < 3L) return(list(colless = NA_real_, sackin = NA_real_))
    tryCatch({
        list(
            colless = castor::tree_imbalance(tree, type = "Colless"),
            sackin  = castor::tree_imbalance(tree, type = "Sackin")
        )
    }, error = function(e) {
        warning("Tree statistics calculation failed: ", e$message)
        list(colless = NA_real_, sackin = NA_real_)
    })
}

## Hommola's method helper function
hommola <- function(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df) {
    # long format symbiont pairwise distances
    dist_df <- tibble::as_tibble(as.matrix(i_symbiont_subtree_dist), rownames = "Symbiont_1") %>%
        tidyr::pivot_longer(-Symbiont_1, names_to = "Symbiont_2", values_to = "Symbiont_Distance")

    # remove self comparisons
    dist_df <- dist_df[dist_df$Symbiont_1 != dist_df$Symbiont_2, ]

    # look up the host for each symbiont in the pair
    dist_df$Host_1 <- i_host_to_symbiont_df$Host[
        match(dist_df$Symbiont_1, i_host_to_symbiont_df$Symbiont)
    ]
    dist_df$Host_2 <- i_host_to_symbiont_df$Host[
        match(dist_df$Symbiont_2, i_host_to_symbiont_df$Symbiont)
    ]

    # long format host pairwise distances
    host_dist <- tibble::as_tibble(as.matrix(i_host_subtree_dist), rownames = "Host_1") %>%
        tidyr::pivot_longer(-Host_1, names_to = "Host_2", values_to = "Host_Distance")
    
    # Join the host distances to the symbiont pairs (Host_1, Host_2)
    dist_df <- suppressMessages(dplyr::left_join(dist_df, host_dist))

    # Drop pairs where host distance could not be looked up (host outside subtree)
    dist_df <- dist_df[!is.na(dist_df$Host_Distance), ]

    # Calculate the correlation between host and symbiont distances
    if (nrow(dist_df) == 0L || length(unique(dist_df$Host_Distance)) == 1L) {
        return(NA_real_)
    }
    cor(dist_df$Host_Distance, dist_df$Symbiont_Distance, method = "pearson")
}

## Hommola Workflow wrapper function
hommola_wf <- function(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df, permutations, seed) {
    hommola_res <- hommola(i_host_subtree_dist, i_symbiont_subtree_dist, i_host_to_symbiont_df)

    if (!is.na(seed)) set.seed(seed)
    permutation_mat <- replicate(permutations, sample(seq_len(nrow(i_host_to_symbiont_df))))

    # list of permuted association data frames
    i_host_to_symbiont_permutations_df <- lapply(seq_len(permutations), function(i) {
        data.frame(
            Host = i_host_to_symbiont_df[, "Host"],
            Symbiont = i_host_to_symbiont_df[permutation_mat[, i], "Symbiont"]
        )
    })

    # calculate hommola for each permuted association
    hommola_res_perms <- vector("list", permutations)
    for (i in seq_len(permutations)) {
        hommola_res_perms[[i]] <- hommola(
            i_host_subtree_dist, 
            i_symbiont_subtree_dist, 
            i_host_to_symbiont_permutations_df[[i]]
            )
    }

    # calculate p-value
    p_value <- (sum(unlist(hommola_res_perms) >= hommola_res) + 1) / (permutations + 1)

    list(Hommola_r = hommola_res,
         Hommola_p = p_value)
}

## PACo Workflow wrapper function

### raw Procrustes residuals helper function
.paco_ss <- function(host_coords, symbiont_coords, assoc_mat) {
   ## build match matrix for symbiont to host coordinates
   links <- which(assoc_mat > 0, arr.ind = TRUE)
   host_pts <- host_coords[links[, 1L], , drop = FALSE]
   sym_pts <- symbiont_coords[links[, 2L], , drop = FALSE]
   ## procrustes residuals
   vegan::procrustes(host_pts, sym_pts)$ss
}

### Vegan PACo wrapper function
paco_wf <- function(host_tree, sym_tree, assoc_mat, permutations, seed) {

    set.seed(seed)
    
    host_d <- cophenetic(host_tree)
    sym_d <- cophenetic(sym_tree)

    ## PCoA with Cailliez correction for negative eigenvalues
    host_pcoa <- vegan::wcmdscale(host_d, eig = TRUE, add = "cailliez")
    sym_pcoa <- vegan::wcmdscale(sym_d, eig = TRUE, add = "cailliez")

    ## keep only the dimensions with positive eigenvalues to avoid noise
    host_coords <- host_pcoa$points[, host_pcoa$eig > 0, drop = FALSE]
    sym_coords <- sym_pcoa$points[, sym_pcoa$eig > 0, drop = FALSE]

    ## ensure row names align with association matrix
    rownames(host_coords) <- rownames(host_d)
    rownames(sym_coords) <- rownames(sym_d)
    host_coords <- host_coords[rownames(assoc_mat), , drop = FALSE]
    sym_coords  <- sym_coords[colnames(assoc_mat), , drop = FALSE]

    ## observed PACo statistic
    obs_ss <- .paco_ss(host_coords, sym_coords, assoc_mat)
    ## permuted PACo statistics
    perm_ss <- vapply(seq_len(permutations), function(i) {
        perm_mat <- assoc_mat
        rownames(perm_mat) <- sample(rownames(assoc_mat))
        perm_mat <- perm_mat[rownames(assoc_mat), , drop = FALSE]
        .paco_ss(host_coords, sym_coords, perm_mat)
    }, numeric(1L))

    ## lower ss = stronger signal, so p = prop. of perms <= obs
    p_value <- (sum(perm_ss <= obs_ss) + 1L) / (permutations + 1L)

    ## return results in same format 
    list(gof = list(ss = obs_ss, p = p_value))
}

# Main codiv function
codiv <- function(Host_tree,
                 Symbiont_tree,
                 Host_to_Symbiont_df,
                 min_hosts = 3L,
                 min_symbionts = 3L,
                 max_symbionts = 500L,
                 span_fraction = 0.1,
                 permutations = 99L,
                 Save_fp,
                 seed = 42L,
                 cores = 4L,
                 methods = c("hommola", "paco", "parafit"),
                 focus_host = character(0),
                 subtree_features = TRUE,
                 continue = TRUE,
                 rarefy_depth = 0L,
                 verbose = TRUE) {
    
    # Input validation
    check_inputs(Host_tree, 
                Symbiont_tree, 
                Host_to_Symbiont_df, 
                min_hosts, 
                min_symbionts, 
                span_fraction, 
                permutations
                )

    # Node labels are required for subtree extraction, so we add them if missing
    if (!is.null(Symbiont_tree$node.label) &&
        anyDuplicated(Symbiont_tree$node.label)) {
            if (verbose) message("Prepending node labels with unique IDs.")
            Symbiont_tree$node.label <- paste0("Node", 
                seq_len(Symbiont_tree$Nnode), "_", Symbiont_tree$node.label)
        }
        if (is.null(Symbiont_tree)) {
            if (verbose) message("Adding unique node labels to symbiont tree.")
            Symbiont_tree$node.label <- paste0("Node_", seq_len(Symbiont_tree$Nnode))
    }

    # Build Subtree list for symbiont tree
    n_nodes <- ape::Nnode(Symbiont_tree)
    subtree_list <- lapply(seq_len(n_nodes), function(node)
        castor::get_subtree_at_node(Symbiont_tree, node)$subtree)
    names(subtree_list) <- Symbiont_tree$node.label

    span_list <- vapply(subtree_list, function(subtree) {
        castor::get_tree_span(subtree, as_edge_count = FALSE)$max_distance
    }, numeric(1))

    sym_tip_counts <- castor::count_tips_per_node(Symbiont_tree)
    host_tip_counts <- vapply(subtree_list, function(subtree)
        length(unique(Host_to_Symbiont_df$Host[
            match(subtree$tip.label, Host_to_Symbiont_df$Symbiont)
        ])), integer(1))
    
    # Filter subtrees based on criteria
    Symbiont_df <- data.frame(
        Node_Label = Symbiont_tree$node.label,
        node = seq_len(n_nodes),
        Subtree_Span = span_list,
        Subtree_Symbiont_Tips = sym_tip_counts,
        Subtree_Host_Tips = host_tip_counts,
        stringsAsFactors = FALSE
    )

    Span_CutOff <- span_fraction * max(Symbiont_df$Subtree_Span)
    n1 <- nrow(Symbiont_df)
    ## filter 1: span
    Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Span <= Span_CutOff, ]
    n2 <- nrow(Symbiont_df)
    if (verbose) message(sprintf("Filtered %d subtrees based on span cutoff (%.2f%% of max span). Remaining: %d",
                                 n1 - n2, span_fraction * 100, n2))
    ## filter 2: min symbiont tips
    Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Symbiont_Tips >= as.integer(min_symbionts), ]
    n3 <- nrow(Symbiont_df)
    if (verbose) message(sprintf("Filtered %d subtrees based on minimum symbiont tips (>= %d). Remaining: %d",
                                 n2 - n3, as.integer(min_symbionts), n3))
    ## filter 3: max symbiont tips
    Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Symbiont_Tips <= as.integer(max_symbionts), ]
    n3b <- nrow(Symbiont_df)
    if (verbose) message(sprintf("Filtered %d subtrees exceeding max symbiont tips (<= %d). Remaining: %d",
                                 n3 - n3b, as.integer(max_symbionts), n3b))
    ## filter 4: min host tips
    Symbiont_df <- Symbiont_df[Symbiont_df$Subtree_Host_Tips >= as.integer(min_hosts), ]
    n4 <- nrow(Symbiont_df)
    if (verbose) message(sprintf("Filtered %d subtrees based on minimum host tips (>= %d). Remaining: %d",
                                 n3b - n4, as.integer(min_hosts), n4))

    Symbiont_df <- Symbiont_df[order(Symbiont_df$Subtree_Symbiont_Tips, Symbiont_df$Subtree_Host_Tips), ]

    # column schema for results
    basic_vars <- c("Node_ID", "Symbiont_Tree", "N_Symbionts", "Symbiont_Colless", "Symbiont_Sackin", "Host_Tree", "N_Hosts", "Host_Colless", "Host_Sackin")
    Hommola_vars <- if ("hommola" %in% methods) 
        c("Hommola_r", "Hommola_pvalue", "Rarefied_Hommola_r", "Rarefied_Hommola_pvalue") else c()
    PACo_vars <- if ("paco" %in% methods) 
        c("Rarefied_PACo_ss", "Rarefied_PACo_pvalue") else c()
    Parafit_vars <- if ("parafit" %in% methods)
        c("Rarefied_ParaFitGlobal", "Rarefied_ParaFit_pvalue") else c()
    subtree_vars <- if (subtree_features) 
        c("Rarefied_TreeDistance", "Rarefied_SharedPhylogenticInfo", 
            "Rarefied_DifferentPhylogeneticInfo", "Rarefied_NyeSimilarity",
      "Rarefied_JaccardRobinsonFoulds","Rarefied_MatchingSplitDistance",
      "Rarefied_MatchingSplitInfoDistance","Rarefied_MutualClusteringInfo") else c()
    
    focus_hosts_vars <- if (length(focus_host)) paste0(focus_host, "_PRESENT") else c()

    cols_vars <- c(basic_vars, Hommola_vars, PACo_vars, Parafit_vars, subtree_vars, focus_hosts_vars)

    # resume: skip completed subtrees
    nodes_to_scan <- Symbiont_df$Node_Label
    if (continue && file.exists(Save_fp)) {
        prev <- suppressWarnings(
            read.table(Save_fp, sep = "\t", header = TRUE, fill = TRUE))
        if ("Hommola_pvalue" %in% colnames(prev)) {
            done <- prev$Node_ID[!is.na(prev$Hommola_pvalue)]
            nodes_to_scan <- nodes_to_scan[!nodes_to_scan %in% done]
            if (verbose) message(sprintf("Resuming from previous run. Skipping %d completed subtrees. Remaining: %d", length(done), length(nodes_to_scan)))

        }
    }


    ## write skeleton only if no existing results — preserves completed data for resume
    if (!file.exists(Save_fp)) {
        skel <- data.frame(matrix(nrow = nrow(Symbiont_df), ncol = length(cols_vars)), stringsAsFactors = FALSE)
        colnames(skel) <- cols_vars
        skel$Node_ID <- Symbiont_df$Node_Label
        write.table(skel, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    ## Per Node temp directory for intermediate results
    tmp_dir <- file.path(dirname(normalizePath(Save_fp)),
                    paste0(".codiv_tmp_", basename(Save_fp)))

    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

    # worker function to process each subtree
    ## fast_dist defined inside codiv() so mclapply forked workers can see it
    fast_dist <- function(tree) {
        n <- length(tree$tip.label)
        ## only_clades = 1:Ntips restricts to tip-tip distances only (faster)
        ## returns n x n matrix; castor indexes tips 1..Ntips
        mat <- castor::get_all_pairwise_distances(tree,
                    only_clades = seq_len(n),
                    as_edge_counts = FALSE)
        rownames(mat) <- colnames(mat) <- tree$tip.label
        as.dist(mat)
    }

    process_node <- function(i) {

        tmp_fp <- file.path(tmp_dir, paste0(make.names(i), ".rds"))
        if (file.exists(tmp_fp)) return(invisible(NULL))

        row <- setNames(vector("list", length(cols_vars)), cols_vars)
        row$Node_ID <- i

        tryCatch({

            i_sym_sub <- subtree_list[[i]]
            i_h2s_df <- Host_to_Symbiont_df[
                Host_to_Symbiont_df$Symbiont %in% i_sym_sub$tip.label,
            ]
            i_host_sub <- ape::keep.tip(Host_tree, unique(i_h2s_df$Host))
            i_h2s_mat <- host_symbiont_links(i_h2s_df)

            i_sym_dist <- fast_dist(i_sym_sub)
            i_host_dist <- fast_dist(i_host_sub)

            ## ── Rarefaction: equalise symbiont sampling across hosts ──────────
            ## For each host, draw rarefy_n symbionts without replacement.
            ## rarefy_n = min per-host count in this subtree when rarefy_depth==0.
            per_host_n <- table(i_h2s_df$Host)
            rarefy_n   <- if (rarefy_depth > 0L) as.integer(rarefy_depth) else as.integer(min(per_host_n))
            rarefy_n   <- max(1L, rarefy_n)  ## safety floor

            node_rng <- seed + which(nodes_to_scan == i)
            set.seed(node_rng)
            i_rar_syms <- unlist(lapply(names(per_host_n), function(h) {
                syms <- i_h2s_df$Symbiont[i_h2s_df$Host == h]
                if (length(syms) <= rarefy_n) syms else sample(syms, rarefy_n)
            }))
            ## B: pre-filter to labels that actually exist in the symbiont tree
            ## (catches name mismatches between association file and tree tip labels)
            i_rar_syms <- intersect(i_rar_syms, i_sym_sub$tip.label)

            i_rar_sym <- ape::keep.tip(i_sym_sub, i_rar_syms)
            i_rar_h2s <- i_h2s_df[i_h2s_df$Symbiont %in% i_rar_syms, ]
            i_rar_mat <- host_symbiont_links(i_rar_h2s)

            ## A+B guard: check keep.tip output and host count before proceeding
            if (is.null(i_rar_sym) ||
                length(i_rar_sym$tip.label) < 2L ||
                length(unique(i_rar_h2s$Host)) < 2L) {
                row$Symbiont_Tree <- ape::write.tree(i_sym_sub)
                row$N_Symbionts   <- length(i_sym_sub$tip.label)
                row$Host_Tree     <- ape::write.tree(i_host_sub)
                row$N_Hosts       <- length(i_host_sub$tip.label)
                row <- lapply(row, function(x) if (is.null(x) || length(x) == 0L) NA else x[1L])
                saveRDS(as.data.frame(row, stringsAsFactors = FALSE), tmp_fp)
                return(invisible(NULL))
            }

            i_rar_sym_d <- fast_dist(i_rar_sym)

            ## basic subtree features
            row$Symbiont_Tree <- ape::write.tree(i_sym_sub)
            row$N_Symbionts <- length(i_sym_sub$tip.label)
            row$Host_Tree <- ape::write.tree(i_host_sub)
            row$N_Hosts <- length(i_host_sub$tip.label)

            ss <- safe_tree_stats(i_sym_sub); hs <- safe_tree_stats(i_host_sub)
            row$Symbiont_Colless <- ss$colless; row$Symbiont_Sackin <- ss$sackin
            row$Host_Colless <- hs$colless; row$Host_Sackin <- hs$sackin

            ## Hommola's method
            if ("hommola" %in% methods) {
                tryCatch({
                    node_seed <- seed + which(nodes_to_scan == i)
                    hr  <- hommola_wf(i_host_dist, i_sym_dist,  i_h2s_df,  permutations, node_seed)
                    rhr <- hommola_wf(i_host_dist, i_rar_sym_d, i_rar_h2s, permutations, node_seed)
                    row$Hommola_r              <- hr$Hommola_r
                    row$Hommola_pvalue         <- hr$Hommola_p
                    row$Rarefied_Hommola_r     <- rhr$Hommola_r
                    row$Rarefied_Hommola_pvalue <- rhr$Hommola_p
                }, error = function(e) warning("Hommola failed for node ", i, ": ", e$message))
            }

            ## PACo method
            if ("paco" %in% methods) {
                tryCatch({
                    pr <- paco_wf(i_host_sub, i_rar_sym, i_rar_mat, permutations, seed + which(nodes_to_scan == i))
                    row$Rarefied_PACo_ss      <- pr$gof$ss
                    row$Rarefied_PACo_pvalue  <- pr$gof$p
                }, error = function(e) warning("PACo failed for node ", i, ": ", e$message))
            }

            ## Parafit method
            if ("parafit" %in% methods) {
                tryCatch({
                    ## guard: parafit requires >= 2 hosts and >= 2 symbionts in the matrix
                    if (nrow(i_rar_mat) >= 2L && ncol(i_rar_mat) >= 2L) {
                        pfr <- ape::parafit(i_host_dist, i_rar_sym_d,
                                            i_rar_mat, nperm = permutations, test.links = FALSE,
                                            seed = seed + which(nodes_to_scan == i),
                                            correction = "cailliez", silent = TRUE)
                        row$Rarefied_ParaFitGlobal  <- pfr$ParaFitGlobal
                        row$Rarefied_ParaFit_pvalue <- pfr$p.global
                    }
                }, error = function(e) warning("ParaFit failed for node ", i, ": ", e$message))
            }

            ## Focus host presence/absence
            if (length(focus_host)) {
                present <- vapply(focus_host,
                                    function(x) x %in% i_h2s_df$Host, logical(1))
                for (k in seq_along(focus_host)) {
                    row[[focus_hosts_vars[k]]] <- present[[k]]
                }
            }

            ## Subtree topology features (TreeDist via host-mirrored symbiont tree)
            ##
            ## Since symbiont tips are MAG IDs and host tips are species names they
            ## share no labels, so direct tip-matching is impossible.  Instead we
            ## build a "host-mirrored" symbiont tree: keep one representative symbiont
            ## per host (the one with the greatest mean patristic distance to all
            ## other symbionts in that host, i.e. the most phylogenetically distinct),
            ## rename its tip to the host name, and prune the host tree to the same
            ## host set.  Both trees then share tip labels and TreeDist metrics are
            ## meaningful: they capture whether the symbiont phylogeny mirrors the
            ## host phylogeny after collapsing within-host diversity to one point.
            if (subtree_features) {
                tryCatch({
                    rar_mat_full <- as.matrix(i_rar_sym_d)

                    ## pick one representative symbiont per host
                    rep_syms <- vapply(names(per_host_n), function(h) {
                        syms_h <- i_rar_h2s$Symbiont[i_rar_h2s$Host == h]
                        syms_h <- intersect(syms_h, rownames(rar_mat_full))
                        if (length(syms_h) == 1L) return(syms_h)
                        ## most phylogenetically distinct = highest mean dist to all others
                        mean_dists <- rowMeans(rar_mat_full[syms_h, , drop = FALSE])
                        syms_h[which.max(mean_dists)]
                    }, character(1))

                    ## build host-mirrored symbiont tree: keep representatives, rename tips
                    mirror_sym  <- ape::keep.tip(i_rar_sym, rep_syms)
                    mirror_sym$tip.label <- names(per_host_n)[match(mirror_sym$tip.label, rep_syms)]

                    ## prune host tree to the same host set
                    mirror_host <- ape::keep.tip(i_host_sub, mirror_sym$tip.label)

                    if (length(mirror_sym$tip.label) >= 2L && length(mirror_host$tip.label) >= 2L) {
                        ## resolve any polytomies before TreeDist (requires binary trees)
                        mirror_sym  <- ape::multi2di(mirror_sym)
                        mirror_host <- ape::multi2di(mirror_host)

                        row$Rarefied_TreeDistance              <- TreeDist::TreeDistance(mirror_sym, mirror_host)
                        row$Rarefied_SharedPhylogenticInfo     <- TreeDist::SharedPhylogeneticInfo(mirror_sym, mirror_host)
                        row$Rarefied_DifferentPhylogeneticInfo <- TreeDist::DifferentPhylogeneticInfo(mirror_sym, mirror_host)
                        row$Rarefied_NyeSimilarity             <- TreeDist::NyeSimilarity(mirror_sym, mirror_host)
                        row$Rarefied_JaccardRobinsonFoulds     <- TreeDist::JaccardRobinsonFoulds(mirror_sym, mirror_host)
                        row$Rarefied_MatchingSplitDistance     <- TreeDist::MatchingSplitDistance(mirror_sym, mirror_host)
                        row$Rarefied_MatchingSplitInfoDistance <- TreeDist::MatchingSplitInfoDistance(mirror_sym, mirror_host)
                        row$Rarefied_MutualClusteringInfo      <- TreeDist::MutualClusteringInfo(mirror_sym, mirror_host)
                    }
                }, error = function(e) {
                    warning("TreeDist calculation failed for node ", i, ": ", e$message)
                })
            }
        }, error = function(e) {
            row$Node_ID <<- i
            row$Symbiont_Tree <<- paste0("ERROR: ", e$message)
        })

        ## coerce every list element to scalar NA if NULL/length-0
        ## so as.data.frame never sees differing row counts
        row <- lapply(row, function(x) {
            if (is.null(x) || length(x) == 0L) NA else x[1L]
        })
        saveRDS(as.data.frame(row, stringsAsFactors = FALSE), tmp_fp)
        invisible(NULL)
    }

    ## dispatch processing across cores
    if (verbose) message(sprintf("Processing %d subtrees across %d cores...", length(nodes_to_scan), cores))

    ## DEBUG: run first node serially in main process so errors print fully.
    ## Remove or set debug_mode <- FALSE once working.
    debug_mode <- Sys.getenv("CODIV_DEBUG", unset = "0") == "1"
    if (debug_mode) {
        message("[DEBUG] Running first node serially to expose errors...")
        ## strip tryCatch for raw error
        i <- nodes_to_scan[1]
        i_sym_sub   <- subtree_list[[i]]
        i_h2s_df    <- Host_to_Symbiont_df[Host_to_Symbiont_df$Symbiont %in% i_sym_sub$tip.label, ]
        i_host_sub  <- ape::keep.tip(Host_tree, unique(i_h2s_df$Host))
        i_h2s_mat   <- host_symbiont_links(i_h2s_df)
        i_sym_dist  <- fast_dist(i_sym_sub)
        i_host_dist <- fast_dist(i_host_sub)
        per_host_n  <- table(i_h2s_df$Host)
        rarefy_n    <- if (rarefy_depth > 0L) as.integer(rarefy_depth) else as.integer(min(per_host_n))
        set.seed(seed + 1L)
        i_rar_syms  <- unlist(lapply(names(per_host_n), function(h) {
            syms <- i_h2s_df$Symbiont[i_h2s_df$Host == h]; if (length(syms) <= rarefy_n) syms else sample(syms, rarefy_n)
        }))
        i_rar_sym   <- ape::keep.tip(i_sym_sub, i_rar_syms)
        i_rar_h2s   <- i_h2s_df[i_h2s_df$Symbiont %in% i_rar_syms, ]
        i_rar_mat   <- host_symbiont_links(i_rar_h2s)
        i_rar_sym_d <- fast_dist(i_rar_sym)
        node_seed   <- seed + 1L
        message("[DEBUG] Setup OK. Testing hommola...")
        hr  <- hommola_wf(i_host_dist, i_sym_dist,  i_h2s_df,  10L, node_seed)
        rhr <- hommola_wf(i_host_dist, i_rar_sym_d, i_rar_h2s, 10L, node_seed)
        message(sprintf("[DEBUG] Hommola OK: r=%.3f p=%.3f  Rarefied: r=%.3f p=%.3f", hr$Hommola_r, hr$Hommola_p, rhr$Hommola_r, rhr$Hommola_p))
        message("[DEBUG] Testing PACo...")
        pr  <- paco_wf(i_host_sub, i_rar_sym, i_rar_mat, 10L, node_seed)
        message(sprintf("[DEBUG] PACo OK: ss=%.4f p=%.3f", pr$gof$ss, pr$gof$p))
        message("[DEBUG] Testing ParaFit...")
        pfr <- ape::parafit(i_host_dist, i_rar_sym_d, i_rar_mat,
                            nperm = 10L, test.links = FALSE,
                            correction = "cailliez", silent = TRUE)
        message(sprintf("[DEBUG] ParaFit OK: global=%.4f p=%.3f", pfr$ParaFitGlobal, pfr$p.global))
        message("[DEBUG] All methods passed for first node. Proceeding with full parallel run.")
    }

    parallel::mclapply(
        nodes_to_scan, 
        process_node, 
        mc.cores = cores,
        mc.preschedule = FALSE
    )

    ## aggregate results
    if (verbose) message("Aggregating results...")


    ## if all nodes were skipped (resume), read existing TSV directly
    if (length(nodes_to_scan) == 0L && file.exists(Save_fp)) {
        if (verbose) message("All nodes already complete — reading from: ", Save_fp)
        Results_df <- read.table(Save_fp, sep = "\t", header = TRUE,
                                 stringsAsFactors = FALSE, fill = TRUE)
        attr(Results_df, "tmp_dir") <- tmp_dir
        return(invisible(Results_df))
    }
    Results_df <- data.frame(
        matrix(nrow = nrow(Symbiont_df), ncol = length(cols_vars)),
        stringsAsFactors = FALSE)
    colnames(Results_df) <- cols_vars
    Results_df$Node_ID <- Symbiont_df$Node_Label

    for (fp in list.files(tmp_dir, pattern = "\\.rds$", full.names = TRUE)) {
        row <- tryCatch(readRDS(fp), error = function(e) NULL)
        if (is.null(row)) next
        idx <- which(Results_df$Node_ID == row$Node_ID[1])
        if (length(idx) == 1L) {
            Results_df[idx, names(row)] <- row
        }
    }

    write.table(Results_df, Save_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    ## tmp_dir kept intentionally — cleaned up by main() after FDR completes
    ## so a re-run can skip node processing and go straight to FDR

    if (verbose) message("Done. Results saved to: ", Save_fp)
    attr(Results_df, "tmp_dir") <- tmp_dir  ## pass path to caller for cleanup
    invisible(Results_df)
}

# =============================================================================
# Empirical Null Model via Full Traversal on Permuted Host Trees
#
# Strategy (mirrors the notebook's approach)
# -------------------------------------------
# Repeat the ENTIRE codiv() traversal N times, each time with host tip labels
# randomly permuted on the host tree. This breaks all real co-diversification
# signal while preserving the structural properties of both trees and the
# marginal distribution of symbionts across hosts.
#
# For each permuted traversal, count the number of nodes that pass a combined
# significance threshold (p <= null_alpha AND stat >= null_r_threshold for the
# chosen method). The resulting N-length count vector is the empirical null
# distribution of "how many significant nodes you expect by chance".
#
# The observed count from the real data is then compared against this null to
# produce an empirical p-value:
#   p_empirical = (# null traversals with count >= observed_count + 1) / (N + 1)
#
# This is output as a summary TSV (<output>_null_model_summary.tsv) and a
# per-permutation counts file (<output>_null_model_counts.tsv).
#
# Parameters
# ----------
# results_df          : data frame returned by codiv() on the real data
# Host_tree           : host phylo object
# Symbiont_tree       : symbiont phylo object
# Host_to_Symbiont_df : full association data frame
# n_traversals        : number of permuted traversals (recommend >= 100)
# methods             : which methods to count hits for (auto-detected if NULL)
# alpha               : p-value threshold for counting a node as a hit
# r_threshold         : minimum stat value threshold (e.g. Hommola r >= 0.75);
#                       set to 0 to count any significant node regardless of effect size
# output_prefix       : file path prefix for summary and counts TSVs
# min_hosts, min_symbionts, max_symbionts, span_fraction : same filters as codiv()
# permutations        : inner permutations per node (keep low, e.g. 99, for speed)
# seed                : RNG seed for first traversal; each traversal uses seed+i
# cores               : parallel cores passed to codiv()
# verbose             : print progress
#
# Returns
# -------
# A list with:
#   $null_counts   : data frame of per-traversal hit counts per method
#   $observed      : named vector of observed hit counts per method
#   $empirical_p   : named vector of empirical p-values per method
#   $summary       : tidy summary data frame
# =============================================================================

null_model_traversal <- function(results_df,
                                  Host_tree,
                                  Symbiont_tree,
                                  Host_to_Symbiont_df,
                                  n_traversals      = 100L,
                                  methods           = NULL,
                                  alpha             = 0.05,
                                  r_threshold       = 0,
                                  output_prefix     = NULL,
                                  min_hosts         = 3L,
                                  min_symbionts     = 3L,
                                  max_symbionts     = 500L,
                                  span_fraction     = 0.1,
                                  permutations      = 99L,
                                  seed              = 42L,
                                  cores             = 4L,
                                  verbose           = TRUE) {

    ## ── method schema: which column holds the stat and which p-value column ──
    method_schema <- list(
        hommola = list(stat_col = "Rarefied_Hommola_r",      p_col = "Rarefied_Hommola_pvalue",  higher_is_good = TRUE),
        paco    = list(stat_col = "Rarefied_PACo_ss",         p_col = "Rarefied_PACo_pvalue",      higher_is_good = FALSE),
        parafit = list(stat_col = "Rarefied_ParaFitGlobal",   p_col = "Rarefied_ParaFit_pvalue",   higher_is_good = TRUE)
    )

    ## auto-detect methods present in results_df
    if (is.null(methods)) {
        methods <- names(method_schema)[vapply(names(method_schema), function(m) {
            cols <- c(method_schema[[m]]$stat_col, method_schema[[m]]$p_col)
            all(cols %in% colnames(results_df)) && any(!is.na(results_df[[method_schema[[m]]$p_col]]))
        }, logical(1))]
    }

    if (length(methods) == 0L)
        stop("No recognised method columns found in results_df. Run codiv() first.")

    if (verbose)
        message(sprintf("[NullModel] Methods detected: %s", paste(methods, collapse = ", ")))

    ## ── helper: count significant hits for one results data frame ─────────────
    count_hits <- function(df, m) {
        schema    <- method_schema[[m]]
        p_vec     <- suppressWarnings(as.numeric(df[[schema$p_col]]))
        stat_vec  <- suppressWarnings(as.numeric(df[[schema$stat_col]]))
        sig_p     <- !is.na(p_vec)    & p_vec    <= alpha
        if (r_threshold != 0) {
            if (schema$higher_is_good)
                sig_stat <- !is.na(stat_vec) & stat_vec >= r_threshold
            else
                sig_stat <- !is.na(stat_vec) & stat_vec <= r_threshold
            return(sum(sig_p & sig_stat, na.rm = TRUE))
        }
        sum(sig_p, na.rm = TRUE)
    }

    ## ── observed hit counts from real data ────────────────────────────────────
    observed <- setNames(
        vapply(methods, function(m) count_hits(results_df, m), integer(1)),
        methods
    )

    if (verbose) {
        message("[NullModel] Observed hit counts (raw inner p-values from real data):")
        for (m in methods)
            message(sprintf("  [%s] %d significant nodes (alpha=%.2f)", m, observed[[m]], alpha))
        if (all(observed == 0L))
            message("[NullModel] WARNING: observed hit count is 0 for all methods. ",
                    "This can happen if the inner permutation p-values are all non-significant. ",
                    "Consider lowering --null_alpha or checking that codiv() ran correctly.")
        message(sprintf("[NullModel] Running %d permuted traversals...", n_traversals))
    }

    ## ── build a temp directory for permuted runs ──────────────────────────────
    tmp_null_dir <- if (!is.null(output_prefix))
        paste0(output_prefix, "_null_tmp")
    else
        tempfile("codiv_null_")
    dir.create(tmp_null_dir, showWarnings = FALSE, recursive = TRUE)

    host_tips <- Host_tree$tip.label

    ## ── run one permuted traversal ─────────────────────────────────────────────
    run_one_perm <- function(perm_idx) {

        perm_tsv <- file.path(tmp_null_dir, sprintf("perm_%04d.tsv", perm_idx))

        ## skip if already done (allows resuming interrupted null runs)
        if (file.exists(perm_tsv)) {
            if (verbose) message(sprintf("[NullModel] Permutation %d already done, skipping.", perm_idx))
            df <- tryCatch(
                read.table(perm_tsv, sep = "\t", header = TRUE,
                           stringsAsFactors = FALSE, fill = TRUE),
                error = function(e) NULL)
            if (!is.null(df)) return(vapply(methods, function(m) count_hits(df, m), integer(1)))
        }

        set.seed(seed + perm_idx)
        permuted_tips              <- sample(host_tips)
        permuted_host_tree         <- Host_tree
        permuted_host_tree$tip.label <- permuted_tips

        ## remap association table to permuted labels
        tip_map <- setNames(permuted_tips, host_tips)
        perm_assoc <- Host_to_Symbiont_df
        perm_assoc$Host <- tip_map[perm_assoc$Host]

        ## run the full traversal (no inner FDR, no subtree features — keep it fast)
        perm_results <- tryCatch(
            codiv(
                Host_tree           = permuted_host_tree,
                Symbiont_tree       = Symbiont_tree,
                Host_to_Symbiont_df = perm_assoc,
                min_hosts           = min_hosts,
                min_symbionts       = min_symbionts,
                max_symbionts       = max_symbionts,
                span_fraction       = span_fraction,
                permutations        = permutations,
                Save_fp             = perm_tsv,
                seed                = seed + perm_idx,
                cores               = cores,
                methods             = methods,
                focus_host          = character(0),
                subtree_features    = FALSE,   ## skip tree-distance features for speed
                continue            = TRUE,
                rarefy_depth        = 0L,      ## auto within each permuted subtree
                verbose             = FALSE
            ),
            error = function(e) {
                warning(sprintf("[NullModel] Permutation %d failed: %s", perm_idx, e$message))
                NULL
            }
        )

        if (is.null(perm_results))
            return(setNames(rep(NA_integer_, length(methods)), methods))

        vapply(methods, function(m) count_hits(perm_results, m), integer(1))
    }

    ## run traversals sequentially (parallel at the node level inside codiv())
    null_counts_list <- lapply(seq_len(n_traversals), function(i) {
        if (verbose) message(sprintf("[NullModel] Permutation %d / %d ...", i, n_traversals))
        run_one_perm(i)
    })

    ## ── reshape into data frame ───────────────────────────────────────────────
    ## Guard: if all traversals failed, null_counts_list is all NA vectors
    n_valid <- sum(vapply(null_counts_list, function(x) !any(is.na(x)), logical(1)))
    if (n_valid == 0L)
        stop("[NullModel] All permuted traversals failed or returned NA. ",
             "Check that codiv() runs successfully on your data with --null_traversals 1 first.")
    if (n_valid < n_traversals && verbose)
        warning(sprintf("[NullModel] %d / %d traversals failed and will appear as NA in counts.",
                        n_traversals - n_valid, n_traversals))

    null_counts_df <- as.data.frame(do.call(rbind, null_counts_list),
                                    stringsAsFactors = FALSE)
    null_counts_df <- cbind(data.frame(permutation = seq_len(n_traversals)), null_counts_df)
    colnames(null_counts_df) <- c("permutation", methods)

    ## ── empirical p-values ────────────────────────────────────────────────────
    empirical_p <- setNames(
        vapply(methods, function(m) {
            null_vec <- null_counts_df[[m]]
            obs      <- observed[[m]]
            (sum(null_vec >= obs, na.rm = TRUE) + 1L) / (sum(!is.na(null_vec)) + 1L)
        }, numeric(1)),
        methods
    )

    ## ── summary data frame ────────────────────────────────────────────────────
    summary_df <- data.frame(
        method          = methods,
        observed_hits   = observed[methods],
        null_mean       = vapply(methods, function(m) mean(null_counts_df[[m]], na.rm = TRUE), numeric(1)),
        null_sd         = vapply(methods, function(m) sd(null_counts_df[[m]],   na.rm = TRUE), numeric(1)),
        null_median     = vapply(methods, function(m) median(null_counts_df[[m]], na.rm = TRUE), numeric(1)),
        null_95pct      = vapply(methods, function(m) quantile(null_counts_df[[m]], 0.95, na.rm = TRUE), numeric(1)),
        empirical_p     = empirical_p[methods],
        alpha_threshold = alpha,
        r_threshold     = r_threshold,
        n_traversals    = n_traversals,
        n_valid_traversals = n_valid,
        stringsAsFactors = FALSE,
        row.names        = NULL
    )

    if (verbose) {
        message("\n[NullModel] Summary:")
        for (i in seq_len(nrow(summary_df))) {
            r <- summary_df[i, ]
            message(sprintf(
                "  [%s] observed=%d  null mean=%.1f (sd=%.1f)  95th pct=%.0f  empirical p=%.4f",
                r$method, r$observed_hits, r$null_mean, r$null_sd, r$null_95pct, r$empirical_p
            ))
        }
    }

    ## ── write output files ────────────────────────────────────────────────────
    ## Derive prefix from tmp_null_dir as fallback so files are ALWAYS written,
    ## even if output_prefix was NULL (e.g. called interactively without a path).
    effective_prefix <- if (!is.null(output_prefix)) {
        output_prefix
    } else {
        ## place files next to the temp dir using its path stem
        sub("_null_tmp$", "", tmp_null_dir)
    }

    ## ensure the output directory exists
    out_dir <- dirname(effective_prefix)
    if (!dir.exists(out_dir))
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    counts_fp  <- paste0(effective_prefix, "_null_model_counts.tsv")
    summary_fp <- paste0(effective_prefix, "_null_model_summary.tsv")

    tryCatch({
        write.table(null_counts_df, counts_fp,  sep = "\t", quote = FALSE, row.names = FALSE)
        if (verbose) message(sprintf("[NullModel] Per-permutation counts written to: %s", counts_fp))
    }, error = function(e) {
        warning(sprintf("[NullModel] Failed to write counts file: %s\n  Error: %s", counts_fp, e$message))
    })

    tryCatch({
        write.table(summary_df, summary_fp, sep = "\t", quote = FALSE, row.names = FALSE)
        if (verbose) message(sprintf("[NullModel] Summary written to: %s", summary_fp))
    }, error = function(e) {
        warning(sprintf("[NullModel] Failed to write summary file: %s\n  Error: %s", summary_fp, e$message))
    })

    ## clean up temp directory
    unlink(tmp_null_dir, recursive = TRUE)

    invisible(list(
        null_counts  = null_counts_df,
        observed     = observed,
        empirical_p  = empirical_p,
        summary      = summary_df
    ))
}

# =============================================================================
# Node-wise Occupancy Scan: Parsimony-based Phylogenetic Signal Test
#
# Motivation
# ----------
# The congruence-based tests (Hommola, PACo, ParaFit) require symbionts to
# appear in multiple host species to form informative pairwise comparisons.
# When host-specificity is high — the majority of symbiont lineages are found
# in only one host — those tests have little power because the host distance
# matrix is nearly zero everywhere.
#
# This scan detects a complementary signal: whether symbiont lineages that ARE
# host-specific tend to cluster together on the symbiont tree according to
# their host of origin. It uses ALL symbiont tips, including host-specific ones,
# and is therefore powered by the dominant signal in highly specific datasets.
#
# Statistical approach
# --------------------
# For each internal node v of the symbiont tree:
#   1. Extract the subtree rooted at v.
#   2. Label each tip with its host species (from the association table).
#      Tips found in multiple hosts receive all their host labels (multi-host).
#   3. Compute the PARSIMONY SCORE: the minimum number of host-state transitions
#      required to explain the observed tip labels on the subtree topology,
#      using Fitch parsimony for unordered discrete states.
#      A low score = host labels are phylogenetically clustered = signal.
#   4. Permute tip labels randomly (preserving per-tip host multiplicity) and
#      recompute parsimony for each permutation to build a null distribution.
#   5. p-value = (# null scores <= observed + 1) / (permutations + 1).
#
# Multi-host tips
# ---------------
# A symbiont found in k hosts contributes k labels. In Fitch parsimony this
# is handled naturally: a multi-host tip's state set is the union of all its
# host labels, which never requires a transition to reach any of them.
# This conservatively reduces the parsimony score for multi-host tips,
# which is appropriate — we don't want to penalise shared lineages.
#
# Effect size: Normalised Parsimony Score (NPS)
# ---------------------------------------------
# Raw parsimony depends on subtree size and number of hosts. We normalise:
#   NPS = 1 - (observed_parsimony / max_possible_parsimony)
# where max_possible = (n_tips - 1) assuming maximum scrambling.
# NPS = 1 means perfect clustering; NPS = 0 means no clustering beyond random.
# NPS is comparable across subtrees of different sizes.
#
# Host phylogenetic weighting (optional, host_tree supplied)
# ----------------------------------------------------------
# Standard Fitch parsimony treats all host transitions as equal cost.
# When a host tree is supplied, transition costs are set proportional to
# host phylogenetic distance, making transitions between distantly related
# hosts more costly. This implements a phylogenetically informed version
# of the test and makes the effect size more biologically interpretable.
# When host_tree = NULL, standard unweighted parsimony is used.
#
# Output columns
# --------------
# Node_ID, N_Tips, N_Hosts_Represented, Observed_Parsimony, Max_Parsimony,
# NPS (normalised parsimony score), Parsimony_p, Host_distribution
#
# Parameters
# ----------
# Symbiont_tree       : symbiont phylo object (full tree, not filtered)
# Host_to_Symbiont_df : association data frame (Host, Symbiont columns)
# Host_tree           : optional host phylo for weighted transitions
# permutations        : number of permutations per node (default 999)
# min_tips            : minimum symbiont tips in subtree to test (default 5)
# span_fraction       : max subtree span as fraction of total (default 1.0,
#                       i.e. test all nodes including root)
# seed                : RNG seed
# cores               : parallel cores
# output_fp           : path to write TSV results
# verbose             : print progress
# =============================================================================

occupancy_scan <- function(Symbiont_tree,
                           Host_to_Symbiont_df,
                           Host_tree         = NULL,
                           permutations      = 999L,
                           min_tips          = 5L,
                           span_fraction     = 1.0,
                           seed              = 42L,
                           cores             = 4L,
                           output_fp         = NULL,
                           verbose           = TRUE) {

    ## ── Fitch parsimony (unweighted) ──────────────────────────────────────────
    ## Implements the two-pass Fitch algorithm for unordered discrete states.
    ## state_sets: named list mapping tip names -> character vector of states
    ## Returns integer parsimony score.
    fitch_parsimony <- function(tree, state_sets) {

        n_tips  <- length(tree$tip.label)
        n_nodes <- tree$Nnode
        total   <- n_tips + n_nodes

        ## initialise state sets for all nodes (tips + internals)
        node_states <- vector("list", total)
        for (i in seq_len(n_tips)) {
            nm <- tree$tip.label[i]
            node_states[[i]] <- if (!is.null(state_sets[[nm]])) state_sets[[nm]] else character(0)
        }

        score  <- 0L
        edge   <- tree$edge          ## n_edges x 2 matrix: [parent, child]
        parent <- edge[, 1L]
        child  <- edge[, 2L]

        ## postorder traversal: process children before parents
        node_order <- rev(unique(parent))  ## internal nodes in postorder approx

        for (p in node_order) {
            children_idx <- child[parent == p]
            ## intersect states of all children
            combined <- node_states[[children_idx[1L]]]
            score_add <- 0L
            for (ci in children_idx[-1L]) {
                inter <- intersect(combined, node_states[[ci]])
                if (length(inter) == 0L) {
                    ## union: no shared state -> parsimony cost +1
                    combined  <- union(combined, node_states[[ci]])
                    score_add <- score_add + 1L
                } else {
                    combined <- inter
                }
            }
            score             <- score + score_add
            node_states[[p]]  <- combined
        }
        score
    }

    ## ── Weighted Fitch parsimony (host-distance weighted) ─────────────────────
    ## Uses host phylogenetic distances as transition costs.
    ## Falls back to unweighted if only 1 host represented.
    fitch_weighted <- function(tree, state_sets, host_dist_mat) {

        n_tips  <- length(tree$tip.label)
        n_nodes <- tree$Nnode
        total   <- n_tips + n_nodes

        ## represent each node's state as a named numeric cost vector
        ## cost[h] = minimum cost to reach host h at this node
        all_hosts <- rownames(host_dist_mat)

        tip_costs <- vector("list", total)
        for (i in seq_len(n_tips)) {
            nm     <- tree$tip.label[i]
            states <- if (!is.null(state_sets[[nm]])) state_sets[[nm]] else all_hosts
            costs  <- setNames(rep(Inf, length(all_hosts)), all_hosts)
            costs[states] <- 0
            tip_costs[[i]] <- costs
        }

        score  <- 0
        edge   <- tree$edge
        parent <- edge[, 1L]
        child  <- edge[, 2L]
        node_order <- rev(unique(parent))

        for (p in node_order) {
            children_idx <- child[parent == p]
            ## for each child, compute min-cost extension to every host state
            child_extended <- lapply(children_idx, function(ci) {
                c_costs <- tip_costs[[ci]]
                ## for each host h, cost = min over all states s of (c_costs[s] + dist[s,h])
                vapply(all_hosts, function(h) {
                    min(c_costs + host_dist_mat[, h])
                }, numeric(1))
            })
            ## parent cost = sum of child extensions (independent branches)
            parent_costs <- Reduce("+", child_extended)
            score        <- score + min(parent_costs)
            ## normalise so root cost is 0 for the optimal state
            tip_costs[[p]] <- parent_costs - min(parent_costs)
        }
        score
    }

    ## ── build host distance matrix if host tree supplied ──────────────────────
    host_dist_mat <- NULL
    if (!is.null(Host_tree)) {
        hd  <- cophenetic(Host_tree)
        ## scale to [0,1] so costs are comparable across datasets
        mx  <- max(hd)
        host_dist_mat <- if (mx > 0) hd / mx else hd
    }

    ## ── pre-build tip -> host(s) lookup ───────────────────────────────────────
    ## For multi-host symbionts, state set = all their hosts
    tip_host_map <- split(Host_to_Symbiont_df$Host,
                          Host_to_Symbiont_df$Symbiont)
    tip_host_map <- lapply(tip_host_map, unique)

    ## ── build subtree list ────────────────────────────────────────────────────
    if (!is.null(Symbiont_tree$node.label) &&
        anyDuplicated(Symbiont_tree$node.label)) {
        Symbiont_tree$node.label <- paste0("Node",
            seq_len(Symbiont_tree$Nnode), "_", Symbiont_tree$node.label)
    } else if (is.null(Symbiont_tree$node.label)) {
        Symbiont_tree$node.label <- paste0("Node_", seq_len(Symbiont_tree$Nnode))
    }

    n_nodes     <- ape::Nnode(Symbiont_tree)
    subtree_list <- lapply(seq_len(n_nodes), function(node)
        castor::get_subtree_at_node(Symbiont_tree, node)$subtree)
    names(subtree_list) <- Symbiont_tree$node.label

    ## span filter
    total_span  <- castor::get_tree_span(Symbiont_tree,
                                          as_edge_count = FALSE)$max_distance
    span_cutoff <- span_fraction * total_span

    ## ── filter nodes ──────────────────────────────────────────────────────────
    node_meta <- data.frame(
        Node_ID = Symbiont_tree$node.label,
        stringsAsFactors = FALSE
    )
    node_meta$n_tips <- vapply(subtree_list, function(st)
        length(st$tip.label), integer(1))
    node_meta$span   <- vapply(subtree_list, function(st)
        castor::get_tree_span(st, as_edge_count = FALSE)$max_distance, numeric(1))

    tested <- node_meta[
        node_meta$n_tips >= min_tips &
        node_meta$span   <= span_cutoff, "Node_ID"]

    if (verbose)
        message(sprintf("[Occupancy] Testing %d / %d nodes (min_tips=%d, span_fraction=%.2f)",
                        length(tested), n_nodes, min_tips, span_fraction))

    ## ── worker: test one node ─────────────────────────────────────────────────
    test_node <- function(node_id) {

        st   <- subtree_list[[node_id]]
        tips <- st$tip.label

        ## build state sets for tips in this subtree
        state_sets <- tip_host_map[tips]
        ## tips with no association entry get all hosts (uninformative)
        missing <- vapply(state_sets, is.null, logical(1))
        all_host_states <- unique(Host_to_Symbiont_df$Host)
        state_sets[missing] <- list(all_host_states)

        hosts_rep <- unique(unlist(state_sets))
        n_hosts   <- length(hosts_rep)

        ## need at least 2 host states to have any parsimony cost
        if (n_hosts < 2L)
            return(data.frame(
                Node_ID              = node_id,
                N_Tips               = length(tips),
                N_Hosts_Represented  = n_hosts,
                Observed_Parsimony   = NA_real_,
                Max_Parsimony        = NA_real_,
                NPS                  = NA_real_,
                Parsimony_p          = NA_real_,
                Host_distribution    = paste(
                    names(sort(table(unlist(state_sets)), decreasing=TRUE)),
                    collapse=","),
                stringsAsFactors = FALSE
            ))

        ## observed parsimony
        use_weighted <- !is.null(host_dist_mat) &&
                        all(hosts_rep %in% rownames(host_dist_mat))
        obs_score <- if (use_weighted) {
            sub_hdm <- host_dist_mat[hosts_rep, hosts_rep, drop = FALSE]
            fitch_weighted(st, state_sets, sub_hdm)
        } else {
            fitch_parsimony(st, state_sets)
        }

        ## maximum possible parsimony: every tip has a different host from
        ## its neighbours -> upper bound = n_tips - 1
        max_score <- length(tips) - 1L

        ## NPS: 1 = perfect clustering, 0 = worst case
        nps <- if (max_score > 0) 1 - obs_score / max_score else NA_real_

        ## permutation null: shuffle tip-to-stateSet assignments
        ## preserve the multi-host structure of each tip (shuffle assignments
        ## as whole units so a multi-host tip stays multi-host)
        set.seed(seed + which(tested == node_id))
        tip_indices <- seq_along(tips)
        null_scores <- vapply(seq_len(permutations), function(i) {
            perm_map <- state_sets[sample(tip_indices)]
            names(perm_map) <- tips
            if (use_weighted)
                fitch_weighted(st, perm_map, sub_hdm)
            else
                fitch_parsimony(st, perm_map)
        }, numeric(1))

        ## lower parsimony = stronger clustering = signal
        p_val <- (sum(null_scores <= obs_score) + 1L) / (permutations + 1L)

        ## host distribution string for annotation
        host_tab  <- sort(table(unlist(state_sets)), decreasing = TRUE)
        host_dist_str <- paste(
            paste0(names(host_tab), "=", host_tab), collapse = ";")

        data.frame(
            Node_ID             = node_id,
            N_Tips              = length(tips),
            N_Hosts_Represented = n_hosts,
            Observed_Parsimony  = obs_score,
            Max_Parsimony       = max_score,
            NPS                 = round(nps, 4),
            Parsimony_p         = p_val,
            Host_distribution   = host_dist_str,
            stringsAsFactors    = FALSE
        )
    }

    ## ── parallel dispatch ─────────────────────────────────────────────────────
    if (verbose) message("[Occupancy] Running permutation tests across cores...")

    results_list <- parallel::mclapply(
        tested,
        test_node,
        mc.cores       = cores,
        mc.preschedule = FALSE
    )

    ## guard against failed workers
    failed <- vapply(results_list, inherits, logical(1), "try-error")
    if (any(failed))
        warning(sprintf("[Occupancy] %d nodes failed and will be NA.", sum(failed)))

    results_df <- do.call(rbind, results_list[!failed])
    results_df <- results_df[order(results_df$Parsimony_p,
                                   -results_df$NPS,
                                   na.last = TRUE), ]
    rownames(results_df) <- NULL

    ## ── BH correction across tested nodes ────────────────────────────────────
    ## Standard BH is appropriate here: nodes are not nested in the same way
    ## as the congruence scan (each tip appears in many subtrees but the
    ## parsimony score is computed independently per subtree, so dependence
    ## is positive and BH remains conservative).
    valid    <- !is.na(results_df$Parsimony_p)
    adj_p    <- rep(NA_real_, nrow(results_df))
    adj_p[valid] <- p.adjust(results_df$Parsimony_p[valid], method = "BH")
    results_df$Parsimony_p_BH <- round(adj_p, 4)

    if (verbose) {
        n_sig_raw <- sum(results_df$Parsimony_p    <= 0.05, na.rm = TRUE)
        n_sig_bh  <- sum(results_df$Parsimony_p_BH <= 0.05, na.rm = TRUE)
        message(sprintf(
            "[Occupancy] Significant nodes: %d (raw p<=0.05), %d (BH p<=0.05) out of %d tested",
            n_sig_raw, n_sig_bh, sum(valid)))
        if (n_sig_bh > 0L) {
            top <- head(results_df[!is.na(results_df$Parsimony_p_BH) &
                                    results_df$Parsimony_p_BH <= 0.05, ], 5L)
            for (i in seq_len(nrow(top))) {
                message(sprintf("    %s  NPS=%.3f  p=%.4f  BH=%.4f  n_tips=%d  n_hosts=%d",
                    top$Node_ID[i], top$NPS[i], top$Parsimony_p[i],
                    top$Parsimony_p_BH[i], top$N_Tips[i], top$N_Hosts_Represented[i]))
            }
        }
    }

    if (!is.null(output_fp)) {
        write.table(results_df, output_fp, sep = "\t",
                    quote = FALSE, row.names = FALSE)
        if (verbose) message(sprintf("[Occupancy] Results written to: %s", output_fp))
    }

    invisible(results_df)
}

# Entry point
main <- function() {
    args <- parse_args()

    if (!args$quiet) {
        message(" codiv - Host-Symbiont Codiversification Analysis\n",
                "---------------------------------------------")
        message("Host tree: ", args$host_tree)
        message("Symbiont tree: ", args$symbiont_tree)
        message("Associations file: ", args$associations)
        message("Output file: ", args$output)
        message("Methods: ", args$methods)
        message("Permutations: ", args$permutations)
        message("Cores: ", args$cores)
        message("Seed: ", args$seed)
        if (args$null_traversals > 0L && !args$skip_null)
            message("Null model traversals: ", args$null_traversals)
    }

    # Read input data
    Host_tree <- read_tree(args$host_tree, args$tree_format)
    Symbiont_tree <- read_tree(args$symbiont_tree, args$tree_format)

    ## read associations — auto-detect separator (tab or comma)
    if (!file.exists(args$associations))
        stop("Associations file not found: ", args$associations)
    assoc_first <- readLines(args$associations, n = 1L)
    assoc_sep   <- if (grepl("\t", assoc_first)) "\t" else ","
    assoc_df    <- read.table(args$associations, sep = assoc_sep, header = TRUE,
                                stringsAsFactors = FALSE, check.names = FALSE)
    
    ## normalise column names: trim whitespace, accept "host"/"symbiont" case-insensitively
    colnames(assoc_df) <- trimws(colnames(assoc_df))
    col_map <- setNames(colnames(assoc_df), tolower(colnames(assoc_df)))
    if ("host"     %in% names(col_map)) colnames(assoc_df)[colnames(assoc_df) == col_map[["host"]]]     <- "Host"
    if ("symbiont" %in% names(col_map)) colnames(assoc_df)[colnames(assoc_df) == col_map[["symbiont"]]] <- "Symbiont"
    
    methods     <- trimws(strsplit(args$methods, ",")[[1]])
    focus_hosts <- if (!is.null(args$focus_host) && nchar(args$focus_host) > 0L)
        trimws(strsplit(args$focus_host, ",")[[1]]) else character(0)

    # Ensure output directory exists
    out_dir <- dirname(args$output)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # Run codiv analysis
    results <- codiv(
          Host_tree           = Host_tree,
          Symbiont_tree       = Symbiont_tree,
          Host_to_Symbiont_df = assoc_df,
          min_hosts           = args$min_hosts,
          min_symbionts       = args$min_symbionts,
          max_symbionts       = args$max_symbionts,
          span_fraction       = args$span_fraction,
          permutations        = args$permutations,
          Save_fp             = args$output,
          seed                = args$seed,
          cores               = args$cores,
          methods             = methods,
          focus_host          = focus_hosts,
          subtree_features    = !args$no_subtree_feat,
          continue            = !args$no_continue,
          rarefy_depth        = args$rarefy_depth,
          verbose             = !args$quiet
    )

    ## permutation FDR — run unless --skip_fdr flag set
    if (!args$skip_fdr) {
        if (!args$quiet) message("\nRunning permutation FDR (", args$global_perms, " global permutations)...")
        results <- permutation_fdr(
            results_df          = results,
            Host_tree           = Host_tree,
            Symbiont_tree       = Symbiont_tree,
            Host_to_Symbiont_df = assoc_df,
            methods             = methods,
            global_perms        = as.integer(args$global_perms),
            alpha               = as.numeric(args$fdr_alpha),
            fdr_min_symbionts   = as.integer(args$fdr_min_symbionts),
            seed                = as.integer(args$seed),
            cores               = as.integer(args$cores),
            null_output_fp      = if (!is.null(args$fdr_null_output))
                                    args$fdr_null_output
                                  else
                                    sub("\\.[^.]+$", "_null_dist.tsv", args$output),
            verbose             = !args$quiet
        )
        ## overwrite TSV with FDR columns added
        write.table(results, args$output, sep = "\t", quote = FALSE, row.names = FALSE)
        if (!args$quiet) message("FDR results appended to: ", args$output)
    }

    ## clean up temp RDS directory now that everything is done
    tmp_dir <- attr(results, "tmp_dir")
    if (!is.null(tmp_dir) && dir.exists(tmp_dir)) {
        unlink(tmp_dir, recursive = TRUE)
        if (!args$quiet) message("Temp directory cleaned up.")
    }

    ## ── empirical null model (full traversal on permuted host trees) ──────────
    if (!args$skip_null && args$null_traversals > 0L) {
        if (!args$quiet)
            message(sprintf(
                "\nRunning empirical null model (%d permuted traversals)...",
                args$null_traversals))

        ## strip any trailing extension from output path to build prefix
        out_prefix <- sub("\\.[^.]+$", "", args$output)

        null_model_traversal(
            results_df          = results,
            Host_tree           = Host_tree,
            Symbiont_tree       = Symbiont_tree,
            Host_to_Symbiont_df = assoc_df,
            n_traversals        = as.integer(args$null_traversals),
            methods             = methods,
            alpha               = as.numeric(args$null_alpha),
            r_threshold         = as.numeric(args$null_r_threshold),
            output_prefix       = out_prefix,
            min_hosts           = args$min_hosts,
            min_symbionts       = args$min_symbionts,
            max_symbionts       = args$max_symbionts,
            span_fraction       = args$span_fraction,
            permutations        = args$permutations,
            seed                = args$seed,
            cores               = args$cores,
            verbose             = !args$quiet
        )
    }

    ## ── occupancy scan ────────────────────────────────────────────────────────
    if (args$occupancy) {
        if (!args$quiet)
            message("\nRunning occupancy scan (parsimony-based phylogenetic signal)...")

        occ_fp <- if (!is.null(args$occ_output)) {
            args$occ_output
        } else {
            sub("\\.[^.]+$", "_occupancy.tsv", args$output)
        }

        occupancy_scan(
            Symbiont_tree       = Symbiont_tree,
            Host_to_Symbiont_df = assoc_df,
            Host_tree           = Host_tree,   ## weighted transitions
            permutations        = as.integer(args$occ_permutations),
            min_tips            = as.integer(args$occ_min_tips),
            span_fraction       = 1.0,         ## test all depths
            seed                = as.integer(args$seed),
            cores               = as.integer(args$cores),
            output_fp           = occ_fp,
            verbose             = !args$quiet
        )
    }
}

# =============================================================================
# Generalised permutation-based FDR
#
# Strategy
# --------
# Covers all three methods: hommola (higher r = stronger signal),
# paco (lower ss = stronger signal), parafit (higher GlobalFit = stronger).
#
# For each of B global permutations:
#   1. Shuffle the ENTIRE Host_to_Symbiont_df once — breaks all real signal
#      simultaneously across every node, preserving inter-node correlation
#      structure in the null (critical for nested nodes).
#   2. Re-run the raw statistic (NOT the inner permutation loop) for every
#      tested node under that global shuffle.
#   3. Record, per method, the most extreme statistic seen across all nodes.
#
# The resulting B-length vector is the global null for the best result
# achievable by chance across the whole experiment.
#
# FDR(node_i, method) =
#   (B * mean(global_null_extreme >= stat_i)) /   <- expected false positives
#   max(1, sum(obs_stat >= stat_i))                <- observed positives
#
# "extreme" is directional per method:
#   hommola  : larger r is more extreme  -> use max across nodes
#   paco     : smaller ss is more extreme -> use min across nodes (negate to unify)
#   parafit  : larger GlobalFit is more extreme -> use max across nodes
#
# Parameters
# ----------
# results_df          : data frame returned by codiv()
# Host_tree           : host phylo object
# Symbiont_tree       : symbiont phylo object
# Host_to_Symbiont_df : full association data frame
# subtree_list        : named list of symbiont subtrees; rebuilt if NULL
# methods             : which methods to compute FDR for (auto-detected if NULL)
# global_perms        : number of global permutations (recommend >= 999)
# alpha               : FDR threshold for the printed summary
# seed                : RNG seed
# cores               : parallel cores
# verbose             : print progress
#
# Returns
# -------
# results_df with new columns:
#   Hommola_Perm_FDR    (if hommola results present)
#   PACo_Perm_FDR       (if paco results present)
#   ParaFit_Perm_FDR    (if parafit results present)
#   <Method>_Null_quantile  for each method
# =============================================================================

permutation_fdr <- function(results_df,
                            Host_tree,
                            Symbiont_tree,
                            Host_to_Symbiont_df,
                            subtree_list      = NULL,
                            methods           = NULL,
                            global_perms      = 999L,
                            alpha             = 0.05,
                            fdr_min_symbionts = 10L,   ## exclude tiny nodes from global max
                            seed              = 8675309L,
                            cores             = 4L,
                            null_output_fp    = NULL,  ## path to save raw null distribution
                            verbose           = TRUE) {

  ## ── method schema ─────────────────────────────────────────────────────────
  ## Each entry defines:
  ##   obs_col   : column in results_df holding the observed statistic
  ##   higher_is_extreme : TRUE  -> larger stat = stronger signal (hommola, parafit)
  ##                       FALSE -> smaller stat = stronger signal (paco)
  ##   compute   : function(host_dist, col_sym, col_h2s, host_sub, col_sym_tree)
  ##               -> scalar statistic under a shuffled association table
  method_schema <- list(

    hommola = list(
      obs_col          = "Rarefied_Hommola_r",
      fdr_col          = "Hommola_Perm_FDR",
      q_col            = "Hommola_Null_quantile",
      higher_is_extreme = TRUE,
      compute          = function(host_dist, rar_sym_d, rar_h2s, host_sub, rar_sym) {
        hommola(host_dist, rar_sym_d, rar_h2s)
      }
    ),

    paco = list(
      obs_col          = "Rarefied_PACo_ss",
      fdr_col          = "PACo_Perm_FDR",
      q_col            = "PACo_Null_quantile",
      higher_is_extreme = FALSE,
      compute          = function(host_dist, rar_sym_d, rar_h2s, host_coords, sym_coords) {
        valid_hosts <- intersect(rar_h2s$Host,     rownames(host_coords))
        valid_syms  <- intersect(rar_h2s$Symbiont, rownames(sym_coords))
        rar_h2s_v   <- rar_h2s[rar_h2s$Host %in% valid_hosts &
                                rar_h2s$Symbiont %in% valid_syms, , drop = FALSE]
        if (nrow(rar_h2s_v) == 0L ||
            length(unique(rar_h2s_v$Host)) < 2L ||
            length(unique(rar_h2s_v$Symbiont)) < 2L) return(Inf)
        assoc_mat   <- host_symbiont_links(rar_h2s_v)
        hc <- host_coords[rownames(assoc_mat), , drop = FALSE]
        sc <- sym_coords[ colnames(assoc_mat), , drop = FALSE]
        .paco_ss(hc, sc, assoc_mat)
      }
    ),

    parafit = list(
      obs_col          = "Rarefied_ParaFitGlobal",
      fdr_col          = "ParaFit_Perm_FDR",
      q_col            = "ParaFit_Null_quantile",
      higher_is_extreme = TRUE,
      compute          = function(host_dist, rar_sym_d, rar_h2s, host_coords, sym_coords) {
        valid_hosts <- intersect(rar_h2s$Host,     rownames(as.matrix(host_dist)))
        valid_syms  <- intersect(rar_h2s$Symbiont, rownames(as.matrix(rar_sym_d)))
        rar_h2s_v   <- rar_h2s[rar_h2s$Host %in% valid_hosts &
                                rar_h2s$Symbiont %in% valid_syms, , drop = FALSE]
        if (nrow(rar_h2s_v) == 0L ||
            length(unique(rar_h2s_v$Host)) < 2L ||
            length(unique(rar_h2s_v$Symbiont)) < 2L) return(-Inf)
        hd <- as.matrix(host_dist)[valid_hosts, valid_hosts, drop = FALSE]
        sd <- as.matrix(rar_sym_d)[valid_syms,  valid_syms,  drop = FALSE]
        res <- ape::parafit(
          as.dist(hd), as.dist(sd),
          host_symbiont_links(rar_h2s_v),
          nperm      = 1L,
          test.links = FALSE,
          correction = "cailliez",
          silent     = TRUE
        )
        res$ParaFitGlobal
      }
    )
  )

  ## ── auto-detect which methods are present in results_df ───────────────────
  if (is.null(methods)) {
    methods <- names(method_schema)[vapply(names(method_schema), function(m) {
      col <- method_schema[[m]]$obs_col
      col %in% colnames(results_df) && any(!is.na(results_df[[col]]))
    }, logical(1))]
  }

  if (length(methods) == 0L)
    stop("No recognised method columns found in results_df.")

  if (verbose)
    message(sprintf("Methods detected: %s", paste(methods, collapse = ", ")))

  ## ── tested nodes: must have at least one method result ────────────────────
  obs_cols    <- vapply(methods, function(m) method_schema[[m]]$obs_col, character(1))
  has_any_result <- apply(
    results_df[, obs_cols, drop = FALSE], 1,
    function(row) any(!is.na(row))
  )
  tested_nodes <- results_df$Node_ID[has_any_result]

  if (length(tested_nodes) == 0L)
    stop("No completed nodes found in results_df. Run codiv() first.")

  if (verbose)
    message(sprintf("Nodes to assess: %d  Global permutations: %d",
                    length(tested_nodes), global_perms))

  ## ── rebuild subtree list if not supplied ──────────────────────────────────
  if (is.null(subtree_list)) {
    if (verbose) message("Rebuilding subtree list ...")

    ## apply the same node-label prepending as codiv() so names match Node_IDs
    ## in results_df. Three cases must mirror codiv() exactly:
    ##   1. labels exist and are duplicated  -> prepend "NodeN_"
    ##   2. labels are NULL                  -> assign "Node_N"
    ##   3. labels exist and are unique      -> leave as-is (codiv() does nothing)
    if (!is.null(Symbiont_tree$node.label) &&
        anyDuplicated(Symbiont_tree$node.label)) {
      Symbiont_tree$node.label <- paste0("Node",
        seq_len(Symbiont_tree$Nnode), "_", Symbiont_tree$node.label)
    } else if (is.null(Symbiont_tree$node.label)) {
      Symbiont_tree$node.label <- paste0("Node_", seq_len(Symbiont_tree$Nnode))
    }
    ## case 3: unique labels present — no change needed

    subtree_list <- lapply(seq_len(ape::Nnode(Symbiont_tree)), function(x)
      castor::get_subtree_at_node(Symbiont_tree, x)$subtree)
    names(subtree_list) <- Symbiont_tree$node.label
  }

  ## ── fast_dist defined locally so lapply workers can see it ─────────────
  fast_dist <- function(tree) {
    n   <- length(tree$tip.label)
    mat <- castor::get_all_pairwise_distances(tree,
               only_clades = seq_len(n), as_edge_counts = FALSE)
    rownames(mat) <- colnames(mat) <- tree$tip.label
    as.dist(mat)
  }

  ## ── pre-compute per-node inputs (done once, shared across global perms) ───
  if (verbose) message("Pre-computing node inputs ...")

  node_data <- setNames(
    lapply(tested_nodes, function(node_id) {

      ## NULL if subtree not found (node label mismatch)
      sym_sub  <- subtree_list[[node_id]]
      if (is.null(sym_sub)) return(NULL)

      h2s_df   <- Host_to_Symbiont_df[
        Host_to_Symbiont_df$Symbiont %in% sym_sub$tip.label, ]

      ## guard: keep.tip can return NULL if no hosts match
      valid_hosts <- intersect(unique(h2s_df$Host), Host_tree$tip.label)
      if (length(valid_hosts) < 2L) return(NULL)
      host_sub <- ape::keep.tip(Host_tree, valid_hosts)
      if (is.null(host_sub) || length(host_sub$tip.label) < 2L) return(NULL)

      host_dist <- fast_dist(host_sub)

      ## rarefaction: equalise per-host symbiont counts (auto depth = min per-host)
      per_host_n <- table(h2s_df$Host)
      rarefy_n   <- as.integer(min(per_host_n))   ## always auto in FDR context
      set.seed(seed + which(tested_nodes == node_id))
      rar_syms   <- unlist(lapply(names(per_host_n), function(h) {
          syms <- h2s_df$Symbiont[h2s_df$Host == h]
          if (length(syms) <= rarefy_n) syms else sample(syms, rarefy_n)
      }))
      rar_sym    <- ape::keep.tip(sym_sub, rar_syms)
      if (is.null(rar_sym) || length(rar_sym$tip.label) < 2L) return(NULL)
      rar_sym_d  <- fast_dist(rar_sym)
      rar_h2s    <- h2s_df[h2s_df$Symbiont %in% rar_syms, ]

      ## pre-compute PCoA coords once — workers just do cheap Procrustes
      host_mat   <- as.matrix(host_dist)
      host_pcoa  <- tryCatch(
        vegan::wcmdscale(host_mat, eig = TRUE, add = "cailliez"),
        error = function(e) NULL)
      if (is.null(host_pcoa)) return(NULL)
      host_coords <- host_pcoa$points[, host_pcoa$eig > 0, drop = FALSE]
      rownames(host_coords) <- rownames(host_mat)

      sym_mat    <- as.matrix(rar_sym_d)
      sym_pcoa   <- tryCatch(
        vegan::wcmdscale(sym_mat, eig = TRUE, add = "cailliez"),
        error = function(e) NULL)
      if (is.null(sym_pcoa)) return(NULL)
      sym_coords  <- sym_pcoa$points[, sym_pcoa$eig > 0, drop = FALSE]
      rownames(sym_coords) <- rownames(sym_mat)

      list(
        host_dist   = host_dist,
        host_sub    = host_sub,
        host_coords = host_coords,
        rar_sym     = rar_sym,
        rar_sym_d   = rar_sym_d,
        sym_coords  = sym_coords,
        rar_h2s     = rar_h2s
      )
    }),
    tested_nodes  ## name the list so Filter preserves node IDs
  )

  ## drop nodes where inputs could not be computed
  node_data    <- Filter(Negate(is.null), node_data)
  tested_nodes <- names(node_data)
  if (verbose) message(sprintf("Valid nodes for FDR: %d / %d",
                               length(tested_nodes),
                               nrow(results_df[has_any_result, ])))
  if (length(tested_nodes) == 0L) stop("No valid nodes for FDR after input checks.")

  ## ── one global permutation: returns named vector of extreme stats per method
  ##
  ## IMPORTANT: the global maximum is computed only over nodes with
  ## N_Symbionts >= fdr_min_symbionts. Tiny nodes (e.g. N=3) trivially produce
  ## r=1.0 under random label permutations, which would inflate the global null
  ## to ~1.0 in every permutation and make FDR=1 for all nodes.
  ## Nodes below the threshold still receive FDR scores (compared against the
  ## null built from larger nodes) but do not contribute to the null maximum.
  n_sym_vec  <- setNames(results_df$N_Symbionts, results_df$Node_ID)
  fdr_nodes  <- tested_nodes[
    !is.na(n_sym_vec[tested_nodes]) &
    n_sym_vec[tested_nodes] >= fdr_min_symbionts
  ]
  if (length(fdr_nodes) == 0L) {
    warning(sprintf(
      "[FDR] No nodes with N_Symbionts >= %d. Lowering fdr_min_symbionts or ",
      fdr_min_symbionts,
      "check your data. Using all tested nodes instead."))
    fdr_nodes <- tested_nodes
  }
  if (verbose)
    message(sprintf(
      "[FDR] Nodes contributing to global null (N_Symbionts >= %d): %d / %d",
      fdr_min_symbionts, length(fdr_nodes), length(tested_nodes)))

  run_global_perm <- function(perm_idx) {

    set.seed(seed + perm_idx)

    ## single global shuffle of HOST labels — breaks all real signal while
    ## keeping the symbiont tree and rarefied sym set intact (they are fixed).
    ## Shuffling Symbiont would reference OTUs absent from rar_sym$tip.label
    ## and produce empty node_shuffled after subsetting, making every stat NA.
    shuffled_h2s        <- Host_to_Symbiont_df
    shuffled_h2s$Host   <- sample(shuffled_h2s$Host)

    ## for each method, collect statistic across all nodes then take extreme
    vapply(methods, function(m) {

      schema  <- method_schema[[m]]
      hi      <- schema$higher_is_extreme

      node_stats <- vapply(fdr_nodes, function(node_id) {

        nd <- node_data[[node_id]]

        ## subset global shuffle to symbionts in this node's rarefied sym set
        node_shuffled <- shuffled_h2s[
          shuffled_h2s$Symbiont %in% colnames(as.matrix(nd$rar_sym_d)), ]

        ## need enough diversity to compute the statistic
        if (length(unique(node_shuffled$Host))     < 2L ||
            length(unique(node_shuffled$Symbiont)) < 2L)
          return(if (hi) -Inf else Inf)

        tryCatch(
          schema$compute(nd$host_dist, nd$rar_sym_d,
                         node_shuffled, nd$host_coords, nd$sym_coords),
          error = function(e) if (hi) -Inf else Inf
        )

      }, numeric(1))

      ## most extreme stat across all nodes in this global permutation
      if (hi) max(node_stats, na.rm = TRUE) else min(node_stats, na.rm = TRUE)

    }, numeric(1))  ## returns named numeric(length(methods))
  }

  ## ── dispatch ──────────────────────────────────────────────────────────────
  if (verbose) message("Running global permutations ...")

  perm_results <- parallel::mclapply(
    seq_len(global_perms),
    run_global_perm,
    mc.cores       = cores,
    mc.preschedule = FALSE
  )

  ## Guard: mclapply returns try-error objects for failed workers instead of
  ## numeric vectors. Filter these out before rbind to avoid a corrupt matrix.
  failed_perms <- vapply(perm_results, inherits, logical(1), "try-error")
  if (any(failed_perms)) {
    warning(sprintf(
      "[FDR] %d / %d global permutations failed and will be excluded.",
      sum(failed_perms), global_perms))
    perm_results <- perm_results[!failed_perms]
  }
  if (length(perm_results) == 0L)
    stop("[FDR] All global permutations failed. Cannot compute FDR.")

  ## perm_results is a list of named numeric vectors; reshape to matrix
  ## rows = successful permutations, cols = methods
  null_mat <- do.call(rbind, perm_results)
  colnames(null_mat) <- methods
  n_valid_perms <- nrow(null_mat)   ## may be < global_perms if some failed

  ## ── save raw null distribution to disk ────────────────────────────────────
  ## This is the 999 x n_methods matrix of per-permutation global maxima.
  ## Without saving it, the null is lost after this function returns and
  ## FDR results cannot be recomputed or inspected post-hoc.
  null_df <- as.data.frame(null_mat)
  null_df <- cbind(data.frame(permutation = seq_len(n_valid_perms)), null_df)
  if (!is.null(null_output_fp)) {
    write.table(null_df, null_output_fp, sep = "\t", quote = FALSE, row.names = FALSE)
    if (verbose) message(sprintf("[FDR] Raw null distribution saved to: %s", null_output_fp))
  } else {
    ## default: save alongside the results file, derived from Save_fp if possible
    ## Caller can also pass null_output_fp explicitly to control the path
    if (verbose) message("[FDR] Tip: pass null_output_fp to permutation_fdr() to save the null distribution.")
  }
  ## always attach to returned object as attribute so caller can access it
  ## even if not written to disk
  attr(null_df, "fdr_min_symbionts") <- fdr_min_symbionts
  attr(null_df, "n_fdr_nodes")       <- length(fdr_nodes)

  ## ── compute FDR and null quantile per method per node ─────────────────────
  for (m in methods) {

    schema  <- method_schema[[m]]
    hi      <- schema$higher_is_extreme
    obs_col <- schema$obs_col
    obs_vec <- as.numeric(results_df[[obs_col]])
    names(obs_vec) <- results_df$Node_ID

    null_extreme <- null_mat[, m]   ## length global_perms

    fdr_col <- vapply(results_df$Node_ID, function(node_id) {

      s_i <- obs_vec[node_id]
      if (is.na(s_i)) return(NA_real_)

      ## Storey-style FDR:
      ##   FDR(i) = E[false positives at threshold s_i] / observed positives at s_i
      ##          = (# null draws >= s_i / n_valid_perms) / (# obs >= s_i / n_nodes)
      ## Simplified to:
      ##   = (# null draws as extreme as s_i / n_valid_perms)
      ##     / (# observed as extreme as s_i / n_tested_nodes)
      ## Use n_valid_perms (actual successes) not global_perms (requested),
      ## so failed workers do not inflate the FDR numerator.
      n_tested <- length(obs_vec[!is.na(obs_vec)])
      if (hi) {
        prop_null <- sum(null_extreme >= s_i, na.rm = TRUE) / n_valid_perms
        prop_obs  <- sum(obs_vec      >= s_i, na.rm = TRUE) / max(1L, n_tested)
      } else {
        prop_null <- sum(null_extreme <= s_i, na.rm = TRUE) / n_valid_perms
        prop_obs  <- sum(obs_vec      <= s_i, na.rm = TRUE) / max(1L, n_tested)
      }

      min(1.0, prop_null / max(.Machine$double.eps, prop_obs))

    }, numeric(1))

    q_col <- vapply(results_df$Node_ID, function(node_id) {

      s_i <- obs_vec[node_id]
      if (is.na(s_i)) return(NA_real_)

      ## fraction of valid null draws WEAKER than this node's stat
      ## -> 1.0 = node beats 100% of null (best possible)
      if (hi) mean(null_extreme <= s_i, na.rm = TRUE)
      else    mean(null_extreme >= s_i, na.rm = TRUE)

    }, numeric(1))

    results_df[[schema$fdr_col]] <- round(fdr_col, 4)
    results_df[[schema$q_col]]   <- round(q_col,   4)
  }

  ## ── summary ───────────────────────────────────────────────────────────────
  if (verbose) {
    for (m in methods) {
      fdr_col <- results_df[[method_schema[[m]]$fdr_col]]
      n_sig   <- sum(fdr_col <= alpha, na.rm = TRUE)
      message(sprintf("[%s] Nodes at FDR <= %.2f: %d / %d",
                      m, alpha, n_sig, length(tested_nodes)))

      if (n_sig > 0L) {
        obs_col <- method_schema[[m]]$obs_col
        hi      <- method_schema[[m]]$higher_is_extreme
        top     <- results_df[!is.na(fdr_col) & fdr_col <= alpha, ]
        top     <- top[order(top[[method_schema[[m]]$fdr_col]],
                             if (hi) -top[[obs_col]] else top[[obs_col]]), ]
        for (k in seq_len(min(5L, nrow(top)))) {
          message(sprintf("    %s  stat=%.4f  FDR=%.4f",
                          top$Node_ID[k],
                          top[[obs_col]][k],
                          top[[method_schema[[m]]$fdr_col]][k]))
        }
      }
    }
  }

  attr(results_df, "null_distribution") <- null_df
  return(results_df)
}

## Fast patristic distance — castor C++ backend, returns dist object
## Orders of magnitude faster than adephylo::distTips for large trees
## NOTE: also defined inside codiv() and permutation_fdr() for mclapply worker visibility
fast_dist <- function(tree) {
    n <- length(tree$tip.label)
    mat <- castor::get_all_pairwise_distances(tree,
                only_clades = seq_len(n),
                as_edge_counts = FALSE)
    rownames(mat) <- colnames(mat) <- tree$tip.label
    as.dist(mat)
}

main()