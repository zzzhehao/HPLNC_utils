#' Consensus Clustering
#' 
#' Perform consensus clustering (CC) with any base method.
#' 
#' @param source_matrix A matrix. Multivariate dataset on which to perform the CC.
#' @param k Numeric. Hyperparameter. Number of clusters to produce. Will be passed to [cutree()] and if \code{base_method = "kmean"}, also to [kmeans()].
#' @param reps Number of subsampling to build consensus matrix.
#' @param subset Proportion of source_matrix to subsample in each replication.
#' @param base_method Base clustering method. Either \code{"hclust"} for hierarchical clustering (\code{hclust(method = "ward.D2")}) or \code{"kmeans"} for k-means clustering (\code{kmeans()}). Will be ignored if \code{base_call} is provided.
#' @param base_call A function call for base clustering method, available variables are the subsampled matrix \code{m}, calculated distance matrix \code{dist.matrix}, and everything passed to \code{...}.
#' @param dist_method Distance calculation method shortcut to generate the calls. Will be ignored if \code{base_method = "kmeans"} or if \code{dist_call} is provided. Accept \code{"euclidean"}, \code{"manhattan"}, \code{"cor_pearson"}, and \code{"cor_spearman"}.
#' @param dist_call A function call for calculating distance that is passed to base clustering method. Available variables are the subsampled matrix \code{m}, and everything passed to \code{...}.
#' @param seed Numeric vector of length equal to or greater than \code{reps}. The seed will be sequentially applied to each subsampling. Passing the \code{seed} slot product here will ensure reproducibility of existing results.
#' @param ... Additional arguments for custom function calls. Either for \code{base_call} or \code{dist_call}. You need to add \code{...} in your function call to receive the arguments.
#' 
#' @export
consensus_clustering <- new_class(
    "consensus_clustering",
    properties = list(
        k = class_integer,
        reps = class_integer,
        subset = class_numeric,
        dist_method = class_character,
        dist_call = class_call,
        seed = class_numeric,
        base_method = class_character, 
        base_call = class_call, 
        summary_call = class_call,
        sample_name = class_character,
        assignment = class_data.frame,
        silhouette_width_avg = class_numeric,
        source_matrix = class_numeric,
        cc_matrix = class_numeric,
        M = class_numeric,
        I = class_numeric,
        heatmap = new_S3_class("ggplot"),
        silhouette_plot = new_S3_class("ggplot"),
        summary_clustering = class_any
    ),
    validator = function(self) {
        if (!is.matrix(self@source_matrix)) {
            "@source_matrix must be a matrix."
        }
    },
    constructor = function(
        source_matrix, 
        k, 
        reps = 500, 
        subset = 0.8, 
        base_method = "hclust",
        base_call = NULL,
        dist_method = "cor_pearson",
        dist_call = NULL,
        seed = NULL,
        ...
    ) {  
        # input
        samples <- rownames(source_matrix)
        n <- length(samples)
        
        M <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))
        I <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))

        matrices <- list(M = M, I = I)
      
        # generate seeds
        if (is.null(seed)) {
            seed <- round(runif(reps, 0, 9999), 0)
        } else {
            if (length(seed) < reps) {stop("Seed has to be equal to greater than reps.")}
        }
      
        # distance call
        if (dist_method == "euclidean") {
            dist_call <- expr(dist(m))
        } else if (dist_method == "manhattan") {
            dist_call <- expr(dist(m, method = "manhattan"))
        } else if (dist_method == "cor_pearson") {
            dist_call <- expr(as.dist(1 - cor(t(m), method = "pearson"))) 
        } else if (dist_method == "cor_spearman") {
            dist_call <- expr(as.dist(1 - cor(t(m), method = "spearman"))) 
        } else if (!is.null(dist_call)) {
            cli::cli_alert_warning("{.val {dist_method}} is not a short cut method, using specified expression: {.val {dist_call}}")
        } else if (is.null(dist_call)) {
            cli::cli_abort("{.val {dist_method}} is not a short cut method, please specify the expression for distance calculation.")
        }
      
        # clustering call
        if (is.null(base_call)) {
            if (base_method == "hclust") {
                base_call <- expr(hclust(dist.matrix, method = "ward.D2"))
            } else if (base_method == "kmeans") {
                base_call <- expr(kmeans(m, centers = k))
                cli::cli_alert_warning("k-mean clustering using strictly euclidean distance, anything else will be ignored.")
            }
        }
        
        # core loop
        pb.sub <- cli::cli_progress_bar(total = reps, format = paste0("Permuting iteration", " {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}"), clear = T)
        matrices.perm <- reduce(seq(reps), .init = matrices, .f = \(matrices, r){
            set.seed(seed[r])
            sub.idx <- sample(1:n, size = floor(n * subset), replace = FALSE)

            m <- source_matrix[sub.idx, ]
            dist.matrix <- eval_tidy(dist_call)
            cluster.obj <- eval_tidy(base_call)
            cl <- cutree(cluster.obj, k = k)

            names(cl) <- rownames(source_matrix)[sub.idx]
            connectivity <- outer(cl, cl, "==")
            
            sub.names <- names(cl)
            matrices[["M"]][sub.names, sub.names] <- matrices[["M"]][sub.names, sub.names] + connectivity
            matrices[["I"]][sub.names, sub.names] <- matrices[["I"]][sub.names, sub.names] + 1
            cli::cli_progress_update(id = pb.sub)

            return(matrices)
        })
        cli::cli_progress_done(pb.sub)

        consensus.matrix <- matrices.perm[["M"]] / matrices.perm[["I"]]
        consensus.matrix[is.na(consensus.matrix)] <- 0
      
        # Here merge with summary function
        summary_call <- expr(hclust(as.dist(1 - matrix), method = "ward.D2"))
        summary.ls <- consensus_clustering_summarize(consensus.matrix, k, source_matrix, summary_call)
      
        # silhouette 
        silhouette_dist_call <- do.call(substitute, list(dist_call, list(m = quote(source_matrix))))
        sil <- cluster::silhouette(as.numeric(summary.ls$cl$cc_cluster), eval_tidy(silhouette_dist_call)) %>% as.data.frame()
        sil_data <- sil %>%
            arrange(cluster, desc(sil_width)) %>%
            mutate(id = row_number())
        avg_width <- mean(sil_data$sil_width)

        # silhouette plot
        silhouette_plot <- ggplot(sil_data, aes(x = id, y = sil_width, fill = cluster)) +
            geom_col(width = 1) +  
            # scale_fill_brewer(palette = "Set1") + 
            geom_hline(yintercept = avg_width, linetype = "dashed", color = "red") +
            labs(
                title = paste("Silhouette Plot for k =", k),
                subtitle = paste("Average Silhouette Width =", round(avg_width, 2)),
                x = "Samples",
                y = "Silhouette Width",
                fill = "Cluster"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_blank(), 
                axis.ticks.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        
        new_object(
            S7_object(),
            source_matrix = source_matrix,
            cc_matrix = consensus.matrix,
            k = as.integer(k),
            reps = as.integer(reps),
            subset = subset,
            dist_method = dist_method,
            dist_call = dist_call,
            seed = seed,
            base_method = base_method,
            base_call = base_call,
            summary_call = summary_call, 
            M = matrices.perm[["M"]],
            I = matrices.perm[["I"]],
            sample_name = samples,
            heatmap = summary.ls$heatmap,
            assignment = summary.ls$cl,
            summary_clustering = summary.ls$fc,
            silhouette_plot = silhouette_plot,
            silhouette_width_avg = avg_width
        )
    }
)


#' Generate Heatmap with Dendrogram in ggplot2 system.
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @import ggplot2 
gg_heatmap_hclust <- function(
    hc.obj,
    matrix,
    heat_name = "Heat value",
    title = "",
    font_size = 8,
    dendro.stretch = c(0, 0.6),
    palette = NULL,
    zeroaswhite = T,
    grid = T,
    plot_labs = list(labs(x = NULL, y = NULL)),
    tip_length = 0,
    themes = list(theme()),
    heat_lab = NULL,
    args = list()
){
    dend_order <- hc.obj$labels[hc.obj$order]
    heatmap_data <- matrix %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column("row_id") %>% 
        pivot_longer(-row_id, names_to = "col_id", values_to = "heat") %>% 
        mutate(
            row_id = factor(row_id, levels = dend_order),
            col_id = factor(col_id, levels = dend_order) 
        )

    p.dendro <- hc.obj %>% 
        as.dendrogram() %>%
        ggdendro::dendro_data() %>%
        ggdendro::segment() %>%
        ggplot() +
        geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) + 
        geom_segment(aes(x = 0, y = x, xend = tip_length, yend = x)) +
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = dendro.stretch) + 
        theme_void() + 
        theme(
            plot.margin = margin(0, 0, 0, 0.01),
            aspect.ratio = 5
        ) 

    p.heatmap <- heatmap_data %>%
        ggplot(aes(x = col_id, y = row_id, fill = heat, data_id = row_id, tooltip = row_id))
        # ggplot(aes(x = col_id, y = row_id, fill = heat, data_id = row_id))

    if (grid) {
        p.heatmap <- p.heatmap + ggiraph::geom_tile_interactive(color = "grey85", linewidth = 0.1) 
    }

    if (is.null(palette)) {
        palette <- rcartocolor::carto_pal(n = 15, name = "Sunset")
    }

    if (zeroaswhite) {
        clr <- c("white", palette)
        palette.length <- length(palette) + 1
    } else {
        clr <- palette
        palette.length <- length(palette)
    }

    if (is.null(heat_lab)) {
        heat_lab <- c("0", round(max(matrix), 0))
    }

    p.heatmap <- p.heatmap +
        scale_fill_gradientn(
            name = heat_name,
            colors = clr, 
            values = c(0, seq(0.01, max(matrix), length.out = palette.length)),
            breaks = c(0, max(matrix)),
            labels = heat_lab,
            guide = guide_colorbar(
                title.position = "top",   
                title.hjust = 0.5,        
                label.position = "bottom",
                ticks = F,
                draw.ulim = F,
                draw.llim = F
            ),
            limits = c(0, max(matrix)),
        ) +
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_discrete(expand = c(0, 0), position = "right") + 
        coord_fixed() +
        labs(
            title = title,
            x = NULL,
            y = NULL
        ) + 
        plot_labs +
        theme_minimal() + 
        theme(
            legend.position = "bottom",
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = font_size, margin = margin(l = 5), hjust = 0), 
            plot.margin = margin(l = 0, r = 5, t = 5, b = 5)
        ) +
        themes +
        args

    # merge
    pp <- p.dendro + p.heatmap + patchwork::plot_layout(widths = c(1, 5))
    return(pp)
}

#' Summarize Consensus Clustering Result 
#' 
#' @param mv.data Multivariate dataset used for [consensus_clustering()].
#' @param matrix Consensus matrix produced by [consensus_clustering()].
#' @param k Numerical. Number of groups (k) to be passed to \code{summary.call} and [stats::cutree()].
#' @param summary.call Final clustering call to group samples according to \code{matrix}.
#' @param heat_name Character. Title of the legend for the tile color.
#' @param title Character or an expression. Title of the heatmap. Will be passed to [gg_heatmap_hclust()].
#' @param plot_labs [ggplot2::labs()] call wrapped in a list. Will be passed to [gg_heatmap_hclust()].
#' @param heat_lab A character vector in length of two. Labels for the min and max value of the color tile legend. Will be passed to [gg_heatmap_hclust()].
#' @param ... Further argument passed to [gg_heatmap_hclust()].
#' @importFrom rlang eval_tidy
consensus_clustering_summarize <- function(
    matrix,
    k,
    mv.data,
    summary.call,
    heat_name = "Consensus index",
    title = expr(paste0("Consensus clustering (k = ", k, ")")),
    plot_labs = list(labs(subtitle = "Method: hierarchical clustering (ward.D2)")),
    heat_lab = c("Least similar", "Most similar"),
    ...
){
    final.c <- eval_tidy(summary.call)
    cl <- cutree(final.c, k) %>% as.factor() %>% data.frame(label = names(.), cc_cluster = .)
    heatmap <- gg_heatmap_hclust(
        final.c, 
        matrix,
        heat_name = heat_name,
        title = ifelse(is.character(title), title, eval_tidy(title)),
        plot_labs = plot_labs,
        heat_lab = heat_lab,
        ...
    )
    
    cluster::silhouette(k, mv.data)
  
    res <- list(k = k, matrix = matrix, fc = final.c, cl = cl, heatmap = heatmap)
    return(res)
}

#' Perform Consensus Clustering for a Range of K
#' 
#' @importFrom rlang eval_tidy
#' @importFrom purrr map
#' @importFrom purrr map_vec
#' @import ggplot2
consensus_clustering_mk <- function(
    source_matrix, 
    kmin, 
    kmax, 
    ...
){
    # core loop for CC
    pb.cc <- cli::cli_progress_bar(total = length(kmin:kmax), format = paste0("Running consensus clustering ", " {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}"), clear = T)
    consensus.matrix.ls <- map(kmin:kmax, \(k){
        cc_obj <- consensus_clustering(source_matrix, k, ...)
        cli::cli_progress_update(id = pb.cc)
        return(cc_obj)
    })
    cli::cli_progress_done(id = pb.cc)
    names(consensus.matrix.ls) <- kmin:kmax

    # consensus CDF
    cdf.palette <- rcartocolor::carto_pal(n = length(kmin:kmax), name = "Fall")

    cdf_fun.ls <- map(consensus.matrix.ls, \(m) {
        matrix <- m@cc_matrix
        i <- match(m@k, kmin:kmax)
        stat_function(fun = ecdf(matrix[lower.tri(matrix)]), size = 1, color = cdf.palette[[i]])
    })

    cdf <- ggplot() +
        cdf_fun.ls +
        theme_minimal() +
        theme(
            aspect.ratio = 1,
            panel.border = element_rect(color = "black", fill = "transparent", linewidth = 1),
            panel.grid.minor.x = element_blank()
        ) +
        labs(title = "Consensus CDF", x = "Consensus index", y = "CDF")

    # Delta area
    delta.area <- map_vec(consensus.matrix.ls, \(m) {
        matrix <- m@cc_matrix
        1 - mean(as.vector(matrix[lower.tri(matrix)]))
    })
    delta.area.dec <- c(delta.area[1], diff(delta.area))

    delta <- delta.area.dec %>% 
        data.frame(dK = ., x = 1:length(.)) %>%
        ggplot() +
        geom_line(aes(x = x, y = dK), linewidth = 0.6) +
        geom_point(aes(x = x, y = dK), size = 3.5, color = "white", fill = "white") +
        geom_point(aes(x = x, y = dK), size = 2.5) +
        theme_minimal() +
        theme(
            aspect.ratio = 1,
            panel.border = element_rect(color = "black", fill = "transparent", linewidth = 1),
            panel.grid.minor.x = element_blank()
        ) +
        labs(title = "Delta area", x = "k", y = "Δ(K)")

    res <- list(matrix.ls = consensus.matrix.ls, cdf = cdf, delta = delta)
    return(res)
}

#' @importFrom purrr walk
#' @importFrom ggplot2 ggsave
create_consensus_clustering_report <- function(
    consensus.matrix.ls,
    result_path = "docs/result/maldi/consensus_clustering"
) {
    suppressMessages({suppressWarnings({
        # Save pngs
        dir.create(result_path)
        cdf_path <- paste0(result_path, "/consensus_cdf.png")
        ggsave(cdf_path, consensus.matrix.ls$cdf, width = 6, height = 6)
        delta_path <- paste0(result_path, "/consensus_delta_area.png")
        ggsave(delta_path, consensus.matrix.ls$delta, width = 6, height = 6)
        pb.ggsave <- cli::cli_progress_bar(total = length(consensus.matrix.ls$matrix.ls), format = paste0("Saving consensus clustering result in image file ... ", " {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}"), clear = T)
        walk(consensus.matrix.ls$matrix.ls, \(obj){
            ggsave(paste0(result_path, "/dendro_heatmap_k_", obj[["k"]], ".png"), obj[["heatmap"]], width = 15, height = 10)
            cli::cli_progress_update(id = pb.ggsave)
        })

        cli::cli_process_done(id = pb.ggsave)
    })})

    cat("## Consensus CDF\n\n")
    cat("![](/", cdf_path ,")\n\\pagebreak\n\n", sep = "")
    cat("## Delta Area Curve of consensus clustering\n\n")
    cat("![](/", delta_path ,")\n\n\\pagebreak\n\n", sep = "")

    pb.report <- cli::cli_progress_bar(total = length(consensus.matrix.ls$matrix.ls), format = paste0("Writing consensus clustering report ... ", " {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}"), clear = T)
    walk(names(consensus.matrix.ls$matrix.ls), \(k){
        cat("## Consensus clustering k = ", k, "\n\n")
        cat("![](/", paste0(result_path, "/dendro_heatmap_k_", k, ".png"), ")\n\n\\pagebreak", sep = "")
        cli::cli_progress_update(id = pb.report)
    })
    cli::cli_process_done(id = pb.report)
}