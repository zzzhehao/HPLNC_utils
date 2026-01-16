#' @include class.R
NULL

#' Consensus Clustering
#'
#' @description Perform consensus clustering (CC) with any base method.
#'
#' @name cc
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
#' @seealso [ccmk()]
#' @seealso [consensus_clustering_summarize()]
#' 
#' @return A S7 object of the class \code{cc}, with following properties:
#' \describe{
#'   \item{\code{k}}{Number of clusters. \code{k} used for cc.}
#'   \item{\code{reps}}{Number of iterations}
#'   \item{\code{subset}}{Subset proportion of the original matrix at each subsampling}
#'   \item{\code{dist_method}}{Method shortcut used for distance calculation}
#'   \item{\code{dist_call}}{Actual call for distance calculation}
#'   \item{\code{seed}}{Seeds used for each subsampling}
#'   \item{\code{base_method}}{Method shortcut used for base clustering}
#'   \item{\code{base_call}}{Actual call for base clustering}
#'   \item{\code{summary_call}}{Actual call for summarizing clustering}
#'   \item{\code{sample_name}}{Sample names extracted from the \code{source_matrix}}
#'   \item{\code{assignment}}{Clustering result, cluster assignment of each sample}
#'   \item{\code{silhouette_width_avg}}{Average silhouette width}
#'   \item{\code{source_matrix}}{Source matrix}
#'   \item{\code{cc_matrix}}{Generated consensus matrix}
#'   \item{\code{M}}{A matrix composed by number of times each sample pair clustered together}
#'   \item{\code{I}}{A matrix composed by number of times each samples are subsampled together}
#'   \item{\code{heatmap}}{A [ggplot2::ggplot] object, heatmap of the result}
#'   \item{\code{silhouette_plot}}{A [ggplot2::ggplot] object, silhouette width plot}
#'   \item{\code{summary_clustering}}{Final summarizing clustering result}
#' }
#' @export
cc <- new_class(
    "cc",
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
        source_matrix = class_any,
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
        verbose = T,
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
            seed <- round(runif(reps, 0, reps*1000), 0)
        } else {
            if (length(seed) < reps) {
                stop("Seed has to be equal to greater than reps.")
            }
        }

        # distance call
        if (is.null(dist_call)) {
            if (dist_method == "euclidean") {
                dist_call <- expr(dist(m))
            } else if (dist_method == "manhattan") {
                dist_call <- expr(dist(m, method = "manhattan"))
            } else if (dist_method == "cor_pearson") {
                dist_call <- expr(as.dist(1 - cor(t(m), method = "pearson")))
            } else if (dist_method == "cor_spearman") {
                dist_call <- expr(as.dist(1 - cor(t(m), method = "spearman")))
            } else if (!is.null(dist_call)) {
                cli::cli_alert_warning(
                    "{.val {dist_method}} is not a short cut method, using specified expression: {.val {dist_call}}"
                )
            } else if (is.null(dist_call)) {
                cli::cli_abort(
                    "{.val {dist_method}} is not a short cut method, please specify the expression for distance calculation."
                )
            }
        }

        # clustering call
        if (is.null(base_call)) {
            if (base_method == "hclust") {
                base_call <- expr(hclust(dist.matrix, method = "ward.D2"))
            } else if (base_method == "kmeans") {
                base_call <- expr(kmeans(m, centers = k))
                cli::cli_alert_warning(
                    "k-mean clustering using strictly euclidean distance, anything else will be ignored."
                )
            }
        }

        # core loop
        if (verbose) {
            pb.sub <- cli::cli_progress_bar(
                total = reps,
                format = paste0(
                    "Resampling and clustering {cli::pb_bar} | {cli::pb_percent} ETA: {cli::pb_eta}"
                ),
                clear = T
            )
        }
        matrices.perm <- reduce(seq(reps), .init = matrices, .f = \(matrices, r) {
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

            if (verbose) {
                cli::cli_progress_update(id = pb.sub)
            }

            return(matrices)
        })

        if (verbose) {
            cli::cli_progress_done(pb.sub)
        }

        consensus.matrix <- matrices.perm[["M"]] / matrices.perm[["I"]]
        consensus.matrix[is.na(consensus.matrix)] <- 0

        # Here merge with summary function
        summary_call <- expr(hclust(as.dist(1 - matrix), method = "ward.D2"))
        summary.ls <- consensus_clustering_summarize(
            consensus.matrix, k, summary_call, 
            plot_labs = list(labs(caption = paste(
                paste("distance call:", expr_text(dist_call)),
                paste("base call:", expr_text(base_call)),
                paste("summary call:", expr_text(summary_call)),
                sep = "\n"
            ))))

        # silhouette
        silhouette_dist_call <- do.call(
            substitute,
            list(dist_call, list(m = quote(source_matrix)))
        )
        sil <- cluster::silhouette(
            as.numeric(summary.ls$cl$cc_cluster),
            eval_tidy(silhouette_dist_call)
        ) %>%
            as.data.frame()
        sil_data <- sil %>%
            arrange(cluster, desc(sil_width)) %>%
            mutate(id = row_number())
        avg_width <- mean(sil_data$sil_width)

        # silhouette plot
        silhouette_plot <- ggplot(
            sil_data,
            aes(x = id, y = sil_width, fill = cluster)
        ) +
            geom_col(width = 1) +
            # scale_fill_brewer(palette = "Set1") +
            geom_hline(
                yintercept = avg_width,
                linetype = "dashed",
                color = "red"
            ) +
            labs(
                title = paste("Silhouette Plot for k =", k),
                subtitle = paste(
                    "Average Silhouette Width =",
                    round(avg_width, 2)
                ),
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

        obj <- new_object(
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

        class(obj) <- c("cc", class(obj))
        return(obj)
    }
)

#' Generate Heatmap with Dendrogram in ggplot2 system.
#' 
#' @param hc.obj A hierarchical clustering object produced by [stats::hclust()].
#' @param matrix The same matrix used for generating \code{hc.obj}, will be used for plotting the heatmap. The values must be normalized in range \[0, 1\] if \code{direction = 1}, or centralized at 0 with range \[-1, 1\] if \code{direction = 2}.
#' @param heat_name Character. Name of the parameter indicated by the heat, i.e. name of the values stored in the matrix.
#' @param title Character. Name of the plot.
#' @param font_size Numeric. Font size of the sample label.
#' @param dendro.stretch A numeric vector in length of 2. Controls vertical stretching of the dendrogram to align with heatmap. Will be passed to \code{expand} argument in [ggplot2::scale_y_continuous()].
#' @param direction Numeric. Whether heatmap shows correlation in one direction (the level of correlation) or two (both the direction and the level of correlation). 
#' @param palette A character vector of a color palette used for filling tiles in heatmap. 
#' @param zeroaswhite Logical. Whether the 0 will be replaced by white.
#' @param grid Logical. Whether to show the grid in heatmap.
#' @param plot_labs [ggplot2::labs()] wrapped in a list, will be passed directly to [ggplot2::ggplot()].
#' @param tip_length Numeric. Length of elongation of the dendrogram tips.
#' @param themes [ggplot2::theme()] wrapped in a list, will be passed directly to [ggplot2::ggplot()].
#' @param heat_lab A character vector in length of 2. Labels to show in the legend for the minimal and maximal value.
#' @param args Further ggplot functions to be added to the plot. 
#' 
#' @seealso [stats::hclust()]
#' @seealso [ggplot2::ggplot()]
#' 
#' @return A ggplot object.
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @import ggplot2
#' @export
gg_heatmap_hclust <- function(
    hc.obj,
    matrix,
    heat_name = "Heat value",
    title = "",
    font_size = 8,
    dendro.stretch = c(0, 0.6),
    direction = 1,
    palette = NULL,
    zeroaswhite = T,
    grid = T,
    plot_labs = list(labs(x = NULL, y = NULL)),
    tip_length = 0,
    themes = list(theme()),
    heat_lab = NULL,
    args = NULL,
    narrow = F
) {
    # check direction
    if (!direction %in% c(1, 2)) cli::cli_abort("{.var direction} must be either 1 or 2.")

    dend_order <- hc.obj$labels[hc.obj$order]
    heatmap_data <- matrix %>%
        as.data.frame() %>%
        tibble::rownames_to_column("row_id") %>%
        pivot_longer(-row_id, names_to = "col_id", values_to = "heat") %>%
        mutate(
            row_id = factor(row_id, levels = dend_order),
            col_id = factor(col_id, levels = dend_order)
        )

    if (!narrow) {
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
    } else {
        p.dendro <- hc.obj %>%
            as.dendrogram() %>%
            ggdendro::dendro_data() %>%
            ggdendro::segment() %>%
            ggplot() +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
            geom_segment(aes(y = 0, x = x, yend = tip_length, xend = x)) +
            scale_x_continuous(expand = dendro.stretch) +
            scale_y_continuous(expand = c(0, 0)) +
            theme_void() +
            theme(
                plot.margin = margin(l = 0, r = 0, t = 0, b = 0.1),
                aspect.ratio = 0.2
            )
    }

    p.heatmap <- heatmap_data %>%
        ggplot(aes(
            x = col_id,
            y = row_id,
            fill = heat,
            data_id = row_id,
            tooltip = row_id
        ))

    if (grid) {
        p.heatmap <- p.heatmap +
            geom_tile(
                color = "grey85", 
                linewidth = 0.1
            )
    }

    if (is.null(palette)) {
        palette <- rcartocolor::carto_pal(n = 15, name = "Sunset")
    }

    if (zeroaswhite) {
        if (direction == 1) {
            clr <- c("white", palette)
            palette.length <- length(palette) + 1
        } else if (direction == 2) {
            palette.length <- length(palette)
            if (palette.length %% 2 == 0) {
                clr <- c(palette[1:(palette.length/2)], "white", palette[(palette.length/2+1):palette.length])
            } else {
                clr <- palette
                clr[ceiling(palette.length/2)] <- "white"
            }
            palette.length <- length(palette) + 1
        }
    } else {
        clr <- palette
        palette.length <- length(palette)
    }

    if (is.null(heat_lab)) {
        heat_lab <- c("0", round(max(matrix), 0))
    }

    if (direction == 1) {
        p.heatmap <- p.heatmap +
            scale_fill_gradientn(
                name = heat_name,
                colors = clr,
                values = seq(0, 1, length.out = palette.length),
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
            )
    } else if (direction == 2) {
        p.heatmap <- p.heatmap +
            scale_fill_gradientn(
                name = heat_name,
                colors = clr,
                values = c(0, seq(-1, 1, length.out = palette.length-1)) %>% sort(),
                breaks = c(-1, 1),
                labels = heat_lab,
                guide = guide_colorbar(
                    title.position = "top",
                    title.hjust = 0.5,
                    label.position = "bottom",
                    ticks = F,
                    draw.ulim = F,
                    draw.llim = F
                ),
                limits = c(-1, 1),
            )
    }
    p.heatmap <- p.heatmap +
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
            axis.text.y = element_text(
                size = font_size,
                margin = margin(l = 5),
                hjust = 0
            ),
            plot.margin = margin(l = 0, r = 5, t = 5, b = 5)
        ) +
        themes +
        args

    if (narrow) {
        p.heatmap <- p.heatmap +
            theme(
                plot.margin = margin(l = 0, r = 5, t = 2.5, b = 5)
            ) +
            labs(title = NULL)
    }

    # merge
    if (narrow) {
        pp <- (p.dendro / p.heatmap) +
            patchwork::plot_layout(heights = c(1, 5))
    } else {
        pp <- p.dendro + p.heatmap + patchwork::plot_layout(widths = c(1, 5))
    }
    return(pp)
}

#' Summarize Consensus Clustering Result
#' 
#' WIP, plug in evaluation to generate supporting values.
#'
#' @param matrix Consensus matrix produced by [cc()].
#' @param k Numerical. Number of groups (k) to be passed to \code{summary_call} and [stats::cutree()].
#' @param summary_call Final clustering call to group samples according to \code{matrix}.
#' @param heat_name Character. Title of the legend for the tile color.
#' @param title Character or an expression. Title of the heatmap. Will be passed to [gg_heatmap_hclust()].
#' @param plot_labs [ggplot2::labs()] call wrapped in a list. Will be passed to [gg_heatmap_hclust()].
#' @param heat_lab A character vector in length of two. Labels for the min and max value of the color tile legend. Will be passed to [gg_heatmap_hclust()].
#' @param ... Further argument passed to [gg_heatmap_hclust()].
#' 
#' @seealso [gg_heatmap_hclust()]
#' @seealso [consensus_clustering()]
#' 
#' @return A list of summary result.
#' \describe{
#'   \item{\code{k}}{The number of clusters.}
#'   \item{\code{matrix}}{The consensus matrix input.}
#'   \item{\code{fc}}{Final clustering result, an \code{\link[stats:hclust]{hclust}} object.}
#'   \item{\code{cl}}{A data frame containing the cluster assignments (Clustering result).}
#'   \item{\code{heatmap}}{The consensus clustering heatmap, a \code{ggplot} object.}
#' }
#' 
#' @importFrom rlang eval_tidy
#' @export
consensus_clustering_summarize <- function(
    matrix,
    k,
    summary_call,
    heat_name = "Consensus index",
    title = expr(paste0("Consensus clustering (k = ", k, ")")),
    heat_lab = c("Least similar", "Most similar"),
    ...
) {
    final.c <- eval_tidy(summary_call)
    cl <- cutree(final.c, k) %>%
        as.factor() %>%
        data.frame(label = names(.), cc_cluster = .)
    heatmap <- gg_heatmap_hclust(
        final.c,
        matrix,
        heat_name = heat_name,
        title = ifelse(is.character(title), title, eval_tidy(title)),
        heat_lab = heat_lab,
        ...
    )

    res <- list(
        k = k,
        matrix = matrix,
        fc = final.c,
        cl = cl,
        heatmap = heatmap
    )
    return(res)
}

plot.cc <- function(x, y, ...) {
    x@heatmap
    invisible(NULL)
}

#' Perform Consensus Clustering for a Range of K
#' 
#' @param source_matrix Multivariate matrix to perform [cc()] on.
#' @param kmin Integer. Minimum \code{k} of the range.
#' @param kmax Integer. Maximum \code{k} of the range.
#' @param ... Further argument passed to [cc()].
#' 
#' @seealso [cc()]
#' 
#' @return A S7 object of the class \code{ccmk}. 
#'
#' @importFrom rlang eval_tidy
#' @importFrom purrr map
#' @importFrom purrr map_vec
#' @import ggplot2
#' @export
ccmk <- new_class(
    "ccmk",
    properties = list(
        kmin = class_integer,
        kmax = class_integer,
        length = class_integer,
        source_matrix = class_numeric,
        cdf = new_S3_class("ggplot"),
        delta = new_S3_class("ggplot"),
        silhouette_width = class_data.frame,
        best_local_k = new_property(getter = function(self) {
            max_width <- max(self@silhouette_width$silhouette_width)
            k <- self@silhouette_width[self@silhouette_width$silhouette_width == max_width,]$k
            return(list(k = k, width = max_width))
        }),
        silhouette_width_plot = new_S3_class("ggplot"),
        cc_obj = class_list,
        cc_meta = class_list
    ),
    validator = function(self) {
        if (!is.matrix(self@source_matrix)) {
            "@source_matrix must be a matrix."
        }
    },
    constructor = function(
        source_matrix,
        kmin,
        kmax,
        ...
    ) {
        # core loop for CC
        pb.cc <- cli::cli_progress_bar(
            total = length(kmin:kmax),
            format = paste0(
                "Running consensus clustering k = {cli::pb_current+kmin-1} {cli::pb_bar} | ETA: {cli::pb_eta}"
            ),
            clear = T
        )
        consensus.matrix.ls <- map(kmin:kmax, \(k) {
            cc_obj <- cc(source_matrix, k, verbose = F, ...)
            cli::cli_progress_update(id = pb.cc)
            return(cc_obj)
        })
        cli::cli_progress_done(id = pb.cc)
        names(consensus.matrix.ls) <- kmin:kmax

        # consensus CDF
        cdf.palette <- rcartocolor::carto_pal(
            n = length(kmin:kmax),
            name = "Fall"
        )

        cdf_fun.ls <- map(consensus.matrix.ls, \(m) {
            matrix <- m@cc_matrix
            i <- match(m@k, kmin:kmax)
            stat_function(
                fun = ecdf(matrix[lower.tri(matrix)]),
                size = 1,
                color = cdf.palette[[i]]
            )
        })

        cdf <- ggplot() +
            cdf_fun.ls +
            theme_monochrome(aspect.ratio = 0.6) +
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
            geom_point(
                aes(x = x, y = dK),
                size = 3.5,
                color = "white",
                fill = "white"
            ) +
            geom_point(aes(x = x, y = dK), size = 2.5) +
            theme_minimal() +
            theme(
                aspect.ratio = 0.6,
                panel.border = element_rect(
                    color = "black",
                    fill = "transparent",
                    linewidth = 1
                ),
                panel.grid.minor.x = element_blank()
            ) +
            labs(title = "Delta area", x = "k", y = "Δ(K)")

        silhouette_width <- map_vec(consensus.matrix.ls, \(m) {
            m@silhouette_width_avg
        }) %>%
            data.frame(k = as.numeric(names(.)), silhouette_width = .)

        silhouette_width.plot <- silhouette_width %>%
            ggplot() +
            geom_line(aes(x = k, y = silhouette_width), linewidth = 0.6) +
            geom_point(aes(x = k, y = silhouette_width), color = "white", fill = "white", size = 3.5) +
            geom_point(aes(x = k, y = silhouette_width), size = 2.5) +
            scale_x_continuous(breaks = seq(0, kmax, by = 5)) +
            theme_monochrome(aspect.ratio = 0.6)

        meta <- list(
            reps = consensus.matrix.ls[[1]]@reps,
            subset = consensus.matrix.ls[[1]]@subset,
            dist_call = consensus.matrix.ls[[1]]@dist_call,
            base_call = consensus.matrix.ls[[1]]@base_call,
            summary_call = consensus.matrix.ls[[1]]@summary_call,
            seed = consensus.matrix.ls[[1]]@seed
        )

        new_object(
            S7_object(),
            kmin = as.integer(kmin),
            kmax = as.integer(kmax),
            length = as.integer(length(kmin:kmax)),
            source_matrix = source_matrix,
            cdf = cdf,
            delta = delta,
            silhouette_width = silhouette_width,
            silhouette_width_plot = silhouette_width.plot,
            cc_obj = consensus.matrix.ls,
            cc_meta = meta
        )
    }
)

#' Generate Markdown Report for \code{ccmk} Object
#' @name md_report-ccmk
#' @rdname md_report
#' 
#' @param report_path Path to export the generated figures.
#' 
#' @importFrom purrr walk
#' @importFrom ggplot2 ggsave
method(md_report, ccmk) <- function(x, report_path = read_config("report_folder"), result_path = read_config("result_folder"), report_name = "") {

    result_path <- paste(result_path, "maldi", "consensus_clustering", sep = "/")
    date_code <- lubridate::today() %>% format("%Y%m%d")
    report_name <- ifelse(report_name == "", date_code, paste0(date_code, "_", report_name))
    # Save pngs
    dir.create(result_path)
    cdf_path <- paste0(result_path, "/", report_name, "_consensus_cdf.svg")
    ggsave(cdf_path, x@cdf, width = 6, height = 4)
    delta_path <- paste0(result_path, "/", report_name, "_consensus_delta_area.svg")
    ggsave(delta_path, x@delta, width = 6, height = 4)
    sw_path <- paste0(result_path, "/", report_name, "_silhouette_width.svg")
    ggsave(sw_path, x@silhouette_width_plot, width = 6, height = 4)
    pb.ggsave <- cli::cli_progress_bar(
        total = length(x@cc_obj),
        format = paste0(
            "Saving consensus clustering result in image file ... ",
            " {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}"
        ),
        clear = T
    )
    walk(x@cc_obj, \(obj) {
        ggsave(
            paste0(
                result_path,
                "/", report_name, 
                "_dendro_heatmap_k_",
                obj@k,
                ".svg"
            ),
            obj@heatmap,
            width = 15,
            height = 10
        )
        cli::cli_progress_update(id = pb.ggsave)
    })

    cli::cli_process_done(id = pb.ggsave)

    sec_meta <- paste(
        "## Metadata",
        sep = "\n\n"
    )

    sec_meta.tbl <- knitr::kable(
        data.frame(
            Parameter = c("K min", "K max", "Resampling iteration", "Resampling subset", "distance call", "base method call", "summary call"),
            Value = c(x@kmin, x@kmax, x@cc_meta$reps, x@cc_meta$subset, expr_text(x@cc_meta$dist_call), expr_text(x@cc_meta$base_call), expr_text(x@cc_meta$summary_call))),
        format = "latex",
        align = c("l", "l")
    )

    sec_kinfo <- paste(
        "## Consensus CDF",
        paste0("![](/", cdf_path, ")"),
        "## Delta Area Curve of consensus clustering",
        paste0("![](/", delta_path, ")"),
        "## Average Silhouette Width Curve",
        paste0("![](/", sw_path, ")"),
        sep = "\n\n"
    )

    pb.report <- cli::cli_progress_bar(
        total = length(x@cc_obj),
        format = paste0(
            "Writing consensus clustering report ... ",
            " {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}"
        ),
        clear = T
    )
    sec_run_image <- paste0(map_vec(names(x@cc_obj), \(k) {
        runinfo <- paste0(
            "\n## Consensus clustering k = ", k, "\n\n",
            "![](/",
            paste0(result_path, "/", report_name, "_dendro_heatmap_k_", k, ".svg"),
            ")\n\n\\pagebreak\n\n"
        )
        cli::cli_progress_update(id = pb.report)
        return(runinfo)
    }))
    cli::cli_process_done(id = pb.report)

    return(list(sec_meta = sec_meta, sec_meta.tbl = sec_meta.tbl, sec_kinfo = sec_kinfo, sec_run_image = sec_run_image))
}

#' Ensemble Consensus Clustering
#'
#' @description Combines multiple consensus clustering objects (cc) into a single ensemble result. 
#' @name cc_ensemble
#' @return A S7 object of the class \code{cc_ensemble}, with following properties:
#' \describe{
#'   \item{\code{cc_list}}{A list of [cc] class object used for assembly}
#'   \item{\code{weights}}{A numeric vector representing weight of each cluster in assembly}
#'   \item{\code{weight.info}}{Additional information about weighting}
#'   \item{\code{k}}{Number of cluster (\code{k}) used for final assembly clustering}
#'   \item{\code{type}}{Type of assembly approach}
#'   \item{\code{method}}{Method of assembly}
#'   \item{\code{ensemble_call}}{Actual function call of the final clustering}
#'   \item{\code{ensemble_matrix}}{A consensus matrix produced by assembling [cc] in \code{cc_list}}
#'   \item{\code{source_matrix}}{The original multivariate matrix used to produce all [cc] object}
#'   \item{\code{ensemble_clustering}}{Final clustering result}
#'   \item{\code{assignment}}{Cluster assignment of each sample}
#'   \item{\code{heatmap}}{A [ggplot] object, the heatmap of the clustering result}
#' }
NULL

cc_ensemble <- new_class(
    "cc_ensemble", 
    properties = list(
        cc_list = class_list,
        weights = class_numeric,
        weight.info = class_character,
        k = class_integer,
        type = class_character,
        method = class_character,
        ensemble_call = class_call,
        ensemble_matrix = class_numeric,
        source_matrix = class_any,
        ensemble_clustering = class_any,
        assignment = class_data.frame,
        heatmap = new_S3_class("ggplot")
    ),
    validator = function(self) {
        if (!is.matrix(self@source_matrix)) {
            "@source_matrix must be a matrix."
        }
    }
)

#' @rdname cc_ensemble
#' @description
#' \code{cc_ensemble_cspa} uses a CSPA-like approach. 
#' 
#' @param cc_list A list of \code{cc} objects. Ideally, these are the "best k" results from different methods (e.g., list(pearson_k5, spearman_k6)).
#' @param weights Numeric vector. Weights for each method. If NULL, equal weights are used.
#' @param k Numeric. The number of clusters to cut the final ensemble into. If NULL, it suggests a k based on the ensemble structure (optional extension).
#' @param method Character. How to combine the matrices. "mean" (arithmetic mean) or "max" (takes the max value, useful if you want to link samples if *any* method finds them similar).
#' @param ensemble_call A function call for the final clustering of assembled consensus matrix. Take \code{matrix} as input.
#' @param ... Additional arguments passed to [consensus_clustering_summarize()], such as \code{heat_name} or \code{plot_labs}.
#' 
#' @seealso [cc()]
#' 
#' @export
cc_ensemble_cspa <- function(
    cc_list = list(),
    weights = NULL,
    k = NULL,
    method = "mean", 
    ensemble_call = expr(hclust(as.dist(1 - matrix), method = "ward.D2")),
    ...
){
    # check type and samples
    ref_sample_name <- cc_list[[1]]@sample_name
    check <- map_dfr(cc_list, \(c){
        class.check <- inherits(c, "HPLNC::cc")
        sample.check <- identical(sort(c@sample_name), sort(ref_sample_name))
        return(data.frame(class.check, sample.check))
    }) 

    if (any(!check$class.check)) {
        cli_abort("All objects in {.arg cc_list} must be a {.cls HPLNC::cc} objects.")
    }
    if (any(!check$sample.check)) {
        cli_abort("All objects must contain the same set of samples.")
    }

    # extract cc matrix
    cc_matrices <- map(cc_list, \(c){
        c@cc_matrix[ref_sample_name, ref_sample_name]
    })
    cc_length <- length(cc_matrices)
    cc_width <- map_vec(cc_list, \(c){
        c@silhouette_width_avg
    })
    cc_k <- map_vec(cc_list, ~.x@k)

    if (any(is.null(weights), weights == "equal")) {
        weights <- rep(1, cc_length)
        weight.info <- "Equal weight."
    } else if (length(weights) == 1 & is.character(weights)) {
        if (weights == "sil-1") { # linear
            weight.info <- "Automatic weight calculation based on silhouette width: linear."
            sum.width <- sum(cc_width)
            weights <- cc_width/sum.width
        } else if (weights == "sil-2") { # square
            weight.info <- "Automatic weight calculation based on silhouette width: square power."
            sum.width <- sum(cc_width^2)
            weights <- cc_width^2/sum.width
        } else if (weights == "sil-3") { # cubic
            weight.info <- "Automatic weight calculation based on silhouette width: cubic power."
            sum.width <- sum(cc_width^3)
            weights <- cc_width^3/sum.width
        }
    } else if (length(weights) != cc_length) {
        cli_abort("Length of {.arg weight} must be the same of {.arg cc_list}.")
    } else {
        weight.info <- "Specified."
    }

    if (any(is.null(k), k == "auto")) {
        k <- sum(cc_k*weights) %>% floor()
    }

    if (method == "mean") {
        weights <- weights/sum(weights)
        ensemble <- reduce2(cc_matrices, weights, \(matrix.ac, matrix.n, weight){
            return(matrix.ac + matrix.n * weight)
        }, .init = matrix(rep(0, length(cc_matrices[[1]])), nrow = nrow(cc_matrices[[1]])))
    } else if (method == "max") {
        # Max Logic: Take the highest probability observed across methods
        # (Useful for "Fuzzy Union" - if any method says they are similar, they are)
        if (!is.null(weights)) {
            cli::cli_alert_warning("{.var method} is set to max, weights of individual cc matrix are not considered.")
            weight.info <- "Not considered because method = 'max'"
        }
        ensemble <- purrr::reduce(cc_matrices, \(m1, m2) pmax(m1, m2))
         
    } else {
        cli::cli_abort("Unknown method: {.var {method}}")
    }
    
    # Fix dimnames lost during reduce
    dimnames(ensemble) <- list(ref_sample_name, ref_sample_name)

    summary <- consensus_clustering_summarize(
        matrix = ensemble,
        k = k,
        summary_call = ensemble_call,
        ...
    )

    obj <- cc_ensemble(
        cc_list = cc_list,
        weights = weights,
        weight.info = weight.info,
        k = as.integer(k),
        type = "cspa",
        method = method,
        ensemble_call = ensemble_call,
        ensemble_matrix = ensemble,
        source_matrix = cc_list[[1]]@source_matrix,
        ensemble_clustering = summary$fc,
        assignment = summary$cl,
        heatmap = summary$heatmap
    )

    class(obj) <- c("cc_ensemble", class(obj))
    return(obj)
}

#' Generate Markdown Report for \code{cc_ensemble} Object
#' @rdname md_report
#' @name md_report-cc_ensemble
#' @param report_name A run name to put on the report and file name.
method(md_report, cc_ensemble) <- function(
    x, 
    result_path = read_config("result_folder"), 
    report_name = ""
) {
    result_path <- paste(result_path, "maldi", "consensus_clustering", sep = "/")
    date_code <- lubridate::today() %>% format("%Y%m%d")
    report_name <- ifelse(report_name == "", date_code, paste0(date_code, "_", report_name))

    # Save heatmap svg
    dir.create(result_path)
    heatmap_path <- paste0(result_path, "/", report_name, "_ensemble_", x@type, ".svg")
    ggsave(heatmap_path, x@heatmap, width = 15, height = 10)

    sec_heatmap <- paste0("![](/", heatmap_path, ")")


    # General metadata
    sec_meta.tblgen <- data.frame(
        Parameter = c("Number of CCs", "Ensemble type", "Method", "Summary k", "Ensemble call", "Weighting info"),
        Value = c(length(x@cc_list), x@type, x@method, x@k, expr_text(x@ensemble_call), x@weight.info)) %>%
        knitr::kable(
        format = "latex",
        align = c("l", "l")
    ) %>%
        kableExtra::kable_styling(font_size = 10, latex_options = "HOLD_position")

    # Run meta
    run.meta1 <- map_dfr(x@cc_list, \(c){
        data.frame(
            k = c@k,
            reps = c@reps,
            subset = c@subset,
            `distance call` = expr_text(c@dist_call),
            `base call` = expr_text(c@base_call)
        )
    }) %>%
        mutate(cc = 1:nrow(.), .before = everything())

    run.meta2 <- map_dfr(x@cc_list, \(c){
        data.frame(
            `summary_call` = expr_text(c@summary_call),
            `silhouette width` = c@silhouette_width_avg
        )
    }) %>%
        bind_cols(data.frame(weight = round(x@weights, 2))) %>%
        mutate(cc = 1:nrow(.), .before = everything())

    sec_meta.tbl1 <- knitr::kable(run.meta1, "latex", align = rep("l", ncol(run.meta1))) %>%
        kableExtra::kable_styling(font_size = 10, latex_options = "HOLD_position")

    sec_meta.tbl2 <- knitr::kable(run.meta2, "latex", align = rep("l", ncol(run.meta2))) %>%
        kableExtra::kable_styling(font_size = 10, latex_options = "HOLD_position")

    # Assignment table
    sec_res <- knitr::kable(
        x@assignment, 
        "latex", 
        align = c("l", "l"), 
        col.names = c("Label", "Cluster"), 
        row.names = FALSE,
        longtable = T
    ) %>% kableExtra::kable_styling(font_size = 10, latex_options = c("repeat_header"))

    return(list(sec_meta.tblgen = sec_meta.tblgen, sec_meta.tbl1 = sec_meta.tbl1, sec_meta.tbl2 = sec_meta.tbl2, sec_res = sec_res, sec_heatmap = sec_heatmap))
}

# S3 method
plot.cc_ensemble <- function(x, y, ...) {
    x@heatmap
    invisible(NULL)
}

#### Confidence calculation ####

#' Evaluation of a Consensus Clustering or their Ensemble Run
#' 
#' 
cc_eval <- new_class(
    "cc_eval",
    properties = list(
        sample_stats = class_data.frame,
        cluster_stats = class_data.frame,
        ic_matrix = class_numeric
    ),
    constructor = function(mat, assignment) {
        
        clusters <- assignment$cc_cluster
        diag(mat) <- 1 # self = M(i, i) = 1 
        
        # check alignment
        if(!identical(rownames(mat), assignment$label)) {
            assignment <- assignment[match(rownames(mat), assignment$label), ]
        }
        
        unique_clusters <- sort(unique(clusters))
        n_samples <- nrow(mat)
        
        # item consensus
        suppressMessages({
            ic_matrix <- map(unique_clusters, function(k) {
                in_cluster <- clusters == k
                idx <- which(in_cluster)
                
                col_subset <- mat[, idx, drop = FALSE]
                row_sums <- rowSums(col_subset) 
                
                # denominator Monti et al. 2003 equation 4
                denoms <- rep(length(idx), n_samples)
                denoms[in_cluster] <- length(idx) - 1
                
                # numerator Monti et al. 2003 equation 4, M(i, j) where i is not equal to j
                nums <- row_sums
                nums[in_cluster] <- nums[in_cluster] - 1 # exclude self (=1)
                
                # singleton 
                res <- nums / denoms
                res[denoms == 0] <- 0
                return(res)
            }) %>% 
                bind_cols() %>%
                as.matrix()
        })
        
        colnames(ic_matrix) <- paste0("ic_cl_", unique_clusters)
        rownames(ic_matrix) <- rownames(mat)
        
        # iterate through row indices
        stats_df <- map_dfr(seq_len(n_samples), function(i) {
            
            # Current row context
            current_cl <- clusters[i]
            row_vals <- ic_matrix[i, ]
            
            # Construct the column name for the self-cluster (e.g., "ic_cl_3")
            # This is safer than using numeric indices
            self_col_name <- paste0("ic_cl_", current_cl)
            
            # Get self IC
            ic_self <- row_vals[[self_col_name]]
            
            # Get others: subset vector by excluding the self name
            other_vals <- row_vals[names(row_vals) != self_col_name]
            
            if (length(other_vals) > 0) {
                # Find max among others
                ic_2nd <- max(other_vals)
                
                # Extract cluster ID from the name of the max value
                best_other_name <- names(other_vals)[which.max(other_vals)]
                assignment_2nd <- str_extract(best_other_name, "[0-9]+$")
            } else {
                # Handle k=1 case
                ic_2nd <- 0
                assignment_2nd <- NA_character_
            }
            
            list(
                ic_self = ic_self, 
                ic_2nd = ic_2nd, 
                assignment_2nd = assignment_2nd
            )
        })

        # compile statistics
        sample_stats <- tibble(
            label = rownames(mat),
            assignment = clusters
        ) %>%
            bind_cols(stats_df) %>% 
            bind_cols(as_tibble(ic_matrix)) %>%
            mutate(
                # Confidence score calculation
                confidence_score = (ic_self - ic_2nd) / if_else(ic_self == 0, 1, ic_self)
            ) %>%
            dplyr::select(label, assignment, ic_self, confidence_score, assignment_2nd, ic_2nd, everything())

        
        cluster_stats <- sample_stats %>%
            group_by(assignment) %>%
            summarise(
                cluster_consensus = mean(ic_self), # Monti et al. 2003 equation 3
                n_samples = n(),
                mean_confidence_score = mean(confidence_score)
            ) %>% 
            rename(cluster = assignment)
        
        obj <- new_object(
            S7_object(),
            sample_stats = sample_stats,
            cluster_stats = cluster_stats,
            ic_matrix = ic_matrix
        )

        class(obj) <- c("cc_eval", class(obj))
        return(obj)
    }
)

# Evaluation method to create cc_eval

#' @usage evaluate(cc)
#' @usage evaluate(cc_ensemble)
#' 
#' @details For class [cc] and [cc_ensemble], the function will calculate confidence statistics for samples and clusters.
#' 
#' @name evaluate
#' @export
method(evaluate, cc_ensemble) <- function(x) {
    mat <- x@ensemble_matrix
    assignment <- x@assignment
    
    eval.res <- cc_eval(mat, assignment)
    return(eval.res)
}

#' @export
method(evaluate, cc) <- function(x) {
    mat <- x@cc_matrix
    assignment <- x@assignment

    eval.res <- cc_eval(mat, assignment)
    return(eval.res)
}

#' Print Quality Evaluation Result of \code{cc} or \code{cc_ensemble} Object.
#' 
#' @rdname report
#' @name report-cc_eval
#' @usage report(cc_eval)
#' 
#' @export
method(report, cc_eval) <- function(x) {
    cat("\n-- Cluster stats --\n\n")
    x@cluster_stats %>%
        as_tibble() %>%
        mutate(
            n = round(n_samples, 0),
            .keep = "unused"
        ) %>%
        insight::format_table() %>%
        print()

    cat("\n\n-- Singleton --\n\n")
    if (any(x@cluster_stats$n_samples == 1)) {
        sgt <- paste0(x@cluster_stats$cluster[x@cluster_stats$n_samples == 1], collapse = ", ")
        cli::cli_alert_warning("{sgt} are singleton clusters.")
    } else {
        cat("No singleton cluster found, all good.")
    }

        cat("\n\n-- Sample stats --\n\n")
        
    # sample stats
    x@sample_stats %>% 
        as_tibble() %>%
        mutate(
            label = case_when(
                nchar(label) > 20 ~ str_sub(label, 1L, 15L) %>% paste("..."),
                .default = label
            ),
        ) %>% 
        dplyr::select(label, assignment, ic_self, confidence_score, assignment_2nd, ic_2nd) %>%
        insight::format_table() %>%
        print()
    
}

method(md_report, cc_eval) <- function(x) {
    sec_cluster_stats <- x@cluster_stats %>%
        as_tibble() %>%
        mutate(
            n = round(n_samples, 0),
            .keep = "unused"
        ) %>%
        knitr::kable(., "latex", align = rep("l", ncol(.))) %>%
        kableExtra::kable_styling(font_size = 10, latex_options = "HOLD_position")

    if (any(x@cluster_stats$n_samples == 1)) {
        sgt <- paste0(x@cluster_stats$cluster[x@cluster_stats$n_samples == 1], collapse = ", ")
        sec_singleton <- paste(sgt, "are singleton clusters.")
    } else {
        sec_singleton <- "No singleton cluster found, all good." 
    }

    sec_sample_stats <- map(x@cluster_stats$cluster, \(c){
        tbl <- x@sample_stats %>% 
        filter(assignment == c) %>%
        as_tibble() %>%
        mutate(
            label = case_when(
                nchar(label) > 50 ~ str_sub(label, 1L, 50L) %>% paste("..."),
                .default = label
            ),
            ic_self = round(ic_self, 2),
            confidence_score = round(confidence_score, 2),
            ic_2nd = round(ic_2nd, 2)
        ) %>% 
        dplyr::select(label, assignment, ic_self, confidence_score, assignment_2nd, ic_2nd) %>%
        knitr::kable(., "latex", align = rep("l", ncol(.)), longtable = T, booktabs = T) %>%
        kableExtra::kable_styling(font_size = 8, latex_options = c("HOLD_position", "repeat_header", "stripped"))

        # table title
        ttl <- paste0("Cluster ", c)
        return(list(tbl = tbl, ttl = ttl))
    })

    rn <- rownames(x@ic_matrix)
    sec_ic <- x@ic_matrix %>%
        round(5) %>%
        as_tibble() %>%
        mutate(Label = case_when(
                nchar(rn) > 30 ~ str_sub(rn, 1L, 30L) %>% paste("..."),
                .default = rn
            )) %>%
        dplyr::select(Label, everything()) %>%
        knitr::kable(
            "latex",
            longtable = T,
            booktabs = T
        ) %>%
        kableExtra::kable_styling(
            font_size = 6,
            latex_options = c("HOLD_position", "repeat_header", "stripped")
        )

    return(list(
        sec_cluster_stats = sec_cluster_stats,
        sec_singleton = sec_singleton,
        sec_sample_stats = sec_sample_stats,
        sec_ic = sec_ic
    ))
}