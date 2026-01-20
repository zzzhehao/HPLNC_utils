#' Wrapper for Quarto Render Species Summary
#' 
#' @importFrom lubridate today
render_species_summary <- function() {
    
    # Define your paths and filenames
    input_file <- "docs/template/species_summary_template.qmd"
    output_directory <- "docs/result/delim/"
    output_filename <- paste0("Species_Summary_", format(today(), "%Y%m%d"), ".pdf")

    # Construct the correct Quarto command
    command <- paste0(
        "quarto render ", input_file,
        " -t pdf --profile dev",
        " --output-dir ", output_directory,  # <-- Use --output-dir for the path
        " --output ", output_filename,      # <-- Use --output for the filename
        " --no-clean"
    )
    system(command)
}

#' Generate Species Summary
#' @import ggtext
#' @import ggplot2
#' @import dplyr
summary_sp <- function(df, metadata.station, metadata.specimen, species) {
    palette <- c(paletteer::paletteer_d("nationalparkcolors::Badlands")[3], "#d5d5d5")
    names(palette) <- c("TRUE", "FALSE")

    fill.palette <- c(paletteer::paletteer_d("nationalparkcolors::Badlands")[3], "#dfdfdf00")
    names(fill.palette) <- c("TRUE", "FALSE")

    barplot <- df %>%
        mutate(color = ifelse(gensp_morpho_ZH == species, "TRUE", "FALSE")) %>%
        ggplot() +
            geom_col( # Count bar
                aes(x = station, y = n, fill = color),
                width = 0.8,
                show.legend = FALSE
            ) +
            labs(x = "", y = "", fill = "Morphotypes", title = species, subtitle = paste0("Using morphocheck result from ", as.character(today()))) +
            scale_fill_manual(values = palette) +
            ggnewscale::new_scale_fill() +
            # Site symbol
            geom_point( 
                aes(x = station, y = -4, color = area_abbr_ZH, fill = area_abbr_ZH, shape = area_abbr_ZH), 
                size = 5, 
                show.legend = F,
            ) +
            scale_fill_manual(values = paletteer::paletteer_d("calecopal::eschscholzia")) +
            scale_color_manual(values = paletteer::paletteer_d("calecopal::eschscholzia")) +
            scale_shape_manual(values = c(21:24)) +
            # bar bottom loc_name
            geom_text( 
                aes(x = station, y = -4, label = station), color = "black", size = 3, show.legend = F
            ) +
            # bar basis lab
            geom_text(aes(x = station, y = -1.5, label = area_abbr_ZH), size = 3) + # Omit
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.07)), breaks = c((1:5)*4)) +
            theme(
                plot.title = element_markdown(),
                legend.text = element_markdown(size = 7.5),
                legend.position = "right",
                legend.key.size = unit(0.4, "cm"),
                legend.justification = c(1, 0),
                aspect.ratio = 0.7,
                plot.background = element_rect(fill = "white"),
                panel.background = element_rect(fill = "white"),
                panel.grid.major.x = element_line(color = "#CCC"),
                axis.ticks = element_blank(),
                legend.background = element_rect(fill = "white", color = "black"),
                axis.text.y.left = element_blank()
            ) +
            coord_flip()

    station.HL <- df %>% 
        filter(gensp_morpho_ZH == species) %>% 
        pull(station)

    station.occur <- metadata.station %>% 
        mutate(
            occur = case_when(
                station %in% station.HL ~ "TRUE", 
                .default = "FALSE"
            )
        )

    map <- plot_basemap_bathy(resolution = 10) +
        ggnewscale::new_scale_fill() +
        ggnewscale::new_scale_color() +
        geom_point(
            data = station.occur,
            aes(x = longStartDec, y = latStartDec, color = occur, fill = occur),
            size = 3,
            shape = 21
        ) + 
        scale_color_manual(values = palette) +
        scale_fill_manual(values = fill.palette)

    display.clean <- function(col) {
        ifelse(is.na(col), "", as.character(col))
    }

    specimen <- metadata.specimen %>%
        dplyr::filter(gensp_morpho_ZH == species) %>% 
        dplyr::select(c("DZMB2HH", "voucher", "station", "depthStart", "area_abbr_ZH", "n", "sex_ZH", "gear", "remark_morpho_ZH", "plan", "ON_morpho_ZH")) %>%
        mutate(
            voucher = display.clean(voucher),
            sex_ZH = display.clean(sex_ZH),
            remark_morpho_ZH = display.clean(remark_morpho_ZH),
            ON_morpho_ZH = display.clean(ON_morpho_ZH)
        ) %>%
        rename(
            sex = sex_ZH,
            area = area_abbr_ZH,
            remark = remark_morpho_ZH,
            status = plan,
            depth = depthStart,
            ON = ON_morpho_ZH
        ) %>%
        arrange(area, station, n, voucher)

    return(list(sp = species, barplot = barplot, map = map, specimen = specimen))
}

#' Create Species Summary Report in Image and PDF
#' 
#' @import kableExtra
#' @import ggpubr
#' @import purrr
#' @import dplyr
create_summary_sp <- function() {
    suppressMessages({suppressWarnings({
        metadata <- db_pull("metadata.Specimen.Haploniscidae")
        metadata.ebs <- db_pull("metadata.Station")
        # morphocheck <- load_morphocheck()
    
        metadata.merge <- left_join(metadata, metadata.ebs %>% dplyr::select(c("station", "latStartDec", "longStartDec", "depthZone", "area_abbr_ZH")), by = "station") 
            # left_join(morphocheck %>% dplyr::select(-c("voucher")), by = "DZMB2HH")
    
        # count species per station
        st.sp.summary <- metadata.merge %>%
            group_by(area_abbr_ZH, station, gensp_morpho_ZH) %>%
            summarize(n = sum(n))
    
        # format table so that stations from the same region stay together in plot
        station.level.desc <- st.sp.summary$station %>% unique()
        species.label.asc <- as.character(st.sp.summary$gensp_morpho_ZH) %>% unique() %>% sort()
        st.sp.summary <- st.sp.summary %>%
            mutate(
                station = factor(station, levels = station.level.desc),
                gensp_morpho_ZH = factor(gensp_morpho_ZH, levels = species.label.asc)
            )
    
        # generate summary
        summaryBysp <- map(st.sp.summary %>% .$gensp_morpho_ZH %>% levels(), \(sp){
            summary <- summary_sp(st.sp.summary, metadata.ebs, metadata.merge, sp)
        })
    
        # save barplot and map in picture
        walk(summaryBysp, \(summary){
            # extract information
            combplot <- ggarrange(summary[["barplot"]], summary[["map"]], ncol = 2, nrow = 1)
            sp <- summary[["sp"]]
    
            # save barplot and map as image
            path <- paste0("docs/result/delim/morpho/", sp, ".pdf")
            ggsave(path, combplot, width = 10, height = 5)
        })
    })})

    cat("\n \\pagebreak \n\n")

    # generate pdf content
    walk(summaryBysp, \(summary){
        sp <- summary[["sp"]]
        cat("## ", summary[["sp"]], "\n\n")
        cat("![](", paste0("/docs/result/delim/morpho/", sp, ".pdf"), ")")
        print(
            kable(summary[["specimen"]]) %>%
                kable_styling(font_size = 9, latex_options = "striped") %>%
                column_spec(9, width = "10em") %>%
                column_spec(10, width = "6em")
        )
        cat("\n \\pagebreak \n\n")
    })
}

#' Wrapper for Quarto Render Consensue Clustering Result Report
#' 
#' @importFrom lubridate today
#' @importFrom quarto quarto_render
render_cc_report <- function() {
    
    # Define paths
    input_file <- "docs/template/consensus_clustering_template.qmd"
    output_directory <- "docs/result/maldi/"
    output_filename <- paste0("Consensus_Clustering_result_", format(today(), "%Y%m%d"), ".pdf")

    # Construct the correct Quarto command
    command <- paste0(
        "quarto render ", input_file,
        " -t pdf",
        " --output-dir ", output_directory,  # <-- Use --output-dir for the path
        " --output ", output_filename,      # <-- Use --output for the filename
        " --no-clean --profile dev"
    )
    system(command)
}


render_cc_ensemble_report <- function(
    yaml_path, name, source_path, cc_expr, ...,
    input_file = "docs/template/cc_ensemble_report.qmd",
    output_directory = "docs/result/maldi/"
) {
    report_name <- name
    config_yaml <- write_params_yaml(yaml_path, name, source_path = source_path, cc_expr = cc_expr, ...)

    if (!is.null(report_name)) {report_name <- paste0("_", report_name, "_")}
    output_filename <- paste0("Ensemble_Consensus_Clustering_report", report_name, format(today(), "%Y%m%d"), ".pdf")
    command <- paste0(
        "quarto render ", input_file,
        " -t pdf",
        " --output-dir ", output_directory,
        " --output ", output_filename,
        " --execute-params ", config_yaml,
        " --no-clean --profile dev"
    )
    system(command)
}

render_dist_geoxgen <- function() {
    # Define your paths and filenames
    input_file <- "docs/template/dist-geoxgen_template.qmd"
    output_directory <- "docs/result/delim/"
    output_filename <- paste0("dist_geoxgen_report_", format(today(), "%Y%m%d"), ".pdf")

    # Construct the correct Quarto command
    command <- paste0(
        "quarto render ", input_file,
        " -t pdf --profile dev",
        " --output-dir ", output_directory,  # <-- Use --output-dir for the path
        " --output ", output_filename,      # <-- Use --output for the filename
        " --no-clean"
    )
    system(command)
}

#' Create Geo x Genentic Distance Analysis Report
#' @importFrom patchwork plot_layout
create_dist_geoxgen_report <- function() {
    dist.analysis <- readRDS("data/cache/dist_geoxgen.rds")
    walk(dist.analysis, \(g){
        complot <- (g$gengeo.cor.plot + g$geoxgen.plot) +
            patchwork::plot_layout(widths = c(1, 3)) &
            theme(
                legend.position = "bottom",
                plot.tag = element_text(face = "bold", size = 14)
            )

        output.file.cp <- paste0("docs/result/dist/", g$group_description, ".pdf")
        ggsave(output.file.cp, complot, width = 10, height = 5)
        
        output.file.crlg <- paste0("docs/result/dist/", g$group_description, "_correlogram.pdf")
        ggsave(output.file.crlg, g$correlogram.plot, width = 10, height = 5)
    })

    walk(dist.analysis, \(g){
        cat("## ", g$group_description, "\n\n")
        cat("![](", paste0("/docs/result/dist/", g$group_description, ".pdf"), ")\n\n")
        cat("\n**Metadata**\n\n")
        print(
            kable(g$metadata) %>% 
                kable_styling(font_size = 9, latex_options = "striped")
        )

        cat("\n**Geodesics distance**\n\n")
        print(
            kable(g$geo.dist %>% round(2)) 
        )

        cat("\n**p-distance**\n\n")
        print(
            kable(g$seq.dist.p %>% round(5)) 
        )

        cat("\n**Spatial correlation**\n\n")
        cat("Mantel test:\n\n")
        cat("- z.stat:", g$mantel$z.stat, "\n")
        cat("- p:", g$mantel$p, "\n\n")
        cat("![](", paste0("/docs/result/dist/", g$group_description, "_correlogram.pdf"), ")\n\n")

        cat("\n \\pagebreak \n\n")
    })
}
