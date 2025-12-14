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
        " -t pdf",
        " --output-dir ", output_directory,  # <-- Use --output-dir for the path
        " --output ", output_filename,      # <-- Use --output for the filename
        " --no-clean"
    )
    system(command)
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
            path <- paste0("docs/result/delim/morpho/", sp, ".png")
            ggsave(path, combplot, width = 10, height = 5)
        })
    })})

    cat("\n \\pagebreak \n\n")

    # generate pdf content
    walk(summaryBysp, \(summary){
        sp <- summary[["sp"]]
        cat("## ", summary[["sp"]], "\n\n")
        cat("![](", paste0("/docs/result/delim/morpho/", sp, ".png"), ")")
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
