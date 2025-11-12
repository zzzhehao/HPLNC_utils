#' Wrapper for Quarto Render Species Summary
#' 
#' @importFrom lubridate today
render_species_summary <- function() {
    
    # Define your paths and filenames
    input_file <- "docs/form_template/species_summary_template.qmd"
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