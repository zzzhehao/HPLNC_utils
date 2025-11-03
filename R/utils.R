#' Generate path to given date's asset folder in log (is this really necessary?)
date2assets <- function(date) {
    paste0("docs/logs/assets/", str_pad(as.character(date), 4, "left", "0"))
}

#' Database query of Haploniscid specimen using DZMB-2-HH ID
ID.inquery <- function(ID) {
    DBpullTable("metadata.Specimen.Haploniscidae") %>% filter(DZMB2HH %in% ID) %>% t()
}

#' Clear Temporary File in R (internal)
clearRCache <- function() {
    unlink(tempdir(), recursive = TRUE)
    dir.create(tempdir())
}
