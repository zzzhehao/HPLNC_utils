format_metadata.Specimen.Haploniscidae <- function(tbl, formatting = T) {
    if (!formatting) { # way to escape
        return(tbl)
    }

    print("Formatting...")

    tbl %>% mutate(
        # formatting only
        station = factor(station),
        gear = factor(gear),
        taxa1 = factor(taxa1),
        taxa2 = factor(taxa2),
        taxa3 = factor(taxa3),
        genus_morphology = factor(genus_morphology),
        species_morphology = factor(species_morphology),
        sex = factor(sex),
        cf = factor(cf),
        fixedFirst = factor(fixedFirst),
        fixedNow = factor(fixedNow),
        storage = factor(storage),
        cruise = factor(cruise),
        sampleUnit = factor(sampleUnit),
        sampleSpec = factor(sampleSpec),

        # calculation
        gensp_morphology = factor(paste(genus_morphology, species_morphology)),
        voucher_valid = ifelse(!is.na(voucher), paste0(voucher_prefix, str_pad(voucher, 3, "left", "0")), NA)
    ) %>%
        return()
} 

format_metadata.Sequence.NCBI <- function(tbl, formatting = T, format = 1) {
    if (!formatting) { # way to escape
        return(tbl)
    }
    
    print("Formatting...")
    # Relevant sequence only
    if (format == 1) {
        tbl.format <- tbl %>% 
            dplyr::filter(c_filter_1 & c_filter_2)
    }
    # Custom columns only
    if (format == 2) {
        tbl.format <- tbl %>% 
            dplyr::select(starts_with("c_"), -contains("filter"))
    }

    return(tbl.format)
}

format_metadata.Station <- function(tbl, formatting = T) {
    tbl.format <- tbl %>%
        mutate(
            expedition = factor(expedition), 
            gear = factor(gear),
            depthZone = factor(depthZone),
            region = factor(region),
            area = factor(area)
        )
    return(tbl.format)
}