#' Format Specimen Metadata
#' 
#' @param format numerical. Formatting options. 0 for not including morphocheck result. Default to 1.
#' 
#' @import dplyr
format_metadata.Specimen.Haploniscidae <- function(tbl, formatting = T, format = 1) {
    if (!formatting) { # way to escape
        return(tbl)
    }

    print("Formatting...")

    tbl.fmt <- tbl %>% mutate(
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
    )
    if (format > 0) {
        morphometa <- load_morphocheck()
        tbl.fmt <- left_join(tbl.fmt, morphometa %>% dplyr::select(c("DZMB2HH", "sex_ZH", "stage_ZH", "n", "gensp_morpho_ZH")), by = "DZMB2HH")
    }
    return(tbl.fmt)
} 

#' Format NCBI Sequence Metadata
#' @import dplyr
format_metadata.Sequence.NCBI <- function(tbl, formatting = T, format = 0) {
    if (!formatting) { # way to escape
        return(tbl)
    }
    
    print("Formatting...")
    # Relevant sequence only
    if (format == 0) {
        tbl.format <- tbl %>% 
            dplyr::filter(c_filter_1 & c_filter_2)
    }
    # Custom columns only
    if (format == 1) {
        tbl.format <- tbl %>% 
            dplyr::select(starts_with("c_"), -contains("filter"))
    }

    return(tbl.format)
}

#' Format Station Metadata
#' @import dplyr
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

#' Format morphocheck result
#' @import stringr
#' @import dplyr
format_morphocheck <- function(tbl) {
    tbl %>% mutate(
        ON_morpho_ZH = gsub("\\.", "", ON_morpho_ZH),
        sex_ZH = factor(sex_ZH),
        stage_ZH = factor(stage_ZH),
        n = ifelse(is.na(n), 1, n),
        gen_morpho_ZH = factor(gen_morpho_ZH),
        sp_morpho_ZH = factor(sp_morpho_ZH),
        DZMB2HH = as.numeric(DZMB2HH),
        gensp_morpho_ZH = case_when(
            is.na(gen_morpho_ZH) ~ "Haploniscidae sp.",
            !is.na(gen_morpho_ZH) & is.na(sp_morpho_ZH) ~ paste(gen_morpho_ZH, "sp."),
            str_detect(sp_morpho_ZH, "unicornis|aduncus") ~ "Haploniscus unicornis complex",
            .default = paste(gen_morpho_ZH, sp_morpho_ZH)
        )
    )
}