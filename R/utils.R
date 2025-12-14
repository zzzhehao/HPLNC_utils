
#' Create assets folder for the day
#' @importFrom stringr str_pad
date2assets <- function(date) {
    paste0("docs/logs/assets/", str_pad(as.character(date), 4, "left", "0"))
}

#' Query specimen info with DZMB2HH ID
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
inquery_id <- function(ID) {
    db_pull("metadata.Specimen.Haploniscidae") %>% filter(DZMB2HH %in% ID) %>% t()
}

#' Query specimen info with voucher
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
inquery_voucher <- function(ID) {
    db_pull("metadata.Specimen.Haploniscidae") %>% filter(voucher %in% ID) %>% t()
}

#' Clear temporary cache file of R
clear_r_cache <- function() {
    unlink(tempdir(), recursive = TRUE)
    dir.create(tempdir())
}

#' Initialize daily log
#' @importFrom lubridate today
create_log <- function() {
    date <- format(today(), "%m%d")
    dir.create(date2assets(date))
    qmd.path <- paste0("docs/logs/", date, ".qmd")
    if (qmd.path %in% list.files("docs/logs", "\\.qmd$", full.names = T)) {stop("Log file already exists.")}
    template <- paste0(
'---
title: "', date, ' Log"
date: last-modified
abstract: ""
execute: 
  freeze: true
---

```{r}
devtools::load_all()
library(tidyverse)
```

## ')
    write(template, paste0("docs/logs/", date, ".qmd"))
}

#' Load latest morphocheck result
#' 
#' @import readxl
load_morphocheck <- function() {
    morphocheck <- readxl::read_xlsx("docs/notes/assets/sandbox/ZHH_morphocheck.xlsx", range = readxl::cell_cols("A:R"))
    file_path <- paste0("data/metadata/archive/morphocheck_", as.character(format(today(), "%Y%m%d")), ".csv") # make archive
    write.table(morphocheck, file_path, sep = ";", row.names = F)
    morphorcheck <- read.table(file_path, sep = ";", header = T)
    morphocheck <- format_morphocheck(morphocheck)
    return(morphocheck)
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

#' Map a vector's values using a lookup table
#'
#' @param input_vector The vector of original values
#' @param lut_df The lookup data frame
#' @param key_col The unquoted column name in `lut_df` to match against.
#' @param value_col The unquoted column name in `lut_df` to get the new values from.
#'
#' @import dplyr
#' @import rlang
#' @return A new vector with the mapped values, in the same order as the input.
map_values <- function(input_vector, lut_df, key_col, value_col) {
    key_col_string <- rlang::as_name(rlang::enquo(key_col))

    # Create an unnamed tibble and set the name afterwards
    temp_df <- tibble(input_vector) %>%
        setNames(key_col_string)

    temp_df %>%
        left_join(lut_df, by = key_col_string) %>%
        pull({{ value_col }})
}

#' @param input_value The original value to look up in LUT.
#' @param lut_df The lookup data frame
#' @param key_col The unquoted column name in `lut_df` to match against.
#' @param value_col The unquoted column name in `lut_df` to get the new values from.
map_value <- function(input_value, lut_df, key_col, value_col) {
    key_col_string <- deparse(substitute(key_col))
    value_col_string <- deparse(substitute(value_col))

    match_index <- match(input_value, lut_df[[key_col_string]])

    if (!is.na(match_index)) {
        return(lut_df[[value_col_string]][match_index])
    } else {
        return(NA) 
    }
}

#' Find Best Match in A Vector of Strings using Automatic Gap Detection
#' 
#' @author Zhehao Hu
#' 
#' @importFrom stringdist stringdistmatrix
#' 
#' @param pattern A pattern to look for.
#' @param strings A vector of strings to look for the pattern.
#' @param index Logical. Whether to return index of the best matching item in the string vector (TRUE) or return the best matching item value (FALSE).
#' @param silent Logical. Set to TRUE to disable result printing.
find_best_match <- function(pattern, strings, index = F, silent = F) {
    name.dist <- stringdistmatrix(strings, pattern) %>% as.numeric() 
    name.dist.sorted <- name.dist %>% sort()
    weights <- 1/log(name.dist.sorted[-1]+1) # weight gap significance decreasingly while upper value of the gap increases
    maxGapIndex.sorted <- diff(name.dist.sorted)*weights %>% which.max()
    threshold <- mean(name.dist.sorted[c(maxGapIndex.sorted, maxGapIndex.sorted+1)]) # identify the threshold of the gap
    pattern.match <- strings[name.dist < threshold]

    cat("\nPattern given as:", pattern, "\nMatched", pattern.match, "\n\n")

    if (index) {
        return(which(name.dist < threshold))
    } else {
        return(pattern.match)
    }
}

#' Create Configuration File
#' 
#' Several functions in this package requires user-specific configuration to work. Setting up a configuration file could reduce repetitive manual input in workflow.
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_abort
#' @importFrom yaml write_yaml
create_config <- function() {
    if (length(list.files(".", "HPLNC_config.yaml"))>0) {
        cli_alert_info(paste0("Configuration file already exists: ", list.files(".", "HPLNC_config.yaml")[[1]], ", overwrite?"))
        resp <- readline("Y/n >>> ")
        if (resp != "Y") {cli_abort("Aborted. Please modify the existing configuration.")}
    }
    config <- list(
        "ABGD_executable" = "",
        "cache_folder" = ""
    )
    write_yaml(config, "HPLNC_config.yaml")
}

#' @importFrom yaml read_yaml
#' @importFrom cli cli_abort
read_config <- function(field, silent = F) {
    config <- read_yaml("HPLNC_config.yaml")
    if (!field %in% names(config)) {
        if (silent) {
            return(NULL)
        } else {
            cli_abort(paste0("Field not found in configuration file: ", field))}
    } else {
        return(config[[field]])
    }
}

#' Load cache files
#' 
#' Load cache file according to file name under the cache folder specified by user configuration.
load_cache <- function(name) {
    cache_folder <- read_config("cache_folder", T)
    if (is.null(cache_folder) | cache_folder == "") {
        cache_folder <- "cache"
    }
    target_file <- list.files(cache_folder, pattern = name, full.names = T, recursive = T)
    if (length(target_file) > 1) {cli::cli_alert_warning(paste0("Found more than one cache file under pattern '", name, "', loading '", target_file[[1]], "'."))}
    if (length(target_file) == 0) {return(NULL)}
    return(readRDS(target_file[[1]]))
}