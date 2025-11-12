<<<<<<< HEAD
#' Generate path to given date's asset folder in log (is this really necessary?)
=======
#' Create assets folder for the day
#' @importFrom stringr str_pad
>>>>>>> cc18640 (new util function and update document)
date2assets <- function(date) {
    paste0("docs/logs/assets/", str_pad(as.character(date), 4, "left", "0"))
}

#' Query specimen info with DZMB2HH ID
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
ID.inquery <- function(ID) {
    DBpullTable("metadata.Specimen.Haploniscidae") %>% filter(DZMB2HH %in% ID) %>% t()
}

#' Query specimen info with voucher
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
voucher.inquery <- function(ID) {
    DBpullTable("metadata.Specimen.Haploniscidae") %>% filter(voucher %in% ID) %>% t()
}

#' Clear temporary cache file of R
clearRCache <- function() {
    unlink(tempdir(), recursive = TRUE)
    dir.create(tempdir())
}

#' Initialize daily log
#' @importFrom lubridate today
create_log <- function() {
    date <- format(today(), "%m%d")
    dir.create(date2assets(date))
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
    morphocheck <- readxl::read_xlsx("docs/notes/assets/sandbox/ZHH_morphocheck.xlsx", range = readxl::cell_cols("A:O"))
    file_path <- paste0("data/metadata/archive/morphocheck_", as.character(format(today(), "%Y%m%d"))) # make archive
    write.table(morphocheck, file_path, sep = ";", row.names = F)
    morphorcheck <- read.table(file_path, sep = ";", header = T)
    morphocheck <- format_morphocheck(morphocheck)
    return(morphocheck)
}

#' Generate Species Summary
#' @import ggtext
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

    map <- bathy.basemap(resolution = 10) +
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
create_summary_sp <- function() {
    suppressMessages({suppressWarnings({
        metadata <- DBpullTable("metadata.Specimen.Haploniscidae")
        metadata.ebs <- DBpullTable("metadata.Station")
        morphocheck <- load_morphocheck()
    
        metadata.merge <- left_join(metadata, metadata.ebs %>% dplyr::select(c("station", "latStartDec", "longStartDec", "depthZone", "area_abbr_ZH")), by = "station") %>%
            left_join(morphocheck %>% dplyr::select(-c("voucher")), by = "DZMB2HH")
    
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
