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
