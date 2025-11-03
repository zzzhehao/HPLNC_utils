DBsnapshot <- function(cleaned = F, formatted = F) {

    library(tidyverse)
    library(DBI)

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    dir <- paste0("data/metadata/latest/", format(Sys.time(), "%Y%m%d_%H%M%S"), "/")
    dir.create(dir)

    table <- dbListTables(HPLNCdb)
    table <- table[table != "signature"]

    walk(table, \(tbl.name) {
        tbl <- tbl(HPLNCdb, tbl.name) %>% collect()
    
        write.table(
            tbl, 
            paste0(
                dir,
                tbl.name, 
                ".csv"), 
            sep = ";", 
            row.names = F, 
            col.names = T)
    })

    signature <- tbl(HPLNCdb, "signature") %>% collect()
    write.table(signature, paste0(dir, "signature.csv"), sep = ";", row.names = F, col.names = T)
    dbDisconnect(HPLNCdb)
}

DBlastRequest <- function() {
    library(tidyverse)
    library(DBI)

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    lastRequest <- tbl(HPLNCdb, "signature") %>% collect() %>% arrange(desc(requestId)) %>% pull(requestId) %>% .[1]

    return(lastRequest)
    dbDisconnect(HPLNCdb)
}

DBpullTable <- function(table, cleaned = T, formatting = T, format = NULL) {
    require(tidyverse)
    require(DBI)

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    if (table == "show_table") {
        print(dbListTables(HPLNCdb))
        return("Show tables only.")
    }
    tbl <- tbl(HPLNCdb, table) %>% collect() 

    if (cleaned) {
        if ("del" %in% colnames(tbl)) {
            tbl <- tbl %>%
                dplyr::filter(is.na(del)) %>% 
                dplyr::select(-c("del"))
        }
    }

    if (formatting) {
        tbl_formatter <- paste0("format_", table) 
        if (tbl_formatter %in% ls()) { # check availability of formatting function
            if (is.null(format)) {
                format <- ""
            } else {
                format <- paste0(", format = ", format)
            }
            tbl <- paste0(tbl_formatter, "(tbl, formatting = T", format, ")") %>%
                rlang::parse_expr() %>%
                rlang::eval_tidy()
        }
    }
    return(tbl)
}

