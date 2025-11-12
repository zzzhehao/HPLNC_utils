#' Generate Snapshot from All Tables
#' @import DBI
#' @import tidyverse
#' @import RSQLite
DBsnapshot <- function() {

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

#' Return Signature of the Last Request
#' @import DBI
#' @import tidyverse
#' @import RSQLite
DBlastRequest <- function() {
    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    lastRequest <- tbl(HPLNCdb, "signature") %>% collect() %>% arrange(desc(requestId)) %>% pull(requestId) %>% .[1]

    return(lastRequest)
    dbDisconnect(HPLNCdb)
}

#' Pull Table
#' @param table Table name. Type "show_table" to display all available table names in database.
#' @param cleaned Clean table, rows marked in `del` column will be dropped.
#' @param formatting Format table. Utilize table-specific formatting function to format the table if available. Formatting functions are always named under the rule `format_` + table name. 
#' @param format Argument to pass over to formatting functions. If available different formatting could be chosen.
#' 
#' @import tidyverse
#' @import DBI
#' @import RSQLite
DBpullTable <- function(table, cleaned = T, formatting = T, format = NULL) {
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

