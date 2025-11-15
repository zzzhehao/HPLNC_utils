#' Generate Signature for Database Change
#' 
#' @param table Table name subjected to change
#' @param type Change type
#' @param request Request ID
#' @param con Connection to database as produced by `DBI::dbConnect()`
#' @param time Time of change. Default to current time.
#' @param signature.table Table name to write signature. 
#' @import DBI
#' @import RSQLite
db_sign <- function(table, type, msg, request = request.id, con = dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite"), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), signature.table = "signature") {
    if (!type %in% c("Initial", "Append", "Amend", "Del", "Hard", "Manual")) {
        stop(paste0("Unaccepted type: ", type))
    }
    
    signature <- data.frame(time = time, table = table, type = type, msg = msg, requestId = request)

    if (!signature.table %in% dbListTables(con)) {
        signature <- data.frame(time = time, table = table, type = type, msg = msg, requestId = request) %>% filter(!is.na(time))
        dbWriteTable(con, signature.table, signature, overwrite = T)
    } else {
        dbWriteTable(con, signature.table, signature, append = T)
    }
}

#' Execute YAML Change Request
#' 
#' @import yaml
#' @import rlist
#' @import glue
#' @import rlang
#' @import tidyverse
#' @import DBI
#' 
db_exe_request <- function(request.filename) {
    request.file <- paste0("data/metadata/request/", request.filename, ".yaml")
    
    requests <- yaml.load_file(request.file) %>% 
        flatten() %>% 
        list.map(data.frame(.)) %>% 
        bind_rows() %>% 
        tibble() %>% 
        mutate(
            requestId = paste0(request.filename, "_", str_pad(1:nrow(.), 2, pad = "0")),
            status = ifelse(
                requestId %in% (dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite") %>% tbl("signature") %>% collect() %>% pull(requestId)),
                "complete", 
                "pending"
            )
        )
    
    walk(1:nrow(requests), \(x) {
        request <- requests[x,]

        if (request$status == "complete") {
            return()
        }
    
        request.id <- paste0(request.filename, "_", str_pad(x, 2, pad = "0"))
    
        table <- request$table
        type <- request$type
        filter <- request$filter
        if (!is.null(filter)) {
            filter_expr <- rlang::parse_expr(filter)
        }
        msg <- request$msg
        amend_col <- request$amend_col
        amend_mutate <- request$amend_mutate
        hard_execute <- request$hard_execute
            
        if (!type %in% c("Del", "Amend", "Initial", "Append", "Hard")) {
            stop(paste0("Unknown type: ", type))
        }
    
        HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")
        
        target <- tbl(HPLNCdb, table) %>% collect()
    
        if (type == "Del") {
            target.updated <- target %>% 
                mutate(
                    del = ifelse(eval_tidy(filter_expr, data = .), 
                    paste(request.id), 
                    del)
                )
        }

        if (type == "Amend") {
            target.markedAsAmend <- target %>% 
                mutate(
                    amend = ifelse(eval_tidy(filter_expr, data = .), 
                    paste(request.id), 
                    amend)
                )
            
            target.updated <- 
                paste0(
                    "target.markedAsAmend %>% mutate(", amend_col, " = case_when(", # mutate + case_when to manipulate
                    filter, " ~ ", # filter expression to select target
                    amend_mutate, # mutate expression
                    ", .default = ", amend_col, "))", # set default for unchanged
                    collapse = "") %>% 
                rlang::parse_expr() %>% 
                eval_tidy()   
            
        }

        if (type == "Hard") {
            target.updated <-
                hard_execute %>%
                rlang::parse_expr() %>%
                eval_tidy()

        }

        print(t(request))
        cat("\nParsed Expression:\n\n")

        if (!is.null(filter)) {
            cat("Filter\n")
            print(rlang::parse_expr(filter))
        }
        if (!is.null(amend_mutate)) {
            cat("\nMutate\n")
            print(rlang::parse_expr(amend_mutate))
        }
        if (!is.null(hard_execute)) {
            cat("\nHard\n")
            print(rlang::parse_expr(hard_execute))
            cat("\n")
        }
    
        exec <- readline("Execute request? Type `y` to continue, type 'v' to preview. \n>>> ")
    
        if (exec != "y") {
            if (exec == "v") {
                View(target.updated)
                stop("Preview result")
            } else {
                stop("Unknown response, abort.")
            }
        }
        
        db_sign(table, type, msg, request.id)
        dbWriteTable(HPLNCdb, table, target.updated, overwrite = T)
    })

    print("All request has been completed.")
    dbDisconnect(HPLNCdb)
}

#' Pull Raw Table
#' 
#' @param table Table name.
#' @import tidyverse
#' @import DBI
db_pull_raw <- function(table) {
    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    tbl <- tbl(HPLNCdb, table) %>% collect() 
    dbDisconnect(HPLNCdb)
    return(tbl)
}

#' Write Table
#' 
#' @param tbl Table object to write.
#' @param table Table name.
#' @param msg Change message.
#' @param request.id Request ID associated to this change.
#' 
#' @import tidyverse
#' @import DBI
db_write <- function(tbl, table, msg, request.id) {
    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    db_sign(table, "Manual", msg, request.id)
    dbWriteTable(HPLNCdb, table, tbl, overwrite = T)
    dbDisconnect(HPLNCdb)
}

#' Generate Snapshot from All Tables
#' @import DBI
#' @import tidyverse
#' @import RSQLite
db_snapshot <- function() {

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
db_show_last_request <- function() {
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
db_pull <- function(table, cleaned = T, formatting = T, format = NULL) {
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
        if (tbl_formatter %in% ls("package:HPLNC")) { # check availability of formatting function
            print(paste0("Found formatter: ", tbl_formatter))
            if (is.null(format)) {
                format <- ""
            } else {
                format <- paste0(", format = ", format)
            }
            print(paste0("Executing: ", tbl_formatter, "(tbl, formatting = T", format, ")"))
            tbl <- paste0(tbl_formatter, "(tbl, formatting = T", format, ")") %>%
                rlang::parse_expr() %>%
                rlang::eval_tidy()
        }
    }
    return(tbl)
}

