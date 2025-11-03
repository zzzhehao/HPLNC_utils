DBchange_sign <- function(table, type, msg, request = request.id, con = dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite"), time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), signature.table = "signature") {
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

DBexecute <- function(request.filename) {

    library(yaml)
    library(rlist)
    library(glue)
    library(rlang)
    library(tidyverse)
    library(DBI)
    
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
        
        DBchange_sign(table, type, msg, request.id)
        dbWriteTable(HPLNCdb, table, target.updated, overwrite = T)
    })

    print("All request has been completed.")
    dbDisconnect(HPLNCdb)
}

DBmanual_pull <- function(table) {
    require(tidyverse)
    require(DBI)

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    tbl <- tbl(HPLNCdb, table) %>% collect() 
    dbDisconnect(HPLNCdb)
    return(tbl)
}

DBmanual_write <- function(tbl, table, msg, request.id) {
    require(tidyverse)
    require(DBI)

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")

    DBchange_sign(table, "Manual", msg, request.id)
    dbWriteTable(HPLNCdb, table, tbl, overwrite = T)
    dbDisconnect(HPLNCdb)
}