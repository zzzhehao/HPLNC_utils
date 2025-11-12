#' Extract INSD Metadata from `xml` File Downloaded From NCBI Genebank
#' 
#' @param xml File path to `xml` file
#' @return A data.frame of flattened tabular metadata.
#' @import rlist
#' @import tidyverse
#' @import XML
#' @import xml2
extract_INSD_metadata <- function(xml) {
    suppressMessages({
        ncbi.xml <- read_xml(xml)
        ncbi.xml.parsed <- xmlParse(ncbi.xml)
        ncbi.ls <- xmlToList(ncbi.xml.parsed)

        feature.table <- map(ncbi.ls, \(seq) {
        ele.ls <- seq$`INSDSeq_feature-table`
        feat.keys <- ele.ls %>% list.select(INSDFeature_key)
        map(feat.keys, \(key) {
            feat.ls <- ele.ls %>%
            list.filter(INSDFeature_key == key) %>%
            list.select(INSDFeature_quals)
            feat.ls.fl <- map(feat.ls, \(feat.ls) {
            map(feat.ls[[1]], \(qua.entry) {
                list(qua.entry$INSDQualifier_value) %>% `names<-`(paste(key, qua.entry$INSDQualifier_name, sep = "_"))
            }) %>% flatten()
            }) %>% flatten()
        }) %>% flatten()
        }) %>% bind_rows()

        ref.table <- map(ncbi.ls, \(seq) {
            ref.ls <- seq$INSDSeq_references
            if (length(ref.ls) > 1) { # if there are more than one reference, only one of them should be the final reference.
                ref.entry <- ref.ls %>%
                list.filter(str_detect(INSDReference_journal, "Submitted")) %>%
                .[1]
                if ( # if there is no submitted reference, or the submitted reference is the direct submission before manuscript publication, then just take one, as the direct submission will override the publication title.
                    any(c(
                        is.null(unlist(ref.entry)),
                        ref.entry$INSDReference$INSDReference_title == "Direct Submission"
                ))) { 
                    ref.entry <- ref.ls %>%
                        list.filter(INSDReference_reference == ref.ls %>%
                        list.select(INSDReference_reference) %>%
                        .[[1]])
                }
            } else { # there is only one reference
                ref.entry <- ref.ls %>%
                    list.filter(INSDReference_reference == ref.ls %>%
                    list.select(INSDReference_reference) %>%
                    .[[1]])
            }
            data.frame(
                Ref_authors = paste(ref.entry$INSDReference$INSDReference_authors, collapse = ";"),
                Ref_title = ref.entry$INSDReference$INSDReference_title,
                Ref_journal = ref.entry$INSDReference$INSDReference_journal
            ) %>% flatten()
        }) %>% bind_rows()

        primary.info <- ncbi.ls %>%
        list.select(
            `INSDSeq_primary-accession`,
            `INSDSeq_length`,
            `INSDSeq_definition`,
            `INSDSeq_organism`,
            `INSDSeq_taxonomy`
        ) %>%
        bind_rows()

        INSDMetadata <-
        bind_cols(primary.info, feature.table, ref.table) %>%
        dplyr::rename("INSDSeq_primary-accession" = "`INSDSeq_primary-accession`") %>%
        select(-contains("db_xref"))
    })
    return(INSDMetadata)
}

#' Wrapper Function to Read All `XML` File and Update Information to Databse.
#' 
#' @param XML_folder Path to folder containing `xml` files. All `xml` files in the folder will be loaded.
#' @param regenerate logical. Whether to overwrite the table in database (i.e. regenerating the table). Default to FALSE, which will only add metadata that has not been added in the table.
#' @import tidyverse
#' @import DBI
#' @return Nothing. 
update_INSD_metadata <- function(XML_folder = "data/sequence/INSD", regenerate = F) {
    filter_c <- function(metadata.seq){
        metadata.seq %>%
            mutate(
                c_filter_1 = ifelse(str_detect(INSDSeq_definition, "\\bcytochrome\\b|\\b18S\\b|\\b28S\\b|\\b16S\\b", F), T, F),
                c_filter_2 = ifelse(str_detect(INSDSeq_definition, "unverified|UNVERIFIED", T), T, F) # verified
            )
    }
    
    extract_voucher <- function(metadata.seq){
        labels <- metadata.seq %>%
            mutate(
                c_voucher = ifelse(is.na(source_specimen_voucher), source_isolate, source_specimen_voucher),
                c_organism_id = case_when(
                    is.na(c_voucher) ~ source_organism,  # prevent literal NA for those don't have any specimen voucher info, overides next case, which includes NA cases
                    str_detect(source_organism, paste0("\\b", c_voucher, "\\b"), T) ~ paste(source_organism, c_voucher), 
                    .default = source_organism
                )
            ) %>%
            dplyr::select(c("INSDSeq_primary-accession", "c_voucher", "c_organism_id"))
        left_join(metadata.seq, labels, by = "INSDSeq_primary-accession")
    }
    
    extract_gene <- function(metadata.seq){
        gene_ex <- metadata.seq %>% 
            mutate(
                c_gene = case_when(
                    str_detect(INSDSeq_definition, "\\bcytochrome\\b") ~ "COI",
                    str_detect(INSDSeq_definition, "\\b18S\\b") ~ "18S",
                    str_detect(INSDSeq_definition, "\\b28S\\b") ~ "28S",
                    str_detect(INSDSeq_definition, "\\b16S\\b") ~ "16S",
                )) %>%
            dplyr::select(c("INSDSeq_primary-accession", "c_gene"))
        left_join(metadata.seq, gene_ex, by = "INSDSeq_primary-accession")
    }

    INSD.files <- list.files(XML_folder, pattern = "\\b.xml$")

    paste0("Loading following XML files: \n\n", paste(INSD.files, collapse = "\n"), collapse = "") %>% cat()
    insd <- extract_INSD_metadata(paste(XML_folder, INSD.files[1], sep = "/"))

    tbl.INSD.raw <- map(INSD.files, \(xml) {
        extract_INSD_metadata(paste(XML_folder, xml, sep = "/"))
    }) %>% bind_rows()

    # Filtering, extracting steps as documented in log 0829.qmd
    tbl.INSD <- tbl.INSD.raw %>%
        filter_c() %>%
        extract_voucher() %>%
        extract_gene()

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")
    if (!"metadata.Sequence.NCBI" %in% dbListTables(HPLNCdb)) { # Initialize
        print("No INSD metadata found, initializing ... ")
        dbWriteTable(HPLNCdb, "metadata.Sequence.NCBI", tbl.INSD)
        return()
    } else {
        print("Updating INSD metadata ... ")
        if (regenerate) {
            print("Overwriting INSD metadata ... (regenerate = T)")
            dbWriteTable(HPLNCdb, "metadata.Sequence.NCBI", tbl.INSD, overwrite = T)
            DBchange_sign(
                "metadata.Sequence.NCBI", 
                "Hard", 
                "INSD metadata of NCBI sequences has been regenerated.", 
                Sys.Date() %>% gsub("-", "", .) %>% paste0("INSD"))
            print("INSD metadata table has been successfully regenerated.")
            dbDisconnect(HPLNCdb)
            return()
        }
        INSD.registered <- tbl(HPLNCdb, "metadata.Sequence.NCBI") %>% collect()
        INSD.new <- tbl.INSD %>% 
            dplyr::filter(!`INSDSeq_primary-accession` %in% INSD.registered$`INSDSeq_primary-accession`)
        INSD.updated <- bind_rows(INSD.registered, INSD.new) 
        dbWriteTable(HPLNCdb, "metadata.Sequence.NCBI", INSD.updated, overwrite = T)
        DBchange_sign(
            "metadata.Sequence.NCBI", 
            "Append", 
            "INSD metadata of NCBI sequences has been updated.", 
            Sys.Date() %>% gsub("-", "", .) %>% paste0("INSD"))
        print("INSD metadata table has been successfully updated.")
        dbDisconnect(HPLNCdb)
        return()
    }
}
