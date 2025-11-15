#' Create a Summary Table Mapping Sequence Label/Accession Number with Vouchered Animals.
#' 
#' @description
#' This will summarize the sequence labels in the way that the sequnces originated from the same vouchered animals could be linked together, which is essential for concatenate multilocus alignment. 
#' 
#' @return Nothing. A data.frame will be written into the database. 
#' @param regenerate logical. Whether to overwrite existing table in the database. Otherwise only the rows with vouchered animal not yet in the database will be added. 
#' @import DBI
#' @import dplyr
#' @import RSQLite
#' @import tidyr
#' @import stringr
update_sequence_map <- function(regenerate = F){
    # Generate wide table
    metadata.seq <- db_pull("metadata.Sequence.NCBI", F, F)
    sequence.map.ncbi <- metadata.seq %>% filter(!is.na(c_gene)) %>% pivot_wider(id_cols = c_organism_id, names_from = c_gene, names_prefix = "c_gene_", values_from = `INSDSeq_primary-accession`)

    # Combine with vps names
    allvps <- db_pull("metadata.Specimen.Haploniscidae") %>% pull(voucher) %>% sort() %>% as.character() %>% str_pad(3, "left", "0")
    vps.map <- data.frame(c_organism_id = allvps) %>%
    mutate(
        c_gene_COI = paste0(c_organism_id, "_COI"),
        c_gene_16S = paste0(c_organism_id, "_16S"),
        c_gene_18S = paste0(c_organism_id, "_18S"),
        c_gene_28S = paste0(c_organism_id, "_28S"),
    ) 

    sequence.map <- bind_rows(sequence.map.ncbi, vps.map)
    # generate label without special character for traking the mrbayes workflow
    sequence.map <- sequence.map %>% 
        mutate(
            id = 1:nrow(sequence.map),
            c_organism_label = c_organism_id %>% 
                gsub(" ", "_", .) %>% 
                gsub("\\.", "", .) %>% 
                gsub(":", "-", .) %>% 
                gsub("^Haploniscidae_sp_[A-Z]{3,5}[0-9]{3,4}-[0-9]{2}_", "", .), # CCZ batch with ultimate long name
            c_organism_label = ifelse(nchar(c_organism_label)>50, id, c_organism_label)
        )
    # print(sequence.map$c_organism_label)

    HPLNCdb <- dbConnect(RSQLite::SQLite(), "data/database/HPLNCdb.sqlite")
    if (!"sequence.map" %in% dbListTables(HPLNCdb)) { # Initialize
        print("No sequence map found, initializing ... ")
        dbWriteTable(HPLNCdb, "sequence.map", sequence.map)
        dbDisconnect(HPLNCdb)
        return()
    } else {
        if (regenerate) {
            print("Regenerating sequence.map table ... ")
            dbWriteTable(HPLNCdb, "sequence.map", sequence.map, overwrite = T)
            dbDisconnect(HPLNCdb)
            return()
        }
        sequence.map.registered <- db_pull("sequence.map", F, F)

        sequence.map.new <- sequence.map %>% 
            dplyr::filter(!c_organism_id %in% sequence.map.registered$c_organism_id)
        sequence.map.updated <- bind_rows(sequence.map.registered, sequence.map.new) 
        dbWriteTable(HPLNCdb, "sequence.map", sequence.map.updated, overwrite = T)
        db_sign(
            "sequence.map", 
            "Append", 
            "Sequence map (accession number of each specimen voucher) has been updated.", 
            Sys.Date() %>% gsub("-", "", .) %>% paste0("sequenceMap"))
        print("Sequence map has been successfully updated.")
        dbDisconnect(HPLNCdb)
        return()
    }
}