update_sequence_map <- function(regenerate = F){
    source("functions/DBsnapshot.R")

    require(tidyverse)
    require(DBI)

    # Generate wide table
    metadata.seq <- DBpullTable("metadata.Sequence.NCBI", F, F)
    sequence.map.ncbi <- metadata.seq %>% pivot_wider(id_cols = c_organism_id, names_from = c_gene, values_from = `INSDSeq_primary-accession`)

    # Combine with vps names
    allvps <- DBpullTable("metadata.Specimen.Haploniscidae") %>% pull(voucher) %>% sort()
    vps.map <- data.frame(c_organism_id = allvps) %>%
    mutate(
        COI = paste0(c_organism_id, "_COI"),
        `18S` = paste0(c_organism_id, "_18S"),
        `28S` = paste0(c_organism_id, "_28S"),
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
        sequence.map.registered <- DBpullTable("sequence.map", F, F)

        sequence.map.new <- sequence.map %>% 
            dplyr::filter(!c_organism_id %in% sequence.map.registered$c_organism_id)
        sequence.map.updated <- bind_rows(sequence.map.registered, sequence.map.new) 
        dbWriteTable(HPLNCdb, "sequence.map", sequence.map.updated, overwrite = T)
        DBchange_sign(
            "sequence.map", 
            "Append", 
            "Sequence map (accession number of each specimen voucher) has been updated.", 
            Sys.Date() %>% gsub("-", "", .) %>% paste0("sequenceMap"))
        print("Sequence map has been successfully updated.")
        dbDisconnect(HPLNCdb)
        return()
    }
}