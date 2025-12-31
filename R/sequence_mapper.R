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
    metadata.seq <- db_pull("metadata.Sequence.NCBI", T, T)
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
    metadata.spec <- db_pull("metadata.Specimen.Haploniscidae") %>%
        mutate(
            label_genetics = paste(gensp_morpho_ZH, DZMB2HH, sep = "_") %>%
                gsub("\\.", "", .) %>% gsub(" ", "_", .)
        )
    sequence.map <- sequence.map %>% 
        mutate(
            id = 1:nrow(sequence.map),
            c_organism_label = c_organism_id %>% 
                gsub(" ", "_", .) %>% 
                gsub("\\.", "", .) %>% 
                gsub(":", "-", .) %>% 
                gsub("^Haploniscidae_sp_[A-Z]{3,5}[0-9]{3,4}-[0-9]{2}_", "", .), # CCZ batch with ultra long name
            c_organism_label = case_when(
                str_detect(c_organism_id, "^[A-Z]{3}[0-9]{3}$") ~ {
                    input_key <- str_sub(c_organism_id, 4, 6) %>% as.numeric()
                    map_vec(input_key, ~map_value(.x, metadata.spec, voucher, label_genetics)) %>% paste0(., "_", c_organism_id)
                },
                .default = c_organism_label
            )
        )
    

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


#' Generate Sequence Metadata LUT
#' 
#' @import dplyr
#' @import purrr
#' @import tidyr
generate_sequence_LUT <- function() {
    specmeta <- db_pull("metadata.Specimen.Haploniscidae")
    sequence.map <- db_pull("sequence.map", F, T)

    # Construct c_organims_id and c_organism_label to match with sequence.map
    specmeta.label <- specmeta %>% 
        dplyr::select(c("voucher", "gensp_morpho_ZH")) %>% 
        dplyr::filter(!is.na(voucher)) %>% 
        mutate( 
            c_organism_id = paste0("ZHH", str_pad(as.character(voucher), 3, "left", "0")),  
            c_organism_label = gensp_morpho_ZH,
            .keep = "none")

    # match sequence.map
    sequence.map.LUT <- sequence.map %>%
        rows_update(specmeta.label, by = "c_organism_id") %>%
        # pivot longer
        pivot_longer(starts_with("c_gene_"), names_to = "gene", values_to = "identifier") %>% 
        # choose desired columns
        dplyr::select(c("identifier", "c_organism_id", "c_organism_label"))

    # Extract IDivA subset
    sequence.map.LUT.IDivA <- sequence.map.LUT %>% dplyr::filter(str_detect(identifier, "^[0-9]{3}_.{3}$"))
    # multiplies IDivA subset for two voucher prefixes
    sequence.map.LUT <- list(
        sequence.map.LUT %>% dplyr::filter(str_detect(identifier, "^[0-9]{3}_.{3}$", T)),
        sequence.map.LUT.IDivA %>% mutate(identifier = paste0("VPS", identifier)),
        sequence.map.LUT.IDivA %>% mutate(identifier = paste0("ZHH", identifier))
    ) %>%
        purrr::reduce(bind_rows) %>%
        filter(!is.na(identifier)) %>% 
        mutate(c_organism_label = case_when(str_detect(identifier, "^[A-Z]{3}[0-9]{3}_") ~ paste(c_organism_label, gsub("_.{3}$", "", identifier), sep = "_"), .default = c_organism_label))
    return(sequence.map.LUT)
}

generate_orglab_LUT <- function(tip.taxa) {
    seqmap <- db_pull("sequence.map")
    seqmap <- seqmap %>% 
        mutate(c_organism_label = as.numeric(c_organism_label)) %>%
        filter(!is.na(c_organism_label)) %>%
        dplyr::select(c_organism_label, c_organism_id)
    seqmap %>% 
        filter(c_organism_label %in% tip.taxa & is.na(as.numeric(seqmap$c_organism_id)))
}
