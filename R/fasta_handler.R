#' Clean Sequence Labels in Alignment
#' 
#' @description
#' This only applies to DZMB collection where the sequences are managed using field ID. The gene will be determined by the file name (whether it's COI, 18S, or 28S) and a gene suffix will be added to the voucher (e.g. VPS101_COI, ZHH003_18S). 
#' 
#' @param fasta_path Path to FASTA file
#' @import ape
#' @import stringr
#' @return A new FASTA file will be written in the same directory. The path to the new file will be returned.
clean_label <- function(fasta_path){
    genes <- c("COI", "18S", "28S")
    gene <- genes[str_detect(basename(fasta_path), paste(genes, sep = "|"))]
    aln <- ape::read.FASTA(fasta_path) 
    names(aln) <- names(aln) %>%
        ifelse(
            str_detect(., "^VPS[0-9]{3}"), 
            str_extract(., "^VPS[0-9]{3}") %>% paste0(., "_", gene), 
            .) %>%
        ifelse(str_detect(., "\\.1$"), gsub("\\.1$", "", .), .)
    new_fasta_path <- paste0(dirname(fasta_path), "/", gsub("\\.fasta$", "", basename(fasta_path)), "_cleaned.fasta")
    ape::write.FASTA(aln, new_fasta_path)
    return(new_fasta_path)
}

#' Run Gblocks to Trim the Alignment
#' 
#' @description
#' This is a wrapper function to call Gblocks in CLI. For details about the program please visit: https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html
#' 
#' @param fasta_path File path of the alignment FASTA file.
#' @param gblocks_path File path of the Gblocks excutable.
#' @param t Type of sequence. Gblocks command line parameter, will be passed over as `-t`.
#' @param b1 Minimum number of sequence for a conserved position (>= 50%+1 of the number of sequences, default 50%+1 of the number of sequences). Gblocks command line parameter, will be passed over as `-b1`. Also accepts `min` to choose the lowest value, which maximize the selected number of positions (i.e. least aggressive in trimming).
#' @param b2 Minumum number of sequences for a flank position (>= b1, default 85% of the number of sequences). Gblocks command line parameter, will be passed over as `-b2`. Also accepts `min` to choose the lowest value, which maximize the selected number of positions (i.e. least aggressive in trimming).
#' @param b3 Maximum number of contiguous nonconserved positions (gblocks default to 8, here default to 20). Bigger values increase the selected number of positions.
#' @param b4 Minimum length of a block after gap cleaning (gblocks default to 10, minimum value 2, here default to 2). Bigger values decrease the selected number of position. 
#' @param b5 Allowe gap positions. `n` for absolute no gap positions allowed in the final alignment. `h` for only positions where 50% or more of the sequens have a gap are treated as gap position and are eliminated. `a` for all gap positions can be selected, no positions are eliminated becaused of the gaps.
#' @param args Other arguments that are passed directly to CLI command. 
#' 
#' @import ape
#' 
#' @return Trimmed alignment will be written in new FASTA file under the same directory with suffix `_gbtrimmed.fasta`. The path to the new file will be returned.
gblocks <- function(fasta_path, gblocks_path = "/Users/hu_zhehao/Desktop/Biology-tools/Gblocks_0.91b/Gblocks", t="d", b1="min", b2="min", b3=20, b4=2, b5="a", args=""){
    MSA.dirname <- dirname(fasta_path)
    genes <- c("18S", "28S")
    gene <- genes[str_detect(basename(fasta_path), paste(genes, sep = "|"))]

    b1 <- as.character(b1)
    b2 <- as.character(b2)

    if (b1 == "min") {
        nseq <- length(ape::read.FASTA(fasta_path))
        b1 <- ceiling(nseq/2)+1
    } else {
        b1 <- as.numeric(b1)
        if (is.na(b1)) {stop("Argument b1 only accepts number or `min`")}
    }
    if (b2 == "min") {
        b2 <- b1
    } else {
        b2 <- as.numeric(b2)
        if (is.na(b2)) {stop("Argument b2 only accepts number or `min`")}
    }

    parameters <- paste0("-t=", t, " -b1=", b1, " -b2=", b2, " -b3=", b3, " -b4=", b4, " -b5=", b5, " ", args)
    cat("Parameter used: ", parameters, sep = "")
    gblocks.cmd <- paste(
        gblocks_path,
        fasta_path, parameters) 
    
    gblocks.cmd.basename <- paste0(MSA.dirname,"/gblocksCMD_", gene, ".sh")
    write(gblocks.cmd, gblocks.cmd.basename)
    
    system(paste0("sh ", gblocks.cmd.basename))
    
    files <- list.files(MSA.dirname, full.names = T)
    geneFile <- files[str_detect(files, gene)]
    gbFile <- geneFile[str_detect(geneFile, pattern = "fasta-gb$")]
    if (length(gbFile) > 1) {stop("Multiple fasta-gb files.")}
    
    trimAln <- ape::read.FASTA(gbFile)
    
    new_fasta_name <- gsub(".fasta-gb$", "_gbtrimmed.fasta", gbFile)
    
    ape::write.FASTA(trimAln, new_fasta_name)
    return(new_fasta_name)
}

#' Convert NEXUS alignment to FASTA alignment
#' @param nexus_path Path to NEXUS file.
#' @param fasta_path Path to FASTA file.
#' @param safe_name logical. Whether to modify the name to remove "_".
#' @param remove_outgroup Vector of sequence labels to remove. Default to `NULL` without removing any sequence. This happens before modifying sequence labels if `safe_name` is set to TRUE.
#' @return The FASTA file will be written in the same directory. The path will be returned.
nexus2fasta <- function(nexus_path, fasta_path=NULL, safe_name=F, remove_outgroup=NULL) {
    if (is.null(fasta_path)) {fasta_path <- gsub("\\.nexus$", ".fasta", nexus_path)}
    aln <- ape::read.nexus.data(nexus_path) %>% ape::as.DNAbin()

    if (!is.null(remove_outgroup)) {
        aln <- aln[!names(aln) %in% remove_outgroup]
    }

    if (safe_name) {
        name <- names(aln) %>% gsub("_", "", .)
        names(aln) <- name
    }
    ape::write.FASTA(aln, fasta_path)
    return(fasta_path)
}