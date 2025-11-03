# log 0829.qmd
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
    return()
}

gblocks <- function(fasta_path, t="d", b1="min", b2="min", b3=20, b4=2, b5="a", args=""){
    MSA.dirname <- dirname(fasta_path)
    genes <- c("18S", "28S")
    gene <- genes[str_detect(basename(fasta_path), paste(genes, sep = "|"))]

    if (b1 == "min") {
        nseq <- length(ape::read.FASTA(fasta_path))
        b1 <- ceiling(nseq/2)+1
    }
    if (b2 == "min") {
        b2 <- b1
    }

    # Gblocks Parameters
    # b1 Minimum number of sequences for a conserved position (>= 50%+1 of number of sequences, def 50%+1)
    # b2 minumum number of sequences for a flank position (>= b1, def 85%)
    # b3 maximum number of contiguous nonconserved positions (def 8)
    # b4 minimum length of a block (def 10, >= 2)
    # b5 allowed gap positions (n = None, h = with half, a = all)

    # https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html#command_line

    parameters <- paste0("-t=", t, " -b1=", b1, " -b2=", b2, " -b3=", b3, " -b4=", b4, " -b5=", b5, " ", args)
    cat("Parameter used: ", parameters, sep = "")
    gblocks.cmd <- paste(
        "/Users/hu_zhehao/Desktop/Biology-tools/Gblocks_0.91b/Gblocks", 
        fasta_path, parameters) 
    
    gblocks.cmd.basename <- paste0(MSA.dirname,"/gblocksCMD_", gene, ".sh")
    write(gblocks.cmd, gblocks.cmd.basename)
    
    system(paste0("sh ", gblocks.cmd.basename))
    
    files <- list.files(MSA.dirname, full.names = T)
    geneFile <- files[str_detect(files, gene)]
    gbFile <- geneFile[str_detect(geneFile, pattern = "fasta-gb$")]
    if (length(gbFile) > 1) {stop("Multiple fasta-gb files.")}
    
    trimAln <- read.FASTA(gbFile)
    
    new_fasta_name <- gsub(".fasta-gb$", "_gbtrimmed.fasta", gbFile)
    
    ape::write.FASTA(trimAln, new_fasta_name)
}

nexus2fasta <- function(nexus_path, fasta_path=NULL, secure_name=F, remove_outgroup=NULL) {
    if (is.null(fasta_path)) {fasta_path <- gsub("\\.nexus$", ".fasta", nexus_path)}
    aln <- ape::read.nexus.data(nexus_path) %>% ape::as.DNAbin()
    if (secure_name) {
        name <- names(aln) %>% gsub("_", "", .)
        names(aln) <- name
    }
    if (!is.null(remove_outgroup)) {
        aln <- aln[!names(aln) %in% remove_outgroup]
    }
    ape::write.FASTA(aln, fasta_path)
    return(fasta_path)
}