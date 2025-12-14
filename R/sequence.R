#' Clean Sequence Labels in Alignment
#' 
#' @description
#' This only applies to DZMB collection where the sequences are managed using field ID. The gene will be determined by the file name (whether it's COI, 18S, or 28S) and a gene suffix will be added to the voucher (e.g. VPS101_COI, ZHH003_18S). 
#' 
#' @param fasta_path Path to FASTA file
#' @importFrom ape read.FASTA
#' @importFrom ape write.FASTA
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract
#' @return A new FASTA file will be written in the same directory. The path to the new file will be returned.
#' @export
clean_label <- function(fasta_path, noprefix = F){
    genes <- c("COI", "18S", "28S", "16S")
    gene <- genes[str_detect(basename(fasta_path), paste(genes, sep = "|"))]
    aln <- ape::read.FASTA(fasta_path) 
    names(aln) <- names(aln) %>%
        # VPS formatting
        ifelse(
            str_detect(., "^VPS[0-9]{3}"), 
            str_extract(., "^VPS[0-9]{3}") %>% paste0(., "_", gene), 
            .) %>%
        # ZHH formatting
        ifelse(
            str_detect(., "^ZHH[0-9]{3}"), 
            str_extract(., "^ZHH[0-9]{3}") %>% paste0(., "_", gene), 
            .) %>%
        # NCBI accession number formatting (to primary)
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
#' @importFrom stringr str_detect
#' @importFrom ape read.FASTA
#' @importFrom ape write.FASTA
#' 
#' @return Trimmed alignment will be written in new FASTA file under the same directory with suffix `_gbtrimmed.fasta`. The path to the new file will be returned.
#' @export
exe_gblocks <- function(fasta_path, gblocks_path = "/Users/hu_zhehao/Desktop/Biology-tools/Gblocks_0.91b/Gblocks", t="d", b1="min", b2="min", b3=20, b4=2, b5="a", args=""){
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
#' @export
convert_nexus2fasta <- function(nexus_path, fasta_path=NULL, safe_name=F, remove_outgroup=NULL) {
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

#' Write MRCA Prior into Alignment Nexus File for BEAST2
#' 
#' BEAST2 produces full bifurcation ultrametric tree required to perform GMYC. MRCA Prior enforces certain taxa to be monophyly in the resulted tree, which serves as a method of defining the outgroups. The MRCA Prior information could be written in a nexus block and recognized in BEAUti and BEAST2.
#' 
#' @param nexus A path to the alignment nexus file.
#' @param monophyly_groups A list defining monophyly groups. Each element of this list defines one monophyly group and must be a list itself containing the following two components:
#'   \describe{
#'     \item{\code{name}}{A character string (length 1) specifying the group's name.}
#'     \item{\code{labels}}{A character vector of taxa belonging to that group, this should be identical to the ones in the alginment.}}
#' @importFrom purrr map_vec
#' @export
beast2_write_mrca <- function(nexus, monophyly_groups) {
    MRCA <- map_vec(monophyly_groups, \(g){
        name <- g[["name"]]
        labels <- g[["labels"]]
        group.block <- paste("taxset", name, "=", labels, collapse = " ") %>% paste0(., ";")
        return(group.block)
    }) %>%
        paste("BEGIN sets;", ., "END;", sep = "\n")   
    write(MRCA, nexus, append = T, sep = "\n")
}

#' Perform Automatic Barcoding Gap Detection (ABGD)
#' 
#' @param alignment Path to the FASTA file.
#' @param exe Path to the ABGD executable.
#' @param output_folder Path to output folder.
#' @param p Minimal a priori value (Pmin, default to 0.001). 
#' @param P Maximal a priori value (Pmax, default to 0.1). 
#' @param n Number of steps in \[Pmin, Pmax\] (default to 20)
#' @param d Distance to use (0: Kimura-2P, 1: Jukes-Cantor, 2: Tamura-Nei, 3:simple distance; default to 3). 
#' @param X Mininmum Slope Increase (default to 1.5)
#' @return Path to the output folder. 
#' @export
#' @seealso [abgd_read_result()]
#' @seealso [abgd_read_result_all()]
#' @seealso [abgd_summary()]
#' @seealso [abgd_delim()]
#' @return The output folder where ABGD CLI program will export the result
abgd_execute <- function(
    alignment, exe, output_folder = NULL, 
    a = T, p = 0.001, P = 0.1, n = 20, d = 3, X = 1.5, args = NULL
) {
    if (is.null(output_folder)) {
        output_folder <- paste0(dirname(alignment), "/ABGD")
    } 
    if (a) {
        args <- paste("-a", args)
    }

    abgd_cmd <- paste(
        exe,
        "-p", p,
        "-P", P,
        "-n", n,
        "-d", d,
        "-X", X,
        "-o", output_folder,
        args,
        alignment,
        collapse = " "
    ) %>% gsub("  ", " ", .)

    cli::cli_alert_info(paste0("Executing ABGD command: \n\n", abgd_cmd))
    system(abgd_cmd)

    # ABGD program produce the distmat.txt always at the working directory. Clean this mess up.
    distmat <- list.files(getwd(), paste0(basename(alignment) %>% gsub("\\.fasta$", "", .), "_distmat\\.txt$"), full.names = T)
    cat(distmat)
    distmat_output <- paste(output_folder, basename(distmat), sep = "/")
    cat(distmat_output)
    file.rename(distmat, distmat_output)

    return(output_folder)
}

#' Read ABGD Result
#' 
#' @description
#' I extracted the code from `{delimtools}` function `abgd_tbl` as it just does what this function should. I omit the execution part from that function because it doesn't provide full control of the program, especially the two parameters `p` and `P` which might be sensitive to datasets, and leaves `X` controllable, which should be left as default most of the time except the datasets are very noisy. 
#' @importFrom tidyr separate_longer_delim
#' @export
#' @seealso [abgd_execute()]
#' @seealso [abgd_read_result_all()]
#' @seealso [abgd_summary()]
#' @seealso [abgd_delim()]
#' @return A data frame of the delimitation result
abgd_read_result <- function(path) {
    delim <- readr::read_delim(path, delim = ";", col_names = c("labels"), col_types = "c")
    delim <- delim %>% 
        mutate(
            taxa = gsub("^id: ", "", X2),
            group = gsub("^Group\\[ ", "", labels) %>% 
                gsub(" \\].*$", "", .) %>%
                as.factor(),
            .keep = "none"
        ) %>%
        separate_longer_delim(taxa, " ")
    return(delim)
}

#' Read All ABGD Partition Results
#' 
#' @description
#' A bulk loading version of `abgd_read_result`.
#' 
#' @export
#' @seealso [abgd_execute()]
#' @seealso [abgd_read_result()]
#' @seealso [abgd_summary()]
#' @seealso [abgd_delim()]
#' @return A named list of summary for all partitions
abgd_read_result_all <- function(output_folder) {
    part.files <- list.files(output_folder, pattern = regex("part(init)?\\.[0-9]+\\.txt$"), full.names = T)
    cvs.file <- list.files(output_folder, pattern = regex("res.cvs$"), full.names = T)
    res <- read.table(cvs.file, T, "\t")
    
    abgd.summary <- map(part.files, \(path){
        partition <- ifelse(str_detect(path, "partinit"), "initial", "recursive")
        index <- str_extract(path, "[0-9]+\\.txt$") %>% gsub("\\.txt", "", .) %>% as.numeric()
        list(
            path = path,
            partition_type = partition, 
            index = index,
            group_count = length(levels(abgd_read_result(path)$group)),
            result = abgd_read_result(path)
        )
    })
    return(abgd.summary)
}

#' Generate Summary for ABGD Result
#' 
#' @param output_folder Path to the ABGD output folder.
#' @param silent Logical. Whether show plot or not.
#' 
#' @return A list of length 2.
#'   \describe{
#'     \item{\code{plot}}{A plot showing number of groups fround with each partition will be shown.}
#'     \item{\code{summary}}{A data.frame of number of groups and the prior intraspecific divergence of each partition will be returned.}
#' 
#' @importFrom rlist list.select
#' @export
#' @seealso [abgd_execute()]
#' @seealso [abgd_read_result()]
#' @seealso [abgd_read_result_all()]
#' @seealso [abgd_delim()]
#' @return A named list. \describe{
#'   \item{\code{plot}}{A ggplot object, the summary plot}
#'   \item{\code{summary}}{A data frame of the partition summary}
#' }
abgd_summary <- function(output_folder, silent = F) {
    abgd.res <- abgd_read_result_all(output_folder)
    cvs.file <- list.files(output_folder, pattern = regex("res.cvs$"), full.names = T)
    res <- read.table(cvs.file, T, "\t")
    abgd.summary <- abgd.res %>% 
        list.select(partition_type, index, group_count) %>% 
        bind_rows() %>% 
        arrange(partition_type, index) %>%
        pivot_wider(names_from = partition_type, values_from = group_count) %>%
        mutate(
            fine_group = recursive - initial,
            fine_group_prop = round(fine_group / recursive, 2)
        )
    
    abgd.summary <- bind_cols(abgd.summary, res[1:nrow(abgd.summary),]["prior"])
    summary.plot <- 
        abgd.summary %>% 
        filter(initial >= 2) %>%
        ggplot() +
        geom_point(aes(x = index, y = initial), shape = 1) +
        geom_point(aes(x = index, y = recursive), shape = 16) +
        geom_text(
            aes(x = index, y = max(abgd.summary$recursive)*1.1, label = round(prior, 4)),
            angle = 45, hjust = 0, size = 3.4
        ) +
        scale_x_continuous(n.breaks = nrow(abgd.summary %>% 
            filter(initial >= 2)), expand = expansion(add = c(1, 1.5))) +
        scale_y_continuous(expand = expansion(mult = c(0.2, 0.5))) +
        theme_minimal() +
        theme(
            aspect.ratio = 0.5,
            axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
            axis.line.y.left = element_line(color = "black", linewidth = 0.5),
            panel.grid.minor.x = element_blank()
        ) +
        labs(x = "Partition", y = "Groups found", title = "ABGD Number of Groups", subtitle = "Solid point: recursive partition; hollow point: initial partition. \nNumber above points: prior intraspecific divergence")
    if (!silent) {print(summary.plot)}
    return(list(plot = summary.plot, summary = abgd.summary))
}

#' Perform Species Delimitation with ABGD Result
#' 
#' @param output_folder Output folder of the ABGD result.
#' @param indices A numerical vector indicating the indices of the partitions to be evaluated.
#' @param label A character vector. Name to added to the column name of the output data frame. The length must match with \code{indices}.
#' 
#' @importFrom purrr map2
#' @importFrom rlist list.filter
#' @export
#' @seealso [abgd_execute()]
#' @seealso [abgd_read_result()]
#' @seealso [abgd_read_result_all()]
#' @seealso [abgd_summary()]
#' @return A data frame of the ABGD result.
abgd_delim <- function(output_folder, indices, label = NULL) {
    abgd.res <- abgd_read_result_all(output_folder)
    abgd.res.select <- abgd.res %>%
        list.filter(partition_type == "recursive") %>%
        list.filter(index %in% indices)
    if (is.null(label) & length(indices) > 1) {cli::cli_abort("When more than one partition results are evaluated, respective number of labels to differentiate each partition are required.")}
    if (is.null(label)) {
        label <- c("abgd")
    } else {
        label <- paste0("abgd_", label)
    }
    list.order <- list.select(abgd.res.select, index) %>% unlist()
    
    abgd.delim.res <- abgd.res.select %>%
        map2(label[match(list.order, indices)], \(res, lab) {
            tibble(
                taxa = res$result$taxa,
                !!lab := res$result$group
            )
        }) %>%
        reduce(~left_join(.x, .y, by = "taxa"))
    return(abgd.delim.res)
}

#' Perform GMYC on the Tree Generated by BEAST2
#' 
#' @param contree File path to the maximum credibility tree produced by BEAST2
#' @param prefix Character. Prefix to be added to the column name of the output data frame
#' @param method Character. Indicates whether to perform the single threshold (\code{single}) or the multiple threshold (\code{multiple}) GMYC, or \code{both}.
#' 
#' @importFrom ape is.ultrametric
#' @importFrom ape is.binary
#' @seealso [beast2_write_mrca()]
#' @export
#' @return A data frame of the delimitation result
gmyc_eval <- function(contree, prefix = NULL, method = "both") {
    contree <- ape::read.nexus(contree)
    if (!is.binary(contree) | !is.ultrametric(contree)) {cli::cli_abort("The tree is not binary or ultrametric.")}
    if (method == "both") {method <- c("single", "multiple")}
    res <- map(method, \(m){
        if (!is.null(prefix)) {
            col_name <- paste0(prefix, "_gmyc_", m) %>% gsub("^_", "", .)
        } else {
            col_name <- paste0("gmyc_", m) %>% gsub("^_", "", .)
        }
        splits::gmyc(contree, method = m) %>% splits::spec.list() %>% rename(identifier = sample_name, !!col_name := GMYC_spec)
    }) %>% 
        reduce(~left_join(.x, .y, by = "identifier")) %>%
        dplyr::select(identifier, everything())
    return(res)
}

#' Test for Best Fit Model in IQ-TREE
#' 
#' @description
#' This function calls local IQ-TREE CLI excutable to perform the modelFinder.
#' 
#' @param fasta_path Path to fasta file of the alignment
#' @param run_name Character. A custom name to be added to output files
#' @param output_folder Path to generate all output files from IQ-TREE
#' @param iqtree_exe Path to the IQ-TREE excutable
#' 
#' @return Path to the output folder
#' @export
test_sub_model <- function(
    fasta_path,
    run_name,
    output_folder,
    iqtree_exe = "/Users/hu_zhehao/Desktop/Biology-tools/iqtree-2.3.6-macOS-arm/bin/iqtree2"
) {
    output_folder <- file.path(output_folder, run_name)
    dir.create(output_folder, recursive = T)
    fasta_path_new <- file.path(output_folder, basename(fasta_path))
    file.copy(fasta_path, fasta_path_new, overwrite = T)
    cmd.cd <- paste0("cd ", output_folder)
    cmd.iqtree <- paste0(
        iqtree_exe,
        " -s '", 
        basename(fasta_path),
        "' -st DNA -m MF -pre '", 
        run_name,
        "' --redo"
    )
    cmd.file <- file.path(output_folder, paste0(run_name, "_modelTest.sh"))
    write(cmd.cd, cmd.file)
    write(cmd.iqtree, cmd.file, append = T)
    system(paste0("sh ", cmd.file))
    return(output_folder)
}

#' Exclude Sequences from Alignment by Filtering on Metadata
#' 
#' @description
#' Known issue: this function only filter the sequence based on single filtering criteria. For complex filter, take the source code and write for yourself. 
#' 
#' 
#' @param aln Alignment in DNAbin
#' @param field Character. Field (column) in metadata to apply the filter
#' @param pattern Character. Regex pattern to look for. Matching result will be excluded.
#' @param retain Numeric. Number of sequence to randomly exclude from exclusion (i.e. retain). Useful when excluding a group with large amount of sequence but keeping them represented in alignment. If bigger than number of sequences matching filtering criteria, no sequence will be excluded.
#' @param seed Will be passed to [set.seed]. Controls random selection of sequences when \code{retain} is not 0.
#' @param metadata A data frame. The metadata. Must incldue column given in \code{field}
#' @export
#' @return An alignment in DNAbin
#' @importFrom dplyr slice
exclude_by_meta <- function(aln, field, pattern, retain = 0, seed = NULL, metadata = db_pull("metadata.Sequence.NCBI", silent = T)) {
    exclude.meta <- metadata %>% 
        filter(str_detect(.data[[field]], pattern)) %>% 
        filter(`INSDSeq_primary-accession` %in% names(aln))
    if (retain != 0) {
        if (!is.null(seed)) {
            set.seed(seed)
        }
        if (nrow(exclude.meta) < retain) {
            avlb <- nrow(exclude.meta)
            cli::cli_alert_warning("Only {avlb} sequence matching exclusion criteria, but `retain` was set to {retain}, retaining all.")
            return(aln)
        }
        exclude.meta <- exclude.meta %>%
            slice(-sample(1:n(), retain))
    }
    include.index <- is.na(match(names(aln), exclude.meta$`INSDSeq_primary-accession`))
    aln_reduced <- aln[include.index]
    return(aln_reduced)
}

#' Split Concatenated Alignment Back into Single-locus Alignment
#' 
#' @description
#' Why should I do this? You might ask. Well, to get the actual sequences of each locus used in subsequent analysis to perform the model test.
#' 
#' @param alignment Concatenated alignment
#' @param loci_name Names of the loci
#' @param loci_length Lengths of each loci
#' @param output_folder Output folder
#' 
#' @seealso [concatenate_alignment()]
split_concat_alignment <- function(
    alignment,
    loci_name,
    loci_length,
    output_folder
) {
    partition <- data.frame(
        name = loci_name,
        length = loci_length
    ) %>%
        mutate(
            end = cumsum(length),
            start = lag(end, default = 0)+1
        )

    paths <- map_vec(1:length(loci_name), \(i){
        alignment.subset <- alignment[,partition[i,]$start:partition[i,]$end]

        # drop empty taxa
        char_matrix <- as.character(alignment.subset)
        all_missing <- map_lgl(1:nrow(char_matrix), \(i){
            all(char_matrix[i,] %in% c("-", "N", "n", "?"))
        })

        alignment.active <- alignment.subset[!all_missing, ]

        path <- file.path(output_folder, paste0("reextracted_", loci_name[[i]], ".fasta"))
        write.FASTA(alignment.active, path)
        return(path)
    })
    return(paths)
}

#' Concatenate Multiple Alignment into Multi-locus Alignment
#' 
#' Take individual alignment in fasta file, pair them based on a sequence mapping table, perform sequence filtering based on metadata ([exclude_by_meta()]), glue the paired sequences into concatenated sequence, construct concatenated alignment based on their sequence availability at certain loci.
#' 
#' @param fasta_path A character vector of path to the fasta file of each individual alignment. The order of the path indicates their order in the alignment.
#' @param names.alignment A character vector. Name of the loci. Must have the same length as \code{fasta_path}.
#' @param sequence.map A data frame to map sequences to individuals. Must have columns specified in \code{identifier} and \code{accs_col}
#' @param identifier Character. Name of the column in \code{sequence.map} corresponds to the organism's ID.
#' @param accs_col A character vector. Names of the column in \code{sequence.map} where accession numbers or sequence names of each organism are stored. Sequence names in the alignment of \code{fasta_path} will be matched against these.
#' @param include_strategy Either to only include taxa that has \code{all} sequences in the loci specified in \code{include_strategy_loci}, or taxa that has sequence in \code{any} of the loci.
#' @param include_strategy_loci Names of the loci used to check whether a taxa should be included or not.
#' @param output_folder A path to the output folder.
#' @param exclude_options A list of arguments passed to [exclude_by_meta()], except the first argument. Will used to exclude sequences. A list of the argument list can be provided to perform multiple steps of filtering.
#' 
#' @export
#' @importFrom purrr map_dfr
#' @importFrom dplyr distinct
#' @importFrom dplyr if_any
#' @importFrom dplyr if_all
#' @importFrom tidyselect all_of
#' 
#' @return Concatenated alignment will be written into a fasta file. The function returns a named list with information of the concatenated alignment: \describe{
#'   \item{\code{alignment}}{The acutal concatenated alignment object.}
#'   \item{\code{concat_fasta_path}}{Path to the concatenated alignment file}
#'   \item{\code{loci_fasta_path}}{Path to the fasta file of re-extracted alignment of each locus}
#'   \item{\code{taxa}}{Names of the taxa included in the alignment}
#'   \item{\code{loci_partition}}{A data frame. Storing the names of the loci and their length}
#' }
#' 
#' @seealso [exclude_by_meta()]
#' @seealso [split_concat_alignment()]
concatenate_alignment <- function(
    fasta_path,
    names.alignments,
    sequence.map,
    identifier,
    accs_col,
    include_strategy = "all",
    include_strategy_loci,
    output_folder,
    run_name = NULL,
    exclude_options = NULL
) {
    alignments.ls <- map(fasta_path, ~ape::read.FASTA(.x))

    if (pluck_depth(exclude_options) == 3) {
        alignments.ls <- reduce(
            exclude_options, 
            .init = alignments.ls, 
            .f = \(aln.ls, args){
                aln.ls %>%
                    map(\(aln){
                        rlang::exec(exclude_by_meta, aln = aln, !!!args)
                    })
            }
        )
    } else if (pluck_depth(exclude_options) == 2) {
        alignments.ls <- alignments.ls %>%
            map(\(aln){
                rlang::exec(exclude_by_meta, aln = aln, !!!exclude_options)
            })
    }

    bp.alignments <- map_vec(alignments.ls, \(x){length(x[[1]])})

    # extract sequence name/accession number
    names.seqs <- map(alignments.ls, ~names(.x))
    names(names.seqs) <- names.alignments

    # extarct sequence map of interest
    map.present <- map_dfr(1:length(names.alignments), \(i){
        sequence.map %>% 
            filter(.data[[accs_col[i]]] %in% names.seqs[[names.alignments[i]]])
    }) %>% distinct()

    # extract and pair sequences
    concat.seq.ls <- map(1:nrow(map.present), \(i){
        identifier <- map.present[[identifier]][i]
        
        seqpair <- map(1:length(names.alignments), \(gi){
            accNo <- map.present[[accs_col[gi]]][i]
            if (!is.na(accNo)) {
                if (!accNo %in% names.seqs[[gi]]) {
                    # accs No in seq map but not in alignment
                    seq <- ape::as.DNAbin(list(identifier = rep("N", bp.alignments[gi])))
                    names(seq) <- identifier
                    seqs.available <- F
                } else {
                    # accs No in seq map and in alignment
                    seq <- alignments.ls[[gi]][which(names.seqs[[gi]] == accNo)]
                    names(seq) <- identifier
                    seqs.available <- T
                }
            } else {
                # no accs No in seq map
                seq <- ape::as.DNAbin(list(identifier = rep("N", bp.alignments[gi])))
                    names(seq) <- identifier
                    seqs.available <- F
            }
            return(list(seq = seq, seqs.available = seqs.available))
        })

        seqs.available <- seqpair %>% map_vec(~.x$seqs.available)
        
        concat.seq <- seqpair %>% map(~.x$seq) %>% map(~as.matrix.DNAbin(.x)) %>% reduce(~cbind(.x, .y))
        return(list(concat.seq = concat.seq, seqs.available = seqs.available))
    })

    # glue sequence to concatenated alignment
    concat.alignment <- 
        do.call(rbind, rlist::list.filter(concat.seq.ls %>% map(~.x$concat.seq), !is.null(.)))

    # table for sequnce availability for later use
    concat.seq.comp <- 
        concat.seq.ls %>% 
        map_dfr(\(x){
            df <- data.frame(x$seqs.available) %>% t() %>% as.data.frame()
            colnames(df) <- names.alignments
            rownames(df) <- NULL
            return(df)
        }) %>% 
        bind_cols(c_organism_label = map.present$c_organism_label)


    if (include_strategy == "all") {
        # include strategy: only those with sequence available in specified columns 
        include.taxa <- concat.seq.comp %>% 
            filter(if_all(all_of(c(include_strategy_loci)), ~. == T))
    } else if (include_strategy == "any") {
        # include strategy: only those with sequence available in specified columns 
        include.taxa <- concat.seq.comp %>% 
            filter(if_any(all_of(c(include_strategy_loci)), ~. == T))
    } else {
        # include strategy: who's this
        cli::cli_abort("{.var include_strategy} must be either `all` or `any`, not {include_strategy}.")
    }
    
    # subset alignment
    concat.alignment <- concat.alignment[dimnames(concat.alignment)[[1]] %in% include.taxa$c_organism_label,]

    if (!is.null(run_name)) {
        concat.path <- file.path(output_folder, paste0("concat_", run_name, ".nexus"))
    } else {
        concat.path <- file.path(output_folder, "concat.nexus")
    }
    dir.create(concat.path, recursive = T)
    write.FASTA(concat.alignment, concat.path)

    # Re-extract alignment from each locus for model testing
    alnbyhand.path <- split_concat_alignment(concat.alignment, names.alignments, bp.alignments, output_folder)

    res <- list(
        alignment = concat.alignment,
        concat_fasta_path = concat.path,
        loci_fasta_path = alnbyhand.path,
        taxa = dimnames(concat.alignment)[[1]],
        loci_parition = data.frame(name = names.alignments, length = bp.alignments)
    )
}

