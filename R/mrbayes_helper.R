#' @title Helper Functions to Generate MrBayes Block in Nexus
#' 
#' @name mrbayes
#' 
#' @description
#' These functions help generate MrBayes command to put in nexus block that can be written with the alignment together. \code{mrbayes_mcmc_cmd}, \code{mrbayes_partition_cmd}, \code{mrbayes_lset_cmd}, \code{mrbayes_prset_cmd}, \code{mrbayes_sum_cmd} produce individual part of the command. Their output can be passed into \code{mrbayes_block} to package them as a complete nexus block. Also see [prepare_mrbayes_nexus()] for a more powerful wrapper.
#' 
#' @return A character string. 
#' @seealso [prepare_mrbayes_nexus()]
NULL

#' @rdname mrbayes
#' @title Generate MCMC command for MrBayes Block. 
#' 
#' @param gen Number of MCMC generations
#' @param nchains Number of MCMC chains
#' @param nruns Number of parallel runs
#' @param stopval Threshold value to stop the run
#' 
#' @importFrom purrr map_vec
#' @importFrom dplyr lag
#' 
#' @export
mrbayes_mcmc_cmd <- function(
    gen = 10000000,
    nchains = 4,
    nruns = 4,
    stopval = 0.01,
    ...
) {
    # default list
    defaults <- list(
        ngen = as.integer(gen), 
        nchains = as.integer(nchains),
        nruns = as.integer(nruns),
        relburnin = "yes",
        samplefreq = 1000,
        printfreq = 5000,
        diagnfreq = 10000,
        checkfreq = 10000,
        savebrlens = "yes",
        starttree = "random", 
        stoprule = "yes", 
        stopval = stopval
    )

    pars <- modifyList(defaults, list(...))

    pars_cmd <- map_vec(1:length(pars), \(i){paste0(names(pars)[i], "=", as.character(pars[[i]]))}) %>% paste(collapse = " ")

    cmd <- paste0("mcmc ", pars_cmd, ";")
    return(cmd)
}

#' @rdname mrbayes
#' @export
#' @param loci_length A numeric vector. Length of each locus.
#' @param loci_name A character vector. Name of each locus.
#' @param protein_coding_loci A character vector of the names of the protein coding locus. The locus will be partitioned further based on codon position to allow different substitution rate at each position.
mrbayes_partition_cmd <- function(
    loci_length,
    loci_name,
    protein_coding_loci
) {
    # Checks
    if (length(loci_length) != length(loci_name)) {cli::cli_abort("{.var loci_length} must have the same length as {.var loci_name}")}
    if (any(!protein_coding_loci %in% loci_name)) {cli::cli_abort("{.var protein_coding_loci} must be elements in {.var loci_name}")}

    partition.scheme <- data.frame(
        name = loci_name,
        length = loci_length
    ) %>%
        mutate(
            protein_coding = ifelse(name %in% protein_coding_loci, T, F),
            end = cumsum(length),
            start = lag(end, default = 0) + 1
        )

    partition <- map(1:nrow(partition.scheme), \(i){
        if (partition.scheme[i,]$protein_coding) {
            partition.names <- paste0(partition.scheme[i,]$name, c("_pos1", "_pos2", "_pos3"))
            p1 <- paste0("charset ", partition.names[1], " = ", partition.scheme[i,]$start, "-", partition.scheme[i,]$end, ".\\3;")
            p2 <- paste0("charset ", partition.names[2], " = ", partition.scheme[i,]$start+1, "-", partition.scheme[i,]$end, "-.\\3;")
            p3 <- paste0("charset ", partition.names[3], " = ", partition.scheme[i,]$start+2, "-", partition.scheme[i,]$end, ".\\3;")
            return(list(name = partition.names, charset = c(p1, p2, p3)))
        } else {
            partition.names <- partition.scheme[i,]$name
            p <- paste0("charset ", partition.scheme[i,]$name, " = ", partition.scheme[i,]$start, "-", partition.scheme[i,]$end, ";")
            return(list(name = partition.names, charset = c(p)))
        }
    })

    charset.lines <- map(partition, ~.x$charset) %>% unlist() %>% paste(collapse = "\n")
    partition.names <- map(partition, ~.x$name) %>% unlist()

    part.def <- paste0(
        "partition conc = ", 
        length(map(partition, ~.x$charset) %>% unlist()),
        ": ",
        paste(partition.names, collapse = ", "),
        ";"
    )

    cmd <- paste(
        charset.lines,
        part.def,
        "set partition = conc;",
        sep = "\n"
    )

    return(cmd)
}

#' @rdname mrbayes
#' @export
#' @param nst Value passed to \code{lset nst}
#' @param rates Value passed to \code{lset rates}
#' @param applyto A numeric vector, indicating the index of the partition to apply the \code{lset} settings to. Default to \code{NULL}, which applies to the whole alignment.
mrbayes_lset_cmd <- function(
    nst,
    rates,
    ...,
    applyto = NULL
) {
    if (is.null(applyto)) {
        applyto <- ""
    } else {
        applyto <- paste0(" applyto=(", paste(applyto, collapse = ","), ")")
    }

    pars <- list(...)
    if (length(pars)==0) {
       pars_cmd <- ""
    } else {
        pars_cmd <- map_vec(1:length(pars), \(i){paste0(names(pars)[i], "=", as.character(pars[[i]]))}) %>% paste(collapse = " ")
    }

    cmd <- paste0("lset", applyto, " nst=", nst, " rates=", rates, " ", pars_cmd, ";")
    return(cmd)
}

#' @rdname mrbayes
#' @export
mrbayes_prset_cmd <- function(
    ...,
    applyto = NULL
) {
    if (is.null(applyto)) {
        applyto <- ""
    } else {
        applyto <- paste0(" applyto=(", paste(applyto, collapse = ","), ")")
    }

    pars <- list(...)
    if (length(pars)==0) {
       return(invisible(""))
    } else {
        pars_cmd <- map_vec(1:length(pars), \(i){paste0(names(pars)[i], "=", as.character(pars[[i]]))}) %>% paste(collapse = " ")
    }

    cmd <- paste0("prset", applyto, " ", pars_cmd, ";")
    return(cmd)
}

#' @rdname mrbayes
#' @export
#' @param cmd Either \code{"sump"} or \code{"sumt"}, indicating the generating command. Will be ignored when \code{...} is missing.
mrbayes_sum_cmd <- function(
    ...,
    cmd
) {
    pars <- list(...)
    if (length(pars)==0) {
        return("sump;\nsumt;")
    }
    pars_cmd <- map_vec(1:length(pars), \(i){paste0(names(pars)[i], "=", as.character(pars[[i]]))}) %>% paste(collapse = " ")

    sum_cmd <- paste0(cmd, " ", pars_cmd, ";")
    return(sum_cmd)
}

#' @rdname mrbayes
#' @export
mrbayes_block <- function(
    ...
) {
    cmds <- c(...)
    paste(
        "BEGIN MrBayes;",
        "set autoclose=yes;",
        paste(cmds, collapse = "\n\n"),
        "END;",
        sep = "\n\n"
    )
}

#' @name mrbayes-args
#' @rdname mrbayes
#' @param ... Additional arguments. \describe{
#'   \item{[mrbayes_mcmc_cmd]}{Further parameters to put in mcmc command}
#'   \item{[mrbayes_lset_cmd]}{Further parameters to put in lset command}
#'   \item{[mrbayes_prset_cmd]}{Further parameters to put in prset command}
#'   \item{[mrbayes_sum_cmd]}{Further parameters to put in sumt or sump command, when given, \code{cmd} must also be provided to specify whether sumt or sump is being generated}
#'   \item{[mrbayes_block]}{Commands to be included in the MrBayes block}
#' }
NULL

#' Helper-Wrapper for Generating Nexus File for MrBayes
#' 
#' @description
#' Takes arguments for MrBayes helper construct nexus file with mrbayes block, and export. 
#' 
#' @seealso [mrbayes]
#' 
#' @export
#' @importFrom purrr pluck_depth
#' @importFrom purrr map_vec
#' @importFrom stringr str_detect
#' @importFrom ape write.nexus.data
#' 
#' @param alignment Alignment in DNAbin
#' @param output_path Path of the output nexus file
#' @param mcmc_options A named list of arguments to be passed to [mrbayes_mcmc_cmd()]
#' @param partition logical. If a partition is to be applied to alignment. 
#' @param loci_length A numeric vector. Length of each locus. Will be passed to [mrbayes_partition_cmd()]
#' @param loci_name A character vector. Name of each locus. Will be passed to [mrbayes_partition_cmd()]
#' @param protein_coding_loci A character vector of the names of the protein coding locus. The locus will be partitioned further based on codon position to allow different substitution rate at each position. Will be passed to [mrbayes_partition_cmd()]
#' @param lset_options A named list of arguments to be passed to [mrbayes_lset_cmd()]
#' @param prset_option A named list of arguments to be passed to [mrbayes_prset_cmd()]
#' @param sum_options A named list of arguments to be passed to [mrbayes_prset_cmd()]
#' 
#' @return The file path of written nexus file. 
prepare_mrbayes_nexus <- function(
    alignment,
    output_path,
    mcmc_options,
    partition = T,
    loci_length,
    loci_name,
    protein_coding_loci,
    lset_options,
    prset_options = NULL,
    sum_options = NULL
){
    mcmc <- do.call(mrbayes_mcmc_cmd, mcmc_options)
    if (partition) {
        partition_arg <- mrbayes_partition_cmd(loci_length, loci_name, protein_coding_loci)
    } else {
        partition_arg <- ""
    }

    if (pluck_depth(lset_options) == 3) {
        lset <- map_vec(lset_options, ~do.call(mrbayes_lset_cmd, .x))
    } else if (pluck_depth(lset_options) == 2) {
        lset <- do.call(mrbayes_lset_cmd, lset_options)
    }

    if (pluck_depth(prset_options) == 3) {
        prset <- map_vec(prset_options, ~do.call(mrbayes_prset_cmd, .x))
    } else if (pluck_depth(prset_options) == 2) {
        prset <- do.call(mrbayes_prset_cmd, prset_options)
    } else if (is.null(prset_options)) {
        prset <- ""
    }

    if (pluck_depth(sum_options) == 3) {
        sum_cmd <- map_vec(sum_options, ~do.call(mrbayes_sum_cmd, .x))
    } else if (pluck_depth(sum_options) == 2) {
        sum_cmd <- do.call(mrbayes_sum_cmd, sum_options)
    } else if (is.null(sum_options)) {
        sum_cmd <- mrbayes_sum_cmd()
    }

    arg_ls <- list(
        mcmc,
        lset,
        prset,
        sum_cmd
    )

    if (partition) {
        arg_ls <- c(partition_arg, arg_ls)
    }

    mrbayes_block <- do.call(mrbayes_block, arg_ls)

    dir.create(dirname(output_path))
    file_output_path <- ifelse(str_detect(output_path, "\\.nexus$"), output_path, paste0(output_path, ".nexus"))
    write.nexus.data(alignment, file_output_path, interleaved = F)
    write(mrbayes_block, file_output_path, append = T)
    return(file_output_path)
}
