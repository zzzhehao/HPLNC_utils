#' Load Raw Bruker Flex Spectra
#' 
#' \code{maldi_load_spectra_raw} loads raw MALDI-TOF Spectra from Bruker Flex files or cache if possible.
#' 
#' @details
#' MALDI-TOF Raw Spectra are read internally by \code{MALDIquantForeign::importBrukerFlex(removeEmptySpectra = TRUE)}
#' 
#' @param raw_path Path to folder containing subfolder of each measurement.
#' @param cache_folder Path to output cache. Default to \code{NULL}. If set to \code{NULL}, it will try to use user configuration (see \code{create_config}), if user configuration is not set, it will create a 'cache/' folder under current working directory. 
#' @param return Logical. Whether to return the imported MassSpectrum object.
#' 
#' @importFrom MALDIquantForeign importBrukerFlex
#' 
#' 
#' @export
maldi_load_spectra_raw <- function(raw_path = NULL, cache_folder = NULL, return = T) {
    if (is.null(cache_folder)) {
        if (is.null(read_config("cache_folder"))) {
            cache_folder <- "cache"
        } else {
            cache_folder <- paste(read_config("cache_folder"), "maldi", sep = "/")
        }
    }
    suppressWarnings({dir.create(cache_folder)})
    rawSpectra <- importBrukerFlex(raw_path, removeEmptySpectra = T)
    saveRDS(rawSpectra, paste(cache_folder, "rawSpectra.rds", sep = "/"))
    if (return) {
        return(rawSpectra)
    } else {
        cli::cli_alert_success(paste0("Raw spectra are now cached under ' ", paste(cache_folder, "rawSpectra.rds", sep = "/"), "'."))
    }
}

#' Measurement Version Management 
#' 
#' @details 
#' \code{maldi_set_eval_list} relies on custom naming convention that named the MALDI-TOF measurments under the rule \code{<MALDI-TOF ID>_<attempt>} (e.g. VPSM001-1, VPSM003-3), which enables the function to effectively identify the attempts on the same sample and provide an interactive step for manually choosing attempts that are not the last attempt of the same samples to be evaluated. Default workflow always evaluates the last (non-empty) spectrum.
#' 
#' @param rawSpectra Raw Spectra. 
#' @param last_attempt Logical. Whether evaluate the last attempt of all samples. Interactive step will be introduced if set to \code{FALSE}.
#' @param sep Seperator used to split the measurement name into voucher and attempt. See details. 
#' 
#' @return A dataframe with four columns: \describe{
#' \item{\code{maldi_voucher}: sample voucher}
#' \item{\code{total_attemt}: Total measurements of the sample}
#' \item{\code{eval}: Which measurement to be evaluated}
#' \item{\code{run_name}: Measrument run name to be evaluated}
#' }
#' 
#' @importFrom MALDIquant metaData
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @export
maldi_set_eval_list <- function(rawSpectra, last_attempt = T, sep = "_") {
    # extract sample names (voucher + attempt run)
    MALDI.sample.names <- unlist(rawSpectra) %>% 
        purrr::map(., \(x) metaData(x)$sampleName) %>% 
        unlist() %>% 
        unique()
    # extract run info
        run.info <- data.frame(run_name = MALDI.sample.names) %>% 
            mutate(
                maldi_voucher = str_split_fixed(run_name, "_", 2)[,1],
                attempt = str_split_fixed(run_name, "_", 2)[,2]
            )
    # summarize total number of attempts for each sample
    run.info.summary <- run.info %>%
        group_by(maldi_voucher) %>%
        summarize(total_attempt = max(attempt))
    
    run.eval.temp <- run.info.summary %>%
        mutate(eval = total_attempt)
    if (!last_attempt) {
        input.path <- paste(data_path, "run.eval.csv", sep = "/")
        write.table(run.eval.temp, input.path, sep = ";", row.names = F)
        resp <- readline(paste("Please specify attempt for each sample to be evaluated at '", input.path, "', defaults are the last attempts of each sample. Type `y` to continue after saving your changes. Type anything else to abort."))
        if (resp != "y") {cli::cli_abort("Abort.")}
        run.eval <- read.table(paste(data_path, "run.eval.csv", sep = "/"), sep = ";", header = T) 
    } else {
        run.eval <- run.eval.temp
    }
    run.eval <- run.eval %>% 
        mutate(run_name = paste(maldi_voucher, eval, sep = "_"))
    return(run.eval)
}

#' Write Metadata into Spectra Object (burn the metadata in)
#' 
#' @param rawSpectra A MassSpectrum class object to be modified. Will be ignored if \code{use_cache} is set to \code{TRUE}.
#' @param use_cache Logical. Whether to load \code{rawSpectra} from the cache file. This could prevent from loading large object into your global environment. This option is only available when cache folder is specified by user configuration. See [create_config()].
#' @param run.eval A dataframe produced by [maldi_set_eval_list()]. If set to \code{NULL}, the function with default option will applied to \code{rawSpectra} to produce required dataframe.
#' @param metadata A dataframe of the metadata.
#' @param maldi_voucher_col Column name in \code{metadata} that contains the unique voucher used as MALDI voucher (as in [maldi_set_eval_list()])
#' @param metadata.extraction.list A list of one or two vectors of characters. The first vector should contain the columns name of the \code{metadata} that needs to be burnt into the MassSpectrum object. The second vector should contain the name of the metadata field burnt into the MassSpectrum object. Items from two vectors are paired by index. MALDI unique voucher specified in \code{maldi_voucher_col} should not be included and will be automatically added to the MassSpectrum object under the name \code{maldi_id}. If \code{inherit.name} is set to \code{TRUE} only one vector is required and the second vector will be ignored.
#' @param inherit.name Logical. Whether to inherit the name of the metadata in \code{metadata} in the MassSpectrum object.
#' @param write_cache Logical. Whether to write cache file.
#' @param cache_folder A path to cache folder. Default to \code{NULL} and the function will try to fetch the user configuration. If no path is set, the cache will be written to cache folder under current working directory.
#' @param write_metadata_field Logical. Whether to document the metadata fields written into the MassSpectrum object in this step. Must be \code{TRUE} if use [maldi_detect_peak()].
#' @param return Logical. Wether to return the modified object.
#' 
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom dplyr %>%
#' @importFrom purrr map_lgl
#' @importFrom MALDIquant metaData
#' @importFrom MALDIquant `metaData<-`
#' @importFrom purrr map
#' @importFrom stringr str_split_fixed
#' @importFrom purrr reduce2
#' @export
maldi_burn_metadata <- function(
    rawSpectra = NULL,
    use_cache = F,
    metadata,
    maldi_voucher_col,
    metadata.map,
    inherit.name = T,
    write_cache = T,
    run.eval = NULL,
    cache_folder = NULL,
    write_metadata_field = T,
    return = T
) {
    if (use_cache) {
        cache <- load_cache("rawSpectra\\.rds")
        if (!is.null(cache)) {
            cat("Loading raw spectra from the cache ... ")
            rawSpectra <- cache
        } else {
            if (is.null(rawSpectra)) {
                cli::cli_abort("Cache not found, and raw spectra are not provided.")
            }
        }
    }

    if (is.null(run.eval)) {
        run.eval <- maldi_set_eval_list(rawSpectra)
    }

    # check metadata mapping list
    if (!inherit.name & length(metadata.map)<2) {
        cli::cli_abort("metadata.map doesn't have enough vectors, and inherit.name is set to FALSE")
    }

    # combine metadata with the sample run name by specified column of unique maldi voucher
    referenceTable <- left_join(
        metadata %>% rename(
            "maldi_id" = maldi_voucher_col
        ), 
        run.eval %>% dplyr::select(c("maldi_voucher", "run_name")) %>% rename("maldi_id" = "maldi_voucher"),
        by = "maldi_id"
    )
    
    # create filter mask to drop those not to be evalutated
    filter_index <- map_lgl(rawSpectra, ~metaData(.x)$sampleName %in% run.eval$run_name)

    if (inherit.name) {
        metadata.map[[2]] <- metadata.map[[1]]
    }

    # add maldi id field to the head
    metadata.map <- list(c("maldi_id", metadata.map[[1]]), c("maldi_id", metadata.map[[2]]))

    # extract metadata and write into MassSpectrum object
    rawSpectra.rich <- map(rawSpectra, \(s){
        sampleName <- metaData(s)$sampleName
        maldi_voucher <- sampleName %>% str_split_fixed("_", 2) %>% .[1,1]
        idx <- which(maldi_voucher == referenceTable$maldi_id)
        if (length(idx) != 1) { 
            cli::cli_abort("Found multiple matches for voucher: ", maldi_voucher)
        }
    
        s <- reduce2(metadata.map[[1]], metadata.map[[2]], .init = s, .f = function(current_s, scr_col, dest_col) {
            value <- as.character(referenceTable[[scr_col]][[idx]])
            expr <- rlang::expr({
                metaData(current_s)[[!!dest_col]] <- !!value
                current_s
            })
            current_s <- rlang::eval_tidy(expr)
            return(current_s)
        }) 
        return(s)
    }) %>% .[filter_index]

    if (write_cache) {
        if (is.null(cache_folder)) {
            cache_folder <- read_config("cache_folder", T)
            if (is.null(cache_folder)) {
                cache_folder <- "cache"
            }
        }
        rr.path <- paste(cache_folder, "maldi/rawSpectra.rich.rds", sep = "/")
        saveRDS(rawSpectra.rich, file = rr.path)
        cli::cli_alert_success(paste0("\nRaw spectra with metadata are cached at '", rr.path, "'", sep = ""))

        rt.path <- paste(cache_folder, "maldi/referenceTable.rds", sep = "/")
        saveRDS(referenceTable, file = rt.path)
        cli::cli_alert_success(paste0("\nReference table is cached at '", rt.path, "'", sep = ""))
    }

    if (write_metadata_field) {
        if (is.null(cache_folder)) {
            cache_folder <- read_config("cache_folder", T)
            if (is.null(cache_folder)) {
                cache_folder <- "cache"
            }
        }
        bm.path <- paste(cache_folder, "maldi/burnt_metadata.rds", sep = "/")
        saveRDS(metadata.map[[2]], file = bm.path)
        cli::cli_alert_success(paste0("\nWritten metadata fields are cached at '", bm.path, "'", sep = ""))
    }

    if (return) {
        return(rawSpectra.rich)
    } else {
        if (!write_cache) {
            cli::cli_alert_danger("Both return and write_cache are set to FALSE. Nothing happens.")
        }
    }
}

#' Process Mass Spectra for Normalized, Baseline Removed Spectra Ready for Comparison
#' 
#' @param rawSpectra.rich A MassSpectrum class object with burnt-in metadata. Will be ignored if \code{use_cache} is set to \code{TRUE}.
#' @param use_cache Logical. Whether to load \code{rawSpectra.rich} from the cache file. This could prevent from loading large object into your global environment. This option is only available when cache folder is specified by user configuration. See [create_config()].
#' @param report Logical. Whether to produce pdf report of the spectra for each processing step.
#' @param report_folder A path to a folder to export the pdf report. If set to \code{NULL}, the function will try to fetch the path from user configuration. 
#' @param trim.min Numeric. Minimal mass to be included. Unit in Da.
#' @param trim.max Numeric. Maximal mass to be included. Unit in Da.
#' @param transform_method Character. Method for intensity transformation on the mass spectrum. Will be passed to [MALDIquant::transformIntensity()].
#' @param smooth_method Character. Method for smoothing mass spectrum. Will be passed to [MALDIquant::smoothIntensity()]
#' @param smooth_hws Numerical. Half window size applied to mass spectrum smoothing. Will be passed to [MALDIquant::smoothIntensity()]
#' @param baseline_method Character. Method for baseline removing of the mass spectrum. Will be passed to [MALDIquant::removeBaseline()]
#' @param SNIP_iteration Numerical. Iterations to be used in baseline removal with SNIP algorithm. Details see [MALDIquant::estimateBaseline()]. Will be passed to [MALDIquant::removeBaseline()]
#' @param calibrate_method Character. Method used for mass spectrum intensity calibration. Will be passed to [MALDIquant::calibrateIntensity()]
#' @param average_method Character. Method used for calculating average spectrum from each sample. Will be passed to [MALDIquant::averageMassSpectra()]
#' @param pause Logical. Whether to pause after every report is generated. This enables to stop the analysis in the middle to adjust parameters.
#' @param write_cache Logical. Whether to write cache file.
#' @param cache_folder A path to cache folder. Default to \code{NULL} and the function will try to fetch the user configuration. If no path is set, the cache will be written to cache folder under current working directory.
#' @param return Logical. Wether to return the modified object.
#' 
#' @importFrom purrr walk
#' @importFrom purrr map_chr
#' @importFrom MALDIquant trim
#' @importFrom MALDIquant transformIntensity
#' @importFrom MALDIquant smoothIntensity
#' @importFrom MALDIquant estimateBaseline
#' @importFrom MALDIquant removeBaseline
#' @importFrom MALDIquant calibrateIntensity
#' @importFrom MALDIquant averageMassSpectra
#' @importFrom MALDIquant metaData
#' @importFrom MALDIquant `metaData<-`
#' @importFrom cli cli_progress_update
#' @importFrom cli cli_progress_step
#' @importFrom cli cli_progress_done
#' 
#' @export
maldi_process_spectra <- function(
    rawSpectra.rich = NULL,
    use_cache = F,
    report = F,
    report_folder = NULL,
    trim.min = 2000,
    trim.max = 20000,
    transform_method = "sqrt",
    smooth_method = "SavitzkyGolay",
    smooth_hws = 10,
    baseline_method = "SNIP",
    SNIP_iteration = 20,
    calibrate_method = "TIC",
    average_method = "mean",
    pause = F,
    write_cache = T,
    cache_folder = NULL,
    return = T
) {
    if (use_cache) {
        cache <- load_cache("rawSpectra\\.rich\\.rds")
        if (!is.null(cache)) {
            cat("Loading raw spectra from the cache ... ")
            rawSpectra.rich <- cache
        } else {
            if (is.null(rawSpectra.rich)) {
                cli::cli_abort("Cache not found, and rich raw spectra (with metadata) are not provided.")
            }
        }
    }

    pause <- function() {
        if (pause) {
            resp <- readline("Spectrum processing paused. Type `n` to abort, type anything else to continue.")
            if (resp == "n") {
                cli::cli_abort("User aborted.")
            }
        }
    }
    
    if (is.null(report_folder)) {
        report_folder <- read_config("report_folder", T)
        if (is.null(report_folder)) {
            report_folder <- "report"
        }
        suppressWarnings({dir.create(report_folder)})
    }
    
    cli_progress_step("Loading raw spectra", paste0("Raw mass spectrum report has been exported at '", paste(report_folder, "rawSpectra.pdf", sep = "/"), "'."), spinner = TRUE)

    # generate report for raw spectra
    if (report) {
        suppressMessages({
            pdf(file = paste(report_folder, "1 Raw Spectra.pdf", sep = "/"), height = pdfHeight, width = pdfWidth) 
            walk(seq(along = rawSpectra.rich), \(i) {
                plot(rawSpectra.rich[[i]], main = names[i])
            })
            dev.off()
        })
        cli_progress_update()
    }
    
    # trim spectra
    cli_progress_step("Trimming spectra", spinner = TRUE)
    trimmedSpectra <- trim(rawSpectra.rich, c(trim.min, trim.max))
    # square
    cli_progress_step("Transforming spectra", spinner = TRUE)
    transformedSpectra <- transformIntensity(trimmedSpectra, method = transform_method)
    # smooth 
    cli_progress_step("Smoothing spectra", spinner = TRUE)
    smoothedSpectra <- smoothIntensity(transformedSpectra, method = smooth_method, halfWindowSize = smooth_hws)
    
    if (report) {
        cli_progress_step("Writing baseline estimation report", paste0("Smoothed spectrum with estimated baseline report has been exported at '", paste(report_folder, "2 Spectra_estimated_baseline.pdf", sep = "/")), spinner = TRUE)
        suppressMessages({
            pdf(file = paste(report_folder, "2 Spectra_estimated_baseline.pdf", sep = "/"), height = pdfHeight, width = pdfWidth)
            
            maxIntensity <- max(sapply(smoothedSpectra, function(x)max(intensity(x))))
    
            walk(seq(along = smoothedSpectra), \(i){
                baseline <- estimateBaseline(smoothedSpectra[[i]], method="SNIP", iterations=20)
                baselineth <- estimateBaseline(smoothedSpectra[[i]], method="TopHat")
                par(mar = c(7, 4, 4, 2) + 0.1)
                plot(smoothedSpectra[[i]], ylim=c(0, maxIntensity))
                lines(baseline, col="red")
                lines(baselineth, col="blue")
                mtext(text = "Red: baseline estimated by Statistics-sensitive Non-linear Iterative Peak-clipping algorithm (SNIP); Blue: baseline estimated by Tophat algorithm by van Herk 1996.", side = 1, line = 5, cex = 0.8, col = "black")
            })
    
            dev.off()
        })
        cli_progress_update()
        pause()
    }
    
    # Apply baseline removal
    cli_progress_step("Removing basline", spinner = TRUE)
    if (baseline_method == "SNIP") {
        baselineCorrectedSpectra <- removeBaseline(smoothedSpectra, method="SNIP", iterations = SNIP_iteration)
    } else {
        baselineCorrectedSpectra <- removeBaseline(smoothedSpectra, method = baseline_method)
    }
    
    if (report) {
        cli_progress_step("Writing baseline corrected spectra report", paste0("Baseline corrected spectra report has been exported at '", paste(report_folder, "3 Baseline_corrected_spectra.pdf", sep = "/")), spinner = TRUE)
        suppressMessages({
            pdf(file = paste(report_folder, "3 Baseline_corrected_spectra.pdf", sep = "/"), height = pdfHeight, width = pdfWidth)
    
            walk(seq(along = baselineCorrectedSpectra), \(i){
                plot(baselineCorrectedSpectra[[i]])
            })
    
            dev.off()
        })
        cli_progress_update()
        pause()
    }
    
    # normalize spectrum
    cli_progress_step("Normalizing spectra", spinner = TRUE)
    normalizedSpectra <- calibrateIntensity(baselineCorrectedSpectra, method = calibrate_method)
    
    
    if (report) {
        cli_progress_step("Writing normalized spectra report", paste0("Normalized spectra report has been exported at '", paste(report_folder, "4 Normalized_spectra.pdf", sep = "/")), spinner = TRUE)
        suppressMessages({
            pdf(file = paste(report_folder, "4 Normalized_spectra.pdf", sep = "/"), height = pdfHeight, width = pdfWidth)
    
            walk(seq(along = normalizedSpectra), \(i){
                plot(normalizedSpectra[[i]])
            })
    
            dev.off()
        })
        cli_progress_update()
        pause()
    }
    
    # average spectrum by sample name
    cli_progress_step("Averaging spectra", spinner = TRUE)
    samplenames <- map_chr(rawSpectra.rich, ~metaData(.x)$sampleName)
    averagedSpectra <- averageMassSpectra(normalizedSpectra, labels = samplenames, method = average_method)
    
    if (report) {
        cli_progress_step("Writing averaged spectra report", paste0("Averaged spectra report has been exported at '", paste(report_folder, "5 Average_spectra.pdf", sep = "/")), spinner = TRUE)
        suppressMessages({
            pdf(file = paste(report_folder, "5 Average_spectra.pdf", sep = "/"), height = pdfHeight, width = pdfWidth)
    
            walk(seq(along = averagedSpectra), \(i){
                plot(averagedSpectra[[i]])
            })
    
            dev.off()
        })
        cli_progress_update()
        pause()
    }
    
    if (write_cache) {
        cli_progress_step("Writing cache file", spinner = TRUE)
        if (is.null(cache_folder)) {
            cache_folder <- read_config("cache_folder", T)
            if (is.null(cache_folder)) {
                cache_folder <- "cache"
            }
        }
        as.path <- paste(cache_folder, "maldi/averagedSpectra.rds", sep = "/")
        saveRDS(averagedSpectra, file = as.path)
        cli::cli_alert_success(paste0("\nAveraged spectra are cached at '", as.path, "'", sep = ""))
    }
    
    cli_progress_done()

    if (return) {
        return(averagedSpectra)
    } else {
        if (!write_cache) {
            cli::cli_alert_danger("Both return and write_cache are set to FALSE. Nothing happens.")
        }
    }
}

#' Detect Peaks in Mass Spectra and Convert into Peak Matrix
#' 
#' 
#' @importFrom cli cli_progress_step
#' @importFrom cli cli_progress_done
#' @importFrom cli cli_progress_update
#' @importFrom cli cli_progress_bar
#' @importFrom cli cli_abort
#' @importFrom MALDIquant detectPeaks
#' @importFrom MALDIquant binPeaks
#' @importFrom MALDIquant intensityMatrix
#' @importFrom MALDIquant metaData
#' @importFrom MALDIquant `metaData<-`
#' @importFrom purrr map_dfr
#' @importFrom dplyr %>%
#' @importFrom dplyr all_of
#' @importFrom tidyr unite
maldi_detect_peak <- function(
    averagedSpectra = NULL,
    use_cache = F,
    SNR = 8,
    hws = 8,
    method = "MAD",
    tolerance = 0.002,
    bin_method = "strict", 
    min_peak = 25,
    metadata_field = NULL,
    lab_field = NULL,
    transform_method = "hellinger",
    write_cache = T,
    cache_folder = NULL,
    return = T
) {
    if (use_cache) {
        cli_progress_step("Loading averaged spectra from the cache ...", "Cached averaged spectra loaded.", spinner = T)
        cache <- load_cache("averagedSpectra\\.rds")
        if (!is.null(cache)) {
            averagedSpectra <- cache
            cli_progress_done()
        } else {
            if (is.null(averagedSpectra)) {
                stop("Cache not found, and averaged spectra are not provided.")
            }
        }
    }

    peaks <- detectPeaks(averagedSpectra, halfWindowSize = hws, SNR = SNR, method = method)

    binned_peaks <- peaks
    binning_progress <- vector()
    number_peaks_progress <- 0
    cli_progress_bar("Binning peaks recursively. Iterations:")
    repeat { 
        binned_peaks <-  binPeaks(binned_peaks, tolerance=tolerance, method=bin_method)
        binned_peak_matrix <- as.matrix(intensityMatrix(binned_peaks))
        number_peaks <- ncol(binned_peak_matrix)
        if (number_peaks==number_peaks_progress) {
            break
        }
        number_peaks_progress <- number_peaks
        binning_progress <- append(binning_progress,number_peaks)
        cli_progress_update()
    }
    cli_progress_done()

    peakmatrix <- intensityMatrix(binned_peaks)
    peakmatrix[is.na(peakmatrix)] <- 0 

    if (is.null(metadata_field)) {
        cache <- load_cache("burnt_metadata\\.rds")
        if (is.null(cache)) {
            cli_abort("Cached burnt metadata fields not found. Please check cache files. Metadata has to be cached by `maldi_burn_metadata(write_metadata_field = TRUE)`")
        } else {
            metadata_field <- cache
        }
    }

    if (is.null(lab_field)) {
        lab_field <- metadata_field
    }

    binned_metadata <- binned_peaks %>%
        # extract metadata
        map_dfr(function(x) {
            md <- metaData(x)
            md[metadata_field]
        }) %>%
        # glue to label
        unite("label", all_of(lab_field), sep = "_", remove = F, na.rm = FALSE)

    # quality control: retain only samples with number of peaks exceeding set threshold
    peak_count <- rowSums(peakmatrix != 0)
    idx <- which(peak_count > min_peak)
    peakmatrix <- peakmatrix[idx,]
    binned_metadata <- binned_metadata[idx,]


    row.names(peakmatrix) <- binned_metadata$label

    if (!is.null(transform_method)) {
        peakmatrix <- vegan::decostand(peakmatrix, transform_method)
    }

    res <- list(
        matrix = peakmatrix, 
        metadata = binned_metadata,
        transformation = transform_method
    )

    if (write_cache) {
        cli_progress_step("Writing cache file", spinner = TRUE)
        if (is.null(cache_folder)) {
            cache_folder <- read_config("cache_folder", T)
            if (is.null(cache_folder)) {
                cache_folder <- "cache"
            }
        }
        res.path <- paste(cache_folder, "maldi/peakmatrix.rds", sep = "/")
        saveRDS(res, file = res.path)
        cli::cli_alert_success(paste0("\nPeak matrix and metadata are cached at '", res.path, "'", sep = ""))
    }
    
    cli_progress_done()

    if (return) {
        return(res)
    } else {
        if (!write_cache) {
            cli::cli_alert_danger("Both return and write_cache are set to FALSE. Nothing happens.")
        }
    }
}

