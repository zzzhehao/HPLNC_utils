#' Snapshot the MrBayes Result During the Analysis
#' 
#' @description
#' This function allows a visualization of a consensus tree with currently sampled tree without stopping the MrBayes process.
#' 
#' @param assets_date Date of the asset folder where the MrBayes analysis process is located in. This is a project-specific data structure adaptation. It just simplifies the location of the target files. 
#' @param outgroup Sequence labels of the outgroup sequences. It will use automatic gap detection to find the sequences so it allows a certain degree of ambiguity in the name. 
#' @param subfolder Subfolder where the MrBayes process is located under the asset folder of the date. This is a project-specific data structure adaptation. It just simplifies the location of the target files. 
#' @param treeviz logical. Whether to visualize the tree. Default to TRUE.
#' @param xlim.factor Factor to expand the x-axis of the tree. Increase this value if long labels are omitted because they exceed the plotting area.
#' 
#' @return A ggplot of visualized tree if `treeviz` is set to TRUE. A tidytree "treedata" class object if `treeviz` is set to FALSE.
#' 
#' @import ggtree
#' @import treeio
#' @import ggtext
#' @import aplot
#' @import ggnewscale
#' @import grid
#' @import gridExtra
#' @import ggpubr
midwayTreeViz <- function(assets_date, outgroup, subfolder = "", treeviz = T, xlim.factor = 1.6) {
    suppressMessages({
        require(ggtree)
        require(treeio)
        require(ggtext)
        require(aplot) 
        require(ggnewscale)
        require(grid)
        require(gridExtra)
        require(ggpubr)
    })

    assets_path <- paste("docs/logs/assets", assets_date, subfolder, sep = "/") %>% gsub("/$", "", .)
    intercept.folder <- paste(assets_path, "preview_tree", sep = "/")
    suppressWarnings({dir.create(intercept.folder)})

    mrbayes.files <- list.files(assets_path, "run[0-9]{1}\\.p$|run[0-9]{1}\\.t$")
    nexus.file <- str_split(mrbayes.files[[1]], "\\.", simplify = T) %>% as.character() %>% .[str_detect(., "run[0-9]{1}|\\bp\\b|\\bt\\b", negate = T)] %>% paste(collapse = ".")
    intercept.paths <- paste(intercept.folder, c(mrbayes.files, nexus.file), sep = "/")

    file.copy(paste(assets_path, c(mrbayes.files, nexus.file), sep = "/"), intercept.paths, overwrite = T)

    nexus_nomcmc.file <- intercept.paths[[length(intercept.paths)]]

    # remove mcmc command
    lines <- readLines(nexus_nomcmc.file)
    mcmc_lines <- grep('^mcmc', lines)
    if (length(mcmc_lines) > 0) {
    lines <- lines[-mcmc_lines]
    }
    writeLines(lines, nexus_nomcmc.file)

    # add end statement for tree files
    tree.files <- grep("\\.t$", intercept.paths, value = T)
    walk(tree.files, \(file) {
        write("end;", file, append = T)
    })

    # Manual consensus tree in MrBayes
    cmd <- paste0("exe ", nexus.file)
    w <- readline(paste("Run command` ", cmd, " `at MrBayes console at", intercept.folder, "\n(Type anything to continue, `n` to abort.)\n >>> "))
    if (w == "n") {stop("Abort.")}


    # Tree Viz
    contree <- list.files(intercept.folder, pattern = "\\.con\\.tre$", full.names = T)
    tree <- treeio::read.mrbayes(contree)
    if (treeviz == F) {return(tree)}
        
    # Set node no. as numeric for later reroot operation
    tree@data$node <- as.numeric(tree@data$node)

    # Find exact name of outgroup (automatic gap detection)
    name.dist <- stringdist::stringdistmatrix(tree@phylo$tip.label, outgroup) %>% as.numeric() 
    name.dist.sorted <- name.dist %>% sort()
    weights <- 1/log(name.dist.sorted[-1]+1) # weight gap significance decreasingly while upper value of the gap increases
    maxGapIndex.sorted <- diff(name.dist.sorted)*weights %>% which.max()
    threshold <- mean(name.dist.sorted[c(maxGapIndex.sorted, maxGapIndex.sorted+1)]) # identify the threshold of the gap
    outgroup.match <- tree@phylo$tip.label[name.dist < threshold]

    cat("\nOutgroup given as:", outgroup, "\nMatched", outgroup.match, "\n\n")
    tree.rt <- ape::root(tree, outgroup.match, resolve.root = T, edgelabel = F)

    tree.rt@phylo$tip.label <- tree@phylo$tip.label

    tree.rt.sc <- rescale_tree(tree.rt, "length_mean")
    xmax <- max(c(tree.rt.sc@phylo$edge.length, tree@phylo$edge.length))
    tree <- list(tree = tree.rt.sc, xmax = xmax)

    tree.viz <- ggtree(tree[["tree"]], layout="rectangular") +
        geom_tiplab(
            parse = F,
            nudge_x = 0.003, 
            align = F,
            size = 3.8) + 
        geom_nodelab(
            aes(label = round(as.numeric(prob), 2)), 
            nudge_x = -0.01,
            nudge_y = 0.2, 
            hjust = 1, 
            size = 3.5) +
        geom_rootedge(rootedge = 0.02) +
        xlim(c(0, tree[["xmax"]]*xlim.factor)) +
        geom_treescale(x = 0, y = 0)
    return(tree.viz)
}

#' Visualize Consensus Tree Produced From MrBayes Analysis
#' 
#' @param assets_date Date of the asset folder where the MrBayes analysis process is located in. This is a project-specific data structure adaptation. It just simplifies the location of the target files. 
#' @param outgroup Sequence labels of the outgroup sequences. It will use automatic gap detection to find the sequences so it allows a certain degree of ambiguity in the name. 
#' @param subfolder Subfolder where the MrBayes process is located under the asset folder of the date. This is a project-specific data structure adaptation. It just simplifies the location of the target files. 
#' @param xlim.factor Factor to expand the x-axis of the tree. Increase this value if long labels are omitted because they exceed the plotting area.
#' 
#' @import ggtree
#' @import treeio
#' @import ggtext
#' @import aplot
#' @import ggnewscale
#' @import grid
#' @import gridExtra
#' @import ggpubr
#' 
#' @return A ggplot of visualized tree.
visualize_tree <- function(assets_date, outgroup, subfolder = "", xlim.factor = 1.5) {
    assets_path <- paste("docs/logs/assets", assets_date, subfolder, sep = "/") %>% gsub("/$", "", .)
    contree <- list.files(assets_path, pattern = "\\.con\\.tre$", full.names = T)
    if (length(contree) > 1) {stop("Multiple `.con.tre` files.")}
    tree <- treeio::read.mrbayes(contree)

    # Set node no. as numeric for later reroot operation
    tree@data$node <- as.numeric(tree@data$node)

    # Find exact name of outgroup (automatic gap detection)
    name.dist <- stringdist::stringdistmatrix(tree@phylo$tip.label, outgroup) %>% as.numeric() 
    name.dist.sorted <- name.dist %>% sort()
    weights <- 1/log(name.dist.sorted[-1]+1) # weight gap significance decreasingly while upper value of the gap increases
    maxGapIndex.sorted <- diff(name.dist.sorted)*weights %>% which.max()
    threshold <- mean(name.dist.sorted[c(maxGapIndex.sorted, maxGapIndex.sorted+1)]) # identify the threshold of the gap
    outgroup.match <- tree@phylo$tip.label[name.dist < threshold]

    cat("\nOutgroup given as:", outgroup, "\nMatched", outgroup.match, "\n\n")
    tree.rt <- ape::root(tree, outgroup.match, resolve.root = T, edgelabel = F)

# relableling

    tree.rt@phylo$tip.label <- tree@phylo$tip.label
    sequence.map <- DBpullTable("sequence.map", F, F)
    LUT.label <- sequence.map$c_organism_label
    names(LUT.label) <- sequence.map$COI

    # Fill VPS labels

    specmeta <- DBpullTable("metadata.Specimen.Haploniscidae")

    specmeta.label <- specmeta %>% 
        mutate(label = paste(
            voucher, gensp_morphology
        ))

    LUT.vpslabel <- specmeta.label$label
    names(LUT.vpslabel) <- specmeta.label$voucher

    LUT.label[str_detect(LUT.label, "^VPS")] <- LUT.vpslabel[LUT.label[str_detect(LUT.label, "^VPS")]]

    tree.rt@phylo$tip.label <- LUT.label[tree.rt@phylo$tip.label]

    # tree manipulation

    tree.rt.sc <- rescale_tree(tree.rt, "length_mean")
    xmax <- max(c(tree.rt.sc@phylo$edge.length, tree@phylo$edge.length))
    tree <- list(tree = tree.rt.sc, xmax = xmax)

    tree.viz <- ggtree(tree[["tree"]], layout="rectangular") +
        geom_tiplab(
            parse = F,
            nudge_x = 0.003, 
            align = F,
            size = 2.5) + 
        geom_nodelab(
            aes(label = round(as.numeric(prob), 2)), 
            nudge_x = -0.01,
            nudge_y = 0.2, 
            hjust = 1, 
            size = 3.5) +
        geom_rootedge(rootedge = 0.02) +
        xlim(c(0, tree[["xmax"]]*xlim.factor)) +
        geom_treescale(x = 0, y = 0)
    return(tree.viz)
}