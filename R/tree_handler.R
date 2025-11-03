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

visualize_tree <- function(assets_date, outgroup, subfolder = "", xfac = 1.5) {
    assets_path <- paste("docs/logs/assets", assets_date, subfolder, sep = "/") %>% gsub("/$", "", .)
    contree <- list.files(assets_path, pattern = "\\.con\\.tre$", full.names = T)
    if (length(contree) > 1) {stop("Multiple `.con.tre` files.")}
    tree <- treeio::read.mrbayes(contree)

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
        xlim(c(0, tree[["xmax"]]*xfac)) +
        geom_treescale(x = 0, y = 0)
    return(tree.viz)
}