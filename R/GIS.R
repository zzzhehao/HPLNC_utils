#' Convert DMS or DDM Coordinates to DD Coordinates. 
#' 
#' @author Zhehao Hu
#' @description
#' Also handles various notations and direction letter just as they are written in the cruise reports, which is never ready to use. 
#' @param dms DMS or DDM coordinates. This is a vectorized function, which means you can either pass one character that is latitude or longitude coordinates, or you can pass a vector of those. But not both latitude and longitude at the same time, it only handles one vector at a time. 
#' @return Converted decimal degree coordinates, ready to use. 
#' 
dms_to_decimal <- function(dms) {
    dms <- as.vector(dms)
    dms <- gsub(' ', '', dms)
    dms <- gsub('′', "'", dms)
    dms <- gsub('º', '°', dms)
    dms <- gsub('d', '°', dms)
    
    with_notation <- str_detect(dms, regex("N|S|W|E|-", ignore_case = T))
    neg <- ifelse(with_notation, str_detect(dms, regex("S|W|-", ignore_case = T)), FALSE)
    dms <- gsub(regex('[[:alpha:]]|-', ignore_case = T), "", dms)
  
    # Split the DMS string into components
    parts <- stringr::str_split(dms, "°|'", simplify = TRUE)
    
    # Extract degrees, minutes, and seconds
    degrees <- as.numeric(parts[,1])
    minutes <- as.numeric(parts[,2]) %>% 
        ifelse(is.na(.), 0, .)
    seconds <- as.numeric(parts[,3]) %>% 
        ifelse(is.na(.), 0, .)
    
    # Calculate decimal degrees
    decimal <- round(degrees + (minutes / 60) + (seconds / 3600), digits = 7)
    decimal <- ifelse(neg, -decimal, decimal)
    
    return(decimal)
}

#' Format DD Coordinates with Letters.
#' 
#' @author Zhehao Hu
#' @description
#' Format decimal degree coordinates notated with letters to clean numerical coordinates. West and south with minus.
#' @param dec DD coordinates. This is a vectorized function, which means you can either pass one character that is latitude or longitude coordinates, or you can pass a vector of those. But not both latitude and longitude at the same time, it only handles one vector at a time. 
#' 
dec_format <- function(dec) {
    dec <- as.vector(dec)
    neg <- grepl(regex("S|W", ignore_case = T), dec)
    dec <- gsub(regex('[[:alpha:]]', ignore_case = T), "", dec) %>% as.numeric()
    decimal <- ifelse(neg, -dec, dec)
}

#' Add checkerboard to a ggplot
#' 
#' @param scale_only Logical. If TRUE no checkerboard will be drawn, but the plot will be scaled (cropped or expaneded) to the area according to giving boundaries. 
#' @return Updated ggplot objects.
#' @export
geom_checkerboard <- function(min.lat, max.lat, lat.interval.count, min.long, max.long, long.interval.count, scale_only = F) {

    # Function modified from Jan Boyer: https://stackoverflow.com/a/54896719/24065650
    
    interval.lat <- (max.lat - min.lat) / lat.interval.count
    interval.long <- (max.long - min.long) / long.interval.count
    
    #define a constant (scaled to map size) width for grid rectangles
    #may need to adjust number to get width you prefer
    grid.width <- (max.lat - min.lat)/90
    
    is.even <- function(x) x %% 2 == 0 #function to test if number is even
    
    #dataframe of longitude rectangles
    rects.long <- data.frame(
        x_start = rep(seq(min.long, max.long - interval.long,by = interval.long), 2))
    rects.long$x_end <- rects.long$x_start + interval.long
    rects.long$y_start <- c(
        rep(min.lat, nrow(rects.long)/2),
        rep(max.lat - grid.width, nrow(rects.long)/2))
    rects.long$y_end <- c(
        rep(min.lat + grid.width, nrow(rects.long)/2),
        rep(max.lat, nrow(rects.long)/2))
    
    rects.long$color <- if(is.even(nrow(rects.long)/2)) { #even/odd test
        rep(c("black", "white"), nrow(rects.long)/2) #pattern for even
    } else {
        rep(c(rep(c("black", "white"), (nrow(rects.long) - 2)/4),"black"), 2)} #odd
    
    #dataframe of latitude rectangles
    rects.lat <- data.frame(y_start = rep(seq(min.lat, max.lat - interval.lat,by = interval.lat), 2))
    rects.lat$y_end <- rects.lat$y_start + interval.lat
    rects.lat$x_start <- c(
        rep(min.long, nrow(rects.lat)/2), 
        rep(max.long - grid.width, nrow(rects.lat)/2))
    rects.lat$x_end <- c(
        rep(min.long + grid.width, nrow(rects.lat)/2),
        rep(max.long, nrow(rects.lat)/2))
    
    rects.lat$color <- if(is.even(nrow(rects.lat)/2)) { #even/odd test
        rep(c("black", "white"), nrow(rects.lat)/2) #pattern for even
    } else {
        rep(c(rep(c("black", "white"), (nrow(rects.lat) - 2)/4),"black"), 2)} #odd
    
    #combine latitude and longitude grid
    rects.grid <- rbind(rects.lat, rects.long)
    
    #split into black dataframe and white dataframe
    rects.black <- rects.grid[rects.grid$color == "black",]
    rects.white <- rects.grid[rects.grid$color == "white",]
    
    #define axis breaks to match grid
    axis.breaks.x <- seq(min.long, max.long, interval.long)
    axis.breaks.y <- seq(min.lat, max.lat, interval.lat)
    
    lat_labeller <- function(x) {
        case_when(
            x > 0 ~ paste0(abs(x), "° N"),
            x < 0 ~ paste0(abs(x), "° S"),
            x == 0 ~ "0°"
        )
    }
    
    lon_labeller <- function(x) {
        case_when(
            x > 0 ~ paste0(abs(x), "° E"),
            x < 0 ~ paste0(abs(x), "° W"),
            x == 0 ~ "0°"
        )
    }

    if (scale_only) {
        return(list(
            scale_y_continuous(
                breaks = axis.breaks.y, 
                sec.axis = dup_axis(), 
                expand = c(0, 0), 
                labels = lat_labeller,
                limits = c(min.lat, max.lat)
            ),
            scale_x_continuous(
                breaks = axis.breaks.x, 
                sec.axis = dup_axis(), 
                expand = c(0, 0), 
                labels = lon_labeller,
                limits = c(min.long, max.long)
            ),
            theme_minimal(),
            theme(
                axis.title = element_blank()
            )))
    }
    
    list(
        scale_y_continuous(
            breaks = axis.breaks.y, 
            sec.axis = dup_axis(), 
            expand = c(0, 0), 
            labels = lat_labeller,
            limits = c(min.lat, max.lat)
        ),
        scale_x_continuous(
            breaks = axis.breaks.x, 
            sec.axis = dup_axis(), 
            expand = c(0, 0), 
            labels = lon_labeller,
            limits = c(min.long, max.long)
        ),
        geom_rect(
            data = rects.white, 
            inherit.aes = FALSE, 
            aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end), color = "black", fill = "white"),
        geom_rect(
            data = rects.black, 
            inherit.aes = FALSE,
            aes(xmin = x_start, xmax = x_end, ymin = y_start, ymax = y_end), color = "black", fill = "black"),
        theme_minimal(),
        theme(
            axis.title = element_blank()
        )
    ) %>% return()
}

#' Fetch bathymetry map from NOAA. 
#' 
#' @param cache Logical. Default to TRUE, the function will generate cache for the exact same boundary for future reuse. 
#' @import readr
fetch.NOAA.bathy <- function(boundaries, resolution = 5, cache = T) {
    require(marmap)
    require(scales)

    filename <- paste0("data/cache/bathy_", paste0(boundaries, resolution, collapse = ""), ".rds")

    if (filename %in% list.files("data/cache", full.names = T)) {
        bathy.df <- readRDS(filename)
    } else {
        bathy <- getNOAA.bathy(
            lat1 = boundaries[[1]], #min lat
            lon1 = boundaries[[2]], #min long
            lat2 = boundaries[[3]], #max lat
            lon2 = boundaries[[4]], #max long
            resolution = resolution
        )
        bathy.df <- fortify.bathy(bathy)
        if (cache) {
            write_rds(bathy.df, filename)
        }
    }
    return(bathy.df)
}

#' Fetch Land Elevation Data
#' 
#' @param boundaries A list of boundaries (min.lat, min.long, max.lat, max.long).
#' @param crs Coorinate reference system.
#' @param z Zoom level for elevation.
#' @param cache Logical. Default to TRUE, the function will generate cache for the exact same boundary for future reuse. 
#' 
#' @import elevatr
#' @import terra
#' @import sf
#' @import giscoR
#' @import scales
#' @import raster
#' @import dplyr
#' @import readr
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' 
#' @export
#' 
fetch.land.elv <- function(boundaries, crs = 4326, z = 3, cache = T) {
    filename <- paste0("data/cache/landelev_", paste0(boundaries, z, collapse = ""), ".rds")
    if (filename %in% list.files("data/cache", full.names = T)) {
        land_elev.df <- readRDS(filename)
    } else {
        land_sf <- giscoR::gisco_get_coastallines()
        bbox <- sf::st_bbox(
            c(xmin = boundaries[[2]], ymin = boundaries[[1]], xmax = boundaries[[4]], ymax = boundaries[[3]]), 
            crs = crs)
        
        land_sf_crop <- sf::st_crop(land_sf, bbox)
    
        land_transformed <- sf::st_cast(sf::st_transform(land_sf_crop, crs = crs), "MULTIPOLYGON")
    
        land_elev <- elevatr::get_elev_raster(
            locations = land_transformed, 
            z = z, 
            clip = "locations")
    
        land_elev.df <- as.data.frame(land_elev, xy = T) %>% na.omit() 
            
        colnames(land_elev.df)[3] <- "Elevation"
        land_elev.df <- land_elev.df[land_elev.df$Elevation >= 0,]

        if (cache) {
            write_rds(land_elev.df, filename)
        }
    }

    return(land_elev.df)
}

#' Create a Bathymetric and Topographic Basemap
#'
#' (Add your description here...)
#'
#' @param boundaries A list of boundaries (min.lat, min.long, max.lat, max.long). Default values are suitable for North Atlantic.
#' @param vertical.tiles Number of vertical checkerboard tiles. Will be ignored if `checkerboard` if False.
#' @param horizontal.tiles Number of horizontal checkerboard tiles. Will be ignored if `checkerboard` if False.
#' @param bathy.clr Logical. Add bathymetry colors?
#' @param land Logical. Add land elevation colors?
#' @param checkerboard Logical. Add checkerboard scales?
#' @param crs The coordinate reference system.
#' @param resolution Resolution for \code{getNOAA.bathy}.
#' @param z Zoom level for \code{elevatr::get_elev_raster}.
#'
#' @import ggplot2
#' @import marmap
#' @import ggnewscale
#' @import scales
#' @import sf
#'
#' @export
#'
bathy.basemap <- function(
    boundaries = list(0, -100, 75, 20), 
    vertical.tiles = 5,
    horizontal.tiles = 5,
    checkerboard = T,
    bathy.clr = T,
    land = T,
    crs = 4326,
    resolution = 5,
    z = 3,
    cache = T
) {
    require(ggplot2)

    filename <- paste0("data/cache/bathybasemap_", paste0(boundaries, z, as.numeric(bathy.clr), resolution, crs, collapse = ""), ".rds")

    if (filename %in% list.files("data/cache", full.names = T)) {
        # cat("Found cache map.\n")
        bathy <- readRDS(filename)
    } else {
        bathy <- ggplot() +
            geom_tile(
                data = fetch.NOAA.bathy(boundaries, resolution = resolution), aes(x, y, fill = z))

        if (bathy.clr) {
            bathy <- bathy +
                scale_fill_gradientn(
                    colors = c("#01040b", "#071e33", "#396573", "#598f92"),
                    values = scales::rescale(c(-8000, -5500, -2500, 0)),
                    breaks = c(-2000, -4000, -6000, -8000),
                    limits = c(-10000, 0), 
                    name = "Depth (m)",
                    na.value = "#2f3031"
                ) +
                ggnewscale::new_scale_fill()
        }

        if (land) {
            bathy <- bathy +
                geom_tile(
                    data = fetch.land.elv(boundaries, crs, z), 
                    aes(x = x, y = y, fill = Elevation), 
                    show.legend = F
                ) +
                scale_fill_gradientn(
                    colors = c("#2f3031", "#ad843d", "#ffd56a"),
                    values = scales::rescale(c(0, 1000, 2000)),
                    name = "Altitude (m)"
                ) +
                new_scale_fill()
        }
        bathy <- bathy +
            coord_sf(crs = st_crs(crs)) +
            labs(x = "", y = "")
        if (cache) {
            write_rds(bathy, filename)
        }
    }
    

    if (checkerboard) {
        bathy <- bathy +
        geom_checkerboard(
            boundaries[[1]], 
            boundaries[[3]], 
            vertical.tiles, 
            boundaries[[2]], 
            boundaries[[4]], 
            horizontal.tiles, 
            scale_only = F
        )
    }

    bathy <- bathy +
        theme(
            legend.position = "inside",
            legend.justification = c(0, 0),
            legend.position.inside = c(0.01, 0.01),
            legend.box = "vertical",
            legend.box.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
            legend.key.size = unit(0.5, 'cm'),
            legend.text = element_text(size = 8)
        )

    return(bathy)
}