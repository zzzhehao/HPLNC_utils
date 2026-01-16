#' Convert DMS or DDM Coordinates to DD Coordinates. 
#' 
#' @author Zhehao Hu
#' @description
#' Also handles various notations and direction letter just as they are written in the cruise reports, which is never ready to use. 
#' @param dms DMS or DDM coordinates. This is a vectorized function, which means you can either pass one character that is latitude or longitude coordinates, or you can pass a vector of those. But not both latitude and longitude at the same time, it only handles one vector at a time. 
#' @return Converted decimal degree coordinates, ready to use. 
#' 
gis_convert_dms2dd <- function(dms) {
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
gis_format_DD <- function(dec) {
    dec <- as.vector(dec)
    neg <- grepl(regex("S|W", ignore_case = T), dec)
    dec <- gsub(regex('[[:alpha:]]', ignore_case = T), "", dec) %>% as.numeric()
    decimal <- ifelse(neg, -dec, dec)
}

