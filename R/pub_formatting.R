#' Remove Descriptions on `ggplot` Plots for Publication
#' @import ggplot2
pub_nodescription <- function(ggobj){
    ggobj <- ggobj +
        theme(
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            plot.caption = element_blank()
        )
    return(ggobj)
}