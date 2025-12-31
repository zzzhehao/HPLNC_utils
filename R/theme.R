#' HPLNC Graphic Themes
#' @param ... Further argument to be passed to [ggplot2::theme()]
#' @name themes
NULL

#' @rdname themes
#' @export
theme_monochrome <- function(...){
    list(
        theme_minimal(),
        theme(
            panel.grid.minor.x = element_blank(),
            panel.border = element_rect(
				color = "black",
				fill = "transparent",
				linewidth = 1
			),
        ),
        theme(...)
    )
}

#' @rdname themes
#' @export
theme_monochrome_open <- function(...){
    list(
        theme_minimal(),
        theme(
            axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
            axis.line.y.left = element_line(color = "black", linewidth = 0.5),
            panel.grid.minor.x = element_blank()
        ),
        theme(...)
    )
}

#' @rdname themes
#' @export
theme_pub <- function(...){
    list(
        theme(
            ...
        ),
        labs(title = "", subtitle = "", caption = "")
    )
}

