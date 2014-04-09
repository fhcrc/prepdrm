#' Below detection fill color
#' @export
bd_color <- '#AAAAAA'

#' Missing data fill color
#' @export
na_color <- '#BBBBBB'

#' Palette for filling cells
#' @export
fill_colors <- c(bd_color, bd_color, brewer.pal(8, 'Blues')[-(1:4)])

#' Percent mutant fill scale
#' @export
fill_scale <- scale_fill_manual("% mutant", values=fill_colors, na.value=na_color, drop=FALSE)

#' Text colors - dark and light, depending on fill
text_colors <- c('#333333', '#EEEEEE')

#' Text scale- dark and light, depending on fill
#' @export
text_scale <- scale_color_manual(values=text_colors, guide=FALSE, na.value='#EEEEEE')

#' Palette for filling cells colored by p-value
#' @export
p_value_fill <- scale_fill_manual('p-value\n(FDR corrected)',
                                  values=rev(c(bd_color, brewer.pal(4, 'OrRd')[-1])),
                                  drop=FALSE,
                                  na.value=na_color)

#' Categorize p values
#' @param p P values
#' @return factor
#' @export
categorize_p_values <- function(p) {
  cut(p,
      c(0, 0.001, 0.01, 0.05, 1),
      labels=c('< 0.001', '< 0.01', '< 0.05', '> 0.1'),
      include.lowest=TRUE)
}
