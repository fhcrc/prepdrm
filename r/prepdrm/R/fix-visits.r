#' Fix the 'visit' column in a data frame
#' @param df Data frame, with column "visit" normalized
#' @export
fix_visits <- function(df) {
  df$is_sc3 <- df$visit == 'SC+3'
  df$visit[df$is_sc3] <- 'SC+1'

  # Drop +
  df$visit[df$visit == 'SC+1'] <- 'SC1'
  df$visit <- tolower(df$visit)
  df
}
