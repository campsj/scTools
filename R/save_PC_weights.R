#' Save Principal component weights
#'
#' Save heatmap in designated folder
#' @param df Dataframe with genes as rownames and PC weights as columns
#' @param PC Principal component
#' @return Save csv files of negative and positive gene weights per principal component
#' @export

save_PC_weights <- function(df, PC) {
  df %>%
    rownames_to_column("gene") %>%
    arrange_(PC) %>%
    select_("gene", PC) %>%
    filter_(interp(~ X > 0, X = as.name(PC))) %>%
    mutate_(PC_weight = interp(~ abs(X), X = as.name(PC))) %>%
    select_("gene", "PC_weight") %>%
    column_to_rownames("gene") %>%
    write.csv(paste(PC, "positive_weights.csv", sep = "_"))
  df %>%
    rownames_to_column("gene") %>%
    arrange_(PC) %>%
    select_("gene", PC) %>%
    filter_(interp(~ X < 0, X = as.name(PC))) %>%
    mutate_(PC_weight = interp(~ abs(X), X = as.name(PC))) %>%
    select_("gene", "PC_weight") %>%
    column_to_rownames("gene") %>%
    write.csv(paste(PC, "negative_weights.csv", sep = "_"))
}
