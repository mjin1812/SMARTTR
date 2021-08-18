# Quality check functions


#_____________________ For mouse objects _________________________


#' Detect atlas regions that only show up in a single slice object
#' @description Quality check function to make sure that the regions included for analysis show up
#' in more than 1 slice object, otherwise the user can remove them from the mouse object.
#'
#' The user should run [normalize_cell_counts()] and [get_cell_table()] functions prior to using this function.
#' If the user has run [split_hipp_DV()] the function will automatically account for this
#'
#' @param m
#' @param remove
#' @param output
#'
#' @return
#' @export
#'
#' @examples
detect_single_slice_regions <- function(m, remove = FALSE, output = TRUE){

  # $hipp_DV_normalized_counts
  if (!is.null(m$hipp_DV_normalized_counts)){

    # Do the same thing for hippocampus and append to the log files

  }
}



#_____________________ For experiment objects _________________________

## remove any normalized regions counts that are more than 2 SDs (user-defined) higher
# than their cohort mean

remove_outlier_counts <- function(sd = 2){


}
