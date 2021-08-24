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

remove_outlier_counts <- function(e, by = c("group", "sex"), sd = 2){
}






#' Check if there are enough mice per group
#' @description Check if there are enough mice per group across all regions
#' if the normalized counts data sets are split by specified grouping variables.
#' This function also automatically keeps only the common regions that are found across all comparison groups.
#'
#' @param e experiment object
#' @param by (str, default = c("group", "sex")) The mice attributes used to group the datasets into comparison groups.
#' @param min_n (int, default = 4) The minimum number of mice in each group for region comparisons.
#' @param remove (bool, TRUE) Remove any regions in the combined normalized count dataframes that don't have enough n to do a comparison on.
#' These regions are removed across all comparison groups.
#' @param log (bool, TRUE) Save the regions that don't have enough n into a .csv file in the output folder.
#' @return e experiment object
#' @export
#' @examples e <- enough_mice_per_group(e, by = c("group", "sex"), min_n = 4, remove = TRUE, log = TRUE)
#'
enough_mice_per_group<- function(e, by = c("group", "sex"), min_n = 4, remove = TRUE, log = TRUE){

  # Get the channels for an experiment
  e_info <- attr(e, "info")

  # correct mismatched attributes typed by users
  by <- match_m_attr(by)

  for (channel in e_info$channels){

    ## Get mouse counts
    mouse_count <- e$combined_normalized_counts[[channel]] %>%
      dplyr::select(all_of(c(by, "mouse_ID", "acronym"))) %>% dplyr::distinct() %>%
      dplyr::group_by(across(all_of(c(by, "acronym")))) %>% dplyr::count()


    # Create a new dataframe of just the regions and groups below the min thresh
    below_thresh <- dplyr::filter(mouse_count, n < min_n) %>%
      dplyr::bind_cols(tibble::tibble(channel = channel), .)


    if (exists('below_thresh_combine')){
      below_thresh_combine <- dplyr::bind_rows(below_thresh_combine, below_thresh)
    } else{
      below_thresh_combine <- below_thresh
    }

    # Remove regions that are below the threshold
    if (remove && dim(below_thresh)[1] > 0){
      to_remove <- which(e$combined_normalized_counts[[channel]]$acronym %in% below_thresh$acronym)
      e$combined_normalized_counts[[channel]] <- e$combined_normalized_counts[[channel]][-to_remove,]
      message("Removed regions for ", channel, " channel.")
    }

    ## Keep only common regions found across all groups
    common_regions <- get_common_regions(mouse_count, by)
    to_keep <- which(e$combined_normalized_counts[[channel]]$acronym %in% common_regions$acronym)
    e$combined_normalized_counts[[channel]] <- e$combined_normalized_counts[[channel]][to_keep,]

  }
     ## If there are any regions in any groups below the threshold, print them
  if (dim(below_thresh_combine)[1] > 0){
    message("There were regions below the minimum n:")
    print(below_thresh_combine, n = length(below_thresh_combine$group))

    if (log){
      out_path <- file.path(e_info$output_path, "regions_below_N_thresh.csv")
      write.csv(below_thresh_combine, file = out_path)
      message(paste0("Saved regions below the threshold at location: ", out_path))
    }
  } else{
    message("There were no regions below the minimum threshold.")
  }

  return(e)
}







#_______________ Internal functions _______________

#' Get common regions that are found across all the comparison groups.
#' @description If a single comparison group doesn't have any mice that have brain regions
#' contained in the other groups, no statistical comparison can be made. This region
#' should be removed.
#'
#' @param mouse_count
#' @param by
#'
#' @return unique acronyms that are the intersection across all groups
#'
#' @examples
get_common_regions <- function(mouse_count, by){
  # keeps track of number of group combinations
  combs <- 1
  for (attrib in by){
    combs <-combs*mouse_count[[attrib]] %>% unique() %>% length()
  }

  common_regions <-  mouse_count %>% dplyr::group_by(across(all_of(c("acronym")))) %>%
    dplyr::count() %>% dplyr::rename(group_n = n)

  common_regions <- common_regions[common_regions$group_n == combs, "acronym"]
  return(common_regions)
}






















