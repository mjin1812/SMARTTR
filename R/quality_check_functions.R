# Quality check functions


#_____________________ For mouse objects _________________________

#' Detect atlas regions that only show up in a single slice object within a mouse.
#' @description Quality check function to make sure that the regions included for analysis show up
#' in more than 1 slice object, otherwise the user can remove exceptions from the mouse object.
#'
#' Regions counts derived from only one image may be less accurate. This function can generate a log of these regions
#' so the user can qualitatively evaluate the raw data. Users also have to option of removing these regions automatically
#' from normalized_counts dataframe.
#'
#'
#' The user should run [normalize_cell_counts()] and [get_cell_table()] functions prior to using this function.
#' If the user has run [split_hipp_DV()] with the option to merge, this function will account for for the dorsal and
#' ventral hippocampal counts separately.
#'
#' @param m
#' @param remove (bool, FALSE) Remove any regions in the normalized counts table that
#' @param log (bool, TRUE) Save the regions that don't have enough n into a .csv file in the output folder.
#'
#' @return
#' @export
#'
#' @examples
detect_single_slice_regions <- function(m, remove = FALSE, log = TRUE){
  # $hipp_DV_normalized_counts
  if (!is.null(m$hipp_DV_normalized_counts)){
    # Do the same thing for hippocampus and append to the log files
  }
}

#_____________________ For experiment objects _________________________

#' Detect, log, and remove outlier counts. This function
#' removes any normalized regions counts that are more than `n_sd` standard deviations (default = 2) higher
#' than their cohort mean.
#'
#' @param e experiment object
#' @param by (str, default = c("group", "sex")) The mice attributes used to group the datasets into comparison groups.
#' @param n_sd (int, default = 2). Number of standards deviations above and below which categorizes outliers.
#' @param remove (bool, default = FALSE) Remove all outlier rows from the combined normalized counts dataframe in the experiment object.
#' @param log (bool, default = TRUE) Save the logged outlier regions into a csv file in the output folder.
#'
#' @return e experiment object. Outlier counts in the experiment object are removed if remove = TRUE.
#' @export
#' @examples e <- find_outlier_counts(e, by = c("group","sex"), n_sd = 2, remove = FALSE, log = TRUE)
#'
find_outlier_counts <- function(e, by = c("group", "sex"), n_sd = 2, remove = FALSE, log = TRUE){

  # Get the channels for an experiment
  e_info <- attr(e, "info")

  # correct mismatched attributes typed by users
  by <- match_m_attr(by)

  for (channel in e_info$channels){

    # Get the mean and sd of each group
    stats_df <- e$combined_normalized_counts[[channel]] %>%
      dplyr::group_by(across(all_of(c(by, "acronym")))) %>%
      dplyr::summarise(mean_norm_counts = mean(normalized.count.by.volume),
                       sd_norm_counts = sd(normalized.count.by.volume))

    joined_df <- e$combined_normalized_counts[[channel]] %>% dplyr::inner_join(stats_df, c(by, "acronym")) %>%
      dplyr::mutate(outlier = ifelse(normalized.count.by.volume > mean_norm_counts + sd_norm_counts*n_sd | normalized.count.by.volume < mean_norm_counts - sd_norm_counts*n_sd, TRUE, FALSE))

    # Create a list of outliers to print
    if (isTRUE(log)){
      log_df <- joined_df %>% dplyr::filter(outlier) %>%
        dplyr::select(-c(mean_norm_counts,sd_norm_counts)) %>%
        dplyr::arrange(across(all_of(by)), acronym, mouse_ID)

      if (length(log_df$mouse_ID) > 0 ){
        message("There were ", length(log_df$mouse_ID),
                " outliers found. Outliers were based on within group mean and standard deviation.")
        out_path <- file.path(e_info$output_path, paste0("region_count_outliers_", channel,".csv"))
        write.csv(log_df, file = out_path)
        message(paste0("Saved regions outliers dataframe at:\n", out_path))
      } else {
        message("No regions outliers were found within each group for channel ", channel)
      }
    }

    if (isTRUE(remove)){
      e$combined_normalized_counts[[channel]] <- joined_df %>% dplyr::filter(!outlier) %>%
        dplyr::select(-c(mean_norm_counts,sd_norm_counts, outlier))
    }
  }
  return(e)
}


#' Check if there are enough mice per analysis subgroup across all regions.
#' if the normalized counts data sets are split by specified grouping variables.
#' This function also automatically keeps only the common regions that are found across all comparison groups.
#'
#' @param e experiment object
#' @param by (str, default = c("group", "sex")) The mice attributes used to group the datasets into comparison groups.
#' @param min_n (int, default = 5) The minimum number of mice in each group for region comparisons.
#' @param remove (bool, TRUE) Remove any regions in the combined normalized count dataframes that don't have enough n to do a comparison on.
#' These regions are removed across all comparison groups.
#' @param log (bool, TRUE) Save the regions that don't have enough n into a '.csv' file in the output folder.
#' @return e experiment object
#' @export
#' @examples e <- enough_mice_per_group(e, by = c("group", "sex"), min_n = 4, remove = TRUE, log = TRUE)
#'
enough_mice_per_group <- function(e, by = c("group", "sex"), min_n = 5, remove = TRUE, log = TRUE){

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
    n_comb <- get_common_regions(mouse_count, by)
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

  combs <- mouse_count %>% dplyr::ungroup() %>% dplyr::select(dplyr::all_of(by)) %>% dplyr::distinct()
  n_combs <- dim(combs)[1]

  common_regions <-  mouse_count %>% dplyr::group_by(across(all_of(c("acronym")))) %>%
    dplyr::count() %>% dplyr::rename(group_n = n)

  common_regions <- common_regions[common_regions$group_n == n_combs, "acronym"]
  return(common_regions)
}






















