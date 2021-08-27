
##_____________________ Analysis functions ___________________________


#Calculate % colabeled over cfos and eyfp. Specify which ROIS you want


#' Get the percentage of colabelled cells over cfos or eyfp channels
#' @description This analysis will only include common regions that are included in both the colabelled and
#' cfos or eyfp channels. The colabelled percentage of individual animals will be extracted. These can be
#'
#' @param e
#' @param channel (str, default = "cfos") The channel used as denominator in percentage counts.
#' @param save_table (bool, default = TRUE) Whether to save the output table as a csv in the experiment object output folder.
#' @param roi (str, default = NULL) Whether to generate colabelled percentages for only specific regions of interest. The default is
#' to assess across all rois.
#' @return e experiment object with colabel percent table
#' @export
#'
#' @examples
get_percent_colabel <-function(e, channel = "cfos", save_table = TRUE, roi = NULL){
  # # correct mismatched attributes typed by users
  # by <- match_m_attr(by)

  e_info <- attr(e, "info")

  # Get common regions that are present across both channels
  common_reg <- e$combined_normalized_counts[["colabel"]]$acronym %>% unique() %>%
    intersect( e$combined_normalized_counts[[channel]]$acronym )


  # Check if user only wants a specific roi in common regions
  if (!is.null(roi)){
    if (!all(roi %in% common_reg)){
      message(paste0("The roi specified is not a common region found in both channels. ", "Common regions include:\n"))
      for (reg in common_reg){
        message(reg)
      }
      stop("Set the roi to one of these acronyms.")
    } else{
      common_reg <- roi
    }
  }

  # names to join by
  by <- names(colabel_group_counts) %>% setdiff(c("area.mm2", "volume.mm3", "normalized.count.by.area", "normalized.count.by.volume"))

  # Get the colabel counts
  colabel_counts <-  e$combined_normalized_counts$colabel %>% dplyr::filter(acronym %in% common_reg) %>%
    dplyr::select(tidyselect::all_of(by))

  # Get the channel counts
  channel_counts <-  e$combined_normalized_counts[[channel]] %>% dplyr::filter(acronym %in% common_reg) %>%
    dplyr::select(tidyselect::all_of(by))


  # Create master counts table of percentage colabels
  # Inner join deletes any mice with na values for counts

  var <- sym(paste0("count.", channel)) # symbol for mutate function

  master_counts <- colabel_counts %>%
    dplyr::inner_join(channel_counts, suffix = c(".colabel", paste0(".", channel)),
                      by = setdiff(by, "count")) %>%
    dplyr::mutate(colabel_percentage = count.colabel/{{ var }} * 100)

  if (save_table){
      out_path <- file.path(e_info$output_path, paste0("colabel_percentage_", channel,".csv"))
      write.csv(master_counts, file = out_path)
      message(paste0("Saved colabel count percentages at location: ", out_path))
  }

  e$colabel_percent[[channel]] <- master_counts
  return(e)
}





#' get_correlations
#'
#' @param e experiment object
#' @param by The variable names (columns) to filter the dataframe selection.
#' @param values The value of the variables to filter the combined normalized counts df by.
#' @param channels The channels to correlate.
#' @param method Correlational method to pass to [corrr::correlate()] function.
#' @param ,,, additional parameters to pass to the [corrr:correlate()] function aside from method, x, and y.
#'
#' @return
#' @export
#'
#' @examples
get_correlations <- function(e, by = c("sex", "group"), values = c("female", "AD"),
                             channels = c("cfos", "eyfp", "colabel"),  method = "pearson", ...){


  corr_list <- list()
  names(corr_list) <- names(channels)

  for (channel in channels){

    df_channel <-  e$combined_normalized_counts[[channel]]

    for (k in 1:length(by)){
      # Convert the variable name into a symbol
      var <- rlang::sym(by[k])
      df_channel <- df_channel %>% dplyr::filter(!!var == values[k])
    }

    df_channel <- df_channel %>%  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
      tidyr::pivot_wider(names_from = acronym, values_from = normalized.count.by.volume)

    # Rearrange the correlations to be in "anatomical order"
    anatomical.order <- c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")
    common.regions <-  df_channel %>% dplyr::select(-all_of(c('mouse_ID', by))) %>% names()
    common.regions.ordered <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure) %>%
      purrr::map(intersect,y=common.regions) %>% unlist()


    # Select the order of the columns, perform the correlations, ignore the mouse_ID and group columns
    corrs <- df_channel %>% dplyr::select(all_of(c(common.regions.ordered))) %>%
      corrr::correlate(method = method, ...)
    # %>% dplyr::arrange(-dplyr::row_number())

    corr_list[[channel]] <- corrs
  }

  e$corr_list <- corr_list

  return(e)
}



##_____________________ Plotting functions ___________________________



#' plot_colabel_percent
#'
#' @param e
#' @param by
#' @param roi
#'
#' @return
#' @export
#'
#' @examples
plot_colabel_percent <- function(e, by = c("group", "sex"), roi = NULL){




}



# Generate correlational heatmaps



#' Plot correlation heatmaps
#'
#' @param e
#' @param by
#' @param channel
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
plot_heatmaps <- function(e, by = c("group", "sex"), channel = c("cfos", "eyfp", "colabel"),
                            colors = c("red", "green", "blue")){



  #
  #
  # cross_region_fractions <- cross_region_fractions %>%
  #   mutate(text = signif(total_TP_fraction, digits = 3))
  #
  # p <- ggplot(cross_region_fractions, aes(train_region,
  #                                         test_region,
  #                                         fill = total_TP_fraction,
  #                                         text = text)) +
  #   geom_tile() +
  #   geom_text(aes(label = text)) +
  #   scale_fill_gradient(low = "white", high = "#136aec") +
  #   labs(x = "Train region", y = "Test region", fill = "TP\n fraction\n retained")





}


# install.packages("Hmisc")
#
# library(palmerpenguins)
# library(dplyr)
# penguins_cor <- palmerpenguins::penguins %>%
#   dplyr::select(bill_length_mm, bill_depth_mm, flipper_length_mm) %>%
#   correlate()
#


#
# # Reformat the data into wide format
#
# pivoted <- e$combined_normalized_counts$cfos %>% dplyr::filter(sex == "female", group == "AD") %>%
#   dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
#   tidyr::pivot_wider(  names_from = acronym, values_from = normalized.count.by.volume)
#
#
# # Perform the correlations, ignore the mouse_ID and group columns
# corrs <- pivoted %>% dplyr::select(-(mouse_ID:group)) %>% corrr::correlate(method = "pearson", use ="pairwise.complete.obs")
#
#
# # plot heatmap using the intrinsic heatmap function
# corrs %>% shave() %>% rplot()
#
# # plot a network plot using
#   network_plot(min_cor = .2)




#
#
#   p <- ggplot(cross_region_fractions, aes(train_region,
#                                           test_region,
#                                           fill = total_FP_fraction,
#                                           text = text)) +
#     geom_tile() +
#     geom_text(aes(label = text)) +
#     scale_fill_gradient(low = "white", high = "#cc3337") +
#     labs(x = "Train region", y = "Test region", fill = "FP\n fraction\n retained")












