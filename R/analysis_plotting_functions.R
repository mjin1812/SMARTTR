#' @importFrom ggplot2 ggplot aes theme element_text unit geom_tile geom_text scale_fill_gradient2 labs ggsave
NULL

##_____________________ Analysis functions ___________________________


#Calculate percentage colabeled over cfos and eyfp. Specify which ROIS you want


#' Get the percentage of colabelled cells over cfos or eyfp channels
#' @description This analysis will only include common regions that are included in both the colabelled and
#' cfos or eyfp channels. The colabelled percentage of individual animals will be extracted. These can be
#'
#' @param e
#' @param channel (str, default = "cfos") The channel used as denominator in fraction counts.
#' @param save_table (bool, default = TRUE) Whether to save the output table as a csv in the experiment object output folder.
#' @param roi (str, default = NULL) Whether to generate colabelled percentages for only specific regions of interest. The default is
#' @param individual (bool, default = TRUE) Whether the data should include individual mouse colabelled. If FALSE the colabel percentages are
#' averaged together to assess across all rois.
#' @return e experiment object with colabel percentage table
#' @export
#'
#' @examples
get_percent_colabel <-function(e, channel = "eyfp", save_table = TRUE, roi = "dDG", individual = TRUE){
  # # correct mismatched attributes typed by users
  # by <- match_m_attr(by)

  e_info <- attr(e, "info")

  # Get common regions that are present across both channels
  common_reg <- e$combined_normalized_counts[["colabel"]]$acronym %>% unique() %>%
    intersect( e$combined_normalized_counts[[channel]]$acronym )


  # Check if user only wants a specific roi in common regions
  if (!is.null(roi)){
    if (!all(roi %in% common_reg)){
      message(paste0("The roi specified is not a common region found in both channels. ",
                     "Common regions include:\n"))
      for (reg in common_reg){
        message(reg)
      }
      stop("Set the roi to one of these acronyms.")
    } else{
      common_reg <- roi
    }
  }

  # names to join by
  by <- names(e$combined_normalized_counts$colabel) %>% setdiff(c("area.mm2", "volume.mm3",
                                                                  "normalized.count.by.area",
                                                                  "normalized.count.by.volume"))

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


  # Get the grouping names
  attr_group <- names(master_counts)[1: (which(names(master_counts) == "acronym")-1)]

  # Sort the dataframe by these grouping names
  master_counts <- master_counts %>% dplyr::group_by(across(all_of(attr_group))) %>%
    dplyr::arrange(.by_group = TRUE)

  ## Group by and summarize
  if(!individual){
    # Get the mean of the individual groups
    attr_group_by <- c(attr_group[2:length(attr_group)], "acronym", "name")
    col_to_summarize <- c("count.colabel", "count.eyfp", "colabel_percentage")
    master_counts <- master_counts %>% dplyr::group_by(across(all_of(attr_group_by))) %>%
      dplyr::summarise(across(all_of(col_to_summarize), mean),
                       colabel_percentage.sd = mean(colabel_percentage),
                       colabel_percentage.sem = mean(colabel_percentage)/sqrt(n()),
                       n = n()) %>%
      dplyr::rename_with(function(x){paste0(x,".mean")}, all_of(col_to_summarize))
    indiv_or_avg <- "average"
  } else {
    indiv_or_avg <- "individual"
  }

  if (save_table){
      out_path <- file.path(e_info$output_path, paste0("colabel_percentage_", channel, "_", indiv_or_avg, "_.csv"))
      write.csv(master_counts, file = out_path)
      message(paste0("Saved colabel count percentages at location: ", out_path))
  }

  e$colabel_percent[[channel]][[indiv_or_avg]] <- master_counts
  return(e)
}


#' get_correlations
#'
#' @param e experiment object
#' @param by The variable names (columns) to filter the dataframe selection.
#' @param values The value of the variables to filter the combined normalized counts df by.
#' @param channels The channels to correlate.
#' @param p_adjust_method (bool or str, default = "BH") Benjamini-Hochberg method is recommended.
#'  Apply the named method to control for the inflated false discovery rate or FWER. Set to FALSE
#'  to keep "raw" p values.
#' @param alpha The alpha level for p adjustment.
#'
#' @return e experiment object. The experiment object now has a named `correlation_list` object stored in it.
#' The name of the correlation object is the concatenation of the variable values separated by a "_".
#' This name allows for unambiguous identification of different group correlations analysis in the future.
#' @export
#'
#' @examples
get_correlations <- function(e, by = c("sex", "group"), values = c("female", "AD"),
                             channels = c("cfos", "eyfp", "colabel"),  p_adjust_method = "BH", alpha = 0.05){
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
    df_corr <- df_channel %>% dplyr::select(all_of(c(common.regions.ordered))) %>%
      as.matrix() %>% Hmisc::rcorr()

    rows <- rownames(df_corr$P)
    cols <- colnames(df_corr$P)

    # Remove (set to NA) the comparisons that are duplicates

    suppressWarnings(
    for (r in 1:length(rows)){
      na_col <- which(is.na(df_corr$P[r,]))
      df_corr$P[r, na_col:length(cols)] <- NA
    })

    # adjust the p-value for false discovery rate or FWER
    if (!isFALSE(p_adjust_method)){
      # Calculate without removing NAs
      df_corr$P <- df_corr$P %>% p.adjust(method = p_adjust_method) %>%
        matrix(nrow = length(rows), ncol= length(cols), dim = list(rows, cols))
      df_corr$sig <- df_corr$P <= alpha
    }

    corr_list[[channel]] <- df_corr
  }

  # Values title
  values_title <- paste0(values, collapse = "_")
  e$correlation_list[[values_title]] <- structure(corr_list,
                                     class = "correlation_list",
                                     group_by = by,
                                     values = values)
  return(e)
}



#' correlation_diff_permutation
#'
#' Note that these correlation lists must have the same number of channels to compare
#'
#' @param e
#' @param correlation_list_name_1
#' @param correlation_list_name_2
#' @param n_shuffle
#' @param alpha
#' @param p_adjust_method
#' @param ... additional parameters to [RcppAlgos::permuteGeneral()] aside from n, m, Parallel and repetition
#' @return
#' @export
#'
#' @examples
correlation_diff_permutation <- function(e, correlation_list_name_1 = "female_AD",
                                         correlation_list_name_2 = "male_AD",
                                         n_shuffle = 1000,
                                         alpha = 0.05,
                                         p_adjust_method = "BH",
                                         ...){


  # Return the correlations list data showing the grouping and the values
  attr_group_1 <- attributes(e$correlation_list[[correlation_list_name_1]])
  attr_group_2 <- attributes(e$correlation_list[[correlation_list_name_2]])

  # Get overlapping regions between the two correlational datasets
  channels <- attr_group_1$names
  for (channel in channels) {

    # Obtain combined cell count table for this channel
    df_channel <- e$combined_normalized_counts[[channel]]

    # Get cell count table for the two groups
    df_channel_group_1 <- filter_df_by_char_params(df_channel, attr_group_1$group_by, attr_group_1$values)
    df_channel_group_2 <- filter_df_by_char_params(df_channel, attr_group_2$group_by, attr_group_2$values)

    # pivot longer to a horizontal df of regions with a column for the mouse_ID and group order
    df_channel_group_1  <- df_channel_group_1  %>%  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
      tidyr::pivot_wider(names_from = acronym, values_from = normalized.count.by.volume) %>%
      dplyr::mutate(corr_group = correlation_list_name_1) %>% dplyr::relocate(corr_group, .before = 2) %>%
      dplyr::select(-all_of(attr_group_1$group_by))

    df_channel_group_2  <- df_channel_group_2  %>%  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
      tidyr::pivot_wider(names_from = acronym, values_from = normalized.count.by.volume) %>%
      dplyr::mutate(corr_group = correlation_list_name_2) %>% dplyr::relocate(corr_group, .before = 2) %>%
      dplyr::select(-all_of(attr_group_2$group_by))

    # Get common regions between each group
    common_regions_btwn_groups <- intersect(names(df_channel_group_1), names(df_channel_group_2))

    # Sort names into anatomical order
    common_regions_btwn_groups <- sort_anatomical_order(common_regions_btwn_groups)

    # Select the common regions in anatomical order across the two group dataframes
    df_channel_group_1 <- df_channel_group_1 %>% dplyr::select(mouse_ID, corr_group, all_of(common_regions_btwn_groups))
    df_channel_group_2 <- df_channel_group_2 %>% dplyr::select(mouse_ID, corr_group, all_of(common_regions_btwn_groups))

    # Bind group dfs together
    df_channel_groups <- dplyr::bind_rows(df_channel_group_1, df_channel_group_2)

    # Get group region pairwise correlational differences
    group_1_corr <- df_channel_group_1 %>% dplyr::select(-c(mouse_ID:corr_group)) %>%
      as.matrix() %>% Hmisc::rcorr()
    group_2_corr <- df_channel_group_2 %>% dplyr::select(-c(mouse_ID:corr_group)) %>%
      as.matrix() %>% Hmisc::rcorr()
    test_statistic <- group_1_corr$r - group_2_corr$r

    # # Get an array of distribution of correlation differences
    test_statistic_distributions <- permute_corr_diff_distrib(df_channel_groups,
                                                         correlation_list_name_1 = correlation_list_name_1,
                                                         correlation_list_name_2 = correlation_list_name_2,
                                                         n_shuffle = n_shuffle, ...)
    # For each pairwise distribution, sort the values
    test_statistic_distributions <- apply(test_statistic_distributions, 1:2, sort)

    # calculate the p-value of the permutation
    p_matrix <- matrix(nrow = length(common_regions_btwn_groups),
                       ncol = length(common_regions_btwn_groups),
                       dimnames = dimnames(test_statistic))

    l_reg <- length(common_regions_btwn_groups)
    for (i in 1:l_reg){
      for (j in 1:l_reg){
        null_distrib <- test_statistic_distributions[i,j] %>% unlist()
        p_matrix[i,j] <-  (sum(abs(null_distrib) >= abs(test_statistic[i,j])) + 1) / (n_shuffle + 1)
        if(j>=i){
          p_matrix[i,j:l_reg] <- NA
        }
      }
    }

    # adjust the p-value for false discovery rate or FWER
    if (!isFALSE(p_adjust_method)){
      # Calculate without removing NAs
      p_matrix <- p_matrix %>% p.adjust(method = p_adjust_method) %>%
        matrix(nrow = l_reg, ncol= l_reg, dimnames = dimnames(test_statistic))
    }

    # Store the results in the experiment object
    comparison <- paste(correlation_list_name_1,"vs",correlation_list_name_2, sep = "_")
    e$permutation_p_matrix[[comparison]][[channel]] <- list(p_val = p_matrix, sig = p_matrix < alpha)
  }
  return(e)
}










# df_corr %>% purrr::map(tibble::as_tibble) %>%
#   purrr::map(dplyr::mutate(row = names(.), .before = TRUE))
#
#
# # df_corr$r %>% dplyr::mutate(row = names(.), .before)
# #
# %>% purrr::map(dplyr::mutate, row = names(.data))


#
#     corrs <- df_channel %>% dplyr::select(all_of(c(common.regions.ordered))) %>%
#       corrr::correlate(method = method, ...)
#     # %>% dplyr::arrange(-dplyr::row_number())

#
# #Calculate while removing NAs
#
# indices <- which(!is.na(df_corr$P))
# p_adjust <- df_corr$P[indices] %>% p.adjust(method = "BH")
#
# df_corr_p_adjust <-  rep(NA, length = length(rows)^2)
# df_corr_p_adjust[indices] <- p_adjust
# df_corr_p_adjust <- df_corr_p_adjust %>%
#   matrix(nrow = length(rows), ncol= length(cols), dim = list(rows, cols))

























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
# plot_colabel_percent <- function(e, by = c("group", "sex"), channel = " cfos", roi = NULL,
#                                  colors = c("red", "blue"),
#                                  pattern = "stripes"){
#
#
#
#   test <- e$colabel_percent[[channel]] %>% dplyr::group_by(across(all_of(c(by, "acronym")))) %>%
#     dplyr::summarise(colabel_percent = mean(colabel_percentage),
#                      SD = sd(colabel_percentage),
#                      SE = SD/sqrt(length(colabel_percentage)))
#
#
#
#   # Check if user only wants a specific roi in common regions
#     if (!is.null(roi)){
#       if (!all(roi %in% unique(e$colabel_percent[[channel]]$acronym))){
#         message(paste0("The roi specified is not a common region found in both channels. ", "Common regions include:\n"))
#         for (reg in unique(e$colabel_percent[[channel]]$acronym)){
#           message(reg)
#         }
#         stop("Set the roi to one of these acronyms.")
#       } else{
#         test <- test %>% dplyr::filter(acronym %in% roi)
#
#       }
#     }
#
#   ggplot(test, mapping = aes(x=acronym, y=colabel_percent, fill=group)) +
#     geom_bar(stat="identity", position = "dodge")
# }

  ## ggplot barplot

#
#   plt <- ggplot(e$colabel_percent[[channel]]) +
#     ggpattern::geom_col_pattern(by, colabel)
#
#   ggplot(test, ,mapping = aes(x = group, y = colabel_percentage, fill = sex)) +
#     geom_bar(stat="identity", position = "dodge")
#
#   ggplot(data = experiment, mapping = aes(x=date, y=car_count, fill=site)) +
#     geom_bar(stat="identity", position = "dodge")
#


#
# df <- data.frame(level = c("a", "b", "c", 'd'), outcome = c(2.3, 1.9, 3.2, 1))
#
# ggplot(df) +
#   geom_col_pattern(
#     aes(level, outcome, pattern_fill = level),
#     pattern = 'stripe',
#     fill    = 'white',
#     colour  = 'black'
#   ) +
#   theme_bw(18) +
#   theme(legend.position = 'none') +
#   labs(
#     title    = "ggpattern::geom_pattern_col()",
#     subtitle = "pattern = 'stripe'"
#   ) +
#   coord_fixed(ratio = 1/2)
# Generate correlation heatmaps


#' Plot correlation heatmaps
#'
#' @param e experiment object that contains a named correlation_list object generated by [SMARTR::get_correlations()]
#' @param correlation_list_name The name of the correlation object generated by [SMARTR::get_correlations()]
#' @param colors Hexadecimal code for the colors corresponding to the channels attribute of the correlation_list.
#' @param title Title of the plot. If NULL, the grouping variable names will be taken from the correlation_list object annd used as the title.
#' @param height Height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save the correlation heatmap plots into the figures subdirectory of the
#'  the experiment object output folder.

#'
#' @return e experiment object
#' @export
#'
#' @examples

plot_correlation_heatmaps <- function(e, correlation_list_name = "female_AD", colors = c("#be0000", "#00782e", "#0d7983"), save_plot = TRUE,
                                      title = NULL, height = 10, width = 10, image_ext = ".png"){


  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # Get the attributes for plotting
  if (is.null(e$correlation_list[[correlation_list_name]])){
    stop(paste0("Your experiment object doesn't contain the specified correlation_list object to plot! ",
                "\nRun the function, get_correlations() first."))
  }

  cl_attr <- attributes(e$correlation_list[[correlation_list_name]])
  channels <- cl_attr$names
  names(colors) <- channels

  if (is.null(title)){
    cl_attr$value[1] <- stringr::str_to_title(cl_attr$value[1])
    title <- paste(cl_attr$value, collapse = " ")
  }

  # Create plotting theme for the heatmap
  theme.hm <- ggplot2::theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 8),
                    axis.text.y = element_text(vjust = 0.5, size = 8),
                    plot.title = element_text(hjust = 0.5, size = 36),
                    axis.title = element_text(size = 22),
                    legend.text = element_text(size = 22),
                    legend.key.height = unit(100, "points"),
                    legend.title = element_text(size = 22))

  for (channel in channels){
   # Turn into tibble
   corr_df <-  e$correlation_list[[correlation_list_name]][[channel]] %>% purrr::map(tibble::as_tibble)
   val_names <- names(corr_df)

   for (k in 1:length(corr_df)){
     corr_df[[k]] <- tibble::add_column(corr_df[[k]], row_acronym = names(corr_df[[k]]), .before = TRUE ) %>%
       tidyr::pivot_longer(!row_acronym, names_to = "col_acronym", values_to = val_names[k])
     corr_df[[k]]$row_acronym <- factor(corr_df[[k]]$row_acronym, levels = unique(corr_df[[k]]$row_acronym)) # to keep anatomical level order
     corr_df[[k]]$col_acronym <- factor(corr_df[[k]]$col_acronym, levels = unique(corr_df[[k]]$col_acronym)) # to keep the anatomical level order
   }

   # Combined the correlation plots into one dataframe and add a column if there is a significant comparison
   df <- corr_df$r %>% dplyr::left_join(corr_df$n, by = c("row_acronym", "col_acronym")) %>%
     dplyr::left_join(corr_df$P, by = c("row_acronym", "col_acronym")) %>%
     dplyr::left_join(corr_df$sig, by = c("row_acronym", "col_acronym")) %>%
     dplyr::mutate(sig_text = dplyr::if_else(sig == TRUE, "*", ""))


  # Generate a correlation heatmap in anatomical order
  p <-  ggplot(df, aes(row_acronym, col_acronym, fill = r)) +
        geom_tile() +
        geom_text(aes(label = sig_text), size=8, color = "yellow") +
        scale_fill_gradient2(low = "#4f4f4f",mid = "#ffffff", high = colors[[channel]],
                         aesthetics = c("color","fill"), na.value = "grey50")+
        labs(title = title, x = "Brain Region", y = "Brain Region") +
    theme.hm

  if(save_plot){
    # Plot the heatmap
    quartz(width = width, height = height)
    print(p)

    # Create figure directory if it doesn't already exists
    output_dir <-  file.path(attr(e, "info")$output_path, "figures")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    image_file <- file.path(output_dir, paste0("heatmap_", str_replace(title, " ", "_"), "_", channel, image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
    dev.off()
  }

  e$correlation_heatmaps[[correlation_list_name]][[channel]] <- p
  }

  return(e)
}




#' Create a Volcano plot.
#' NOTE THAT THE INPUT FOR THIS FUNCTION NEEDS TO BE REWORKED TO BE BE A DF OF STATISTICALLY DIFFERENT
#' REGION CORRELATIONS BETWEEN TWO COMPARISON GROUPS.
#'
#' @param e
#' @param colors Hexadecimal code for the color of the plot..
#' @param title Title of the plot.
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save the correlation heatmap plots into the figures subdirectory of the
#' @return e experiment object
#' @export
#' @examples
volcano_plot <- function(e, color = c("#be0000"), save_plot = TRUE,
                         title = NULL, height = 10, width = 10, image_ext = ".png"){
}









#________________ Internal Analysis functions _______________


filter_df_by_char_params <- function(df, by, values){
  for (k in 1:length(by)){
    # Convert the variable name into a symbol
    var <- rlang::sym(by[k])
    df <- df %>% dplyr::filter(!!var == values[k])
  }
  return(df)
}


# Sort the dataframes columns to be in anatomical order

sort_anatomical_order <- function(common_regions){

  anatomical.order <- c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")
  common_regions <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure) %>%
    purrr::map(intersect,y=common_regions) %>% unlist()


}





#' Generate array of null distribution of region pairwise correlation differences.
#' @param df
#' @param correlation_list_name_1
#' @param correlation_list_name_2
#' @param n_shuffle
#' @param ...
#'
#' @return
#'
#' @examples
permute_corr_diff_distrib <- function(df, correlation_list_name_1, correlation_list_name_2, n_shuffle = n_shuffle, ...){


  # Create a 3D matrix to hold the correlation distributions
  n_reg <- length(names(df)) - 2
  region_names <- names(df)[3:length(names(df))]
  corr_diff_matrix <- array(dim= c(n_reg, n_reg, n_shuffle))
  dimnames(corr_diff_matrix) <- list(region_names, region_names, 1:n_shuffle)

  # Get permutation sampling combinations
  n_mice <- length(df$mouse_ID)
  possible_perm_combs <-  RcppAlgos::permuteGeneral(n_mice, m = n_mice, Parallel = TRUE, repetition = FALSE, ...)
  sampled_perm_combs <- possible_perm_combs[sample(factorial(n_mice), size = n_shuffle, replace = TRUE),]


  for (n in 1:n_shuffle){
    # reorder the group labels based on the permutation analysis
    comb <- sampled_perm_combs[n, ]
    df$mouse_ID <- df$mouse_ID[comb]
    df$corr_group <- df$corr_group[comb]

    # create matrices as input for rcorr
    matrix_list <-  df %>% dplyr::select(-c(mouse_ID)) %>% dplyr::group_by(corr_group) %>%
      dplyr::group_map(as.matrix, .keep = TRUE)
    element_1_name <- matrix_list[[1]][,"corr_group"] %>% unique()

    # Reorder the matrix in case element 1 doesn't correspond to correlation_list_name_1
    if (element_1_name != correlation_list_name_1){
      matrix_list <- list(matrix_list[[2]], matrix_list[[1]])
    }
    names(matrix_list) <- c(correlation_list_name_1, correlation_list_name_2)

    # calculate R coefficients for each region
    correlations_list <- vector(mode = "list", length = 2)
    names(correlations_list) <- c(correlation_list_name_1, correlation_list_name_2)

    correlations_list[[correlation_list_name_1]] <- matrix_list[[correlation_list_name_1]][,-1] %>% Hmisc::rcorr()
    correlations_list[[correlation_list_name_2]] <- matrix_list[[correlation_list_name_2]][,-1] %>% Hmisc::rcorr()

    # subtract R coefficient differences
    corr_diff_matrix[,,n] <- correlations_list[[correlation_list_name_1]]$r - correlations_list[[correlation_list_name_2]]$r
  }
  return(corr_diff_matrix)

}


  #
  # cross_region_fractions <- cross_region_fractions %>%
  #   mutate(text = signif(total_TP_fraction, digits = 3))x
  #
  # p <- ggplot(cross_region_fractions, aes(train_region,
  #                                         test_region,
  #                                         fill = total_TP_fraction,
  #                                         text = text)) +
  #   geom_tile() +
  #   geom_text(aes(label = text)) +
  #   scale_fill_gradient(low = "white", high = "#136aec") +
  #   labs(x = "Train region", y = "Test region", fill = "TP\n fraction\n retained")


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












