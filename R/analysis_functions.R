
##_____________________ Analysis functions ___________________________

#' Get the percentage of colabelled cells over either cfos or eyfp channels.
#' @description This analysis will only include common regions that are included in both the colabelled and
#' cfos or eyfp channels. The colabelled percentage of individual animals will be calculated with the option to export the data.
#' @param e experiment object
#' @param by (str) Attribute names to group by, e.g. by = c("group", "sex"). These will generate analysis subgroups that are
#' averaged together to assess across all rois.
#' @param colabel_channel (str, default = "colabel") The channel used as the numerator in fraction counts. The string 'colabel' in the pipeline
#' refers to colocalized 'eyfp' and 'cfos' channels. For other colocalized channels, import the channel using [SMARTTR::import_segmentation_custom()]
#' or your own customized import channel.
#' @param channel (str, default = "eyfp") The channel used as denominator in fraction counts.
#' @param save_table (bool, default = TRUE) Whether to save the output table as a csv in the experiment object output folder.
#' @param rois (str, default = NULL) Whether to generate colabelled percentages for only specific regions of interest, e.g. rois = c("HY", "DG").
#' Child regions of specified rois will also be searched for.
#' @param individual (bool, default = FALSE) Whether the data should include individual mouse colabelled percentages rather than the average.
#' If FALSE the colabel percentages are averaged across all analysis subgroups determined by the `by` parameter
#' @return e experiment object with colabelled percentage table stored in it.
#' @export
#' @examples
#' \dontrun{
#' e <- get_percent_colabel(e, c("group", "sex"), channel = "eyfp"))
#' }
get_percent_colabel <- function(e, by, colabel_channel = "colabel",
                               channel = "eyfp", save_table = TRUE, rois = NULL, individual = TRUE){
  e_info <- attr(e, "info")

  # Get common regions that are present across both channels
  common_reg <- e$combined_normalized_counts[[colabel_channel]]$acronym %>% unique() %>%
    intersect( e$combined_normalized_counts[[channel]]$acronym )

  if (!is.null(rois)){
    common_reg <- rois_intersect_region_list(common_reg, rois)
  }

  if (is.null(by)){ # Join by all attributes in the dataset
    by <- names(e$combined_normalized_counts$colabel) %>% setdiff(c("area.mm2", "volume.mm3",
                                                                    "normalized.count.by.area",
                                                                    "normalized.count.by.volume"))
  } else { # join by user defined attributes
    by <- c("mouse_ID", by, "acronym", "name", "count")
  }
  colabel_counts <-  e$combined_normalized_counts[[colabel_channel]] %>% dplyr::filter(.data$acronym %in% common_reg) %>%
    dplyr::select(all_of(by))
  channel_counts <-  e$combined_normalized_counts[[channel]] %>% dplyr::filter(.data$acronym %in% common_reg) %>%
    dplyr::select(all_of(by))

  # Create master counts table of percentage colabels
  # Inner join deletes any mice with na values for counts
  var <- rlang::sym(paste0("count.", channel)) # symbol for mutate
  master_counts <- colabel_counts %>%
    dplyr::inner_join(channel_counts, suffix = c(".colabel", paste0(".", channel)),
                      by = setdiff(by, "count")) %>%
    dplyr::mutate(colabel_percentage = .data$count.colabel/{{ var }} * 100)

  attr_group <- names(master_counts)[1: (which(names(master_counts) == "acronym")-1)]

  # Sort the dataframe by these grouping names
  master_counts <- master_counts %>% dplyr::group_by(across(all_of(attr_group))) %>%
    dplyr::arrange(.by_group = TRUE)

  if(!individual){
    # Get the mean of the individual groups
    attr_group_by <- c(attr_group[2:length(attr_group)], "acronym", "name")
    col_to_summarize <- c("count.colabel", paste0("count.", channel), "colabel_percentage")
    master_counts <- master_counts %>%
      dplyr::group_by(across(all_of(attr_group_by))) %>%
      dplyr::summarise(n = dplyr::n(),
                       colabel_percentage.sd = sd(.data$colabel_percentage),
                       colabel_percentage.sem = .data$colabel_percentage.sd/sqrt(.data$n),
                       across(all_of(col_to_summarize), mean)) %>%
      dplyr::rename_with(function(x){paste0(x,".mean")}, all_of(col_to_summarize)) %>%
      dplyr::relocate(all_of(c("colabel_percentage.sd", "colabel_percentage.sem", "n")), .after = last_col())
    indiv_or_avg <- "average"
  } else {
    indiv_or_avg <- "individual"
  }
  if (save_table){
    # Create table directory if it doesn't already exists
    output_dir <-  file.path(attr(e, "info")$output_path, "tables")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
      out_path <- file.path(output_dir, paste0(make.names(colabel_channel),"_percentage_", channel, "_", indiv_or_avg, ".csv"))
      write.csv(master_counts, file = out_path)
      message(paste0("Saved ", colabel_channel, " count percentages at location: ", out_path))
  }
  e$colabel_percent[[channel]][[indiv_or_avg]] <- master_counts
  return(e)
}


#' Get regional cross correlations and their p-values in a correlation list object.
#' @description This analysis will get regional cross correlations based on cell counts normalized by region volume.
#' @param e experiment object
#' @param by (str) Attribute names to group by, e.g. c("sex", "group")
#' @param values (str) The respective values of the attributes entered for the `by` parameter to generate a specific analysis group,
#' e.g.values = c("female", "AD"). Length must be the same as `by`.
#' @param channels (str, channels =  c("cfos", "eyfp", "colabel") The channels to process.
#' @param p_adjust_method (bool or str, default = "none") This parameter is fed into the p.adjust function. Options: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")
#'  Apply the named method to control for the inflated false discovery rate or family wise error rate (FWER). Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param region_order (list, default = NULL)  optional list with the first element named "acronym" supplying a vector as region acronyms and the second element named
#' "order"  supplying an vector of integers determining numerical order, e.g. 1, 1, 2, 2.
#' @param alpha (num, default = 0.05) The alpha level for significance applied AFTER p-adjustment.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param method (str, default = "pearson", options = c("pearson", "spearman")) Specifies the type of correlations to compute. Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements,
#' @param anatomical.order (str, default = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The order of the regions to plot in the heatmaps.
#' using midranks for ties. See also [Hmisc::rcorr()]
#' @return e experiment object. The experiment object now has a named `correlation_list` object stored in it.
#' The name of the correlation object is the concatenation of the variable values separated by a "_".
#' This name allows for unambiguous identification of different analysis subgroups in the future.
#' @export
#' @examples
#' \dontrun{
#' e <- get_correlations(e, by = c("sex", "group"), values = c("female", "AD"),
#' channels = c("cfos", "eyfp", "colabel"), alpha = 0.05)
#' }
#' @seealso [Hmisc::rcorr()]
get_correlations <- function(e, by, values,
                             channels = c("cfos", "eyfp", "colabel"),
                             p_adjust_method = "none",
                             alpha = 0.05,
                             ontology = "allen",
                             method = "pearson",
                             anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                             region_order = NULL){
  corr_list <- list()
  names(corr_list) <- names(channels)

  if (length(by) != length(values)){
    warning("Your values vector seems too long. Please double check!")
    values <- values[1:length(by)]
  }

  for (channel in channels){
    df_channel <-  e$combined_normalized_counts[[channel]] %>%
      filter_df_by_char_params(by, values) %>% dplyr::distinct()

    df_channel <- df_channel %>%  dplyr::select(any_of("mouse_ID"):any_of("acronym"), "normalized.count.by.volume") %>%
      tidyr::pivot_wider(names_from = any_of("acronym"), values_from = any_of("normalized.count.by.volume"), values_fill = NA)

    # Rearrange the correlations to be in "anatomical order"
    common.regions <-  df_channel %>% dplyr::ungroup() %>% dplyr::select(-any_of(c('mouse_ID', SMARTTR::attr2match$mouse))) %>% names()

    if (is.null(region_order)){
      if (tolower(ontology) == "allen"){
        common.regions.ordered <- anatomical.order %>% purrr::map(get.sub.structure) %>%
          purrr::map(intersect, y=common.regions) %>% unlist()
      } else {
        common.regions.ordered <- anatomical.order %>% purrr::map(get.sub.structure.custom, ontology = ontology) %>%
          purrr::map(intersect, y=common.regions) %>% unlist()
      }
      } else if (is.list(region_order)){
      common.regions.ordered <- region_order$acronym[region_order$order]
      common.regions.ordered <- intersect(common.regions.ordered , common.regions)
    } else {
      stop("You did not supply a list object for the variable `region_order`."
      )
    }
    df_channel <- df_channel %>% dplyr::ungroup() %>% dplyr::select(all_of(common.regions.ordered)) %>%
      as.matrix()

    df_corr <-  try_correlate(df_channel, type = method)

    if (!isFALSE(p_adjust_method)){
      lowertri <- df_corr$P %>%  lower.tri(diag = FALSE)
      lower_p <- df_corr$P[lowertri]
      uppertri <- df_corr$P %>%  upper.tri(diag = FALSE)
      upper_p <- df_corr$P[uppertri]

      lower_p <- lower_p %>% stats::p.adjust(method = p_adjust_method)
      df_corr$P[lowertri] <- lower_p
      upper_p <- upper_p %>% stats::p.adjust(method = p_adjust_method)
      df_corr$P[uppertri] <- upper_p
    }
    df_corr$sig <- df_corr$P <= alpha
    diag(df_corr$sig) <- FALSE
    corr_list[[channel]] <- df_corr
  }
  values_title <- paste0(values, collapse = "_")
  e$correlation_list[[values_title]] <- structure(corr_list,
                                     class = "correlation_list",
                                     group_by = by,
                                     values = values)
  return(e)
}


#' Get regional cross correlations and their p-values in a correlation list object.
#' @description This analysis will get regional cross correlations across left and right brain regions based on cell counts normalized by region volume.
#' This is essentially the same functionality as `get_correlations()` except it supports division of regions by hemisphere.
#' @param e experiment object
#' @param by (str) Attribute names to group by, e.g. c("sex", "group")
#' @param values (str) The respective values of the attributes entered for the `by` parameter to generate a specific analysis group,
#' e.g.values = c("female", "AD"). Length must be the same as `by`.
#' @param channels (str, channels =  c("cfos", "eyfp", "colabel") The channels to process.
#' @param p_adjust_method (bool or str, default = "none") This parameter is fed into the p.adjust function. Options: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")
#'  Apply the named method to control for the inflated false discovery rate or family wise error rate (FWER). Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param region_order (list, default = NULL)  optional list with the first element named "acronym" supplying a vector as region acronyms and the second element named
#' "order"  supplying an vector of integers determining numerical order, e.g. 1, 1, 2, 2.
#' @param alpha (num, default = 0.05) The alpha level for significance applied AFTER p-adjustment.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param method (str, default = "pearson", options = c("pearson", "spearman")) Specifies the type of correlations to compute. Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements,
#' @param anatomical.order (str, default = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The order of the regions to plot in the heatmaps.
#' using midranks for ties. See also [Hmisc::rcorr()]
#' @return e experiment object. The experiment object now has a named `correlation_list` object stored in it.
#' The name of the correlation object is the concatenation of the variable values separated by a "_".
#' This name allows for unambiguous identification of different analysis subgroups in the future.
#' @examples
#' \dontrun{
#' e <- get_correlations(e, by = c("sex", "group"), values = c("female", "AD"),
#' channels = c("cfos", "eyfp", "colabel"), alpha = 0.05)
#' }
#' @seealso [Hmisc::rcorr()]
get_correlations_hemi <- function(e, by, values,
                             channels = c("cfos", "eyfp", "colabel"),
                             p_adjust_method = "none",
                             alpha = 0.05,
                             ontology = "allen",
                             method = "pearson",
                             anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                             region_order = NULL){
  corr_list <- list()
  names(corr_list) <- names(channels)

  if (length(by) != length(values)){
    warning("Your values vector seems too long. Please double check!")
    values <- values[1:length(by)]
  }

  for (channel in channels){
    df_channel <-  e$combined_normalized_counts[[channel]] %>%
      filter_df_by_char_params(by, values) %>% dplyr::distinct()

    df_channel <- df_channel %>% mutate(acronym =  ifelse(hemisphere=="left", paste0(acronym, "_L"), paste0(acronym, "_R"))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of("hemisphere"))


    df_channel <- df_channel %>%  dplyr::select(any_of("mouse_ID"):any_of("acronym"), "normalized.count.by.volume") %>%
      tidyr::pivot_wider(names_from = any_of("acronym"), values_from = any_of("normalized.count.by.volume"), values_fill = NA)


    # Rearrange the correlations to be in "anatomical order"
    common.regions.hemi <-  df_channel %>% dplyr::ungroup() %>% dplyr::select(-any_of(c('mouse_ID', SMARTTR::attr2match$mouse))) %>% names()
    common.regions <- common.regions.hemi %>% stringr::str_replace_all("_L", "") %>% stringr::str_replace_all("_R", "") %>% unique()

    if (is.null(region_order)){
      if (tolower(ontology) == "allen"){
        common.regions.ordered <- anatomical.order %>% purrr::map(get.sub.structure) %>%
          purrr::map(intersect, y=common.regions) %>% unlist()
      } else {
        common.regions.ordered <- anatomical.order %>% purrr::map(get.sub.structure.custom, ontology = ontology) %>%
          purrr::map(intersect, y=common.regions) %>% unlist()
      }
    } else if (is.list(region_order)){
      common.regions.ordered <- region_order$acronym[region_order$order]
      common.regions.ordered <- intersect(common.regions.ordered , common.regions)
    } else {
      stop("You did not supply a list object for the variable `region_order`."
      )
    }


    ############### left off here
    common.regions.ordered.hemi <- c(paste(common.regions.ordered, "_L", sep =""), paste(common.regions.ordered, "_R", sep = ""))


    df_channel <- df_channel %>% dplyr::ungroup() %>% dplyr::select(any_of(common.regions.ordered.hemi)) %>%
      as.matrix()

    df_corr <-  try_correlate(df_channel, type = method)

    if (!isFALSE(p_adjust_method)){
      lowertri <- df_corr$P %>%  lower.tri(diag = FALSE)
      lower_p <- df_corr$P[lowertri]
      uppertri <- df_corr$P %>%  upper.tri(diag = FALSE)
      upper_p <- df_corr$P[uppertri]

      lower_p <- lower_p %>% stats::p.adjust(method = p_adjust_method)
      df_corr$P[lowertri] <- lower_p
      upper_p <- upper_p %>% stats::p.adjust(method = p_adjust_method)
      df_corr$P[uppertri] <- upper_p
    }
    df_corr$sig <- df_corr$P <= alpha
    diag(df_corr$sig) <- FALSE
    corr_list[[channel]] <- df_corr
  }
  values_title <- paste0(values, collapse = "_")
  e$correlation_list[[values_title]] <- structure(corr_list,
                                                  class = "correlation_list",
                                                  group_by = by,
                                                  values = values)
  return(e)
}



#' This function performs a permutation analysis to compare the region pairwise correlation coefficients between two different analysis groups.
#'
#' @description The data from two different analysis groups are compared by specifying the
#' `correlation_list_name_1` and `correlation_list_name_2` parameters. Note that both of these analysis groups must have the same
#' number of channels to compare. The functions `get_correlations()` needs to have been run for each of these analysis groups prior to
#' running this function. The test statistics used is the pearson values of those in correlation_list_name_2 subtracted from corresponding Pearson values in correlation_list_name_1.
#'
#' @param e experiment object
#' @param correlation_list_name_1 (str) The name of the correlation list object used as the first group for comparison.
#' @param correlation_list_name_2 (str) The name of the correlation list object used as the second group for comparison.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param n_shuffle (int, default = 1000) The number of permutation shuffles.
#' @param method (str, default = "pearson", options = c("pearson", "spearman")) Specifies the type of correlations to compute.
#' Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements, using midranks for ties. See also [Hmisc::rcorr()]
#' @param seed (int, default = 5) Random seed for future replication.
#' @param p_adjust_method (bool or str, default = "none")
#'  Apply the named method to control for the inflated false discovery rate or FWER. Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param alpha (float, default = 0.05) The alpha cutoff for significance between region pairwise correlation differences
#' @param progressbar (bool, default = TRUE) Display a progress bar for the processing time of the permutation.
#' @param ... additional parameters to pass to `permute_corr_diff_distrib()`
#'
#' @return e experiment object. The experiment object now has a list called `permutation_p_matrix` stored in it. Elements of this `permutation_p_matrix` list are
#' the outputs of different permutation comparison analyses. These elements are named by the groups that were compared.
#' @export
#' @seealso [SMARTTR::get_correlations()]
#' @examples
#' \dontrun{e <- correlation_diff_permutation(sundowning,
#'                                             correlation_list_name_1 = "female_AD",
#'                                             correlation_list_name_2 = "female_control",
#'                                             channels = c("cfos", "eyfp", "colabel"),
#'                                             n_shuffles = 1000,
#'                                             seed = 5,
#'                                             p_adjust_method = "none"
#'                                             alpha = 0.001
#'                                             )
#'}
#'
correlation_diff_permutation <- function(e,
                                         correlation_list_name_1,
                                         correlation_list_name_2,
                                         channels = c("cfos", "eyfp", "colabel"),
                                         n_shuffle = 1000,
                                         method = "pearson",
                                         seed = 5,
                                         p_adjust_method = "none",
                                         alpha = 0.05,
                                         progressbar = TRUE,
                                         ...){

  attr_group_1 <- attributes(e$correlation_list[[correlation_list_name_1]])
  attr_group_2 <- attributes(e$correlation_list[[correlation_list_name_2]])

  p_matrix_list <- vector(mode = "list", length = length(channels))
  names(p_matrix_list) <- channels

  for (channel in channels) {
    df_channel <- e$combined_normalized_counts[[channel]]
    df_channel_group_1 <- filter_df_by_char_params(df_channel, attr_group_1$group_by, attr_group_1$values)
    df_channel_group_2 <- filter_df_by_char_params(df_channel, attr_group_2$group_by, attr_group_2$values)

    df_channel_group_1  <- df_channel_group_1  %>%  dplyr::select(any_of("mouse_ID"):any_of("acronym"), "normalized.count.by.volume") %>%
      tidyr::pivot_wider(names_from = any_of("acronym"), values_from = any_of("normalized.count.by.volume")) %>%
      dplyr::mutate(corr_group = correlation_list_name_1) %>% dplyr::relocate(any_of("corr_group"), .before = 2) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of(SMARTTR::attr2match$mouse))

    df_channel_group_2  <- df_channel_group_2  %>%  dplyr::select(any_of("mouse_ID"):any_of("acronym"), "normalized.count.by.volume") %>%
      tidyr::pivot_wider(names_from = any_of("acronym"), values_from = any_of("normalized.count.by.volume")) %>%
      dplyr::mutate(corr_group = correlation_list_name_2) %>% dplyr::relocate(any_of("corr_group"), .before = 2) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of(SMARTTR::attr2match$mouse))

    common_regions_btwn_groups <- intersect(names(df_channel_group_1), names(df_channel_group_2))
    df_channel_group_1 <- df_channel_group_1 %>% dplyr::select(any_of(c("mouse_ID", "corr_group")), all_of(common_regions_btwn_groups)) %>%
      dplyr::mutate(mouse_ID = as.character(.data$mouse_ID))
    df_channel_group_2 <- df_channel_group_2 %>% dplyr::select(any_of(c("mouse_ID", "corr_group")), all_of(common_regions_btwn_groups)) %>%
      dplyr::mutate(mouse_ID = as.character(.data$mouse_ID))
    df_channel_groups <- dplyr::bind_rows(df_channel_group_1, df_channel_group_2)

    group_1_corr <- df_channel_group_1 %>% dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>% try_correlate(type = method)
    group_2_corr <- df_channel_group_2 %>%  dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>% try_correlate(type = method)
    test_statistic <- group_2_corr$r - group_1_corr$r

    suppressWarnings(
    test_statistic_distributions <- permute_corr_diff_distrib(df_channel_groups,
                                                         correlation_list_name_1 = correlation_list_name_1,
                                                         correlation_list_name_2 = correlation_list_name_2,
                                                         n_shuffle = n_shuffle,
                                                         seed = seed, method = method,  progressbar = progressbar, ...)
    )

    # message("dim before sort", dim(test_statistic_distributions))
    # For each pairwise distribution, sort the values
    # test_statistic_distributions <- apply(test_statistic_distributions, 1:2, sort)
    # message("dim after sort", dim(test_statistic_distributions))
    # do.call(sort, 1:2)
    # test_statistic_distributions <- test_statistic_distributions %>% aperm(c(1, 2, 3))
    # aperm(c(3, 2, 1))

    p_matrix <- matrix(nrow = dim(test_statistic)[1],
                       ncol = dim(test_statistic)[1],
                       dimnames = dimnames(test_statistic))
    l_reg <- dimnames(test_statistic)[1] %>% unlist()
    for (i in 1:length(l_reg)){
      for (j in 1:length(l_reg)){
        null_distrib <- test_statistic_distributions[l_reg[i],l_reg[j],] %>% unlist()
        p_matrix[l_reg[i],l_reg[j]] <-  (sum(abs(null_distrib) >= abs(test_statistic[l_reg[i],l_reg[j]])) + 1) / (n_shuffle + 1)

        if(j>=i){
          p_matrix[i,j:length(l_reg)] <- NA
        }
      }
    }

    if (!isFALSE(p_adjust_method)){
      p_matrix <- p_matrix %>% stats::p.adjust(method = p_adjust_method) %>%
        matrix(nrow = length(l_reg), ncol= length(l_reg), dimnames = dimnames(test_statistic))
    }
    p_matrix_list[[channel]] <- list(p_val = p_matrix,
                                     sig = p_matrix < alpha,
                                     test_statistic = test_statistic,
                                     group_1_pearson = group_1_corr$r,
                                     group_2_pearson = group_2_corr$r,
                                     permutation_null_distribution = test_statistic_distributions,
                                     alpha = alpha)
  }
  comparison <- paste(correlation_list_name_1,"vs",correlation_list_name_2, sep = "_")
  if(is.null(e$permutation_p_matrix)){
    e$permutation_p_matrix <- list()
  }
  e$permutation_p_matrix[[comparison]] <- p_matrix_list
  return(e)
}


#' This function performs a permutation analysis to compare the region pairwise correlation coefficients between two different analysis groups.
#'
#' @description The data from two different analysis groups are compared by specifying the
#' `correlation_list_name_1` and `correlation_list_name_2` parameters. Note that both of these analysis groups must have the same
#' number of channels to compare. The functions `get_correlations()` needs to have been run for each of these analysis groups prior to
#' running this function. The test statistics used is the pearson values of those in correlation_list_name_2 subtracted from corresponding Pearson values in correlation_list_name_1.
#'
#' @param e experiment object
#' @param correlation_list_name_1 (str) The name of the correlation list object used as the first group for comparison.
#' @param correlation_list_name_2 (str) The name of the correlation list object used as the second group for comparison.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param n_shuffle (int, default = 1000) The number of permutation shuffles.
#' @param method (str, default = "pearson", options = c("pearson", "spearman")) Specifies the type of correlations to compute.
#' Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements, using midranks for ties. See also [Hmisc::rcorr()]
#' @param seed (int, default = 5) Random seed for future replication.
#' @param p_adjust_method (bool or str, default = "none")
#'  Apply the named method to control for the inflated false discovery rate or FWER. Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param alpha (float, default = 0.05) The alpha cutoff for significance between region pairwise correlation differences
#' @param progressbar (bool, default = TRUE) Display a progress bar for the processing time of the permutation.
#' @param ... additional parameters to pass to `permute_corr_diff_distrib()`
#'
#' @return e experiment object. The experiment object now has a list called `permutation_p_matrix` stored in it. Elements of this `permutation_p_matrix` list are
#' the outputs of different permutation comparison analyses. These elements are named by the groups that were compared.
#' @export
#' @seealso [SMARTTR::get_correlations()]
#' @examples
#' \dontrun{e <- correlation_diff_permutation(sundowning,
#'                                             correlation_list_name_1 = "female_AD",
#'                                             correlation_list_name_2 = "female_control",
#'                                             channels = c("cfos", "eyfp", "colabel"),
#'                                             n_shuffles = 1000,
#'                                             seed = 5,
#'                                             p_adjust_method = "none"
#'                                             alpha = 0.001
#'                                             )
#'}
#'
correlation_diff_permutation_hemi <- function(e,
                                         correlation_list_name_1,
                                         correlation_list_name_2,
                                         channels = c("cfos", "eyfp", "colabel"),
                                         n_shuffle = 1000,
                                         method = "pearson",
                                         seed = 5,
                                         p_adjust_method = "none",
                                         alpha = 0.05,
                                         progressbar = TRUE,
                                         ...){

  attr_group_1 <- attributes(e$correlation_list[[correlation_list_name_1]])
  attr_group_2 <- attributes(e$correlation_list[[correlation_list_name_2]])

  p_matrix_list <- vector(mode = "list", length = length(channels))
  names(p_matrix_list) <- channels

  for (channel in channels) {
    df_channel <- e$combined_normalized_counts[[channel]]

    df_channel <- df_channel %>% mutate(acronym =  ifelse(hemisphere=="left", paste0(acronym, "_L"), paste0(acronym, "_R"))) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of("hemisphere"))

    df_channel_group_1 <- filter_df_by_char_params(df_channel, attr_group_1$group_by, attr_group_1$values)
    df_channel_group_2 <- filter_df_by_char_params(df_channel, attr_group_2$group_by, attr_group_2$values)

    df_channel_group_1  <- df_channel_group_1  %>%  dplyr::select(any_of("mouse_ID"):any_of("acronym"), "normalized.count.by.volume") %>%
      tidyr::pivot_wider(names_from = any_of("acronym"), values_from = any_of("normalized.count.by.volume")) %>%
      dplyr::mutate(corr_group = correlation_list_name_1) %>% dplyr::relocate(any_of("corr_group"), .before = 2) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of(SMARTTR::attr2match$mouse))

    df_channel_group_2  <- df_channel_group_2  %>%  dplyr::select(any_of("mouse_ID"):any_of("acronym"), "normalized.count.by.volume") %>%
      tidyr::pivot_wider(names_from = any_of("acronym"), values_from = any_of("normalized.count.by.volume")) %>%
      dplyr::mutate(corr_group = correlation_list_name_2) %>% dplyr::relocate(any_of("corr_group"), .before = 2) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of(SMARTTR::attr2match$mouse))

    common_regions_btwn_groups <- intersect(names(df_channel_group_1), names(df_channel_group_2))
    df_channel_group_1 <- df_channel_group_1 %>% dplyr::select(any_of(c("mouse_ID", "corr_group")), all_of(common_regions_btwn_groups)) %>%
      dplyr::mutate(mouse_ID = as.character(.data$mouse_ID))
    df_channel_group_2 <- df_channel_group_2 %>% dplyr::select(any_of(c("mouse_ID", "corr_group")), all_of(common_regions_btwn_groups)) %>%
      dplyr::mutate(mouse_ID = as.character(.data$mouse_ID))
    df_channel_groups <- dplyr::bind_rows(df_channel_group_1, df_channel_group_2)

    group_1_corr <- df_channel_group_1 %>% dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>% try_correlate(type = method)
    group_2_corr <- df_channel_group_2 %>%  dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>% try_correlate(type = method)
    test_statistic <- group_2_corr$r - group_1_corr$r

    suppressWarnings(
    test_statistic_distributions <- permute_corr_diff_distrib(df_channel_groups,
                                                         correlation_list_name_1 = correlation_list_name_1,
                                                         correlation_list_name_2 = correlation_list_name_2,
                                                         n_shuffle = n_shuffle,
                                                         seed = seed, method = method,  progressbar = progressbar, ...)
    )

    # message("dim before sort", dim(test_statistic_distributions))
    # For each pairwise distribution, sort the values
    # test_statistic_distributions <- apply(test_statistic_distributions, 1:2, sort)
    # message("dim after sort", dim(test_statistic_distributions))
    # do.call(sort, 1:2)
    # test_statistic_distributions <- test_statistic_distributions %>% aperm(c(1, 2, 3))
    # aperm(c(3, 2, 1))

    p_matrix <- matrix(nrow = dim(test_statistic)[1],
                       ncol = dim(test_statistic)[1],
                       dimnames = dimnames(test_statistic))
    l_reg <- dimnames(test_statistic)[1] %>% unlist()
    for (i in 1:length(l_reg)){
      for (j in 1:length(l_reg)){
        null_distrib <- test_statistic_distributions[l_reg[i],l_reg[j],] %>% unlist()
        p_matrix[l_reg[i],l_reg[j]] <-  (sum(abs(null_distrib) >= abs(test_statistic[l_reg[i],l_reg[j]])) + 1) / (n_shuffle + 1)

        if(j>=i){
          p_matrix[i,j:length(l_reg)] <- NA
        }
      }
    }

    if (!isFALSE(p_adjust_method)){
      p_matrix <- p_matrix %>% stats::p.adjust(method = p_adjust_method) %>%
        matrix(nrow = length(l_reg), ncol= length(l_reg), dimnames = dimnames(test_statistic))
    }
    p_matrix_list[[channel]] <- list(p_val = p_matrix,
                                     sig = p_matrix < alpha,
                                     test_statistic = test_statistic,
                                     group_1_pearson = group_1_corr$r,
                                     group_2_pearson = group_2_corr$r,
                                     permutation_null_distribution = test_statistic_distributions,
                                     alpha = alpha)
  }
  comparison <- paste(correlation_list_name_1,"vs",correlation_list_name_2, sep = "_")
  if(is.null(e$permutation_p_matrix)){
    e$permutation_p_matrix <- list()
  }
  e$permutation_p_matrix[[comparison]] <- p_matrix_list
  return(e)
}



#' Export the permutation results as a csv file. This automatically saves into the tables folder.
#'
#' @param e experiment object
#' @param permutation_groups (str vec, default = "all") Default is to export all analyses run. Supply names of the specific analyses found under `e$permutation_p_matrix %>% names()` if you want to export specific groups.
#' @param channels (str, channels =  c("cfos", "eyfp", "colabel") The channels to process.
#' @param ontology (str, default = "allen") Set to "unified" for Kim Lab unified ontology
#' @param filter_significant (bool, default = TRUE) If FALSE, this keeps all comparisons. Otherwise exports only the most different and significant permutations.
#' @export
#' @return NULL
#' @examples
#' \dontrun{
#' e <- export_permutation_results(e, permutation_groups = "all", filter_significant = TRUE)
#' }
#'
export_permutation_results <- function(e,
                                       permutation_groups = "all",
                                       channels = c("cfos"),
                                       ontology = "allen",
                                       filter_significant =  TRUE){
  if (permutation_groups == "all") {
    permutation_groups <- e$permutation_p_matrix %>% names()
  } else{
    # Check that permutation_groups is in the names vector for e$permutation_p_matrix
    if (!all(permutation_groups %in% names(e$permutation_p_matrix))){
      stop("The permutation groups you have supplied are not all in the experiment object. Please fix the names supplied or run the permutation analysis if you have not done so.")
    }
  }

  for (channel in channels) {
    for (pg in permutation_groups){
      regions <- colnames(e$permutation_p_matrix[[pg]][[channel]]$sig)
      sig_df <- e$permutation_p_matrix[[pg]][[channel]]$sig %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "sig") %>%
        tidyr::drop_na()
      pval_df <-e$permutation_p_matrix[[pg]][[channel]]$p_val %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "p_val") %>%
        tidyr::drop_na()
      test_statistic <- e$permutation_p_matrix[[pg]][[channel]]$test_statistic
      test_statistic[upper.tri(test_statistic, diag = TRUE)] <- NA
      test_statistic <- test_statistic %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "test_statistic") %>%
        tidyr::drop_na()

      permutation_results <- dplyr::left_join(test_statistic, pval_df, by = c("R1", "R2"), keep = FALSE) %>%
        dplyr::left_join(sig_df, by = c("R1", "R2"), keep = FALSE)

      if (ontology == "allen"){
        permutation_results <- permutation_results %>% dplyr::mutate(name1 = name.from.acronym(.data$R1),
                                                                     name2 = name.from.acronym(.data$R2),
                                                                     .before = test_statistic)
      } else {
        permutation_results <- permutation_results %>% dplyr::mutate(name1 = name.from.acronym.custom(.data$R1, ontology = ontology),
                                                                     name2 = name.from.acronym.custom(.data$R2, ontology = ontology),
                                                                     .before = test_statistic)
      }
      if (filter_significant) {
        permutation_results <- permutation_results %>% dplyr::filter(.data$sig)
      }
      output_dir <-  file.path(attr(e, "info")$output_path, "tables")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      if (filter_significant) {
      utils::write.csv(permutation_results, file.path(output_dir, paste0("permutation_results", "_", pg, "_", channel, "_significant.csv")), row.names = FALSE)
      } else {
        utils::write.csv(permutation_results, file.path(output_dir, paste0("permutation_results", "_", pg, "_", channel, ".csv")), row.names = FALSE)
      }
    }
  }
  return(NULL)
}


#' Export the permutation results as a csv file. This automatically saves into the tables folder. Meant for
#' datasets with hemispheric information.
#'
#' @param e experiment object
#' @param permutation_groups (str vec, default = "all") Default is to export all analyses run. Supply names of the specific analyses found under `e$permutation_p_matrix %>% names()` if you want to export specific groups.
#' @param channels (str, channels =  c("cfos", "eyfp", "colabel") The channels to process.
#' @param ontology (str, default = "allen") Set to "unified" for Kim Lab unified ontology
#' @param filter_significant (bool, default = TRUE) If FALSE, this keeps all comparisons. Otherwise exports only the most different and significant permutations.
#' @export
#' @return NULL
#' @examples
#' \dontrun{
#' e <- export_permutation_results(e, permutation_groups = "all", filter_significant = TRUE)
#' }
#'
export_permutation_results_hemi <- function(e,
                                       permutation_groups = "all",
                                       channels = c("cfos"),
                                       ontology = "allen",
                                       filter_significant =  TRUE){
  if (permutation_groups == "all") {
    permutation_groups <- e$permutation_p_matrix %>% names()
  } else{
    # Check that permutation_groups is in the names vector for e$permutation_p_matrix
    if (!all(permutation_groups %in% names(e$permutation_p_matrix))){
      stop("The permutation groups you have supplied are not all in the experiment object. Please fix the names supplied or run the permutation analysis if you have not done so.")
    }
  }

  for (channel in channels) {
    for (pg in permutation_groups){
      regions <- colnames(e$permutation_p_matrix[[pg]][[channel]]$sig)
      sig_df <- e$permutation_p_matrix[[pg]][[channel]]$sig %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "sig") %>%
        tidyr::drop_na()
      pval_df <-e$permutation_p_matrix[[pg]][[channel]]$p_val %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "p_val") %>%
        tidyr::drop_na()
      test_statistic <- e$permutation_p_matrix[[pg]][[channel]]$test_statistic
      test_statistic[upper.tri(test_statistic, diag = TRUE)] <- NA
      test_statistic <- test_statistic %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "test_statistic") %>%
        tidyr::drop_na()

      permutation_results <- dplyr::left_join(test_statistic, pval_df, by = c("R1", "R2"), keep = FALSE) %>%
        dplyr::left_join(sig_df, by = c("R1", "R2"), keep = FALSE)

      if (ontology == "allen"){
        permutation_results <- permutation_results %>% dplyr::mutate(name1 = name.from.acronym({function (x)return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}(.data$R1)),
                                                                     name2 = name.from.acronym({function (x)return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}(.data$R2)),
                                                                     .before = test_statistic)
      } else {
        permutation_results <- permutation_results %>% dplyr::mutate(name1 = name.from.acronym.custom({function (x)return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}(.data$R1), ontology = ontology),
                                                                     name2 = name.from.acronym.custom({function (x)return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}(.data$R2), ontology = ontology),
                                                                     .before = test_statistic)
      }
      if (filter_significant) {
        permutation_results <- permutation_results %>% dplyr::filter(.data$sig)
      }
      # Create table directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "tables")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }

      if (filter_significant) {
      utils::write.csv(permutation_results, file.path(output_dir, paste0("permutation_results", "_", pg, "_", channel, "_significant.csv")), row.names = FALSE)
      } else {
        utils::write.csv(permutation_results, file.path(output_dir, paste0("permutation_results", "_", pg, "_", channel, ".csv")), row.names = FALSE)
      }
    }
  }
  return(NULL)
}

#' Create a permutation difference network
#' This function requires [SMARTTR::export_permutation_results()] or  [SMARTTR::export_permutation_results_hemi()] to have been run for the specific permutation analysis of interest, with the `filter_significant`
#' parameter set to FALSE, so that users can set their own p-value threshold for filtering out edges.
#' @param e experiment object
#' @param permutation_group (str) Supply names of the specific analyses found under `e$permutation_p_matrix %>% names()`
#' @param channel (str, default = "cfos")
#' @param ontology (str, default = "allen") Set to "unified" for Kim Lab unified ontology
#' @param p_value_thresh (numeric, default = 0.001) User chosen p-value for thresholding edges to include in the network.
#' @param anatomical.order (vec, default = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"))
#' The default super region acronym list that groups all subregions in the dataset.
#' @param community_detection_function (function, default = tidygraph::group_fast_greedy)
#' \href{https://tidygraph.data-imaginist.com/reference/group_graph.html}{tidygraph grouping function} supporting
#' various algorithms for community detection
#' @param ... Additional arguments for the community_detection_function chosen
#' @return e experiment object
#' @export
#' @examples
#' \dontrun{
#' e <- create_perm_diff_network(e, permutation_group = "Control_vs_Isoflurane", channel = "cfos", ontology = "unified", p_value_thresh = 0.001, weights = NULL)
#' }
#'
create_perm_diff_network <- function(e,
                                     permutation_group,
                                     channel = "cfos",
                                     p_value_thresh = 0.001,
                                     ontology = "allen",
                                     anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                                     community_detection_function = tidygraph::group_fast_greedy,
                                     ...){
  output_dir <-  file.path(attr(e, "info")$output_path, "tables")
  perm_df <- read.csv(file.path(output_dir, paste0("permutation_results", "_", permutation_group, "_", channel, ".csv")))
  perm_df <- perm_df %>% dplyr::filter(.data$p_val <= p_value_thresh)
  perm_df <- perm_df %>% dplyr::mutate(sign = ifelse(.data$test_statistic >= 0, "pos", "neg"),
                                       weight = abs(.data$test_statistic)) %>%
    dplyr::rename(from = .data$R1,
                  to = .data$R2,
                  p.value = .data$p_val)

  perm_network <- tidygraph::tbl_graph(
    edges = perm_df,
    directed = FALSE)
  perm_network <- perm_network %>% activate(nodes) %>%
    dplyr::mutate(degree = tidygraph::centrality_degree()) %>%
    dplyr::filter(.data$degree>0)
  # Get unique common regions
  acronyms <- perm_network %>% activate(nodes) %>% tibble::as_tibble()
  acronyms <- acronyms$name
  # get the parent super region
  super.region <- acronyms %>% strip_hemi_suffix()
  if (tolower(ontology) == "allen"){
    for (sup.region in anatomical.order){
      super.region[super.region %in% get.sub.structure(sup.region)] <- sup.region
    }
  } else {
    for (sup.region in anatomical.order){
      super.region[super.region %in% get.sub.structure.custom(sup.region, ontology = ontology)] <- sup.region
    }
  }
  perm_network <- perm_network %>% activate(nodes) %>%
    dplyr::mutate(super.region =  factor(super.region, levels = intersect(anatomical.order, unique(super.region))),
                  btw = tidygraph::centrality_betweenness(weights = NULL, directed = FALSE),
                  hub_score = tidygraph::centrality_hub(scale = TRUE))  # Kleinberg's hubness centrality score
  perm_network <- perm_network %>% activate(nodes) %>%
    mutate(community = as.factor(do.call(community_detection_function, list(...))))
  if(is.null(e$permutation_diff_network)){
    e$permutation_diff_network <- list()
  }
  if(is.null(e$permutation_diff_network[[permutation_group]])){
    e$permutation_diff_network[[permutation_group]] <- list()
  }
  e$permutation_diff_network[[permutation_group]][[channel]] <- perm_network
  return(e)
}

#' Create graph objects for plotting different analysis subgroups.
#'
#' @param e experiment object
#' @param correlation_list_name (str) Name of the correlation list object used to generate the networks.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param proportional_thresh (float, default = NULL) Takes precedent over the `alph`a and the `pearson_thresh` parameters. Input the desired edge proportion (i.e., edge density) as your desired sparsity constraint.
#' @param alpha (float, default = 0.05) The significance threshold for including brain regions in the network. if NULL or NA,
#' this threshold is not applied.
#' @param pearson_thresh (float, default = 0) The pearson correlation coefficient threshold to apply for filtering out
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param filter_isolates (logical, default = FALSE) Whether to filter out the number of isolated (zero degree) nodes from the network. Default is to retain them.
#' @param anatomical.order (vec, c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The default super region acronym list that groups all subregions in the dataset.
#'
#' @return e experiment object. This object now has a new added element called `networks.` This is a list storing a
#' graph object per channel for each network analysis run.
#' The name of each network (`network_name`) is the same as the `correlation_list_name`
#' used to generate the network. This `network_name` is fed as a parameter into the
#' [SMARTTR::plot_networks()] function.
#' @export
#' @examples
#' \dontrun{
#' e sundowning <- create_networks(sundowning, correlation_list_name = "female_control", alpha = 0.05)
#' }
#' @seealso [SMARTTR::plot_networks()]

create_networks <- function(e,
                            correlation_list_name,
                            channels = c("cfos", "eyfp", "colabel"),
                            proportional_thresh = NULL,
                            alpha = 0.05,
                            pearson_thresh = 0,
                            ontology = "allen",
                            anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                            filter_isolates = FALSE){

  networks <- vector(mode = "list", length = length(channels))
  names(networks) <- channels
  for (channel in channels){

    #__________________ Create Edge Tibble _____________
    corr_df <-  e$correlation_list[[correlation_list_name]][[channel]] %>% purrr::map(tibble::as_tibble)
    val_names <- names(corr_df)

    for (k in 1:length(corr_df)){
      corr_df[[k]] <- tibble::add_column(corr_df[[k]], row_acronym = names(corr_df[[k]]), .before = TRUE) %>%
        tidyr::pivot_longer(!any_of("row_acronym"), names_to = "col_acronym", values_to = val_names[k])
    }
    df <- corr_df$r %>%
      dplyr::left_join(corr_df$n, by = c("row_acronym", "col_acronym")) %>%
      dplyr::left_join(corr_df$P, by = c("row_acronym", "col_acronym")) %>%
      dplyr::left_join(corr_df$sig, by = c("row_acronym", "col_acronym")) %>%
      tidyr::drop_na()
    edges <- df %>% dplyr::mutate(sign = ifelse(.data$r >= 0, "pos", "neg"),
                                  weight = abs(.data$r)) %>%
       dplyr::rename(from = any_of("row_acronym"),
                     to = any_of("col_acronym"),
                     p.value = any_of("P"))

    # ___________Create Node tibble ________________
    acronyms <-  e$correlation_list[[correlation_list_name]][[channel]]$r %>% rownames()
    super.region <- acronyms %>% strip_hemi_suffix()
    if (tolower(ontology) == "allen"){
      for (sup.region in anatomical.order){
        super.region[super.region %in% get.sub.structure(sup.region)] <- sup.region
      }
    } else {
      for (sup.region in anatomical.order){
        super.region[super.region %in% get.sub.structure.custom(sup.region, ontology = ontology)] <- sup.region
      }
    }
    nodes <- tibble::tibble(name = acronyms, super.region = super.region)

    # _____________ Create the network ________________
    network <- tidygraph::tbl_graph(nodes = nodes,
                                    edges = edges,
                                    directed = FALSE)
    network <- network %>%
      tidygraph::convert(tidygraph::to_simple) %>%
      activate(edges) %>%
      mutate(r = purrr::map_dbl(.orig_data, ~ .x[1,]$r),
             n = purrr::map_int(.orig_data, ~ .x[1,]$n),
             p.value = purrr::map_dbl(.orig_data, ~ .x[1,]$p.value),
             sig = purrr::map_lgl(.orig_data, ~ .x[1,]$sig),
             sign = purrr::map_chr(.orig_data, ~ .x[1,]$sign),
             weight = purrr::map_dbl(.orig_data, ~ .x[1,]$weight)) %>%
      tidygraph::select(-.tidygraph_edge_index,  -.orig_data) %>%
      activate(nodes) %>%
      tidygraph::select(-.tidygraph_node_index)

    if (!is.na(proportional_thresh) && !is.null(proportional_thresh)){
      edges <- network %>% activate(edges) %>% dplyr::arrange(dplyr::desc(.data$weight)) %>% dplyr::as_tibble()
      row <- round(igraph::gsize(network)*proportional_thresh)
      pearson_thresh <- edges$weight[[row]]
      network <- network %>% activate(edges) %>% dplyr::filter(.data$weight >= pearson_thresh)

    } else{
      if(!is.na(alpha) && !is.null(alpha)){
        network <- network %>% activate(edges) %>% dplyr::filter(.data$p.value < alpha)
      }
      if(!is.na(pearson_thresh) && !is.null(pearson_thresh)){
        network <- network %>% activate(edges) %>% dplyr::filter(.data$weight > pearson_thresh)
      }
    }

    network <- network %>% activate(nodes) %>%
      dplyr::mutate(super.region = factor(.data$super.region, levels = unique(.data$super.region)),
                    degree = tidygraph::centrality_degree(),
                    triangles = igraph::count_triangles(.),
                    # clust.coef = ifelse(degree < 2, 0 , 2*triangles/(degree*(degree - 1))),
                    clust.coef = igraph::transitivity(., type = "local", isolates = "zero"))
    if (filter_isolates){
      network <- network %>% activate(nodes) %>%
        dplyr::filter(.data$degree > 0)
    }
    # Add distance, efficiency, and btw metrics
    d <- igraph::distances(network, weights = NA)
    d[which(!is.finite(d))] <- NA
    diag(d) <- NA
    network <- network %>%
      dplyr::mutate(avg.dist = rowSums(d, na.rm = TRUE)/(n()-1),
             efficiency = rowSums(d^(-1),na.rm = TRUE)/(n()-1),
             btw = tidygraph::centrality_betweenness(weights = NULL, directed = FALSE),
             group = as.factor(tidygraph::group_walktrap()))
    networks[[channel]] <- network
  }
  e$networks[[correlation_list_name]] <- networks
  return(e)
}

#' Summarize multiple networks.
#' calculate network statistics for each network. This is not meant to summarize networks created using `create_joined_networks`.
#' @param e experiment object
#' @param network_names (str) The names of the networks to generate summary tables for, e.g. network_names = c("female_AD", "female_control")
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param save_stats (bool, default = TRUE) Save the summary stats as a csv file in the output folder. Note that the clustering calculated is an average of the local vertex clustering.
#' @param save_degree_distribution (bool, default = TRUE) Save the network degree distributions (frequencies of each degree) across each comparison group as a csv file.
#' @param save_betweenness_distribution (bool, default = TRUE) Save the betweenness distribution and summary as a csv.
#' @param save_efficiency_distribution (bool, default = TRUE) Save the efficiency distribution and summary as a csv.
#' @return e experiment object
#' @export
#' @examples
#' \dontrun{
#' e <- get_network_statistics(e,  network_names = c("female_AD", "female_control"),
#' channels = c("cfos", "eyfp", "colabel"), save_stats = TRUE, save_degree_distribution = TRUE)
#' }
summarise_networks <- function(e,
                               network_names,
                               channels = c("cfos", "eyfp", "colabel"),
                               save_stats = TRUE,
                               save_degree_distribution = TRUE,
                               save_betweenness_distribution = TRUE,
                               save_efficiency_distribution = TRUE){

  # Check all the networks named are made
  if(!all(network_names %in% names(e$networks))){
    stop("The network names that you specified are not all stored in your experiment object.",
         "Please check the names of the networks you want to summarize!")
  }

  # Check that they all contain the channels to compare with each other
  for (network in network_names){
    if(!all(channels %in% names(e$networks[[network]]))){
      stop("At least one of your networks does not contain data for the channel(s) you want to summarize.",
           "Please double check and change your channels parameter.")
    }
  }

  for (channel in channels){
    nodes_df_list <- list()
    edges_df_list <- list()
    for (network in network_names){
      nodes_df_list[[network]] <-  e$networks[[network]][[channel]] %>% tidygraph::activate(nodes) %>%
        tibble::as_tibble() %>% dplyr::mutate(group = network)
      edges_df_list[[network]] <-  e$networks[[network]][[channel]] %>% tidygraph::activate(edges) %>%
        tibble::as_tibble() %>% dplyr::mutate(group = network)
    }
    nodes_df <- dplyr::bind_rows(nodes_df_list) %>% dplyr::mutate(group=factor(.data$group, levels = network_names))
    edges_df <- dplyr::bind_rows(edges_df_list) %>% dplyr::mutate(group=factor(.data$group, levels = network_names))

    network_stats_df <- nodes_df %>% dplyr::group_by_at("group") %>%
      dplyr::summarise_if(is.numeric,list(~mean(.),~sd(.), ~sem(.))) %>%
      dplyr::left_join(nodes_df %>% dplyr::group_by_at("group") %>% dplyr::summarise(n.nodes = n())) %>%
      dplyr::left_join(edges_df %>% dplyr::group_by_at("group") %>% dplyr::summarise(n.edges = n())) %>%
      dplyr::mutate(edge.density = .data$n.edges/(.data$n.nodes*(.data$n.nodes-1)/2)) %>%
      dplyr::rename_all(stringr::str_replace,pattern = "_",replacement=".")

    #make table of degree frequency (useful for plotting degree histogram outside of R)
    degree_distribution_df <- nodes_df %>% dplyr::group_by_at(c("group", "degree")) %>% dplyr::count(.data$degree)

    # Make a dataframe of degree (useful)
    degree_distribution <- nodes_df %>% dplyr::select(any_of(c("name", "group", "degree"))) %>%
      dplyr::group_by_at("group") %>% dplyr::arrange(.data$group, dplyr::desc(.data$degree)) %>%
      dplyr::mutate(name = factor(.data$name, levels=.data$name))

    #make table of betweenness frequency (useful for plotting degree histogram outside of R)
    betweenness_distribution_df <- nodes_df %>% dplyr::group_by_at(c("group", "btw")) %>% dplyr::count(.data$btw)

    # Make a dataframe of betweenness (useful)
    betweenness_distribution <- nodes_df %>% dplyr::select(any_of(c("name", "group", "btw"))) %>%
      dplyr::group_by_at("group") %>% dplyr::arrange(.data$group, dplyr::desc(.data$btw)) %>%
      dplyr::mutate(name = factor(.data$name, levels=.data$name))

    #make table of efficiency frequency (useful for plotting degree histogram outside of R)
    efficiency_distribution_df <- nodes_df %>% dplyr::group_by_at(c("group", "efficiency")) %>% dplyr::count(.data$efficiency)

    # Make a dataframe of efficiency (useful)
    efficiency_distribution <- nodes_df %>% dplyr::select(any_of(c("name", "group", "efficiency"))) %>%
      dplyr::group_by_at("group") %>% dplyr::arrange(.data$group, dplyr::desc(.data$efficiency)) %>%
      dplyr::mutate(name = factor(.data$name, levels=.data$name))

    # Create table directory if it doesn't already exists
    output_dir <-  file.path(attr(e, "info")$output_path, "tables")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }

    if(save_stats){
      write.csv(network_stats_df,  file.path(output_dir, paste0("summary_networks_stats_", channel, ".csv")))
    }

    if(save_degree_distribution){
      write.csv(degree_distribution,  file.path(output_dir, paste0("networks_degree_distributions_", channel, ".csv")))
      write.csv(degree_distribution_df,  file.path(output_dir, paste0("networks_degree_distributions_summary_", channel, ".csv")))
    }

    if(save_betweenness_distribution){
      write.csv(betweenness_distribution,  file.path(output_dir, paste0("networks_betweenness_distributions_", channel, ".csv")))
      write.csv( betweenness_distribution_df,  file.path(output_dir, paste0("networks_betweenness_distributions_summary_", channel, ".csv")))
    }

    if(save_efficiency_distribution){
      write.csv(efficiency_distribution,  file.path(output_dir, paste0("networks_efficiency_distributions_", channel, ".csv")))
      write.csv(efficiency_distribution_df,  file.path(output_dir, paste0("networks_efficiency_distributions_summary_", channel, ".csv")))
    }

    # Store the network summary data into channels
    e$networks_summaries[[channel]] <- list(networks_nodes = nodes_df,
         networks_edges = edges_df,
         networks_stats = network_stats_df,
         networks_degree_distrib = degree_distribution_df)
  }
  return(e)
}



#' Create a joined network to visualize overlapping connections with the same outer joined node set.
#'
#' @param e experiment object
#' @param correlation_list_names (str vec) character vector of the two correlation lists used to include in a joined network
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param proportional_thresh (float, default = NULL) Takes precedent over the `alpha` and the `pearson_thresh` parameters. Input the desired edge proportion (i.e., edge density) as your desired sparsity constraint.
#' @param alpha (float, default = 0.05) The significance threshold for including brain regions in the network. if NULL or NA,
#' this threshold is not applied.
#' @param pearson_thresh (float, default = 0.8) The pearson correlation coefficient threshold to apply for filtering out
#' @param proportional_thresh2 (NULL) If not NULL, this gives the option of filtering the second network by a different proportional threshold from the first.
#' @param alpha2 (NULL) If not NULL, this gives the option of filtering the second network by a different alpha from the first. The `alpha` parameter will then be used as the threshold for network 1.
#' @param pearson_thresh2 (NULL) If not NULL, this gives the option of filtering the second network by a different pearson threshold from the first network.
#' The `pearson_thresh` parameter will then be used as the threshold for network 1.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param filter_isolates (logical, default = TRUE) Whether to filter out the number of isolated (zero degree) nodes from the network. Default is to retain them.
#' @param export_overlapping_edges (bool, default  = TRUE) Whether to export the overlapping edges between the two networks as a csv into the `table` directory.
#' @param anatomical.order (vec, c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The default super region acronym list that groups all subregions in the dataset.
#' @param proportional_sort
#'
#' @return e experiment object. This object now has a new added element called `networks.` This is a list storing a
#' graph object per channel for each network analysis run. The name of each network (`network_name`) is the same as the `correlation_list_name`
#' used to generate the network. This `network_name` is fed as a parameter into the
#' [SMARTTR::plot_networks()] function.
#' @export
#' @examples
#' \dontrun{sundowning <- create_joined_networks(sundowning,
#' correlation_list_name = "female_control", alpha = 0.05)
#' }
#' @seealso [SMARTTR::plot_networks()]

create_joined_networks <- function(e,
                                   correlation_list_names = c("male_agg", "female_non"),
                                   channels = "cfos",
                                   ontology = "unified",
                                   alpha = 0.001,
                                   pearson_thresh = 0.9,
                                   proportional_thresh = NULL,
                                   proportional_sort = "pearson",
                                   alpha2 = NULL,
                                   pearson_thresh2 = NULL,
                                   proportional_thresh2 = NULL,
                                   filter_isolates = TRUE,
                                   anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                                   export_overlapping_edges = TRUE){

  # List to store the networks
  networks <- vector(mode = "list", length = length(channels))
  nodes_joined <- vector(mode = "list")
  edges_joined <- vector(mode = "list")
  names(networks) <- channels

  for (channel in channels){
    for (correlation_list_name in correlation_list_names){
      #__________________ Create Edge Tibble _____________
      corr_df <-  e$correlation_list[[correlation_list_name]][[channel]] %>% purrr::map(tibble::as_tibble)
      val_names <- names(corr_df)

      for (k in 1:length(corr_df)){
        corr_df[[k]] <- tibble::add_column(corr_df[[k]], row_acronym = names(corr_df[[k]]), .before = TRUE) %>%
          tidyr::pivot_longer(!any_of("row_acronym"), names_to = "col_acronym", values_to = val_names[k])
      }

      # Combined the correlation plots into one dataframe
      df <- corr_df$r %>% dplyr::left_join(corr_df$n, by = c("row_acronym", "col_acronym")) %>%
        dplyr::left_join(corr_df$P, by = c("row_acronym", "col_acronym")) %>%
        dplyr::left_join(corr_df$sig, by = c("row_acronym", "col_acronym")) %>%
        tidyr::drop_na()

      # Edges dataframe
      edges_joined[[correlation_list_name]] <- df %>% dplyr::mutate(sign = ifelse(.data$r >= 0, "pos", "neg"),
                                                                    weight = abs(.data$r)) %>%
        dplyr::rename(from = any_of("row_acronym"),
                      to = any_of("col_acronym"),
                      p.value = any_of("P")) %>%
        tibble::add_column(network = correlation_list_name)

      acronyms <-  e$correlation_list[[correlation_list_name]][[channel]]$r %>% rownames()
      # get the parent super region
      super.region <- acronyms %>% strip_hemi_suffix()

      if (tolower(ontology) == "allen"){
        for (sup.region in anatomical.order){
          super.region[super.region %in% get.sub.structure(sup.region)] <- sup.region
        }
      } else {
        for (sup.region in anatomical.order){
          super.region[super.region %in% get.sub.structure.custom(sup.region, ontology = ontology)] <- sup.region
        }
      }
      nodes_joined[[correlation_list_name]] <- tibble::tibble(name = acronyms, super.region = super.region)
    }

    network1 <- tidygraph::tbl_graph(nodes = nodes_joined[[1]],
                                     edges = edges_joined[[1]],
                                     directed = FALSE)

    network2 <- tidygraph::tbl_graph(nodes = nodes_joined[[2]],
                                     edges = edges_joined[[2]],
                                     directed = FALSE)
    network1 <- network1 %>%
      tidygraph::convert(tidygraph::to_simple) %>%
      tidygraph::activate(edges) %>%
      mutate(r = purrr::map_dbl(.orig_data, ~ .x[1,]$r),
             n = purrr::map_int(.orig_data, ~ .x[1,]$n),
             p.value = purrr::map_dbl(.orig_data, ~ .x[1,]$p.value),
             sig = purrr::map_lgl(.orig_data, ~ .x[1,]$sig),
             sign = purrr::map_chr(.orig_data, ~ .x[1,]$sign),
             weight = purrr::map_dbl(.orig_data, ~ .x[1,]$weight),
             network = purrr::map_chr(.orig_data, ~ .x[1,]$network)) %>%
      tidygraph::select(-.tidygraph_edge_index,  -.orig_data) %>%
      tidygraph::activate(nodes) %>%
      tidygraph::select(-.tidygraph_node_index)

    network2 <- network2 %>%
      tidygraph::convert(tidygraph::to_simple) %>%
      tidygraph::activate(edges) %>%
      mutate(r = purrr::map_dbl(.orig_data, ~ .x[1,]$r),
             n = purrr::map_int(.orig_data, ~ .x[1,]$n),
             p.value = purrr::map_dbl(.orig_data, ~ .x[1,]$p.value),
             sig = purrr::map_lgl(.orig_data, ~ .x[1,]$sig),
             sign = purrr::map_chr(.orig_data, ~ .x[1,]$sign),
             weight = purrr::map_dbl(.orig_data, ~ .x[1,]$weight),
             network = purrr::map_chr(.orig_data, ~ .x[1,]$network)) %>%
      tidygraph::select(-.tidygraph_edge_index,  -.orig_data) %>%
      tidygraph::activate(nodes) %>%
      tidygraph::select(-.tidygraph_node_index)

    if (!is.na(proportional_thresh) && !is.null(proportional_thresh)){

      if ("pearson" == tolower(proportional_sort)){

        edges <- network1 %>% activate(edges) %>% dplyr::arrange(dplyr::desc(.data$weight)) %>% dplyr::as_tibble()
        row <- round(igraph::gsize(network1)*proportional_thresh)
        pearson_thresh <- edges$weight[[row]]
        network1 <- network1 %>% activate(edges) %>% dplyr::filter(.data$p.value >= pearson_thresh)
        print(edges$p.value[[row]])

        # network 2
        edges <- network2 %>% activate(edges) %>% dplyr::arrange(dplyr::desc(.data$weight)) %>% dplyr::as_tibble()
        if (!is.na(proportional_thresh2) && !is.null(proportional_thresh2)){
          row <- round(igraph::gsize(network2)*proportional_thresh2)
          } else{
          row <- round(igraph::gsize(network2)*proportional_thresh)
          }
        pearson_thresh <- edges$weight[[row]]
        network2 <- network2 %>% activate(edges) %>% dplyr::filter(.data$weight >= pearson_thresh)
        print(edges$p.value[[row]])

      } else if ("alpha" == tolower(proportional_sort)) {
        edges <- network1 %>% activate(edges) %>% dplyr::arrange(.data$p.value) %>% dplyr::as_tibble()
        row <- round(igraph::gsize(network1)*proportional_thresh)
        alpha_thresh <- edges$p.value[[row]]
        network1 <- network1 %>% activate(edges) %>% dplyr::filter(.data$p.value <= alpha_thresh)
        # print(edges$p.value)
        print(alpha_thresh)

        edges <- network2 %>% activate(edges) %>% dplyr::arrange(.data$p.value) %>% dplyr::as_tibble()
        if (!is.na(proportional_thresh2) && !is.null(proportional_thresh2)){
          row <- round(igraph::gsize(network2)*proportional_thresh2)
        } else{
          row <- round(igraph::gsize(network2)*proportional_thresh)
        }

        alpha_thresh <- edges$p.value[[row]]
        network2 <-network2 %>% activate(edges) %>% dplyr::filter(.data$p.value <= alpha_thresh)
        print(alpha_thresh)
    }
      } else{
      # filter by alpha
      if(!is.na(alpha) && !is.null(alpha)){
        network1 <- network1 %>% tidygraph::activate(edges) %>% dplyr::filter(.data$p.value < alpha)
        if(!is.na(alpha2) && !is.null(alpha2)){
          network2 <- network2 %>% tidygraph::activate(edges) %>% dplyr::filter(.data$p.value < alpha2)
        } else{
          network2 <- network2 %>% tidygraph::activate(edges) %>% dplyr::filter(.data$p.value < alpha)
        }
      }

      if(!is.na(pearson_thresh) && !is.null(pearson_thresh)){
        network1 <- network1 %>% tidygraph::activate(edges) %>% dplyr::filter(.data$weight > pearson_thresh)
        if (!is.na(pearson_thresh2) && !is.null(pearson_thresh2)){
          network2 <- network2 %>% tidygraph::activate(edges) %>% dplyr::filter(.data$weight > pearson_thresh2)
        } else {
          network2 <- network2 %>% tidygraph::activate(edges) %>% dplyr::filter(.data$weight > pearson_thresh)
        }
      }
    }

    network <- tidygraph::graph_join(network1, network2)
    network <- network %>% tidygraph::activate(nodes) %>%
      dplyr::mutate(super.region = factor(.data$super.region, levels = unique(.data$super.region)),
                    degree = tidygraph::centrality_degree(mode="all"))

    if (filter_isolates){
      network <- network %>% activate(nodes) %>%
        dplyr::filter(.data$degree > 0)
    }

    networks[[channel]] <- network

    if (export_overlapping_edges) {
      joined_network_name <- paste(correlation_list_names, collapse = "_")
      networks[[channel]]%>% tidygraph::activate(nodes) %>% tibble::as_tibble() -> nodes_df
      networks[[channel]] %>% tidygraph::activate(edges) %>% tibble::as_tibble() -> edges_df
      edges_df <- edges_df %>% mutate(from = nodes_df$name[edges_df$from],
                                      to =   nodes_df$name[edges_df$to])
      edges_df_p1 <- edges_df %>% dplyr::filter(.data$network == correlation_list_names[1])
      edges_df_p2 <- edges_df %>% dplyr::filter(.data$network == correlation_list_names[2])
      mutual_edges <- edges_df_p1 %>% dplyr::inner_join(edges_df_p2, by = c("from", "to"),
                                                        suffix = c(paste0(".", correlation_list_names[1]), paste0(".", correlation_list_names[2])))

      output_dir <-  file.path(attr(e, "info")$output_path, "tables")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      file_path <- file.path(output_dir, paste0("joined_networks_mutual_edges_",joined_network_name, "_", channel,".csv"))
      utils::write.csv(mutual_edges, file = file_path)
      message(paste0("Saved mutual edge list at: ", file_path))
    }
  }

  e$networks[[paste(correlation_list_names, collapse = "_")]] <- networks
  return(e)
}

#' Implement rewiring algorithms to current empirical networks to randomize certain network properties.
#'
#' Not that this keeps other characteristics constant (such as preserved degree sequence). These null networks can them be used to compare
#' against and normalize the empirical networks. Currently this essentially erases edge metrics and treats networks like binary graphs.
#' Edge weights are not used in calculating
#' network topology metrics.
#' @param e experiment object
#' @param network_name (str) Name of the network
#' @param channels (str)  Vector of channels to process
#' @param method (str, default = "ms") "ms" implements Maslov-Sneppen rewiring approach (annuls all network properties except for network size, connection density, and degree distribution).
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param n_rewires (int, default = 10000) The number of rewires for randomization for "ms" rewiring implementation. Recommended to be the larger of either 10,000 or 10*No. edges in a graph.
#' @param n_networks (int, default = 100) The number of random networks to create
#' @param seed (int, default = 5) Random seed for future replication.
#' @param return_graphs (logical, default = FALSE) if TRUE, returns a list organized by channel containing a sublist, with each element containing a tidygraph object. This must be FALSE if you want to run
#' you want to summarize the null network statistics with [SMARTTR::summarize_null_networks()]
#' @return Summary table of rewired network properties of all nodes showing the average of all randomized network properties generated.
#' @export
#' @examples
#' \dontrun{
#' summary_table <- rewire_network(e, network_name = "network1", channels = "cfos",
#' n_rewire =  igraph::gsize(e$networks$network1$cfos)*100, n_networks = 100)
#' }
rewire_network <- function(e,
                           network_name,
                           channels = "cfos",
                           method = "ms",
                           ontology = "allen",
                           n_rewires = 10000,
                           n_networks = 100,
                           return_graphs = FALSE,
                           seed = 5){

  # null_graphs <- vector(mode = "list", length = length(channels))
  null_nodes <- vector(mode = "list", length = length(channels))
  names(null_nodes) <- channels

  if (return_graphs){
    null_graphs <- vector(mode = "list", length = length(channels))
    names(null_graphs) <- channels
  }

  set.seed(seed)

  # Network check
  if (is.null(e$networks[[network_name]])){
    stop("The network_name you provided does not exist. Please run create_networks() before implementing this function.")
  }

  for (channel in channels) {
    #  Channel check
    if (is.null(e$networks[[network_name]][[channel]])){
      stop("The network for this channel does not exist. Please run create_networks() before implementing this function")
    }

    if (method == "ms") {
      g <- e$networks[[network_name]][[channel]] %>% maslov_sneppen_rewire(n_rewires = n_rewires)
      g <- g %>% activate(nodes) %>% dplyr::mutate(iter = 1) %>% dplyr::select(-any_of("group"))
      null_nodes[[channel]] <- g %>% dplyr::as_tibble()

      if (return_graphs){
        graph_list <- vector(mode = "list", length = n_networks)
        graph_list[[1]] <- g
      }

      for (i in 2:n_networks){
        g <- e$networks[[network_name]][[channel]] %>% maslov_sneppen_rewire(n_rewires = n_rewires)
        g <- g %>% activate(nodes) %>% dplyr::mutate(iter = i) %>% dplyr::select(-any_of("group"))
        null_node <- g  %>% dplyr::as_tibble()
        null_nodes[[channel]] <- dplyr::bind_rows(null_nodes[[channel]], null_node)

        if (return_graphs){
          graph_list[[i]] <- g
        }
      }
    }
    if (return_graphs) {
      null_graphs[[channel]] <- graph_list
    }
  }
  if (return_graphs){
    return(null_graphs)
  } else{
    return(null_nodes)
  }
}


#' Summarize the parameters of the rewired null networks generated by [SMARTTR::rewire_network()]
#'
#' @param null_nodes_list a list of output summary tables (1 per network) of rewired network properties of all nodes from [SMARTTR::rewire_network()].
#' @param network_names (str vec) Name of the networks that were rewired (in the same order as the list). `null_nodes_list` and `network_names` must have
#' the same length.
#' @param channel (str)   channel to process
#'
#' @return a list of length 2. The first element is named `global_summary` and contains a table of global summary statistics.
#' The second element is named `node_summary`, and contained per node statistics averaged from multiple null networks.
#' @export
#' @examples
#' \dontrun{
#' rewire_summary <- rewire_network(e, "Context", channels = "eyfp", return_graphs = FALSE)
#' summarized_null_networks <- summarize_null_networks(rewire_summary,
#' network_names = "Context", channel = "eyfp")
#' }
summarize_null_networks <- function(null_nodes_list,
                                    network_names = NULL,
                                    channel = "cfos"){
  summaries <- vector(mode="list", length = 2)
  names(summaries) <- c("node_summary", "global_summary")
  for (n in 1:length(null_nodes_list)){
    if (n == 1){
    null_node_summary <- null_nodes_list[[n]][[channel]] %>% dplyr::group_by_at(c("name", "super.region")) %>%
      dplyr::summarise(null_degree = mean(.data$degree, na.rm = TRUE),
                null_triangles = mean(.data$triangles, na.rm = TRUE),
                null_clust.coef = mean(.data$clust.coef, na.rm = TRUE),
                null_dist = mean(.data$avg.dist, na.rm = TRUE),
                null_efficiency = mean(.data$efficiency, na.rm = TRUE),
                null_btw = mean(.data$btw, na.rm = TRUE)) %>%
      dplyr::mutate(group = network_names[[n]], .before = any_of("null_degree"))
    } else {
      toadd <- null_nodes_list[[n]][[channel]] %>% dplyr::group_by_at(c("name", "super.region")) %>%
        summarise(null_degree = mean(.data$degree, na.rm = TRUE),
                  null_triangles = mean(.data$triangles, na.rm = TRUE),
                  null_clust.coef = mean(.data$clust.coef, na.rm = TRUE),
                  null_dist = mean(.data$avg.dist, na.rm = TRUE),
                  null_efficiency = mean(.data$efficiency, na.rm = TRUE),
                  null_btw = mean(.data$btw, na.rm = TRUE)) %>%
        dplyr::mutate(group = network_names[[n]], .before = any_of("null_degree"))
      null_node_summary <- null_node_summary  %>% dplyr::bind_rows(toadd)
    }
  }
  ## Global level normalization
  null_global_summary <- null_node_summary %>% dplyr::ungroup() %>%  dplyr::group_by_at("group") %>%
    dplyr::summarise_if(is.numeric,list(~mean(.),~sd(.), ~sem(.)))
  summaries[["node_summary"]] <- null_node_summary
  summaries[["global_summary"]] <- null_global_summary
  return(summaries)
}


#________________ Internal Analysis functions _______________

#' @noRd
filter_df_by_char_params <- function(df, by, values){
  for (k in 1:length(by)){
    # Convert the variable name into a symbol
    var <- rlang::sym(by[k])
    df <- df %>% dplyr::filter(!!var == values[k])
    # df <- df %>% dplyr::filter(dplyr::all_of(by[k]) == values[k])
  }
  return(df)
}

# Sort the dataframes columns to be in anatomical order
#' @noRd
sort_anatomical_order <- function(common_regions, ontology ="allen",
                                  anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")){

  if (tolower(ontology) == "allen"){
    common_regions <- anatomical.order %>% purrr::map(get.sub.structure) %>%
      purrr::map(intersect,y=common_regions) %>% unlist()
  } else {
    common_regions <- anatomical.order %>% purrr::map(get.sub.structure.custom, ontology=ontology) %>%
      purrr::map(intersect,y=common_regions) %>% unlist()
  }
  return(common_regions)
}



#' Internal algorithm for maslov-sneppen rewiring
#'
#' This algorithm is not appropriate if you would like to take weights into account
#' See Maslov & Sneppen (2002)"Specificity and stability in topology of protein networks"
#' @param network tidygraph graph object
#' @param n_rewires number of rewires.
#' @return a tidygraph graph
#' @noRd

maslov_sneppen_rewire <- function(network, n_rewires = 10000){

  g <- network %>%
    igraph::rewire(igraph::keeping_degseq(loops = FALSE, niter = n_rewires)) %>%
    tidygraph::as_tbl_graph()
  g <-  g %>% tidygraph::activate(edges) %>% dplyr::select(any_of(c("from", "to")))
  # Recalculate the nodal metrics
  g <- g %>% tidygraph::activate(nodes) %>% dplyr::mutate(degree = tidygraph::centrality_degree(),
                                                          triangles = igraph::count_triangles(.),
                                                          clust.coef = igraph::transitivity(., type = "local", isolates = "zero"))

  # Add distance, efficiency, and btw metrics
  d <- igraph::distances(g, weights = NA)
  d[which(!is.finite(d))] <- NA
  diag(d) <- NA
  g <- g %>%
    dplyr::mutate(avg.dist = rowSums(d, na.rm = TRUE)/(n()-1),
                  efficiency = rowSums(d^(-1),na.rm = TRUE)/(n()-1),
                  btw = tidygraph::centrality_betweenness(weights = NULL,
                                                          directed = FALSE))
  return(g)

}


#' Generate array of null distribution of region pairwise correlation differences.
#' @param df dataframe
#' @param correlation_list_name_1 permutation group 1
#' @param correlation_list_name_2 permutation group 2
#' @param n_shuffle number of shuffles
#' @param seed random seed for replication
#' @param method (str, default = "pearson", options = c("pearson", "spearman")) Specifies the type of correlations to compute.
#' Spearman correlations are the Pearson linear correlations computed on the ranks of non-missing elements, using midranks for ties. See also [Hmisc::rcorr()]
#' @param progressbar (bool, default = TRUE) Display a progress bar for the processing time of the permutation.
#' @return a matrix storing the permutation with the dimensions no. regions x no. regions x n_shuffle
#' @noRd
permute_corr_diff_distrib <- function(df, correlation_list_name_1, correlation_list_name_2,
                                      n_shuffle = n_shuffle, seed = 5, method = "pearson", progressbar = TRUE){
  set.seed(seed)
  # Create a 3D matrix to hold the correlation distributions
  df <- df %>% dplyr::select(any_of("corr_group"), dplyr::where(is.numeric))
  n_reg <- length(names(df)) - 1
  region_names <- names(df)[2:length(names(df))]
  corr_diff_matrix <- array(dim= c(n_reg, n_reg, n_shuffle))
  dimnames(corr_diff_matrix) <- list(region_names, region_names, 1:n_shuffle)

  # Get original group order
  corr_groups <- df$corr_group

  if (progressbar){
    pb <- txtProgressBar(min = 0, max = n_shuffle, initial = 0, style = 3)
  }

  for (n in 1:n_shuffle){
    # Shuffle the group labels
    df$corr_group <- sample(corr_groups, replace = FALSE)
    matrix_list <-  df %>% dplyr::group_by_at("corr_group") %>%
      dplyr::group_map(as.matrix, .keep = TRUE)
    element_1_name <- matrix_list[[1]][,"corr_group"] %>% unique()

    # Reorder the matrix in case element 1 doesn't correspond to correlation_list_name_1
    if (element_1_name != correlation_list_name_1){
      matrix_list <- list(matrix_list[[2]], matrix_list[[1]])
    }
    names(matrix_list) <- c(correlation_list_name_1, correlation_list_name_2)

    correlations_list <- vector(mode = "list", length = 2)
    names(correlations_list) <- c(correlation_list_name_1, correlation_list_name_2)
    correlations_list[[correlation_list_name_1]] <- matrix_list[[correlation_list_name_1]][,-1] %>% try_correlate(type = method)
    correlations_list[[correlation_list_name_2]] <- matrix_list[[correlation_list_name_2]][,-1] %>% try_correlate(type = method)
    corr_diff_matrix[,,n] <- correlations_list[[correlation_list_name_2]]$r - correlations_list[[correlation_list_name_1]]$r

    if (progressbar){
      setTxtProgressBar(pb,n)
    }
  }
  return(corr_diff_matrix)
}

#' Get  a list of intersecting regions to a list of common regions
#' @param common_reg A comprehensive list of all regions (and all existing subregions) that the brain area in the `rois` list will be compared against.
#' @param rois A list of rois whose regions and subregions will be compared the the `common_reg` list
#' @return common_reg A list of the rois (or any of its subregions) that intersected with the common_reg list.
#' @noRd
rois_intersect_region_list <- function(common_reg, rois){
  # Get rois of all the child regions
  child_r <- get.acronym.child(rois)
  while (length(child_r) > 0){
    rois <- c(rois, child_r)
    child_r <- get.acronym.child(child_r) %>% na.omit()
  }
  message("Checking also for child regions of specified roi(s). Only rois and child rois that overlap with common
            regions of both channels will be used.")

  if (all(!rois %in% common_reg)){
    message(paste0("None of the roi(s) specified or their child regions is a common region found in both channels. ",
                   "Common regions include:\n"))
    for (reg in common_reg){
      message(reg)
    }
    stop("Set the rois argument to a subset of these acronyms.")
  } else{
    common_reg <- intersect(rois, common_reg)
  }
  return(common_reg)
}


#' Try to correlate
#' @param df_channel dataframe of normalized counts, with rows as animals and columns as regions
#' @return  a list with elements r, the matrix of correlations, n the matrix of number of observations used in analyzing each
#' pair of variables, and P, the asymptotic P-values.
#' @noRd
try_correlate <- function(df_channel, type = "pearson"){
  tryCatch({
    # Code that may throw an error
    df_corr <- df_channel %>% Hmisc::rcorr(type = type)
    return(df_corr)
  },
  error = function(err) {
    message(c("One or more of your brain regions has a n below the recommended value."))
    message("\nHere's the original error message:")
    message(err)
    message(c("\nCalculating pearsons, using alternative method but we ",
              "recommend increasing the sample size or excluding that region due to insufficient data"))

    df_channel %>% dim() -> shape
    df_corr <- vector(mode = "list")
    df_corr$r <- matrix(nrow = shape[2], ncol = shape[2])
    rownames(df_corr$r) <- colnames(df_channel)
    colnames(df_corr$r) <- colnames(df_channel)
    df_corr$n <-  df_corr$r
    df_corr$P <- df_corr$r
    for(r1 in 1:shape[2]){
      for(r2 in 1:r1){
        df_corr$n[r1,r2] <- sum(df_channel[,r1] & df_channel[,r2], na.rm = TRUE)
        ct <- stats::cor.test(df_channel[,r1], df_channel[,r2], method = type, exact = FALSE)
        df_corr$P[r1,r2] <- ct$p.value
        df_corr$r[r1,r2] <- ct$estimate
      }
    }
    # df_corr is lower triangle.
    # need to mirror to upper tri
    upper_tri <- upper.tri(df_corr$r)
    df_corr$P[upper_tri] <- t(df_corr$P)[upper_tri]
    df_corr$r[upper_tri] <- t(df_corr$r)[upper_tri]
    df_corr$n[upper_tri] <- t(df_corr$n)[upper_tri]
    return(df_corr)
  })
}



strip_hemi_suffix <- function (x){return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}
