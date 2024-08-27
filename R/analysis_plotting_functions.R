#' @importFrom ggplot2 ggplot aes theme element_text unit geom_tile geom_text scale_fill_gradient2 labs ggsave guide_legend xlab ylab xlim ylim
NULL

#' @importFrom tidyselect all_of
NULL

#' @importFrom dplyr n mutate summarize summarise across arrange group_by
NULL

#' @importFrom tidygraph activate convert
NULL

#' @importFrom tidyr pivot_longer pivot_wider
NULL

##_____________________ Analysis functions ___________________________

#' Get the percentage of colabelled cells over either cfos or eyfp channels.
#' @description This analysis will only include common regions that are included in both the colabelled and
#' cfos or eyfp channels. The colabelled percentage of individual animals will be calculated with the option to export the data.
#' @param e experiment object
#' @param by (str) Attribute names to group by, e.g. by = c("group", "sex"). These will generate analysis subgroups that are
#' averaged together to assess across all rois.
#' @param colabel_channel (str, default = "colabel") The channel used as the numerator in fraction counts. The string 'colabel' in the pipeline
#' refers to colocalized 'eyfp' and 'cfos' channels. For other colocalized channels, import the channel using [SMARTR::import_segmentation_custom()]
#' or your own customized import channel.
#' @param channel (str, default = "eyfp") The channel used as denominator in fraction counts.
#' @param save_table (bool, default = TRUE) Whether to save the output table as a csv in the experiment object output folder.
#' @param rois (str, default = NULL) Whether to generate colabelled percentages for only specific regions of interest, e.g. rois = c("HY", "DG").
#' Child regions of specified rois will also be searched for.
#' @param individual (bool, default = FALSE) Whether the data should include individual mouse colabelled percentages rather than the average.
#' If FALSE the colabel percentages are averaged across all analysis subgroups determined by the `by` parameter
#' @return e experiment object with colabelled percentage table stored in it.
#' @export
#' @examples e <- get_percent_colabel(e, c("group", "sex", channel = "eyfp"))

get_percent_colabel <-function(e, by, colabel_channel = "colabel",
                               channel = "eyfp", save_table = TRUE, rois = NULL, individual = TRUE){

  # correct mismatched attributes typed by users
  e_info <- attr(e, "info")

  # Get common regions that are present across both channels
  common_reg <- e$combined_normalized_counts[[colabel_channel]]$acronym %>% unique() %>%
    intersect( e$combined_normalized_counts[[channel]]$acronym )

  # Check if user only wants a specific roi in common regions
  if (!is.null(rois)){
    common_reg <- rois_intersect_region_list(common_reg, rois)
  }

  if (is.null(by)){
    # Join by all attributes in the dataset
    by <- names(e$combined_normalized_counts$colabel) %>% setdiff(c("area.mm2", "volume.mm3",
                                                                    "normalized.count.by.area",
                                                                    "normalized.count.by.volume"))
  } else {
    # join by user defined attributes
    by <- c("mouse_ID", by, "acronym", "name", "count")
  }

  # Get the colabel counts
  colabel_counts <-  e$combined_normalized_counts[[colabel_channel]] %>% dplyr::filter(acronym %in% common_reg) %>%
    dplyr::select(all_of(by))

  # Get the channel counts
  channel_counts <-  e$combined_normalized_counts[[channel]] %>% dplyr::filter(acronym %in% common_reg) %>%
    dplyr::select(all_of(by))

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
    col_to_summarize <- c("count.colabel", paste0("count.", channel), "colabel_percentage")

    master_counts <- master_counts %>%
      dplyr::group_by(across(all_of(attr_group_by))) %>%
      dplyr::summarise(n = dplyr::n(),
                       colabel_percentage.sd = sd(colabel_percentage),
                       colabel_percentage.sem = colabel_percentage.sd/sqrt(n),
                       across(all_of(col_to_summarize), mean)) %>%
      dplyr::rename_with(function(x){paste0(x,".mean")}, all_of(col_to_summarize)) %>%
      dplyr::relocate(colabel_percentage.sd, colabel_percentage.sem, n, .after = last_col())

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
#'
#' @param e experiment object
#' @param by (str) Attribute names to group by, e.g. c("sex", "group")
#' @param values (str) The respective values of the attributes entered for the `by` parameter to generate a specific analysis group,
#' e.g.values = c("female", "AD").
#' @param channels (str, channels =  c("cfos", "eyfp", "colabel") The channels to process.
#' @param p_adjust_method (bool or str, default = "none") This parameter is fed into the p.adjust function. Options: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",  "fdr", "none")
#'  Apply the named method to control for the inflated false discovery rate or family wise error rate (FWER). Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param region_order (list, default = NULL)  optional list with the first element named "acronym" supplying a vector as region acronyms and the second element named
#' "order"  supplying an vector of integers determining numerical order, e.g. 1, 1, 2, 2.
#' @param alpha (num, default = 0.05) The alpha level for significance applied AFTER p-adjustment.
#' @param anatomical.order (vec, c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The default super region acronym list that sorts all subregions in the dataset for grouped
#' positioning in the correlation matrix. Supercedes `region_order` parameter.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#'
#' @return e experiment object. The experiment object now has a named `correlation_list` object stored in it.
#' The name of the correlation object is the concatenation of the variable values separated by a "_".
#' This name allows for unambiguous identification of different analysis subgroups in the future.
#' @export
#' @examples e <- get_correlations(e, by = c("sex", "group"), values = c("female", "AD"),
#' channels = c("cfos", "eyfp", "colabel"),  p_adjust_method = "BH", alpha = 0.05)
#' @seealso [Hmisc::rcorr()]
get_correlations <- function(e, by, values,
                             channels = c("cfos", "eyfp", "colabel"),
                             p_adjust_method = "none",
                             alpha = 0.05,
                             ontology = "allen",
                             anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                             region_order = NULL){
  corr_list <- list()
  names(corr_list) <- names(channels)

  for (channel in channels){

    # filter by parameters
    df_channel <-  e$combined_normalized_counts[[channel]] %>%
      filter_df_by_char_params(by, values) %>% dplyr::distinct()

    # Pivot wider
    df_channel <- df_channel %>%  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
      tidyr::pivot_wider(names_from = acronym, values_from = normalized.count.by.volume, values_fill = NA)

    # Rearrange the correlations to be in "anatomical order"
    common.regions <-  df_channel %>% dplyr::ungroup() %>% dplyr::select(-any_of(c('mouse_ID', attr2match$mouse))) %>% names()

    if (is.null(region_order)){
      if (tolower(ontology) == "allen"){
        common.regions.ordered <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure) %>%
          purrr::map(intersect, y=common.regions) %>% unlist()
      } else {
        # anatomical.order <- c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")
        common.regions.ordered <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure.custom, ontology = ontology) %>%
          purrr::map(intersect, y=common.regions) %>% unlist()
      }
      } else if (is.list(region_order)){
      common.regions.ordered <- region_order$acronym[region_order$order]
      common.regions.ordered <- intersect(common.regions.ordered , common.regions)
    } else {
      stop("You did not supply a list object for the variable `region_order`."
      )
    }

    df_channel <- df_channel %>% dplyr::ungroup() %>% dplyr::select(all_of(c(common.regions.ordered))) %>%
      as.matrix()

    # Get minimum number of mice that we have data for across all regions
    # n <- colSums(!is.na(df_channel)) %>% min()

    # Try-catch statement
    df_corr <-  try_correlate(df_channel)

    if (!isFALSE(p_adjust_method)){
      lowertri <- df_corr$P %>%  lower.tri(diag = FALSE)
      lower_p <- df_corr$P[lowertri]
      uppertri <- df_corr$P %>%  upper.tri(diag = FALSE)
      upper_p <- df_corr$P[uppertri]

      lower_p <- lower_p %>% p.adjust(method = p_adjust_method)
      df_corr$P[lowertri] <- lower_p
      upper_p <- upper_p %>% p.adjust(method = p_adjust_method)
      df_corr$P[uppertri] <- upper_p
    }

    df_corr$sig <- df_corr$P <= alpha
    diag(df_corr$sig) <- FALSE
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

#' This function performs a permutation analysis to compare the region pairwise correlation coefficients between two different analysis groups.
#'
#' @description The data from two different analysis groups are compared by specifying the
#' `correlation_list_name_1` and `correlation_list_name_2` parameters. Note that both of these analysis groups must have the same
#' number of channels to compare. The functions `get_correlations()` needs to have been run for each of these analysis groups prior to
#' running this function. The test statistics used is the pearson values of those in correlation_list_name_2 subtracted from corresponding Pearson values in correlation_list_name_1.
#' @param e experiment object
#' @param correlation_list_name_1 (str) The name of the correlation list object used as the first group for comparison.
#' @param correlation_list_name_2 (str) The name of the correlation list object used as the second group for comparison.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param n_shuffle (int, default = 1000) The number of permutation shuffles.
#' @param seed (int, default = 5) Random seed for future replication.
#' @param p_adjust_method (bool or str, default = "BH") Benjamini-Hochberg method is recommended.
#'  Apply the named method to control for the inflated false discovery rate or FWER. Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param alpha (float, default = 0.05) The alpha cutoff for significance between region pairwise correlation differences
# @param ... additional parameters to [RcppAlgos::permuteGeneral()] aside from n, m, Parallel and repetition
#' @return e experiment object. The experiment object now has a list called `permutation_p_matrix` stored in it. Elements of this `permutation_p_matrix` list are
#' the outputs of different permutation comparison analyses. These elements are named by the groups that were compared.
#' @export
#' @seealso [SMARTR::get_correlations()]
#' @examples e <- correlation_diff_permutation(sundowning,
#'                                             correlation_list_name_1 = "female_AD",
#'                                             correlation_list_name_2 = "female_control",
#'                                             channels = c("cfos", "eyfp", "colabel"),
#'                                             n_shuffles = 1000,
#'                                             seed = 5,
#'                                             p_adjust_method = FALSE
#'                                             alpha = 0.001,
#'                                             )
#'
correlation_diff_permutation <- function(e,
                                         correlation_list_name_1,
                                         correlation_list_name_2,
                                         channels = c("cfos", "eyfp", "colabel"),
                                         n_shuffle = 1000,
                                         seed = 5,
                                         p_adjust_method = "BH",
                                         alpha = 0.05,
                                         ...){

  # Return the correlations list data showing the grouping and the values
  attr_group_1 <- attributes(e$correlation_list[[correlation_list_name_1]])
  attr_group_2 <- attributes(e$correlation_list[[correlation_list_name_2]])

  # # Get overlapping regions between the two correlational datasets
  p_matrix_list <- vector(mode = "list", length = length(channels))
  names(p_matrix_list) <- channels

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
      dplyr::ungroup() %>%
      dplyr::select(-any_of(attr2match$mouse))

    df_channel_group_2  <- df_channel_group_2  %>%  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
      tidyr::pivot_wider(names_from = acronym, values_from = normalized.count.by.volume) %>%
      dplyr::mutate(corr_group = correlation_list_name_2) %>% dplyr::relocate(corr_group, .before = 2) %>%
      dplyr::ungroup() %>%
      dplyr::select(-any_of(attr2match$mouse))

    # Get common regions between each group
    common_regions_btwn_groups <- intersect(names(df_channel_group_1), names(df_channel_group_2))

    # Sort names into anatomical order
    # common_regions_btwn_groups <- sort_anatomical_order(common_regions_btwn_groups)

    # Select the common regions in anatomical order across the two group dataframes
    df_channel_group_1 <- df_channel_group_1 %>% dplyr::select(mouse_ID, corr_group, all_of(common_regions_btwn_groups)) %>% dplyr::mutate(mouse_ID = as.character(mouse_ID))

    df_channel_group_2 <- df_channel_group_2 %>% dplyr::select(mouse_ID, corr_group, all_of(common_regions_btwn_groups)) %>% dplyr::mutate(mouse_ID = as.character(mouse_ID))

    # Bind group dfs together
    df_channel_groups <- dplyr::bind_rows(df_channel_group_1, df_channel_group_2)

    # df_channel_groups <- dplyr::bind_cols(df_channel_groups %>% dplyr::select(c(mouse_ID, corr_group)),
    #                                       df_channel_groups %>% dplyr::select(dplyr::where(is.numeric)))

    # Get group region pairwise correlational differences
    group_1_corr <- df_channel_group_1 %>% dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>% try_correlate()
    group_2_corr <- df_channel_group_2 %>%  dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>% try_correlate()
    test_statistic <- group_2_corr$r - group_1_corr$r

    # Get an array of distribution of correlation differences

    suppressWarnings(
    test_statistic_distributions <- permute_corr_diff_distrib(df_channel_groups,
                                                         correlation_list_name_1 = correlation_list_name_1,
                                                         correlation_list_name_2 = correlation_list_name_2,
                                                         n_shuffle = n_shuffle,
                                                         seed = seed, ...)
    )

    # message("dim before sort", dim(test_statistic_distributions))
    # For each pairwise distribution, sort the values
    # test_statistic_distributions <- apply(test_statistic_distributions, 1:2, sort)
    # message("dim after sort", dim(test_statistic_distributions))
    # do.call(sort, 1:2)
    # test_statistic_distributions <- test_statistic_distributions %>% aperm(c(1, 2, 3))
    # aperm(c(3, 2, 1))

    # calculate the p-value of the permutation
    p_matrix <- matrix(nrow = dim(test_statistic)[1],
                       ncol = dim(test_statistic)[1],
                       dimnames = dimnames(test_statistic))

    # Change to names just to be precise?
    l_reg <- dimnames(test_statistic)[1] %>% unlist()
    for (i in 1:length(l_reg)){
      for (j in 1:length(l_reg)){
        # print(l_reg[i])
        # print(l_reg[j])
        # test_statistic_distributions %>% dim() %>% print()
        null_distrib <- test_statistic_distributions[l_reg[i],l_reg[j],] %>% unlist()
        p_matrix[l_reg[i],l_reg[j]] <-  (sum(abs(null_distrib) >= abs(test_statistic[l_reg[i],l_reg[j]])) + 1) / (n_shuffle + 1)

        if(j>=i){
          p_matrix[i,j:length(l_reg)] <- NA
        }
      }
    }

    # adjust the p-value for false discovery rate or FWER
    if (!isFALSE(p_adjust_method)){
      # Calculate without removing NAs
      p_matrix <- p_matrix %>% p.adjust(method = p_adjust_method) %>%
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
  # Store the results in the experiment object
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
        drop_na()


      pval_df <-e$permutation_p_matrix[[pg]][[channel]]$p_val %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "p_val") %>%
        drop_na()


      test_statistic <- e$permutation_p_matrix[[pg]][[channel]]$test_statistic
      test_statistic[upper.tri(test_statistic, diag = TRUE)] <- NA

      test_statistic <- test_statistic %>%
        tibble::as_tibble(rownames="R1") %>%
        tidyr::pivot_longer(cols = all_of(regions), names_to = "R2", values_to = "test_statistic") %>%
        drop_na()

      permutation_results <- dplyr::left_join(test_statistic, pval_df, by = c("R1", "R2"), keep = FALSE) %>%
        dplyr::left_join(sig_df, by = c("R1", "R2"), keep = FALSE)


      if (ontology == "allen"){
        permutation_results <- permutation_results %>% dplyr::mutate(name1 = name.from.acronym(R1),
                                                                     name2 = name.from.acronym(R2),
                                                                     .before = test_statistic)
      } else (
        permutation_results <- permutation_results %>% dplyr::mutate(name1 = name.from.acronym.custom(R1, ontology = ontology),
                                                                     name2 = name.from.acronym.custom(R2, ontology = ontology),
                                                                     .before = test_statistic)
      )

      if (filter_significant) {
        permutation_results <- permutation_results %>% dplyr::filter(sig)
      }

      # Create table directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "tables")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }

      write.csv(permutation_results, file.path(output_dir, paste0("permutation_results", "_", pg, "_", channel, ".csv")), row.names = FALSE)
    }
  }
}






#' Create %>% graph objects for plotting different analysis subgroups.
#'
#' @param e experiment object
#' @param correlation_list_name (str) Name of the correlation list object used to generate the networks.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param alpha (float, default = 0.05) The significance threshold for including brain regions in the network. if NULL or NA,
#' this threshold is not applied.
#' @param pearson_thresh (float, default = 0.8) The pearson correlation coefficient threshold to apply for filtering out
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param anatomical.order (vec, c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The default super region acronym list that groups all subregions in the dataset.
#' @return e experiment object. This object now has a new added element called `networks.` This is a list storing a
#' graph object per channel for each network analysis run.
#' The name of each network (`network_name`) is the same as the `correlation_list_name`
#' used to generate the network. This `network_name` is fed as a parameter into the
#' [SMARTR::plot_network()] function.
#' @export
#' @examples e sundowning <- create_networks(sundowning, correlation_list_name = "female_control", alpha = 0.05)
#' @seealso [SMARTR::plot_networks()]

create_networks <- function(e,
                            correlation_list_name,
                            channels = c("cfos", "eyfp", "colabel"),
                            alpha = 0.05,
                            ontology = "allen",
                            anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB"),
                            pearson_thresh = 0.8){

  # List to store the networks
  networks <- vector(mode = "list", length = length(channels))
  names(networks) <- channels

  for (channel in channels){

    #__________________ Create Edge Tibble _____________
    corr_df <-  e$correlation_list[[correlation_list_name]][[channel]] %>% purrr::map(tibble::as_tibble)
    val_names <- names(corr_df)

    for (k in 1:length(corr_df)){
      corr_df[[k]] <- tibble::add_column(corr_df[[k]], row_acronym = names(corr_df[[k]]), .before = TRUE) %>%
        tidyr::pivot_longer(!row_acronym, names_to = "col_acronym", values_to = val_names[k])
    }

    # Combined the correlation plots into one dataframe
    df <- corr_df$r %>% dplyr::left_join(corr_df$n, by = c("row_acronym", "col_acronym")) %>%
      dplyr::left_join(corr_df$P, by = c("row_acronym", "col_acronym")) %>%
      dplyr::left_join(corr_df$sig, by = c("row_acronym", "col_acronym")) %>%
      tidyr::drop_na()

    # Edges dataframe
    edges <- df %>% dplyr::mutate(sign = ifelse(r >= 0, "pos", "neg"),
                                  weight = abs(r)) %>%
       dplyr::rename(from = row_acronym,
             to = col_acronym,
             p.value = P)

    # ___________Create Node tibble ________________


    # Get unique common regions
    acronyms <-  e$correlation_list[[correlation_list_name]][[channel]]$r %>% rownames()

    # get the parent super region
    super.region <- acronyms

    if (tolower(ontology) == "allen"){
      for (sup.region in anatomical.order){
        super.region[super.region %in% SMARTR::get.sub.structure(sup.region)] <- sup.region
      }
    } else {
      for (sup.region in anatomical.order){
        super.region[super.region %in% SMARTR::get.sub.structure.custom(sup.region, ontology = ontology)] <- sup.region
      }
    }


    # get the color
    # color <- wholebrain::color.from.acronym(super.region)
    # nodes <- tibble::tibble(name = acronyms, super.region = super.region, color = color)
    nodes <- tibble::tibble(name = acronyms, super.region = super.region)

    # _____________ Create the network ________________
    network <- tidygraph::tbl_graph(nodes = nodes,
                                    edges = edges,
                                    directed = FALSE)

    # Remove self-loops and parallel edges
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

    # filter by alpha
    if(!is.na(alpha) && !is.null(alpha)){
      network <- network %>% activate(edges) %>% dplyr::filter(p.value < alpha)
    }

    # Filter by pearson correlation coefficient
    if(!is.na(pearson_thresh) && !is.null(pearson_thresh)){
      network <- network %>% activate(edges) %>% dplyr::filter(weight > pearson_thresh)
    }

    network <- network %>% activate(nodes) %>%
      dplyr::mutate(super.region = factor(super.region,levels = unique(super.region)),
                    degree = tidygraph::centrality_degree(),
                    triangles = igraph::count_triangles(.),
                    clust.coef = ifelse(degree < 2, 0 , 2*triangles/(degree*(degree - 1)))) %>%
      dplyr::filter(degree > 0)

    # Add distance, efficiency, and btw metrics
    d <- igraph::distances(network, weights = NA)
    d[which(!is.finite(d))] <- NA
    diag(d) <- NA
    network <- network %>%
      dplyr::mutate(avg.dist = rowSums(d, na.rm = TRUE)/(n()-1),
             efficiency = rowSums(d^(-1),na.rm = TRUE)/(n()-1),
             btw = tidygraph::centrality_betweenness(weights = rep(1,nrow(tidygraph::.E())),
                                          directed = FALSE),
             group = as.factor(tidygraph::group_walktrap()))
    networks[[channel]] <- network
  }

  e$networks[[correlation_list_name]] <- networks
  return(e)
}

#' Summarize multiple networks. Create summary dataframes of for multiple networks and
#' calculate network statistics for each network.
#' @param e experiment object
#' @param network_names (str) The names of the networks to generate summary tables for, e.g. network_names = c("female_AD", "female_control")
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param save_stats (bool, default = TRUE) Save the summary stats as a csv file in the output folder
#' @param save_degree_distribution (bool, default = TRUE) Save the network degree distributions (frequencies of each degree) across each comparison group as a csv file.
#' @param save_betweenness_distribution (bool, default = TRUE) Save the betweenness distribution and summary as a csv.
#' @param save_efficiency_distribution (bool, default = TRUE) Save the efficiency distribution and summary as a csv.
#' @return e experiment object
#' @export
#' @examples e <- get_network_statistics(e,  network_names = c("female_AD", "female_control"),
#' channels = c("cfos", "eyfp", "colabel"), save_stats = TRUE, save_degree_distribution = TRUE)
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

  # Create summary tibble
  for (channel in channels){
    nodes_df_list <- list()
    edges_df_list <- list()
    for (network in network_names){

      nodes_df_list[[network]] <-  e$networks[[network]][[channel]] %>% tidygraph::activate(nodes) %>%
        tibble::as_tibble() %>% dplyr::mutate(group = network)
      edges_df_list[[network]] <-  e$networks[[network]][[channel]] %>% tidygraph::activate(edges) %>%
        tibble::as_tibble() %>% dplyr::mutate(group = network)
    }
    nodes_df <- dplyr::bind_rows(nodes_df_list) %>% dplyr::mutate(group=factor(group, levels = network_names))
    edges_df <- dplyr::bind_rows(edges_df_list) %>% dplyr::mutate(group=factor(group, levels = network_names))

    network_stats_df <- nodes_df %>% dplyr::group_by(group) %>%
      dplyr::summarise_if(is.numeric,list(~mean(.),~sd(.), ~sem(.))) %>%
      dplyr::left_join(nodes_df %>% dplyr::group_by(group) %>% dplyr::summarise(n.nodes = n())) %>%
      dplyr::left_join(edges_df %>% dplyr::group_by(group) %>% dplyr::summarise(n.edges = n())) %>%
      dplyr::mutate(edge.density = n.edges/(2*n.nodes*(n.nodes-1))) %>%
      dplyr::rename_all(stringr::str_replace,pattern = "_",replacement=".")

    #make table of degree frequency (useful for plotting degree histogram outside of R)
    degree_distribution_df <- nodes_df %>% dplyr::group_by(group, degree) %>% dplyr::count(degree)


    # Make a dataframe of degree (useful)
    degree_distribution <- nodes_df %>% dplyr::select(name, group, degree) %>%
      dplyr::group_by(group) %>% dplyr::arrange(group, dplyr::desc(degree)) %>%
      dplyr::mutate(name = factor(name, levels=name))

    #make table of betweenness frequency (useful for plotting degree histogram outside of R)
    betweenness_distribution_df <- nodes_df %>% dplyr::group_by(group, btw) %>% dplyr::count(btw)


    # Make a dataframe of betweenness (useful)
    betweenness_distribution <- nodes_df %>% dplyr::select(name, group, btw) %>%
      dplyr::group_by(group) %>% dplyr::arrange(group, dplyr::desc(btw)) %>%
      dplyr::mutate(name = factor(name, levels=name))

    #make table of efficiency frequency (useful for plotting degree histogram outside of R)
    efficiency_distribution_df <- nodes_df %>% dplyr::group_by(group, efficiency) %>% dplyr::count(efficiency)


    # Make a dataframe of efficiency (useful)
    efficiency_distribution <- nodes_df %>% dplyr::select(name, group, efficiency) %>%
      dplyr::group_by(group) %>% dplyr::arrange(group, dplyr::desc(efficiency)) %>%
      dplyr::mutate(name = factor(name, levels=name))

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
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param alpha (float, default = 0.05) The significance threshold for including brain regions in the network. if NULL or NA,
#' this threshold is not applied.
#' @param pearson_thresh (float, default = 0.8) The pearson correlation coefficient threshold to apply for filtering out
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param correlation_list_names (str vec) character vector of the two correlation lists used to include in a joined network
#' @param export_overlapping_edges (bool, default  = TRUE) Whether to export the overlapping edges between the two networks as a csv into the `table` directory.
#' @param anatomical.order (vec, c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")) The default super region acronym list that groups all subregions in the dataset.
#'
#' @return e experiment object. This object now has a new added element called `networks.` This is a list storing a
#' graph object per channel for each network analysis run. The name of each network (`network_name`) is the same as the `correlation_list_name`
#' used to generate the network. This `network_name` is fed as a parameter into the
#' [SMARTR::plot_network()] function.
#' @export
#' @examples e sundowning <- create_networks(sundowning, correlation_list_name = "female_control", alpha = 0.05)
#' @seealso [SMARTR::plot_networks()]

create_joined_networks <- function(e,
                                   correlation_list_names = c("male_agg", "female_non"),
                                   channels = "cfos",
                                   ontology = "unified",
                                   alpha = 0.001,
                                   pearson_thresh = 0.9,
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
          tidyr::pivot_longer(!row_acronym, names_to = "col_acronym", values_to = val_names[k])
      }

      # Combined the correlation plots into one dataframe
      df <- corr_df$r %>% dplyr::left_join(corr_df$n, by = c("row_acronym", "col_acronym")) %>%
        dplyr::left_join(corr_df$P, by = c("row_acronym", "col_acronym")) %>%
        dplyr::left_join(corr_df$sig, by = c("row_acronym", "col_acronym")) %>%
        tidyr::drop_na()

      # Edges dataframe
      edges_joined[[correlation_list_name]] <- df %>% dplyr::mutate(sign = ifelse(r >= 0, "pos", "neg"),
                                                                    weight = abs(r)) %>%
        dplyr::rename(from = row_acronym,
                      to = col_acronym,
                      p.value = P) %>%
        tibble::add_column(network = correlation_list_name)
      # ___________Create Node tibble ________________


      # Get unique common regions
      acronyms <-  e$correlation_list[[correlation_list_name]][[channel]]$r %>% rownames()

      # get the parent super region
      super.region <- acronyms

      if (tolower(ontology) == "allen"){
        for (sup.region in anatomical.order){
          super.region[super.region %in% SMARTR::get.sub.structure(sup.region)] <- sup.region
        }
      } else {
        for (sup.region in anatomical.order){
          super.region[super.region %in% SMARTR::get.sub.structure.custom(sup.region, ontology = ontology)] <- sup.region
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



    # Remove self-loops and parallel edges
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

    network <- tidygraph::graph_join(network1, network2)

    # filter by alpha
    if(!is.na(alpha) && !is.null(alpha)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::filter(p.value < alpha)
    }

    # Filter by pearson correlation coefficient
    if(!is.na(pearson_thresh) && !is.null(pearson_thresh)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::filter(weight > pearson_thresh)
    }

    network <- network %>% tidygraph::activate(nodes) %>%
      dplyr::mutate(super.region = factor(super.region,levels = unique(super.region)),
                    degree = tidygraph::centrality_degree(mode="all")) %>%
      dplyr::filter(degree>0)

    networks[[channel]] <- network

    # Whether to export the list of overlapping edges
    if (export_overlapping_edges) {
      joined_network_name <- paste(correlation_list_names, collapse = "_")
      e$networks[[joined_network_name]][[channel]] %>% activate(nodes) %>% as_tibble() -> nodes_df
      e$networks[[joined_network_name]][[channel]] %>% activate(edges) %>% as_tibble() -> edges_df
      edges_df <- edges_df %>% mutate(from = nodes_df$name[edges_df$from],
                                      to =   nodes_df$name[edges_df$to])
      edges_df_p1 <- edges_df %>% dplyr::filter(network == correlation_list_names[1])
      edges_df_p2 <- edges_df %>% dplyr::filter(network == correlation_list_names[2])
      mutual_edges <- edges_df_p1 %>% dplyr::inner_join(edges_df_p2, by = c("from", "to"), suffix = c(paste0(".", correlation_list_names[1]), paste0(".", correlation_list_names[2])))

      # Create table directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "tables")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      file_path <- file.path(output_dir, paste0("joined_networks_mutual_edges_",joined_network_name, "_", channel,".csv"))
      write.csv(mutual_edges, file = file_path)
      message(paste0("Saved mutual edge list at: ", file_path))
    }
  }

  e$networks[[paste(correlation_list_names, collapse = "_")]] <- networks
  return(e)
}



##_____________________ Plotting functions ___________________________

#' This function allows for plotting of colabelled cells over either the "cfos" or "eyfp" channels. And allows for specification
#' of specific brain regions to plot. Two different mouse attributes can be used as categorical variables to map to either the color or
#' pattern aesthetics of the bar plot, e.g. sex and experimental group.
#' The color aesthetic takes precedence over the pattern aesthetic so if you only want to use one mouse attribute, for plotting
#' set it to the `color_mapping` parameter and set the `pattern_mapping` parameter to NULL.
#'
#' @param e experiment object
#' @param channel (str, default = "eyfp") The channel used as denominator in fraction counts.
#' @param rois
#' @param color_mapping
#' @param colors
#' @param pattern_mapping
#' @param patterns
#' @param error_bar (str, c("sd", "sem)) options for which type of error bar to display, standard deviation or standard error of the mean.
#' @return p Plot handle to the figure
#' @export
#' @examples plot_percentage_colabel
plot_percent_colabel <- function(e,
                                 colabel_channel = "colabel",
                                 channel = "eyfp",
                                 rois = c("AAA", "dDG", "HY"),
                                 color_mapping = "sex",
                                 colors = c("#952899", "#358a9c"),
                                 pattern_mapping = NULL,
                                 patterns = c("gray100", 'hs_fdiagonal', "hs_horizontal", "gray90", "hs_vertical"),
                                 error_bar = "sem",
                                 ylim = c(0, 100),
                                 plot_individual = TRUE,
                                 height = 8,
                                 width = 8,
                                 print_plot = TRUE,
                                 save_plot = TRUE,
                                 image_ext = ".png"){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # Re-running the get_percent_colabel function for the color &/or pattern mapping
  if(is.null(pattern_mapping)){
    by <- match_m_attr(color_mapping)
  } else{
    by <- match_m_attr(c(color_mapping, pattern_mapping))
  }

  # Get the average and individual percentages
  e <- get_percent_colabel(e, by = by, colabel_channel = colabel_channel, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
  e <- get_percent_colabel(e, by = by, colabel_channel = colabel_channel, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)

  # Get rois of all the child regions
  child_r <- SMARTR::get.acronym.child(rois)
  while (length(child_r) > 0){
    rois <- c(rois, child_r)
    child_r <- SMARTR::get.acronym.child(child_r) %>% na.omit()
  }

  # # Check that a colabel percentage dataframe exists for this channel
  # if(is.null(e$colabel_percent[[channel]])){
  #   message("The percent colabelled data does not exist in your experiment object for this channel. Autorunning the function now...")
  #   e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
  #   e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)
  # }

  # Make sure the by (if not NULL) is equal to the by parameter fed to get_percent_colabel() used to generate the dataframes
  # attrib <- attributes(e$colabel_percent[[channel]]$average)
  # stored_by <- names(attrib$groups) %>% setdiff(c("acronym", ".rows"))
  #
  # if (!is.null(by) && all(by %in% stored_by)){
  #   loop <- TRUE
  #   while(loop){
  #     message("The 'by' parameter used here doesn't include attributes input to the 'by' parameter that was used in get_percent_colabel().\n",
  #             "The groups that you intend to compare may be different.")
  #
  #     inp <- readline("Would you like to rerun the get_percent_colabel() function using the new 'by' parameters for grouping?: Y/N?" )
  #     if (inp=="Y" || inp=="y") {
  #       e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
  #       e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)
  #       loop <- FALSE
  #     } else if ( inp=="N" || inp == "n") {
  #       stop("Cannot continue with plotting using the current function parameters.")
  #     }
  #   }
  # } else{
  #   by <- stored_by
  # }

  # Error bar check
  if (error_bar == "sem"){
    error_var <- sym("colabel_percentage.sem")
  } else if (error_bar == "sd") {
    error_var <- sym("colabel_percentage.sd")
  } else {
    stop("You did not supply a valid option for the error_bar. Valid options are 'sem' and 'sd'.")
  }

  # Set up custom bar plot theme
  bar_theme <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     strip.background = element_rect(fill = "white"),
                     strip.text =  element_text(size = 20, color = "black"),
                     plot.title = element_text(hjust = 0.5, size = 36,  color = "black"),
                     axis.text.y = element_text(size = 20,  color = "black"),
                     axis.title = element_text(size = 20, color = "black"),
                     legend.text = element_text(size = 20,  color = "black"),
                     legend.title = element_text(size = 20, color = "black"),
                     # axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

  # Create  variable names for mapping
  color_var <- sym(color_mapping)

  # Use pattern as a variable if not NULL
  if (!is.null(pattern_mapping)){

    # Check if ggpattern is installed
    if (!requireNamespace("ggpattern", quietly = TRUE)) {
      message("Package \"ggpattern\" (>= 0.2.0) is needed for the pattern mapping to work. Installing it now.",
           call. = FALSE)
      install.packages('ggpattern')
    }

    pattern_var <- sym(pattern_mapping)

    # Create the plot
    p <- ggplot(e$colabel_percent[[channel]]$average %>% dplyr::filter(acronym %in% rois),
                aes(x = interaction(!!pattern_var, !!color_var),
                    y = colabel_percentage.mean,
                    fill = !!color_var)) +
      ggpattern::geom_bar_pattern(aes(pattern_type = !!pattern_var),
                                  stat = "identity",
                                  color = "black",
                                  pattern_color = "black",
                                  pattern_fill = "black",
                                  pattern = "magick",
                                  position = position_dodge(width = .5),
                                  show.legend = TRUE,
                                  width = 0.4) +
      ggplot2::scale_fill_manual(values = colors) +
      ggpattern::scale_pattern_type_manual(values = patterns) +
      geom_errorbar(aes(ymin = colabel_percentage.mean - !!error_var,
                        ymax = colabel_percentage.mean + !!error_var),
                    width=.2,
                    position = position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', size=0.5, linetype='solid')) +
      facet_wrap(~acronym,
                 strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = paste0("co-labeled / ", channel, "+ cells (%)")) +
      bar_theme

    if (plot_individual){
      df <- e$colabel_percent[[channel]]$individual %>% dplyr::filter(acronym %in% rois)
      p <- p + geom_jitter(data = df,
                           aes(x = interaction(!!pattern_var, !!color_var),
                               y = colabel_percentage),
                           size = 2,
                           width = 0.1)
    }
  } else{

    # Create the plot
    p <- ggplot(e$colabel_percent[[channel]]$average %>% dplyr::filter(acronym %in% rois),
                aes(x = !!color_var,
                    y = colabel_percentage.mean,
                    fill = !!color_var)) +
      ggplot2::geom_bar(stat='identity',
                        position = position_dodge(width = .5),
                        color = "black",
                        show.legend = TRUE,
                        width = 0.4) +
      scale_fill_manual(values = colors) +
      geom_errorbar(aes(ymin = colabel_percentage.mean - !!error_var,
                        ymax = colabel_percentage.mean + !!error_var),
                    width=.2,
                    position = position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', size=0.5, linetype='solid')) +
      facet_wrap(~acronym,
                 strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = paste0("co-labeled / ", channel, "+ cells (%)")) +
      bar_theme

    if (plot_individual){
      df <- e$colabel_percent[[channel]]$individual %>% dplyr::filter(acronym %in% rois)
      p <- p + geom_jitter(data = df,
                           aes(x = !!color_var,
                               y = colabel_percentage),
                           size = 2,
                           width = 0.1)
    }
  }

  if (print_plot){
    quartz()
    print(p)
  }

  if(save_plot){
    # Plot the plot
    quartz(width = width, height = height)
    print(p)
    # Create figure directory if it doesn't already exists
    output_dir <-  file.path(attr(e, "info")$output_path, "figures")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    image_file <- file.path(output_dir,
                            paste0("Colabeled_over_", channel, "_plot_by_", paste0(by, collapse = "_"), image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
  }
  return(p)
}


#' This function allows for plotting of normalized cell counts by area across specific regions to plot. Two different mouse attributes can be used as categorical variables to map to either the color or
#' pattern aesthetics of the bar plot, e.g. sex and experimental group.
#' The color aesthetic takes precedence over the pattern aesthetic so if you only want to use one mouse attribute, for plotting
#' set it to the `color_mapping` parameter and set the `pattern_mapping` parameter to NULL.
#'
#' @param e experiment object
#' @param channel (str, default = "eyfp") The channel used as denominator in fraction counts.
#' @param rois (vec,  default = c("AAA", "dDG", "HY")) Allen acronyms of the ROIS that the user would like to plot
#' @param color_mapping (str, default = "group") The variable name that maps subgroups  you would like to graphically distinguish through colors.
#' @param colors (vec, default = c("#952899", "#358a9c")) A vector of hexadecimal color codes for each subgroup distinguished by the color mapping variable.
#' @param pattern_mapping (str, default = NULL) variable name that maps subgroups  you would like to graphically distinguish through bar patterns.Set to NULL if not in use.
#' @param patterns  (default = c("gray100", 'hs_fdiagonal', "hs_horizontal", "gray90", "hs_vertical") Available patterns in ggpattern package to map to subgroups distinguished by the pattern mapping variable.
#' @param ylab (str, default =  bquote('Cell counts '('cells/mm'^3))) unit of measurement
#' @param ylim
#' @param plot_individua (bool) whether to plot individual points
#' @param height (numeric) height of the plot in inches to save as.
#' @param width (numeric) height of the plot in inches to save as.
#' @param print_plot
#' @param save_plot (book) Whether to save the plot in the experiment figures folder.
#' @param image_ext (str, default = ".png") Extension of the output file
#' @param error_bar (str, c("sd", "sem)) options for which type of error bar to display, standard deviation or standard error of the mean.
#'
#' @return p Plot handle to the figure
#' @export
#' @examples plot_percentage_colabel
plot_cell_counts <- function(e,
                             channel = "eyfp",
                             rois = c("AAA", "dDG", "HY"),
                             color_mapping = "group",
                             colors = c("#952899", "#358a9c"),
                             pattern_mapping = NULL,
                             patterns = c("gray100", 'hs_fdiagonal', "hs_horizontal", "gray90", "hs_vertical"),
                             ylab = bquote('Cell counts '('cells/mm'^3)),
                             error_bar = "sem",
                             ylim = c(0, 100),
                             plot_individual = TRUE,
                             height = 8,
                             width = 8,
                             print_plot = TRUE,
                             save_plot = TRUE,
                             image_ext = ".png"){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # Re-running the get_percent_colabel function for the color &/or pattern mapping
  if(is.null(pattern_mapping)){
    by <- match_m_attr(color_mapping)
  } else{
    by <- match_m_attr(c(color_mapping, pattern_mapping))
  }

  if(is.null(e$combined_normalized_counts)){
    warning("You have not run `combine_cell_counts()` yet. Running this automatically based on your color_mapping and pattern_mappings parameters.")
   e <- combine_cell_counts(e, by = by)
  }
  # Get the average and individual percentages
  # e <- get_percent_colabel(e, by = by, colabel_channel = colabel_channel, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
  # e <- get_percent_colabel(e, by = by, colabel_channel = colabel_channel, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)

  # Get rois of all the child regions
  child_r <- SMARTR::get.acronym.child(rois)
  while (length(child_r) > 0){
    rois <- c(rois, child_r)
    child_r <- SMARTR::get.acronym.child(child_r) %>% na.omit()
  }


  df <- e$combined_normalized_counts[[channel]] %>% dplyr::group_by(group, acronym, name) %>%
    dplyr::summarise(normalized.count.mean = mean(normalized.count.by.volume),
                     normalized.count.sem = sd(normalized.count.by.volume)/sqrt(n()),
                     normalized.count.sd = sd(normalized.count.by.volume),
                     n = n())

  df <- df %>% dplyr::filter(acronym %in% rois)
  df_indiv <- e$combined_normalized_counts[[channel]]  %>% dplyr::filter(acronym %in% rois)


  # Error bar check
  if (error_bar == "sem"){
    error_var <- sym("normalized.count.sem")
  } else if (error_bar == "sd") {
    error_var <- sym("normalized.count.sem")
  } else {
    stop("You did not supply a valid option for the error_bar. Valid options are 'sem' and 'sd'.")
  }

  # Set up custom bar plot theme
  bar_theme <- theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     strip.background = element_rect(fill = "white"),
                     strip.text =  element_text(size = 20, color = "black"),
                     plot.title = element_text(hjust = 0.5, size = 36,  color = "black"),
                     axis.text.y = element_text(size = 20,  color = "black"),
                     axis.title = element_text(size = 20, color = "black"),
                     legend.text = element_text(size = 20,  color = "black"),
                     legend.title = element_text(size = 20, color = "black"),
                     # axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

  # Create  variable names for mapping
  color_var <- sym(color_mapping)

  # Use pattern as a variable if not NULL
  if (!is.null(pattern_mapping)){

    # Check if ggpattern is installed
    if (!requireNamespace("ggpattern", quietly = TRUE)) {
      message("Package \"ggpattern\" (>= 0.2.0) is needed for the pattern mapping to work. Installing it now.",
              call. = FALSE)
      install.packages('ggpattern')
    }

    pattern_var <- sym(pattern_mapping)

    # Create the plot
    p <- ggplot(df,
                aes(x = interaction(!!pattern_var, !!color_var),
                    y = normalized.count.mean,
                    fill = !!color_var)) +
      ggpattern::geom_bar_pattern(aes(pattern_type = !!pattern_var),
                                  stat = "identity",
                                  color = "black",
                                  pattern_color = "black",
                                  pattern_fill = "black",
                                  pattern = "magick",
                                  position = position_dodge(width = .5),
                                  show.legend = TRUE,
                                  width = 0.4) +
      ggplot2::scale_fill_manual(values = colors) +
      ggpattern::scale_pattern_type_manual(values = patterns) +
      geom_errorbar(aes(ymin = normalized.count.mean - !!error_var,
                        ymax = normalized.count.mean + !!error_var),
                    width=.2,
                    position = position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', size=0.5, linetype='solid')) +
      facet_wrap(~acronym,
                 strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = ylab) +
      bar_theme

    if (plot_individual){
      p <- p + geom_jitter(data = df_indiv,
                           aes(x = interaction(!!pattern_var, !!color_var),
                               y = count),
                           size = 2,
                           width = 0.1)
    }
  } else{

    # Create the plot
    p <- ggplot(df,
                aes(x = !!color_var,
                    y = normalized.count.mean,
                    fill = !!color_var)) +
      ggplot2::geom_bar(stat='identity',
                        position = position_dodge(width = .5),
                        color = "black",
                        show.legend = TRUE,
                        width = 0.4) +
      scale_fill_manual(values = colors) +
      geom_errorbar(aes(ymin = normalized.count.mean - !!error_var,
                        ymax = normalized.count.mean + !!error_var),
                    width=.2,
                    position = position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', size=0.5, linetype='solid')) +
      facet_wrap(~acronym,
                 strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = ylab) +
      bar_theme

    if (plot_individual){
      p <- p + geom_jitter(data = df_indiv,
                           aes(x = !!color_var,
                               y = normalized.count.by.volume),
                           size = 2,
                           width = 0.1)
    }
  }

  if (print_plot){
    quartz()
    print(p)
  }

  if(save_plot){
    # Plot the plot
    quartz(width = width, height = height)
    print(p)
    # Create figure directory if it doesn't already exists
    output_dir <-  file.path(attr(e, "info")$output_path, "figures")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    image_file <- file.path(output_dir,
                            paste0("region_count_", channel, "_plot_by_", paste0(by, collapse = "_"), image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
  }
  return(p)
}






#' Plot normalized cell counts
#' @description Plot the cell counts normalized by volume for a given channel
#'
#' @param e experiment object
#' @param channels (str, default = c("cfos", "eyfp", "colabel"))
#' @param by (str) Attribute names to group by, e.g. c("sex", "group")
#' @param values (list) A list with a length the number of groups desired for plotting. Each element of the list is a vector in the order of the
#' the respective values for the attributes entered for the `by` parameter to generate a specific analysis group. Each vector should be unique to generate a
#' uniquely colored bar.
#' e.g.values = c("female", "AD").
#' @param colors (str, default = c("white", "lightblue")) Hexadecimal codes corresponding to the groups (respectively) to plot. The length of this vector should be the length of the list.
#' @param title (str, default = NULL) An optional title for the plot
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param print_plot (bool, default = TRUE) Whether to display the plot (in addition to saving the plot)
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of
#'  the experiment object output folder.
#' @param limits (c(0,100000)) Range of the normalized cell counts.
#' @param flip_axis plot cell counts on x-axis rather than y-axis.
#' @param unit_label (str, default = bquote('Cell counts '('cells/mm'^3))) Default unit label for the graphs
#' @param region_label_angle (int, default = 50) Angle of the region label.
#' @param label_text_size  (int, default = 8) Size of the text element for the region labels.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param legend.position
#' @param facet_background_color (default = NULL) Set to a hexadecimal string, e.g."#FFFFFF", when you want to shade the background of the graph. Defaults to no background when NULL.
#' @param anatomical.order (default = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU","TH", "HY", "MB", "HB", "CB")) Default way to group subregions into super regions order
#' @param legend.justification
#' @param legend.position.inside
#' @param legend.direction

#' The exact acronyms used will be ontology dependent so make sure that these match the `ontology` parameter
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @example
#' p_list <- plot_normalized_counts(e, channels = "cfos", by = c("sex", "group"), values = list(c("female", "non"), c("female", "agg")), colors = c("white", "lightblue"))
#'
plot_normalized_counts <- function(e,
                                   channels = c("cfos", "eyfp", "colabel"),
                                   by = c("sex", "group"),
                                   values = list(c("female", "non"),
                                                 c("female", "agg"),
                                                 c("female", "control"),
                                                 c("male", "agg"),
                                                 c("male", "control")),
                                   # groups = c("Context", "Shock"),
                                   colors = c("white", "lightblue", "black", "red", "green"),
                                   ontology = "allen",
                                   # regions_to_remove = c("CTX", "grey", "MB", "TH", "HY"),
                                   title = NULL,
                                   unit_label = bquote('Cell counts '('cells/mm'^3)),
                                   anatomical.order = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU",
                                                         "TH", "HY", "MB", "HB", "CB"),
                                   height = 7,
                                   width = 20,
                                   region_label_angle = 50,
                                   label_text_size = 8,
                                   print_plot = TRUE,
                                   save_plot = TRUE,
                                   flip_axis = FALSE,
                                   legend.justification = c(0, 0),
                                   legend.position = "inside",
                                   legend.position.inside = c(0.05, 0.6),
                                   legend.direction = "vertical",
                                   limits = c(0, 100000),
                                   facet_background_color =  NULL,
                                   image_ext = ".pdf"){
  # check os and set graphics window
  if (get_os() != "osx") {
    quartz <- X11
  }
  # create plot list containing number of plots equal to the number of channels
  p_list <- vector(mode = 'list', length = length(channels))
  labels <- purrr::map(values, function(x){paste(x, collapse = "_")}) %>% unlist()

  names(p_list) <- channels
  for (k in 1:length(channels)) {

    channel <- channels[[k]]
    channel_counts <- NULL
    for (g in 1:length(values)){
      if (is.null(channel_counts)){
        # filter by parameters
        channel_counts <-  e$combined_normalized_counts[[channel]] %>% dplyr::ungroup() %>%
          filter_df_by_char_params(by, values[[g]]) %>% dplyr::distinct() %>%
          dplyr::select(dplyr::all_of(c(by, "mouse_ID", "name", "acronym", "normalized.count.by.volume"))) %>%
          dplyr::group_by(across(all_of(c(by, "acronym", "name" )))) %>%
          dplyr::summarise(mean_normalized_counts = mean(normalized.count.by.volume),
                           sem = sd(normalized.count.by.volume, na.rm=TRUE)/sqrt(n()))
      } else {
        to_bind <-  e$combined_normalized_counts[[channel]] %>% dplyr::ungroup() %>%
          filter_df_by_char_params(by, values[[g]]) %>% dplyr::distinct() %>%
          dplyr::select(dplyr::all_of(c(by, "mouse_ID", "name", "acronym", "normalized.count.by.volume"))) %>%
          dplyr::group_by(across(all_of(c(by, "acronym", "name" )))) %>%
          dplyr::summarise(mean_normalized_counts = mean(normalized.count.by.volume),
                           sem = sd(normalized.count.by.volume, na.rm=TRUE)/sqrt(n()))
        channel_counts <- dplyr::bind_rows(channel_counts, to_bind)
      }
    }

    # # Filter regions meeting minimum n of both groups
    # region_counts <- channel_counts %>% group_by(acronym, name) %>% summarize(n = sum(!is.na(sem))) %>%
    #   filter(n == length(values))
    # channel_counts <- channel_counts %>% filter(acronym %in% region_counts$acronym)

    # if (!(isFALSE(regions_to_remove) || is.null(regions_to_remove) || is.na(regions_to_remove))){
    #   # remove parent regions (if any are specified) containing residual counts which are better represented in subregions
    #   channel_counts <- channel_counts[-c(which(channel_counts$acronym %in% regions_to_remove)),]
    #   # gather parent regions and assign these regions to subregions containing counts
    # }

    if (tolower(ontology) == "allen") {
      regions.ordered <- anatomical.order %>%
        purrr::map(SMARTR::get.sub.structure) %>%
        unlist()
    } else {
      regions.ordered <- anatomical.order %>%
        purrr::map(SMARTR::get.sub.structure.custom, ontology = ontology) %>%
        unlist()
    }

    common.regions.ordered <- channel_counts$name[unique(match(regions.ordered, channel_counts$acronym))]
    channel_counts$name <- factor(channel_counts$name, levels = common.regions.ordered)

    if (tolower(ontology) == "allen") {
    channel_counts <- channel_counts %>%
      mutate(parent = get.sup.structure(acronym, matching.string = anatomical.order)) %>%
      na.omit()
    channel_counts$parent <- factor(channel_counts$parent, levels = anatomical.order)
    } else {
      channel_counts <- channel_counts %>%
        mutate(parent = get.sup.structure.custom(acronym, ontology = ontology, matching.string = anatomical.order)) %>%
        na.omit()
      channel_counts$parent <- factor(channel_counts$parent, levels = anatomical.order)
    }


    # Annoying work around for the dynamic references of column names with a string vector
    dynamic_ref <-  rep(".data[[by[[1]]]]", times = length(by))
    for (b in 1:length(by)){
      dynamic_ref[[b]] <- dynamic_ref[[b]] %>% gsub("1", as.character(b), .)
    }
    dynamic_ref <- paste(dynamic_ref, collapse = ", ")
    text <- paste( "channel_counts <- channel_counts %>% mutate(unique_groups=paste(", dynamic_ref, ", sep='_' ))", collapse = "")
    eval((parse(text=text)))

    # plot normalized cell counts, grouped by specified groups and ordered by parent region
    channel_counts$unique_groups <- factor(channel_counts$unique_groups, levels = rev(labels))

    if (flip_axis) {
      p <- channel_counts %>%
        ggplot(aes(x = mean_normalized_counts, y = name,
                   fill = unique_groups), color = "black") +
        geom_col(position = position_dodge(0.8), width = 0.8, color = "black") +
        geom_errorbar(aes(xmin = mean_normalized_counts - sem,
                          xmax = mean_normalized_counts + sem, y = name),
                      position = position_dodge(0.8),
                      width = 0.5,
                      color = "black") +
        labs(title = title,
             x = unit_label,
             y = "",
             fill = "Group") +
        scale_x_continuous(expand = c(0,0), limits = limits) +
        scale_fill_manual(values=c(colors)) +
        facet_grid(parent~., scales = "free", space = "free_y", switch = "y") +
        theme_bw() +
        theme(plot.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank()) +
        theme(axis.line = element_line(color = 'black')) +
        theme(legend.justification = legend.justification,
              legend.position = legend.position,
              legend.position.inside = legend.position.inside,
              legend.direction = legend.direction) +
        theme(axis.text.y = element_text(angle = region_label_angle,
                                         hjust = 1,
                                         size = label_text_size,
                                         color = "black"),
              axis.text.x = element_text(color = "black")) +
        theme(strip.text.y = element_text(angle = 0, margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")),
              strip.placement = "outside",
              strip.background = element_rect(color = "black",
                                              fill = "lightblue"),
              strip.switch.pad.grid = unit(0.1, "in"))
    } else if (isFALSE(flip_axis)) {
      p <- channel_counts %>%
        ggplot(aes(y = mean_normalized_counts, x = name,
                   fill = unique_groups), color = "black") +
        geom_col(position = position_dodge(0.8), width = 0.8, color = "black") +
        geom_errorbar(aes(ymin = mean_normalized_counts - sem,
                          ymax = mean_normalized_counts + sem, x = name),
                      position = position_dodge(0.8),
                      width = 0.5,
                      color = "black") +
        labs(title = title,
             y = unit_label,
             x = "",
             fill = "Group") +
        scale_y_continuous(expand = c(0,0), limits = limits) +
        scale_fill_manual(values=colors, labels = labels) +

        facet_grid(~parent, scales = "free_x", space = "free_x", switch = "x") +
        theme_bw() +
        theme(
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank()) +
        theme(axis.line = element_line(color = 'black')) +
        theme(legend.justification = legend.justification,
              legend.position = legend.position,
              legend.position.inside = legend.position.inside,
              legend.direction = legend.direction) +
        theme(axis.text.x = element_text(angle = region_label_angle, hjust = 1, size = label_text_size, color = "black"),
              axis.text.y = element_text(color = "black")) +
        theme(strip.text.x = element_text(angle = 0, margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")),
              strip.placement = "outside",
              strip.background = element_rect(color = "black",
                                              fill = "lightblue")) +
        theme(plot.margin = margin(1,1.5,0,1.5, "cm"))
    }

    if(!is.null(facet_background_color)){
      p <- p +
      theme(panel.background = element_rect(fill = facet_background_color, color = "black"))
    }

    if (print_plot) {
      quartz()
      print(p)
    }

    if (save_plot) {
      # create figures directory if not already created
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }

      # save in figures
      ggsave(p, filename = paste0(channels[k], "_normalized_counts", image_ext),
             path = file.path(attr(e, 'info')$output_path, "figures"), width = width, height = height, limitsize = FALSE)


    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

#' Plot correlation heatmaps
#'
#' @param e experiment object. Must contain a named correlation_list object generated by [SMARTR::get_correlations()]
#' @param correlation_list_name (str) The name of the correlation object generated by [SMARTR::get_correlations()]
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Must exist in the channels attribute of the correlation_list.
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal code for the colors corresponding channels parameter. Color values can also be
#' input compatible with ggplot2 plotting functions.
#' @param print_plot (bool, default = TRUE) Print the plot as graphics windows.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#' the experiment object output folder.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param plot_title (str, default = NULL) If NULL, the `correlation_list_name` will used as the title with underscores removed.
#' @param height (int) Height of the plot in inches.
#' @param sig_color (str, default = "yellow") Color of the significance symbol in R
#' @param width (int) Width of the plot in inches.
#' @param theme.hm Option to use custom ggplot2 theme if the user wants. See default values as example.
#' @param sig_size (default = 7) Point size for significance symbol.
#' @param sig_nudge_y (default = -0.7) Relative amount to nudge the significance symbols in the y direction to center over each square.
#' @param anatomical.order (default = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU","TH", "HY", "MB", "HB", "CB")) Default way to group subregions into super regions order
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples plot_correlation_heatmaps(e, correlation_list_name = "female_AD")
#' @seealso [SMARTR::get_correlations()]

plot_correlation_heatmaps <- function(e, correlation_list_name,
                                      channels = c('cfos', 'eyfp', 'colabel'),
                                      colors = c("#be0000", "#00782e", "#f09b08"),
                                      sig_color = "yellow",
                                      sig_nudge_y = -0.7,
                                      sig_size = 7,
                                      ontology = "allen",
                                      anatomical.order = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU",
                                                           "TH", "HY", "MB", "HB", "CB"),
                                      print_plot = TRUE, save_plot = TRUE, image_ext = ".png",
                                      plot_title = NULL, height = 10, width = 10,
                                      theme.hm = ggplot2::theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 8),
                                                                axis.text.y = element_text(vjust = 0.5, size = 8),
                                                                plot.title = element_text(hjust = 0.5, size = 36),
                                                                axis.title = element_text(size = 18),
                                                                legend.text = element_text(size = 22),
                                                                legend.key.height = unit(100, "points"),
                                                                legend.title = element_text(size = 22),
                                                                panel.spacing = unit(0.2, "lines"),
                                                                strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
                                                                strip.text.y = element_text(angle = 270,  hjust = 0.5, vjust = 0.5, size = 10),
                                                                strip.placement = "outside",
                                                                strip.background = element_rect(color = "black", fill = "lightblue"))
                                      ){

  # Detect the OS and set quartz (as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # Get the attributes for plotting
  if (is.null(e$correlation_list[[correlation_list_name]])){
    stop(paste0("Your experiment object doesn't contain the specified correlation_list object to plot! ",
                "\nRun the function, get_correlations() first."))
  }

  cl_attr <- attributes(e$correlation_list[[correlation_list_name]])
  # channels <- cl_attr$names
  names(colors) <- channels

  if (is.null(plot_title)){
    cl_attr$values <- stringr::str_to_title(cl_attr$values)
    plot_title <- paste(cl_attr$value, collapse = " ")
  }

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (channel in channels){
   # Turn into tibblee$
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

   # Get row and column parent super regions
  df <- df %>% dplyr::mutate(row_parent = get.super.regions(row_acronym, anatomical.order = anatomical.order, ontology = ontology), .before  = row_acronym) %>%
    dplyr::mutate(row_parent = factor(row_parent, levels = anatomical.order))
  df <- df %>% dplyr::mutate(col_parent = get.super.regions(col_acronym, anatomical.order = anatomical.order, ontology = ontology), .before  = row_acronym) %>%
    dplyr::mutate(col_parent = factor(col_parent, levels = rev(anatomical.order)))


  # Generate a correlation heatmap in anatomical order
  n_facet <- df$row_parent %>% unique() %>% length()
  p <-  ggplot(df, aes(row_acronym, col_acronym, fill = r)) +
        facet_grid(col_parent ~ row_parent,
                   space = "free",
                   margins = FALSE,
                   scales = "free") +
          geom_tile() +
        geom_text(aes(label = sig_text), size=sig_size, position = position_nudge(y = sig_nudge_y), color = sig_color) +
        scale_fill_gradient2(low = "#4f4f4f",mid = "#ffffff", high = colors[[channel]],
                         aesthetics = c("color","fill"), na.value = "grey50",
                         limits=c(-1, 1)) +
        labs(title = plot_title, x = "Brain Region", y = "Brain Region") +
    theme.hm


  if (print_plot){
    quartz()
    print(p)
  }

  if(save_plot){
    # Plot the heatmap
    quartz(width = width, height = height)
    print(p)

    # Create figure directory if it doesn't already exists
    output_dir <-  file.path(attr(e, "info")$output_path, "figures")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    image_file <- file.path(output_dir, paste0("heatmap_", stringr::str_replace(plot_title, " ", "_"), "_", channel, image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
  }
    # Store the plot handle
    p_list[[channel]] <- p
  }
  return(p_list)
}



#' Plot the results of the permutation histogram used to determine the p-value of the pairwise region comparison
#'
#'
# plot_permutation_histogram <- function(e, roi = "dDG"){
# }





#' Create a Volcano plot.
#'
#' Plot the correlation difference between two comparison groups into a volcano plot. The function
#' [SMARTR::correlation_diff_permutation()] must be run first in order to generate results to plot.
#'
#' @param e experiment object
#' @param permutation_comparison The name of the correlation group comparisons to plot.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) channels to plot.
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal code for the colors corresponding to the
#' channels attribute of the correlation_list. Color values can also be input compatible with ggplot2 plotting functions.
#' @param print_plot (bool, default = TRUE) Print the plot as graphics windows.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#' the experiment object output folder.
#' @param image_ext (default = ".png") image extension to save the plot as.
#' @param title Title of the plot.
#' @param height height of the plot in inches.
#' @param ylim
#' @param plt_theme (default = NULL) Add a [ggplot2::theme()] to the plot. If NULL, the default is taken..
#' @param point_size (default = 1) Size of the plotted points.
#' @param width width of the plot in inches.
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples volcano_plot(e, permutation_comparison = "female_AD_vs_male_AD", channels = c("cfos", "eyfp", "colabel"),
#' colors =  c("#be0000", "#00782e", "#f09b08"), save_plot = TRUE, title = NULL, ylim = c(0, 3), height = 8,
#' width = 10, print_plot = TRUE, image_ext = ".png")


volcano_plot <- function(e,
                         permutation_comparison = "female_AD_vs_male_AD",
                         channels = c("cfos", "eyfp", "colabel"),
                         colors =  c("#be0000", "#00782e", "#f09b08"),
                         save_plot = TRUE,
                         title = NULL,
                         ylim = c(0, 3),
                         height = 8,
                         width = 10,
                         print_plot = TRUE,
                         plt_theme = NULL,
                         point_size = 1,
                         image_ext = ".png"){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  ## plotting theme
  if (is.null(plt_theme)){
    plt_theme <- ggplot2::theme_classic() + theme(text = element_text(size = 22),
                                                   line = element_line(size = 1),
                                                   plot.title = element_text(hjust = 0.5, size = 36),
                                                   axis.ticks.length = unit(5.5,"points"),
                                                   axis.text.x = element_text(colour = "black"),
                                                   axis.text.y = element_text(colour = "black"))
  }

  if (is.null(title)){
    title <-  permutation_comparison
  }

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){
    # Get alpha levels
    alpha <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$alpha

    # get the p-values, significance, and corr_diff
    p_vals <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$p_val %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    corr_diffs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$test_statistic %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    sigs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$sig %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")

    # Pivot long
    p_vals <- p_vals %>% tidyr::pivot_longer(col = -"rowreg", values_drop_na = TRUE,
                                      values_to = "p_val", names_to = "colreg")
    corr_diffs <- corr_diffs %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                              values_to = "corr_diff", names_to = "colreg")
    sigs <- sigs %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                  values_to = "sig", names_to = "colreg")

    # Combine into master df for plotting
    df <- p_vals %>% dplyr::inner_join(corr_diffs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(sigs, by = c("rowreg", "colreg"))


    p <- ggplot(df, aes(x = corr_diff, y = -log10(p_val))) +
      geom_point(size = point_size) +
      geom_point(data = subset(df, sig > 0 & corr_diff <= -1 | sig > 0 & corr_diff >= 1), color = colors[k]) +
      geom_vline(xintercept = c(-1, 1), color = colors[k], size = 1) +
      geom_hline(yintercept = -log10(alpha), color = colors[k], size = 1) +
      ggplot2::coord_cartesian(xlim = c(-2.1, 2.1)) +
      # xlim(c(-2.2, 2.2)) +
      # scale_x_reverse() +
      ylim(ylim) +
      labs(title = title, x = "Correlation Difference", y = "-log(p-value)") +
      plt_theme

    if (print_plot){
      quartz()
      print(p)
    }

    if(save_plot){
      # Plot the plot
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }

      image_file <- file.path(output_dir,
                              paste0("volcano_plot_", str_replace(title, " ", "_"), "_", channels[k], image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }

  return(p_list)
}



#' Create a parallel coordinate plot
#' @description Plot the correlation difference between two comparison groups into a parallel coordinate plot. The function
#' [SMARTR::correlation_diff_permutation()] must be run first in order to generate results to plot.
#'
#' @param e experiment object
#' @param permutation_comparison The name of the correlation group comparisons to plot.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) channels to plot
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal codes corresponding to the channels (respectively) to plot.
#' @param x_label_group_1 (str, NULL) The label for the first group in the permutation analysis. Note: this is to customize the graph labels. It does not reverse the group order.
#' @param x_label_group_2 (str, NULL) The label for the second group in the permutaiton analysis. Note: this is to customize the graph labels. It does not reverse the group order.
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param print_plot (bool, default = TRUE) Whether to display the plot (in addition to saving the plot)
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param reverse_group_order (bool, default = TRUE) Reverse the order of the groups on the x-axis.
#' @param image_ext (default = ".png") image extension to save the plot as.
#' @param plt_theme (default = NULL) Add a [ggplot2::theme()] to the plot. If NULL, the default is taken.
#' @param force (default =1) Force of the text repel between text labels.
#' @param nudge_x (vec, default = 2:5) a vector determining the jitter between labels.
#' @param label_size (default = 30) Default font size for region labels.
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @example
parallel_coordinate_plot <- function(e,
                                     permutation_comparison = "AD_vs_control",
                                     channels = c("cfos", "eyfp", "colabel"),
                                     colors =  c("#be0000", "#00782e", "#f09b08"),
                                     x_label_group_1 = NULL,
                                     x_label_group_2 = NULL,
                                     height = 10,
                                     width = 10,
                                     print_plot = TRUE,
                                     save_plot = TRUE,
                                     reverse_group_order= FALSE,
                                     force = 1,
                                     plt_theme = NULL,
                                     label_size = 30,
                                     image_ext = ".png",
                                     nudge_x = 2:5
                                     ){


  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # parse group names
  if(is.null(x_label_group_1)){
    group_1 <- stringr::str_split(permutation_comparison, "_vs_", simplify = TRUE)[,1] %>%
      stringr::str_split("_", simplify = TRUE) %>% paste(collapse = " ")
  } else {
    group_1 <- x_label_group_1
  }
  if(is.null(x_label_group_2)){
    group_2 <- stringr::str_split(permutation_comparison, "_vs_", simplify = TRUE)[,2] %>%
      stringr::str_split("_", simplify = TRUE) %>% paste(collapse = " ")
  } else {
    group_2 <- x_label_group_2
  }

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    # Get alpha
    alpha <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$alpha

    # get the p-values, significance, and corr_diff
    p_vals <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$p_val %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    corr_diffs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$test_statistic %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    sigs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$sig %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    group_1_pearson <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$group_1_pearson %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    group_2_pearson <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$group_2_pearson %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")


    # Pivot long
    p_vals <- p_vals %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                      values_to = "p_val", names_to = "colreg")
    corr_diffs <- corr_diffs %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                              values_to = "corr_diff", names_to = "colreg")
    sigs <- sigs %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                  values_to = "sig", names_to = "colreg")
    group_1_pearson <- group_1_pearson %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                  values_to = group_1, names_to = "colreg")
    group_2_pearson <- group_2_pearson %>% tidyr::pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                  values_to = group_2, names_to = "colreg")

    # Combine into master df for plotting
    df <- p_vals %>% dplyr::inner_join(corr_diffs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(sigs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(group_1_pearson, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(group_2_pearson, by = c("rowreg", "colreg")) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(c(group_1, group_2)), names_to = "group", values_to = "corr")



    df <- df %>% dplyr::filter(sig, abs(corr_diff) >= 1) %>%
      dplyr::mutate(group = factor(group, levels = c(group_1, group_2)),
             nudge = ifelse(group == group_1, -0.1, 0.1)) %>%
      dplyr::arrange(group, corr_diff) %>%
      mutate(text = paste(rowreg, colreg, sep = "."),
             group_plot = paste(rowreg, colreg, sep = "."))

    if (isTRUE(reverse_group_order)){
      df <- df %>% dplyr::mutate(group = fct_rev(group),
                                 nudge = ifelse(group == group_2, -0.1, 0.1))
    }

    # Error handling
    tryCatch({
      df[seq(2, nrow(df)/2, by = 2),]$text <- ""
      df[seq(nrow(df)/2+1, nrow(df), by = 2),]$text <- ""

    }, error = function(e) {
      message("could not divide by 2. Skipping")
    })

    # Plotting theme
    if (is.null(plt_theme)){
      plt_theme <- ggplot2::theme_classic() +
        theme(text = element_text(size = 22),
              line = element_line(size = 1),
              plot.title = element_text(hjust = 0.5, size = 36),
              axis.ticks.length = unit(5.5, "points"),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black")
        )
    }

    # Create parallel coordinate plot
    p <- ggplot(df, aes(x = group, y = corr, group = group_plot)) +
      ggplot2::geom_line(alpha = 0.5, color = colors[k], size = 3) +
      ggplot2::geom_point(size = 4, alpha = 0.5, color = colors[k]) +
      ggrepel::geom_text_repel(aes(label = text),
                      size = label_size,
                      color = colors[k], direction = "y",
                      force = force,
                      ylim = c(-1, 1),
                      segment.alpha = 0.3,
                      nudge_x = dplyr::pull(df, nudge)*nudge_x, max.iter = 20000) +
      ggplot2::geom_hline(yintercept = 0,linetype=2,size=1.2) +
      xlab("Group") + ylab("Correlation") +
      expand_limits(y=c(-1,1)) + plt_theme


    if (print_plot){
      quartz(height = height, width = width)
      print(p)
    }

    if(save_plot){
      # Plot
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0(permutation_comparison, "_parallel_coordinate_plot_",
                                                 channels[k], "_", image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

#' Plot the networks stored in an experiment object
#'
#' @param e experiment object
#' @param network_name (str, default = "AD")
#' @param title (str, default = NULL) Title of network plot
#' @param channels (str, default = c("cfos", "eyfp", "colabel"))
#' @param height Height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.
#'  the experiment object output folder.
#' @param edge_color (str, default = "firebrick") Color of the network edges.
#' @param degree_scale_limit (vec, default = c(1,10)) Scale limit for degree size
#' @param network_radius
#' @param graph_theme (default = NULL) Add a [ggraph::theme()] to the network graph. If NULL, the default is taken.
#' @param label_size (default = 5) Default font size for network region labels.
#' @param label_offset (default = 0.15) Distance of label from nodes.
#' @param region_legend (default = TRUE) Boolean determining whether or not to show the region legend categorizing subregions into their largest parent region. Only works well if the Allen ontology is used for the dataset.
#' @param correlation_edge_width_limit
#' @param edge_type (default = "arc") "arc" or "diagonal".
#' Can also be a hexadecimal color code written as a string.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
plot_networks <- function(e,
                          network_name = "AD",
                          title = NULL,
                          channels = c("cfos", "eyfp", "colabel"),
                          edge_color = "firebrick",
                          height = 15,
                          width = 15,
                          edge_type = "arc",
                          region_legend = TRUE,
                          degree_scale_limit = c(1,10),
                          correlation_edge_width_limit = c(0.8,1),
                          image_ext = ".png",
                          network_radius = 300,
                          print_plot = TRUE,
                          graph_theme = NULL,
                          label_size = 5,
                          label_offset = 0.15,
                          save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (channel in channels){

    network <- e$networks[[network_name]][[channel]]

    # _______________ Plot the network ________________________________
    if (is.null(graph_theme)){
      graph_theme <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5,size = 28),
                                                     legend.text = element_text(size = 15),
                                                     legend.title = element_text(size = 15))
    }

    if (edge_type == "diagonal"){
      p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
        ggraph::geom_edge_diagonal(aes(color = sign, width = abs(weight)),
                                   edge_alpha = 0.6, n = 1000)

    } else if (edge_type == "arc"){
      p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
        ggraph::geom_edge_arc(aes(color = sign, width = abs(weight)),
                              edge_alpha = 0.6, n = 1000)

    }
   p <- p + ggraph::geom_node_point(aes(size = degree,
                                    color = super.region),
                                show.legend = TRUE) +
      ggraph::geom_node_text(aes(x = (sqrt(x^2+y^2)+label_offset)*cos(atan(y/x))*sign(x),
                         y = abs((sqrt(x^2+y^2)+label_offset)*sin(atan(y/x)))*sign(y),
                         angle = atan(y/x)*180/pi,
                         label = name),
                     repel = FALSE, color = "grey25",
                     size = label_size,
                     show.legend = NA) +
      ggraph::scale_edge_color_manual(values = c(pos = edge_color,
                                                 neg = "grey20"),
                                      labels = c(pos = "Positive", neg = "Negative"),
                                      name = "Correlation",
                                      guide = guide_legend(order = 1)) +


      ggraph::scale_edge_width(limits=correlation_edge_width_limit,range = c(1,3),name = "Correlation Strength",
                       guide = guide_legend(order = 3)) +
      ggraph::scale_color_viridis(name = "Anatomical Region",
                                  discrete = TRUE,
                                  option = "D",
                                  guide = guide_legend(override.aes = list(size=5), order=4)) +
      ggplot2::scale_size(limits = degree_scale_limit, name="Degree",range=c(4,10),
                 guide = guide_legend(order = 2)) +
      ggplot2::coord_equal() + graph_theme

    if (is.null(title)){
      title <- paste(network_name, channel)
    }
    p <-  p + ggplot2::ggtitle(title)

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("network_", network_name, "_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    # Store the plot handle
    p_list[[channel]] <- p
  }
  return(p_list)
}


#' Plot the degree distributions
#' @description
#' Plot a stacked bar plot of the degree distributions.
#'
#' @param e experiment object
#' @param channels (str, default = c("cfos", "eyfp")) Channels to plot.
#' @param color_palettes (str, default = c("reds", "greens")) Color palettes from [grDevices::hcl.colors] that are used to for plotting networks for each channel, respectively.
#' @param colors_manual (str, default = NULL ) Manually choose the hexadecimal color codes to create a custom color palette, e.g. colors_manual = c("#660000", "#FF0000", "#FF6666").
#' Warning: this color will be applied to all channels. It's recommended to set the channels parameter to a single channel if this parameter is used.
#' @param title (str, default = "my_title")
#'  the experiment object output folder.
#' @param labels (e.g. labels = c(network1_name = "network 1 label", network2_name = "network 2 label))
#' The legend labels to correspond with your respective network names.
#' @param height (int, default = 15) Height of the plot in inches.
#' @param width (int, default = 15) Width of the plot in inches.
#' @param xlim (vec, default = c(0,20)) axes limits x-axis
#' @param ylim (vec, default = c(0,15))axes limits of y-axis
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param theme.gg (default = NULL) Option to use custom ggplot2 theme if the user wants
#' @param image_ext (default = ".png") image extension to the plot as.
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_degree_distributions <- function(e,
                                      channels = c("cfos", "eyfp"),
                                      color_palettes = c("reds", "greens"),
                                      colors_manual = NULL,
                                      labels = c(female_AD = "female_AD_label",
                                                 female_control = "female_control_label") ,
                                      title = "my_title",
                                      height = 15,
                                      width = 15,
                                      xlim = c(0,20),
                                      ylim = c(0,15),
                                      image_ext = ".png",
                                      print_plot = TRUE,
                                      theme.gg = NULL,
                                      save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  if (is.null(theme.gg)){
    theme.gg <- ggplot2::theme_classic() +
                theme(text = element_text(size = 22),
                      line = element_line(size = 1),
                      plot.title = element_text(hjust = 0.5, size = 36),
                      axis.ticks.length = unit(5.5,"points"))
  }


  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()


    if (!is.null(colors_manual)){
      colors <- colors_manual
    } else {
      colors <- grDevices::hcl.colors(n_groups, color_palettes[[k]])
    }

    p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_nodes,
                    aes(degree)) +
      ggplot2::geom_bar(aes(fill = group), color="black") +
      ggplot2::scale_fill_manual(values = colors[1:n_groups], name = "Group",
                        labels = labels) +
      scale_x_continuous(breaks = 1:20) +
      xlab("Degree") +
      ylab("Frequency") +
      xlim(xlim) +
      ylim(ylim) +
      theme.gg


    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }


    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_degree_distributions_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}



#' Plot the mean degree of the networks in a barplot. Error bars are plotted as SEM.
#' @param e experiment object
#' @param labels (str) The legend labels to correspond with your network names, e.g. labels = c(network1_name = "network 1 label", network2_name = "network 2 label).
#' These are the same network names used in the function [`SMARTR::summarise_networks()`].
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param color_palettes (str, default = c("reds", "greens")) Color palettes from [grDevices::hcl.colors] that are used to for plotting networks for each channel, respectively.
#' @param colors_manual (str, default = NULL ) Manually choose the hexadecimal color codes to create a custom color palette, e.g. colors_manual = c("#660000", "#FF0000", "#FF6666").
#' Warning: this will be applied to all channels. It's recommended to set the channels parameter to a single channel if this parameter is used.
#' @param title (str, default = "my_title) plot title
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#'  the experiment object output folder.
#' @param label_angle (int, default = 60)
#' @param height (int, default = 10) Height of the plot in inches.
#' @param width (int, default = 10) Width of the plot in inches.
#' @param ylim (vec, default = c(0,10)) Axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_mean_degree <- function(e,
                             color_palettes = c("reds", "greens"),
                             colors_manual = NULL,
                             channels = c("cfos", "eyfp"),
                             labels = c("AD" = "AD_label", "control" = "control_label") ,
                             title = "my_title",
                             height = 10,
                             width = 10,
                             label_angle = 60,
                             rev_x_scale = FALSE,
                             ylim = c(0, 70),
                             image_ext = ".png",
                             print_plot = TRUE,
                             save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.small.xh <- ggplot2::theme_classic() +
    theme(text = element_text(size = 22),
          line = element_line(size = 1),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
          axis.ticks.length = unit(5.5,"points"))

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()


    if (!is.null(colors_manual)){
      colors <- colors_manual
    } else {
      colors <- grDevices::hcl.colors(n_groups, color_palettes[[k]])
    }

    # Reverse the scale of the categorical variable when plotting
    if (rev_x_scale){
      p <- ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = reorder(group, dplyr::desc(group)), degree.mean))
    } else{

      p <- ggplot(e$networks_summaries[[channel]]$networks_stats, aes(group, degree.mean))
    }


    p <- p +
      ggplot2::geom_col(aes(fill = group), color = "black", size = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = degree.mean - (degree.sd/sqrt(n.nodes)),
                        ymax = degree.mean + (degree.sd/sqrt(n.nodes))),
                        width = 0.2, size = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") + ylab("Mean Degree") +
      ylim(ylim) +
      theme.small.xh

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_degree_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}



#' Plot mean clustering coefficient
#' @description
#' Plot the mean clustering coefficients of the networks in a barplot. Error bars are plotted as SEM.
#'
#' @param e experiment object
#' @param labels (e.g. labels = c(network1_name = "network 1 label", network2_name = "network 2 label)) The legend labels to correspond with your network names.
#' @param title (str, default = "my_title) plot title
#' @param color_palettes (str, default = c("reds", "greens")) Color palettes from [grDevices::hcl.colors] that are used to for plotting networks characteristics for each channel, respectively.
#' @param colors_manual (str, default = NULL ) Manually choose the hexadecimal color codes to create a custom color palette, e.g. colors_manual = c("#660000", "#FF0000", "#FF6666").
#' Warning: this will be applied to all channels. It's recommended to set the channels parameter to a single channel if this parameter is used.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 10) Height of the plot in inches.
#' @param width (int, default = 10) Width of the plot in inches.
#' @param ylim (vec, default = c(0,10)) Axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param label_angle (int, default = 60)
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_mean_clust_coeff <- function(e,
                             color_palettes = c("reds", "greens"),
                             colors_manual = NULL,
                             channels = c("cfos", "eyfp"),
                             labels = c("AD" = "AD_label", "control" = "control_label") ,
                             title = "my_title",
                             height = 10,
                             width = 10,
                             label_angle = 60,
                             rev_x_scale = FALSE,
                             ylim = c(0, 0.7),
                             image_ext = ".png",
                             print_plot = TRUE,
                             save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.small.xh <- ggplot2::theme_classic() +
    theme(text = element_text(size = 22),
          line = element_line(size = 1),
          axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.ticks.length = unit(5.5,"points"))

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()


    if (!is.null(colors_manual)){
      colors <- colors_manual
    } else {
      colors <- grDevices::hcl.colors(n_groups, color_palettes[[k]])
    }

    # Reverse the scale of the categorical variable when plotting
    if (rev_x_scale){
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(group, dplyr::desc(group)), clust.coef.mean))
    } else{

      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(group, clust.coef.mean))
    }

    p <- p +
      ggplot2::geom_col(aes(fill = group), color = "black", size = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = clust.coef.mean - (clust.coef.sd/sqrt(n.nodes)),
                                 ymax = clust.coef.mean + (clust.coef.sd/sqrt(n.nodes))),
                             width = 0.2, size = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") +
      ylab("Mean Clustering Coefficient") +
      ylim(ylim) +
      theme.small.xh


    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_clust_coeff_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

#' Plot mean global efficiency
#' @description
#' Plot the mean global efficiency of the networks in a barplot. Error bars are plotted as SEM.
#'
#' @param e experiment object
#' @param labels The labels to correspond with your network names.
#' @param title (str, default = "my_title) plot title
#' @param color_palettes (str, default = c("reds", "greens")) Color palettes from [grDevices::hcl.colors] that are used to for plotting networks characteristics for each channel, respectively.
#' @param colors_manual (str, default = NULL ) Manually choose the hexadecimal color codes to create a custom color palette, e.g. colors_manual = c("#660000", "#FF0000", "#FF6666").
#' Warning: this will be applied to all channels. It's recommended to set the channels parameter to a single channel if this parameter is used.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 10) Height of the plot in inches.
#' @param width (int, default = 10) Width of the plot in inches.
#' @param label_angle (int, default = 60)
#' @param ylim (vec, default = c(0,10)) Axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_mean_global_effic <- function(e,
                                  color_palettes = c("reds", "greens"),
                                  colors_manual = NULL,
                                  channels = c("cfos", "eyfp"),
                                  labels = c("AD" = "AD_label", "control" = "control_label") ,
                                  title = "my_title",
                                  height = 10,
                                  width = 10,
                                  label_angle = 60,
                                  rev_x_scale = FALSE,
                                  ylim = c(0, 0.7),
                                  image_ext = ".png",
                                  print_plot = TRUE,
                                  save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.small.xh <- ggplot2::theme_classic() +
    theme(text = element_text(size = 22),
          line = element_line(size = 1),
          axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.ticks.length = unit(5.5,"points"))

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()


    if (!is.null(colors_manual)){
      colors <- colors_manual
    } else {
      colors <- grDevices::hcl.colors(n_groups, color_palettes[[k]])
    }

    # Reverse the scale of the categorical variable when plotting
    if (rev_x_scale){
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(group, dplyr::desc(group)), efficiency.mean))
    } else{
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(group, efficiency.mean))
    }

    p <- p +
      ggplot2::geom_col(aes(fill = group), color = "black", size = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = efficiency.mean - (efficiency.sd/sqrt(n.nodes)),
                                 ymax = efficiency.mean + (efficiency.sd/sqrt(n.nodes))),
                             width = 0.2, size = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") +
      ylab("Global Efficiency") +
      ylim(ylim) +
      theme.small.xh

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_global_efficiency_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

#' Plot mean betweenness centrality
#' @description
#' Plot the mean betweenness centrality of the networks in a barplot. Error bars are plotted as SEM.
#' @param e experiment object
#' @param labels The labels to correspond with your network names.
#' @param title (str, default = "my_title) plot title
#' @param color_palettes (str, default = c("reds", "greens")) Color palettes from [grDevices::hcl.colors] that are used to for plotting networks characteristics for each channel, respectively.
#' @param colors_manual (str, default = NULL ) Manually choose the hexadecimal color codes to create a custom color palette, e.g. colors_manual = c("#660000", "#FF0000", "#FF6666").
#' Warning: this will be applied to all channels. It's recommended to set the channels parameter to a single channel if this parameter is used.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 10) Height of the plot in inches.
#' @param width (int, default = 10) Width of the plot in inches.
#' @param label_angle (int, default = 60)
#' @param ylim (vec, default = c(0,10)) Axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_mean_between_centrality <- function(e,
                                   color_palettes = c("reds", "greens"),
                                   colors_manual = NULL,
                                   channels = c("cfos", "eyfp"),
                                   labels = c("AD" = "AD_label", "control" = "control_label") ,
                                   title = "my_title",
                                   height = 10,
                                   width = 10,
                                   label_angle = 60,
                                   rev_x_scale = FALSE,
                                   ylim = c(0, 50),
                                   image_ext = ".png",
                                   print_plot = TRUE,
                                   save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.small.xh <- ggplot2::theme_classic() +
    theme(text = element_text(size = 22),
          line = element_line(size = 1),
          axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.ticks.length = unit(5.5,"points"))

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_stats$group %>% length()
    if (!is.null(colors_manual)){
      colors <- colors_manual
    } else {
      colors <- grDevices::hcl.colors(n_groups, color_palettes[[k]])
    }

    # Reverse the scale of the categorical variable when plotting
    if (rev_x_scale){
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(group, dplyr::desc(group)), btw.mean))
    } else{
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(group, btw.mean))
    }

    p <- p +
      ggplot2::geom_col(aes(fill = group), color = "black", size = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = btw.mean - (btw.sd/sqrt(n.nodes)),
                                 ymax = btw.mean + (btw.sd/sqrt(n.nodes))),
                             width = 0.2, size = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") +
      ylab("Mean Betweenness Centrality") +
      ylim(ylim) +
      theme.small.xh

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_betweenness_centrality_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}




#' Plot the degree distributions across regions
#' @description
#' Bar plot of degree per region in descending magnitude
#'
#' @param e experiment object
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 15) Height of the plot in inches.
#' @param width (int, default = 20) Width of the plot in inches.
#' @param ylim (vec, default = c(0,15))axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param colors (str, default = ) String vector of hexadecimal color codes corresponding to to each channel plotted.
#' @param network (str, default = "AD") Which network to plot the degree distribution across regions
#' @param title (str, default = "")
#' @param sort_super_region (bool, default = FALSE) Whether to divide into subfacets based on which parent region
#' @param region_label_angle (int, default = 60) Angle of region labels.
#' @param label_text_size (int, default = 12) Font size of region labels.
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_degree_regions <- function(e,
                                channels = c("cfos", "eyfp"),
                                colors = c("red", "green"),
                                network = "AD",
                                title = "",
                                height = 10,
                                width = 20,
                                ylim = c(0, 15),
                                sort_super_region = FALSE,
                                region_label_angle = 60,
                                label_text_size = 12,
                                image_ext = ".png",
                                print_plot = TRUE,
                                save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.bar <-  ggplot2::theme_classic() +
    theme(axis.text.x = element_text(angle = region_label_angle, hjust = 1, size = label_text_size, color = "black"),
          axis.text.y = element_text(color = "black", size = 20),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 24),
          strip.clip = "off",
          strip.text.x = element_text(angle = 0, margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")),
          strip.placement = "outside",
          strip.background = element_rect(color = "black",
                                          fill = "lightblue"),
          plot.margin = margin(1,1.5,0,1.5, "cm"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()

    if (sort_super_region) {
      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(group == network) %>%
        dplyr::arrange(super.region,  dplyr::desc(degree)) %>%
        dplyr::mutate(name = factor(name, levels = name))

      p <- df %>%
        ggplot2::ggplot(aes(name, degree)) +
        ggplot2::geom_col(fill = colors[[k]], color = "black") +
        facet_grid(~super.region, scales = "free_x", space = "free_x", switch = "x")  +
        xlab("Brain Region") + ylab("Degree") +
        ylim(ylim)  + theme.bar

    } else {

      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(group == network) %>%
        dplyr::arrange(dplyr::desc(degree)) %>%
        dplyr::mutate(name = factor(name, levels = name))

      p <- df %>%
        ggplot2::ggplot(aes(name, degree)) +
        ggplot2::geom_col(fill = colors[[k]], color = "black") +
        xlab("Brain Region") + ylab("Degree") +
        ylim(ylim) +
        theme.bar
    }

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }


    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_degree_per_region_", channel, "_",network, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}



#' Plot the betweenness distributions across regions
#' @description
#' Bar plot of betweenness per region in descending magnitude
#'
#' @param e experiment object
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 15) Height of the plot in inches.
#' @param width (int, default = 20) Width of the plot in inches.
#' @param ylim (vec, default = c(0,15))axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param colors (str, default = ) String vector of hexadecimal color codes corresponding to to each channel plotted.
#' @param network (str, default = "AD") Which network to plot the betweenness distribution across regions
#' @param title (str, default = "")
#' @param sort_super_region (bool, default = FALSE) Whether to divide into subfacets based on which parent region
#' @param region_label_angle (int, default = 60) Angle of region labels.
#' @param label_text_size (int, default = 12) Font size of region labels.
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_betweenness_regions <- function(e,
                                     channels = c("cfos", "eyfp"),
                                     colors = c("red", "green"),
                                     network = "AD",
                                     title = "",
                                     height = 10,
                                     width = 20,
                                     ylim = c(0, 15),
                                     sort_super_region = FALSE,
                                     region_label_angle = 60,
                                     label_text_size = 12,
                                     image_ext = ".png",
                                     print_plot = TRUE,
                                     save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.bar <-  ggplot2::theme_classic() +
    theme(axis.text.x = element_text(angle = region_label_angle, hjust = 1, size = label_text_size, color = "black"),
          axis.text.y = element_text(color = "black", size = 20),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 24),
          strip.clip = "off",
          strip.text.x = element_text(angle = 0, margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")),
          strip.placement = "outside",
          strip.background = element_rect(color = "black",
                                          fill = "lightblue"),
          plot.margin = margin(1,1.5,0,1.5, "cm"),
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()

    if (sort_super_region) {
      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(group == network) %>%
        dplyr::arrange(super.region,  dplyr::desc(btw)) %>%
        dplyr::mutate(name = factor(name, levels = name))

      p <- df %>%
        ggplot2::ggplot(aes(name, btw)) +
        ggplot2::geom_col(fill = colors[[k]], color = "black") +
        facet_grid(~super.region, scales = "free_x", space = "free_x", switch = "x")  +
        xlab("Brain Region") + ylab("Betweenness") +
        ylim(ylim)  + theme.bar

    } else {

      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(group == network) %>%
        dplyr::arrange(dplyr::desc(btw)) %>%
        dplyr::mutate(name = factor(name, levels = name))

      p <- df %>%
        ggplot2::ggplot(aes(name, btw)) +
        ggplot2::geom_col(fill = colors[[k]], color = "black") +
        xlab("Brain Region") + ylab("Betweenness") +
        ylim(ylim) +
        theme.bar
    }

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_betweenness_per_region_", channel, "_", network, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}



#________________ Internal Analysis functions _______________


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

sort_anatomical_order <- function(common_regions, ontology ="allen",
                                  anatomical.order = c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")){

  if (tolower(ontology) == "allen"){
    common_regions <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure) %>%
      purrr::map(intersect,y=common_regions) %>% unlist()
  } else {
    common_regions <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure.custom, ontology=ontology) %>%
      purrr::map(intersect,y=common_regions) %>% unlist()
  }
  return(common_regions)

  }


#' Plot the networks stored in an experiment object
#' @param e experiment object
#' @param correlation_list_names (str vec) character vector of the two correlation lists used to include in a joined network, e.g., `correlation_list_names = c("male_agg", "female_non")`
#' @param title (str, default = NULL) Title of network plot
#' @param channels (str, default = c("cfos", "eyfp", "colabel"))
#' @param height Height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.
#'  the experiment object output folder.
#' @param absolute_weight (bool, default = TRUE) Whether to plot absolute weights. If TRUE, the edge_colors and edge_colors_label should not contain values for positive and negative correlations.
#' @param edge_color (str, default =  c(male_agg_pos = "Positive male", male_agg_neg = "Negative male", female_non_pos = "Positive female", female_non_neg = "Negative female")) Color of the network edges as a named vector.
#' @param degree_scale_limit (vec, default = c(1,10)) Scale limit for degree size
#' @param network_radius
#' @param graph_theme (default = NULL) Add a [ggraph::theme()] to the network graph. If NULL, the default is taken.
#' @param label_size (default = 5) Default font size for network region labels.
#' @param label_offset (default = 0.15) Distance of label from nodes.
#' @param region_legend (default = TRUE) Boolean determining whether or not to show the region legend categorizing subregions into their largest parent region. Only works well if the Allen ontology is used for the dataset.
#' @param correlation_edge_width_limit
#' Can also be a hexadecimal color code written as a string.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
plot_joined_networks <- function(e,
                                 correlation_list_names = c("male_agg", "female_non"),
                                 title = NULL,
                                 channels = "cfos",
                                 absolute_weight = TRUE,
                                 edge_colors = c(male_agg_pos =  "#06537f",
                                                 male_agg_neg = "#526c7a",
                                                 female_non_pos = "#C70039",
                                                 female_non_neg = "#71585f"),
                                 edge_color_labels = c(male_agg_pos = "Positive male",
                                                       male_agg_neg = "Negative male",
                                                       female_non_pos = "Positive female",
                                                       female_non_neg = "Negative female"),
                                 height = 15,
                                 width = 15,
                                 region_legend = TRUE,
                                 degree_scale_limit = c(1,10),
                                 correlation_edge_width_limit = c(0.8,1),
                                 image_ext = ".png",
                                 network_radius = 300,
                                 print_plot = TRUE,
                                 graph_theme = NULL,
                                 reverse_plot_edge_order = FALSE,
                                 transparent_edge_group1 = TRUE,
                                 transparent_edge_group2 = FALSE,
                                 label_size = 5,
                                 label_offset = 0.15,
                                 save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (channel in channels){

    joined_network_name <- paste(correlation_list_names, collapse = "_")
    network <- e$networks[[joined_network_name]][[channel]]
    p1 <-  correlation_list_names[1]
    p2 <-  correlation_list_names[2]


    # _______________ Plot the network ________________________________
    if (is.null(graph_theme)){
      graph_theme <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5,size = 28),
                                                   legend.text = element_text(size = 15),
                                                   legend.title = element_text(size = 15)
      )

    }

    if (isTRUE(absolute_weight)) {
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(color = factor(network))
    } else {
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(color=paste(network, sign, sep = "_"))
    }

    if (isTRUE(reverse_plot_edge_order)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(color = forcats::fct_rev(color))
    }

    if (isTRUE(transparent_edge_group1) &&  isTRUE(transparent_edge_group2)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = 0)
    } else if (isFALSE(transparent_edge_group1) && isTRUE(transparent_edge_group2)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = if_else(network == p2, 0.0, 0.6))
    } else if (isTRUE(transparent_edge_group1) && isFALSE(transparent_edge_group2)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = if_else(network == p1, 0.0, 0.6))
    } else {
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = 0.6)
    }


    p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
      # ggraph::geom_edge_diagonal(aes(color = color, width = abs(weight)),
      #                            edge_alpha = 0.6, n = 1000) +
      ggraph::geom_edge_arc(aes(color = color, width = abs(weight), edge_alpha = edge_alpha),
                            n = 1000) +
      ggraph::geom_node_point(aes(size = degree,
                                  color = super.region),
                              show.legend = TRUE) +
      ggraph::geom_node_text(aes(x = (sqrt(x^2+y^2)+label_offset)*cos(atan(y/x))*sign(x),
                                 y = abs((sqrt(x^2+y^2)+label_offset)*sin(atan(y/x)))*sign(y),
                                 angle = atan(y/x)*180/pi,
                                 label = name),
                             repel = FALSE, color = "grey25",
                             size = label_size,
                             show.legend = NA) +
      ggraph::scale_edge_color_manual(values = edge_colors,
                                      labels = edge_color_labels,
                                      name = "Correlation",
                                      guide = guide_legend(order = 1)) +
      ggraph::scale_edge_width(limits=correlation_edge_width_limit,range = c(1,4),name = "Correlation Strength",
                               guide = guide_legend(order = 3)) +
      ggraph::scale_color_viridis(name = "Anatomical Region",
                                  discrete = TRUE,
                                  option = "D",
                                  guide = guide_legend(override.aes = list(size=5), order=4)) +
      ggplot2::scale_size(limits = degree_scale_limit, name="Degree",range=c(4,10),
                          guide = guide_legend(order = 2)) +
      ggplot2::coord_equal() + graph_theme

    if (is.null(title)){
      title <- paste(network_name, channel)
    }
    p <-  p + ggplot2::ggtitle(title)

    if (print_plot){
      quartz(width = width, height = height)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      quartz(width = width, height = height)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("network_", joined_network_name, "_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    # Store the plot handle
    p_list[[channel]] <- p
  }
  return(p_list)
}






#' Generate array of null distribution of region pairwise correlation differences.
#' @param df
#' @param correlation_list_name_1
#' @param correlation_list_name_2
#' @param n_shuffle
#' @param seed random seed for replication
#' @param ...
#' @return
#'
#' @examples
permute_corr_diff_distrib <- function(df, correlation_list_name_1, correlation_list_name_2,
                                      n_shuffle = n_shuffle, seed = 5, ...){


  # Set the random seed
  set.seed(seed)

  # Create a 3D matrix to hold the correlation distributions
  n_reg <- length(names(df)) - 2
  region_names <- names(df)[3:length(names(df))]
  corr_diff_matrix <- array(dim= c(n_reg, n_reg, n_shuffle))
  dimnames(corr_diff_matrix) <- list(region_names, region_names, 1:n_shuffle)

  # Get permutation sampling combinations
  # n_mice <- length(df$mouse_ID)
  # possible_perm_combs <-  RcppAlgos::permuteGeneral(n_mice, m = n_mice, Parallel = TRUE, repetition = FALSE, ...)
  # sampled_perm_combs <- possible_perm_combs[sample(factorial(n_mice), size = n_shuffle, replace = TRUE),]

  # Get original group order
  corr_groups <- df$corr_group

  for (n in 1:n_shuffle){

    # reorder the group labels based on the permutation analysis
    # comb <- sampled_perm_combs[n, ]
    # df$mouse_ID <- df$mouse_ID[comb]
    # df$corr_group <- df$corr_group[comb]

    # Shuffle the group labels
    df$corr_group <- sample(corr_groups, replace = FALSE)

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

    correlations_list[[correlation_list_name_1]] <- matrix_list[[correlation_list_name_1]][,-1] %>% try_correlate()
    correlations_list[[correlation_list_name_2]] <- matrix_list[[correlation_list_name_2]][,-1] %>% try_correlate()

    # subtract R coefficient differences
    corr_diff_matrix[,,n] <- correlations_list[[correlation_list_name_2]]$r - correlations_list[[correlation_list_name_1]]$r
  }
  return(corr_diff_matrix)

}

#' Get  a list of intersecting regions to a list of common regions
#' @param common_reg A comprehensive list of all regions (and all existing subregions) that the brain area in the `rois` list will be compared against.
#' @param rois A list of rois whose regions and subregions will be compared the the `common_reg` list
#' @return common_reg A list of the rois (or any of its subregions) that intersected with the common_reg list.
#' @examples
#'
rois_intersect_region_list <- function(common_reg, rois){
  # Get rois of all the child regions
  child_r <- SMARTR::get.acronym.child(rois)
  while (length(child_r) > 0){
    rois <- c(rois, child_r)
    child_r <- SMARTR::get.acronym.child(child_r) %>% na.omit()
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

#' Standard error function
#' @param x (vec)
#' @return
#' @export
#' @examples
sem <- function(x){
  sd(x) / sqrt(length(x))
}



#' Try to correlate
#'
#' @param df_channel
#'
#' @return
#'
#' @examples
try_correlate <- function(df_channel){
  tryCatch({
    # Code that may throw an error
    df_corr <- df_channel %>% Hmisc::rcorr()
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
        ct <- cor.test(df_channel[,r1], df_channel[,r2])
        df_corr$P[r1,r2] <- ct$p.value
        df_corr$r[r1,r2] <- ct$estimate
      }
    }
    # df_corr is lower triangle.
    # need to ,mirror to upper tri
    upper_tri <- upper.tri(df_corr$r)
    df_corr$P[upper_tri] <- t(df_corr$P)[upper_tri]
    df_corr$r[upper_tri] <- t(df_corr$r)[upper_tri]
    df_corr$n[upper_tri] <- t(df_corr$n)[upper_tri]
    return(df_corr)
  })
}


