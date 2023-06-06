#' @importFrom ggplot2 ggplot aes theme element_text unit geom_tile geom_text scale_fill_gradient2 labs ggsave guide_legend xlab ylab xlim ylim
NULL

#' @importFrom tidyselect all_of
NULL

#' @importFrom dplyr n mutate summarize summarise across arrange group_by
NULL

#' @importFrom tidygraph activate
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
    col_to_summarize <- c("count.colabel", "count.eyfp", "colabel_percentage")

    master_counts <- master_counts %>%
      dplyr::group_by(across(all_of(attr_group_by))) %>%
      dplyr::summarise(n = dplyr::n(),
                       colabel_percentage.sd = sd(colabel_percentage),
                       colabel_percentage.sem = colabel_percentage.sd/sqrt(n),
                       across(all_of(col_to_summarize), mean),) %>%
      dplyr::rename_with(function(x){paste0(x,".mean")}, all_of(col_to_summarize)) %>%
      dplyr::relocate(colabel_percentage.sd, colabel_percentage.sem, n, .after = last_col())

    indiv_or_avg <- "average"
  } else {
    indiv_or_avg <- "individual"
  }

  if (save_table){
      out_path <- file.path(e_info$output_path, paste0(make.names(colabel_channel),"_percentage_", channel, "_", indiv_or_avg, ".csv"))
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
#' e.g.values = c("female", "AD").
#' @param channels (str, channels =  c("cfos", "eyfp", "colabel") The channels to process.
#' @param p_adjust_method (bool or str, default = "BH") Benjamini-Hochberg method is recommended.
#'  Apply the named method to control for the inflated false discovery rate or family wise error rate (FWER). Set to FALSE or "none"
#'  to keep "raw" p values. See also [stats::p.adjust()] for the correction options.
#' @param alpha (num, default = 0.05) The alpha level for significance applied AFTER p-adjustment.
#' @return e experiment object. The experiment object now has a named `correlation_list` object stored in it.
#' The name of the correlation object is the concatenation of the variable values separated by a "_".
#' This name allows for unambiguous identification of different analysis subgroups in the future.
#' @export
#' @examples e <- get_correlations(e, by = c("sex", "group"), values = c("female", "AD"),
#' channels = c("cfos", "eyfp", "colabel"),  p_adjust_method = "BH", alpha = 0.05)
#' @seealso [Hmisc::rcorr()]
get_correlations <- function(e, by, values,
                             channels = c("cfos", "eyfp", "colabel"),  p_adjust_method = "BH", alpha = 0.05){
  corr_list <- list()
  names(corr_list) <- names(channels)

  for (channel in channels){

    # filter by parameters
    df_channel <-  e$combined_normalized_counts[[channel]] %>%
      filter_df_by_char_params(by, values) %>% dplyr::distinct()

    # Pivot wider
    df_channel <- df_channel %>%  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
      tidyr::pivot_wider(names_from = acronym, values_from = normalized.count.by.volume)

    # Rearrange the correlations to be in "anatomical order"
    anatomical.order <- c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")
    common.regions <-  df_channel %>% dplyr::select(-all_of(c('mouse_ID', by))) %>% names()
    common.regions.ordered <- anatomical.order %>% purrr::map(SMARTR::get.sub.structure) %>%
      purrr::map(intersect, y=common.regions) %>% unlist()


    df_channel <- df_channel %>% dplyr::select(all_of(c(common.regions.ordered))) %>%
      as.matrix()

    # Get minimum number of mice that we have data for across all regions
    n <- colSums(!is.na(df_channel)) %>% min()

    if (n > 4){
      # Select the order of the columns, perform the correlations, ignore the mouse_ID and group columns
      df_corr <- df_channel %>% Hmisc::rcorr()
    } else {
      message(c("One or more of your brain regions has a n below 5. \nCalculating pearsons, but we ",
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
    }

    # Adjust the correlation matrix (don't include p-values for correlation against self)
    lowertri <- df_corr$P %>%  lower.tri(diag = FALSE)
    df_corr$P[!lowertri] <- NA

    lowertri <- df_corr$P %>%  lower.tri(diag = TRUE)
    df_corr$n[!lowertri] <- NA
    df_corr$r[!lowertri] <- NA

    # adjust the p-value for false discovery rate or FWER
    if (!isFALSE(p_adjust_method)){
      # Calculate without removing NAs
      rows <- rownames(df_corr$P)
      cols <- colnames(df_corr$P)
      df_corr$P <- df_corr$P %>% p.adjust(method = p_adjust_method) %>%
        matrix(nrow = length(rows), ncol= length(cols), dim = list(rows, cols))
    }

    df_corr$sig <- df_corr$P <= alpha
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
#' running this function.
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

    # Get an array of distribution of correlation differences
    test_statistic_distributions <- permute_corr_diff_distrib(df_channel_groups,
                                                         correlation_list_name_1 = correlation_list_name_1,
                                                         correlation_list_name_2 = correlation_list_name_2,
                                                         n_shuffle = n_shuffle,
                                                         seed = seed, ...)
    # For each pairwise distribution, sort the values
    test_statistic_distributions <- apply(test_statistic_distributions, 1:2, sort)
    # do.call(sort, 1:2)
    # test_statistic_distributions <- test_statistic_distributions %>% aperm(c(1, 2, 3))
    # aperm(c(3, 2, 1))

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

#' Create %>% graph objects for plotting different analysis subgroups.
#' @param e experiment object
#' @param correlation_list_name (str) Name of the correlation list object used to generate the networks.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param alpha (float, default = 0.05) The significance cutoff for including brain regions in the network.
#'
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
                            alpha = 0.05){

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
    edges <- df %>% dplyr::mutate(weight = abs(r),
                           sign = ifelse(r >= 0, "pos", "neg")) %>%
       dplyr::rename(from = row_acronym,
             to = col_acronym,
             p.value = P)

    # ___________Create Node tibble ________________
    anatomical.order <- c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")

    # Get unique common regions
    acronyms <-  e$correlation_list[[correlation_list_name]][[channel]]$r %>% rownames()

    # get the parent super region
    super.region <- acronyms
    for (sup.region in anatomical.order){
      super.region[super.region %in% SMARTR::get.sub.structure(sup.region)] <- sup.region
    }

    # get the color
    color <- wholebrain::color.from.acronym(super.region)
    nodes <- tibble::tibble(name = acronyms, super.region = super.region, color = color)

    # _____________ Create the network ________________
    network <- tidygraph::tbl_graph(nodes = nodes,
                                    edges = edges,
                                    directed = FALSE)

    # filter by alpha
    network <- network %>% activate(edges) %>% dplyr::filter(p.value < alpha)
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
#' @param save_betweenness_distribution (bool, default = TRUE) Save the betweenness distribution.
#' @return e experiment object
#' @export
#' @examples e <- get_network_statistics(e,  network_names = c("female_AD", "female_control"),
#' channels = c("cfos", "eyfp", "colabel"), save_stats = TRUE, save_degree_distribution = TRUE)
summarise_networks <- function(e,
                               network_names,
                               channels = c("cfos", "eyfp", "colabel"),
                               save_stats = TRUE,
                               save_degree_distribution = TRUE,
                               save_betweenness_distribution = TRUE){

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

  # Get output path
  outpath <- attr(e, "info")$output_path

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
      dplyr::summarise_if(is.numeric,list(~mean(.),~sd(.))) %>%
      dplyr::left_join(nodes_df %>% dplyr::group_by(group) %>% dplyr::summarise(n.nodes = n())) %>%
      dplyr::left_join(edges_df %>% dplyr::group_by(group) %>% dplyr::summarise(n.edges = n())) %>%
      dplyr::mutate(edge.density = n.edges/(2*n.nodes*(n.nodes-1))) %>%
      dplyr::rename_all(stringr::str_replace,pattern = "_",replacement=".")

    #make table of degree frequency (useful for plotting degree histogram outside of R)
    degree_distribution_df <- nodes_df %>% dplyr::group_by(group, degree) %>% dplyr::count(degree)

    # Make a dataframe of betweenness (useful)
    betweenness_distribution <- nodes_df %>% dplyr::select(name, group, btw) %>%
      dplyr::group_by(group) %>% dplyr::arrange(group, dplyr::desc(btw)) %>%
      dplyr::mutate(name = factor(name, levels=name))

    if(save_stats){
      write.csv(network_stats_df,  file.path(outpath, paste0("summary_networks_stats_", channel, ".csv")))
    }

    if(save_degree_distribution){
      write.csv(degree_distribution_df,  file.path(outpath, paste0("networks_degree_distributions_", channel, ".csv")))
    }

    if(save_betweenness_distribution){
      write.csv( betweenness_distribution,  file.path(outpath, paste0("networks_betweenness_distributions_", channel, ".csv")))
    }

    # Store the network summary data into channels
    e$networks_summaries[[channel]] <- list(networks_nodes = nodes_df,
         networks_edges = edges_df,
         networks_stats = network_stats_df,
         networks_degree_distrib = degree_distribution_df)
  }
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
    p <-ggplot(e$colabel_percent[[channel]]$average %>% dplyr::filter(acronym %in% rois),
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
#' @param width (int) Width of the plot in inches.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples plot_correlation_heatmaps(e, correlation_list_name = "female_AD") # No return value
#' @seealso [SMARTR::get_correlations()]

plot_correlation_heatmaps <- function(e, correlation_list_name ,
                                      channels = c('cfos', 'eyfp', 'colabel'),
                                      colors = c("#be0000", "#00782e", "#f09b08"),
                                      print_plot = TRUE, save_plot = TRUE, image_ext = ".png",
                                      plot_title = NULL, height = 10, width = 10){

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

  # Create plotting theme for the heatmap
  theme.hm <- ggplot2::theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 8),
                    axis.text.y = element_text(vjust = 0.5, size = 8),
                    plot.title = element_text(hjust = 0.5, size = 36),
                    axis.title = element_text(size = 18),
                    legend.text = element_text(size = 22),
                    legend.key.height = unit(100, "points"),
                    legend.title = element_text(size = 22))

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

  # Generate a correlation heatmap in anatomical order
  p <-  ggplot(df, aes(row_acronym, col_acronym, fill = r)) +
        geom_tile() +
        geom_text(aes(label = sig_text), size=8, color = "yellow") +
        scale_fill_gradient2(low = "#4f4f4f",mid = "#ffffff", high = colors[[channel]],
                         aesthetics = c("color","fill"), na.value = "grey50")+
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
    image_file <- file.path(output_dir, paste0("heatmap_", str_replace(plot_title, " ", "_"), "_", channel, image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
  }
    # Store the plot handle
    p_list[[channel]] <- p
  }
}


###########################################################################################
#  Breakpoint 4.7.2022
###########################################################################################

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
#' @param width width of the plot in inches.
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
                         image_ext = ".png"){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  ## plotting theme
  plot_theme <- ggplot2::theme_classic() + theme(text = element_text(size = 22),
                                            line = element_line(size = 1),
                                            plot.title = element_text(hjust = 0.5, size = 36),
                                            axis.ticks.length = unit(5.5,"points"))
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
      geom_point() +
      geom_point(data = subset(df, sig > 0 & corr_diff <= -1 | sig > 0 & corr_diff >= 1), color = colors[k]) +
      geom_vline(xintercept = c(-1, 1), color = colors[k], size = 1) +
      geom_hline(yintercept = -log10(alpha), color = colors[k], size = 1) +
      xlim(c(-2, 2)) +
      ylim(ylim) +
      labs(title = title, x = "Correlation Difference", y = "-log(p-value)") +
      plot_theme

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
}



#' Plot normalized cell counts
#' @description Plot the cell counts normalized by volume for a given channel
#' @param e experiment object
#' @param channels (str, default = c("cfos", "eyfp", "colabel"))
#' @param colors (str, default = c()) Hexadecimal codes corresponding to the channels (respectively) to plot.
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param print_plot (bool, default = TRUE) Whether to display the plot (in addition to saving the plot)
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @example
plot_normalized_counts <- function(e,
                                   channels = c("cfos", "eyfp", "colabel"),
                                   title = NULL,
                                   colors = c("#FFFFFF", "lightblue"),
                                   height = 7,
                                   width = 20,
                                   print_plot = TRUE,
                                   save_plot = TRUE) {


  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (channel in channels) {

  # select normalized counts for the given channel and generate mean and sem stats by region
    channel_counts <- lh$combined_normalized_counts[[channel]] %>%
      select(group, mouse_ID, name, acronym, normalized.count.by.volume) %>%
      group_by(group, acronym, name) %>%
      summarise(n = n(),
                mean_normalized_counts = mean(normalized.count.by.volume),
                sem = sd(normalized.count.by.volume)/sqrt(n))

    # generate cell counts plot for the given channel (with standard error bars, slanted region names and clean theme)

    p <- channel_counts %>%
      ggplot(aes(y = mean_normalized_counts, x = name,
                 fill = group), color = "black") +
      geom_col(position = position_dodge(0.8), width = 0.8, color = "black") +
      geom_errorbar(aes(ymin = mean_normalized_counts - sem,
                        ymax = mean_normalized_counts + sem, x = name),
                    position = position_dodge(0.8),
                    width = 0.5,
                    color = "black") +
      labs(title = title,
           y = bquote('Normalized cell counts '('cells/mm'^3)),
           x = "",
           fill = "Group") +
      scale_y_continuous(expand = c(0,0)) +
      scale_fill_manual(values=colors) +
      theme_bw() +
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 50, hjust = 1))

    p <-  p + ggplot2::ggtitle(title)

    # print plot if indicated

    if (print_plot) {
      quartz()
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

    # Save the plot handle
    p_list[[channels[channel]]] <- p

  }

}


#' Create a parallel coordinate plot
#' @description Plot the correlation difference between two comparison groups into a parallel coordinate plot. The function
#' [SMARTR::correlation_diff_permutation()] must be run first in order to generate results to plot.
#' @param e experiment object
#' @param permutation_comparison The name of the correlation group comparisons to plot.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) channels to plot
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal codes corresponding to the channels (respectively) to plot.
#' @param x_label_group_1 (str, NULL) The label for the first group in the permutation analysis.
#' @param x_label_group_2 (str, NULL) The label for the second group in the permutaiton analysis.
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param print_plot (bool, default = TRUE) Whether to display the plot (in addition to saving the plot)
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param image_ext (default = ".png") image extension to save the plot as.
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
                                     image_ext = ".png"){


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
      dplyr::mutate(group = factor(group, levels = c(group_2, group_1)),
             nudge = ifelse(group == group_1, 0.1, -0.1)) %>%
      dplyr::arrange(group, corr_diff) %>%
      mutate(text = paste(rowreg, colreg, sep = "."),
             group_plot = paste(rowreg, colreg, sep = "."))


    df[seq(2, nrow(df)/2, by = 2),]$text <- ""
    df[seq(nrow(df)/2+1, nrow(df), by = 2),]$text <- ""

    # Plotting theme
    theme.small.xh <- ggplot2::theme_classic() +
      theme(text = element_text(size = 22), line = element_line(size = 1),
            plot.title = element_text(hjust = 0.5, size = 36), axis.ticks.length = unit(5.5, "points"))

    # Create parallel coordinate plot
    p <- ggplot(df, aes(x = group, y = corr, group = group_plot)) +
      ggplot2::geom_line(alpha = 0.5, color = colors[k], size = 3) +
      ggplot2::geom_point(size = 4, alpha = 0.5, color = colors[k]) +
      ggrepel::geom_text_repel(aes(label = text),
                      color = colors[k], direction = "y", force = 1,
                      ylim = c(-1, 1),
                      segment.alpha = 0.3,
                      nudge_x = dplyr::pull(df, nudge)*1:5, max.iter = 20000) +
      ggplot2::geom_hline(yintercept = 0,linetype=2,size=1.2) +
      xlab("Group") + ylab("Correlation") +
      expand_limits(y=c(-1,1)) + theme.small.xh


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
}

#' Plot the networks stored in an experiment object
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
                          image_ext = ".png",
                          print_plot = TRUE,
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
    theme.network <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5,size = 28),
                                           legend.text = element_text(size = 15),
                                           legend.title = element_text(size = 15))

    p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
      ggraph::geom_edge_diagonal(aes(color = sign, width = abs(weight)),
                                 edge_alpha = 0.6, n = 1000) +
      ggraph::geom_node_point(aes(size = degree,
                                  color = super.region)) +
      ggraph::geom_node_text(aes(x = (sqrt(x^2+y^2)+0.1)*cos(atan(y/x))*sign(x),
                         y = abs((sqrt(x^2+y^2)+0.1)*sin(atan(y/x)))*sign(y),
                         angle = atan(y/x)*180/pi,
                         label = name),
                     repel = FALSE, color = "grey25",
                     size = 5) +
      ggraph::scale_edge_color_manual(values = c(pos = edge_color,
                                         neg = "grey20"),
                              labels = c(pos = "Positive", neg = "Negative"),
                              name = "Correlation",
                              guide = guide_legend(order = 1,
                                                   override.aes = list(size = 4))) +
      ggraph::scale_color_viridis(name = "Anatomical Region",
                          discrete = TRUE,
                          option = "D",
                          guide = guide_legend(override.aes = list(size=5), order=4)) +
      ggraph::scale_edge_width(limits=c(0.8,1),range = c(1,3),name = "Correlation Strength",
                       guide = guide_legend(order = 3)) +
      ggplot2::scale_size(limits = c(1,6),name="Degree",range=c(4,10),
                 guide = guide_legend(order = 2)) +
      ggplot2::coord_equal() + theme.network

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
}


#' Plot the degree distributions
#' @description
#' Plot a stacked bar plot of the degree distributions.
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
#' @param image_ext (default = ".png") image extension to the plot as.
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
                                      save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.small.xh <- ggplot2::theme_classic() +
    theme(text = element_text(size = 22),
          line = element_line(size = 1),
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
      image_file <- file.path(output_dir, paste0("networks_degree_distributions_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
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
}







#' Plot mean betweenness centrality
#' @description
#' Plot the mean betweenness centrality of the networks in a barplot. Error bars are plotted as SEM.
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
#' @param title
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_degree_regions <- function(e,
                                     channels = c("cfos", "eyfp"),
                                     colors = c("red", "green"),
                                     network = "AD",
                                     title = "my_title",
                                     height = 10,
                                     width = 20,
                                     ylim = c(0, 15),
                                     image_ext = ".png",
                                     print_plot = TRUE,
                                     save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.bar <-  ggplot2::theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 12),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 24))

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()

    p <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(group == network) %>%
      dplyr::arrange(dplyr::desc(degree)) %>%
      dplyr::mutate(name = factor(name, levels = name)) %>%
      ggplot2::ggplot(aes(name, degree)) +
      ggplot2::geom_col(fill = colors[[k]], color = "black") +
      xlab("Brain Region") + ylab("Degree") +
      ylim(ylim) +
      theme.bar


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
      image_file <- file.path(output_dir, paste0("networks_degree_per_region_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
}



#' Plot the betweenness distributions across regions
#' @description
#' Bar plot of betweenness per region in descending magnitude
#'
#' @param e experiment object
#' @param network (str) Which network to plot the betweenness distribution across regions
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 15) Height of the plot in inches.
#' @param width (int, default = 20) Width of the plot in inches.
#' @param ylim (vec, default = c(0,15))axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = TRUE) Whether to print the plot as an output.s
#' @param colors (str, default = ) String vector of hexadecimal color codes corresponding to to each channel plotted.
#' @param title (str, default = "my_title")
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples

plot_betweenness_regions <- function(e,
                                     network,
                                     channels = c("cfos", "eyfp"),
                                     colors = c("red", "green"),
                                     title = "my_title",
                                     height = 10,
                                     width = 20,
                                     ylim = c(0, 60),
                                     image_ext = ".png",
                                     print_plot = TRUE,
                                     save_plot = TRUE){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  theme.bar <-  ggplot2::theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 12),
          axis.text.y = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 36),
          axis.title = element_text(size = 24),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 24))

  # List to store the returned plot handles
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){

    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()

    p <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(group == network) %>%
      dplyr::arrange(dplyr::desc(btw)) %>%
      dplyr::mutate(name = factor(name, levels = name)) %>%
      ggplot2::ggplot(aes(name, btw)) +
      ggplot2::geom_col(fill = colors[[k]], color = "black") +
      xlab("Brain Region") + ylab("Betweenness") +
      ylim(ylim) +
      theme.bar

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
      image_file <- file.path(output_dir, paste0("networks_betweenness_per_region_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    # Store the plot handle
    p_list[[channels[k]]] <- p
  }
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
#' @param seed random seed for replication
#' @param ...
#' @return
#'
#' @examples
permute_corr_diff_distrib <- function(df, correlation_list_name_1, correlation_list_name_2, n_shuffle = n_shuffle, seed = 5, ...){


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

    correlations_list[[correlation_list_name_1]] <- matrix_list[[correlation_list_name_1]][,-1] %>% Hmisc::rcorr()
    correlations_list[[correlation_list_name_2]] <- matrix_list[[correlation_list_name_2]][,-1] %>% Hmisc::rcorr()

    # subtract R coefficient differences
    corr_diff_matrix[,,n] <- correlations_list[[correlation_list_name_1]]$r - correlations_list[[correlation_list_name_2]]$r
  }
  return(corr_diff_matrix)

}






#' Title
#'
#' @param common_reg A comprehensive list of all regions (and all existing subregions) that the brain area in the `rois` list will be compared against.
#' @param rois A list of rois whose regions and subregions will be compared the the `common_reg` list
#'
#' @return common_reg A list of the rois (or any of its subregions) that intersected with the common_reg list.

#'
#' @examples
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


