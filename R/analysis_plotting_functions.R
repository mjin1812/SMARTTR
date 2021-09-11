#' @importFrom ggplot2 ggplot aes theme element_text unit geom_tile geom_text scale_fill_gradient2 labs ggsave guide_legend
NULL

#' @importFrom tidyselect all_of
NULL

#' @importFrom dplyr n mutate
NULL

#' @importFrom tidygraph activate
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
#' @param by (str, default = NULL) Attribute names to group by, e.g. by = c("group", "sex")
#' averaged together to assess across all rois.
#' @return e experiment object with colabel percentage table stored in it.
#' @export
#'
#' @examples
get_percent_colabel <-function(e, by = NULL, channel = "eyfp", save_table = TRUE, rois = "dDG", individual = TRUE){
  # # correct mismatched attributes typed by users
  # by <- match_m_attr(by)

  e_info <- attr(e, "info")

  # Get common regions that are present across both channels
  common_reg <- e$combined_normalized_counts[["colabel"]]$acronym %>% unique() %>%
    intersect( e$combined_normalized_counts[[channel]]$acronym )


  # Check if user only wants a specific roi in common regions
  if (!is.null(rois)){
    if (!all(rois %in% common_reg)){
      message(paste0("The roi(s) specified is not a common region found in both channels. ",
                     "Common regions include:\n"))
      for (reg in common_reg){
        message(reg)
      }
      stop("Set the rois argument to a subset of these acronyms.")
    } else{
      common_reg <- rois
    }
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
  colabel_counts <-  e$combined_normalized_counts$colabel %>% dplyr::filter(acronym %in% common_reg) %>%
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

    # filter by parameters
    df_channel <-  e$combined_normalized_counts[[channel]] %>%
      filter_df_by_char_params(by, values)

    # Pivot wider
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



#' correlation_diff_permutation
#' Note that these correlation lists must have the same number of channels to compare
#'
#' @param e experiment object
#' @param correlation_list_name_1 (default = "AD") The name of the correlation list object used as the first group for comparison.
#' @param correlation_list_name_2 (default = "control") The name of the correlation list object used as the second group for comparison.
#' @param n_shuffle (int, default = 1000) The number of permutation shuffles.
#' @param alpha (float, default = 0.05) The alpha cutoff for significance between region pairwise correlation differences
#' @param p_adjust_method (bool or str, default = "BH") Benjamini-Hochberg method is recommended.
#'  Apply the named method to control for the inflated false discovery rate or FWER. Set to FALSE
#'  to keep "raw" p values.
#' @param seed (int, default = 5) Random seed for future replication.
#' @param ... additional parameters to [RcppAlgos::permuteGeneral()] aside from n, m, Parallel and repetition
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @return e experiment object. The experiment object now has a list called `permutation_p_matrix` stored in it. Elements of this `permutation_p_matrix` list are
#' the outputs of different permutation comparison analyses. These elements are named by the groups that were compared.
#' @export
#' @examples
correlation_diff_permutation <- function(e,
                                         correlation_list_name_1 = "AD",
                                         correlation_list_name_2 = "control",
                                         n_shuffle = 1000,
                                         alpha = 0.05,
                                         p_adjust_method = "BH",
                                         seed = 5,
                                         channels = c("cfos", "eyfp", "colabel"),
                                         ...){


  # Return the correlations list data showing the grouping and the values
  attr_group_1 <- attributes(e$correlation_list[[correlation_list_name_1]])
  attr_group_2 <- attributes(e$correlation_list[[correlation_list_name_2]])

  # # Get overlapping regions between the two correlational datasets
  # channels <- attr_group_1$names

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

    # # Get an array of distribution of correlation differences
    test_statistic_distributions <- permute_corr_diff_distrib(df_channel_groups,
                                                         correlation_list_name_1 = correlation_list_name_1,
                                                         correlation_list_name_2 = correlation_list_name_2,
                                                         n_shuffle = n_shuffle,
                                                         seed = seed, ...)
    # For each pairwise distribution, sort the values
    test_statistic_distributions <- apply(test_statistic_distributions, 1:2, sort) %>% aperm( c(3, 2, 1))

    # calculate the p-value of the permutation
    p_matrix <- matrix(nrow = length(common_regions_btwn_groups),
                       ncol = length(common_regions_btwn_groups),
                       dimnames = dimnames(test_statistic))

    l_reg <- length(common_regions_btwn_groups)
    for (i in 1:l_reg){
      for (j in 1:l_reg){
        null_distrib <- test_statistic_distributions[i,j,]
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
                                     permutation_null_distribution = test_statistic_distributions)
  }
  # Store the results in the experiment object
  comparison <- paste(correlation_list_name_1,"vs",correlation_list_name_2, sep = "_")

  if(is.null(e$permutation_p_matrix)){
    e$permutation_p_matrix <- list()
  }
  e$permutation_p_matrix[[comparison]] <- p_matrix_list
  return(e)
}




#' Create network graph objects for plotting
#' @param e experiment object
#' @param correlation_list_name (str, default = "AD") Name of the correlation list object used to generate the networks.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param edge_color (str, default = "firebrick")  The name of the color or a hexadecimal code used to plot the edges.
#' @param alpha (float, default = 0.05) The significance cutoff for including brain regions in the network.
#'
#' @return e experiment object. This object now has a `networks` list with each element of the list storing network data. The name of
#' each network (`network_name`)is the same as the `correlation_list_name` used to generate the network. This `network_name` is fed into
#' [SMARTR::plot_network()] as a parameter.
#' @export
#'
#' @examples
create_networks <- function(e,
                            correlation_list_name = "AD",
                            channels = c("cfos", "eyfp", "colabel"),
                            edge_color = "firebrick",
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
                                    edges = edges, directed = FALSE)

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
      mutate(avg.dist = rowSums(d, na.rm = TRUE)/(n()-1),
             efficiency = rowSums(d^(-1),na.rm = TRUE)/(n()-1),
             btw = tidygraph::centrality_betweenness(weights = rep(1,nrow(.E())),
                                          directed = FALSE),
             group = as.factor(tidygraph::group_walktrap()))
    networks[[channel]] <- network
  }

  e$networks[[correlation_list_name]] <- networks
  return(e)
}

##_____________________ Plotting functions ___________________________


## THis needs to be cleaned up so that patterns are optional

#' plot_colabel_percent
#'
#' @param e experiment object
#' @param channel
#' @param rois
#' @param color_mapping
#' @param colors
#' @param pattern_mapping
#' @param patterns
#' @param error_bar (str, c("sd", "sem)) options for which type of error bar to display, standard deviation or standard error of the mean.
#' @param by
#'
#' @return
#' @export
#'
#' @examples
plot_percent_colabel <- function(e, by = NULL, channel = "cfos",  rois = c("AAA", "dDG", "HY"),
                                 color_mapping = "sex",
                                 colors = c("#952899", "#358a9c"),
                                 pattern_mapping = "group",
                                 patterns = c("gray100", 'hs_fdiagonal', "hs_horizontal", "gray90", "hs_vertical"),
                                 error_bar = "sem", ylim = c(0, 100),
                                 plot_individual = TRUE,
                                 height = 10, width = 10,
                                 print_plot = TRUE,
                                 save_plot = TRUE,
                                 image_ext = ".png"){

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

  # Check if ggpattern is installed
  if (!requireNamespace("ggpattern", quietly = TRUE)) {
    stop("Package \"ggpattern\" (>= 0.2.0) is needed for this function to work. Please install it now.",
         call. = FALSE)
  }

  # Check that a colabel percentage dataframe exists for this channel
  if(is.null(e$colabel_percent[[channel]])){
    message("The percent colabelled data does not exist in your experiment object for this channel. Autorunning the function now...")
    e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
    e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)
  }

  # Make sure the by (if not NULL) is equal to the by parameter fed to get_percent_colabel() used to generate the dataframes
  attrib <- attributes(e$colabel_percent[[channel]]$average)
  stored_by <- names(attrib$groups) %>% setdiff(c("acronym", ".rows"))

  if (!is.null(by) && !setequal(stored_by, by)){
    loop <- TRUE
    while(loop){
      message("The 'by' parameter used here is not identical to the 'by' parameter that was used in get_percent_colabel().\n",
              "The groups that you intend to compare may be different.")

      inp <- readline("Would you like to rerun the get_percent_colabel() function using the new 'by' parameters for grouping?: Y/N?" )
      if (inp=="Y" || inp=="y") {

        e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
        e <- get_percent_colabel(e, by = by, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)
        loop <- FALSE
      } else if ( inp=="N" || inp == "n") {
        stop("Cannot continue with plotting using the current function parameters.")
      }
    }
  } else{
    by <- stored_by
  }

  # If rois aren't supplied, set as the rois in the dataset
  if (is.null(rois)){
    rois <- unique(e$colabel_percent[[channel]]$average$acronym)
  }

  # Error bar check
  if (error_bar == "sem"){
    error_var <- sym("colabel_percentage.sd")
  } else if (error_bar == "sd") {
    error_var <- sym("colabel_percentage.sem")
  } else {
    stop("You did not supply a valid option for the error_bar. Valid options are 'sem' and 'sd'.")
  }

  # Create  variable names for mapping
  color_var <- sym(color_mapping)
  pattern_var <- sym(pattern_mapping)

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

  # Create the plot
  p <- ggplot(e$colabel_percent[[channel]]$average, aes(x = interaction(!!pattern_var, !!color_var),
                                                        y = colabel_percentage.mean,
                                                        fill = !!color_var)) +
    geom_bar_pattern(aes(pattern_type = !!pattern_var),
                     stat = "identity",
                     color = "black",
                     pattern_color = "black",
                     pattern_fill = "black",
                     pattern = "magick",
                     position = position_dodge(width = .5),
                     show.legend = TRUE,
                     width = 0.4) +
    scale_fill_manual(values = colors) +
    scale_pattern_type_manual(values = patterns) +
    geom_errorbar(aes(ymin = colabel_percentage.mean - !!error_var,
                      ymax = colabel_percentage.mean + !!error_var),
                  width=.2,
                  position = position_dodge(width = .5)) +
    geom_hline(yintercept = 0, element_line(colour = 'black', size=0.5, linetype='solid')) +
    facet_wrap(~acronym,
               strip.position = "bottom") +
    coord_cartesian(ylim = ylim) +
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
    dev.off()
  }
}


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
#' @param store store the plot object in the experiment objects

#'
#' @return e experiment object
#' @export
#'
#' @examples

plot_correlation_heatmaps <- function(e, correlation_list_name = "female_AD", colors = c("#be0000", "#00782e", "#f09b08"),
                                      print_plot = TRUE, save_plot = TRUE,
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
    image_file <- file.path(output_dir, paste0("heatmap_", str_replace(title, " ", "_"), "_", channel, image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
    dev.off()
    }
  }
}




#' Create a Volcano plot.
#'
#' Plot the correlation difference between two comparison groups into a volcano plot. The function
#' [SMARTR::correlation_diff_permutation()] must be run first in order to generate results to plot.
#'
#' @param e
#' @param permutation_comparison The name of the correlation group comparisons to plot.
#' @param title Title of the plot.
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save the correlation heatmap plots into the figures subdirectory of the
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) channels to plot
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal codes corresponding to the channels (respectively) to plot.
#' @param alpha (float, default = 0.05) The alpha value used to plot horizontal line. Point above the line are considered significant.
#' @param store (bool, default = FALSE) Store the plot object into your experiment.
#'
#' @return e experiment object
#' @export
#' @examples
volcano_plot <- function(e,
                         permutation_comparison = "female_AD_vs_male_AD",
                         channels = c("cfos", "eyfp", "colabel"),
                         colors =  c("#be0000", "#00782e", "#f09b08"),
                         alpha = 0.05,
                         save_plot = TRUE,
                         title = NULL,
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

  for (k in 1:length(channels)){
    # get the p-values, significance, and corr_diff
    p_vals <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$p_val %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    corr_diffs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$test_statistic %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")
    sigs <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$sig %>%
      tibble::as_tibble(rownames = NA) %>% tibble::rownames_to_column(var = "rowreg")


    # Pivot long
    p_vals <- p_vals %>% pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                      values_to = "p_val", names_to = "colreg")
    corr_diffs <- corr_diffs %>% pivot_longer(col = - "rowreg", values_drop_na = TRUE,
                                              values_to = "corr_diff", names_to = "colreg")
    sigs <- sigs %>% pivot_longer(col = - "rowreg", values_drop_na = TRUE,
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
      ylim(c(0, 3)) +
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
      dev.off()
    }
  }
}





#' Plot the networks stored in an experiment object
#'
#' @param e experiment object
#' @param network_name (str, default = "AD")
#' @param title
#' @param channels
#' @param height
#' @param width
#' @param image_ext
#' @param save_plot
#'
#' @return
#' @export
#'
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

  for (channel in channels){

    network <- e$networks[[network_name]][[channel]]

    # _______________ Plot the network ________________________________
    theme.network <- theme_graph() + theme(plot.title = element_text(hjust = 0.5,size = 28),
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
    p <-  p + ggtitle(title)

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
      image_file <- file.path(output_dir, paste0("network_", network_name, "_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
      dev.off()
    }
  }

}




# plot_permutation_histogram <- function(e, roi = "dDG"){

#
# }

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






