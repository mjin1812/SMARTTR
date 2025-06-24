

##_____________________ Plotting functions ___________________________

#' This function allows for plotting of colabelled cells over either the "cfos" or "eyfp" channels.
#'
#' Allows for specification of specific brain regions to plot. Two different mouse attributes can be used as categorical variables to map to either the color or
#' pattern aesthetics of the bar plot, e.g. sex and experimental group.
#' The color aesthetic takes precedence over the pattern aesthetic so if you only want to use one mouse attribute, for plotting
#' set it to the `color_mapping` parameter and set the `pattern_mapping` parameter to NULL.
#'
#' @param e experiment object
#' @param colabel_channel (str, default = "colabel") The channel used as the numerator in fraction counts.
#' @param channel (str, default = "eyfp") The channel used as denominator in fraction counts.
#' @param rois character vector of region acronyms, e.g. c("AAA", "DG)
#' @param color_mapping (str, default = "sex") The name of the categorical variable (e.g., "sex", "age", etc.) to map to the color aesthetic of the bar plot.
#' @param colors (str) character vector of the color values desired for the groups.
#' @param pattern_mapping (str, default = "sex") The name of the categorical variable (e.g., "sex", "age", etc.) to map to the pattern aesthetic of the bar plot.
#' @param patterns (str, default = c("gray100", 'hs_fdiagonal', "hs_horizontal", "gray90", "hs_vertical"), Pattern types to define subgroups.
#' @param error_bar (str, c("sd", "sem)) options for which type of error bar to display, standard deviation or standard error of the mean.
#' @param ylim  (default = c(0,100)) The range of the y-axis
#' @param plot_individual (boo) Whether or not to plot multiple
#' @param height (default = 8) height of graphics devices in inches
#' @param width (default = 8) height of graphics device in inches
#' @param print_plot (bool, default = FALSE) whether or not to print the plot for just for display
#' @param save_plot (bool, default = TRUE) whether or not to save the plor
#' @param image_ext (default = ".png") extension determining the image type to save as
#' @return p Plot handle to the figure
#' @export
#' @examples
#' \dontrun{
#'  plot_percentage_colabel
#' }
plot_percent_colabel <- function(e,
                                 colabel_channel = "colabel",
                                 channel = "eyfp",
                                 rois = c("AAA", "dDG", "HY"),
                                 color_mapping = "sex",
                                 colors = c("#952899", "#358a9c"),
                                 pattern_mapping = NULL,
                                 patterns = c("gray100", 'hs_fdiagonal', "hs_horizontal
                                              ", "gray90", "hs_vertical"),
                                 error_bar = "sem",
                                 ylim = c(0, 100),
                                 plot_individual = TRUE,
                                 height = 8,
                                 width = 8,
                                 print_plot = FALSE,
                                 save_plot = TRUE,
                                 image_ext = ".png"){

  if(is.null(pattern_mapping)){
    by <- match_m_attr(color_mapping)
  } else{
    by <- match_m_attr(c(color_mapping, pattern_mapping))
  }

  e <- get_percent_colabel(e, by = by, colabel_channel = colabel_channel, channel = channel, save_table = FALSE, rois = rois, individual = FALSE)
  e <- get_percent_colabel(e, by = by, colabel_channel = colabel_channel, channel = channel, save_table = FALSE, rois = rois, individual = TRUE)

  child_r <- get.acronym.child(rois)
  while (length(child_r) > 0){
    rois <- c(rois, child_r)
    child_r <- get.acronym.child(child_r) %>% na.omit()
  }

  if (error_bar == "sem"){
    error_var <- rlang::sym("colabel_percentage.sem")
  } else if (error_bar == "sd") {
    error_var <-rlang::sym("colabel_percentage.sd")
  } else {
    stop("You did not supply a valid option for the error_bar. Valid options are 'sem' and 'sd'.")
  }

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
                     axis.line.y = element_line(colour = 'black', linewidth = 0.5, linetype ='solid'))

  color_var <- rlang::sym(color_mapping)

  if (!is.null(pattern_mapping)){
    if (!requireNamespace("ggpattern", quietly = TRUE)) {
      stop(
        "Package \"ggpattern\" (>= 0.2.0) is needed for the pattern mapping parameter to work. Please install it or set `pattern_mapping` to NULL.",
        call. = FALSE)
    }
    pattern_var <- rlang::sym(pattern_mapping)

    # Create the plot
    p <- ggplot(e$colabel_percent[[channel]]$average %>% dplyr::filter(.data$acronym %in% rois),
                aes(x = interaction(!!pattern_var, !!color_var),
                    y = .data$colabel_percentage.mean,
                    fill = !!color_var)) +
      ggpattern::geom_bar_pattern(aes(pattern_type = !!pattern_var),
                                  stat = "identity",
                                  color = "black",
                                  pattern_color = "black",
                                  pattern_fill = "black",
                                  pattern = "magick",
                                  position = ggplot2::position_dodge(width = .5),
                                  show.legend = TRUE,
                                  width = 0.4) +
      ggplot2::scale_fill_manual(values = colors) +
      ggpattern::scale_pattern_type_manual(values = patterns) +
      ggplot2::geom_errorbar(aes(ymin = .data$colabel_percentage.mean - !!error_var,
                                 ymax = .data$colabel_percentage.mean + !!error_var),
                             width=.2,
                             position = ggplot2::position_dodge(width = .5)) +
      ggplot2::geom_hline(yintercept = 0, element_line(colour = 'black', linewidth = 0.5, linetype='solid')) +
      ggplot2::facet_wrap(~.data$acronym,
                          strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = paste0("co-labeled / ", channel, "+ cells (%)")) +
      bar_theme

    if (plot_individual){
      df <- e$colabel_percent[[channel]]$individual %>% dplyr::filter(.data$acronym %in% rois)
      p <- p + geom_jitter(data = df,
                           aes(x = interaction(!!pattern_var, !!color_var),
                               y = .data$colabel_percentage),
                           size = 2,
                           width = 0.1)
    }
  } else{
    p <- ggplot(e$colabel_percent[[channel]]$average %>% dplyr::filter(.data$acronym %in% rois),
                aes(x = !!color_var,
                    y = .data$colabel_percentage.mean,
                    fill = !!color_var)) +
      ggplot2::geom_bar(stat='identity',
                        position = ggplot2::position_dodge(width = .5),
                        color = "black",
                        show.legend = TRUE,
                        width = 0.4) +
      scale_fill_manual(values = colors) +
      geom_errorbar(aes(ymin = .data$colabel_percentage.mean - !!error_var,
                        ymax = .data$colabel_percentage.mean + !!error_var),
                    width=.2,
                    position = ggplot2::position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', linewidth = 0.5, linetype='solid')) +
      ggplot2::facet_wrap(~.data$acronym,
                          strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = paste0("co-labeled / ", channel, "+ cells (%)")) +
      bar_theme

    if (plot_individual){
      df <- e$colabel_percent[[channel]]$individual %>% dplyr::filter(.data$acronym %in% rois)
      p <- p + geom_jitter(data = df,
                           aes(x = !!color_var,
                               y = .data$colabel_percentage),
                           size = 2,
                           width = 0.1)
    }
  }
  if (print_plot){
    dev.new(noRStudioGD=TRUE)
    print(p)
  }

  if(save_plot){
    dev.new(width = width, height = height, noRStudioGD=TRUE)
    print(p)
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


#' This function allows for plotting of normalized cell counts by area across specific regions to plot.
#'
#' Two different mouse attributes can be used as categorical variables to map to either the color or
#' pattern aesthetics of the bar plot, e.g. sex and experimental group.
#' The color aesthetic takes precedence over the pattern aesthetic so if you only want to use one mouse attribute, for plotting
#' set it to the `color_mapping` parameter and set the `pattern_mapping` parameter to NULL.
#'
#' Set to be deprecated in favor of plot_normalized_counts()
#'
#' @param e experiment object
#' @param channel (str, default = "eyfp") The channel used as denominator in fraction counts.
#' @param rois (vec,  default = c("AAA", "dDG", "HY")) Allen acronyms of the ROIS that the user would like to plot
#' @param color_mapping (str, default = "group") The variable name that maps subgroups  you would like to graphically distinguish through colors.
#' @param colors (vec, default = c("#952899", "#358a9c")) A vector of hexadecimal color codes for each subgroup distinguished by the color mapping variable.
#' @param pattern_mapping (str, default = NULL) variable name that maps subgroups  you would like to graphically distinguish through bar patterns. Set to NULL if not in use. Should not be the same mapping as color_mapping.
#' @param patterns  (default = c("gray100", 'hs_fdiagonal', "hs_horizontal", "gray90", "hs_vertical") Available patterns in ggpattern package to map to subgroups distinguished by the pattern mapping variable.
#' @param ylab (str, default =  bquote('Cell counts '('cells/mm'^3))) unit of measurement
#' @param ylim (vec, default = c(0,100)) Y-axis range.
#' @param plot_individual (bool) whether to plot individual points
#' @param height (numeric) height of the plot in inches to save as.
#' @param width (numeric) height of the plot in inches to save as.
#' @param print_plot (bool) Whether to to print the plot.
#' @param save_plot (bool) Whether to save the plot in the experiment figures folder.
#' @param image_ext (str, default = ".png") Extension of the output file
#' @param error_bar (str, c("sd", "sem)) options for which type of error bar to display, standard deviation or standard error of the mean.
#' @return p Plot handle to the figure
#' @noRd
#' @examples
#' \dontrun{
#' # To add example
#' }
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
                             print_plot = FALSE,
                             save_plot = TRUE,
                             image_ext = ".png"){

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

  child_r <- get.acronym.child(rois)
  while (length(child_r) > 0){
    rois <- c(rois, child_r)
    child_r <- get.acronym.child(child_r) %>% na.omit()
  }


  df <- e$combined_normalized_counts[[channel]] %>% dplyr::group_by_at(c("group", "acronym", "name")) %>%
    dplyr::summarise(normalized.count.mean = mean(.data$normalized.count.by.volume),
                     normalized.count.sem = sd(.data$normalized.count.by.volume)/sqrt(n()),
                     normalized.count.sd = sd(.data$normalized.count.by.volume),
                     n = n())

  df <- df %>% dplyr::filter(.data$acronym %in% rois)
  df_indiv <- e$combined_normalized_counts[[channel]]  %>% dplyr::filter(.data$acronym %in% rois)

  if (error_bar == "sem"){
    error_var <- rlang::sym("normalized.count.sem")
  } else if (error_bar == "sd") {
    error_var <- rlang::sym("normalized.count.sem")
  } else {
    stop("You did not supply a valid option for the error_bar. Valid options are 'sem' and 'sd'.")
  }

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
                     axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

  color_var <- rlang::sym(color_mapping)
  if (!is.null(pattern_mapping)){
    if (!requireNamespace("ggpattern", quietly = TRUE)) {
      message("Package \"ggpattern\" (>= 0.2.0) is needed for the pattern mapping to work. Installing it now.",
              call. = FALSE)
      utils::install.packages('ggpattern')
    }
    pattern_var <- rlang::sym(pattern_mapping)
    p <- ggplot(df,
                aes(x = interaction(!!pattern_var, !!color_var),
                    y = .data$normalized.count.mean,
                    fill = !!color_var)) +
      ggpattern::geom_bar_pattern(aes(pattern_type = !!pattern_var),
                                  stat = "identity",
                                  color = "black",
                                  pattern_color = "black",
                                  pattern_fill = "black",
                                  pattern = "magick",
                                  position = ggplot2::position_dodge(width = .5),
                                  show.legend = TRUE,
                                  width = 0.4) +
      ggplot2::scale_fill_manual(values = colors) +
      ggpattern::scale_pattern_type_manual(values = patterns) +
      geom_errorbar(aes(ymin = .data$normalized.count.mean - !!error_var,
                        ymax = .data$normalized.count.mean + !!error_var),
                    width=.2,
                    position = ggplot2::position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', linewidth = 0.5, linetype='solid')) +
      ggplot2::facet_wrap(~.data$acronym,
                          strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = ylab) +
      bar_theme

    if (plot_individual){
      p <- p + geom_jitter(data = df_indiv,
                           aes(x = interaction(!!pattern_var, !!color_var),
                               y = .data$count),
                           size = 2,
                           width = 0.1)
    }
  } else{
    p <- ggplot(df,
                aes(x = !!color_var,
                    y = .data$normalized.count.mean,
                    fill = !!color_var)) +
      ggplot2::geom_bar(stat='identity',
                        position = ggplot2::position_dodge(width = .5),
                        color = "black",
                        show.legend = TRUE,
                        width = 0.4) +
      scale_fill_manual(values = colors) +
      geom_errorbar(aes(ymin = .data$normalized.count.mean - !!error_var,
                        ymax = .data$normalized.count.mean + !!error_var),
                    width=.2,
                    position = ggplot2::position_dodge(width = .5)) +
      geom_hline(yintercept = 0, element_line(colour = 'black', linewidth = 0.5, linetype='solid')) +
      ggplot2::facet_wrap(~.data$acronym,
                          strip.position = "bottom") +
      ylim(ylim) +
      geom_text(aes(x=c(2.5), y= 0.4, label=c("|")),
                vjust=1.2, size=3) +
      labs(y = ylab) +
      bar_theme
    if (plot_individual){
      p <- p + geom_jitter(data = df_indiv,
                           aes(x = !!color_var,
                               y = .data$normalized.count.by.volume),
                           size = 2,
                           width = 0.1)
    }
  }

  if (print_plot){
    dev.new(noRStudioGD=TRUE)
    print(p)
  }

  if(save_plot){
    # Plot the plot
    dev.new(width = width, height = height, noRStudioGD=TRUE)
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
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param title (str, default = NULL) An optional title for the plot
#' @param unit_label (str, default = bquote('Cell counts '('cells/mm'^3))) Default unit label for the graphs
#' @param anatomical.order (default = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU","TH", "HY", "MB", "HB", "CB")) Default way to group subregions into super regions order
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param print_plot (bool, default = FALSE) Whether to display the plot (in addition to saving the plot)
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of
#'  the experiment object output folder.
#' @param flip_axis plot cell counts on x-axis rather than y-axis.
#' @param reverse_colors (bool, default = FALSE) Whether to reverse the color order. This may depend on the order in which you entered the `colors` parameter
#' @param limits (c(0,100000)) Range of the normalized cell counts.
#' @param facet_background_color (default = NULL) Set to a hexadecimal string, e.g."#FFFFFF", when you want to shade the background of the graph. Defaults to no background when NULL.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param strip_background_colors (default = "lightblue) Enter custom codes to control the strip background colors, e.g. c(Isocortex = "#5571a9", OLF = "#64bdc4",
#' HPF = "#d2875b", CTXsp = "#87a3db", CNU = "#466496", TH = "#7e72af", HY = "#8e7960",  MB = "#d796c8", HB = "#646464"). If more than one color is used, you must install the package ggh4x.
#' @param plot_theme (ggplot2 theme object) This allows for fine tuning the aesthetics of the figure. Default parameters shown:
#' ggplot2::theme(plot.background = element_blank(),
#'              panel.grid.major = element_blank(),
#'               panel.grid.minor = element_blank(),
#'               panel.border = element_blank(),
#'               axis.line = element_line(color = 'black'),
#'               legend.justification = c(0, 0),
#'               legend.position = "inside",
#'               legend.position.inside = c(0.05, 0.6),
#'               legend.direction = "vertical",
#'               axis.text.y = element_text(angle = 50,
#'                                          hjust = 1,
#'                                          size = 8,
#'                                          color = "black"),
#'               axis.text.x = element_text(color = "black"),
#'               strip.text.y = element_text(angle = 0,
#'                                           margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")),
#'               strip.placement = "outside",
#'               strip.switch.pad.grid = unit(0.1, "in"))
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p_list <- plot_normalized_counts(e, channels = "cfos", by = c("sex", "group"),
#' values = list(c("female", "non"), c("female", "agg")), colors = c("white", "lightblue"))
#' }
plot_normalized_counts <- function(e,
                                   channels = c("cfos", "eyfp", "colabel"),
                                   by = c("sex", "group"),
                                   values = list(c("female", "non"),
                                                 c("female", "agg"),
                                                 c("female", "control"),
                                                 c("male", "agg"),
                                                 c("male", "control")),
                                   colors = c("white", "lightblue", "black", "red", "green"),
                                   ontology = "allen",
                                   title = NULL,
                                   unit_label = bquote('Cell counts '('cells/mm'^3)),
                                   anatomical.order = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU",
                                                        "TH", "HY", "MB", "HB", "CB"),
                                   height = 7,
                                   width = 20,
                                   print_plot = FALSE,
                                   save_plot = TRUE,
                                   flip_axis = FALSE,
                                   reverse_colors = FALSE,
                                   limits = c(0, 100000),
                                   facet_background_color =  NULL,
                                   strip_background_colors =  "lightblue",
                                   plot_theme = ggplot2::theme(plot.background = element_blank(),
                                                               panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank(),
                                                               panel.border = element_blank(),
                                                               axis.line = element_line(color = 'black'),
                                                               legend.justification = c(0, 0),
                                                               legend.position = "inside",
                                                               legend.position.inside = c(0.05, 0.6),
                                                               legend.direction = "vertical",
                                                               axis.text.y = element_text(angle = 50,
                                                                                          hjust = 1,
                                                                                          size = 8,
                                                                                          color = "black"),
                                                               axis.text.x = element_text(color = "black"),
                                                               strip.text.y = element_text(angle = 0,
                                                                                           margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt")),
                                                               strip.placement = "outside",
                                                               strip.switch.pad.grid = unit(0.1, "in")),
                                   image_ext = ".pdf"){
  p_list <- vector(mode = 'list', length = length(channels))
  labels <- purrr::map(values, function(x){paste(x, collapse = "_")}) %>% unlist()

  if (length(strip_background_colors) > 1){
    if (!requireNamespace("ggh4x", quietly = TRUE)) {
      message("Package \"ggh4x\" is needed to plot more than 1 strip background color. Installing it now.",
              call. = FALSE)
      utils::install.packages('ggh4x')
    }

    background_x <- vector(mode = "list", length = length(strip_background_colors))
    for (c in 1:length(strip_background_colors)){
      background_x[[c]] <- element_rect(fill = strip_background_colors[[c]], color = "black")
    }
  }

  names(p_list) <- channels
  for (k in 1:length(channels)) {
    channel <- channels[[k]]
    channel_counts <- NULL
    for (g in 1:length(values)){
      if (is.null(channel_counts)){
        channel_counts <-  e$combined_normalized_counts[[channel]] %>% dplyr::ungroup() %>%
          filter_df_by_char_params(by, values[[g]]) %>% dplyr::distinct() %>%
          dplyr::select(dplyr::all_of(c(by, "mouse_ID", "name", "acronym", "normalized.count.by.volume"))) %>%
          dplyr::group_by(across(all_of(c(by, "acronym", "name" )))) %>%
          dplyr::summarise(mean_normalized_counts = mean(.data$normalized.count.by.volume),
                           sem = sd(.data$normalized.count.by.volume, na.rm=TRUE)/sqrt(n()))
      } else {
        to_bind <-  e$combined_normalized_counts[[channel]] %>% dplyr::ungroup() %>%
          filter_df_by_char_params(by, values[[g]]) %>% dplyr::distinct() %>%
          dplyr::select(dplyr::all_of(c(by, "mouse_ID", "name", "acronym", "normalized.count.by.volume"))) %>%
          dplyr::group_by(across(all_of(c(by, "acronym", "name" )))) %>%
          dplyr::summarise(mean_normalized_counts = mean(.data$normalized.count.by.volume),
                           sem = sd(.data$normalized.count.by.volume, na.rm=TRUE)/sqrt(n()))
        channel_counts <- dplyr::bind_rows(channel_counts, to_bind)
      }
    }
    if (tolower(ontology) == "allen") {
      regions.ordered <- anatomical.order %>%
        purrr::map(get.sub.structure) %>%
        unlist()
    } else {
      regions.ordered <- anatomical.order %>%
        purrr::map(get.sub.structure.custom, ontology = ontology) %>%
        unlist()
    }
    common.regions.ordered <- channel_counts$name[unique(match(regions.ordered, channel_counts$acronym))]
    channel_counts$name <- factor(channel_counts$name, levels = common.regions.ordered)
    if (tolower(ontology) == "allen") {
      channel_counts <- channel_counts %>%
        mutate(parent = get.sup.structure(.data$acronym, matching.string = anatomical.order)) %>%
        na.omit()
      channel_counts$parent <- factor(channel_counts$parent, levels = anatomical.order)
    } else {
      channel_counts <- channel_counts %>%
        mutate(parent = get.sup.structure.custom(.data$acronym, ontology = ontology, matching.string = anatomical.order)) %>%
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

    if (reverse_colors){
      channel_counts$unique_groups <- factor(channel_counts$unique_groups, levels = rev(labels))
    } else{
      channel_counts$unique_groups <- factor(channel_counts$unique_groups, levels = labels)
    }

    if (flip_axis) {
      p <- channel_counts %>%
        ggplot(aes(x = .data$mean_normalized_counts, y = .data$name,
                   fill = .data$unique_groups), color = "black") +
        geom_col(position = ggplot2::position_dodge(0.8), width = 0.8, color = "black") +
        geom_errorbar(aes(xmin = .data$mean_normalized_counts - .data$sem,
                          xmax = .data$mean_normalized_counts + .data$sem,
                          y = .data$name),
                      position = ggplot2::position_dodge(0.8),
                      width = 0.5,
                      color = "black") +
        labs(title = title,
             x = unit_label,
             y = "",
             fill = "Group") +
        ggplot2::scale_x_continuous(expand = c(0,0), limits = limits) +
        ggplot2::scale_fill_manual(values=colors)
      if (length(strip_background_colors) > 1) {
        p <- p +  ggh4x::facet_grid2(.data$parent~., scales = "free", space = "free_y", switch = "y",
                                     strip = ggh4x::strip_themed(background_y = background_x)) + plot_theme
      } else {
        p <- p + theme(strip.background = element_rect(color = "black",
                                                       fill =  strip_background_colors)) + plot_theme
      }
    } else if (isFALSE(flip_axis)) {
      p <- channel_counts %>%
        ggplot(aes(y = .data$mean_normalized_counts, x = .data$name,
                   fill = .data$unique_groups), color = "black") +
        geom_col(position = ggplot2::position_dodge(0.8), width = 0.8, color = "black") +
        geom_errorbar(aes(ymin = .data$mean_normalized_counts - .data$sem,
                          ymax = .data$mean_normalized_counts + .data$sem, x = .data$name),
                      position = ggplot2::position_dodge(0.8),
                      width = 0.5,
                      color = "black") +
        labs(title = title,
             y = unit_label,
             x = "",
             fill = "Group") +
        scale_y_continuous(expand = c(0,0), limits = limits) +
        scale_fill_manual(values=colors, labels = labels)
      if (length(strip_background_colors) > 1) {
        p <- p +  ggh4x::facet_grid2(.data$parent~., scales = "free", space = "free_y", switch = "y",
                                     strip = ggh4x::strip_themed(background_y = background_x)) + plot_theme
      } else {
        p <- p + theme(strip.background = element_rect(color = "black",
                                                       fill =  strip_background_colors)) + plot_theme
      }
    }
    if(!is.null(facet_background_color)){
      p <- p +
        theme(panel.background = element_rect(fill = facet_background_color, color = "black"))
    }
    if (print_plot) {
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if (save_plot) {
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      ggsave(p, filename = paste0(channels[k], "_normalized_counts", image_ext),
             path = file.path(attr(e, 'info')$output_path, "figures"), width = width, height = height, limitsize = FALSE)
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

#' Plot correlation heatmaps
#'
#' @param e experiment object. Must contain a named correlation_list object generated by [SMARTTR::get_correlations()]
#' @param correlation_list_name (str) The name of the correlation object generated by [SMARTTR::get_correlations()]
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Must exist in the channels attribute of the correlation_list.
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal code for the colors corresponding channels parameter. Color values can also be
#' input compatible with ggplot2 plotting functions. If value is NULL scale_fill_viridis will be used.
#' @param print_plot (bool, default = FALSE) Print the plot as graphics windows.
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
#' @param strip_background_colors (default = "lightblue) Enter custom codes to control the strip background colors, e.g. c(Isocortex = "#5571a9", OLF = "#64bdc4",
#' HPF = "#d2875b", CTXsp = "#87a3db", CNU = "#466496", TH = "#7e72af", HY = "#8e7960",  MB = "#d796c8", HB = "#646464", CB = "#dddd7f").
#' If more than one color is used, you must install the package ggh4x. If a single color is used, it will
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' plot_correlation_heatmaps(e, correlation_list_name = "female_AD")
#' }
#' @seealso [SMARTTR::get_correlations()]

plot_correlation_heatmaps <- function(e,
                                      correlation_list_name,
                                      channels = c('cfos', 'eyfp', 'colabel'),
                                      colors = c("#be0000", "#00782e", "#f09b08"),
                                      sig_color = "yellow",
                                      sig_nudge_y = -0.7,
                                      sig_size = 7,
                                      ontology = "allen",
                                      anatomical.order = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU",
                                                           "TH", "HY", "MB", "HB", "CB"),
                                      print_plot = FALSE,
                                      save_plot = TRUE,
                                      image_ext = ".png",
                                      plot_title = NULL,
                                      height = 10,
                                      width = 10,
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
                                                                strip.placement = "outside"),
                                      # strip_background_colors = c(Isocortex = "#5571a9", OLF = "#64bdc4", HPF = "#d2875b",
                                      #                             CTXsp = "#87a3db", CNU = "#466496", TH = "#7e72af",
                                      #                             HY = "#8e7960",  MB = "#d796c8", HB = "#646464", CB = "#dddd7f")
                                      strip_background_colors = "lightblue"
){
  # Get the attributes for plotting
  if (is.null(e$correlation_list[[correlation_list_name]])){
    stop(paste0("Your experiment object doesn't contain the specified correlation_list object to plot! ",
                "\nRun the function, get_correlations() first."))
  }

  cl_attr <- attributes(e$correlation_list[[correlation_list_name]])
  if (!is.null(colors)){
    names(colors) <- channels
  }

  if (is.null(plot_title)){
    cl_attr$values <- stringr::str_to_title(cl_attr$values)
    plot_title <- paste(cl_attr$value, collapse = " ")
  }

  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (channel in channels){
    corr_df <-  e$correlation_list[[correlation_list_name]][[channel]] %>% purrr::map(tibble::as_tibble)
    val_names <- names(corr_df)

    for (k in 1:length(corr_df)){
      corr_df[[k]] <- tibble::add_column(corr_df[[k]], row_acronym = names(corr_df[[k]]), .before = TRUE ) %>%
        tidyr::pivot_longer(!any_of("row_acronym"), names_to = "col_acronym", values_to = val_names[k])
      corr_df[[k]]$row_acronym <- factor(corr_df[[k]]$row_acronym, levels = unique(corr_df[[k]]$row_acronym))
      corr_df[[k]]$col_acronym <- factor(corr_df[[k]]$col_acronym, levels = unique(corr_df[[k]]$col_acronym))
    }

    df <- corr_df$r %>% dplyr::left_join(corr_df$n, by = c("row_acronym", "col_acronym")) %>%
      dplyr::left_join(corr_df$P, by = c("row_acronym", "col_acronym")) %>%
      dplyr::left_join(corr_df$sig, by = c("row_acronym", "col_acronym")) %>%
      dplyr::mutate(sig_text = dplyr::if_else(.data$sig == TRUE, "*", ""))

    df <- df %>% dplyr::mutate(row_parent = get.super.regions(.data$row_acronym, anatomical.order = anatomical.order,
                                                              ontology = ontology), .before  = any_of("row_acronym")) %>%
      dplyr::mutate(row_parent = factor(.data$row_parent, levels = anatomical.order))
    df <- df %>% dplyr::mutate(col_parent = get.super.regions(.data$col_acronym, anatomical.order = anatomical.order,
                                                              ontology = ontology), .before  = any_of("row_acronym")) %>%
      dplyr::mutate(col_parent = factor(.data$col_parent, levels = rev(anatomical.order)))
    p <-  ggplot(df, aes(.data$row_acronym, .data$col_acronym, fill = .data$r))
    if (length(strip_background_colors) > 1) {
      if (!requireNamespace("ggh4x", quietly = TRUE)) {
        message("Package \"ggh4x\" is needed to plot more than 1 strip background color. Installing it now.",
                call. = FALSE)
        utils::install.packages('ggh4x')
      }
      background_x <- vector(mode = "list", length = length(strip_background_colors))
      for (c in 1:length(strip_background_colors)){
        background_x[[c]] <- element_rect(fill = strip_background_colors[[c]], color = "black")
      }
      p <-  p +
        ggh4x::facet_grid2(col_parent ~ row_parent,
                           space = "free",
                           margins = FALSE,
                           scales = "free",
                           strip = ggh4x::strip_themed(background_y = rev(background_x),
                                                       background_x = background_x)) +
        geom_tile() +
        geom_text(aes(label =  .data$sig_text), size=sig_size, position = position_nudge(y = sig_nudge_y), color = sig_color) +
        labs(title = plot_title, x = "Brain Region", y = "Brain Region") +
        theme.hm
    } else {
      p <- p +
        ggplot2::facet_grid(.data$col_parent ~ .data$row_parent,
                            space = "free",
                            margins = FALSE,
                            scales = "free") +
        geom_tile() +
        geom_text(aes(label = .data$sig_text), size=sig_size, position = position_nudge(y = sig_nudge_y), color = sig_color) +
        labs(title = plot_title, x = "Brain Region", y = "Brain Region") +
        ggplot2::theme(strip.background = element_rect(color = "black", fill = strip_background_colors)) +
        theme.hm
    }
    if (is.null(colors)){
      p <- p + scale_fill_viridis_c(aesthetics = c("color","fill"), option="magma", na.value = "grey50",
                                    limits=c(-1, 1))
    } else {
      p <- p + scale_fill_gradient2(low = "#4f4f4f",mid = "#ffffff", high = colors[[channel]],
                                    aesthetics = c("color","fill"), na.value = "grey50",
                                    limits=c(-1, 1))
    }

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("heatmap_", stringr::str_replace(plot_title, " ", "_"), "_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channel]] <- p
  }
  return(p_list)
}


#' Plot correlation heatmaps for a dataset with hemispheric information
#'
#' @param e experiment object. Must contain a named correlation_list object generated by [SMARTTR::get_correlations()]
#' @param correlation_list_name (str) The name of the correlation object generated by [SMARTTR::get_correlations()]
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Must exist in the channels attribute of the correlation_list.
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal code for the colors corresponding channels parameter. Color values can also be
#' input compatible with ggplot2 plotting functions.
#' @param plot_title (str, default = NULL) If NULL, the `correlation_list_name` will used as the title with underscores removed
#' @param sig_color (str, default = "yellow") Color of the significance symbol in R
#' @param sig_size (default = 7) Point size for significance symbol.
#' @param sig_nudge_y (default = -0.7) Relative amount to nudge the significance symbols in the y direction to center over each square.
#' @param anatomical.order (default = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU","TH", "HY", "MB", "HB", "CB")) Default way to group subregions into super regions order
#' @param strip_background_colors Named vector of hexadecimal strings corresponding to base parent regions.
#' @param height (int) Height of the plot in inches.
#' @param width (int) Width of the plot in inches.
#' @param print_plot (bool, default = FALSE) Print the plot as graphics windows.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#' the experiment object output folder.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param theme.hm Option to use custom ggplot2 theme if the user wants. See default values as example.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
plot_correlation_heatmaps_hemi <- function(e,
                                           correlation_list_name,
                                           channels = c('rcbv'),
                                           colors = NULL, # if NULL, scale_fill_viridis is used
                                           sig_color = "yellow",
                                           sig_nudge_y = -0.7,
                                           sig_size = 7,
                                           ontology = "allen",
                                           anatomical.order = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU",
                                                                "TH", "HY", "MB", "HB", "CB"),
                                           print_plot = TRUE, save_plot = TRUE, image_ext = ".png",
                                           plot_title = NULL, height = 10, width = 10,
                                           strip_background_colors = c(Isocortex = "#70ff71", OLF = "#75aadb", HPF = "#7ed04b",
                                                                       CTXsp = "#8ada87", CNU = "#98d6f9", TH = "#ff7080",
                                                                       HY = "#e64438",  MB = "#ff64ff", HB = "#ff9b88"),

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
                                                                     strip.placement = "outside"
                                           )
){


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
    df <- df %>% dplyr::mutate(row_parent = get.super.regions({function (x)  return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}(row_acronym),
                                                              anatomical.order = anatomical.order, ontology = ontology), .before  = row_acronym) %>%
      dplyr::mutate(row_parent = factor(row_parent, levels = anatomical.order))
    df <- df %>% dplyr::mutate(col_parent = get.super.regions({function (x)  return(str_remove_all(x, "_L") %>% str_remove_all("_R"))}(col_acronym),
                                                              anatomical.order = anatomical.order, ontology = ontology), .before  = row_acronym) %>%
      dplyr::mutate(col_parent = factor(col_parent, levels = rev(anatomical.order)))
    # df <- df %>% dplyr::mutate(sig_text = ifelse(is.finite(sig_text), sig_text, ""))


    # Arrange the custom strip colors
    background_x <- vector(mode = "list", length = length(strip_background_colors))
    for (c in 1:length(strip_background_colors)){
      background_x[[c]] <- element_rect(fill = strip_background_colors[[c]], color = "black")
    }

    # Generate a correlation heatmap in anatomical order
    n_facet <- df$row_parent %>% unique() %>% length()
    p <-  ggplot(df, aes(row_acronym, col_acronym, fill = r)) +
      ggh4x::facet_grid2(col_parent ~ row_parent,
                         space = "free",
                         margins = FALSE,
                         scales = "free",
                         strip = ggh4x::strip_themed(background_y = background_x,
                                                     background_x = background_x)) +
      geom_tile() +
      geom_text(aes(label = sig_text), size=sig_size, position = position_nudge(y = sig_nudge_y), color = sig_color) +
      labs(title = plot_title, x = "Brain Region", y = "Brain Region") +
      theme.hm

    if (is.na(colors)){
      p <- p + scale_fill_viridis_c(aesthetics = c("color","fill"), option="magma", na.value = "grey50",
                                    limits=c(-1, 1))

    } else {
      p <- p + scale_fill_gradient2(low = "#4f4f4f",mid = "#ffffff", high = colors[[channel]],
                                    aesthetics = c("color","fill"), na.value = "grey50",
                                    limits=c(-1, 1))
    }

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      # Plot the heatmap
      dev.new(width = width, height = height, noRStudioGD=TRUE)
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
#' [SMARTTR::correlation_diff_permutation()] must be run first in order to generate results to plot.
#'
#' @param e experiment object
#' @param permutation_comparison The name of the correlation group comparisons to plot.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) channels to plot.
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal code for the colors corresponding to the
#' channels attribute of the correlation_list. Color values can also be input compatible with ggplot2 plotting functions.
#' @param print_plot (bool, default = FALSE) Print the plot as graphics windows.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#' the experiment object output folder.
#' @param image_ext (default = ".png") image extension to save the plot as.
#' @param title Title of the plot.
#' @param height height of the plot in inches.
#' @param ylim (vec, default = c(0,y)) Y-axis range (logarithmic).
#' @param plt_theme (default = NULL) Add a [ggplot2::theme()] to the plot. If NULL, the default is taken..
#' @param point_size (default = 1) Size of the plotted points.
#' @param width width of the plot in inches.
#'
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' volcano_plot(e, permutation_comparison = "female_AD_vs_male_AD",
#' channels = c("cfos", "eyfp", "colabel"), colors =  c("#be0000", "#00782e", "#f09b08"),
#' save_plot = TRUE, title = NULL, ylim = c(0, 3), height = 8,
#' width = 10, print_plot = FALSE, image_ext = ".png")
#' }

volcano_plot <- function(e,
                         permutation_comparison = "female_AD_vs_male_AD",
                         channels = c("cfos", "eyfp", "colabel"),
                         colors =  c("#be0000", "#00782e", "#f09b08"),
                         save_plot = TRUE,
                         title = NULL,
                         ylim = c(0, 3),
                         height = 8,
                         width = 10,
                         print_plot = FALSE,
                         plt_theme = NULL,
                         point_size = 1,
                         image_ext = ".png"){
  ## plotting theme
  if (is.null(plt_theme)){
    plt_theme <- ggplot2::theme_classic() + theme(text = element_text(size = 22),
                                                  line = element_line(linewidth = 1),
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


    p <- ggplot(df, aes(x = .data$corr_diff, y = -log10(.data$p_val))) +
      geom_point(size = point_size) +
      geom_point(data = subset(df, df$sig > 0 & df$corr_diff <= -1 | df$sig > 0 & df$corr_diff >= 1), color = colors[k], size = point_size) +
      geom_vline(xintercept = c(-1, 1), color = colors[k], linewidth = 1) +
      geom_hline(yintercept = -log10(alpha), color = colors[k], linewidth = 1) +
      ggplot2::coord_cartesian(xlim = c(-2.1, 2.1)) +
      ylim(ylim) +
      labs(title = title, x = "Correlation Difference", y = "-log(p-value)") +
      plt_theme

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)

      # Create figure directory if it doesn't already exists
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }

      title <- stringr::str_replace(title, "[\\/:*?\"<>|]", "_") %>%  stringr::str_replace(., " ", "_")
      image_file <- file.path(output_dir,
                              paste0("volcano_plot_", title, "_", channels[k], image_ext))
      print(image_file)
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    p_list[[channels[k]]] <- p
  }

  return(p_list)
}



#' Create a parallel coordinate plot
#' @description Plot the correlation difference between two comparison groups into a parallel coordinate plot. The function
#' [SMARTTR::correlation_diff_permutation()] must be run first in order to generate results to plot.
#'
#' @param e experiment object
#' @param permutation_comparison The name of the correlation group comparisons to plot.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) channels to plot
#' @param colors (str, default = c("#be0000", "#00782e", "#f09b08")) Hexadecimal codes corresponding to the channels (respectively) to plot.
#' @param x_label_group_1 (str, NULL) The label for the first group in the permutation analysis. Note: this is to customize the graph labels. It does not reverse the group order.
#' @param x_label_group_2 (str, NULL) The label for the second group in the permutaiton analysis. Note: this is to customize the graph labels. It does not reverse the group order.
#' @param height height of the plot in inches.
#' @param width width of the plot in inches.
#' @param print_plot (bool, default = FALSE) Whether to display the plot (in addition to saving the plot)
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
#' @examples
#' \dontrun{
#' p_list <- parallel_coordinate_plot(e, permutation_comparison = "Context_vs_Shock",
#'  channels = "cfos", colors ="#be0000", x_label_group1 = "Context", x_label_group_2 = "Shock")
#' }
parallel_coordinate_plot <- function(e,
                                     permutation_comparison = "AD_vs_control",
                                     channels = c("cfos", "eyfp", "colabel"),
                                     colors =  c("#be0000", "#00782e", "#f09b08"),
                                     x_label_group_1 = NULL,
                                     x_label_group_2 = NULL,
                                     height = 10,
                                     width = 10,
                                     print_plot = FALSE,
                                     save_plot = TRUE,
                                     reverse_group_order= FALSE,
                                     force = 1,
                                     plt_theme = NULL,
                                     label_size = 30,
                                     image_ext = ".png",
                                     nudge_x = 2:5
){
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

  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){
    alpha <- e$permutation_p_matrix[[permutation_comparison]][[channels[k]]]$alpha
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

    df <- p_vals %>% dplyr::inner_join(corr_diffs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(sigs, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(group_1_pearson, by = c("rowreg", "colreg")) %>%
      dplyr::inner_join(group_2_pearson, by = c("rowreg", "colreg")) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(c(group_1, group_2)), names_to = "group", values_to = "corr")

    df <- df %>% dplyr::filter(.data$sig, abs(.data$corr_diff) >= 1) %>%
      dplyr::mutate(group = factor(.data$group, levels = c(group_1, group_2)),
                    nudge = ifelse(.data$group == group_1, -0.1, 0.1)) %>%
      dplyr::arrange(.data$group, .data$corr_diff) %>%
      mutate(text = paste(.data$rowreg, .data$colreg, sep = "."),
             group_plot = paste(.data$rowreg, .data$colreg, sep = "."))

    if (isTRUE(reverse_group_order)){
      df <- df %>% dplyr::mutate(group = forcats::fct_rev(.data$group),
                                 nudge = ifelse(.data$group == group_2, -0.1, 0.1))
    }

    tryCatch({
      df[seq(2, nrow(df)/2, by = 2),]$text <- ""
      df[seq(nrow(df)/2+1, nrow(df), by = 2),]$text <- ""

    }, error = function(e) {
      message("could not divide by 2. Skipping")
    })

    if (is.null(plt_theme)){
      plt_theme <- ggplot2::theme_classic() +
        theme(text = element_text(size = 22),
              line = element_line(linewidth = 1),
              plot.title = element_text(hjust = 0.5, size = 36),
              axis.ticks.length = unit(5.5, "points"),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black")
        )
    }

    p <- ggplot(df, aes(x = .data$group, y = .data$corr, group = .data$group_plot)) +
      ggplot2::geom_line(alpha = 0.5, color = colors[k], linewidth = 3) +
      ggplot2::geom_point(size = 4, alpha = 0.5, color = colors[k]) +
      ggrepel::geom_text_repel(aes(label = .data$text),
                               size = label_size,
                               color = colors[k], direction = "y",
                               force = force,
                               ylim = c(-1, 1),
                               segment.alpha = 0.3,
                               nudge_x = dplyr::pull(df, .data$nudge)*nudge_x, max.iter = 20000) +
      ggplot2::geom_hline(yintercept = 0,linetype=2,linewidth=1.2) +
      xlab("Group") + ylab("Correlation") +
      ggplot2::expand_limits(y=c(-1,1)) + plt_theme

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0(permutation_comparison, "_parallel_coordinate_plot_",
                                                 channels[k], "_", image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
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
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.
#'  the experiment object output folder.
#' @param edge_color (str, default = "firebrick") Color of the network edges.
#' @param degree_scale_limit (vec, default = c(1,10)) Scale limit for degree size
#' @param graph_theme (default = NULL) Add a [ggraph::theme_graph()] to the network graph. If NULL, the default is taken.
#' @param label_size (default = 5) Default font size for network region labels.
#' @param label_offset (default = 0.15) Distance of label from nodes.
#' @param region_legend (default = TRUE) Boolean determining whether or not to show the region legend categorizing subregions into their largest parent region. Only works well if the Allen ontology is used for the dataset.
#' @param correlation_edge_width_limit (default = c(0.8,1))
#' @param edge_type (default = "arc") "arc" or "diagonal".
#' @param anatomical.colors (vector, default = NULL) NULL defaults to viridis. A named vector of hexadecimal codes for the anatomical super regions.
#' e.g. anatomical.colors = c(Isocortex = "#5571a9", OLF = "#64bdc4", HPF = "#d2875b", CTXsp = "#87a3db", CNU = "#466496", TH = "#7e72af", HY = "#8e7960",  MB = "#d796c8", HB = "#646464")
#' @param edge_thickness_range (default = c(1,5)) Thickness range of the edges.
#' @param node_size_range (default = c(1,8)) Node size range.
#' Can also be a hexadecimal color code written as a string.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_networks(e, network_name = "Shock", channels = "cfos")
#' }
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
                          anatomical.colors = NULL,
                          correlation_edge_width_limit = c(0.8,1),
                          image_ext = ".png",
                          print_plot = FALSE,
                          graph_theme = NULL,
                          label_size = 5,
                          edge_thickness_range = c(1,5),
                          node_size_range = c(1, 8),
                          label_offset = 0.15,
                          save_plot = TRUE){
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels
  for (channel in channels){
    network <- e$networks[[network_name]][[channel]]
    if (is.null(graph_theme)){
      graph_theme <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5,size = 28),
                                                   legend.text = element_text(size = 15),
                                                   legend.title = element_text(size = 15))
    }

    network <- network %>% activate(nodes) %>% dplyr::arrange(.data$super.region)
    if (edge_type == "diagonal"){
      p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
        ggraph::geom_edge_diagonal(aes(color = .data$sign, width = abs(.data$weight)),
                                   edge_alpha = 0.6, n = 1000)

    } else if (edge_type == "arc"){
      p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
        ggraph::geom_edge_arc(aes(color = .data$sign, width = abs(.data$weight)),
                              edge_alpha = 0.6, n = 1000)

    }
    p <- p + ggraph::geom_node_point(aes(size = .data$degree,
                                         color = .data$super.region),
                                     show.legend = TRUE) +
      ggraph::geom_node_text(aes(x = (sqrt(.data$x^2+.data$y^2)+label_offset)*cos(atan(.data$y/.data$x))*sign(.data$x),
                                 y = abs((sqrt(.data$x^2+.data$y^2)+label_offset)*sin(atan(.data$y/.data$x)))*sign(.data$y),
                                 angle = atan(.data$y/.data$x)*180/pi,
                                 label = .data$name),
                             repel = FALSE, color = "grey25",
                             size = label_size,
                             show.legend = NA) +
      ggraph::scale_edge_color_manual(values = c(pos = edge_color,
                                                 neg = "grey20"),
                                      labels = c(pos = "Positive", neg = "Negative"),
                                      name = "Correlation",
                                      guide = guide_legend(order = 1)) +
      ggraph::scale_edge_width(limits=correlation_edge_width_limit,
                               range = edge_thickness_range, name = "Correlation Strength",
                               guide = guide_legend(order = 3))
    if (is.null(anatomical.colors)){
      p <- p + ggraph::scale_color_viridis(name = "Anatomical Region",
                                           discrete = TRUE,
                                           option = "D",
                                           guide = guide_legend(override.aes = list(size=max(node_size_range)), order=4)) +
        ggplot2::scale_size(limits = degree_scale_limit, name="Degree", range = node_size_range,
                            guide = guide_legend(order = 2)) +
        ggplot2::coord_equal() + graph_theme
    } else{
      p <- p + ggplot2::scale_color_manual(name = "Anatomical Region",
                                           values = anatomical.colors,
                                           guide = guide_legend(override.aes = list(size=max(node_size_range)), order=4)) +
        ggplot2::scale_size(limits = degree_scale_limit, name="Degree", range = node_size_range,
                            guide = guide_legend(order = 2)) +
        ggplot2::coord_equal() + graph_theme
    }
    if (is.null(title)){
      title <- paste(network_name, channel)
    }
    p <-  p + ggplot2::ggtitle(title)
    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }
    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("network_", network_name, "_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channel]] <- p
  }
  return(p_list)
}

#' Plot the networks stored in an experiment object
#'
#' @param e experiment object
#' @param correlation_list_names (str vec) character vector of the two correlation lists used to include in a joined network, e.g., `correlation_list_names = c("male_agg", "female_non")`
#' @param title (str, default = NULL) Title of network plot
#' @param channels (str, default = c("cfos", "eyfp", "colabel"))
#' @param height Height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.
#'  the experiment object output folder.
#' @param absolute_weight (bool, default = TRUE) Whether to plot absolute weights. If TRUE, the edge_colors and edge_colors_label should not contain values for positive and negative correlations.
#' @param degree_scale_limit (vec, default = c(1,10)) Scale limit for degree size
#' @param graph_theme (default = NULL) Add a [ggraph::theme_graph()] to the network graph. If NULL, the default is taken.
#' @param label_size (default = 5) Default font size for network region labels.
#' @param label_offset (default = 0.15) Distance of label from nodes.
#' @param region_legend (default = TRUE) Boolean determining whether or not to show the region legend categorizing subregions into their largest parent region. Only works well if the Allen ontology is used for the dataset.
#' @param correlation_edge_width_limit (default = c(0.8, 1)) Range for the width size of the edges.
#' @param edge_colors (vec) vector of hexidecimal codes as strings. Assign a group name to the vector element. e.g. c(male_agg_pos =  "#06537f",
#' male_agg_neg = "#526c7a", female_non_pos = "#C70039", female_non_neg = "#71585f")
#' @param edge_color_labels (vec) vector of edge labels as strings. e.g. c(male_agg_pos = "Positive male",
#' male_agg_neg = "Negative male", female_non_pos = "Positive female", female_non_neg = "Negative female")
#' @param transparent_edge_group1 (bool) logical to render edges transparent
#' @param transparent_edge_group2 (bool) logical to render edges transparent
#' @param edge_thickness_range (default = c(1,5))
#' @param node_size_range (default = c(1, 8))
#' @param anatomical.colors (NuLL or vec, default = NULL) Colors for the parent region as a named vector of hexadecimal regions.
#' Can also be a hexadecimal color code written as a string.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p_list <- plot_joined_networks(anesthesia, correlation_list_names = c("male_agg", "female_non"),
#' channels = "cfos", edge_colors = c(male_agg =  "#06537f", female_non = "#C70039"),
#' edge_color_labels = c(male_agg = "Male aggressor", female_non = "Female non-aggressor"),
#' degree_scale_limit = c(1,45), correlation_edge_width_limit = c(0.8, 1.0),
#' height = 30, width = 30, label_size = 13, label_offset = 0.08, image_ext = ".png")
#' }
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
                                 print_plot = FALSE,
                                 graph_theme = NULL,
                                 transparent_edge_group1 = TRUE,
                                 transparent_edge_group2 = FALSE,
                                 label_size = 5,
                                 label_offset = 0.15,
                                 edge_thickness_range = c(1,5),
                                 node_size_range = c(1, 8),
                                 anatomical.colors = NULL,
                                 save_plot = TRUE){

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
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(color = factor(.data$network))
    } else {
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(color = paste(.data$network, .data$sign, sep = "_"))
    }

    if (isTRUE(transparent_edge_group1) &&  isTRUE(transparent_edge_group2)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = 0)
    } else if (isFALSE(transparent_edge_group1) && isTRUE(transparent_edge_group2)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = if_else(.data$network == p2, 0.0, 0.6))
    } else if (isTRUE(transparent_edge_group1) && isFALSE(transparent_edge_group2)){
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = if_else(.data$network == p1, 0.0, 0.6))
    } else {
      network <- network %>% tidygraph::activate(edges) %>% dplyr::mutate(edge_alpha = 0.6)
    }

    network <- network %>% activate(nodes) %>% dplyr::arrange(.data$super.region)
    p <- ggraph::ggraph(network, layout = "linear", circular = TRUE) +
      ggraph::geom_edge_arc(aes(color = .data$color, width = abs(.data$weight), edge_alpha = .data$edge_alpha),
                            n = 1000) +
      ggraph::geom_node_point(aes(size = .data$degree,
                                  color = .data$super.region),
                              show.legend = TRUE) +
      ggraph::geom_node_text(aes(x = (sqrt(.data$x^2+.data$y^2)+label_offset)*cos(atan(.data$y/.data$x))*sign(.data$x),
                                 y = abs((sqrt(.data$x^2+.data$y^2)+label_offset)*sin(atan(.data$y/.data$x)))*sign(.data$y),
                                 angle = atan(.data$y/.data$x)*180/pi,
                                 label = .data$name),
                             repel = FALSE, color = "grey25",
                             size = label_size,
                             show.legend = NA) +
      ggraph::scale_edge_color_manual(values = edge_colors,
                                      labels = edge_color_labels,
                                      name = "Correlation",
                                      guide = guide_legend(order = 1)) +
      ggraph::scale_edge_width(limits=correlation_edge_width_limit,range = edge_thickness_range, name = "Correlation Strength",
                               guide = guide_legend(order = 3))
    if (is.null(anatomical.colors)){
      p <- p + ggraph::scale_color_viridis(name = "Anatomical Region",
                                           discrete = TRUE,
                                           option = "D",
                                           guide = guide_legend(override.aes = list(size=max(node_size_range)), order=4)) +
        ggplot2::scale_size(limits = degree_scale_limit, name="Degree",range=node_size_range,
                            guide = guide_legend(order = 2)) +
        ggplot2::coord_equal() + graph_theme
    } else{
      p <- p + ggplot2::scale_color_manual(name = "Anatomical Region",
                                           values = anatomical.colors,
                                           guide = guide_legend(override.aes = list(size=max(node_size_range)), order=4)) +
        ggplot2::scale_size(limits = degree_scale_limit, name="Degree",range=node_size_range,
                            guide = guide_legend(order = 2)) +
        ggplot2::coord_equal() + graph_theme
    }
    if (is.null(title)){
      title <- paste(joined_network_name, channel)
    }
    p <-  p + ggplot2::ggtitle(title)
    if (print_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
    }
    p_list[[channel]] <- p

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("network_", joined_network_name, "_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
  }
  return(p_list)
}



#' Plot a permutation difference network
#' e experiment object
#'
#' @param e experiment object
#' @param permutation_group (str) Supply names of the specific analyses found under `e$permutation_p_matrix %>% names()`
#' @param channel (str, default = "cfos")
#' @param layout (str, default = "stress") equivalent to the
#' \href{https://ggraph.data-imaginist.com/reference/ggraph.html}{ggraph} `layout` parameter.
#' @param seed random seed for replication
#' @param degree_scale_limit (vec, default = c(1,15)) Scale limit for degree size
#' @param node_size_range (default = c(1,8)) Node size range.
#' @param edge_color_pos (str, default = "#04C3C3") Color of the edges representing positive correlation differences
#' @param edge_color_neg (str, default = "#BF1F1F") Color of the edges representing negative correlation differences
#' @param pos_edge_legend_label (str, default = "") Legend label for positive difference edges, e.g. "Iso > Ctrl"
#' @param neg_edge_legend_label (str, default = "") Legend label for negative difference edges, e.g. "Iso < Ctrl"
#' @param title (str, default = "") Plot title.
#' @param correlation_edge_width_limit (default = c(0, 2)) Range of correlation differences in the network.
#' @param edge_thickness_range (default = c(0,6))  Range for the width size of the edges.
#' @param node_text_size (default = 12)  Font size for node text.
#' @param graph_theme (default = NULL) Add a [ggraph::theme_graph()] to the network graph. If NULL, the default is taken.
#' @param mark_hull  (bool, default = TRUE) Plot detected communities using a convex hull polygon.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param height Height of the plot in inches.
#' @param width width of the plot in inches.
#' @param image_ext (default = ".png")
#' @param anatomical.colors (vec, default = c(Isocortex = "#70ff71", OLF = "#75aadb", HPF = "#7ed04b", CTXsp = "#8ada87", CNU = "#98d6f9",
#' TH = "#ff7080", HY = "#e64438",  MB = "#ff64ff", HB = "#ff9b88")) Colors for the parent region as a named vector of hexadecimal regions.
#' @return p A plot handle.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_perm_diff_network(anesthesia,
#' permutation_group = "Control_vs_Isoflurane",
#' channel = "cfos",
#' pos_edge_legend_label = "Iso > Ctrl",
#' neg_edge_legend_label = "Iso < Ctrl",
#' layout = "fr",
#' mark_hull=TRUE,
#' node_text_size = 14,
#' seed = 15,
#' width = 20,
#' height = 20,
#' graph_theme = graph_theme,
#' anatomical.colors = anatomical.colors,
#' image_ext = ".png")
#' }
plot_perm_diff_network <- function(e,
                                   permutation_group,
                                   channel = "cfos",
                                   layout ="stress",
                                   seed = 30,
                                   degree_scale_limit= c(1,15),
                                   node_size_range = c(3,15),
                                   edge_color_pos = "#04C3C3",
                                   edge_color_neg =  "#BF1F1F",
                                   anatomical.colors = c(Isocortex = "#70ff71", OLF = "#75aadb", HPF = "#7ed04b",
                                                         CTXsp = "#8ada87", CNU = "#98d6f9", TH = "#ff7080",
                                                         HY = "#e64438",  MB = "#ff64ff", HB = "#ff9b88"),
                                   pos_edge_legend_label = "",
                                   neg_edge_legend_label = "",
                                   title = "",
                                   correlation_edge_width_limit = c(0,2),
                                   edge_thickness_range = c(0,6),
                                   node_text_size = 12,
                                   graph_theme = NULL,
                                   mark_hull = TRUE,
                                   save_plot = TRUE,
                                   height=15,
                                   width=15,
                                   image_ext = ".png"){

  if (is.null(graph_theme)){
    graph_theme <- ggraph::theme_graph() + theme(plot.title = element_text(hjust = 0.5,size = 28),
                                                 legend.text = element_text(size = 15),
                                                 legend.title = element_text(size = 15))
  }
  perm_network <- e$permutation_diff_network[[permutation_group]][[channel]]
  set.seed(seed)
  p <- ggraph::ggraph(perm_network, layout = layout) +
    ggraph::geom_edge_link(aes(color = sign, width = abs(.data$weight)),
                           edge_alpha = 0.6, n = 1000) +
    ggraph::geom_node_point(aes(size = .data$degree, color = .data$super.region)) +
    ggraph::geom_node_text(aes(label = .data$name), size =  node_text_size, color = "black", repel = TRUE, check_overlap = TRUE) +
    ggplot2::scale_color_manual(name = "Anatomical Region",
                                values = anatomical.colors,
                                guide = guide_legend(override.aes = list(size=max(node_size_range)), order=4)) +
    ggplot2::coord_equal() +
    ggraph::scale_edge_color_manual(values = c(pos = edge_color_pos,
                                               neg = edge_color_neg),
                                    labels = c(pos = pos_edge_legend_label, neg = neg_edge_legend_label),
                                    name = "Functional direction",
                                    guide = guide_legend(order = 1)) +

    ggraph::scale_edge_width(limits=correlation_edge_width_limit,
                             range = edge_thickness_range, name = "Pearson correlation difference",
                             guide = guide_legend(order = 3))  +
    ggplot2::scale_size(limits = degree_scale_limit, name ="Degree",
                        range = node_size_range,
                        guide = guide_legend(order = 2)) +
    ggplot2::ggtitle(title) +
    graph_theme

  if (mark_hull){
    p <- p +  ggforce::geom_mark_hull(aes(.data$x, .data$y, fill = .data$community, color = .data$community), alpha = .1,
                                      concavity = 5, radius = unit(6, "mm"), show.legend = FALSE) +
      ggplot2::scale_fill_viridis_d(option = "D")
    ### add an if/else statement here
    # scale_fill_manual(values = c( "#3c9ea5", "#2f3480", "#9d1c7b", "#ff7f00", "#4daf4a", "#e41a1c","#984ea3",  "#bdd9bf", "#377eb8","#27ae60"))
  }

  if (save_plot){
    dev.new(width = width, height = height, noRStudioGD=TRUE)
    print(p)
    output_dir <-  file.path(attr(e, "info")$output_path, "figures")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }
    image_file <- file.path(output_dir,
                            paste0("permutation_difference_network_", permutation_group,"_mark_hull_", mark_hull, image_ext))
    ggsave(filename = image_file,  width = width, height = height, units = "in")
  }
  return(p)
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
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param theme.gg (default = NULL) Option to use custom ggplot2 theme if the user wants
#' @param image_ext (default = ".png") image extension to the plot as.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_degree_distributions(e, channels = "cfos",
#' labels = c(female_AD = "female_AD_label", female_control = "female_control_label"),
#' title = "my title", image_ext = ".png")
#' }

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
                                      print_plot = FALSE,
                                      theme.gg = NULL,
                                      save_plot = TRUE){
  if (is.null(theme.gg)){
    theme.gg <- ggplot2::theme_classic() +
      theme(text = element_text(size = 22),
            line = element_line(linewidth = 1),
            plot.title = element_text(hjust = 0.5, size = 36),
            axis.ticks.length = unit(5.5,"points"))
  }

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
                         aes(.data$degree)) +
      ggplot2::geom_bar(aes(fill = .data$group), color="black") +
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
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_degree_distributions_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}



#' Plot the mean degree of the networks in a barplot. Error bars are plotted as SEM.
#'
#' @param e experiment object
#' @param labels (str) The legend labels to correspond with your network names, e.g. labels = c(network1_name = "network 1 label", network2_name = "network 2 label).
#' These are the same network names used in the function [SMARTTR::summarise_networks()].
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param color_palettes (str, default = c("reds", "greens")) Color palettes from [grDevices::hcl.colors] that are used to for plotting networks for each channel, respectively.
#' @param colors_manual (str, default = NULL ) Manually choose the hexadecimal color codes to create a custom color palette, e.g. colors_manual = c("#660000", "#FF0000", "#FF6666").
#' Warning: this will be applied to all channels. It's recommended to set the channels parameter to a single channel if this parameter is used.
#' @param title (str, default = "my_title) plot title
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#'  the experiment object output folder.
#' @param label_angle (int, default = 60)
#' @param height (int, default = 10) Height of the plot in inches.
#' @param width (int, default = 10) Width of the plot in inches.
#' @param ylim (vec, default = c(0,10)) Axes limits of y-axis
#' @param theme.gg (default = NULL) Option to use custom ggplot2 theme if the user wants
#' @param image_ext (default = ".png") image extension to the plot as.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_mean_degree(e, colors_manual = c("#660000", "#FF0000"),
#' channels = "cfos", labels = c("AD" = "AD_label", "control" = "control_label"),
#' title = "my title", ylim = c(0,100), image_ext = ".png")
#' }

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
                             theme.gg = NULL,
                             image_ext = ".png",
                             print_plot = FALSE,
                             save_plot = TRUE){
  if (is.null(theme.gg)){
    theme.gg <- ggplot2::theme_classic() +
      theme(text = element_text(size = 22),
            line = element_line(linewidth = 1),
            plot.title = element_text(hjust = 0.5, size = 36),
            axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
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

    if (rev_x_scale){
      p <- ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(.data$group, dplyr::desc(.data$group)), .data$degree.mean))
    } else{

      p <- ggplot(e$networks_summaries[[channel]]$networks_stats, aes(.data$group, .data$degree.mean))
    }
    p <- p +
      ggplot2::geom_col(aes(fill = .data$group), color = "black", linewidth = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = .data$degree.mean - (.data$degree.sd/sqrt(.data$n.nodes)),
                                 ymax = .data$degree.mean + (.data$degree.sd/sqrt(.data$n.nodes))),
                             width = 0.2, linewidth = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") + ylab("Mean Degree") +
      ylim(ylim) +
      theme.gg
    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }
    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }
    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_degree_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}

#' Plot mean clustering coefficient
#' @description
#' Plot the mean clustering coefficients of the networks in a barplot. Error bars are plotted as SEM.
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
#' @param theme.gg (default = NULL) Option to use custom ggplot2 theme if the user wants
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_mean_clust_coeff(e, colors_manual = c("#660000", "#FF0000"),
#' channels = "cfos", labels = c("AD" = "AD_label", "control" = "control_label"),
#' title = "my title", ylim = c(0, 0.7), image_ext = ".png")
#' }

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
                                  theme.gg = NULL,
                                  image_ext = ".png",
                                  print_plot = FALSE,
                                  save_plot = TRUE){
  if (is.null(theme.gg)){
    theme.gg <- ggplot2::theme_classic() +
      theme(text = element_text(size = 22),
            line = element_line(linewidth = 1),
            plot.title = element_text(hjust = 0.5, size = 36),
            axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
            axis.ticks.length = unit(5.5,"points"))
  }

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

    if (rev_x_scale){
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(.data$group, dplyr::desc(.data$group)), .data$clust.coef.mean))
    } else{

      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(.data$group, .data$clust.coef.mean))
    }

    p <- p +
      ggplot2::geom_col(aes(fill = .data$group), color = "black", linewidth = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = .data$clust.coef.mean - (.data$clust.coef.sd/sqrt(.data$n.nodes)),
                                 ymax = .data$clust.coef.mean + (.data$clust.coef.sd/sqrt(.data$n.nodes))),
                             width = 0.2, linewidth = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") +
      ylab("Mean Clustering Coefficient") +
      ylim(ylim) +
      theme.gg

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_clust_coeff_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
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
#' @param theme.gg (default = NULL) Option to use custom ggplot2 theme if the user wants
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_mean_global_effic(e, colors_manual = c("#660000", "#FF0000"),
#' channels = "cfos", labels = c("AD" = "AD_label", "control" = "control_label"),
#' title = "my title", ylim = c(0, 0.7), image_ext = ".png")
#' }

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
                                   theme.gg = NULL,
                                   image_ext = ".png",
                                   print_plot = FALSE,
                                   save_plot = TRUE){
  if (is.null(theme.gg)){
    theme.gg <- ggplot2::theme_classic() +
      theme(text = element_text(size = 22),
            line = element_line(linewidth = 1),
            plot.title = element_text(hjust = 0.5, size = 36),
            axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
            axis.ticks.length = unit(5.5,"points"))
  }

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

    if (rev_x_scale){
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(.data$group, dplyr::desc(.data$group)), .data$efficiency.mean))
    } else{
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(.data$group, .data$efficiency.mean))
    }

    p <- p +
      ggplot2::geom_col(aes(fill = .data$group), color = "black", linewidth = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = .data$efficiency.mean - (.data$efficiency.sd/sqrt(.data$n.nodes)),
                                 ymax = .data$efficiency.mean + (.data$efficiency.sd/sqrt(.data$n.nodes))),
                             width = 0.2, linewidth = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") +
      ylab("Global Efficiency") +
      ylim(ylim) +
      theme.gg

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_global_efficiency_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
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
#' @param label_angle (int, default = 60)
#' @param ylim (vec, default = c(0,10)) Axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.s
#' @param rev_x_scale (bool, default = FALSE) Reveres the scale of the categorical variables
#' @param theme.gg (default = NULL) Option to use custom ggplot2 theme if the user wants
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_mean_between_centrality(e, colors_manual = c("#660000", "#FF0000"),
#' channels = "cfos", labels = c("AD" = "AD_label", "control" = "control_label"),
#' title = "my title", ylim = c(0, 50), image_ext = ".png")
#' }

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
                                         theme.gg = NULL,
                                         image_ext = ".png",
                                         print_plot = FALSE,
                                         save_plot = TRUE){
  if (is.null(theme.gg)){
    theme.gg <- ggplot2::theme_classic() +
      theme(text = element_text(size = 22),
            line = element_line(linewidth = 1),
            plot.title = element_text(hjust = 0.5, size = 36),
            axis.text.x = element_text(angle = label_angle, hjust = 1, color = "black"),
            axis.ticks.length = unit(5.5,"points"))
  }

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
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(x = stats::reorder(.data$group, dplyr::desc(.data$group)), .data$btw.mean))
    } else{
      p <- ggplot2::ggplot(e$networks_summaries[[channel]]$networks_stats, aes(.data$group, .data$btw.mean))
    }

    p <- p +
      ggplot2::geom_col(aes(fill = .data$group), color = "black", linewidth = 1, width = 0.6) +
      ggplot2::geom_errorbar(aes(ymin = .data$btw.mean - (.data$btw.sd/sqrt(.data$n.nodes)),
                                 ymax = .data$btw.mean + (.data$btw.sd/sqrt(.data$n.nodes))),
                             width = 0.2, linewidth = 1) +
      ggplot2::scale_fill_manual(values = colors,
                                 name = "Group",
                                 guide="none") +
      ggplot2::scale_x_discrete(labels = labels) +
      xlab("Group") +
      ylab("Mean Betweenness Centrality") +
      ylim(ylim) +
      theme.gg

    if (!is.null(title)){
      title <- paste(title)
      p <-  p + ggplot2::ggtitle(title)
    }

    if (print_plot){
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)

      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_mean_betweenness_centrality_", channel, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }

    p_list[[channels[k]]] <- p
  }
  return(p_list)
}




#' Plot the degree distributions across regions
#' @description
#' Bar plot of degree per region in descending magnitude
#' @param e experiment object
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) Channels to plot
#' @param height (int, default = 15) Height of the plot in inches.
#' @param width (int, default = 20) Width of the plot in inches.
#' @param ylim (vec, default = c(0,15))axes limits of y-axis
#' @param image_ext (default = ".png") image extension to the plot as.
#' @param save_plot (bool, default = TRUE) Save into the figures subdirectory of the
#'  the experiment object output folder.
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.s
#' @param colors (str, default = ) String vector of hexadecimal color codes corresponding to to each channel plotted.
#' @param network (str, default = "AD") Which network to plot the degree distribution across regions
#' @param filter_isolates (default = TRUE) Avoid plotting isolated nodes (zero value)
#' @param title (str, default = "") Plot title.
#' @param sort_super_region (bool, default = FALSE) Whether to divide into subfacets based on which parent region
#' @param region_label_angle (int, default = 60) Angle of region labels.
#' @param label_text_size (int, default = 12) Font size of region labels.
#' @param theme.bar User option to use their own ggplot theme
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_degree_regions(e, colors = c("#660000", "#FF0000"),
#' channels = "cfos", region_label_angle = 60,
#' ylim = c(0, 15), image_ext = ".png")
#' }

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
                                filter_isolates = TRUE,
                                image_ext = ".png",
                                print_plot = FALSE,
                                save_plot = TRUE,
                                theme.bar = NULL){
  if (is.null(theme.bar)){
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
            plot.margin = ggplot2::margin(1,1.5,0,1.5, "cm"),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank())
  }

  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels

  for (k in 1:length(channels)){
    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()

    if (sort_super_region) {
      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(.data$group == network) %>%
        dplyr::arrange(.data$super.region,  dplyr::desc(.data$degree)) %>%
        dplyr::mutate(name = factor(.data$name, levels = .data$name))

      if (isTRUE(filter_isolates)){
        df <- df %>% dplyr::filter(.data$degree > 0)
      }

      p <- df %>%
        ggplot2::ggplot(aes(.data$name, .data$degree)) +
        ggplot2::geom_col(fill = colors[[k]], color = "black") +
        ggplot2::facet_grid(~.data$super.region, scales = "free_x", space = "free_x", switch = "x")  +
        xlab("Brain Region") + ylab("Degree") +
        ylim(ylim)  + theme.bar

    } else {
      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(.data$group == network) %>%
        dplyr::arrange(dplyr::desc(.data$degree)) %>%
        dplyr::mutate(name = factor(.data$name, levels = .data$name))

      if (isTRUE(filter_isolates)){
        df <- df %>% dplyr::filter(.data$degree > 0)
      }

      p <- df %>%
        ggplot2::ggplot(aes(.data$name, .data$degree)) +
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
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_degree_per_region_", channel, "_",network, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
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
#' @param print_plot (bool, default = FALSE) Whether to print the plot as an output.s
#' @param colors (str, default = ) String vector of hexadecimal color codes corresponding to to each channel plotted.
#' @param network (str, default = "AD") Which network to plot the betweenness distribution across regions
#' @param title (str, default = "")
#' @param sort_super_region (bool, default = FALSE) Whether to divide into subfacets based on which parent region
#' @param region_label_angle (int, default = 60) Angle of region labels.
#' @param label_text_size (int, default = 12) Font size of region labels.
#' @param filter_isolates (default = TRUE) Avoid plotting isolated nodes (zero value)
#' @param theme.bar User option to use their own ggplot theme
#'  the experiment object output folder.
#' @return p_list A list the same length as the number of channels, with each element containing a plot handle for that channel.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_betweenness_regions(e, colors ="#660000", channels = "cfos", region_label_angle = 60,
#' ylim = c(0, 50), image_ext = ".png")
#' }
plot_betweenness_regions <- function(e,
                                     channels = c("cfos", "eyfp"),
                                     colors = c("red", "green"),
                                     network = "AD",
                                     title = "",
                                     height = 10,
                                     width = 20,
                                     ylim = c(0, 15),
                                     filter_isolates = TRUE,
                                     sort_super_region = FALSE,
                                     region_label_angle = 60,
                                     label_text_size = 12,
                                     image_ext = ".png",
                                     print_plot = FALSE,
                                     save_plot = TRUE,
                                     theme.bar = NULL){
  if (is.null(theme.bar)){
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
            plot.margin = ggplot2::margin(1,1.5,0,1.5, "cm"),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank())
  }
  p_list <- vector(mode='list', length = length(channels))
  names(p_list) <- channels
  for (k in 1:length(channels)){
    channel <- channels[[k]]
    n_groups <-  e$networks_summaries[[channel]]$networks_degree_distrib$group %>%
      unique() %>% length()

    if (sort_super_region) {
      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(.data$group == network) %>%
        dplyr::arrange(.data$super.region,  dplyr::desc(.data$btw)) %>%
        dplyr::mutate(name = factor(.data$name, levels = .data$name))

      if (isTRUE(filter_isolates)){
        df <- df %>% dplyr::filter(.data$degree > 0)
      }

      p <- df %>%
        ggplot2::ggplot(aes(.data$name, .data$btw)) +
        ggplot2::geom_col(fill = colors[[k]], color = "black") +
        ggplot2::facet_grid(~.data$super.region, scales = "free_x", space = "free_x", switch = "x")  +
        xlab("Brain Region") + ylab("Betweenness") +
        ylim(ylim)  + theme.bar

    } else {

      df <- e$networks_summaries[[channel]]$networks_nodes %>% dplyr::filter(.data$group == network) %>%
        dplyr::arrange(dplyr::desc(.data$btw)) %>%
        dplyr::mutate(name = factor(.data$name, levels = .data$name))

      if (isTRUE(filter_isolates)){
        df <- df %>% dplyr::filter(.data$degree > 0)
      }

      p <- df %>%
        ggplot2::ggplot(aes(.data$name, .data$btw)) +
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
      dev.new(noRStudioGD=TRUE)
      print(p)
    }

    if(save_plot){
      dev.new(width = width, height = height, noRStudioGD=TRUE)
      print(p)
      output_dir <-  file.path(attr(e, "info")$output_path, "figures")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }
      image_file <- file.path(output_dir, paste0("networks_betweenness_per_region_", channel, "_", network, image_ext))
      ggsave(filename = image_file,  width = width, height = height, units = "in")
    }
    p_list[[channels[k]]] <- p
  }
  return(p_list)
}


