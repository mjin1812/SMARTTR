test_that("testing plot_normalized_counts", {
  skip_on_cran()
  skip_on_ci()

  e <- make_analyzed_experiment()
  expect_error(plot_normalized_counts(e,
                         by = "group",
                         values = list("Context", "Shock"),
                         channels = "eyfp",
                         print_plot = FALSE,
                         save_plot = FALSE,
                         flip_axis = FALSE),
                  regexp = NA)

  expect_error(plot_normalized_counts(e,
                         by = "group",
                         values = list("Context", "Shock"),
                         channels = "eyfp",
                         print_plot = FALSE,
                         save_plot = TRUE,
                         height = 10,
                         width = 4,
                         image_ext = ".jpg",
                         flip_axis = TRUE),
               regexp = NA)

  file_path <- file.path(attr(e, "info")$output_path, "figures/eyfp_normalized_counts.jpg")
  expect_true(file.exists(file_path))
})


test_that("plotting and saving correlation heatmaps works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_analyzed_experiment()
  expect_warning(expect_warning(
  p <- plot_correlation_heatmaps(e,
                                 correlation_list_name = "Shock",
                                 channels = "eyfp",
                                 colors ="#00782e",
                                 save_plot = TRUE,
                                 image_ext = ".jpg",
                                 height = 5,
                                 width = 5),
  regexp = "missing values or values outside the scale range"),
  regexp =  "missing values or values outside the scale range")
  expect_s3_class(p$eyfp, "ggplot")
  file_path <- file.path(attr(e, "info")$output_path, "figures/heatmap_Shock_eyfp.jpg" )
  expect_true(file.exists(file_path))
  graphics.off()
})

test_that("plotting a volcano plot works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_analyzed_experiment()
  p <- volcano_plot(e,
                    permutation_comparison = "Context_vs_Shock",
                    channels = "eyfp",
                    colors = "#00782e",
                    save_plot = FALSE,
                    title = NULL,
                    ylim = c(0, 4))
  expect_silent(print(p$eyfp))
})


test_that("plotting a parallel coordinate plot works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_analyzed_experiment()
  p <- parallel_coordinate_plot(e,
                                permutation_comparison = "Context_vs_Shock",
                                channels = "eyfp",
                                colors = "#00782e",
                                save_plot = FALSE,
                                x_label_group_1 = "Context",
                                x_label_group_2 = "Shock",
                                force = 0,
                                label_size = 2,
                                nudge_x = 2:5)

  expect_s3_class(p$eyfp, "ggplot")
})


test_that("plotting networks works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_analyzed_experiment()
  p <- plot_networks(e,
                     network_name = "Shock",
                     channels = "eyfp",
                     edge_color = "#00782e",
                     save_plot = FALSE,
                     degree_scale_limit = c(1,20),
                     correlation_edge_width_limit = c(0.8,1),
                     print_plot = FALSE,
                      graph_theme = NULL,
                      label_size = 1,
                      edge_thickness_range = c(.2,1),
                      node_size_range = c(1, 4),
                      label_offset = 0.15)
 expect_s3_class(p$eyfp, "ggplot")
})


test_that("plotting joined networks works", {
  skip_on_cran()
  skip_on_ci()
  e <- make_analyzed_experiment()
  p <- plot_joined_networks(e,
                            correlation_list_names = c("Context", "Shock"),
                            channels = "eyfp",
                            absolute_weight = TRUE,
                            edge_colors = c(Context =  "#526c7a",
                                            Shock = "#00782e"),
                            edge_color_labels = c(Context = "Context",
                                                  Shock = "Shock"),
                            degree_scale_limit = c(1,30),
                            region_legend = TRUE,
                            correlation_edge_width_limit = c(0.8,1),
                            print_plot = FALSE,
                            graph_theme = NULL,
                            transparent_edge_group1 = FALSE,
                            transparent_edge_group2 = FALSE,
                            label_size = 1,
                            label_offset = 0.1,
                            edge_thickness_range = c(1,2),
                            node_size_range = c(1, 8),
                            anatomical.colors = NULL,
                            save_plot = FALSE)
 expect_s3_class(p$eyfp, "ggplot")
})

