

#Calculate % colabeled over cfos and eyfp. Specify which ROIS you want







##


get_correlations <- function(e, by = c("sex", "group"), values = c("female", "AD"),
                             channels = c("cfos", "eyfp", "colabel"),  method = "pearson"){


  # Reformat the data into wide format
  pivoted <- e$combined_normalized_counts$cfos %>% dplyr::filter(sex == "female", group == "AD") %>%
    dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
    tidyr::pivot_wider(  names_from = acronym, values_from = normalized.count.by.volume)

  # Perform the correlations, ignore the mouse_ID and group columns
  corrs <- pivoted %>% dplyr::select(-(mouse_ID:group)) %>% corrr::correlate(method = method)





}





# Generate correlational heatmaps



plot_heatmaps <- function(e, by = c("group", "sex"), channel = c("cfos", "eyfp", "colabel"),
                            colors = c("red", "green", "blue")){


  Hmisc





}


install.packages("Hmisc")

library(palmerpenguins)
library(dplyr)
penguins_cor <- palmerpenguins::penguins %>%
  dplyr::select(bill_length_mm, bill_depth_mm, flipper_length_mm) %>%
  correlate()




# Reformat the data into wide format

pivoted <- e$combined_normalized_counts$cfos %>% dplyr::filter(sex == "female", group == "AD") %>%
  dplyr::select(mouse_ID:acronym, normalized.count.by.volume) %>%
  tidyr::pivot_wider(  names_from = acronym, values_from = normalized.count.by.volume)


# Perform the correlations, ignore the mouse_ID and group columns
corrs <- pivoted %>% dplyr::select(-(mouse_ID:group)) %>% corrr::correlate(method = "pearson")


# plot heatmap using the intrinsic heatmap function
corrs %>% shave() %>% rplot()

# plot a network plot using
  network_plot(min_cor = .2)






  p <- ggplot(cross_region_fractions, aes(train_region,
                                          test_region,
                                          fill = total_FP_fraction,
                                          text = text)) +
    geom_tile() +
    geom_text(aes(label = text)) +
    scale_fill_gradient(low = "white", high = "#cc3337") +
    labs(x = "Train region", y = "Test region", fill = "FP\n fraction\n retained")












