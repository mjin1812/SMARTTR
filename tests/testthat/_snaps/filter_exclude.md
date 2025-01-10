# filtering the regions of a dataset works

    Code
      e$combined_normalized_counts$cfos %>% head()
    Output
      # A tibble: 6 x 9
      # Groups:   acronym [6]
        mouse_ID group  acronym name  count area.mm2 volume.mm3 normalized.count.by.~1
           <dbl> <chr>  <chr>   <chr> <int>    <dbl>      <dbl>                  <dbl>
      1      829 Conte~ LHA     Late~  1467    2.94     0.0264                    499.
      2      829 Conte~ LPO     Late~   471    1.07     0.00961                   441.
      3      829 Conte~ MPO     Medi~   476    0.651    0.00586                   731.
      4      829 Conte~ PH      Post~  2021    2.18     0.0196                    929.
      5      829 Conte~ PVH     Para~    67    0.162    0.00146                   412.
      6      829 Conte~ STN     Subt~    31    0.272    0.00245                   114.
      # i abbreviated name: 1: normalized.count.by.area
      # i 1 more variable: normalized.count.by.volume <dbl>

