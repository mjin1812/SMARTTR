# filtering the regions of a dataset works

    Code
      e$combined_normalized_counts$cfos %>% head()
    Output
      # A tibble: 6 x 9
      # Groups:   acronym [6]
        mouse_ID group   acronym name                count area.mm2 volume.mm3 normalized.count.by.~1 normalized.count.by.~2
           <dbl> <chr>   <chr>   <chr>               <int>    <dbl>      <dbl>                  <dbl>                  <dbl>
      1      829 Context LHA     Lateral hypothalam~  1467    2.94     0.0264                    499.                 55488.
      2      829 Context LPO     Lateral preoptic a~   471    1.07     0.00961                   441.                 49027.
      3      829 Context MPO     Medial preoptic ar~   476    0.651    0.00586                   731.                 81232.
      4      829 Context PH      Posterior hypothal~  2021    2.18     0.0196                    929.                103167.
      5      829 Context PVH     Paraventricular hy~    67    0.162    0.00146                   412.                 45820.
      6      829 Context STN     Subthalamic nucleus    31    0.272    0.00245                   114.                 12677.
      # i abbreviated names: 1: normalized.count.by.area, 2: normalized.count.by.volume

