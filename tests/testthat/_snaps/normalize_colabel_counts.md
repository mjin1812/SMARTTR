# normalizing colabelled counts works

    Code
      e$combined_normalized_counts$colabel_over_eyfp %>% head()
    Output
      # A tibble: 6 x 9
      # Groups:   acronym [6]
        mouse_ID group  acronym name  count area.mm2 volume.mm3 normalized.count.by.~1
           <dbl> <chr>  <chr>   <chr> <dbl>    <dbl>      <dbl>                  <dbl>
      1      829 Conte~ ACA     Ante~ 0.617     2.74     0.0247                  0.617
      2      829 Conte~ ACB     Nucl~ 0.367     4.13     0.0371                  0.367
      3      829 Conte~ AI      Agra~ 0.272     1.63     0.0147                  0.272
      4      829 Conte~ AUDd    Dors~ 0.267     3.99     0.0359                  0.267
      5      829 Conte~ AUDp    Prim~ 0.384     3.35     0.0301                  0.384
      6      829 Conte~ AUDv    Vent~ 0.295     3.53     0.0318                  0.295
      # i abbreviated name: 1: normalized.count.by.area
      # i 1 more variable: normalized.count.by.volume <dbl>

---

    Code
      e$combined_normalized_counts$colabel_over_cfos %>% head()
    Output
      # A tibble: 6 x 9
      # Groups:   acronym [6]
        mouse_ID group acronym name   count area.mm2 volume.mm3 normalized.count.by.~1
           <dbl> <chr> <chr>   <chr>  <dbl>    <dbl>      <dbl>                  <dbl>
      1      829 Cont~ ACA     Ante~ 0.115      2.74     0.0247                 0.115 
      2      829 Cont~ ACB     Nucl~ 0.0209     4.13     0.0371                 0.0209
      3      829 Cont~ AI      Agra~ 0.0333     1.63     0.0147                 0.0333
      4      829 Cont~ AUDd    Dors~ 0.149      3.99     0.0359                 0.149 
      5      829 Cont~ AUDp    Prim~ 0.154      3.35     0.0301                 0.154 
      6      829 Cont~ AUDv    Vent~ 0.0988     3.53     0.0318                 0.0988
      # i abbreviated name: 1: normalized.count.by.area
      # i 1 more variable: normalized.count.by.volume <dbl>

