# Flexible combining of cell counts

    Code
      e$combined_counts_per_slice_list$cfos %>% head()
    Output
      # A tibble: 6 x 12
      # Groups:   slice_name, AP, right.hemisphere, acronym [6]
        mouse_ID age   sex   group   slice_name    AP right.hemisphere acronym name   
           <dbl> <chr> <chr> <chr>   <chr>      <dbl> <lgl>            <chr>   <chr>  
      1      829 5.41  male  Context 1_1         0.79 FALSE            ACA     Anteri~
      2      829 5.41  male  Context 1_1         0.79 TRUE             ACA     Anteri~
      3      829 5.41  male  Context 1_1         0.79 FALSE            ACB     Nucleu~
      4      829 5.41  male  Context 1_1         0.79 TRUE             ACB     Nucleu~
      5      829 5.41  male  Context 1_1         0.79 FALSE            AI      Agranu~
      6      829 5.41  male  Context 1_1         0.79 TRUE             AI      Agranu~
      # i 3 more variables: count <int>, area.mm2 <dbl>, volume.mm3 <dbl>

---

    Code
      e$combined_normalized_counts$cfos %>% head()
    Output
      # A tibble: 6 x 11
      # Groups:   acronym [6]
        mouse_ID age   sex   group   acronym name            count area.mm2 volume.mm3
           <dbl> <chr> <chr> <chr>   <chr>   <chr>           <int>    <dbl>      <dbl>
      1      829 5.41  male  Context AAA     Anterior amygd~   293    0.704    0.00634
      2      829 5.41  male  Context ACA     Anterior cingu~  3400    2.74     0.0247 
      3      829 5.41  male  Context ACB     Nucleus accumb~  1054    4.13     0.0371 
      4      829 5.41  male  Context AI      Agranular insu~   661    1.63     0.0147 
      5      829 5.41  male  Context APN     Anterior prete~   397    2.07     0.0186 
      6      829 5.41  male  Context ARH     Arcuate hypoth~   276    0.492    0.00443
      # i 2 more variables: normalized.count.by.area <dbl>,
      #   normalized.count.by.volume <dbl>

