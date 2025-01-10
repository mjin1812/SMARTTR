# Getting cell table works

    Code
      m$cell_table$cfos %>% head()
    Output
      # A tibble: 6 x 15
        animal    AP     x     y intensity  area    id color   right.hemisphere     ML
         <dbl> <dbl> <dbl> <dbl>     <dbl> <int> <dbl> <chr>   <lgl>             <dbl>
      1    733 -2.14 3876.  328.      76.8    62   241 #009fac FALSE            -0.967
      2    733 -2.14 3864.  335.      81.8   149   241 #009fac FALSE            -0.976
      3    733 -2.14 3838.  346.      99.1   244   241 #009fac FALSE            -1.01 
      4    733 -2.14 3868.  348.      95.6    84   241 #009fac FALSE            -0.969
      5    733 -2.14 3790.  350.      64.1    53   241 #009fac FALSE            -1.06 
      6    733 -2.14 3528.  360.      74.4    60   532 #009fac FALSE            -1.39 
      # i 5 more variables: DV <dbl>, acronym <chr>, name <chr>, image <chr>,
      #   slice_name <chr>

