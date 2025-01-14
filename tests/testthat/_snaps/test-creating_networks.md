# Test different parameters to create single networks

    Code
      e$networks$Shock$eyfp %>% tidygraph::activate(nodes) %>% tibble::as_tibble()
    Output
      # A tibble: 78 x 9
         name  super.region degree triangles clust.coef avg.dist efficiency    btw
         <chr> <fct>         <dbl>     <dbl>      <dbl>    <dbl>      <dbl>  <dbl>
       1 RSP   Isocortex         6        10      0.667     2.44      0.418   5.14
       2 PTLp  Isocortex         7         9      0.429     2.39      0.431  51.5 
       3 VISC  Isocortex        18        64      0.418     1.99      0.540 148.  
       4 ACA   Isocortex        20        82      0.432     1.91      0.563 234.  
       5 PERI  Isocortex         6         2      0.133     2.45      0.416 124.  
       6 TEa   Isocortex        13        47      0.603     2.23      0.479  43.8 
       7 AI    Isocortex        20        74      0.389     1.92      0.558 189.  
       8 ECT   Isocortex         8        14      0.5       2.60      0.414  96.6 
       9 GU    Isocortex        10        18      0.4       2.36      0.452  94.5 
      10 ILA   Isocortex        10        16      0.356     2.57      0.427  73.7 
      # i 68 more rows
      # i 1 more variable: group <fct>

