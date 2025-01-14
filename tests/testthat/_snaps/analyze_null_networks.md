# test that rewire_network() and summarizing networks works as expected

    Code
      summarytable$eyfp
    Output
      # A tibble: 780 x 9
         name  super.region degree triangles clust.coef avg.dist efficiency    btw
         <chr> <fct>         <dbl>     <dbl>      <dbl>    <dbl>      <dbl>  <dbl>
       1 RSP   Isocortex         5         1     0.1        2.05      0.318  58.2 
       2 PTLp  Isocortex         3         0     0          2.35      0.275  21.5 
       3 VISC  Isocortex         3         0     0          2.58      0.256   9.94
       4 ACA   Isocortex         2         0     0          2.65      0.243   4.5 
       5 PERI  Isocortex         5         0     0          2.01      0.323  79.0 
       6 TEa   Isocortex         8         3     0.107      1.90      0.358 148.  
       7 AI    Isocortex         4         0     0          2.05      0.313  61.2 
       8 ECT   Isocortex         9         2     0.0556     1.73      0.383 226.  
       9 GU    Isocortex         5         2     0.2        1.97      0.328  59.3 
      10 ILA   Isocortex         1         0     0          2.77      0.226   0   
      # i 770 more rows
      # i 1 more variable: iter <dbl>

---

    Code
      summary_list$global_summary
    Output
      # A tibble: 1 x 19
        group null_degree_mean null_triangles_mean null_clust.coef_mean null_dist_mean
        <chr>            <dbl>               <dbl>                <dbl>          <dbl>
      1 Cont~             2.95                 0.7               0.0672           1.76
      # i 14 more variables: null_efficiency_mean <dbl>, null_btw_mean <dbl>,
      #   null_degree_sd <dbl>, null_triangles_sd <dbl>, null_clust.coef_sd <dbl>,
      #   null_dist_sd <dbl>, null_efficiency_sd <dbl>, null_btw_sd <dbl>,
      #   null_degree_sem <dbl>, null_triangles_sem <dbl>, null_clust.coef_sem <dbl>,
      #   null_dist_sem <dbl>, null_efficiency_sem <dbl>, null_btw_sem <dbl>

