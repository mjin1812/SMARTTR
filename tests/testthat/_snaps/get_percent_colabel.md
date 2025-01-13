# You can analyze subsets of rois

    Code
      e$colabel_percent$eyfp$individual %>% head()
    Output
      # A tibble: 6 x 7
      # Groups:   mouse_ID, group [1]
        mouse_ID group acronym name        count.colabel count.eyfp colabel_percentage
           <dbl> <chr> <chr>   <chr>               <int>      <int>              <dbl>
      1      669 Shock dCA1    Field CA1 ~            15        307               4.89
      2      669 Shock dCA3    Field CA3 ~            10        156               6.41
      3      669 Shock dDG     Dentate gy~             9        173               5.20
      4      669 Shock ENT     Entorhinal~            45        392              11.5 
      5      669 Shock SUB     Subiculum              16        236               6.78
      6      669 Shock vCA1    Field CA1 ~            37        368              10.1 

