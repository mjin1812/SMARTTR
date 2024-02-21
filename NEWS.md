# SMARTR 0.0.0.9000
* First release of SMARTR

# SMARTR 1.0.0
* Major change in the way area/volume of regions is calculated to account for bug in underlying atlas in wholebrain that pulls the wrong region registration outlines for certain subregions. Specifically, areas such as the DG-sg and DG-mo are represented twice in the same atlas plate when the hippocampus curves rostral-caudally. Previously only one of the two subregion outlines were used to calculate area, which should now be fixed.
* Change in the `exclude_anatomy`, `get_registered_volumes`, and `normalize_cell_counts` functions to include a `simplify_keywords` parameter which automatically re categorizes subregions to their major parent regions based on specific keywords. 
* Modified the `normalize_cell_counts` and `combine_cell_counts` function to also store counts and areas/volume by slice in the mouse and experiment objects, respectively. This allows for easy evaluation of the raw mapped dataset.
