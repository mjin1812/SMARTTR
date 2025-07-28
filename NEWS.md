# SMARTTR 1.0.1
* To adhere to CRAN naming policy, we are changing our package repository name to SMARTTR. The new package website can be assessed [here](https://mjin1812.github.io/SMARTTR/). Listed below is the previous changelog of SMARTR. 

## SMARTR 0.0.0.9000
* First release of SMARTR

## SMARTR 1.0.0
* Major change in the way area/volume of regions is calculated to account for bug in underlying atlas in wholebrain that pulls the wrong region registration outlines for certain subregions. Specifically, areas such as the DG-sg and DG-mo are represented twice in the same atlas plate when the hippocampus curves rostral-caudally. Previously only one of the two subregion outlines were used to calculate area, which should now be fixed.
* Change in the `exclude_anatomy`, `get_registered_volumes`, and `normalize_cell_counts` functions to include a `simplify_keywords` parameter which automatically re categorizes subregions to their major parent regions based on specific keywords. 
* Modified the `normalize_cell_counts` and `combine_cell_counts` function to also store counts and areas/volume by slice in the mouse and experiment objects, respectively. This allows for easy evaluation of the raw mapped dataset.
* Our various analysis and visualization function are now compatible with the  Kim lab's [unified atlas](https://kimlab.io/brain-map/atlas/). The unified atlas synthesizes label from two of the most highly used atlases --Franklin-Paxinos (FP) and the common coordinate framework (CCF) from the Allen Institute.Check out the various help pages for analysis and plotting functions to see which ones now include`ontology` as a relevant parameter. Simply change the value from `allen` to `unified` to process datasets using this ontology.
* New function `simplify_cell_count` can process external datasets and simplify regions based on keywords. Simplification is based on user indicated atlas (`allen` or `unified`)

## SMARTR 1.0.1
* New function `export_permutation_results` reformats permutation results and saves it nicely for users.
* Tables are now saved into an autocreated `tables` folder in the output directory for organization.
* Added `create_joined_networks` and `plot_joined_networks` as as capability to visualize the overlapping edges in two networks with outer joined nodes. The capability to export a list of overlapping edges between the joined network is included as a parameter in `create_joined_networks`. See the documentation for more details.
* Altered `plot_correlation_heatmaps() function to automatically plot parent regions as subfacets`

## SMARTTR 1.0.2
* Added `rewire_network()` function to implement maslov sneppen rewiring to generate null networks. Added `summarize_null_networks()` function to summarize the output of multiple rewiring iterations.
* Added patch for bug in the `create_joined_networks()` function
* Added filter isolates option for networks. This allows isolated nodes to be retained for global calculations and retains the same node numbers between networks.
* Added `import_mapped_datasets()` as a way to import an externally mapped dataset and `check_ontology_coding()` to check the ontology coding of an imported dataset. A tutorial on how to use these function are documented on our article titled `Importing External Datasets`.
* Added new progress bar display option for the permutation analysis function `correlation_diff_permutation()`
* Added new feature to plot a permutation difference network using `plot_perm_diff_network()`. 
* Added new heatmap plotting options allow for viridis color scheme and new custom faceting colors for parent regions using ggh4x
* Bug fixes in `plot_normalized_counts()`, `volcano_plot()`, `normalize_cell_counts()`
