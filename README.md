
# SMARTR: a mapping, analysis, and visualization package for dual-ensemble datasets across the mouse brain


## Welcome to the SMARTR package!

SMARTR is a package designed for the high-throughput mapping and analysis of dual-labelled ensemble datasets. Genetic tagging strategies such as the ArcCreERT^2^ mouse line allow for investigation of neural ensembles underlying behavior across multiple time points.

![*Tagging two ensembles using the ArcCreERT^2^ mouse line*](man/figures/1.Pipeline_tagging_schematic.jpg)


However, few tools exist to map dual-labelled ensembles, as well as their overlap across the brain to a standardized atlas space. The SMARTR package facilitates this by bridging an optimized [approach](https://osf.io/ynqp7/) for cell-segmentation and colocalization that was developed in FIJI/ImageJ with registration and mapping capabilities from the [wholebrain](https://github.com/tractatus/wholebrain) and [SMART](https://github.com/mjin1812/SMART) packages. 

![*Pipeline overview*](man/figures/2.general_pipeline_schematic.jpg)

Moreover, SMARTR provides a streamlined API for storing metadata related to imaging, animal subject, and experimental parameters and groupings. Data management is intrinsically built into the package, which helps automate the process of combining ensemble datasets across multiple images, animals, and experimental groupings.

Finally, SMARTR provides a set of built-in analysis and visualization functions to conduct network analysis for each ensemble dataset for immediate downstream analysis. SMARTR is designed to facilitate dual-ensemble brain mapping projects by lowering the technical barrier for registration, segmentation, and statistical analysis.

Check out the [Get Started](./articles/SMARTR.html) page to see a more thorough breakdown of the pipeline and get started using an example dataset! 


## Contact
Answers to common questions and issues may be found in the FAQ page. If questions are not otherwise answered through the FAQ or through the tutorial, feel free to contact me at 

