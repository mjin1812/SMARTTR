#' SMARTR: Wholebrain analysis package in development for the Denny Lab
#'
#' @description
#' The base of this pipeline is an R package in development with the working title, `SMARTR`, a self-referential play on
#' a previous package developed as an extension to `wholebrain` called `SMART`.  This package allows for the
#' user-friendly pre-processing of segmentation data generated from ImageJ to a be compatible with the [`wholebrain`] package to generate
#' region-based cell counts that are normalized by volume.
#' It will also provides tools for data analysis based on experimental groupings.
#'
#' @details The package currently allows for easy implementation of the following:
#' 1) Setting up the pipeline by specifying experimentparameters, and save directories.
#' 2) The interactive registration process.
#' 3) Importing raw segmentation data from .txt files generated from ImageJ for multiple channels.
#' 4) Optionally creating a filter for the 'cfos' and 'eyfp' channels to clean segmented counts.
#' 5) Creating a segmentation object that is compatible with `wholebrain` functions.
#' 6) Forward warping and mapping the data onto the standardized mouse atlas.
#' 7) Cleaning the mapped data in all the following ways: + Removing cells that map outside the boundaries of the atlas.
#'    + Omitting regions by a default list of regions to omit.
#'    + Omitting regions by user specified region acronyms.
#'    + Removing Layer 1 cells
#'    + Removing cells from a contralateral hemisphere per slice if the registrations are divided by right and left hemispheres.
#' 8)  Obtaining cell counts normalized by region volume (per mm^2^) and region areas (per mm^2^).
#' @md
#' @docType package
#' @name SMARTR
NULL

