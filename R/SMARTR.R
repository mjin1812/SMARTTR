#' SMARTR: A mapping, analysis, and visualization package for wholebrain dual-ensemble coronal datasets.
#'
#' @description
#' The base of this pipeline is an R package in with the title, `SMARTR`, a self-referential play on
#' a previous package developed as an extension to `wholebrain` called [`SMART`].  This package allows for the
#' user-friendly pre-processing of segmentation data generated from ImageJ to a be compatible with the `wholebrain` package to generate
#' region-based cell counts that are normalized by volume. It will also provides tools for data analysis based on experimental groupings.
#'
#' @details  Object descriptions
#'
#' The data for analysis will be stored in objects that allow for more neat bundling of useful information together.
#'
#' A [slice] object will contain all the data related to registration, segmentation for each channel, and cell counts for a particular image.
#' It will also contain “metadata” about your experimental images, such as what the experimenter-assigned slice ID is,
#' which brain atlas AP coordinate matches best with the given image, and what the path to the image used for registration is.
#' These metadata are stored as the object’s attributes.
#'
#' A [mouse] object is an object that will store multiple slice objects (and therefore all the information in it),
#' and will eventually store the combined cell data and the region cell counts normalized by volume.
#' Like a slice, it will also contain “metadata” about your mouse stored as attributes.





#' An [experiment]An experiment object consists of a list of processed mouse objects with raw data from slices omitted, and experimental attributes stored as a list.
#' It will also contain “metadata” about your experimental personnel and analysis groups stored as attributes.
#'
#' @section The package currently allows for easy implementation of the following steps:
#'
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
#' @docType _PACKAGE
#' @name SMARTR
NULL

