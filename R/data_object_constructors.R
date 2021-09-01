##                Constructor for mouse object
#_________________________________________________________________________

## Hidden constructor for mouse
new_mouse <- function(data = list(slices = list(),
                                  cell_table = NULL,
                                  normalized_counts = NULL),
                      info = list(mouse_ID = 'set ID',
                                  sex = 'female',
                                  age = NULL,
                                  cre_genotype = NULL,
                                  reporter = NULL,
                                  strain = NULL,
                                  experiment = 'create experiment name',
                                  group = 'name experiment group: e.g. control or AD',
                                  drug = NULL,
                                  cohort = 'name mouse cohort',
                                  input_path = 'set input path',
                                  output_path = 'set output path')
){
  stopifnot(is.list(data))
  structure(data,
            class = 'mouse',
            info = info)
}

## Validator (placeholder)
validate_mouse <- function(m){
  data <- unclass(m)
  info <- attributes(m)
  info_names <- names(info$info)

  # Check if the user specified any special arguments into their info path
  valid_attributes <- c("mouse_ID",
                        "sex",
                        "age",
                        "cre_genotype",
                        "reporter",
                        "strain",
                        "experiment",
                        "group",
                        "drug",
                        "cohort",
                        "input_path",
                        "output_path")

  if (!all(info_names %in% valid_attributes)){
    message(paste0("Specified a custom attribute name. The attribute will be stored in the object but it can only be used for documentation purposes.",
                   "\nPlease check if the custom attribute may apply to one of the standard attribute names."))
  }


  return(m)
}

## Helper constructor
#' Create a mouse object
#'
#' @param mouse_ID (str, default = 'set ID') e.g. '1_1'
#' @param sex (str, default = "female")
#' @param age (str, default = NULL)
#' @param cre_genotype (str, default = NULL)
#' @param reporter (str, default = NULL)
#' @param strain (str, default = NULL) e.g. 'B6'
#' @param experiment (str, default = NULL) e.g. 'sundowning'
#' @param group (str, default = NULL) e.g. 'control' or 'AD'
#' @param drug (str, default = NULL) e.g. 'vehicle' or 'ketamine'
#' @param cohort (str, default = NULL)
#' @param input_path (str, default = 'set input path') Path to the mouse data TODO: currently not used. Will include search function to parse out files and individual slice information for a mouse.
#' @param output_path (str, default = 'set output path') Set the path to the folder you want to save the mouse RDATA file to.
#' @param ... additional custom keyword pair attributes you'd like to store
#'
#' @description mouse() constructs an S3 object of class 'mouse'. A mouse object
#' consists of a list of slice objects and attributes stored as a list.
#'
#' The slice objects are added to a mouse object with the function [SMARTR::add_slice()].
#' Each slice is a named element in the mouse object list, with the naming convention dependent on the slice ID and hemisphere attributes of the slice object.
#'
#' If you are processing either a left or right hemisphere, the slice is named with the convention: "slice_ID" appended with its "hemisphere"
#' If the hemisphere attribute is NULL,  i.e. if the whole slice aligns well with a single atlas plate and there is no need to create separate slice objects per hemisphere,
#' then the slice is named with the convention: "slice_ID"
#'
#' @details
#' The mouse attributes can be assigned as arguments to the mouse constructor function. See the parameters listed for the default values for these attributes
#' Note that you are able to add custom attributes as keyword pairs, if you would like to keep track of an additional piece of information. However, this will only serve a descriptive purpose
#' and will not be used for analysis. You may not need to use all mouse attributes but fill out as many are applicable to your experiment.
#'
#' @usage mouse_example <- mouse() # initializes a mouse object
#' @returns A mouse, a colloquial term for an object of class 'mouse'. A 'mouse' object
#' is also a list, with class list.
#' @seealso
#' See also [SMARTR::slice()] for the description of a slice object and it's attributes.
#' See also [SMARTR::add_slice()] for the description of how to add a slice object to a mouse object.
#' @export

mouse <- function(mouse_ID = 'set ID',
                  sex = 'female',
                  age = NULL,
                  cre_genotype = NULL,
                  reporter = NULL,
                  strain = NULL,
                  experiment = NULL,
                  group = NULL,
                  drug = NULL,
                  cohort = NULL,
                  input_path = 'set input path',
                  output_path = 'set output path',
                  ... ){


  info = list(mouse_ID = mouse_ID,
              sex = sex ,
              age = age,
              cre_genotype = cre_genotype,
              reporter = reporter,
              strain = strain,
              experiment = experiment,
              group = group,
              drug = drug,
              cohort = cohort,
              input_path = input_path,
              output_path = output_path,
              ...)

  validate_mouse(new_mouse(info = info))
}


##             Constructor for slice object
#______________________________________________________________________


new_slice <- function(data = list(registration_obj = NULL,           #list per section
                                  raw_segmentation_data = NULL,    #list of length 3(per channel)
                                  segmentation_filters = NULL,     #Each element containing sublist the length of registration_obj
                                  segmentation_obj = NULL,
                                  raw_forward_warped_data = NULL,
                                  forward_warped_data = NULL,      # After cleaning up with excluded regions
                                  areas = NULL),
                      info = list(slice_ID = NA,
                                  coordinate = -1,
                                  atlas_plate = NA,
                                  conversion_factor = 1.0833,   # From pixels to microns
                                  bin = 1,
                                  z_width = 24,            #Measured in um
                                  hemisphere = NULL,
                                  channels = c('eyfp', 'cfos', 'colabel'),
                                  registration_path = 'set registration image path',
                                  # segmentation_path = 'set segmentation image path',    # Segmentation path may not be used
                                  slice_directory = NULL,
                                  regions_excluded = c("layer 1","PIR1","TR1","PAA1","NLOT1","OT1","MOBgl","OV","VLPO","SO",
                                    "BA","TU","MEAav","ME","TMv","PVp","SUMl","SCzo","fiber tracts","VS")

)
){
  stopifnot(is.list(data))
  structure(data,
            class = 'slice',
            info = info)
}

## Validator (placeholder)
validate_slice <- function(s){

  data <- unclass(s)
  info <- attributes(s)
  info_names <- names(info$info)

  # Check if the user specified any special arguments into their info path
  valid_attributes <- c("slice_ID",
                        "coordinate",
                        "atlas_plate",
                        "conversion_factor",
                        "bin",
                        "z_width",
                        "hemisphere",
                        "channels",
                        "registration_path",
                        "slice_directory",
                        "regions_excluded")

  if (!all(info_names %in% valid_attributes)){

    print(info_names %in% valid_attributes)


    print(info_names)
    print(valid_attributes)

    message(paste0("Specified a custom attribute name. The attribute will be stored in the object but it can only be used for documentation purposes.",


                                  "\nPlease check if the custom attribute may apply to one of the standard attribute names."))
  }

  return(s)
}

## Helper constructor
#' Create a slice object
#'
#' @param slice_ID (str, default = NA) Slice name
#' @param coordinate (num, default = -1) Allen mouse brain atlas coordinate aligning with slice.
#' @param atlas_plate (int, default = NA) Atlas place number. TODO: Currently unused
#' @param conversion_factor (num, 1.0833) pixel-to-micron conversion factor
#' @param bin (int, default = 1) Whether the registration image was binned in ImageJ.
#' @param z_width (num, default = 24) The z-stack width in um.
#' @param hemisphere (str, default = NULL) Hemisphere to process. "left", "right" or NULL is legal.
#' @param channels (str, default = c("cfos", "eyfp", "colabel")) The channels to process.
#' @param registration_path (str, default = 'set registration image path') TODO: deprecate this in favor of slice_directory
#' @param slice_directory (str, default = NULL) The directory where slice information is stored. TODO: Change the import and registration functions to only rely on this path.
#' @param regions_excluded (str, default = ("layer 1","PIR1","TR1","PAA1","NLOT1","OT1","MOBgl","OV","VLPO","SO",
#' "BA","TU","MEAav","ME","TMv","PVp","SUMl","SCzo","fiber tracts","VS")) The list of acronyms corresponding to regions to exclude for this slice.
#' @param ... additional custom keyword pair attributes you'd like to store
#'
#' @description slice() constructs an S3 object of class 'slice'. A slice object
#' consists of a list of lists storing information about registration, segmentation, raw mapped data and cleaned mapped data.
#' The object attributes are also stored as a list.
#'
#' @details
#' The slice attributes can be assigned as arguments to the slice constructor function. See the parameters listed for the default values for these attributes
#' Note that you are able to add custom attributes as keyword pairs, if you would like to keep track of an additional piece of information. However, this will only serve a descriptive purpose
#' and will not be used for analysis. You may not need to use all slice attributes but fill out as many are applicable to your experiment.
#'
#' @usage slice_example <- slice() # initializes a slice object
#' @returns A slice, a colloquial term for an object of class 'slice'. A 'slice' object
#' is also a list, with class list.
#' @export


slice <- function(slice_ID = NA,
                  coordinate = -1,
                  atlas_plate = NA,
                  conversion_factor = 1.0833,   # From pixels to microns
                  bin = 1,
                  z_width = 24,            #Measured in um
                  hemisphere = NULL,
                  channels = c('eyfp', 'cfos', 'colabel'),
                  registration_path = 'set registration image path',
                  # segmentation_path = 'set segmentation image path',    # Segmentation path may not be used
                  slice_directory = NULL,
                  regions_excluded = c("layer 1","PIR1","TR1","PAA1","NLOT1","OT1","MOBgl","OV","VLPO","SO",
                                       "BA","TU","MEAav","ME","TMv","PVp","SUMl","SCzo","fiber tracts","VS"),
                  ...){


  info <- list(slice_ID = slice_ID,
               coordinate = coordinate,
               atlas_plate = atlas_plate,
               conversion_factor = conversion_factor,   # From pixels to microns
               bin = bin,
               z_width = z_width,            #Measured in um
               hemisphere = hemisphere,
               channels = channels,
               registration_path = registration_path,
               # segmentation_path = segmentation_path,    # Segmentation path may not be used
               slice_directory = slice_directory,
               regions_excluded = regions_excluded,
               ...)


  validate_slice(new_slice(info = info))
}


##             Constructor for experiment object
#______________________________________________________________________

## Hidden constructor for experiment
new_experiment <- function(data = list(mice = list()),
                           info = list(experiment_name = NULL,
                                       experimenters = NULL,
                                       channels = NULL,
                                       experiment_groups = NULL,
                                       drug_groups = NULL,
                                       sex_groups = NULL,
                                       cohorts = NULL,
                                       strains = NULL,
                                       genotypes = NULL,
                                       reporters = NULL,
                                       ages = NULL,
                                       output_path = 'set output path for your experiment'
                                       )

){
  stopifnot(is.list(data))
  stopifnot(is.list(info))
  structure(data,
            class = 'wb_experiment',
            info = info)
}

## Experiment Validator (placeholder)
validate_experiment <- function(e){
  data <- unclass(e)
  info <- attributes(e)
  info_names <- names(info$info)

  # Check if the user specified any special arguments into their info path
  valid_attributes <- c("experiment_name",
                        "experimenters",
                        "channels",
                        "experiment_groups",
                        "drug_groups",
                        "sex_groups",
                        "cohorts",
                        "strains",
                        "genotypes",
                        "reporters",
                        "ages",
                        "output_path")

  if (!all(info_names %in% valid_attributes)){
    message(paste0("Specified a custom attribute name. The attribute will be stored, but it can only be used for documentation purposes.",
                   "\nPlease check if the custom attribute may apply to one of the standard attribute names."))
  }

  return(e)
}


## Helper constructor
#' Create an experiment object
#' @param experiment_name (str, default = NULL)
#' @param experimenters (str, default = NULL)
#' @param channels (str, default = NULL) Autogenerated with [add_mouse()] function. Will detect all unique channels stored in a mouse object.
#' @param experiment_groups (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param drug_groups (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param sex_groups (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param cohorts (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param strains (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param genotypes (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param reporters (str, default = NULL) Autogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param ages (str, default = NULL). AUtogenerated with [add_mouse()] function. Must exactly match the string values from the mouse objects.
#' @param output_path (str, default = 'set output path for your experiment') Where to save the RData file for your experiment object
#' @param ... additional custom keyword pair attributes you'd like to store
#'
#' @description experiment() constructs an S3 object of class 'wb_experiment'. An experiment object
#' consists of a list of processed mouse objects with raw data from slices omitted, and experimental attributes stored as a list.
#' @details
#'
#' The experimental attributes can be assigned as arguments to the experiment constructor function. See the parameters listed for the default values for these attributes
#' Note that you are able to add custom attributes as keyword pairs, if you would like to keep track of an additional piece of information. However, this will only serve a descriptive purpose
#' and will not be used for analysis. You may not need to use all experimental attributes but fill out as many are applicable to your experiment.
#'
#' @usage my_experiment <- experiment() # construct an experiment object
#' @returns An experiment, a colloquial term for an object of class 'wb_experiment'. An 'experiment' object
#' is also a list, with class list.
#' @seealso
#' See also [SMARTR::mouse()] for the description of a mouse object and it's attributes.
#' @export

experiment <- function(experiment_name = NULL,
                       experimenters = NULL,
                       channels = NULL,
                       experiment_groups = NULL,
                       drug_groups = NULL,
                       sex_groups = NULL,
                       cohorts = NULL,
                       strains = NULL,
                       genotypes = NULL,
                       reporters = NULL,
                       ages = NULL,
                       output_path = 'set output path for your experiment',
                       ...){


  info <- list(experiment_name = experiment_name,
               experimenters = experimenters,
               channels = channels,
               experiment_groups = experiment_groups,
               drug_groups = drug_groups,
               sex_groups = sex_groups,
               cohorts = cohorts,
               strains = strains,
               genotypes = genotypes,
               reporters = reporters,
               ages = ages,
               output_path = output_path,
               ...)

  validate_experiment(new_experiment(info = info))
}















































