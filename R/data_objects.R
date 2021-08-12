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
                                  drug = 'drug group: e.g. vehicle or ketamine',
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
  # values <- unclass(m)
  # m_attrib <- attributes(m)
  return(m)
}

## Helper (placeholder)
#' Create a mouse object
#'
#' @description mouse() constructs an S3 object of class 'mouse'. A mouse object
#' consists of a list of slice objects and attributes stored as a list.
#'
#' @details
#' ## The list of slice object.
#' The slice objects are added to a mouse object with the function [SMARTR::add_slice()].
#' Each slice is a named element in the mouse object list, with the naming convention dependent on the slice ID and hemisphere attributes of the slice object.
#'
#' If you are processing either a left or right hemisphere, the slice is named with the convention: "slice_ID" appended with its "hemisphere"
#' If the hemisphere attribute is NULL,  i.e. if the whole slice aligns well with a single atlas plate and there is no need to create separate slice objects per hemisphere,
#' then the slice is named with the convention: "slice_ID"
#'
#' ## Attributes are strings with initialized (default) values listed below.
#' 1. mouse_ID = 'set ID'
#' 2. sex = "female"
#' 3. age = NULL
#' 4. cre_genotype = NULL
#' 5. reporter = NULL
#' 6. strain = NULL
#' 7. experiment = 'create experiment name'
#' 8. group = 'name experiment group: e.g. control or AD'
#' 9. drug = 'drug group: e.g. vehicle or ketamine',
#' 10. cohort = 'name mouse cohort'
#' 11. input_path = 'set input path'
#' 12. output_path = 'set output path'
#'
#' Note that you may not need to use all of these attributes but fill out as many are applicable to your experiment.
#'
#' @usage mouse_example <- mouse() # initializes a mouse object
#' @returns A mouse, a colloquial term for an object of class 'mouse'. A 'mouse' object
#' is also a list, with class list.
#' @seealso
#' See also [SMARTR::slice()] for the description of a slice object and it's attributes.
#' See also [SMARTR::add_slice()] for the description of how to add a slice object to a mouse object.
#' @export

mouse <- function(... ){
  validate_mouse(new_mouse(...))
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
                                  segmentation_path = 'set segmentation image path',    # Segmentation path may not be used
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
  # values <- unclass(s)
  # m_attrib <- attributes(s)
  return(s)
}

## Helper (placeholder)
#' Create a slice object
#'
#' @description slice() constructs an S3 object of class 'slice'. A slice object
#' consists of a list of lists storing information about registration, segmentation, raw mapped data and cleaned mapped data.
#' The object attributes are also stored as a list.
#'
#' ## Attributes is a list with named initialized (default) values listed below:
#' 1. slice_ID = NA,
#' 2. coordinate = -1,
#' 3. atlas_plate = NA,
#' 4. bin = 1,
#' 5. z_width = 24,            #Measured in um
#' 6. hemisphere = NULL,
#' 7. channels = c('eyfp', 'cfos', 'colabel'),
#' 8. registration_path = 'set registration image path',
#' 9. segmentation_path = 'set segmentation image path',    # path to images used for wholebrain segmentation method (currently unused)
#' 10. slice_directory = NULL                               # Slice directory for importing segmentation data if registration path is not assigned.
#' 11. regions_excluded = c("layer 1","PIR1","TR1","PAA1","NLOT1","OT1","MOBgl","OV","VLPO","SO",
#' "BA","TU","MEAav","ME","TMv","PVp","SUMl","SCzo","fiber tracts","VS")
#'
#' @usage slice_example <- slice() # initializes a slice object
#' @returns A slice, a colloquial term for an object of class 'slice'. A 'slice' object
#' is also a list, with class list.
#' @export


slice <- function(...){
  validate_slice(new_slice(...))
}


##             Constructor for experiment object
#______________________________________________________________________


# Attributes:
# Experiment_groups (Autogenerate)
# drug_groups (autogenerate)
# Output path
#
#













