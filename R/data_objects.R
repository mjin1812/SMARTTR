##                Constructor for mouse object
#_________________________________________________________________________

## Hidden constructor for mouse
new_mouse <- function(data = list(slices = list()
                                  ),
                      info = list(mouse_ID = 'set ID',
                                  experiment = 'create experiment name',
                                  group = 'name experiment group',
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
#' ## Attributes are strings with initialized values listed below:
#' 1. mouse_ID = 'set ID'
#' 2. experiment = 'create experiment name'
#' 3. group = 'name experiment group',
#' 4. cohort = 'name mouse cohort'
#' 5. input_path = 'set input path'
#' 6. output_path = 'set output path'
#'
#' @usage mouse_example <- mouse() # initializes a mouse object
#' @returns A mouse, a colloquial term for an object of class 'mouse'. A 'mouse' object
#' is also a list, with class list.
#' @export

mouse <- function(... ){
  validate_mouse(new_mouse(...))
}


##             Constructor for slice object
#______________________________________________________________________


new_slice <- function(data = list(registration_obj = NULL,           #list per section
                                  raw_segmentation_data = list(),    #list of length 3(per channel)
                                  segmentation_filters = list(),     #Each element containing sublist the length of registration_obj
                                  segmentation_obj = list(),
                                  raw_forward_warped_data = list(),
                                  forward_warped_data = list(),      # After cleaning up with excluded regions
                                  areas = list()),
                      info = list(slice_ID = NA,
                                  coordinate = -1,
                                  atlas_plate = NA,
                                  conversion_factor = 1.0833,
                                  bin = 1,
                                  z_width = 24,            #Measured in um
                                  hemisphere = NULL,
                                  channels = c('eyfp', 'cfos', 'colabel'),
                                  registration_path = 'set registration image path',
                                  segmentation_path = 'set segmentation image path',    # Segmentation path may not be used
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
#' ## Attributes is a list with named initialized values listed below:
#' 1. slice_ID = NA,
#' 2. coordinate = -1,
#' 3. atlas_plate = NA,
#' 4. bin = 1,
#' 5. z_width = 24,            #Measured in um
#' 6. hemisphere = NULL,
#' 7. channels = c('eyfp', 'cfos', 'colabel'),
#' 8. registration_path = 'set registration image path',
#' 9. segmentation_path = 'set segmentation image path',    # Segmentation path may not be used
#' 10. regions_excluded = c("layer 1","PIR1","TR1","PAA1","NLOT1","OT1","MOBgl","OV","VLPO","SO",
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

