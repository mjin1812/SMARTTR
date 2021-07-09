##                Constructor for mouse object
#_________________________________________________________________________

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
slice <- function(...){
  validate_slice(new_slice(...))
}


##             Constructor for experiment object
#______________________________________________________________________

