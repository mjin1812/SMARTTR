
#' @title Add slice to a mouse object
#' @usage m <- add_slice(m, s, replace = FALSE)
#' @param m mouse object
#' @param s slice object
#' @param replace (bool, default = FALSE) Replace a slice already contained in a mouse object.
#' @export

add_slice <- function(m, s, replace = FALSE){

  # Get slice ID & hemisphere
  slice_ID <- attr(s, 'info')$slice_ID
  hemisphere <- attr(s, 'info')$hemisphere

  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }

  # First slice
  if (length(m$slices) < 1){

    m$slices[[1]] <- s
    names(m$slices) <- slice_name

  } else{

    # check in list of previous stored slice names for a match
    stored_names <- names(m$slices)

    # If there is a previous slice that matches, give user warning to set replace to TRUE to overwrite
    match <- FALSE

    for (stored_name in stored_names){

      if (identical(stored_name, slice_name)){

        match <- slice_name
        message(paste0('There was existing data found for slice ', slice_ID,
                       ", ", hemisphere, ' hemisphere\n'))

        if (replace){
          # replace slice
          m$slices[[slice_name]] <- s
          message(paste0('Replaced existing slice data!'))
        } else {
          stop(paste0('If you want to replace a previous slice object,',
                      'then set the "replace" argument to "TRUE".'))
        }
      }
    }

    # If there were no matches, store the slice as a new slice

    if (isFALSE(match)){
      index <- length(m$slices) + 1
      m$slices[[index]] <- s
      names(m$slices)[index] <- slice_name
    }
  }
  return(m)
}



###______________Generic functions ________________###

## Creating generic function for registration
#' Register
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
register <- function(x, ...) {
  UseMethod("register")
}

## Creating generic function for segmentation importation
#' import_segmentation
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
import_segmentation <- function(x, ...){
  UseMethod('import_segmentation')
}

## Creating generic function for filtering raw segmentation data
make_segmentation_filter <- function(x, ...){
  UseMethod('make_segmentation_filter')
}

## Create segmentation object
#' make_segmentation_object
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
make_segmentation_object <- function(x, ...){
  UseMethod('make_segmentation_object')
}


## forward_warp segmentation data to atlas space
#' map_cells_to_atlas
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
map_cells_to_atlas <- function(x, ...){
  UseMethod('map_cells_to_atlas')
}

## excluded designated regions from the mapped cell data table
#' exclude_anatomy
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
exclude_anatomy <- function(x, ...){
  UseMethod('exclude_anatomy')
}

## excluded designated regions from the mapped cell data table
#' get_registered_areas
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
get_registered_areas <- function(x, ...){
  UseMethod('get_registered_areas')
}


# #_________ TESTING NEW SEGMENTATION METHOD _________
#
# ## Creating generic function for segmentation importation
# import_segmentation_weka <- function(x, ...){
#   UseMethod('import_segmentation_weka')
# }

#__________________ Mouse & Slice object methods for generic functions __________________________
##____ Creating method for registration of slice object____
#' Register a slice in a slice object
#'
#' @param s slice object
#' @param filter (list)
#' @param ... additional parameters to pass to the SMART::registration2() function, besides 'input', 'coordinate', 'filter' & 'correspondance'
#'
#' @return s slice object
#' @export
#'
#' @examples s <- register(s)
register.slice <- function(s,
                           filter = NULL,
                           ...){
  # Get slice information
  info <- attributes(s)

  # If the user did not set a filter and there is one in the object,
  # set to the filter that is stored
  if (is.null(filter) && !is.null(s$raw_segmentation_obj$filter)){
    filter <- s$raw_segmentation_obj$filter
  }

  # Check if registration already exists
  if (!is.null(s$registration_obj)){
    message(paste0("\nRegistration data for this slice already exists.\nThe existing registration data will be used."))
  }


  # Run the interactive registration correction loop
  quartz()
  s$registration_obj <-  SMART::registration2(input = info$info$registration_path,
                                              coordinate = info$info$coordinate,
                                              filter = filter,
                                              correspondance = s$registration_obj,
                                              ...)
  return(s)
}

##____ Creating method for registration of mouse object____
#' Register a slice in a mouse object
#'
#' @param m  mouse object
#' @param slice_ID ID of slice (str)
#' @param hemisphere 'left', 'right' (str) or NULL if both hemispheres are included
#' @param filter (list)
#' @param replace (bool, default = FALSE) Replace a registration already contained in a mouse object by resetting to NULL value before registration improvement loop.
#' @param ... additional parameters to pass to the SMART::registration2() function, besides 'input', 'coordinate', 'filter' & 'correspondance'
#'
#' @return m  mouse object
#' @export
#'
#' @examples m <- register(m, slice_ID = '1_10', hemisphere = "left", filter = my_filter)
#' TODO: add popup window option to better visualize the internal structures of the images
register.mouse <- function(m,
                           slice_ID = '1_10',
                           hemisphere = NULL,  # set hemisphere to 'left' or 'right' to specify which slice should be pulled
                           filter = NULL,
                           replace = FALSE,
                           ...){


  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }


  # check in list of previous stored slice names for a match
  stored_names <- names(m$slices)

  # Check if there is a previous slice that matches
  match <- FALSE

  for (stored_name in stored_names){
    if (identical(stored_name, slice_name)){
      match <- slice_name

      # Check if registration data already exists
      if (!is.null(m$slices[[match]]$registration_obj)){

        if (replace){
          # replace the slice
          m$slices[[match]]$registration_obj <- NULL

          # Register from scratch
          m$slices[[match]] <- register(m$slices[[match]],
                                           filter = filter, ...)
        } else{
          message("Improving previous registration for this slice. Set replace = TRUE to start from scratch!...\n")

          if (is.null(filter)){
            filter <- SMARTR::segmentation.object$filter
            filter$resize <-  m$slices[[match]]$registration_obj$resize
          }

          # Register matched slices with existing data
          m$slices[[match]] <- register(m$slices[[match]],
                                        filter = filter, ...)
        }
      } else  {
        message("Starting a new registration!")

        # Register matched slices with no registration data yet
        m$slices[[match]] <- register(m$slices[[match]],
                                      filter = filter, ...)

      }
    }
  }
  if (isFALSE(match)){
    message(paste0("There were no slices matching the name ", slice_name, " found in your mouse!" ))
  }
  return(m)
}



##____ Creating method for importing segmentation data for a slice object____

#' Import segmentation data
#' @description Method for importing segmentation data for a slice object
#' @param s : slice object
#' @param mouse_ID : ID of mouse (str)
#' @param channels : channels to import (str vector)
#'
#' @return s slice object
#' @export
#'
#' @examples s <-  import_segmentation(s, mouse_ID = "255", channels = c("cfos", "eyfp", "colabel"))

#' TODO: modify paths to search for correct files regardless of upper or lower case. Use grep and stringi to match multiple patterns

import_segmentation.slice <- function(s,
                                      mouse_ID = '255',
                                      channels = c('cfos', 'eyfp', 'colabel')){

  # Get slice information
  info <- attributes(s)

  ## Use registration image path to find segmentation data path stem
  path_stem <- dirname(info$info$registration_path)
  section <- info$info$slice_ID

  for (k in 1:length(channels)){
    if (tolower(channels[k]) == 'colabel') {
      coloc.table <- read.delim(file.path(path_stem,
                                          paste0( mouse_ID,'_',section,"_Fast_R_cfos_SpotSegmentation_ColocOnly.txt")))
      s$raw_segmentation_data[[k]] <- coloc.table

    } else{
      if (tolower(channels[k]) == 'cfos'){

        # Create import data paths for cfos
        meas_path <- file.path(path_stem, paste0("M_R_cfos_", mouse_ID,'_',section,".txt" ))
        quant_path <- file.path(path_stem, paste0("Q_R_cfos_",mouse_ID,'_',section,"_",'cfos',".txt" ))

      } else if (tolower(channels[k]) == 'eyfp'){

        # Create import data paths for eyfp
        meas_path <- file.path(path_stem, paste0("M_G_eyfp_", mouse_ID,'_',section,".txt" ))
        quant_path <- file.path(path_stem, paste0("Q_G_eyfp_",mouse_ID,'_',section,"_eYFP",".txt" ))
      }

      print(meas_path)
      meas <- read.csv(meas_path)
      print(quant_path)
      quant <- read.csv(quant_path)
      meas$X2_Pix <- meas$CX..pix./(info$info$bin*info$info$conversion_factor) #create position column to account for binning
      meas$Y2_Pix <- meas$CY..pix./(info$info$bin*info$info$conversion_factor) #same as above
      counts <- cbind(meas, quant) #create combined table
      s$raw_segmentation_data[[k]] <- counts
    }
  }
  # Name the elements of the raw data
  names(s$raw_segmentation_data) <- channels
  return(s)
}


##____ Creating method for importing segmentation data for a mouse object____
#' Import raw segmentation data
#' @description Method for importing segmentation data for a mouse object
#' @param  m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere : 'left', 'right' or NULL (both)
#' @param channels : channels to import (str vector)
#' @param replace : replace existing raw segmentation data (bool)
#' @return m mouse object
#' @export

import_segmentation.mouse <- function(m,
                                      slice_ID = '1_10',
                                      hemisphere = NULL,
                                      channels = c('cfos', 'efyp', 'colabel'),
                                      replace = FALSE){

  # Get mouse ID
  mouse_ID <- attr(m, 'info')$mouse_ID

  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }

  # check in list of previous stored slice names for a match
  stored_names <- names(m$slices)

  # Check if there is a previous slice that matches
  match <- FALSE

  for (stored_name in stored_names){
    if (identical(stored_name, slice_name)){
      match <- slice_name

      # Check if segmentation data already exists
      if (!is.null(m$slices[[match]]$raw_segmentation_data)){

        if (replace){
          # replace the slice's raw segmentation data
          m$slices[[match]]$raw_segmentation_data <- NULL

          # Import segmentation from scratch
          m$slices[[match]] <- import_segmentation(m$slices[[match]],
                                                   mouse_ID = mouse_ID,
                                                   channels = channels)
        } else{
          stop(paste0("There is existing segmentation data for this slice! If you want to overwrite it, set replace to TRUE."))

        }
      } else  {
        # Import segmentation with matched slices with no segmentation data yet
        m$slices[[match]] <- import_segmentation(m$slices[[match]],
                                                 mouse_ID = mouse_ID,
                                                 channels = channels)

      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!"))
  }
  return(m)
}

################################ START of work in progress #########################

##____ Method for making an eyfp filter____
# s: slice object
# channels: channels to filter (str vector), can filter both eyfp & cfos
# params:  (list) Has same length as channels. Each element contains a vector of parameters names used for filtering that channel
# ranges" (list of lists) Has same length as channels. Each element contains a list of ranges for the parameters used for filtering.

# TODO: clean up the name of the make filter function internally


make_segmentation_filter.slice <- function(s,
                                      channels = c('eyfp'),
                                      params =list(c("Vol..unit.","Moment1","Moment2",
                                                     "Moment3","Moment4","Sigma")),
                                      ranges = list(list(c(200,12000),c(3,50),
                                                    c(0,600),c(0,2000),c(0,5),c(20,Inf)))){

  for (k in 1:length(channels)){
    if (channels[k] == 'eyfp'){
      data <- s$raw_segmentation_data$eyfp
    } else if (channels[k] == 'cfos'){
      data <- s$raw_segmentation_data$eyfp
    }
    filter <- make.eyfp.filter(data, params = params[[k]], ranges = ranges[[k]])
    s$segmentation_filters[[k]] <- filter

  message('Processed ', channels[k], " channel")
  }
  names(s$segmentation_filters) <- channels
  return(s)
}


##____ Method for creating filter for making eyfp filter for mouse____
# m : mouse object
# slice_ID : ID of slice (str)
# hemisphere: 'left', 'right' or NULL (both)
# channels: channels to import (str vector)
# params:  (list) Has same length as channels. Each element contains a vector of parameters names used for filtering that channel
# ranges" (list of lists) Has same length as channels. Each element contains a list of ranges for the parameters used for filtering.

make_segmentation_filter.mouse <- function(m,
                                           slice_ID = '1_10',
                                           hemisphere = 'left',
                                           channels = c('eyfp'),
                                           params =list(c("Vol..unit.","Moment1","Moment2","Moment3","Moment4","Sigma")),
                                           ranges = list(list(c(200,12000),c(3,50),
                                                              c(0,600),c(0,2000),c(0,5),c(20,Inf)))){
  # Get mouse ID
  mouse_ID <- attr(m, 'info')$mouse_ID

  # Get slice information
  info  <- list()
  index <- vector(mode = 'integer')

  for (k in 1:length(m$slices)){
    if (tolower(attr(m$slices[[k]], 'info')$slice_ID) == tolower(slice_ID)){     # case insensitive mouse id check
      index <- c(index, k)
      info[[length(index)]] <- attr(m$slices[[k]], 'info')
    }
  }

  # If there are multiple slice matches, match by the hemisphere selected
  if (length(index) > 1 && !is.null(hemisphere)){

    # implement check for length of 2?
    if (tolower(info[[1]]$hemisphere) == tolower(hemisphere)){   # If the first match is correct hemisphere, remove other match
      info[[2]] <- NULL
      index <- index[1]
    } else {                       #If the second match is the correct hemisphere
      info[[1]] <- NULL
      index <- index[2]
    }
  } else {
    message(paste0('There were multiple slice objects with the same ID! \n',
                   'Set the hemisphere argument if you want to specify which slice to create a filter for,\n',
                   'otherwise all matched slices will have a filter created.'))
  }


  for (k in 1:length(index)){
    ## add print statement indicating which hemisphere is being registered
    message(paste0('This is slice ',info[[k]]$slice_ID))
    message(paste0('Creating filter(s) for the ',info[[k]]$hemisphere, ' hemisphere'))

    m$slices[[index[k]]] <- make_segmentation_filter(m$slices[[index[k]]],
                                                     channels = channels,
                                                     params = params,
                                                     ranges = ranges)
  }
  return(m)
}


################################ END of work in progress ##############################


##____Method for creating segmentation object for each slice____

#' Create a segmentation object
#' Creates a wholebrain segmentation object and stores it.
#' @description Make a wholebrain compatible segmentation object for a slice in a slice object
#' @param s : slice object
#' @param mouse_ID : ID of mouse (str)
#' @param channels : (str vec) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
#' @param use_filter : TRUE (bool). Use a filter to create more curated segmentation object from the raw segmentation data.
#'
#' @return s slice object
#' @export
#'
#' @examples
#'

# Note: this function needs a lot of cleaning; for now always create a filter for the
# 'eyfp' channel and and set the 'use_filter' variable to TRUE.

## TODO: CLEAN UP so that when use_filter = FALSE, processing 'colabel' channel doesn't throw an error
## TODO: Figure out how Marcos's colocalization filter works
## TODO: Change the colocalization filtering step so that one approach isn't hardcoded into the pipeline

make_segmentation_object.slice <- function(s,
                                           mouse_ID = '255',
                                           channels = c('cfos', 'eyfp', 'colabel'),
                                           use_filter = TRUE){
  # If use_filter is TRUE, check if a filter exists for either the cfos channel or eyfp
  # Get slice information
  info <- attributes(s)

  ## Use registration image path to find segmentation data path stem
  path_stem <- dirname(info$info$registration_path)
  section <- info$info$slice_ID

  for (k in 1:length(channels)){
    ## Check if a filter currently exists for cfos and eyfp channel OR if a filter exists for colabeles
    if (use_filter && !is.null(s$segmentation_filters[[channels[k]]]) |
        use_filter && channels[k] =='colabel'){

      ## special filtering process for colabel to match with eyfp channel by volume
      ## (MARCOS's approach in 'Helper.R', THIS NEEDS TO BE CHANGED SO IT's NOT HARDCODED INTO THE PIPELINE)
      if (channels[k]=='colabel'){
        eyfp_meas_16bit <- read.csv(file.path(path_stem,
                                              paste0("M_", mouse_ID,'_',section,"_Fast_G_eYFP_LabelImage_16bit.txt")))

        # Create a list of the filtered eyfp cells
        eyfp_counts_filtered <- s$raw_segmentation_data[['eyfp']][-s$segmentation_filters[['eyfp']], ]

        coloc_cell_list <- identify.colabeled.cells.filter(s$raw_segmentation_data[[channels[k]]],
                                                          eyfp_counts_filtered$Nb+1,
                                                          eyfp_meas_16bit$Label)
        coloc.counts <- eyfp_counts_filtered[match(coloc_cell_list, eyfp_counts_filtered$Nb+1),]
        seg.coloc <- SMARTR::segmentation.object
        seg.coloc$soma$x <- coloc.counts$X2_Pix
        seg.coloc$soma$y <- coloc.counts$Y2_Pix
        seg.coloc$soma$area <- coloc.counts$Vol..pix.
        seg.coloc$soma$intensity <- coloc.counts$Mean

        # storing segmentation object in  slice
        s$segmentation_obj[[k]] <- seg

      } else{
        # Using filter to create filtered counts
        counts_filtered <- s$raw_segmentation_data[[channels[k]]][-s$segmentation_filters[[channels[k]]], ]
        seg <- SMARTR::segmentation.object
        seg$soma$x <- counts_filtered$X2_Pix
        seg$soma$y <- counts_filtered$Y2_Pix
        seg$soma$area <- counts_filtered$Vol..pix.
        seg$soma$intensity <- counts_filtered$Mean
      }
    } else{
      # create segmentation object
      seg <- SMARTR::segmentation.object
      seg$soma$x <- s$raw_segmentation_data[[channels[k]]]$X2_Pix
      seg$soma$y <- s$raw_segmentation_data[[channels[k]]]$Y2_Pix
      seg$soma$area <- s$raw_segmentation_data[[channels[k]]]$Vol..pix.
      seg$soma$intensity <- s$raw_segmentation_data[[channels[k]]]$Mean

      # storing segmentation object in  slice
      s$segmentation_obj[[k]] <- seg
    }
  }
  names(s$segmentation_obj) <- channels
  return(s)
}

##____Method for creating segmentation object for each mouse____
#' Creates a wholebrain segmentation object and stores it
#' @description Make a wholebrain compatible segmentation object for a slice in a mouse object
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere: 'left', 'right' or NULL (both)
#' @param channels : (vec) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
#' @param replace : replace existing raw segmentation data (bool)
#' @param use_filter: TRUE (bool). Use a filter to create more curated segmentation object from the raw segmentation data.


make_segmentation_object.mouse <- function(m,
                                           slice_ID = '1_10',
                                           hemisphere = NULL,
                                           channels = c('cfos', 'eyfp', 'colabel'),
                                           replace = FALSE,
                                           use_filter = TRUE){

  # Get mouse ID
  mouse_ID <- attr(m, 'info')$mouse_ID

  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }

  # check in list of previous stored slice names for a match
  stored_names <- names(m$slices)

  # Check if there is a previous slice that matches
  match <- FALSE

  for (stored_name in stored_names){
    if (identical(stored_name, slice_name)){
      match <- slice_name

      # Check if segmentation object data already exists
      if (!is.null(m$slices[[match]]$segmentation_obj)){

        if (replace){
          # replace the slice's segmentation object
          m$slices[[match]]$segmentation_obj <- NULL

          # Make segmentation object from scratch
          m$slices[[match]] <- make_segmentation_object(m$slices[[match]],
                                                        mouse_ID = mouse_ID,
                                                        channels = channels,
                                                        use_filter = use_filter)
        } else{
          stop(paste0("There is an existing segmentation object for this slice! If you want to overwrite it, set replace to TRUE."))
        }
      } else  {

        # Import segmentation with matched slices with no segmentation data yet
        m$slices[[match]] <- make_segmentation_object(m$slices[[match]],
                                                      mouse_ID = mouse_ID,
                                                      channels = channels,
                                                      use_filter = use_filter)
      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!" ))
  }
  return(m)
}

##____Method for forward warping the segmentation data per slice____

#' @title Map cells to atlas for slice object
#' @description Method for forward warping segmentation data to atlas space for a slice object.
#' Requires segmentation objects to be made for channels specified and a registration object.
#' @param s : slice object
#' @param channels : channels to filter (str vector), can filter both eyfp & cfos
#' @param clean : TRUE (bool). Remove cells that don't map to any regions.
#' @param display : TRUE (bool). Display the results of the forward warp for the slice.
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the wholebrain::inspect.registration() function
#' @export

map_cells_to_atlas.slice <- function(s,
                                     channels = c('cfos', 'eyfp', 'colabel'),
                                     clean =  TRUE,
                                     display = TRUE,
                                     ...){

  for (k in 1:length(channels)){
    cell_data <- inspect.registration(s$registration_obj,
                                      s$segmentation_obj[[channels[k]]],
                                      forward.warps = TRUE,
                                      device = display,
                                      ...)
    s$raw_forward_warped_data[[k]] <- cell_data

    if(clean){
      cell_data <- cell_data[!cell_data$id==0,]
      s$forward_warped_data[[k]] <- cell_data
    } else{
      s$forward_warped_data[[k]] <- cell_data
    }
  }

  names(s$raw_forward_warped_data) <- channels
  names(s$forward_warped_data) <- channels

  return(s)
}

##____Method for forward warping the segmentation data per mouse____
#' @title Map cells to atlas for slice object
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere 'left', 'right' or NULL (both)
#' @param channels (vec of str) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
#' @param clean : TRUE (bool). Remove cells that don't map to any regions.
#' @param display : TRUE (bool). Display the results of the forward warp for the slice.display
#' @param replace : FALSE (bool). Replace current forward warped data, both raw and cleaned.
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the wholebrain::inspect.registration() function
#' @return m mouse object
#' @export
#'
#' @examples
#'

map_cells_to_atlas.mouse <- function(m,
                                     slice_ID = "1_10",
                                     hemisphere = 'left',
                                     channels = c('cfos', 'efyp', 'colabel'),
                                     clean =  TRUE,
                                     display = TRUE,
                                     replace = FALSE,
                                     ...){

  # Get mouse ID
  mouse_ID <- attr(m, 'info')$mouse_ID

  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }

  # check in list of previous stored slice names for a match
  stored_names <- names(m$slices)

  # Check if there is a previous slice that matches
  match <- FALSE

  for (stored_name in stored_names){
    if (identical(stored_name, slice_name)){
      match <- slice_name

      # Check if raw forward mapped data already exists
      if (!is.null(m$slices[[match]]$raw_forward_warped_data)){

        if (replace){
          # replace the slice's forward mapped data
          m$slices[[match]]$raw_forward_warped_data <- NULL
          m$slices[[match]]$forward_warped_data <- NULL

          # forward mapped data from scratch
          m$slices[[match]] <- map_cells_to_atlas(m$slices[[match]],
                                                  channels = channels,
                                                  clean = clean,
                                                  display = display,
                                                  ...)

        } else{
          stop(paste0("There is an existing segmentation object for this slice! If you want to overwrite it, set replace to TRUE."))
        }
      } else {
        # forward mapped data from scratch
        m$slices[[match]] <- map_cells_to_atlas(m$slices[[match]],
                                                channels = channels,
                                                clean = clean,
                                                display = display,
                                                ...)
      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!" ))
  }
  return(m)
}


##____Method for excluding regions, layer 1, and hemispheres per slice ____
#' @title Post-processing to exclude and clean up mapped data
#' @description Method for excluding regions, layer 1, and hemispheres, and  per slice
#' @param s: slice object
#' @param channels: channels to filter (str vector), can filter both eyfp & cfos
#' @param clean: TRUE (bool). Remove cells that don't map to any regions.
#' @param exclude_regions: NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
#' @param exclude_hemisphere: TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1: TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)
#' @export

# NOTE: # This function automatically excludes the default regions included in the attribute "regions_excluded" in each slice
# IN ADDITION to the regions added to the 'exclude_regions' parameter. Regions added to the 'exclude_regions' parameter will then be updated in the slice attribute

exclude_anatomy.slice <- function(s,
                                  channels = c('cfos', 'eyfp', 'colabel'),
                                  clean = TRUE,
                                  exclude_regions = NULL,
                                  exclude_hemisphere = TRUE,
                                  exclude_layer_1 = TRUE,
                                  plot_filtered = TRUE){


  # Get slice information & combine exclude regions parameter with default regions excluded;
  # Reassign so that all regions excluded are tracked
  info <- attributes(s)
  exclude_regions <- unique(c(exclude_regions, info$info$regions_excluded))
  attr(s, 'info')$regions_excluded <- exclude_regions

  # is.null check for exclude regions
  if (!is.null(exclude_regions)){

    # Containing all subregions to exclude too
    all_excluded_regions <- c()
    for (r in 1:length(exclude_regions)) {
      all_excluded_regions <- c(all_excluded_regions, exclude_regions[r],
                                wholebrain::get.sub.structure(exclude_regions[r]))
    }
    # May  not be necessary, included just in case
    all_excluded_regions <- unique(all_excluded_regions)
  } else {
    all_excluded_regions <- NULL
  }

  ## Filtering per channel
  for (k in 1:length(channels)){

    dataset <- s$forward_warped_data[[channels[k]]]

    # 1) Filter out regions
    dataset <- dataset[!(dataset$acronym %in% all_excluded_regions), ]

    # 2) Filter out hemisphere
    if (exclude_hemisphere){
      if (info$info$hemisphere == "right"){
        dataset <- dataset[dataset$right.hemisphere,]
      } else if (info$info$hemisphere == "left"){
        dataset <- dataset[!dataset$right.hemisphere,]
      }
    }

    # 3) Filter out layer 1
    if (exclude_layer_1){
      dataset <- dataset[-grep("layer 1",dataset$name, ignore.case = TRUE, value = FALSE),]
    }

    # Filter out cell counts that are out of bounds
    if(clean){
      dataset <- dataset[!dataset$id==0,]
      s$forward_warped_data[[k]] <- dataset
    }

    # Reassign cleaned dataset
    s$forward_warped_data[[channels[k]]] <- dataset

  }

  ## plotting schematic plots to check that the region cell counts are cleaned
  if (plot_filtered){
    for (k in 1:length(channels)){
      quartz()
      dataset <- s$forward_warped_data[[channels[k]]]
      wholebrain::schematic.plot(dataset)
    }
  }
  return(s)
}

##____Method for excluding regions, layer 1, and hemispheres per mouse ____

#' Method for excluding regions, layer 1, and hemispheres per mouse
#' @description Method for excluding cell counts in specified regions, layer 1, out-of-bounds cells counts, and hemispheres per mouse.
#'  PLEASE NOTE that this function automatically excludes the default regions included in the attribute "regions_excluded" in each slice
#'  IN ADDITION to the regions added to the 'exclude_regions' parameter.
#'  Regions added to the 'exclude_regions' parameter will then be updated in the slice attribute to keep track of overall which regions were excluded per slice.
#'
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere : hemisphere of slice
#' @param channels : channels to filter (str vector), can filter both eyfp & cfos
#' @param clean clean: TRUE (bool). Remove cells that don't map to any regions.
#' @param exclude_regions : NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
#' @param exclude_other_hemisphere : TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1 : TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)
#'
#' @return m mouse object
#' @export
#' @examples m <- exclude
exclude_anatomy.mouse <- function(m,
                                  slice_ID = "1_10",
                                  hemisphere = NULL,
                                  channels = c('cfos', 'eyfp', 'colabel'),
                                  clean = TRUE,
                                  exclude_regions = NULL,
                                  exclude_other_hemisphere = TRUE,
                                  exclude_layer_1 = TRUE){

  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
    if (exclude_other_hemisphere){
      stop("Can't exclude contralateral hemisphere without reference hemisphere! i.e., hemisphere = NULL and exclude_other_hemisphere = TRUE. ")
    }

  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }

  # print(slice_name)

  # check in list of previous stored slice names for a match
  stored_names <- names(m$slices)

  # Check if there is a previous slice that matches
  match <- FALSE

  for (stored_name in stored_names){
    if (identical(stored_name, slice_name)){
      match <- slice_name

      # Check if raw forward mapped data already exists
      if (!is.null(m$slices[[match]]$forward_warped_data)){

        # If forward mapped data exists, exclude the designated areas
        m$slices[[match]] <- exclude_anatomy(m$slices[[match]],
                                             channels = channels,
                                             clean = clean,
                                             exclude_regions = exclude_regions,
                                             exclude_hemisphere = exclude_other_hemisphere,
                                             exclude_layer_1 = exclude_layer_1)
      } else {
        stop(paste0("There no forward mapped data to clean up! Run map_cells_to_atlas() to generate a dataset first."))
      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!" ))
  }
  return(m)
}



# ##____Method for getting regional areas per slice ____
#' Method for getting regional areas per slice
#' @description returns a dataframe with columns 'name' (full region name), 'acronym', 'area', 'right.hemisphere
#' @param s : slice object
#' @return s slice object
#' @export
#' @examples s <- get_registered_areas(2)
get_registered_areas.slice <- function(s){

  # Get conversion factor from microns to pixels
    cf <- attr(s, 'info')$conversion_factor
    s$areas <- get.registered.areas(s$forward_warped_data, registration = s$registration_obj, conversion.factor = cf)
    return(s)

}


## ___ Method for getting regional areas per slice
#' Method for getting regional areas per slice in a mouse object
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere : hemisphere of slice
#' @return m mouse object
#' @export
#'
#' @examples m <- get_registered_areas(m, slice_ID = "1_10", hemisphere = "left")
get_registered_areas.mouse <- function(m,
                                       slice_ID,
                                       hemisphere = NULL,
                                       replace = FALSE){

  if (is.null(hemisphere)){
    slice_name <- slice_ID
  } else {
    slice_name <- paste0(slice_ID, "_", hemisphere)
  }

  # check in list of previous stored slice names for a match
  stored_names <- names(m$slices)

  # Check if there is a previous slice that matches
  match <- FALSE

  for (stored_name in stored_names){
    if (identical(stored_name, slice_name)){
      match <- slice_name

      # Check if areas data already exists
      if (!is.null(m$slices[[match]]$areas)){

        # Check if user wants existing data to be overwritten
        if (replace){
          m$slices[[match]] <- get_registered_areas(m$slices[[match]])

        } else {
          stop(paste0("There is an existing areas data for this slice! If you want to overwrite it, set replace to TRUE."))

        }

      } else {
        m$slices[[match]] <- get_registered_areas(m$slices[[match]])
      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!" ))
  }
return(m)
}




#__________________ Mouse object specific functions __________________________




## normalized_cell_counts per mm^3 function





#__________________ Experiment object specific functions __________________________



## Experiment function to automatically extract hippocampus as a region and split DV to a
# parameter to























#__________________ Internal Functions __________________________


## Modification of Marcos' function
get.registered.areas <- function(cell.data.list, registration, conversion.factor = 1){

  # bind the cell data into on dataframe
  cell.data <- do.call("rbind", cell.data.list)

  # create a regions dataframe
  regions.left <- unique(cell.data[!cell.data$right.hemisphere,]$acronym)
  regions.right <- unique(cell.data[cell.data$right.hemisphere,]$acronym)

  regions.left <- data.frame(acronym = regions.left, right.hemisphere = rep(FALSE, length(regions.left)), stringsAsFactors = FALSE)
  regions.right <- data.frame(acronym = regions.right, right.hemisphere = rep(TRUE, length(regions.right)), stringsAsFactors = FALSE)
  regions <- rbind(regions.left, regions.right)

  # Sort region acronyms to alphabetical order
  regions <- regions[order(regions$acronym),]

  # Get region registration information
  region.info <- list()

  for (i in 1:nrow(regions)) {
    region.data <- wholebrain::get.region(regions$acronym[i],registration)
    region.data[,1:4] <- region.data[,1:4]*conversion.factor
    region.info <- c(region.info, list(region.data[region.data$right.hemisphere==regions$right.hemisphere[i],]))
  }

  areas <- data.frame(name=as.character(wholebrain::name.from.acronym(regions$acronym)),
                      acronym=regions$acronym,
                      right.hemisphere=regions$right.hemisphere,
                      area=rep(0,nrow(regions)),
                      stringsAsFactors = FALSE)

  # Use Gauss' area formula (shoelace formula)
  for (i in 1:length(region.info)) {
    r <- region.info[[i]]
    area <- 0
    for (j in 1:(nrow(r)-1)) {
      area <- area + r$xT[j]*r$yT[j+1] - r$xT[j+1]*r$yT[j]
    }
    area <- 0.5*abs(area + r$xT[nrow(r)]*r$yT[1] - r$xT[1]*r$yT[nrow(r)])
    areas$area[i] <- area
  }
  return(areas)
}



