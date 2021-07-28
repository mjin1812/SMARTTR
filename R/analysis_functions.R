
##______________ Package imports ________________

#' @importFrom magrittr %>%
NULL



###______________Generic functions ________________###

## Creating generic function for registration
#' Register (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'

register <- function(x, ...) {
  UseMethod("register")
}

## Creating generic function for segmentation importation
#' import_segmentation (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export

import_segmentation <- function(x, ...){
  UseMethod('import_segmentation')
}

## Creating generic function for filtering raw segmentation data
make_segmentation_filter <- function(x, ...){
  UseMethod('make_segmentation_filter')
}

## Create segmentation object
#' make_segmentation_object (generic funciton)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export

make_segmentation_object <- function(x, ...){
  UseMethod('make_segmentation_object')
}


## forward_warp segmentation data to atlas space
#' map_cells_to_atlas (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export


map_cells_to_atlas <- function(x, ...){
  UseMethod('map_cells_to_atlas')
}

## excluded designated regions from the mapped cell data table
#' exclude_anatomy (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export

exclude_anatomy <- function(x, ...){
  UseMethod('exclude_anatomy')
}

## excluded designated regions from the mapped cell data table
#' get_registered_areas (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export

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
#' @rdname register
#' @param s slice object
#' @param filter (list)
#' @param ... additional parameters to pass to the SMART::registration2() function, besides 'input', 'coordinate', 'filter' & 'correspondance'
#' @return s slice object
#' @examples s <- register(s)
#' @export

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
#' Register a slice in a mouse object. If a slice has been previously registered, the default behavior is to continue modifying the
#' previous registration. Use the `replace` parameter to change this behavior.
#'
#' @rdname register
#' @param m  mouse object
#' @param slice_ID ID of slice (str)
#' @param hemisphere 'left', 'right' (str) or NULL if both hemispheres are included
#' @param filter (list)
#' @param replace (bool, default = FALSE) Replace a registration already contained in a mouse object by resetting to NULL value before registration improvement loop.
#' @param ... additional parameters to pass to the SMART::registration2() function, besides 'input', 'coordinate', 'filter' & 'correspondance'
#' @return m  mouse object
#'
#' @examples m <- register(m, slice_ID = '1_10', hemisphere = "left", filter = my_filter)
#' @export
# TODO: add popup window option to better visualize the internal structures of the images
register.mouse <- function(m,
                           slice_ID = '1_10',
                           hemisphere = NULL,  # set hemisphere to 'left' or 'right' to specify which slice should be pulled
                           filter = NULL,
                           replace = FALSE,
                           ...){


  # get or set output path
  output_path <- attr(m, 'info')$output_path
  if (is.null(output_path)) {
    output_path <- "../"
  } else {
    output_path <- file.path(output_path, "registrations")

    if (!file.exists(output_path)){
      dir.create(output_path)
    }
  }

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
                                        filter = filter,
                                        output.folder = output_path, ...)
        } else{
          message("Improving previous registration for this slice. Set replace = TRUE to start from scratch!...\n")

          if (is.null(filter)){
            filter <- SMARTR::segmentation.object$filter
            filter$resize <-  m$slices[[match]]$registration_obj$resize
          }

          # Register matched slices with existing data
          m$slices[[match]] <- register(m$slices[[match]],
                                        filter = filter,
                                        output.folder = output_path, ...)
        }
      } else  {
        message("Starting a new registration!")

        # Register matched slices with no registration data yet
        m$slices[[match]] <- register(m$slices[[match]],
                                      filter = filter,
                                      output.folder = output_path, ...)

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
#' @rdname import_segmentation
#' @description Method for importing segmentation data for a slice object
#' @param s : slice object
#' @param mouse_ID : ID of mouse (str)
#' @param channels : channels to import (str vector)
#' @return s slice object
#' @examples s <-  import_segmentation(s, mouse_ID = "255", channels = c("cfos", "eyfp", "colabel"))
#' @export

# TODO: modify paths to search for correct files regardless of upper or lower case. Use grep and stringi to match multiple patterns

import_segmentation.slice <- function(s,
                                      mouse_ID = '255',
                                      channels = c('cfos', 'eyfp', 'colabel')){


  # Get slice information
  info <- attributes(s)

  ## Use registration image path to find segmentation data path stem
  # Check if slice has a set directory, otherwise borrow from the registration path
  if (is.null(info$info$slice_directory)){
    if (info$info$registration_path == 'set registration image path'){
      stop(paste0("This function doesn't know where to look to import segmentation data!\n",
                  "Set either the registration_path or the slice_directory in your slice attributes to continue."))
    } else {
      path_stem <- dirname(info$info$registration_path)
    }
  } else{
    path_stem <- info$info$slice_directory
  }

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
#'
#' @rdname import_segmentation
#' @description Method for importing segmentation data for a mouse object
#' @param  m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere : 'left', 'right' or NULL (both)
#' @param channels : channels to import (str vector)
#' @param replace : replace existing raw segmentation data (bool)
#' @return m mouse object
#' @examples m <-  import_segmentation(m, slice_ID = "1_10", channels = c("cfos", "eyfp", "colabel"), replace = FALSE)
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
    filter <- make.filter(data, params = params[[k]], ranges = ranges[[k]])
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


#' Creates a wholebrain segmentation object and stores it.
#'
#' @rdname make_segmentation_object
#' @description Make a wholebrain compatible segmentation object for a slice in a slice object
#' @param s : slice object
#' @param mouse_ID : ID of mouse (str)
#' @param channels : (str vec) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
#' @param use_filter : TRUE (bool). Use a filter to create more curated segmentation object from the raw segmentation data.
#'
#' @return s slice object
#' @examples s <- make_segmentation_object(s, mouse_ID = "255", channels = c("cfos", "eyfp"), use_filter = FALSE)
#' @export

# Note: this function needs cleaning; for now always create a filter for the
# 'eyfp' channel and and set the 'use_filter' variable to TRUE.

# Note: this function needs cleaning; for now always create a filter for the
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
#' @rdname make_segmentation_object
#' @description Make a wholebrain compatible segmentation object for a slice in a mouse object
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere: 'left', 'right' or NULL (both)
#' @param channels : (vec) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
#' @param replace : replace existing raw segmentation data (bool)
#' @param use_filter: TRUE (bool). Use a filter to create more curated segmentation object from the raw segmentation data.
#' @examples m <-  make_segmentation_object(m, slice_ID = '1_9', hemisphere = 'left', channels = c('eyfp', 'cfos'), use_filter = FALSE)
#' @export


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
#' @rdname map_cells_to_atlas
#' @param s : slice object
#' @param channels : channels to filter (str vector), can filter both eyfp & cfos
#' @param clean : TRUE (bool). Remove cells that don't map to any regions.
#' @param display : TRUE (bool). Display the results of the forward warp for the slice.
#' @param mouse_ID : (str) mouse ID
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the wholebrain::inspect.registration() function
#' @examples s <- map_cells_to_atlas(s, channels c('cfos' , 'eyfp', 'colabel'), clean = TRUE, display = TRUE, mouse_ID = "255")
#' @export


map_cells_to_atlas.slice <- function(s,
                                     channels = c('cfos', 'eyfp', 'colabel'),
                                     clean =  TRUE,
                                     display = TRUE,
                                     mouse_ID = NULL,
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

    if (!is.null(mouse_ID)){
      s$forward_warped_data[[k]]$animal <-  mouse_ID
    }


  }

  names(s$raw_forward_warped_data) <- channels
  names(s$forward_warped_data) <- channels

  return(s)
}

##____Method for forward warping the segmentation data per mouse____
#' @title Map cells to atlas for slice within a mouse object
#' @description  Method for forward warping segmentation data to atlas space for a slice within a mouse object.
#' Requires segmentation objects to be made for the channels specified and a registration.
#' @rdname map_cells_to_atlas
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere 'left', 'right' or NULL (both)
#' @param channels (vec of str) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
#' @param clean : TRUE (bool). Remove cells that don't map to any regions.
#' @param display : TRUE (bool). Display the results of the forward warp for the slice.display
#' @param replace : FALSE (bool). Replace current forward warped data, both raw and cleaned.
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the wholebrain::inspect.registration() function
#' @return m mouse object
#' @examples m <- map_cells_to_atlas(m, slice_ID = "1_10", hemisphere = NULL, channels = c("cfos", "eyfp", "colabel"), clean = TRUE, replace = FALSE)
#' @export

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
                                                  mouse_ID = mouse_ID,
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
                                                mouse_ID = mouse_ID,
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
#' @description Method for excluding user specified regions, layer 1, and contralateral hemisphere  per slice.
#' This function automatically excludes the default regions included in the attribute "regions_excluded" in each slice
#' IN ADDITION to the regions added to the 'exclude_regions' parameter.
#' Regions added to the 'exclude_regions' parameter will then be updated in the slice attribute to keep track of what was excluded.
#' @rdname exclude_anatomy
#' @param s: slice object
#' @param channels: channels to filter (str vector), can filter both eyfp & cfos
#' @param clean: TRUE (bool). Remove cells that don't map to any regions.
#' @param exclude_regions: NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
#' @param exclude_hemisphere: TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1: TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)
#' @param plot_filtered : TRUE (bool) pop up window to check the excluded anatomy.
#' @examples s <-  exclude_anatomy(s, channels = c('cfos', 'eyfp', 'colabel'), clean = TRUE, exclude_regions = NULL, exclude_hemisphere = TRUE, exclude_layer_1 = TRUE, plot_filtered = TRUE)
#' @export


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

      # Get rid of rows that contain any NA values
      dataset <- tidyr::drop_na(dataset)
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

#' Method for excluding regions, layer 1, and hemispheres in a slice per mouse
#' @description Method for excluding cell counts in specified regions, layer 1, out-of-bounds cells counts, and hemispheres per mouse.
#' This function automatically excludes the default regions included in the attribute "regions_excluded" in each slice
#' IN ADDITION to the regions added to the 'exclude_regions' parameter.
#' Regions added to the 'exclude_regions' parameter will then be updated in the slice attribute to keep track of what was excluded.
#' @rdname exclude_anatomy
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere : hemisphere of slice
#' @param channels : channels to filter (str vector), can filter both eyfp & cfos
#' @param clean clean: TRUE (bool). Remove cells that don't map to any regions.
#' @param exclude_regions : NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
#' @param exclude_other_hemisphere : TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1 : TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)
#' @return m mouse object
#' @examples m <- exclude_anatomy(m, slice_ID = "1_10", hemisphere = NULL, channels = c('cfos', 'eyfp', 'colabel'), clean = TRUE,
#' exclude_regions = NULL,
#' exclude_other_hemisphere = TRUE,
#' exclude_layer_1 = TRUE
#' @export


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
#' @rdname get_registered_areas
#' @description returns a dataframe with columns 'name' (full region name), 'acronym', 'area' (in microns^2^), 'right.hemisphere'
#' @param s : slice object
#' @return s slice object
#' @examples s <- get_registered_areas(s)
#' @export
get_registered_areas.slice <- function(s){

  # Get conversion factor from microns to pixels
    cf <- attr(s, 'info')$conversion_factor
    s$areas <- get.registered.areas(s$forward_warped_data, registration = s$registration_obj, conversion.factor = cf)
    tidyr::drop_na(s$areas)

    return(s)
}


## ___ Method for getting regional areas per slice
#' Method for getting regional areas per slice in a mouse object
#' @rdname get_registered_areas
#' @param m : mouse object
#' @param slice_ID : ID of slice (str)
#' @param hemisphere : hemisphere of slice
#' @return m mouse object
#' @examples m <- get_registered_areas(m, slice_ID = "1_10", hemisphere = "left", replace = FALSE)
#' @export
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




#' Get cell tables
#' @description This function stores a list in a mouse object that is the length of the channels parameter. Each element of the list is a dataframe
#' containing combined cell counts for one channel across all slices processed for that mouse.
#' By default, if one slices has no dataset processed for that particular channel, that slice will be skipped over.
#' The function will run properly but a warning will be thrown indicating you should go back and
#' generate a mapped dataset for that particular slice and channel.
#' @param m mouse object
#' @param channels  channels : c("cfos", "eyfp") (vec) Channels to process.
#' @return m mouse object
#' @examples m <- get_cell_tables(m, channels = c("cfos"m "eyfp", "colabel"))
#' @export

get_cell_table <- function(m,
                           channels = c("cfos", "eyfp")){

  #A list the length of "channels" combining cell counts across all slices
  #All NA values are removed to avoid an error so just check that each data frame is not NULL, otherwise print an error

  cell_table <- vector(mode = "list", length = length(channels))
  names(cell_table) = channels

  for (channel in channels){
    for (s in m$slices){

      # Get slice information
      s_info <- attr(s, "info")
      hemisphere <- s_info$hemisphere
      slice_ID <-  s_info$slice_ID

      # Omit hemisphere name if there is no hemisphere value in the attributes
      if (is.null(hemisphere)){
        slice_name <- slice_ID
      } else {
        slice_name <- paste0(slice_ID, "_", hemisphere)
      }


      # Check for empty dataframe
      if  (is.null(s$forward_warped_data[[channel]])){
        warning(paste0("\n The forward warped data from slice ", slice_name,
                       ", channel: ", channel, " is missing. Skipping this slice while creating a combined cell table.",
                       "\nIf you want to include it, go back and create the mapped dataset for this slice, then re-run this function! "))
      }

      # Add the slice name to the dataset
      s$forward_warped_data[[channel]]$slice_name <- rep(slice_name, times = length(s$forward_warped_data[[channel]]$animal))

      # Bind the dataset
      cell_table[[channel]] <- rbind(cell_table[[channel]], s$forward_warped_data[[channel]])

    }

    # Add the mouse ID and remove rows with NA values
    # cell_table[[channel]]$<- <- rep(attr(m, "info")$mouse_ID, times = length(cell_table[[channel]]$animal))
    cell_table[[channel]] <- tidyr::drop_na(cell_table[[channel]])
  }

  m$cell_table <-  cell_table
  return(m)
}



## normalized_cell_counts per mm^3 function
#
#' Normalize cell counts per mm^^2 or by mm^3^ (if multiplying by the stack size)
#' @description Run this function after all the slices that you want to process are finished being added
#' and you have run the function get_cell_table() for your mouse.
#' @param m mouse object
#' @param combine_hemispheres bool : Combine normalized cell counts from both hemispheres
#' @param simplify_regions bool : simplify the normalized region counts based on keywords in `simplify_keywords`
#' @param simplify_keywords str vector : c("layer","part","stratum","division"), keywords to search through region names and simplify to parent structure
#' @return
#' @examples m <- normalize_cell_counts(m, combine_hemispheres = TRUE, simplify_regions = TRUE)
#' @export
normalize_cell_counts <- function(m,
                                  combine_hemispheres = TRUE,
                                  simplify_regions = TRUE,
                                  simplify_keywords = c("layer","part","stratum","division")){


  # 1) NULL areas check

    for(s in m$slices){

      # Get slice information
      s_info <- attr(s, "info")
      hemisphere <- s_info$hemisphere
      slice_ID <-  s_info$slice_ID

      # Omit hemisphere name if there is no hemisphere value in the attributes
      if (is.null(hemisphere)){
        slice_name <- slice_ID
      } else {
        slice_name <- paste0(slice_ID, "_", hemisphere)
      }

      ## Check that this slice has run get_registered_areas
      if (is.null(s$areas)){
        stop(paste0("Slice ", slice_name ," in your mouse dataset has an empty areas vector.\n",
                    "Run the function get_registered_areas() for this slice before you can normalized all your cell counts by total area or by volume."))
        }
    }


  # 2)  Getting total areas across all slices, separated by hemisphere

  # rbind all the areas table
  aggregate_areas <- m$slices[[1]]$areas
  if (length(m$slices)  >=  2){
    for (k in 2:length(m$slices)){
      aggregate_areas <- rbind(aggregate_areas, m$slices[[k]]$areas)
    }
  }

  # Drop any NAs
  aggregate_areas <- tidyr::drop_na(aggregate_areas)

  # Use aggregate to sum up areas from the same region across slices
  total_areas <- aggregate(aggregate_areas$area,
                           by = list(acronym = aggregate_areas$acronym, right.hemisphere = aggregate_areas$right.hemisphere),
                           FUN =sum) %>% dplyr::rename(area = x)

  # Get z_width from slice attributes of the first slice
  z_width <- attr(m$slices[[1]], "info")$z_width


  # then use Cell_table to create combined cell counts normalized by area
  normalized_counts <- normalize.registered.areas(m$cell_table, total_areas, z_width = z_width ) #tabulate number of cells in each region
  names(normalized_counts) <- names(m$cell_table)


  # Clean NAs, check for combining hemispheres or simplifying regions
  for (channel in names(normalized_counts)){

    # clean NA values in all channels
    normalized_counts[[channel]] <- tidyr::drop_na(normalized_counts[[channel]])

    # Combine hemispheres
    if (combine_hemispheres){
      normalized_counts[[channel]] <- plyr::ddply(normalized_counts[[channel]], c("acronym", "name"), numcolwise(sum))

      # Recalculated normalized count by area
      normalized_counts[[channel]]$normalized.count.by.area <- normalized_counts[[channel]]$count/normalized_counts[[channel]]$area.mm2

      # Recalculate normalized count by volume
      normalized_counts[[channel]]$normalized.count.by.volume <- normalized_counts[[channel]]$count/normalized_counts[[channel]]$volume.mm3

    }

    # Simplify regions
    if (simplify_regions){

      # normalize the cell counts by specific keyworks
      normalized_counts <- simplify.regions(normalized_counts, keywords = simplify_keywords)
    }
  }

  # Store the normalized counts in the mouse object
  m$normalized_counts <- normalized_counts
  return(m)
}







#__________________ Experiment object specific functions __________________________



## Experiment function to automatically extract hippocampus as a region and split DV to a
# parameter to















#__________________ Internal Functions __________________________


## Modification of Marcos' function
#' Title
#'
#' @param cell.data.list
#' @param registration
#' @param conversion.factor
#'
#' @return
#'
#' @examples
get.registered.areas <- function(cell.data.list, registration, conversion.factor = 1){

  # bind the channel cell data into a single dataframe
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
    region.data <- wholebrain::get.region(regions$acronym[k],registration)
    region.data[,1:4] <- region.data[,1:4]*conversion.factor
    region.info <- c(region.info, list(region.data[region.data$right.hemisphere==regions$right.hemisphere[k],]))
  }

  areas <- data.frame(name=as.character(wholebrain::name.from.acronym(regions$acronym)),
                      acronym=regions$acronym,
                      right.hemisphere=regions$right.hemisphere,
                      area=rep(0,nrow(regions)),
                      stringsAsFactors = FALSE)

  # Use Gauss' area formula (shoelace formula)
  for (i in 1:length(region.info)) {
    r <- region.info[[k]]
    area <- 0
    for (j in 1:(nrow(r)-1)) {
      area <- area + r$xT[j]*r$yT[j+1] - r$xT[j+1]*r$yT[j]
    }
    area <- 0.5*abs(area + r$xT[nrow(r)]*r$yT[1] - r$xT[1]*r$yT[nrow(r)])
    areas$area[k] <- area
  }

  return(areas)
}


#____

## Modification of Marcos' function

#' normalize registered areas
#'
#' @param cell_table cell_tables list from get_cell_table()
#' @param total_areas computed total areas list
#' @param z_width stack width in microns (default = 24)
#'
#' @return
#'
#' @examples
normalize.registered.areas <- function(cell_table, total_areas, z_width = 24){

  # Initializing a normalized cell count list that will store df of  normalized cell counts
  normalized.list <- vector(mode = "list", length = length(cell_table))
  names(normalized.list) <- names(cell_table)

  for(channel in names(cell_table)) {

    # Get cell data for channel
    cell.data <-  cell_table[[channel]]

    # Table tally of how many counts there are in a given region
    counts <- as.data.frame(table(cell.data$right.hemisphere, cell.data$acronym), stringsAsFactors = FALSE)
    names(counts) <- c("right.hemisphere","acronym","count")

    # Merge the total areas with the counts df by acronym and hemisphere
    areas <- merge(counts, total_areas, by=c("acronym","right.hemisphere"), all = TRUE)


    # Normalize count by area mm^2 and volume
    areas$area.mm2 <- areas$area*1e-6
    areas$volume.mm3 <- areas$area*z_width*1e-9

    areas$normalized.count.by.area <- areas$count/areas$area.mm2
    areas$normalized.count.by.volume <- areas$count/areas$volume.mm3

    # Add acronym names
    areas$name <- as.character(name.from.acronym(areas$acronym))


    # store the data into the normalized list
    normalized.list[[channel]] <- areas
  }

  return(normalized.list)
}


#' Modification of Marcos' combine regions function
#'
#' @param normalized_counts list with length = No. channels., each channel element is df of normalized counts
#' @param keywords a vector of keywords with which to simplify the region names
#'
#' @return a list with length = No. channels. Each channel df has simplified parent acronyms and names is summed from original df based on them.
#'
#' @examples
simplify.regions <- function(normalized_counts, keywords = c("layer","part","stratum","division")) {


  # Initialize empty list vector to store the simplified counts
  simplified_counts <- vector(mode = "list", length = length(normalized_counts))
  names(simplified_counts) <- names(normalized_counts)

  for(channel in names(normalized_counts)) {

    # Extract the dataframe for one channel
    simplified_counts_chan <- normalized_counts[[channel]]

    # loop through each keyword
    for (s in keywords) {

      # Look for row indices where the names contain the keywords
      k <- grep(s, simplified_counts_chan$name, value = FALSE, ignore.case = TRUE)

      while(length(k) > 0){

        # Store the parent acronym and full name
        simplified_counts_chan$acronym[k] <- get.acronym.parent(simplified_counts_chan$acronym[k])
        simplified_counts_chan$name[k] <- as.character(name.from.acronym(simplified_counts_chan$acronym[k]))

        # Continue storing storing the parent names until there are no more that match the keyword
        k <- grep(s, simplified_counts_chan$name, value = FALSE, ignore.case = TRUE)
      }
    }

    # Check if the data has been collapsed by hemisphere
    if (is.null(simplified_counts_chan$right.hemisphere)){

      # step to collapse by name/acronym
      simplified_counts_chan <- plyr::ddply(simplified_counts_chan, c("acronym", "name"), numcolwise(sum))
    } else {
      # step to collapse by name/acronym and hemisphere
      simplified_counts_chan <- plyr::ddply(simplified_counts_chan, c("acronym", "name", "right.hemisphere"), numcolwise(sum))

    }

    # Recalculated normalized count by area
    simplified_counts_chan$normalized.count.by.area <- simplified_counts_chan$count/simplified_counts_chan$area.mm2

    # Recalculate normalized count by volume
    simplified_counts_chan$normalized.count.by.volume <- simplified_counts_chan$count/simplified_counts_chan$volume.mm3

    # Store in simplified counts list
    simplified_counts[[channel]] <- simplified_counts_chan
  }

  return(simplified_counts)
}

###################### Unedited OLD functions #################




####################################

make.filter <- function(data, params = c("Vol..unit.","Moment1","Moment2","Moment3","Moment4","Sigma"),
                             ranges = list(c(200,12000),c(3,50),c(0,600),c(0,2000),c(0,5),c(20,Inf))){
  non.cells <- c()
  for (i in 1:length(params)){
    non.cells <- c(non.cells,which(data[,params[k]] < ranges[[k]][1] | data[,params[k]] > ranges[[k]][2]))
  }
  return(unique(non.cells))
}

#######################################





# eyfp.counts  - dataframe
# eyfp.counts.16bit
#
#





# Modified version of Marcos identify.colabelled.cells.filter() functions
#' Get colabelled cells data table
#'
#' @param coloc.table
#' @param eyfp.counts
#' @param eyfp.counts.16bit
#' @param volume
#' @param overlap
#'
#' @return returns a dataframe of colabelled cell counts
#'
#' @examples
get.colabeled.cells <- function(coloc.table, eyfp.counts, eyfp.counts.16bit, volume = 25, overlap = 0.5) {
  # p is position
  # mi is max index
  # mp is max proportion
  # mv is max volume
  # mo is max object number for eyfp


  # Get position index of the volume parameteres in the coloc.table
  p <- seq(4,ncol(coloc.table),3)
  names(p) <- names(coloc.table)[p]  # Helper line to keep track

  # Get column indices of the volumn column with the largest objects overlap
  mi <- max.col(coloc.table[,p])

  # extracting column for objects, volumes, and proportions based on names
  mp_names <- paste0("P", 1:7)
  mv_names <- paste0("V", 1:7)
  mo_names <- paste0("O", 1:7)

  # Get selection indices
  select <- cbind(1:length(coloc.table$X),mi)

  # Extracting max proportion
  mp <- coloc.table[,mp_names][select]
  mv <- coloc.table[,mv_names][select]
  mo <- coloc.table[,mo_names][select]


  # Filter out objects that are smaller than the volume threshold and the less than
  # The proportion overlap threshold
  mo <- mo[mv>=volume & mp>=overlap]



  # Split by object name and value, split by character position
  obj.val.16 <- strsplit(eyfp.counts.16bit$Name,"-") %>% lapply(substr, start = 4, stop = 10) %>% unlist()
  val.16 <-  obj.val.16[seq(2, length(obj.val.16), by = 2)] # Object value
  obj.16 <-  obj.val.16[seq(1, length(obj.val.16), by = 2)] # Object number

  obj_index <- match(mo, val.16) %>% na.omit()
  mot <- obj.16[obj_index] %>% as.integer()

  # mot is matched object number from 16bit measure df that is in mo
  obj.val <- strsplit(eyfp.counts$Name,"-") %>% lapply(substr, start = 4, stop = 10) %>% unlist()
  obj <-  obj.val[seq(1, length(obj.val), by = 2)] # Object number


  # index of the eyfp rows corresponding to matched object name
  index <- match(mot, obj) %>% na.omit()
  coloc.data <- eyfp.counts[index, unique(names(eyfp.counts))] %>% dplyr::distinct()

  return(coloc.data)




  ## EXTRAS
  # ## get the matched object number
  # matched_obj <- mot[match(obj , mot)]

  ## compare the X&Y coordinates of the 16bit and normal based on matched object number
  # val <-  obj.val[seq(2, length(obj.val), by = 2)] # Object value
  # df.16 <- eyfp.counts.16bit[obj_index,] %>% dplyr::distinct()

}


