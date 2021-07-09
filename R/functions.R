# library(SMART)
# library(devtools)
# library(tidyverse)





#' @title Add slice to a mouse object
#' @usage m <- add_slice(m, s, replace = FALSE)
#' @param m mouse object
#' @param s slice object
#' @param replace (bool, default = FALSE) Replace slice already contained in a mouse object.
#' @export

add_slice <- function(m, s, replace = FALSE){
  # Get slice ID & hemisphere
  slice_ID <- attr(s, 'info')$slice_ID
  hemisphere <- attr(s, 'info')$hemisphere

  # First slice
  if (length(m$slices) < 1){
    m$slices[[1]] <- s

  } else{

     # Get slice information
      info  <- list()
      index <- vector(mode = 'integer')

      # Check for previous matches
      for (k in 1:length(m$slices)){
        if (tolower(attr(m$slices[[k]], 'info')$slice_ID) == tolower(slice_ID)){     # case insensitive mouse id check
          index <- c(index, k)
         info[[length(index)]] <- attr(m$slices[[k]], 'info')
       }
     }

      # If this is a new slice, append to end of existing list
      if (length(index) < 1 ){
        index <-  length(m$slices) + 1
        m$slices[[index]] <- s

      } else {

        if (length(index) > 1 && !is.null(hemisphere)){ # If there are multiple slice matches, match by the hemisphere selected

          if (tolower(info[[1]]$hemisphere) == tolower(hemisphere)){   # If the first match is correct hemisphere, remove other match
            info[[2]] <- NULL
            index <- index[1]
          } else {                       #If the second match is the correct hemisphere
            info[[1]] <- NULL
            index <- index[2]
          }


          if (replace){
            ## add print statement indicating which hemisphere is being registered
            message(paste0('Replacing existing data for slice ',info[[1]]$slice_ID,
                           ", ", hemisphere, ' hemisphere'))
            m$slices[[index]] <- s

          } else {
            message(paste0('There is existing data for slice ',info[[1]]$slice_ID,
                           ", ", hemisphere, ' hemisphere\n',
                           'If you want to replace this previous slice object, set the "replace" argument to "TRUE".'))
          }


        } else if (length(index) == 1 &&  !is.null(hemisphere)){

          if (tolower(info[[1]]$hemisphere) == tolower(hemisphere)){

            if (replace){
              ## add print statement indicating which hemisphere is being registered
              message(paste0('Replacing existing data for slice ',info[[1]]$slice_ID,
                             ", ", hemisphere, ' hemisphere'))
              m$slices[[index]] <- s

            } else {
              message(paste0('There is existing data for slice ',info[[1]]$slice_ID,
                             ", ", hemisphere, ' hemisphere\n',
                             'If you want to replace this previous slice object, set the "replace" argument to "TRUE".'))
            }

          } else {
            index <-  length(m$slices) + 1
            m$slices[[index]] <- s
          }

        } else if (length(index) > 1){

             stop(paste0('There are multiple slice objects in this mouse with the same ID! \n',
                          'If you want to replace one of these previous slice objects, set the hemisphere attribute to specify which slice replace\n',
                          'then set the "replace" argument to "TRUE".'))
        }
      }
  }
  return(m)
}





###______________Generic functions ________________###

## Creating generic function for registration
register <- function(x, ...) {
  UseMethod("register")
}

## Creating generic function for segmentation importation
import_segmentation <- function(x, ...){
  UseMethod('import_segmentation')
}

## Creating generic function for filtering raw segmentation data
make_segmentation_filter <- function(x, ...){
  UseMethod('make_segmentation_filter')
}

## Create segmentation object
make_segmentation_object <- function(x, ...){
  UseMethod('make_segmentation_object')
}


## forward_warp segmentation data to atlas space
map_cells_to_atlas <- function(x, ...){
  UseMethod('map_cells_to_atlas')
}

## excluded designated regions from the mapped cell data table
exclude_anatomy <- function(x, ...){
  UseMethod('exclude_anatomy')
}

## excluded designated regions from the mapped cell data table
get_registered_areas <- function(x, ...){
  UseMethod('get_registered_areas')
}



#
#
# #_________ TESTING NEW SEGMENTATION METHOD _________
#
# ## Creating generic function for segmentation importation
# import_segmentation_weka <- function(x, ...){
#   UseMethod('import_segmentation_weka')
# }


###_____________ Methods for generic functions _________________####

##____ Creating method for registration of slice object____
# s: slice object
# filter: (list)

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
    message(paste0("\nRegistration data for this slice already exists.\nThe existing registration data will be used.",
                "\nSet the registration data to NULL if you want to start from scratch.\n"))
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
# m : mouse object
# slice_ID: ID of slice (str)
# hemisphere: 'left', 'right' (str) or NULL (if both hemispheres are included)
# filter: (list)

# TODO: maybe? add looping through all slices option or have that be outside of this specific function
# TODO: add popup window option to better visualize the internal structures of the images

register.mouse <- function(m,
                           slice_ID = '1_10',
                           hemisphere = NULL,  # set hemisphere to 'left' or 'right' to specify which slice should be pulled
                           filter = NULL,
                           ...){
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
            'Try setting the hemisphere argument if you want to specify which slice to use,\n',
             'otherwise all matched slices will be looped through.'))
  }

  for (k in 1:length(index)){

    ## add print statement indicating which hemisphere is being registered

    message(paste0('\nThis is slice ',info[[k]]$slice_ID))
    message(paste0('Registering the ',info[[k]]$hemisphere, ' hemisphere'))

    # Register matched slices
    m$slices[[index[k]]] <- register(m$slices[[index[k]]],
                              filter = filter, ...)
  }
  return(m)
}


##____ Creating method for importing segmentation data for a slice object____
# s : slice object
# mouse_ID: ID of mouse (str)
# channels: channels to import (str vector)

# TODO: modify paths to search for correct files regardless of upper or lower case. Use grep and stringi to match multiple patterns

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

      meas <- read.csv(meas_path)
      quant <- read.csv(quant_path)
      meas$X2_Pix <- meas$CX..pix./info$info$bin #create position column to account for binning
      meas$Y2_Pix <- meas$CY..pix./info$info$bin #same as above
      counts <- cbind(meas, quant) #create combined table
      s$raw_segmentation_data[[k]] <- counts
    }
  }
  # Name the elements of the raw data
  names(s$raw_segmentation_data) <- channels
  return(s)
}


##____ Creating method for importing segmentation data for a mouse object____
# m : mouse object
# slice_ID : ID of slice (str)
# hemisphere: 'left', 'right' or NULL (both)
# channels: channels to import (str vector)

import_segmentation.mouse <- function(m,
                                      slice_ID = '1_10',
                                      hemisphere = NULL,
                                      channels = c('cfos', 'efyp', 'colabel')){

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
                   'Try setting the hemisphere argument if you want to specify which slice to use,\n',
                   'otherwise all matched slices will have segmentation data imported.'))
  }

  for (k in 1:length(index)){

    ## add print statement indicating which hemisphere is being registered
    message(paste0('This is slice ',info[[k]]$slice_ID))
    message(paste0('Importing segmentation data for the ',info[[k]]$hemisphere, ' hemisphere'))

    m$slices[[index[k]]] <- import_segmentation(m$slices[[index[k]]],
                                                mouse_ID = mouse_ID,
                                                channels = channels)
  }

  return(m)
}


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

##____Method for creating segmentation object for each slice____
# s : slice object
# mouse_ID: ID of mouse (str)
# channels : (vec) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
# use_filter: TRUE (bool). Use a filter to create more curated segmentation object from the raw segmentation data.


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
# m : mouse object
# slice_ID : ID of slice (str)
# hemisphere: 'left', 'right' or NULL (both)
# channels : (vec) Channels to process; Note: always set 'eyfp' before 'colabel' in the vector
# use_filter: TRUE (bool). Use a filter to create more curated segmentation object from the raw segmentation data.


make_segmentation_object.mouse <- function(m,
                                           slice_ID = '1_10',
                                           hemisphere = NULL,
                                           channels = c('cfos', 'eyfp', 'colabel'),
                                           use_filter = TRUE){

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
                   'Set the hemisphere argument if you want to specify which slice to create segmentation object for,\n',
                   'otherwise all matched slices will have a segmentation object created.'))
  }

  for (k in 1:length(index)){


    ## add print statement indicating which hemisphere is being registered
    message(paste0('This is slice ',info[[k]]$slice_ID))
    message(paste0('Creating segmentation object(s) for the ',info[[k]]$hemisphere, ' hemisphere'))

    m$slices[[index[k]]] <- make_segmentation_object(m$slices[[index[k]]],
                                                      mouse_ID = mouse_ID,
                                                      channels = channels,
                                                      use_filter = use_filter
                                                      )
  }
  return(m)
}


##____Method for forward warping the segmentation data per slice____

# s: slice object
# channels: channels to filter (str vector), can filter both eyfp & cfos
# clean: TRUE (bool). Remove cells that don't map to any regions.
# display: TRUE (bool). Display the results of the forward warp for the slice.

map_cells_to_atlas.slice <- function(s,
                                     channels = c('cfos', 'eyfp', 'colabel'),
                                     clean =  TRUE,
                                     display = TRUE){

  for (k in 1:length(channels)){
    cell_data <- inspect.registration(s$registration_obj,
                                      s$segmentation_obj[[channels[k]]],
                                      forward.warps = TRUE,
                                      device = display)
    fil
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

map_cells_to_atlas.mouse <- function(m,
                                     slice_ID = "1_10",
                                     hemisphere = 'left',
                                     channels = c('cfos', 'efyp', 'colabel'),
                                     clean =  TRUE,
                                     display = TRUE){

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
                   'Set the hemisphere argument if you want to specify which slice to select,\n',
                   'otherwise all matched slices will be forward warped.'))
  }


  for (k in 1:length(index)){
    ## add print statement indicating which hemisphere is being registered
    message(paste0('This is slice ',info[[k]]$slice_ID))
    message(paste0('Forward warping the ',info[[k]]$hemisphere, ' hemisphere'))

    m$slices[[index[k]]] <- map_cells_to_atlas(m$slices[[index[k]]],
                                              channels = channels,
                                              clean = clean,
                                              display = display)
  }
  return(m)
}



##____Method for excluding regions, layer 1, and hemispheres per slice ____


# s: slice object
# channels: channels to filter (str vector), can filter both eyfp & cfos
# exclude_regions: NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
# exclude_hemisphere: TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
# exclude_layer_1: TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)

# NOTE: # This function automatically excludes the default regions included in the attribute "regions_excluded" in each slice
# IN ADDITION to the regions added to the 'exclude_regions' parameter. Regions added to the 'exclude_regions' parameter will then be updated in the slice attribute

exclude_anatomy.slice <- function(s,
                                  channels = c('cfos', 'eyfp', 'colabel'),
                                  exclude_regions = NULL,
                                  exclude_hemisphere = TRUE,
                                  exclude_layer_1 = TRUE){


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

    # Reassign cleaned dataset
    s$forward_warped_data[[channels[k]]] <- dataset

  }
  return(s)
}

##____Method for excluding regions, layer 1, and hemispheres per mouse ____
# m: mouse object
# slice_ID : ID of slice (str)
# channels: channels to filter (str vector), can filter both eyfp & cfos
# exclude_regions: NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
# exclude_hemisphere: TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
# exclude_layer_1: TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)

# NOTE: # This function automatically excludes the default regions included in the attribute "regions_excluded" in each slice
# IN ADDITION to the regions added to the 'exclude_regions' parameter. Regions added to the 'exclude_regions' parameter will then be updated in the slice attribute

exclude_anatomy.mouse <- function(m,
                                  slice_ID = "1_10",
                                  channels = c('cfos', 'eyfp', 'colabel'),
                                  exclude_regions = NULL,
                                  exclude_hemisphere = TRUE,
                                  exclude_layer_1 = TRUE){

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
                   'Set the hemisphere argument if you want to specify which slice to clean.,\n',
                   'otherwise all matched slices will be cleaned.'))
  }


  for (k in 1:length(index)){
    ## add print statement indicating which hemisphere is being registered
    message(paste0('This is slice ',info[[k]]$slice_ID))
    message(paste0('Cleaning cell counts from the ',info[[k]]$hemisphere, ' hemisphere'))

    m$slices[[index[k]]] <- exclude_anatomy.mouse(m$slices[[index[k]]],
                                                  channels = channels,
                                                  exclude_regions = exclude_regions,
                                                  exclude_hemisphere = exclude_hemisphere,
                                                  exclude_layer_1 = exclude_layer_1)
  }
  return(m)
}





##____Method for getting regional areas per slice ____
# m: mouse object
# slice_ID : ID of slice (str)
# channels: channels to filter (str vector), can filter both eyfp & cfos
# exclude_regions: NULL (str vector); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
# exclude_hemisphere: TRUE (bool); excludes the contralateral hemisphere from one indicated in slice attribute
# exclude_layer_1: TRUE (bool); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)


get_registered_areas.slice <- function(s){



}


## excluded designated regions from the mapped cell data table
get_registered_areas <- function(x, ...){
  UseMethod('get_registered_areas')
}






