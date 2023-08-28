
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

import_segmentation_ij <- function(x, ...){
  UseMethod('import_segmentation_ij')
}


## Creating generic function for custom segmentation importation
#' import_segmentation (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export

import_segmentation_custom <- function(x, ...){
  UseMethod('import_segmentation_custom')
}


## Creating generic function for filtering raw segmentation data
#' make_segmentation_filter (generic function)
#'
#' @param x
#' @param ...
#'
#' @return
#' @export

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

get_registered_volumes <- function(x, ...){
  UseMethod('get_registered_volumes')
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

  # Detect the OS and set quartz( as graphing function)
  if(get_os() != "osx"){
    quartz <- X11
  }

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
#' @param slice_ID (str)  ID of slice
#' @param hemisphere  (str, default = NULL) 'left', 'right' or NULL if both hemispheres are included
#' @param filter (list) Wholebrain filter with parameters.
#' @param replace (bool, default = FALSE) Replace a registration already contained in a mouse object by resetting to NULL value before registration improvement loop.
#' @param ... additional parameters to pass to the SMART::registration2() function, besides 'input', 'coordinate', 'filter' & 'correspondance'
#' @return m  mouse object
#'
#' @examples m <- register(m, slice_ID = '1_10', hemisphere = "left", filter = my_filter)
#' @export
# TODO: add popup window option to better visualize the internal structures of the images
register.mouse <- function(m,
                           slice_ID = NA,
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
#' @rdname import_segmentation_ij
#' @description Method for importing segmentation data for a slice object
#' @param s slice object
#' @param mouse_ID (str) ID of mouse
#' @param channels (str vector, default = NULL) All the channels the user wants to import.
#' If NULL, defaults to the channels stored in the slice object attributes.
#' @return s slice object
#' @examples
#' s <-  import_segmentation(s, mouse_ID = "255") # Defaults to channels stored in slice attributes
#' s <-  import_segmentation(s, mouse_ID = "255", channels = c("cfos", "eyfp", "colabel")) # Specify channels
#' @export
#' @note The designated `colabel` channel name in this pipeline will auto import the output of the batch_3D_MultiColocalization.ijm macro provided in the pre-processing pipeline.
#' If you have a separate method used for detecting colabelled cells, please use a different naming convention for this channel,
#' e.g. "colabel_PV_cfos", and import using a customized import function such as [SMARTR::import_segmentation_custom()].

import_segmentation_ij.slice <- function(s,
                                      mouse_ID = NA,
                                      channels = NULL){


  # Get slice information
  info <- attr(s, "info")

  if (is.null(channels) && !is.null(info$channels)){
    channels <- info$channels
  } else if (is.null(info$channels)){
    stop("You must specify the channels parameter.")
  }

  ## Use registration image path to find segmentation data path stem
  # Check if slice has a set directory, otherwise borrow from the registration path
  if (is.null(info$slice_directory)){
    if (info$registration_path == 'set registration image path'){
      stop(paste0("This function doesn't know where to look to import segmentation data!\n",
                  "Set either the registration_path or the slice_directory in your slice attributes to continue."))
    } else {
      path_stem <- dirname(info$registration_path)
    }
  } else{
    path_stem <- info$slice_directory
  }

  section <- info$slice_ID
  txt_files <-  list.files(path = path_stem, all.files = TRUE, pattern = ".txt")

  for (k in 1:length(channels)){

    if (tolower(channels[k]) == 'colabel') {
      # Find the files necessary for colocalization analysis
      coloc_path <- paste0( mouse_ID,'_',section,"_ColocOnlyTable.txt")
      M_image_A_path <- paste0("M_", mouse_ID, '_', section,"_Coloc_A_C2.txt")
      M_image_B_path <- paste0("M_", mouse_ID, '_', section,"_Coloc_B_C1.txt")


      # Best match helper in case of user name edit errors
      coloc_path <- txt_files[stringdist::amatch(coloc_path, txt_files, maxDist=Inf)]
      M_image_A_path <- txt_files[stringdist::amatch(M_image_A_path, txt_files, maxDist=Inf)]
      M_image_B_path <- txt_files[stringdist::amatch(M_image_B_path, txt_files, maxDist=Inf)]

      message("Imported the following files: \n")
      print(coloc_path)
      print(M_image_A_path)
      print(M_image_B_path)

      # Read
      coloc.table <- read.delim(file.path(path_stem, coloc_path), stringsAsFactors = FALSE)
      image_A_objects <- read.csv(file.path(path_stem,  M_image_A_path), stringsAsFactors = FALSE)
      image_B_objects <- read.csv(file.path(path_stem,  M_image_B_path), stringsAsFactors = FALSE)

      # store the coloc table and the eyfp 16 bit measurements as a combined list
      s$raw_segmentation_data[[k]] <- list(coloc_table = coloc.table,
                                           image_A_objects = image_A_objects,
                                           image_B_objects = image_B_objects)
    }
      else{
        # Create import data paths for custom named channel
        channel <- tolower(channels[k])
        meas_path <- paste0("M_C_", channel, "_", mouse_ID,'_',section,".txt")
        quant_path <- paste0("Q_C_", channel, "_", mouse_ID,'_',section,"_", channel,".txt" )

        # match to exact name in directory
        meas_path <- txt_files[stringdist::amatch(meas_path, txt_files, maxDist=Inf)]
        quant_path <- txt_files[stringdist::amatch(quant_path, txt_files, maxDist=Inf)]
        message("Imported the following files: \n")
        print(meas_path)
        print(quant_path)

        # Read cell counts
        meas <- read.csv( file.path(path_stem, meas_path), stringsAsFactors = FALSE )
        quant <- read.csv(file.path(path_stem, quant_path), stringsAsFactors = FALSE)
        counts <- cbind(meas, quant) #create combined table
        counts <- counts[,unique(names(counts))]
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
#' @rdname import_segmentation_ij
#' @description Method for importing segmentation data for a mouse object
#' @param  m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str)'left', 'right' or NULL
#' @param channels (str vector, default = NULL) channels to import. If NULL, defaults to the channels stored in the slice object attributes.
#' @param replace (bool, default = FALSE) replace existing raw segmentation data
#' @return m mouse object
#' @examples m <-  import_segmentation(m, slice_ID = "1_10", channels = c("cfos", "eyfp", "colabel"), replace = FALSE)
#' @export


import_segmentation_ij.mouse <- function(m,
                                      slice_ID = NA,
                                      hemisphere = NULL,
                                      channels = NULL,
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
          m$slices[[match]] <- import_segmentation_ij(m$slices[[match]],
                                                   mouse_ID = mouse_ID,
                                                   channels = channels)
        } else{
          stop(paste0("There is existing segmentation data for this slice! If you want to overwrite it, set replace to TRUE."))
        }
      } else  {
        # Import segmentation with matched slices with no segmentation data yet
        m$slices[[match]] <- import_segmentation_ij(m$slices[[match]],
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



##____ Creating CUSTOM method for importing segmentation data for a slice object____
#' Import segmentation data
#' @rdname import_segmentation_custom
#' @description Custom method for importing segmentation data for a slice object.
#' This is a flexible method for importing channels aside from `cfos` and `eyfp` that use the same ImageJ segmentation macros common to the
#' Denny Lab. Other labs may use this method. This method works best when following saving segmentation data with the same naming conventions using the
#' \href{https://imagejdocu.list.lu/plugin/stacks/3d_ij_suite/start}{3D Roi Manager} within the
#' \href{dx.doi.org/10.1093/bioinformatics/btt276}{3D ImageJ Suite}.
#' Segmentation data is saved in two'.txt' files which output the results of the Measure 3D and Quantif 3D options of the 3D Roi Manager plugin,respectively.
#' A \href{https://imagejdocu.list.lu/tutorial/working/tutorial_for_3d_roi_manager}{description} of what each of those these options measure
#' is provided in the online documentation.
#'
#' The naming conventions for the ".txt" file storing the Quantif3D results are `Q_*_{channel}_*_{channel}.txt`,
#' where the `*` indicates a wildcard character(s) and {channel} is the channel name without brackets.
#' E.g. "Q_G_eYFP_258_1_1_eYFP.txt"
#'
#'
#' The naming conventions for the ".txt" file storing the Measure 3D are `M_*_{channel}_*.txt` where
#' the `*` indicates a wildcard character(s) and {channel} is the channel name without brackets.
#' E.g "M_G_eYFP_258_1_1.txt".
#'
#' The wildcards characters may be used to store things like date or slice naming information.
#'
#' The locations of these files must be specified in the `slice_directory` attribute of the slice_object.
#' Otherwise, the root folder containing the registration image is searched.
#' This attribute can be stored when initializing the slice object or can be edited afterwards.
#'
#' @param s slice object
#' @param channel (str) channel to import.
#' @param x_col (int, optional) The column index of the x pixel location in the txt file result from Measure 3D.
#' @param y_col (int, optional) The column index of the y pixel location in the txt file result from Measure 3D.
#' @param meas_path (chr, optional, default = NULL).
#' @param quant_path (chr, optional, default = NULL).
#' @return s slice object
#' @examples s <-  import_segmentation_custom(s, mouse_ID = "255", channel = "cfos")
#' @export


import_segmentation_custom.slice <- function(s,
                                             channel,
                                             x_col = NULL,
                                             y_col = NULL,
                                             meas_path = NULL,
                                             quant_path = NULL){

  # If the slice directory doesn't exist, it is set to the root directory of the registration image
  info <- attr(s, "info")
  if (is.null(info$slice_directory)){
    info$slice_directory <-  dirname(info$registration_path)
    attr(s, "info") <- info
  }

  if (is.null(quant_path) && is.null(meas_path)){
      # If file paths don't exist find them based off the root slice directory
    seg_f <- find_segmentation_files(info$slice_directory, channel)
    meas_path <- seg_f[1]
    quant_path <- seg_f[2]

  } else {
    if (!file.exists(quant_path) && !file.exists(meas_path)){
      warning('One of the files supplied in either the meas_path or quant_path arguments does not exist.\n
              Will search in the slice root directory for a matching file.')
      seg_f <- find_segmentation_files(info$slice_directory, channel)
      meas_path <- seg_f[1]
      quant_path <- seg_f[2]
    }
  }

  print(paste0("Importing ", meas_path))
  print(paste0("Importing ", quant_path))
  meas <- read.csv( file.path(info$slice_directory, meas_path), stringsAsFactors = FALSE)
  quant <- read.csv(file.path(info$slice_directory, quant_path), stringsAsFactors = FALSE)


  if (is.null(x_col)){
    meas$X2_Pix <- meas$CX..pix./(info$bin) #create position column to account for binning
  } else{
    meas$X2_Pix <- meas[, x_col]/(info$bin) #create position column to account for binning
  }

  if (is.null(y_col)){
    meas$Y2_Pix <- meas$CY..pix./(info$bin) #same as above
  } else{
    meas$Y2_Pix <- meas[, y_col]/(info$bin) #create position column to account for binning
  }

  counts <- cbind(meas, quant) #create combined table
  counts <- counts[,unique(names(counts))]

  # Add to current data
  n_cur <-  length(names(s$raw_segmentation_data))
  names_cur <- names(s$raw_segmentation_data)
  if (n_cur > 0){
    k <-  n_cur + 1
    s$raw_segmentation_data[[k]] <- counts
    names(s$raw_segmentation_data) <- c(names_cur, channel)
  } else {
    s$raw_segmentation_data[[1]] <- counts
    names(s$raw_segmentation_data) <- channel
  }
  return(s)
}

##____ Creating method for importing segmentation data for a mouse object____
#' Import raw segmentation data
#'
#' @rdname import_segmentation_custom
#' @description Method for custom importation of segmentation data for a mouse object
#' @param  m mouse object
#' @param channel (str) channel to import
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str)'left', 'right' or NULL
#' @param x_col = NULL,
#' @param y_col = NULL,
#' @param meas_path = NULL,
#' @param quant_path = NULL
#' @return m mouse object
#' @examples m <-  import_segmentation(m, slice_ID = "1_10", channels = c("PV"), replace = FALSE)
#' @export


import_segmentation_custom.mouse <- function(m,
                                             channel,
                                             slice_ID = NA,
                                             hemisphere = NULL,
                                             x_col = NULL,
                                             y_col = NULL,
                                             meas_path = NULL,
                                             quant_path = NULL){
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
          warning(paste0("There is existing segmentation data for this slice! Adding on top of data."))
      }

      m$slices[[match]] <- import_segmentation_custom(m$slices[[match]],
                                                      x_col = x_col,
                                                      y_col = y_col,
                                                      meas_path = meas_path,
                                                      quant_path = quant_path,
                                                      channel = channel)
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!"))
  }
  return(m)
}

##____Optional Method for making an eyfp filter____

#' Make a segmentation filter for a slice object
#' @rdname make_segmentation_filter
#' @param s
#' @param channels (str vector, default = "eyfp") Channels to process. If NULL, defaults to the channels stored in the slice object attributes (not recommended).
#' @param params (list) Has same length as channels. Each element contains a vector of parameters names used for filtering that channel
#' @param ranges (list of lists) Has same length as channels. Each element of outer list corresponds the order of channels you want to process.
#' Inner list contains vectors of parameter ranges for that channel.
#'
#' @return s slice object. Vector of indices of cells to remove are stored as the channel filters in the slice object.
#' @export
#'
#' @examples s <- make_segmentation_filter(s, channels = c('eyfp'),
#' params = list(c("Vol..unit.","Moment1","Moment2","Moment3","Moment4","Sigma")),
#' ranges = list(list(c(200, 12000), c(3, 50), c(0, 600), c(0, 2000), c(0, 5), c(20, Inf))))

make_segmentation_filter.slice <- function(s,
                                           channels = "eyfp",
                                           params = list(c("Vol..unit.","Moment1","Moment2", "Moment3","Moment4","Sigma")),
                                           ranges = list(list(c(200,12000),
                                                              c(3,50),
                                                              c(0,600),
                                                              c(0,2000),
                                                              c(0,5),
                                                              c(20,Inf)))){
  info <- attr(s, 'info')
  if (is.null(channels) && !is.null(info$channels)){
    channels <- info$channels
  } else if (is.null(info$channels)){
    stop("You must specify the channels parameter.")
  }

  for (channel in channels){
    if (channel %in% names(s$raw_segmentation_data)){
      data <- s$raw_segmentation_data[[channel]]
    } else {
      stop(paste0("The ", channel, " channel has not been imported yet!"))
    }

    # Make the filter
    filter <- make.filter(data, params = params[[k]], ranges = ranges[[k]])
    s$segmentation_filters[[channel]] <- filter
    message('Processed ', channel, " channel")
  }
  return(s)
}


##____ Method for creating filter for making eyfp filter for mouse____

#' Make a segmentation filter for a slice within a mouse object
#' @rdname make_segmentation_filter
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere 'left', 'right' or NULL (both)
#' @param channels (str vector, default = "eyfp") Channels to process. If NULL, defaults to the channels stored in the slice object attributes (not recommended).
#' @param params  (list) Has same length as channels. Each element contains a vector of parameters names used for filtering that channel
#' @param ranges (list of lists) Has same length as channels. Each element of outer list corresponds the order of channels you want to process.
#' Inner list contains vectors of parameter ranges for that channel.
#' @param replace (bool, default = FALSE) replace existing filters.

#'
#' @return m mouse object.  Vector of indices of cells to remove are stored as the channel filters in the slice object within the mouse.
#' @export
#'
#' @examples m <- make_segmentation_filter(m, slice_ID = '1_10', hemisphere = NULL , channels = c('eyfp'),
#' params = list(c("Vol..unit.","Moment1","Moment2","Moment3","Moment4","Sigma")),
#' ranges = list(list(c(200, 12000), c(3, 50), c(0, 600), c(0, 2000), c(0, 5), c(20, Inf))))
make_segmentation_filter.mouse <- function(m,
                                           slice_ID = NA,
                                           hemisphere = NULL,
                                           channels = c('eyfp'),
                                           params =list(c("Vol..unit.","Moment1","Moment2","Moment3","Moment4","Sigma")),
                                           ranges = list(list(c(200,12000),c(3,50),
                                                              c(0,600),c(0,2000),c(0,5),c(20,Inf))),
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
      if (!is.null(m$slices[[match]]$segmentation_filters)){

        if (replace){
          # replace the slice's raw segmentation data
          m$slices[[match]]$segmentation_filters <- NULL

          # Import segmentation from scratch
          m$slices[[match]] <- make_segmentation_filter(m$slices[[match]],
                                                        channels = channels,
                                                        params = params,
                                                        ranges = ranges)
        } else{
          stop(paste0("A segmentation filter already exists for this slice. Set replace = TRUE to overwrite the current filter(s). "))
        }
      } else  {
        # Import segmentation with matched slices with no segmentation data yet
        m$slices[[match]] <- make_segmentation_filter(m$slices[[match]],
                                                      channels = channels,
                                                      params = params,
                                                      ranges = ranges)
      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!"))
  }
  return(m)
}


##____Method for creating segmentation object for each slice____


#' Creates a wholebrain segmentation object and stores it.
#'
#' @rdname make_segmentation_object
#' @description Make a wholebrain compatible segmentation object for a slice in a slice object
#'
#' @param s slice object
#' @param mouse_ID (str) ID of mouse
#' @param channels (str vector, default = NULL) Channels to process. If NULL, defaults to the channels stored in the slice object attributes.
#' @param use_filter (bool, default = FALSE). Use a filter to create more curated segmentation object from the raw segmentation data.
#' @param ... (optional) additional volume and overlap parameters for get.colabeled.cells().
#'
#' @return s slice object
#' @examples s <- make_segmentation_object(s, mouse_ID = "255", channels = c("cfos", "eyfp"), use_filter = FALSE)
#' @export
#' @note If you are processing the colabel channel, the X and Y positions of colabelled cells are the average of the x,y centroid coordinates of the colabelled objects

make_segmentation_object.slice <- function(s,
                                           mouse_ID = NA,
                                           channels = NULL,
                                           use_filter = FALSE,
                                           ... ){

  # Get slice information
  info <- attr(s, 'info')
  if (is.null(channels) && !is.null(info$channels)){
    channels <- info$channels
  } else if (is.null(info$channels)){
    stop("You must specify the channels parameter.")
  }

  # Check if colabel channel is included (MAY NOT BE NECESSARY)
  if ('colabel' %in% channels){
    # Rearrange the channels vector so that colabel is always processed last
    channels <- c(channels[!channels %in% 'colabel'], channels[channels %in% 'colabel'])
  }

  # Create a segmentation list to store the channel segmentation objects
  segmentation_list  <- vector(mode = 'list', length = length(channels))
  names(segmentation_list) <- names(channels)

  for (channel in channels){
    if (channel == 'colabel'){
      seg.coloc <- get.colabeled.cells(s$raw_segmentation_data$colabel$coloc_table,
                                              s$raw_segmentation_data$colabel$image_A_objects,
                                              s$raw_segmentation_data$colabel$image_B_objects,
                                              ...)
      seg.coloc$soma$x <- seg.coloc$soma$x/(info$bin)
      seg.coloc$soma$y <- seg.coloc$soma$y/(info$bin)
      segmentation_list[[channel]] <- seg.coloc
    } else {
      #Other channels
      if (use_filter){
        if (is.null(s$segmentation_filter[[channel]])){
          warning(paste0("No filter is available for ", channel," even though `use_filter` parameter is TRUE.\nThe channel will be processed without usage of an filter."))
          counts <- s$raw_segmentation_data[[channel]]
        } else {
          counts <- s$raw_segmentation_data[[channel]][-s$segmentation_filters[[channel]], ]
        }
      } else {
        counts <- s$raw_segmentation_data[[channel]]
      }
      seg <- SMARTR::segmentation.object
      seg$soma$x <- counts$CX..pix./(info$bin) #create position column to account for binning
      seg$soma$y <- counts$CY..pix./(info$bin) #same as above
      seg$soma$area <- counts$Vol..pix.
      seg$soma$intensity <- counts$Mean
      # Store into segmentation list
      segmentation_list[[channel]] <- seg
    }
  }
 s$segmentation_obj <- segmentation_list
  return(s)
}

##____Method for creating segmentation object for each mouse____

#' Creates a wholebrain segmentation object and stores it
#' @rdname make_segmentation_object
#' @description Make a wholebrain compatible segmentation object for a slice in a mouse object
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str, default = NULL) 'left', 'right' or NULL (both)
#' @param channels (str vector, default = NULL) Channels to process. If NULL, defaults to the channels stored in the slice object attributes.
#' @param replace (bool, default = FALSE) replace existing raw segmentation data
#' @param use_filter (bool, default = FALSE). Use a filter to create more curated segmentation object from the raw segmentation data.
#' @param ... (optional) additional volume and overlap parameters for get.colabeled.cells().
#' @examples m <-  make_segmentation_object(m, slice_ID = '1_9', hemisphere = 'left', channels = c('eyfp', 'cfos', 'colabel'), use_filter = FALSE)
#' @export
#'
make_segmentation_object.mouse <- function(m,
                                           slice_ID = NA,
                                           hemisphere = NULL,
                                           channels = NULL,
                                           replace = FALSE,
                                           use_filter = FALSE,
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

      # Check if segmentation object data already exists
      if (!is.null(m$slices[[match]]$segmentation_obj)){

        if (replace){
          # replace the slice's segmentation object
          m$slices[[match]]$segmentation_obj <- NULL

          # Make segmentation object from scratch
          m$slices[[match]] <- make_segmentation_object(m$slices[[match]],
                                                        mouse_ID = mouse_ID,
                                                        channels = channels,
                                                        use_filter = use_filter,
                                                        ...)
        } else{
          stop(paste0("There is an existing segmentation object for this slice! If you want to overwrite it, set replace to TRUE."))
        }
      } else  {

        # Import segmentation with matched slices with no segmentation data yet
        m$slices[[match]] <- make_segmentation_object(m$slices[[match]],
                                                      mouse_ID = mouse_ID,
                                                      channels = channels,
                                                      use_filter = use_filter,
                                                      ... )
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
#' @param s slice object
#' @param channels (str vector, default = NULL) Channels to process. If NULL, defaults to the channels stored in the slice object attributes.
#' @param clean (bool, default = TRUE). Remove cells that don't map to any regions.
#' @param display (bool, default = TRUE). Display the results of the forward warp for the slice.
#' @param mouse_ID (str) mouse ID
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the [wholebrain::inspect.registration()] function
#' @examples s <- map_cells_to_atlas(s, channels c('cfos' , 'eyfp', 'colabel'), clean = TRUE, display = TRUE, mouse_ID = "255")
#' @export


map_cells_to_atlas.slice <- function(s,
                                     channels = NULL,
                                     clean =  TRUE,
                                     display = TRUE,
                                     mouse_ID = NULL,
                                     ...){

  # Get slice information
  info <- attr(s, 'info')
  if (is.null(channels) && !is.null(info$channels)){
    channels <- info$channels
  } else if (is.null(info$channels)){
    stop("You must specify the channels parameter.")
  }

  for (channel in channels){
    cell_data <- wholebrain::inspect.registration(s$registration_obj,
                                      s$segmentation_obj[[channel]],
                                      forward.warps = TRUE,
                                      device = display,
                                      ...)
    # s$raw_forward_warped_data[[channel]] <- cell_data # Deprecate this, it's unnecessary

    if(clean){
      cell_data <- cell_data[!cell_data$id==0,]
    }

    if (!is.null(mouse_ID)){
      cell_data$animal <-  mouse_ID
    }

    s$forward_warped_data[[channel]] <- cell_data
  }
  return(s)
}

##____Method for forward warping the segmentation data per mouse____
#' @title Map cells to atlas for slice within a mouse object
#' @description  Method for forward warping segmentation data to atlas space for a slice within a mouse object.
#' Requires segmentation objects to be made for the channels specified and a registration.
#' @rdname map_cells_to_atlas
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str) 'left', 'right' or NULL (both)
#' @param channels (str vector, default = NULL) Channels to process. If NULL, defaults to the channels stored in the slice object attributes.
#' @param clean (bool, default = TRUE). Remove cells that don't map to any regions.
#' @param display (bool, default = TRUE). Display the results of the forward warp for the slice.display
#' @param replace (bool, default = FALSE). Replace current forward warped data, both raw and cleaned.
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the [wholebrain::inspect.registration()] function
#' @return m mouse object
#' @examples m <- map_cells_to_atlas(m, slice_ID = "1_10", hemisphere = NULL, channels = c("cfos", "eyfp", "colabel"), clean = TRUE, replace = FALSE)
#' @export

map_cells_to_atlas.mouse <- function(m,
                                     slice_ID = NA,
                                     hemisphere = NULL,
                                     channels = NULL,
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
          stop(paste0("There is a existing mapped data for this slice! If you want to overwrite it, set replace to TRUE."))
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
#'
#' @param s slice object
#' @param channels (str vector, default = NULL) Channels to process. If NULL, defaults to the channels stored in the slice object attributes.
#' @param clean (bool, default = TRUE ). Remove cells that don't map to any regions. This option is recommended.
#' @param exclude_right_regions (str vector, default = NULL); acronyms of regions you want to exclude from right hemi,in addition to regions that will by default be excluded in the slice attribute 'right_regions_excluded'
#' @param exclude_left_regions (str vector, default = NULL); acronyms of regions you want to exclude from left hemi, in addition to regions that will by default be excluded in the slice attribute 'left_regions_excluded'
#' @param exclude_hemisphere (bool, default = TRUE); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1 (bool, default = TRUE) excludes all counts from layer 1
#' @param include_right_regions (str vector, default = NULL) Acronyms of regions to include from the right hemi; if not NULL, takes precedence over `exclude_right_regions` & all other regions will be excluded. Typically, this is used for slices with poor quality/lots of tears.
#' @param include_left_regions (str vector, default = NULL) Acronyms of regions to include from the light hemi; if not NULL, takes precedence over `exclude_left_regions` & all other regions will be excluded.  Typically, this is used for slices with poor quality/lots of tears.
#' @param plot_filtered (bool, default = TRUE) pop up window to check the excluded anatomy.
#'
#' @examples s <-  exclude_anatomy(s, channels = c('cfos', 'eyfp', 'colabel'), clean = TRUE, exclude_regions = NULL, exclude_hemisphere = TRUE, exclude_layer_1 = TRUE, plot_filtered = TRUE)
#' @export

exclude_anatomy.slice <- function(s,
                                  channels = NULL,
                                  clean = TRUE,
                                  exclude_right_regions = NULL,
                                  exclude_left_regions = NULL,
                                  exclude_hemisphere = TRUE,
                                  exclude_layer_1 = TRUE,
                                  include_right_regions = NULL,
                                  include_left_regions = NULL,
                                  plot_filtered = TRUE){


  # Get slice information & combine exclude regions parameter with default regions excluded;
  info <- attr(s, 'info')
  if (is.null(channels) && !is.null(info$channels)){
    channels <- info$channels
  } else if (is.null(info$channels)){
    stop("You must specify the channels parameter.")
  }

  include_right_regions <- unique(c(include_right_regions, info$right_regions_included))
  include_left_regions <- unique(c(include_left_regions, info$left_regions_included))


  # if (!is.null(include_right_regions)  ){
  if (!(length(include_right_regions) < 1)){
    # The regions to include parameter is being used. Exclude all other regions
    # Concatenate list of all regions
    include_right_regions <-  find_all_subregions(include_right_regions)
    attr(s, 'info')$right_regions_included <- include_right_regions
  } else{
    # Reassign so that all regions excluded on the right side are tracked
    exclude_right_regions <- unique(c(exclude_right_regions, info$right_regions_excluded))
    attr(s, 'info')$right_regions_excluded <- exclude_right_regions
    all_excluded_right_regions <-  find_all_subregions(exclude_right_regions)

  }



  # if (!is.null(include_left_regions)){
  if (!(length(include_left_regions) < 1)){
    # The regions to include parameter is being used. Exclude all other regions
    # Concatenate list of all regions
    include_left_regions <-  find_all_subregions(include_left_regions)
    attr(s, 'info')$left_regions_included <- include_left_regions
  } else{
    # Reassign so that all regions excluded on the left side are tracked
    exclude_left_regions <- unique(c(exclude_left_regions, info$left_regions_excluded))
    attr(s, 'info')$left_regions_excluded <- exclude_left_regions
    all_excluded_left_regions <-  find_all_subregions(exclude_left_regions)
  }




  ## Filtering per channel
  for (channel in channels){
    dataset <- s$forward_warped_data[[channel]]

    # 1) Filter out right and left regions

    # if (!is.null(include_right_regions)){
    if (!(length(include_right_regions)<1)){
      dataset <-  dataset[!dataset$right.hemisphere | (dataset$right.hemisphere & (dataset$acronym %in% include_right_regions)),]
    } else {
      dataset <- dataset[!(dataset$right.hemisphere & (dataset$acronym %in% all_excluded_right_regions)),]
    }

    # 1) Filter out left and left regions
    # if (!is.null(include_left_regions)){
    if (!(length(include_left_regions)<1)){
      dataset <-  dataset[dataset$right.hemisphere | (!dataset$right.hemisphere & (dataset$acronym %in% include_left_regions)),]
    } else {
      dataset <- dataset[!(!dataset$right.hemisphere & (dataset$acronym %in% all_excluded_left_regions)),]
    }


    # 2) Filter out hemisphere
    if (exclude_hemisphere){
      if (info$hemisphere == "right"){
        dataset <- dataset[dataset$right.hemisphere,]
      } else if (info$hemisphere == "left"){
        dataset <- dataset[!dataset$right.hemisphere,]
      }
    }

    # 3) Filter out layer 1 of the Cortex
    if (exclude_layer_1){
      dataset <- dataset[-grep("layer 1",dataset$name, ignore.case = TRUE, value = FALSE),]
    }

    # Filter out cell counts that are out of bounds
    if(clean){
      dataset <- dataset[!dataset$id==0,]
      dataset <- tidyr::drop_na(dataset) # Get rid of rows that contain any NA values
    }

    # Dataset
    s$forward_warped_data[[channel]] <- dataset

    ## plotting schematic plots to check that the region cell counts are cleaned
    if (plot_filtered){
      wholebrain::schematic.plot(dataset, device = TRUE)
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
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str) 'left', 'right' or NULL (both)
#' @param channels (str vector, default = NULL) Channels to process. If NULL, defaults to the channels stored in the slice object attributes.
#' @param clean (bool, default = TRUE ). Remove cells that don't map to any regions.
#' @param exclude_right_regions (str vector, default = NULL); acronyms of regions you want to exclude from right hemi,in addition to regions that will by default be excluded in the slice attribute 'right_regions_excluded'
#' @param exclude_left_regions (str vector, default = NULL); acronyms of regions you want to exclude from left hemi, in addition to regions that will by default be excluded in the slice attribute 'left_regions_excluded'
#' @param exclude_hemisphere (bool, default = TRUE); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1 (bool, default = TRUE); excludes all counts from layer 1 (TEMPORARY, may not be hardcoded in later)
#' @param plot_filtered (bool, default = TRUE) pop up window to check the excluded anatomy.
#' @return m mouse object
#' @examples m <- exclude_anatomy(m, slice_ID = "1_10", hemisphere = NULL, channels = c('cfos', 'eyfp', 'colabel'), clean = TRUE,
#' exclude_regions = NULL,
#' exclude_hemisphere = TRUE,
#' exclude_layer_1 = TRUE
#' @export


 exclude_anatomy.mouse <- function(m,
                                  slice_ID = NA,
                                  hemisphere = NULL,
                                  channels = NULL,
                                  clean = TRUE,
                                  exclude_right_regions = NULL,
                                  exclude_left_regions = NULL,
                                  exclude_hemisphere = FALSE,
                                  exclude_layer_1 = TRUE,
                                  plot_filtered = TRUE){

  # Omit hemisphere name if there is no hemisphere value in the attributes
  if (is.null(hemisphere)){
    slice_name <- slice_ID
    if (exclude_hemisphere){
      stop("Can't exclude contralateral hemisphere without reference hemisphere! i.e., hemisphere = NULL and exclude_hemisphere = TRUE. ")
    }

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
      if (!is.null(m$slices[[match]]$forward_warped_data)){

        # If forward mapped data exists, exclude the designated areas
        m$slices[[match]] <- exclude_anatomy(m$slices[[match]],
                                             channels = channels,
                                             clean = clean,
                                             exclude_right_regions = exclude_right_regions,
                                             exclude_left_regions = exclude_left_regions,
                                             exclude_hemisphere = exclude_hemisphere,
                                             exclude_layer_1 = exclude_layer_1,
                                             plot_filtered = plot_filtered)
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



# ##____Method for getting regional areas and volumes per slice ____
#' Method for getting regional areas and volumes per slice
#' @rdname get_registered_volumes
#' @description Calculate the registered area (microns^2^) and the regional volumes (microns^3^) of all regions contained in a slice.
#' @param s slice object
#' @return s slice object with a stored dataframe with columns 'name' (full region name), 'acronym', 'area' (in microns^2^), 'volume' (in microns^3^) 'right.hemisphere'
#' @examples s <- get_registered_areas(s)
#' @export
#' @md

get_registered_volumes.slice <- function(s){

  # Get conversion factor from pixel to microns
  # get z-width of this slice (in microns)
    cf <- attr(s, 'info')$conversion_factor
    z_width <- attr(s, 'info')$z_width

    s$volumes <- get.registered.areas(s$forward_warped_data, registration = s$registration_obj, conversion.factor = cf) %>%
      dplyr::mutate(volume = area*z_width) %>% tidyr::drop_na()
    return(s)
}


## ___ Method for getting regional areas and volumes per slice
#' Method for getting regional areas and volumes for each slice in a mouse object
#' @rdname get_registered_areas
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str) 'left', 'right' or NULL (both)
#' @return m mouse object
#' @examples m <- get_registered_areas(m, slice_ID = "1_10", hemisphere = "left", replace = FALSE)
#' @export

get_registered_volumes.mouse <- function(m,
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

      # Check if volumes data already exists
      if (!is.null(m$slices[[match]]$volumes)){

        # Check if user wants existing data to be overwritten
        if (replace){
          m$slices[[match]] <- get_registered_volumes(m$slices[[match]])

        } else {
          stop(paste0("There is an existing volume data for this slice! If you want to overwrite it, set replace to TRUE."))
        }
      } else {
        m$slices[[match]] <- get_registered_volumes(m$slices[[match]])
      }
    }
  }
  if (isFALSE(match)){
    stop(paste0("There were no slices matching the name ", slice_name, " found in your mouse!" ))
  }
return(m)
}



#__________________ slice object specific functions __________________________



#' Adjust brain outline.
#' @description  This function takes a slice object and first applies a filter with default settings to
#' the image set as the slice registration path. An interative user loop allows for easy adjustment of the brain threshold
#' since the wholebrain GUI tends to be a bit buggy. This function then returns a filter with the adjusted brain threshold.
#'
#' @param s slice object
#' @param filter (list, default = NULL) If the user passes their own filter list, it will use that instead of the presaved filter list from SMARTR.
#'
#' @return filter (list) wholebrain compatible filter
#' @export
#'
#' @examples filter <- adjust_brain_outline(s); s <- register(s, filter = filter) Adjust the brain threshold, then run register on the slice object
adjust_brain_outline <- function(s, filter = NULL){

  regi_path <- attr(s, "info")$registration_path

  if (is.null(filter)){
    filter <- SMARTR::filter
    filter$brain.threshold <- 10
  } else if (!is.list(filter)){
    stop("You did not supply a valid list for the filter parameter.")
  }

  cat(paste0("Trying default brain threshold of ", filter$brain.threshold, "\n" ))
  filter <- wholebrain::segment(regi_path, filter = filter)
  filter <- filter$filter

  change_done <-FALSE
  while (!change_done) {
    cat(paste0("Your brain threshold is: ", filter$brain.threshold ))
    inp <- readline("Do you want to change your threshold: Y/N?" )
    if (inp=="Y" || inp=="y") {

      filter$brain.threshold <- as.integer(readline("Enter your new brain threshold: "))
      filter <- wholebrain::segment(regi_path, filter = filter)
      filter <- filter$filter

    } else if ( inp=="N" || inp == "n") {
      # exit out of loop
      change_done  <- TRUE
    }
  }

  return(filter)
}



#__________________ Mouse object specific functions __________________________


#' Get cell tables
#' @description This function stores a list in a mouse object that is the length of the channels parameter. Each element of the list is a dataframe
#' containing combined cell counts for one channel across all slices processed for that mouse.
#' By default, if one slice has no dataset processed for that particular channel, that slice will be skipped over.
#' The function will run properly but a warning will be thrown indicating you should go back and
#' generate a mapped dataset for that particular slice and channel.
#' @param m mouse object
#' @param channels (vec, default = c("cfos", "eyfp", "colabel")) Channels to process.
#' @return m mouse object
#' @examples m <- get_cell_tables(m, channels = c("cfos", "eyfp", "colabel"))
#' @export

get_cell_table <- function(m,
                           channels = c("cfos", "eyfp", "colabel")){

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

      # Add the slice name to the dataset dataframe
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
#' Normalize cell counts per mm^2^ or by mm^3^ (if multiplying by the stack size).
#' @description Run this function after all the slices that you want to process are finished being added
#' and you have combined your cell counts with [SMARTR::get_cell_table()]. This functions process all channels
#' where a cell table was made using the latter function.
#' @param m mouse object
#' @param combine_hemispheres (bool, default = TRUE) Combine normalized cell counts from both hemispheres
#' @param simplify_regions (bool, default = TRUE ) simplify the normalized region counts based on keywords in the internal function, `simplify_keywords`
#' @param simplify_keywords (str vec, default =  c("layer","part","stratum","division")). Keywords to search through region names and simplify to parent structure
#' @return
#' @examples m <- normalize_cell_counts(m, combine_hemispheres = TRUE, simplify_regions = TRUE)
#' @export
#' @ms
#'
# TODO: Change the internal implementation of plyr::ddply to a function consistent with dplyr
#
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

      ## Check that this slice has run get_registered_volumes
      if (is.null(s$volumes)){
        stop(paste0("Slice ", slice_name ," in your mouse dataset has an empty areas vector.\n",
                    "Run the function get_registered_volumes() for this slice before you can normalized all your cell counts by total area or by volume."))
        }
    }

  # ____________________________________________________________________________
  # 2)  Getting total areas across all slices, separated by hemisphere

  # rbind all the areas table
  aggregate_volumes <- m$slices[[1]]$volumes
  if (length(m$slices)  >=  2){
    for (k in 2:length(m$slices)){
      aggregate_volumes <- rbind(aggregate_volumes, m$slices[[k]]$volumes)
    }
  }

  # Summarize by acronym and right.hemisphere
  total_volumes <- dplyr::group_by(aggregate_volumes, acronym, right.hemisphere, name) %>%
    dplyr::summarise(area.mm2 = sum(area)*1e-6, volume.mm3 = sum(volume)*1e-9)
  m$total_volumes <- total_volumes

  # ____________________________________________________________________________
  # 3)
  # then use Cell_table to create combined cell counts normalized by area
  normalized_counts <- normalize.registered.areas(m$cell_table, total_volumes) #tabulate number of cells in each region
  names(normalized_counts) <- names(m$cell_table)

  # Clean NAs, check for combining hemispheres or simplifying regions
  for (channel in names(normalized_counts)){

    # clean NA values in all channels
    normalized_counts[[channel]] <- tidyr::drop_na(normalized_counts[[channel]])

    # Combine hemispheres
    if (combine_hemispheres){
      # Collapse the hemisphere data and recalculating normalized counts by area and volume
      normalized_counts[[channel]] <-  dplyr::group_by(normalized_counts[[channel]], acronym, name) %>%
        dplyr::summarise_at(c("count", "area.mm2", "volume.mm3"), sum) %>%
        dplyr::mutate(normalized.count.by.area = count/area.mm2,
                      normalized.count.by.volume = count/volume.mm3)
    }
  }

  # Simplify regions
  if (simplify_regions){
    # normalize the cell counts by specific keyworks
    normalized_counts <- simplify.regions(normalized_counts, keywords = simplify_keywords)
  }

  # Store the normalized counts in the mouse object
  m$normalized_counts <- normalized_counts
  return(m)
}



## Auto split the HPF dataset into Dorsal and Ventral

#' @title Split the hippocampal dataset to dorsal and ventral regions.
#' @description If [normalize_cell_counts()] has already been run,
#' the existing hippocampus data will be deleted and either replaced with new dorsal/ventral hippocampus data when merge = TRUE (recommended) or
#' it will be stored in a separate dataframe if the user wants analyze the hippocampus separately for analysis and the `normalized_counts` dataframe will be unchanged.
#' @param m mouse object
#' @param AP_coord (double) The AP coordinate which separates dorsal hippocampus from ventral hippocampus. Counts anterior to this are considered dorsal, counts posterior to
#' this are considered ventral.
#' @param combine_hemispheres (bool, default = TRUE) Combine normalized cell counts from both hemispheres
#' @param merge (bool, default = TRUE) Recommended setting. Whether to include the division of dHipp and vHipp into the normalized counts dataframe of the mouse or to store as a separate element for isolated analysis.
#' Merging will replace the hippocampal counts currently in the normalized counts table. Not merging will create a separate data element in the mouse. TODO: Currently not functional.
#' @param simplify_regions (bool, default = TRUE ) simplify the normalized region counts based on keywords in the internal function, `simplify_keywords`
#' @param simplify_keywords (str, default =  c("layer","part","stratum","division")). Keywords to search through region names and simplify to parent structure
#' @param rois (str, default = c("DG", "CA1", "CA2", "CA3")) List of the region acronyms to include in the hippocampus
#' @return m mouse object with DV cell counts from the hippocampus included
#' @export
#'
#' @examples
split_hipp_DV <- function(m,
                          AP_coord = -2.7,
                          combine_hemispheres = TRUE,
                          merge = TRUE,
                          simplify_regions = TRUE,
                          simplify_keywords = c("layer","part","stratum","division"),
                          rois = c("DG", "CA1", "CA2", "CA3")
                          ){


  # Channels
  channels <- names(m$cell_table)

  # Loop through all rois and get all subregions
  regions <- c(rois)
  for (roi in rois) {
    regions <- c(regions, SMARTR::get.sub.structure(roi))
  }

  ## Isolating counts from the hippocampus & sorting them by coordinate
  dHipp_cell_table <- vector(mode = "list", length = length(channels))
  vHipp_cell_table <- vector(mode = "list", length = length(channels))
  names(dHipp_cell_table) <- channels
  names(vHipp_cell_table) <- channels

  # process each channel separately
  hipp_AP_coordinates <- c()
  for (channel in channels){

    # Get cell table of just the hippocampus
    hipp_df <- m$cell_table[[channel]] %>% dplyr::filter(acronym %in% regions)

    # Store dorsal counts if there are any, else erase
    if (sum(hipp_df$AP > AP_coord) > 0){
      dHipp_cell_table[[channel]] <- hipp_df[hipp_df$AP > AP_coord, ]
    } else {
      dHipp_cell_table[[channel]] <- NULL
    }

    # Store ventral counts if there are any, else erase
    if (sum(hipp_df$AP <= AP_coord) > 0){
      vHipp_cell_table[[channel]] <- hipp_df[hipp_df$AP <= AP_coord, ]
    } else {
      vHipp_cell_table[[channel]] <- NULL
    }

    # Keep track of the unique AP coordinate (slice identifier)
    hipp_AP_coordinates <- c(hipp_AP_coordinates, hipp_df$AP) %>% unique()

    # Remove existing hippocampus data from the normalized cell counts table if it exists
    if (!is.null(m$normalized_counts) & merge){
      hipp_delete <- which(m$normalized_counts[[channel]]$acronym %in% regions)
      m$normalized_counts[[channel]] <- m$normalized_counts[[channel]][-hipp_delete,]

      message(paste0("Removed existing hippocampus data for ", channel,
                     " channel from the normalized counts dataframe"))
    }

  }

  # Get hippocampal areas
  DV_volumes <- get_hipp_DV_volumes(m, AP_coord = AP_coord, rois = rois, hipp_AP_coordinates)

  # Create combine cell_table and initialize list for normalized hippocampus cell counts
  hipp_cell_tables <- list(dHipp_cell_table,vHipp_cell_table)
  hipp_norm_counts <- vector(mode = "list", length = 2)
  names(hipp_cell_tables) <- c("dorsal", "ventral")
  names(hipp_norm_counts) <- c("dorsal", "ventral")

  # Get normalized cell counts
  for (dv in names(hipp_cell_tables)){
    if(length(hipp_cell_tables[[dv]]) > 0){
      hipp_norm_counts[[dv]] <- normalize.registered.areas(hipp_cell_tables[[dv]], DV_volumes[[dv]])

      # Simplify the subregions regions
      if (simplify_regions){
        hipp_norm_counts[[dv]] <- suppressWarnings(simplify.regions(hipp_norm_counts[[dv]],
                                                                    keywords = simplify_keywords))
        message("Simplified regions.")
      }

      for (channel in names(hipp_norm_counts[[dv]])){

        # clean NA values in all channels
        hipp_norm_counts[[dv]][[channel]] <- tidyr::drop_na(hipp_norm_counts[[dv]][[channel]])

        # Rename the regions to match regions of the dorsal and ventral hippocampus
        if (dv == "dorsal"){
          dorsal_acronyms <- paste0("d", hipp_norm_counts[[dv]][[channel]]$acronym)
          # print(dorsal_acronyms)
          hipp_norm_counts[[dv]][[channel]]$acronym <- dorsal_acronyms
          hipp_norm_counts[[dv]][[channel]]$name <- SMARTR::name.from.acronym(dorsal_acronyms) %>% as.character()

        } else {

          ventral_acronyms <- paste0("v", hipp_norm_counts[[dv]][[channel]]$acronym)
          # print(ventral_acronyms)
          hipp_norm_counts[[dv]][[channel]]$acronym <- ventral_acronyms
          hipp_norm_counts[[dv]][[channel]]$name <- SMARTR::name.from.acronym(ventral_acronyms) %>% as.character()
        }

        # combine hemispheres
        if (combine_hemispheres){
        hipp_norm_counts[[dv]][[channel]] <-  dplyr::group_by(hipp_norm_counts[[dv]][[channel]], acronym, name) %>%
          dplyr::summarise_at(c("count", "area.mm2", "volume.mm3"), sum) %>%
          dplyr::mutate(normalized.count.by.area = count/area.mm2,
                        normalized.count.by.volume = count/volume.mm3)
          message(paste0("Combined hemispheres for ", channel, " channel."))
        }
      }
    } else{
      hipp_norm_counts[[dv]] <- NULL
    }
  }


  suppressWarnings(
  # Whether to merge the information in each channel to the large normalized counts dataframe
  if (merge){
    # merged right.hemisphere check
    if (is.null(m$normalized_counts[[1]]$right.hemisphere) && is.null(hipp_norm_counts[[1]][[1]]$right.hemisphere) ||
        !is.null(m$normalized_counts[[1]]$right.hemisphere) && !is.null(hipp_norm_counts[[1]][[1]]$right.hemisphere)){

      ## continue with merge
      for (dv in names(hipp_norm_counts)){
        for (channel in channels){
          m$normalized_counts[[channel]] <- rbind(m$normalized_counts[[channel]], hipp_norm_counts[[dv]][[channel]]) %>% arrange(acronym)
        }
      }
      message("Merged into main normalized counts dataframe.")

    } else {
      stop(paste0("If you want to merge the dorsal ventral split dataset into the main normalized counts dataframe, the \n",
                  "combine_hemispheres parameter needs to be identical when you run normalize_cell_counts() and split_hipp_DV()."))
    }
  } else {
    m$hipp_DV_normalized_counts <- hipp_norm_counts
    message("Data is now stored in hipp_DV_normalized_counts dataframe. Future analysis of the hippocampus will be stored separately.\n")
  }
  )
  return(m)
}

 #__________________ Experiment object specific functions __________________________


#' @title Combine cell counts across all mice in an experiment into a single dataframe.
#' @description This function also stores the mouse attribute names (not experiment attributes) as columns that will be used as categorical variables to make analysis subgroups.
#' The values of these attributes (`group`, `drug`, `age`) will automatically be converted to a string values for consistency.
#' @param e experiment object
#' @param by (str) names of the experiment attributes (categorical variables) that will be used to create analysis subgroups.
#' @return
#' @export
#' @examples e <- combine_norm_cell_counts(e, by = c('groups', 'sex'))

combine_norm_cell_counts <- function(e, by){

  # Get experiment info
  e_info <- attr(e, "info")

  # Fix close but wrong attribute names
  if (!all(by %in% names(e_info))){
        by <- c(by[by %in% names(e_info)], m2e_attr(by[!by %in% names(e_info)]))
  }

  # initialize list to store the combined dataframes and the attributes
  combined_norm_counts_list <- vector(mode = "list", length(e_info$channels))
  names(combined_norm_counts_list) <- e_info$channels

  for (channel in e_info$channels){
    for (m in 1:length(e$mice)){

      # Get mouse info
      m_info <- attr(e$mice[[m]], "info")
      df <- e$mice[[m]]$normalized_counts[[channel]]
      for (attrib in by){
        # Add column keeping track of the mouse attribute
        add_col <- m_info[[e2m_attr(attrib)]] %>% toString() %>% tibble::tibble()
        names(add_col)  <- e2m_attr(attrib)
        df <- df %>% tibble::add_column(add_col, .before = TRUE)
      }

      # Always add the mouse ID
      df <- df %>% tibble::add_column(mouse_ID = m_info[["mouse_ID"]], .before = TRUE)

      if (m == 1){
        combined_norm_counts <- df
      } else{
        # add mouse dataframe of norm counts to growing combined norm counts table
        combined_norm_counts <- combined_norm_counts %>% dplyr::bind_rows(df)
      }
    }
    combined_norm_counts_list[[channel]] <- combined_norm_counts
  }
  e$combined_normalized_counts <- combined_norm_counts_list
  return(e)
}



#' @title Normalize colabel counts over a designated denominator channel.
#' @description This function can only be run after running [SMARTR::combine_norm_cell_counts()]. It divides the colabelled cell counts by
#' a designated normalization channel to provide a normalized ratio. Please note that the areas and volumes cancel out in this operation.
#' This is not designed to work on multiple hemispheres. Please combine cell counts across multiple hemispheres when you run [SMARTR::normalize_cell_counts()].
#'
#' @param e experiment object
#' @param denominator_channel (str, default = "eyfp") The exact name of the channel used for normalization
#'
#' @return e An experiment object with a new dataframe with the normalized ratios of colabelled counts over the designated denominator counts.
#' Because the volumes and region areas cancel out, the values of count, normalized.count.by.area, and normalized.count.by.volume are all the same.
#' This is to provide a consistent input dataframe into the analysis functions.
#' @export
#' @seealso [SMARTR::combine_norm_cell_counts()] &  [SMARTR::normalize_cell_counts()]
#' @examples e <- normalize_colabel_counts(e, denominator_channel = "eyfp")
normalize_colabel_counts <- function(e, denominator_channel = "eyfp"){


  colabel_normalized <- dplyr::inner_join(e$combined_normalized_counts$colabel,
                                          e$combined_normalized_counts[[denominator_channel]],
                                          by= c("mouse_ID",
                                                "group",
                                                "acronym",
                                                "name"),
                                          unmatched = "drop",
                                          suffix = c(".colabel", ".denom")) %>% dplyr::mutate(count = count.colabel/count.denom,
                                                                                       area.mm2 = area.mm2.denom,
                                                                                       volume.mm3 = volume.mm3.denom,
                                                                                       normalized.count.by.area = count.colabel/count.denom,
                                                                                       normalized.count.by.volume = count.colabel/count.denom)
  colabel_normalized <- colabel_normalized %>% dplyr::select(mouse_ID, group, acronym, name, count, area.mm2, volume.mm3, normalized.count.by.area, normalized.count.by.volume)
  # Returns normalized channel
  norm_chan <- paste0("colabel_over_", denominator_channel)
  e$combined_normalized_counts[[norm_chan]] <- colabel_normalized


  # Update the channel information
  attr(e, "info")$channels <- c( attr(e, "info")$channels, norm_chan)

  return(e)
}




#__________________ Internal Functions __________________________


# Rearrange things to be anatomical order

# anatomical_order<- c("Isocortex","OLF","HPF","CTXsp","CNU","TH","HY","MB","HB","CB")
# anat.order$region <- anatomical.order %>% map(get.sub.structure) %>% map(intersect,y=common.regions) %>% unlist()
# for(super.region in anatomical.order){
# anat.order$super.region[anat.order$region %in% get.sub.structure(super.region)] <- super.region
#}
#' Get aggregate volumes in the hippocampus
#'
#' @param m
#' @param AP_coord
#'
#' @return a list the length of channels, total_volumes
#'
#' @examples
get_hipp_DV_volumes <- function(m, AP_coord = -2.7, rois = c("DG", "CA1", "CA2", "CA3"), hipp_AP_coordinates){

  DV <- c("dorsal", "ventral")

  # Loop through all rois and get all subregions
  regions <- c(rois)
  for (roi in rois) {
    regions <- c(regions, wholebrain::get.sub.structure(roi))
  }

  # Loop through all slices and concatenate ONLY areas from the hippocampus.
  # Split these areas into dorsal and ventral
  aggregate_volumes_dorsal <- data.frame()
  aggregate_volumes_ventral <- data.frame()
  for (slice in m$slices){
    coordinate <- attr(slice, "info")$coordinate
    if (coordinate %in% hipp_AP_coordinates){
      # dorsal
      if (coordinate > AP_coord){
        aggregate_volumes_dorsal <- rbind(aggregate_volumes_dorsal, slice$volumes)
      }
      # Ventral
      if (coordinate <= AP_coord){
        aggregate_volumes_ventral <- rbind(aggregate_volumes_ventral, slice$volumes)
      }
    }
  }

  # Drop any NA values
  aggregate_volumes_dorsal  <- tidyr::drop_na(aggregate_volumes_dorsal)
  aggregate_volumes_ventral <- tidyr::drop_na(aggregate_volumes_ventral)

  # Store total volumes
  total_volumes_hipp <- list("dorsal" = aggregate_volumes_dorsal,
                              "ventral" = aggregate_volumes_ventral)

  for (dv in names(total_volumes_hipp)){
    if (length(total_volumes_hipp[[dv]]) > 0){

      total_volumes_hipp[[dv]] <- total_volumes_hipp[[dv]] %>% dplyr::group_by(acronym, right.hemisphere, name) %>%
        dplyr::summarise(area.mm2 = sum(area)*1e-6, volume.mm3 = sum(volume)*1e-9)

      # Store only regions in the hippocampus
      total_volumes_hipp[[dv]] <- total_volumes_hipp[[dv]][total_volumes_hipp[[dv]]$acronym %in% regions,]

    } else {
      message("There was no volume data found for the ", dv, " hippocampus!")
    }
  }
  return(total_volumes_hipp)
}


## Modification of Marcos' function
#' Get the registered areas
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

  regions.left <- tibble::tibble(acronym = regions.left, right.hemisphere = rep(FALSE, length(regions.left)))
  regions.right <-tibble::tibble(acronym = regions.right, right.hemisphere = rep(TRUE, length(regions.right)))
  regions <- rbind(regions.left, regions.right)

  # Sort region acronyms to alphabetical order
  regions <- regions[order(regions$acronym),]

  # Get region registration information
  region.info <- list()

  for (k in 1:nrow(regions)) {

    region.data <- wholebrain::get.region(regions$acronym[k],registration)
    region.data[,1:4] <- region.data[,1:4]*conversion.factor
    region.info <- c(region.info, list(region.data[region.data$right.hemisphere==regions$right.hemisphere[k],]))
  }

  areas <- tibble::tibble(name=as.character(wholebrain::name.from.acronym(regions$acronym)),
                      acronym=regions$acronym,
                      right.hemisphere=regions$right.hemisphere,
                      area=rep(0,nrow(regions)))

  # Use Gauss' area formula (shoelace formula)
  for (k in 1:length(region.info)) {
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
#' @param total_volumes computed total volumes list
#'
#' @return
#'
#' @examples
normalize.registered.areas <- function(cell_table, total_volumes){

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
    norm_counts <- merge(counts, total_volumes, by=c("acronym","right.hemisphere"), all = TRUE)
    norm_counts$normalized.count.by.area <- norm_counts$count/norm_counts$area.mm2
    norm_counts$normalized.count.by.volume <- norm_counts$count/norm_counts$volume.mm3

    # Add acronym names
    norm_counts$name <- as.character(wholebrain::name.from.acronym(norm_counts$acronym))

    # store the data into the normalized list
    normalized.list[[channel]] <- norm_counts %>% tidyr::drop_na() %>% tibble::as_tibble()
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
        simplified_counts_chan$acronym[k] <- wholebrain::get.acronym.parent(simplified_counts_chan$acronym[k])
        simplified_counts_chan$name[k] <- as.character(wholebrain::name.from.acronym(simplified_counts_chan$acronym[k]))

        # Continue storing storing the parent names until there are no more that match the keyword
        k <- grep(s, simplified_counts_chan$name, value = FALSE, ignore.case = TRUE)
      }
    }

    # Check if the data has been collapsed by hemisphere
    if (!"right.hemisphere" %in% names(simplified_counts_chan)){

      # step to collapse by name/acronym
      # simplified_counts_chan <- plyr::ddply(simplified_counts_chan, c("acronym", "name"), plyr::numcolwise(sum))

      simplified_counts_chan <- dplyr::group_by(simplified_counts_chan, acronym, name) %>%
        dplyr::summarise_at(c("count", "area.mm2", "volume.mm3"), sum) %>%
        dplyr::mutate(normalized.count.by.area = count/area.mm2,
                      normalized.count.by.volume = count/volume.mm3)

    } else {
      # step to collapse by name/acronym and hemisphere
      # simplified_counts_chan <- plyr::ddply(simplified_counts_chan, c("acronym", "name", "right.hemisphere"), plyr::numcolwise(sum))

      simplified_counts_chan <- dplyr::group_by(simplified_counts_chan, acronym, name, right.hemisphere) %>%
        dplyr::summarise_at(c("count", "area.mm2", "volume.mm3"), sum) %>%
        dplyr::mutate(normalized.count.by.area = count/area.mm2,
                      normalized.count.by.volume = count/volume.mm3)
    }

    # # Recalculated normalized count by area
    # simplified_counts_chan$normalized.count.by.area <- simplified_counts_chan$count/simplified_counts_chan$area.mm2
    #
    # # Recalculate normalized count by volume
    # simplified_counts_chan$normalized.count.by.volume <- simplified_counts_chan$count/simplified_counts_chan$volume.mm3

    # Store in simplified counts list
    simplified_counts[[channel]] <- simplified_counts_chan
  }

  return(simplified_counts)
}



#' Make a filter for the cfos or eyfp channel
#'
#' @param data
#' @param params
#' @param ranges
#'
#' @return
#'
#' @examples
make.filter <- function(data,
                        params = c("Vol..unit.","Moment1","Moment2","Moment3","Moment4","Sigma"),
                        ranges = list(c(200,12000),c(3,50),c(0,600),c(0,2000),c(0,5),c(20,Inf))){
  non.cells <- c()
  for (k in 1:length(params)){
    non.cells <- c(non.cells,which(data[,params[k]] <= ranges[[k]][1] | data[,params[k]] >= ranges[[k]][2]))
  }
  return(unique(non.cells))
}




#' Get colabelled cells data table. This is designed specifically to create a segmentation object from the imported raw files that are outputs from the batch_3D_MultiColocalization.ijm macro.
#' @param coloc_table
#' @param image_A_objects
#' @param image_B_objects
#' @param overlap (default = 0.5) Minimum fraction of object volume overlap from image A (Ch2) with object from image B (Ch1). Fraction is of image A objects.
#' @param volume (default = 25) Minimum threshold for colocalized volume in voxels.
#' @param euc_centroid_dist (default = 30) Euclidean threshold in pixels between the centroid coordinates of two flagged overlapping objects. If this distance is exceeded, there will be an error message.
#' @export
#' @return returns segmentation object with x Y coordinates as an average of the centroid coordinates. The area will be replaced with the volume of the image A object.
#' Intensity is also just replaced with the volume of the image A object.
#'
#' @examples
get.colabeled.cells <- function(coloc_table,
                                image_A_objects,
                                image_B_objects,
                                volume = 25,
                                euc_centroid_dist = 30,
                                overlap = 0.5){



  # calculate number of columns
  end_col <- (length(names(coloc_table)) - 2)/3

  # extracting column for objects, volumes, and proportions based on names
  mp_names <- paste0("P", 1:end_col)
  mv_names <- paste0("V", 1:end_col)
  mo_names <- paste0("O", 1:end_col)

  # Get column indices of the column with the largest objects percentage overlap
  mi <- max.col(coloc_table[mp_names], "first")

  # Get selection indices
  select <- cbind(1:length(coloc_table$X),mi)

  # Extracting max proportion
  mp <- coloc_table[,mp_names][select]
  mv <- coloc_table[,mv_names][select]
  mo <- coloc_table[,mo_names][select]

  # Filter out objects that are smaller than the volume threshold and the less than
  # The proportion overlap threshold
  colabel_ind <- which(mv>=volume & mp>=overlap)
  mo <- mo[colabel_ind]
  ma <- coloc_table$X[colabel_ind]

  x_a <- image_A_objects$CX..pix.[ma]
  y_a <- image_A_objects$CY..pix.[ma]
  x_b <- image_B_objects$CX..pix.[mo]
  y_b <- image_B_objects$CY..pix.[mo]

  # Quality check
  if (max(sqrt((x_a-x_b)^2+(y_a-y_b)^2)) > euc_centroid_dist){
    stop("The euclidian distance between the centroid coordinates of the
            overlapping objects seems too high! Please double check the output of the raw colocalization files!")
  }

  x <-  rowMeans(cbind(x_a, x_b))
  y <-  rowMeans(cbind(y_a, y_b))

  # Create a segmentation object
  seg.coloc <- SMARTR::segmentation.object
  seg.coloc$soma$x  <-   x
  seg.coloc$soma$y  <-   y
  seg.coloc$soma$area   <-   image_B_objects$Vol..pix.[mo]
  seg.coloc$soma$intensity  <-   image_B_objects$Vol..pix.[mo] # Dummy
  return(seg.coloc)
}







#' Find segmentation files following the naming conventions of the denny lab given a channel name and a root slice directory
#' @param slice_directory
#' @param channel
#' @return returns a vector of two paths. The first element is the path to the measurement data. The second path is the path to the
#' quantification data.
#' @examples
find_segmentation_files <- function(slice_directory, channel){

  files <-  list.files(slice_directory, pattern = "\\.txt")
  lc_files <- files %>% tolower()
  meas_path <- files[grepl(paste0("m_._", channel, "_."), lc_files)]
  quant_path <- files[grepl(paste0("q_._", channel, "_."), lc_files) & grepl(paste0(channel, ".txt"), lc_files)]

  if (length(meas_path) < 1 || length(quant_path) < 1){
    stop("Could not find the file within the slice directory. Please check your file naming conventions.")
  }
  return(c(meas_path, quant_path))
}


#' find_all_subregions
#'
#' @param regions A string vector of the Allen Mouse Brain Atlas abbreviated
#' regions of regions to exclude
#' @return all_regions. A string vector of all the input & subregions within them to exclude
#' @examples
find_all_subregions <- function(regions){

  # is.null check for exclude regions
  if (!is.null(regions)){
    # Containing all subregions to exclude too
    all_regions <- c()
    for (region in regions) {
      all_regions <- c(all_regions,
                       region,
                       SMARTR::get.sub.structure(region))
    }
    all_regions <- unique(all_regions)
  } else {
    all_regions <- NULL
  }
  return(all_regions)
}

















