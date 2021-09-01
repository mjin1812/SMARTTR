
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
#' @rdname import_segmentation
#' @description Method for importing segmentation data for a slice object
#' @param s slice object
#' @param mouse_ID (str) ID of mouse
#' @param channels (str vector) channels to import
#' @return s slice object
#' @examples s <-  import_segmentation(s, mouse_ID = "255", channels = c("cfos", "eyfp", "colabel"))
#' @export
#' @note Note that in order to import the raw data for co-labeled cells, you need both the txt files that end in "ColocOnly.txt" and "eYFP_16bit.txt".
#' The ImageJ colabelling algorithm calculates overlap between two segmented 3D objects from different channel.
#' One channel serves as the "reference" to later generate X & Y coordinates of colabelled cells in the make_segmentation_object() function.
#' By convention, we use the eYFP channel rather than cfos..
#'
# TODO: modify paths to search for correct files regardless of upper or lower case. Use grep and stringi to match multiple patterns

import_segmentation.slice <- function(s,
                                      mouse_ID = NA,
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
  txt_files <-  list.files(path = path_stem, all.files = TRUE, pattern = ".txt")

  for (k in 1:length(channels)){

    if (tolower(channels[k]) == 'colabel') {


      coloc_path <- paste0( mouse_ID,'_',section,"_Fast_R_cfos_SpotSegmentation_ColocOnly.txt")
      eyfp_16bit_path <- paste0("M_", mouse_ID, '_', section,"_Fast_G_eYFP_16bit.txt")

      # match to exact name in directory
      coloc_path <- txt_files[stringdist::amatch(coloc_path, txt_files, maxDist=Inf)]
      eyfp_16bit_path <- txt_files[stringdist::amatch(eyfp_16bit_path, txt_files, maxDist=Inf)]

    # ### CLEAN THIS UP AFTER NAMING CONVENTIONS ARE STANDARDIZED
    #   if (!file.exists(eyfp_16bit_path)){
    #     eyfp_16bit_path <- file.path(path_stem, paste0("M_", mouse_ID, '_', section, "_Fast_G_eYFP_LabelImage_16bit"))
    #   }

      print(coloc_path)
      print(eyfp_16bit_path)

      # Read
      coloc.table <- read.delim(file.path(path_stem, coloc_path), stringsAsFactors = FALSE)
      eyfp.meas.16bit <- read.csv(file.path(path_stem, eyfp_16bit_path), stringsAsFactors = FALSE)

      # store the coloc table and the eyfp 16 bit measurements as a combined list
      s$raw_segmentation_data[[k]] <- list(coloc_table = coloc.table,
                                           eyfp_counts_16bit = eyfp.meas.16bit)

    }

    else{
      if (tolower(channels[k]) == 'cfos'){

        # Create import data paths for cfos
        meas_path <- paste0("M_R_cfos_", mouse_ID,'_',section,".txt" )
        quant_path <- paste0("Q_R_cfos_",mouse_ID,'_',section,"_",'cfos',".txt" )

        # match to exact name in directory
        meas_path <- txt_files[stringdist::amatch(meas_path, txt_files, maxDist=Inf)]
        quant_path <- txt_files[stringdist::amatch(quant_path, txt_files, maxDist=Inf)]

      } else if (tolower(channels[k]) == 'eyfp'){


        # Create import data paths for eyfp
        meas_path <- paste0("M_G_eyfp_", mouse_ID,'_',section,".txt" )
        quant_path <- paste0("Q_G_eyfp_",mouse_ID,'_',section,"_eYFP",".txt" )

        # match to exact name in directory
        meas_path <- txt_files[stringdist::amatch(meas_path, txt_files, maxDist=Inf)]
        quant_path <- txt_files[stringdist::amatch(quant_path, txt_files, maxDist=Inf)]
      }

      print(meas_path)
      print(quant_path)
      meas <- read.csv( file.path(path_stem, meas_path), stringsAsFactors = FALSE )
      quant <- read.csv(file.path(path_stem, quant_path), stringsAsFactors = FALSE)
      meas$X2_Pix <- meas$CX..pix./(info$info$bin) #create position column to account for binning
      meas$Y2_Pix <- meas$CY..pix./(info$info$bin) #same as above
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
#' @rdname import_segmentation
#' @description Method for importing segmentation data for a mouse object
#' @param  m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str)'left', 'right' or NULL
#' @param channels (str vector) channels to import
#' @param replace (bool, default = FALSE) replace existing raw segmentation data
#' @return m mouse object
#' @examples m <-  import_segmentation(m, slice_ID = "1_10", channels = c("cfos", "eyfp", "colabel"), replace = FALSE)
#' @export


import_segmentation.mouse <- function(m,
                                      slice_ID = NA,
                                      hemisphere = NULL,
                                      channels = c('cfos', 'eyfp', 'colabel'),
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
##____ Method for making an eyfp filter____

#' Make a segmentation filter for a slice object
#' @rdname make_segmentation_filter
#' @param s
#' @param channels (str vector) channels to process, "eyfp" and "cfos" are legal.
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
                                           channels = c('eyfp'),
                                           params = list(c("Vol..unit.","Moment1","Moment2", "Moment3","Moment4","Sigma")),
                                           ranges = list(list(c(200,12000),
                                                              c(3,50),
                                                              c(0,600),
                                                              c(0,2000),
                                                              c(0,5),
                                                              c(20,Inf)))){

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

#' Make a segmentation filter for a slice within a mouse object
#' @rdname make_segmentation_filter
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere 'left', 'right' or NULL (both)
#' @param channels (str vector) Channels to process , "eyfp" and "cfos" are legal.
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
#' @param s slice object
#' @param mouse_ID (str) ID of mouse
#' @param channels (str vec) Channels to process.
#' @param use_filter (bool, default = FALSE). Use a filter to create more curated segmentation object from the raw segmentation data.
#' @param ... additional volume and overlap parameters for get.colabeled.cells().
#'
#' @return s slice object
#' @examples s <- make_segmentation_object(s, mouse_ID = "255", channels = c("cfos", "eyfp"), use_filter = FALSE)
#' @export
#' @note If you are processing the colabel channel, you must always have the raw eyfp segmentation data imported.

make_segmentation_object.slice <- function(s,
                                           mouse_ID = NA,
                                           channels = c('cfos', 'eyfp', 'colabel'),
                                           use_filter = FALSE,
                                           ... ){

  # Get slice information
  info <- attributes(s)

  ## Use registration image path to find segmentation data path stem
  # path_stem <- dirname(info$info$registration_path)


  # Check if colabel channel is included (MAY NOT BE NECESSARY)
  if ('colabel' %in% channels){
    # Rearrange the channels vector so that colabel is always processed last
    channels <- c(channels[!channels %in% 'colabel'], channels[channels %in% 'colabel'])

  }

  # Create a segmentation list to store the channel segmentatoin objects
  segmentation_list  <- vector(mode = 'list', length = length(channels))
  names(segmentation_list) <- names(channels)


  for (channel in channels){

    if (channel == 'colabel'){

      # filter check
      if (use_filter) {
        if (is.null(s$segmentation_filter[['eyfp']])) {
          warning(paste0("No filter is available for eyfp, even though `use_filter` parameter is TRUE.\nThe eyfp channel is used to process the colabel channel,",
                         "\n so the colabel channel will be processed without usage of an eyfp filter."))
          eyfp_counts <- s$raw_segmentation_data[['eyfp']]
        } else {
          eyfp_counts <- s$raw_segmentation_data[['eyfp']][-s$segmentation_filters[['eyfp']], ]
        }
      } else {
        eyfp_counts <- s$raw_segmentation_data[['eyfp']]
      }

      # obtain df of colocalized cell counts
      coloc.counts <- get.colabeled.cells(s$raw_segmentation_data$colabel$coloc_table,
                                        eyfp_counts,
                                        s$raw_segmentation_data$colabel$eyfp_counts_16bit,
                                        ...)

      seg.coloc <- SMARTR::segmentation.object
      seg.coloc$soma$x <- coloc.counts$X2_Pix
      seg.coloc$soma$y <- coloc.counts$Y2_Pix
      seg.coloc$soma$area <- coloc.counts$Vol..pix.
      seg.coloc$soma$intensity <- coloc.counts$Mean


      # store into segmentation list
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
      seg$soma$x <- counts$X2_Pix
      seg$soma$y <- counts$Y2_Pix
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
#' @param hemisphere (str) 'left', 'right' or NULL (both)
#' @param channels (vec) Channels to process.
#' @param replace (bool, default = FALSE) replace existing raw segmentation data
#' @param use_filter (bool, default = FALSE). Use a filter to create more curated segmentation object from the raw segmentation data.
#' @param ... additional volume and overlap parameters for get.colabeled.cells().
#' @examples m <-  make_segmentation_object(m, slice_ID = '1_9', hemisphere = 'left', channels = c('eyfp', 'cfos'), use_filter = FALSE)
#' @export
#'
make_segmentation_object.mouse <- function(m,
                                           slice_ID = NA,
                                           hemisphere = NULL,
                                           channels = c('cfos', 'eyfp', 'colabel'),
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
#' @param channels (vec) Channels to process.
#' @param clean (bool, default = TRUE). Remove cells that don't map to any regions.
#' @param display (bool, default = TRUE). Display the results of the forward warp for the slice.
#' @param mouse_ID (str) mouse ID
#' @param ... additional parameters besides 'registration', 'segmentation', 'forward.warps', and 'device' to pass to the [wholebrain::inspect.registration()] function
#' @examples s <- map_cells_to_atlas(s, channels c('cfos' , 'eyfp', 'colabel'), clean = TRUE, display = TRUE, mouse_ID = "255")
#' @export


map_cells_to_atlas.slice <- function(s,
                                     channels = c('cfos', 'eyfp', 'colabel'),
                                     clean =  TRUE,
                                     display = TRUE,
                                     mouse_ID = NULL,
                                     ...){

  for (k in 1:length(channels)){
    cell_data <- wholebrain::inspect.registration(s$registration_obj,
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
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str) 'left', 'right' or NULL (both)
#' @param channels (vec) Channels to process.
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
#' @param s slice object
#' @param channels (vec) Channels to process.
#' @param clean (bool, default = TRUE ). Remove cells that don't map to any regions.
#' @param exclude_regions (str vector, default = NULL); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
#' @param exclude_hemisphere (bool, default = TRUE); excludes the contralateral hemisphere from one indicated in slice attribute
#' @param exclude_layer_1 (bool, default = TRUE) excludes all counts from layer 1
#' @param plot_filtered (bool, default = TRUE) pop up window to check the excluded anatomy.
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
      dataset <- s$forward_warped_data[[channels[k]]]
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
#' @param channels (vec) Channels to process.
#' @param clean (bool, default = TRUE ). Remove cells that don't map to any regions.
#' @param exclude_regions (str vector, default = NULL); acronyms of regions you want to exclude IN ADDITION to regions that will by default be excluded in the slice attribute 'regions_excluded'
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
                                  channels = c('cfos', 'eyfp', 'colabel'),
                                  clean = TRUE,
                                  exclude_regions = NULL,
                                  exclude_hemisphere = TRUE,
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



# ##____Method for getting regional areas per slice ____
#' Method for getting regional areas per slice
#' @rdname get_registered_areas
#' @description Calculate the registered area (in microns^2^) of all regions contained in a slice.
#' @param s slice object
#' @return s slice object with a stored dataframe with columns 'name' (full region name), 'acronym', 'area' (in microns^2^), 'right.hemisphere'
#' @examples s <- get_registered_areas(s)
#' @export
#' @md

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
#' @param m mouse object
#' @param slice_ID (str) ID of slice
#' @param hemisphere (str) 'left', 'right' or NULL (both)
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
#' @param channels (vec) Channels to process.
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
#' Normalize cell counts per mm^2^ or by mm^3^ (if multiplying by the stack size)
#' @description Run this function after all the slices that you want to process are finished being added
#' and you have run the function get_cell_table() for your mouse.
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

      ## Check that this slice has run get_registered_areas
      if (is.null(s$areas)){
        stop(paste0("Slice ", slice_name ," in your mouse dataset has an empty areas vector.\n",
                    "Run the function get_registered_areas() for this slice before you can normalized all your cell counts by total area or by volume."))
        }
    }

  # ____________________________________________________________________________
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

  # Get z_width from slice attributes of the first slice
  z_width <- attr(m$slices[[1]], "info")$z_width

  # Summarize by acronym and right.hemisphere
  total_volumes <- dplyr::group_by(aggregate_areas, acronym, right.hemisphere, name) %>%
    dplyr::summarise(area.mm2 = sum(area)*1e-6, volume.mm3 = area.mm2*z_width*1e-3)
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

#' Split the hippocampal dataset to dorsal and ventral regions. If [normalize_cell_counts()] has already been run,
#' the existing hippocampus data will be deleted and either replaced with new dorsal/ventral hippocampus data with merge (recommended) or
#' it will be stored in a separate dataframe if the user wants analyze the hippocampus separately for specialized analysis.
#' @param m mouse object
#' @param AP_coord (double) The AP coordinate which separates dorsal hippocampus from ventral hippocampus. Counts anterior to this are considered dorsal, counts posterior to
#' this are considered ventral.
#' @param combine_hemispheres (bool, default = TRUE) Combine normalized cell counts from both hemispheres
#' @param merge (bool, default = TRUE) Recommended setting. Whether to include the division of dHipp and vHipp into the normalized counts dataframe of the mouse or to store as a separate element for isolated analysis.
#' Merging will replace the hippocampal counts currently in the normalized counts table. Not merging will create a separate data element in the mouse. TODO: Currently not functional.
#' @param simplify_regions (bool, default = TRUE ) simplify the normalized region counts based on keywords in the internal function, `simplify_keywords`
#' @param simplify_keywords (str, default =  c("layer","part","stratum","division")). Keywords to search through region names and simplify to parent structure
#' @param rois (str) List of the region acronyms to include in the hippocampus
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
    regions <- c(regions, wholebrain::get.sub.structure(roi))
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
    rois_ind <- which(m$cell_table[[channel]]$acronym %in% regions)
    hipp_df <- m$cell_table[[channel]][rois_ind,]

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
    if (!is.null(m$normalized_counts)){
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
  return(m)
}

#__________________ Experiment object specific functions __________________________


#' Combine cell counts across the mice in an experiment
#'
#' @param e experiment object
#' @param by (str) names of the experiment attributes used to divide up the cell counts. Will be same attributes to group by during analysis.
#'
#' @return
#' @export
#'
#' @examples e <- combine_norm_cell_counts(e, by = c('groups', 'sex'))
combine_norm_cell_counts <- function(e, by = c("groups", "sex")){

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
        add_col <- tibble::tibble(m_info[[e2m_attr(attrib)]])
        names(add_col)  <- e2m_attr(attrib)
        df <- df %>% tibble::add_column(add_col, .before = TRUE)
      }

      # Always add the mouse ID
      add_col <- tibble(m_info[["mouse_ID"]])
      names(add_col)  <- "mouse_ID"
      df <- df %>% tibble::add_column(add_col, .before = TRUE)

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

  aggregate_areas_dorsal <- data.frame()
  aggregate_areas_ventral <- data.frame()
  for (slice in m$slices){
    coordinate <- attr(slice, "info")$coordinate
    if (coordinate %in% hipp_AP_coordinates){
      # dorsal
      if (coordinate > AP_coord){
        aggregate_areas_dorsal <- rbind(aggregate_areas_dorsal, slice$areas)
      }
      # Ventral
      if (coordinate <= AP_coord){
        aggregate_areas_ventral <- rbind(aggregate_areas_ventral, slice$areas)
      }
    }
  }

  # Drop any NA values
  aggregate_areas_dorsal  <- tidyr::drop_na(aggregate_areas_dorsal)
  aggregate_areas_ventral <- tidyr::drop_na(aggregate_areas_ventral)

  # Get z_width from slice attributes of the first slice
  z_width <- attr(m$slices[[1]], "info")$z_width

  # Store total volumes
  total_volumes_hipp <- vector(mode = "list", length = 2)
  names(total_volumes_hipp) <- DV

  # iterator
  k <- 1
  for (areas in list(aggregate_areas_dorsal, aggregate_areas_ventral)){
    if (length(areas) > 0){
      total_volumes_hipp[[k]] <- dplyr::group_by(areas, acronym, right.hemisphere, name) %>%
        dplyr::summarise(area.mm2 = sum(area)*1e-6, volume.mm3 = area.mm2*24*1e-3 )

      # Store only regions in the hippocampus
      total_volumes_hipp[[k]] <- total_volumes_hipp[[k]][total_volumes_hipp[[k]]$acronym %in% regions,]

    } else {
      message("There was no area data found for the ", DV[k], " hippocampus!")
      total_volumes_hipp[[DV[k]]] <- NULL
    }
    k <- k + 1
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
    if (is.null(simplified_counts_chan$right.hemisphere)){

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

  # calculate number of columns
  end_col <- (length(names(coloc.table)) - 2)/3

  # extracting column for objects, volumes, and proportions based on names
  mp_names <- paste0("P", 1:end_col)
  mv_names <- paste0("V", 1:end_col)
  mo_names <- paste0("O", 1:end_col)

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


