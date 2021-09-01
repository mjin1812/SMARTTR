###_____________ General housekeeping functions __________####

#' @title Save experiment data
#' @description Saves experiment object into it's attribute output path as an RDATA file
#' @usage save_experiment(e)
#' @param ... parameter to pass experiment object
#' @param timestamp (bool) save the object with a date tag
#' @export
#' @example e <- save_experiment(e, timestamp = TRUE)

save_experiment <- function(..., timestamp = FALSE){
  info <- attr(..., 'info')

  if (timestamp){
    save(..., file = file.path(info$output_path, paste0(info$experiment_name,'_experiment_', Sys.Date(),'.RDATA')))
  } else {
    save(..., file = file.path(info$output_path, paste0(info$experiment_name,'_experiment','.RDATA')))
  }
}


#' @title Save mouse data
#' @description Saves mouse object into it's attribute output path as an RDATA file
#' @usage save_mouse(m)
#' @param ... parameter to pass mouse object
#' @param timestamp (bool) save the object with a date tag
#' @export
#' @example m <- save_mouse(m, timestamp = TRUE)

save_mouse <- function(..., timestamp = FALSE){
  info <- attr(..., 'info')

  if (timestamp){
    save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'_', Sys.Date(),'.RDATA')))
  } else {
    save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'.RDATA')))
  }
}


#' Print attributes of experiment object
#' @param e experiment object
#' @export

print.wb_experiment <- function(e){
  print(attr(e, 'info'))
}

## Printing methods for mouse and slices
#' Print attributes of mouse object
#' @param m mouse object
#' @export

print.mouse <- function(m){
  print(attr(m, 'info'))
}

#' Print attributes of slice object
#' @param s slice object
#' @export

print.slice <- function(s){
  print(attr(s, 'info'))
}


#' Print attributes of correlation_list object
#' @param s slice object
#' @export

print.correlation_list <- function(cl){
  print(attr(cl, 'info'))
}



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


#' Add mouse object to an experiment
#'
#' @description This function takes an experiment object and mouse object and adds the mouse object to the experiment. The mouse's unprocessed data, including all of it's individual
#' slice information and raw imported segmentation and registration data will be not be added from the mouse object to save space. Any desire to modify this data must be done at the mouse
#' object level before continuing further.
#'
#' This function will also read the individual mouse attributes and automatically populate the experimental attributes that are relevant. For example, the 'group' attribute of
#' a mouse will be read and automatically added to the experiment object's 'experimental_groups' attribute if it is a new unique experimental group name.
#'
#' @param e experiment object
#' @param m mouse object
#' @param replace (bool, default = FALSE) Replace a mouse already contained in an experiment object.
#' @return an experiment object
#' @export
#'
#' @examples
add_mouse <- function(e, m, replace = FALSE){

  # Read the mouse's attributes
  m_info <- attr(m, "info")
  mouse_ID <- m_info$mouse_ID


  # Read the experiment attributes
  e_info <- attr(e, 'info')

  # First mouse stored
  if (length(e$mice) < 1){
    m$slices <- NULL  # Delete slice information
    e$mice[[1]] <- m
    names(e$mice) <-  mouse_ID
  } else{
    # Check list of previously stored slice names
    stored_mice <- names(e$mice)

    # match flag
    match <- FALSE

    for (stored_mouse in stored_mice){
      if (identical(stored_mouse, mouse_ID)){

        match <- stored_mouse
        message(paste0('There was existing data found for mouse ', mouse_ID, '\n'))

        if (replace){
          # replace slice
          m$slices <- NULL # Delete slice information
          e$mice[[mouse_ID]] <- m
          message(paste0('Replaced existing mouse data!'))
        } else {
          stop(paste0('If you want to replace a previous mouse object,',
                      'then set the "replace" argument to "TRUE".'))
        }
      }
    }

    # If there were no matches, store the mouse as a new mouse
    if (isFALSE(match)){
      index <- length(e$mice) + 1
      m$slices <- NULL # Delete slice information
      e$mice[[index]] <- m
      names(e$mice)[index] <- mouse_ID
    }
  }

  # Check with the other experimental attributes and see if there are any new unique ones -> store it
  attr2match <- SMARTR::attr2match

  mouse_attr <- attr(m, 'info')
  exp_attr <- attr(e,"info")

  for (attrib in attr2match$exp){
    exp_attr[[attrib]] <- c(mouse_attr[[e2m_attr(attrib)]], exp_attr[[attrib]]) %>% unique()
  }

  # Detect the number of channels in the mouse
  channels <- names(m$normalized_counts)

  if (!all(channels %in% exp_attr[["channels"]])){
    exp_attr[["channels"]] <- c(exp_attr[["channels"]], channels) %>% unique()
  }

  # Reassign the attributes
  attr(e,"info") <- exp_attr

  return(e)
}











#__________________ Internal Functions __________________________


m2e_attr <- function(attribs){
  attr2match <- SMARTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$exp[stringdist::amatch(attr, attr2match$mouse, maxDist=Inf)])
  }
  return(matched)
}

e2m_attr <- function(attribs){
  attr2match <- SMARTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$mouse[stringdist::amatch(attr, attr2match$exp, maxDist=Inf)])
  }
  return(matched)
}

match_e_attr <- function(attribs){
  attr2match <- SMARTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$exp[stringdist::amatch(attr, attr2match$exp, maxDist=Inf)])
  }
  return(matched)
}

match_m_attr <- function(attribs){
  attr2match <- SMARTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$mouse[stringdist::amatch(attr, attr2match$mouse, maxDist=Inf)])
  }
  return(matched)
}





# Copies of Wholebrain's ontology search functions applied to the SMARTR custom ontology for the hippocampus

#' @export
name.from.acronym <- function (x)
{
  unlist(lapply(x, function(y) {
    if (length(which(SMARTR::ontology$acronym == y)) != 0) {
      return(SMARTR::ontology$name[which(SMARTR::ontology$acronym == y)])
    }
    else {
      return(NA)
    }
  }))
}

#' @export
name.from.id <- function (x)
{
  unlist(lapply(x, function(y) {
    if (length(which(SMARTR::ontology$id == y)) != 0) {
      return(SMARTR::ontology$name[which(SMARTR::ontology$id == y)])
    }
    else {
      return(NA)
    }
  }))
}


#' @export
get.sub.structure <- function (x)
{
  tmp <- get.acronym.child(x)
  if (sum(is.na(tmp)/length(tmp)) != 0) {
    tmp <- x
    return(tmp)
  }
  tmp2 <- tmp
  for (i in tmp) {
    tmp2 <- append(tmp2, get.sub.structure(i))
  }
  return(tmp2)
}



#' @export
get.sup.structure <- function (x, matching.string = c("CTX", "CNU", "IB",
                                 "MB", "HB", "grey", "root", "VS",
                                 "fiber tracts"))
{
  if (x %in% matching.string) {
    return(x)
  }
  tmp <- get.acronym.parent(x)
  if ((tmp %in% matching.string)) {
    tmp <- x
  }
  tmp2 <- tmp
  while (!(tmp %in% matching.string)) {
    tmp2 <- tmp
    tmp <- get.acronym.parent(tmp)
  }
  if (tmp == "root" | tmp == "grey" | tmp == "CH") {
    return(tmp2)
  }
  else {
    return(tmp)
  }
}


#' @export
get.acronym.child <-  function (x)
  {
    ids <- unlist(lapply(x, function(y) {
      if (length(which(SMARTR::ontology$parent == SMARTR::ontology$id[which(SMARTR::ontology$acronym ==
                                                            y)])) != 0) {
        if (y == "root") {
          return("997")
        }
        else {
          return(SMARTR::ontology$id[which(SMARTR::ontology$parent == SMARTR::ontology$id[which(SMARTR::ontology$acronym ==
                                                                          y)])])
        }
      }
      else {
        return(NA)
      }
    }))
    return(acronym.from.id(ids))
  }

#' @export
get.acronym.parent <- function (x)
{
  ids <- unlist(lapply(x, function(y) {
    if (length(which(SMARTR::ontology$acronym == y)) != 0) {
      if (y == "root") {
        return("997")
      }
      else {
        return(SMARTR::ontology$parent[which(SMARTR::ontology$acronym ==
                                       y)])
      }
    }
    else {
      return(NA)
    }
  }))
  return(acronym.from.id(ids))
}

#' @export
acronym.from.id <- function (x)
{
  unlist(lapply(x, function(y) {
    if (length(which(SMARTR::ontology$id == y)) != 0) {
      return(as.character(SMARTR::ontology$acronym[which(SMARTR::ontology$id ==
                                                   y)]))
    }
    else {
      return(NA)
    }
  }))
}
