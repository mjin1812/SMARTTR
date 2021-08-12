###_____________ General housekeeping functions __________####

#' @title Save mouse data
#' @description Saves mouse object into it's attribute output path as an RDATA file
#' @usage save_mouse(m)
#' @param ... parameter to pass mouse object
#' @param timestamp (bool) save the object with a date tag
#' @export
#' @example m <- save_mouse(m)

save_mouse <- function(..., timestamp = FALSE){
  info <- attr(..., 'info')

  if (timestamp){
    save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'_', Sys.Date(),'.RDATA')))
  } else {
    save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'.RDATA')))
  }
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






#__________________ Internal Functions __________________________

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
