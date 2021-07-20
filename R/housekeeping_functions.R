###_____________ General housekeeping functions __________####

#' @title Save mouse data
#' @description Saves mouse object into it's attribute output path as an RDATA file
#' @usage save_mouse(m)
#' @param ... parameter to pass mouse object
#' @export
#' @example m <- save_mouse(m)

save_mouse <- function(...){
  info <- attr(..., 'info')
  save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'.RDATA'))  )
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

