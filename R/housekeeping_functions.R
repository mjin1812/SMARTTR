###_____________ General housekeeping functions __________####

#' Save mouse into it's attribute output path
#'
#' @usage save_mouse(m)
#' @param m mouse object
#' @export

save_mouse <- function(m){
  info <- attr(m, 'info')
  save(m, file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'.RDATA')))
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

