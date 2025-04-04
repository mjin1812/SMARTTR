##______________ Package imports ________________
#' @importFrom ggplot2 ggplot aes theme element_text element_line element_rect element_blank unit labs ggsave guide_legend xlab ylab xlim ylim labs
#' @importFrom ggplot2 geom_col geom_tile geom_text geom_jitter geom_hline geom_vline geom_errorbar geom_point
#' @importFrom ggplot2 scale_y_continuous scale_x_continuous scale_fill_manual scale_fill_gradient2 scale_fill_viridis_c position_nudge
#' @importFrom tidyselect all_of any_of last_col
#' @importFrom dplyr n mutate summarize summarise across arrange group_by if_else desc
#' @importFrom tidygraph activate convert
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom grDevices dev.new
#' @importFrom utils write.csv read.csv setTxtProgressBar txtProgressBar
#' @importFrom magrittr %>%
#' @importFrom stats sd na.omit
#' @importFrom rlang .data

# ## ____________Globals bindings _____________

utils::globalVariables(c(".", "!!", ".orig_data", ".tidygraph_edge_index",
                         ".tidygraph_node_index","nodes", "edges"))

#' @title Save experiment data
#' @description Saves experiment object into it's attribute output path as an RDATA file
#'  save_experiment(e)
#' @param ... parameter to pass experiment object
#' @param timestamp (bool) save the object with a date tag
#' @export
#' @examples
#' \dontrun{
#' e <- save_experiment(e, timestamp = TRUE)
#' }
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
#'  save_mouse(m)
#' @param ... parameter to pass mouse object
#' @param timestamp (bool) save the object with a date tag
#' @export
#' @examples
#' \dontrun{
#' m <- save_mouse(m, timestamp = TRUE)
#' }
save_mouse <- function(..., timestamp = FALSE){
  info <- attr(..., 'info')

  if (timestamp){
    save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'_', Sys.Date(),'.RDATA')))
  } else {
    save(..., file = file.path(info$output_path, paste0('mouse_',info$mouse_ID,'.RDATA')))
  }
}

#' Print attributes of experiment object
#' @param x 	experiment object
#' @param ... further arguments passed to or from other methods.
#' @export
#' @examples
#' e <- experiment()
#' print(e)

print.experiment <- function(x, ...){
  print(attr(x, 'info'))
}

## Printing methods for mouse and slices
#' Print attributes of mouse object
#' @param x 	mouse object
#' @param ... further arguments passed to or from other methods.
#' @export
#' @examples
#' m <- mouse()
#' print(m)

print.mouse <- function(x, ...){
  print(attr(x, 'info'))
}

#' Print attributes of slice object
#' @param x 	slice object
#' @param ... further arguments passed to or from other methods.
#' @export
#' @examples
#' s <- slice()
#' print(s)

print.slice <- function(x, ...){
  print(attr(x, 'info'))
}

#' Print attributes of correlation_list object
#' @param x 	an object used to select a method. A correlation list
#' @param ... further arguments passed to or from other methods.
#' @export
print.correlation_list <- function(x, ...){
  print(attributes(x))
}

#' @title Add slice to a mouse object
#'
#' @param m mouse object
#' @param s slice object
#' @param replace (bool, default = FALSE) Replace a slice already contained in a mouse object.
#' @export
#' @examples
#' m <- mouse()
#' s <- slice()
#' m <- add_slice(m, s, replace = FALSE)
#' m <- add_slice(m, s, replace = TRUE)

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
#' @examples
#' e <- experiment(experiment_name = "example")
#' m <- mouse(mouse_ID = "Mouse1")
#' e <- add_mouse(e, m)
#' e <- add_mouse(e, m, replace = TRUE)
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
      if (identical(stored_mouse, toString(mouse_ID))){

        match <- stored_mouse
        message(paste0('There was existing data found for mouse ', mouse_ID, '\n'))

        if (replace){
          m$slices <- NULL
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
      m$slices <- NULL
      e$mice[[index]] <- m
      names(e$mice)[index] <- mouse_ID
    }
  }
  # Check with the other experimental attributes and see if there are any new unique ones -> store it
  attr2match <- SMARTTR::attr2match
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

#' Import externally mapped datasets into an experiment
#'
#' @description This function takes an experiment object and imports externally mapped datasets.
#'
#' @param e experiment object
#' @param normalized_count_paths (str vec, default = NULL, optional) For importing external datasets ONLY.
#' A character vector, with each element being a path to the respective .csv or .xlsx file in the same order as the `channels` attribute.
#' Must have the same order and length as the `channels` parameter.
#' You can add channels as names to the vector elements to avoid ambiguity and ensure correct importation for the right channel.
#' @param ... additional parameters to pass to either the [readr::read_csv()] function or [readxl::read_excel()]
#' @return e an experiment object with the imported dataset
#' @export
#'
import_mapped_datasets <- function(e, normalized_count_paths, ...){
  exp_attr <- attr(e,"info")
  exp_attr$channels
  if(is.null(names(normalized_count_paths))){
    names(normalized_count_paths)  <- exp_attr$channels
  } else {
    if (!all(names(normalized_count_paths) %in%  exp_attr$channels)){
      stop("The names for the channels you provided in `normalized_count_paths` does not match the `channels`
           attribute of your experiment object. Please check for case differences!")
    }
  }
  combined_normalized_counts <- vector(mode = "list", length = length(names(normalized_count_paths)))
  names(combined_normalized_counts) <- names(normalized_count_paths)
  for (channel in names(normalized_count_paths)){
    combined_normalized_counts[[channel]] <- read_check_file(normalized_count_paths[[channel]], ...) %>%
       dplyr::mutate(across(any_of("counts"):any_of("normalized.count.by.volume"), as.numeric)) %>%
      dplyr::mutate(across(!any_of("counts"):any_of("normalized.count.by.volume"), as.character))
  }

  columns <- colnames(combined_normalized_counts[[channel]])
  grouping_vars <- columns[!columns %in% c("mouse_ID", "acronym", "name", "counts", "volume.mm3", "normalized.count.by.volume")]
  if (length(grouping_vars) > 1){
    groups_df <- combined_normalized_counts[[channel]] %>% dplyr::select(any_of(grouping_vars)) %>% dplyr::group_by_at(grouping_vars[1]) %>% dplyr::summarise_at(grouping_vars[-1], unique)
  } else{
    groups_df <- combined_normalized_counts[[channel]] %>% dplyr::select(any_of(grouping_vars)) %>% dplyr::summarise_at(grouping_vars, unique)
  }
  for (variable in grouping_vars){
    variable_e <- m2e_attr(variable )
    exp_attr[[variable_e]] <-  groups_df[[variable]] %>% unique()
  }
  attr(e,"info") <- exp_attr
  e$combined_normalized_counts <- combined_normalized_counts
  return(e)
}


#' Reset the root path for the folder containing the registration and segmentation data.
#' @description This function takes a mouse object and also a `input_path` as the root folder for that mouse.
#' It then adjusts all the paths for the registration and segmentation data read to be relative to the root folder.
#' This function is especially useful if you have changed the computers you are analyzing and drive mappings may be different.
#' @param m mouse object
#' @param input_path (default = NULL) Reset the root directory of the mouse object.
#' @param print (bool, default = TRUE) Print the changes in the console.
#' @return m a mouse object
#' @export
#' @examples
#' \dontrun{
#' m <- reset_mouse_root(m, input_path = "C:/Users/Documents/Mice/mouse_1/", print = TRUE)
#' }
reset_mouse_root <- function(m, input_path = NULL, print = TRUE){

  if (is.null(input_path)){
    stop("Please provide a new root folder path using the input_path parameter...")
  }

  message("OLD mouse output path: ", attr(m, "info")$output_path, "\n",
          "New mouse output path: ", input_path)
  attr(m, "info")$output_path <- input_path

  root_files <- list.files(path = input_path, pattern = "MAX", recursive = TRUE)


  if (length(root_files) < 1){
    warning("There were no registration files found in the directory set as the input path. Please recheck where your folder is to make sure the input_path provided is correct.")
  } else {
    for (k in 1:length(m$slices)){
      s <- m$slices[[k]]
      s_reg_path <- attr(s, "info")$registration_path
      matched_root_file <- root_files[stringdist::amatch(s_reg_path, root_files, maxDist = Inf)]
      attr(m$slices[[k]], "info")$registration_path <- file.path(input_path, matched_root_file)
      if (print){
        message("Changed ", s_reg_path, " to ", attr(m$slices[[k]], "info")$registration_path)
      }
      if (!is.null(m$slices[[k]]$registration_obj)){

        root_files2 <- root_files %>%
          stringr::str_replace(., "_distorted.tif", "") %>%
          stringr::str_replace(., "_distorted.png", "") %>%
          stringr::str_replace(., "_undistorted.tif", "") %>%
          stringr::str_replace(., "_undistorted.png", "") %>% unique()

        matched_root_file <- root_files2[stringdist::amatch(m$slices[[k]]$registration_obj$outputfile, root_files2, maxDist = Inf)]
        m$slices[[k]]$registration_obj$outputfile <- file.path(input_path, matched_root_file)
      }

    }
  }
  return(m)
}

#__________________ Internal Functions __________________________

#' @noRd
m2e_attr <- function(attribs){
  attr2match <- SMARTTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$exp[stringdist::amatch(attr, attr2match$mouse, maxDist=Inf)])
  }
  return(matched)
}

#' @noRd
e2m_attr <- function(attribs){
  attr2match <- SMARTTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$mouse[stringdist::amatch(attr, attr2match$exp, maxDist=Inf)])
  }
  return(matched)
}

#' @noRd
match_e_attr <- function(attribs){
  attr2match <- SMARTTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$exp[stringdist::amatch(attr, attr2match$exp, maxDist=Inf)])
  }
  return(matched)
}

#' @noRd
match_m_attr <- function(attribs){
  attr2match <- SMARTTR::attr2match
  matched <- c()
  for (attr in attribs){
    matched <- c(matched, attr2match$mouse[stringdist::amatch(attr, attr2match$mouse, maxDist=Inf)])
  }
  return(matched)
}


#' Read a csv or excel file as a tibble. Checks first that the file exists, and that it is a csv or xlsx format.
#' @param x A file path
#' @param ... additional parameters to pass to either the [readr::read_csv] function or [readxl::read_excel]
#' @return tibble dataframe
#' @export
read_check_file <- function(x, ...) {
  # Check if file exists at path
  if (!file.exists(x)){
    stop("File does not exist at the expected location")
  }
  # Check file type
  ext <- tools::file_ext(x)
  file.read <- switch(ext,
         csv = ,
         .csv = readr::read_csv,
         xlsx = ,
         .xlsx = readxl::read_excel,
         stop("Invalid file type. `csv` or `xlsx` are allowed")
  )
  return(file.read(x, ...))
}


#' Check for redundant parent regions included in a list of acronyms in a plate. For example, if all the the subregions for
#' the hypothalamus are represented, the HY should not be included in the list.
#'
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param acronyms (vec) a vector of acronyms to check for possible parents that are redundantly included in the vector.
#' @return A list containing two elements: one vector of unique child acronyms, a vector of the parent regions considered redundant
#' @export
#' @examples
#' check_redundant_parents(acronyms = c("ACA5", "ACA1", "ACA"))
#' check_redundant_parents(acronyms = c("VMHDM", "VMHSh", "VMH"), ontology = "unified")
check_redundant_parents <- function(acronyms, ontology = "allen"){
  child_acronyms <- acronyms
  redundant_parents <- c()
  keep_going <- TRUE
  while (keep_going){

    if (tolower(ontology)=="allen"){
      parent_acronyms <- get.acronym.parent(child_acronyms)
    } else{
      parent_acronyms <- get.acronym.parent.custom(child_acronyms, ontology = ontology)
    }
    intersection <- intersect(parent_acronyms, acronyms)
    if (length(intersection) > 0){
      acronyms <- setdiff(acronyms, intersection)
      redundant_parents <- c(redundant_parents, intersection)
      child_acronyms <- parent_acronyms
    } else{
      keep_going <- FALSE
    }
  }
  redundant_parents <- redundant_parents %>% unique()

  # print(paste0("redundant parents ", redundant_parents)) # This is just for troubleshooting

  return(list(unique_acronyms = setdiff(acronyms, redundant_parents),
       redundant_parents = redundant_parents))
}


#' Simplify dataframe by keywords.
#'
#' @param df (tibble) Must contain columns "acronym" and "name"
#' @param keywords (vec, default = c("layer","part","stratum","division", "leaflet", "Subgeniculate", "island", "Islands", "Fields of Forel", "Cajal", "Darkschewitsch", "Precommissural")) a list of keywords to simplify based on region name.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param dont_fold (vec, default = c("Dorsal part of the lateral geniculate complex", "Ventral posterolateral nucleus of the thalamus, parvicellular part", "
#' Ventral posteromedial nucleus of the thalamus, parvicellular part","Ventral posterolateral nucleus of the thalamus, parvicellular part",
#' "Ventral posteromedial nucleus of the thalamus, parvicellular part","Substantia nigra")) Regions that are exceptions to being folded into their parent regions.
#' @return df
#' @export
#' @examples
#' df <- dplyr::tibble(acronym = c("LGd", "GU4", "dCA1so" ),
#' name = c("Dorsal part of the lateral geniculate complex",
#' "Gustatory areas, layer 4", "Field CA2, stratum oriens (dorsal)"))
#' simplify_by_keywords(df)
#'
#' df <- dplyr::tibble(acronym = c("LPBD", "MPBE", "CPre"),
#' name = c("Lateral parabrachial nucleus, dorsal part",
#' "Medial parabrachial nucleus, external part", "Caudoputamen- rostral extreme"))
#' simplify_by_keywords(df, keywords = c("dorsal part", "external part",
#' "Caudoputamen-"), ontology = "unified")
#'
simplify_by_keywords <- function(df,
                                 keywords = c("layer","part","stratum","division", "leaflet", "Subgeniculate", "island", "Islands", "Fields of Forel", "Cajal", "Darkschewitsch", "Precommissural"),
                                 ontology = "allen",
                                 dont_fold = c("Dorsal part of the lateral geniculate complex",
                                                "Ventral posterolateral nucleus of the thalamus, parvicellular part",
                                                "Ventral posteromedial nucleus of the thalamus, parvicellular part",
                                                "Ventral posterolateral nucleus of the thalamus, parvicellular part",
                                                "Ventral posteromedial nucleus of the thalamus, parvicellular part",
                                                "Substantia nigra")){

  for (s in keywords) {
    # Look for row indices where the names contain the keywords
    k <- grep(s, df$name, value = FALSE, ignore.case = TRUE)

    if (all(!is.na(dont_fold)) && all(!dont_fold=="")){
      k_omit <- grep(paste(dont_fold, collapse = "|"), df$name, value = FALSE, ignore.case = TRUE)
      k <- setdiff(k, k_omit)
    }

    while(length(k) > 0){
      if (tolower(ontology) == "allen"){
        # Store the parent acronym and full name
        df$acronym[k] <- get.acronym.parent(df$acronym[k])
        df$name[k] <- as.character(name.from.acronym(df$acronym[k]))
      } else {
        # df$acronym[k] <- SMARTTR::get.acronym.parent.custom(df$acronym[k], ontology = ontology)
        # df$name[k] <- as.character(SMARTTR::name.from.acronym.custom(df$acronym[k], ontology = ontology))
        ids <- id.from.acronym.custom(df$acronym[k], ontology = ontology)
        parent_ids <- parentid.from.id.custom(ids, ontology=ontology)
        df$acronym[k] <- acronym.from.id.custom(parent_ids, ontology = ontology)
        df$name[k] <- name.from.id.custom(parent_ids, ontology = ontology)
      }

      # Continue storing storing the parent names until there are no more that match the keyword
      k <- grep(s, df$name, value = FALSE, ignore.case = TRUE)
      if (all(!is.na(dont_fold)) && all((!dont_fold==""))){
        k_omit <- grep(paste(dont_fold, collapse = "|"), df$name, value = FALSE, ignore.case = TRUE)
        k <- setdiff(k, k_omit)
      }
    }
  }
  # df <- df %>% dplyr::arrange(acronym)
  return(df)
}


#' Simplify vector of acronyms by keywords.
#'
#' @param vec (vector) Must contain acronyms
#' @param keywords (vec, default = c("layer","part","stratum","division")) a list of keywords to simplify based on region name.
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @param dont_fold (vec, default = c("Dorsal part of the lateral geniculate complex", "Ventral posterolateral nucleus of the thalamus, parvicellular part", "
#' Ventral posteromedial nucleus of the thalamus, parvicellular part","Ventral posterolateral nucleus of the thalamus, parvicellular part",
#' "Ventral posteromedial nucleus of the thalamus, parvicellular part","Substantia nigra")) Regions that are exceptions to being folded into their parent regions.
#'
#' @return df, dataframe as a tibble with included long name and acronyms that are simplified to parents
#' @export
#'
#' @examples
#' acronyms <- c("LGd", "GU4", "dCA1so" )
#' simplify_vec_by_keywords(acronyms)
#'
#' acronyms <- c("LPBD", "MPBE", "CPre")
#' simplify_vec_by_keywords(acronyms,
#' keywords = c("dorsal part", "external part", "Caudoputamen-"), ontology = "unified")
simplify_vec_by_keywords <- function(vec,
                                     keywords = c("layer","part","stratum","division", "leaflet", "Subgeniculate", "island", "Islands", "Fields of Forel", "Cajal", "Darkschewitsch", "Precommissural"),
                                     ontology = "allen",
                                     dont_fold = c("Dorsal part of the lateral geniculate complex",
                                                   "Ventral posterolateral nucleus of the thalamus, parvicellular part",
                                                   "Ventral posteromedial nucleus of the thalamus, parvicellular part",
                                                   "Ventral posterolateral nucleus of the thalamus, parvicellular part",
                                                   "Ventral posteromedial nucleus of the thalamus, parvicellular part",
                                                   "Substantia nigra"
                                     )){
  name <- vector(mode="character", length = length(vec))
  for (v in 1:length(vec)){
    if (tolower(ontology) == "allen"){
      name[v] <- name.from.acronym(vec[v]) %>% as.character()
    } else {
      name[v] <- name.from.acronym.custom(vec[v], ontology = ontology) %>% as.character()
    }
  }

  df <- tibble::tibble(name=name,
                       acronym=vec) %>% tidyr::drop_na()

  # loop through each keyword
  for (s in keywords) {

    # Look for row indices where the names contain the keywords
    k <- grep(s, df$name, value = FALSE, ignore.case = TRUE)
    if (all(!is.na(dont_fold)) && all(!dont_fold=="")){
      k_omit <- grep(paste(dont_fold, collapse = "|"), df$name, value = FALSE, ignore.case = TRUE)
      k <- setdiff(k, k_omit)
    }

    while(length(k) > 0){

      if (tolower(ontology) == "allen"){
      # Store the parent acronym and full name
        df$acronym[k] <- get.acronym.parent(df$acronym[k])
        df$name[k] <- as.character(name.from.acronym(df$acronym[k]))
      } else {
        # df$acronym[k] <- SMARTTR::get.acronym.parent.custom(df$acronym[k], ontology = ontology)
        # df$name[k] <- as.character(SMARTTR::name.from.acronym.custom(df$acronym[k], ontology = ontology))
        ids <- id.from.acronym.custom(df$acronym[k], ontology = ontology)
        parent_ids <- parentid.from.id.custom(ids, ontology=ontology)
        df$acronym[k] <- acronym.from.id.custom(parent_ids, ontology = ontology)
        df$name[k] <- name.from.id.custom(parent_ids, ontology = ontology)

      }

      # Continue storing storing the parent names until there are no more that match the keyword
      k <- grep(s, df$name, value = FALSE, ignore.case = TRUE)
      if (all(!is.na(dont_fold)) && all(!dont_fold=="")){
        k_omit <- grep(paste(dont_fold, collapse = "|"), df$name, value = FALSE, ignore.case = TRUE)
        k <- setdiff(k, k_omit)
      }
    }
  }
  # df <- df %>% dplyr::arrange(acronym)
  return(df)
}

#' Get region ontology name from acronym
#' @description Get whole regional names from acronyms based on a lookup table including SMARTTR's custom ontology for the
#' dorsal ventral split of the hippocampus
#' @param x (str) Regional acronym or vector of regional acronyms
#' @return full ontology names of regions
#' @noRd

name.from.acronym <- function (x){
  x <- na.omit(x)
  unlist(lapply(x, function(y){ SMARTTR::ontology$name[SMARTTR::ontology$acronym %in% y] %>% return()})) %>% return()
}

#' Get region ontology name from ID
#' @description Similar to wholebrain package's search functions to get whole regional names from a numerical ID lookup table
#' including SMARTTR's custom ontology for the  dorsal ventral split of the hippocampus
#' @param x (int) integer ID
#' @return full ontology names of regions
#' @noRd

name.from.id <- function (x){
  x <- na.omit(x)
  unlist(lapply(x, function(y){ SMARTTR::ontology$name[SMARTTR::ontology$id %in% y] %>% return()})) %>% return()
}




#' Get subregion acronyms
#' @description Search function to get ALL substructure acronyms from parent acronym.
#' @param x (str vector) Regional acronyms in a vector
#' @return acroynms of subregions
#' @noRd

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




#' Updated function to return character vector of the super region acronyms after inputting a character vector of acronyms
#'
#' @param acronym (str vec) character vector of acronyms. Cannot be a factor vector.
#' @param anatomical.order (default = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU","TH", "HY", "MB", "HB", "CB")) Default way to group subregions into super regions order
#' @param ontology (str, default = "allen") Region ontology to use. options = "allen" or "unified"
#' @return vector the same length as acronym consisting of the super regions of the acronyms
#' @noRd
get.super.regions <- function(acronym, anatomical.order = c("Isocortex", "OLF", "HPF", "CTXsp", "CNU",
                                                            "TH", "HY", "MB", "HB", "CB"), ontology = "allen"){
  # get the parent super region
  super.region <- acronym %>% as.character()

  if (tolower(ontology) == "allen"){
    for (sup.region in anatomical.order){
      super.region[super.region %in% get.sub.structure(sup.region)] <- sup.region
    }
  } else {
    for (sup.region in anatomical.order){
      super.region[super.region %in% get.sub.structure.custom(sup.region, ontology = ontology)] <- sup.region
    }
  }
  return(super.region)
}

#' Get parent region acronyms
#' @description Function to get the acronym of parent regions in the ontology
#' @param x (str) Regional acronym
#' @param matching.string (str vector, default = c("CTX", "CNU", "IB", "MB", "HB", "grey", "root", "VS", "fiber tracts")) Vector of the basest parent levels to stop at.
#' @return acronyms as a str
#' @noRd

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


#' Get acronyms of child structures
#' @description Function to get the acronym of parent regions in the ontology
#' @param x (str) Regional acronym
#' @return acronyms as str
#' @noRd

get.acronym.child <-  function (x)
  {
    ids <- unlist(lapply(x, function(y) {
      if (length(which(SMARTTR::ontology$parent == SMARTTR::ontology$id[which(SMARTTR::ontology$acronym ==
                                                            y)])) != 0) {
        if (y == "root") {
          return("997")
        }
        else {
          return(SMARTTR::ontology$id[which(SMARTTR::ontology$parent == SMARTTR::ontology$id[which(SMARTTR::ontology$acronym ==
                                                                          y)])])
        }
      }
      else {
        return(NA)
      }
    }))
    return(acronym.from.id(ids))
  }

#' Get parent region acronyms
#' @description Function to get the acronym of parent regions
#' @param x (str) Regional acronym
#' @return acronyms as str
#' @noRd

get.acronym.parent <- function (x)
{
  ids <- unlist(lapply(x, function(y) {
    if (length(which(SMARTTR::ontology$acronym == y)) != 0) {
      if (y == "root") {
        return("997")
      }
      else {
        return(SMARTTR::ontology$parent[which(SMARTTR::ontology$acronym ==
                                       y)])
      }
    }
    else {
      return(NA)
    }
  }))
  return(acronym.from.id(ids))
}

#'  Retrieves acronym from an id
#' @param x region as an id number
#' @return acronyms
#' @noRd
acronym.from.id <- function (x)
{
  unlist(lapply(x, function(y) {
    if (length(which(SMARTTR::ontology$id == y)) != 0) {
      return(as.character(SMARTTR::ontology$acronym[which(SMARTTR::ontology$id ==
                                                   y)]))
    }
    else {
      return(NA)
    }
  }))
}

#' Get current OS
#' @return str of current operating system
#' @noRd

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
}


#' Code to evaluate an expression and
#'
#' @param expr expresion to evaluation
#' @param cpu n cpus to run with
#' @param elapsed time before timeout
#' @return either outcome of evaluated expression or timeout
#' @noRd
with_timeout <- function(expr, cpu=1, elapsed=1){
  expr <- substitute(expr)
  envir <- parent.frame()
  setTimeLimit(cpu = cpu, elapsed = elapsed, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  eval(expr, envir = envir)
}

# List of standard tiny regions to remove
# Oculomotor nucleus
# Trochlear nucleus
# Lateral terminal nucleus of the accessory optic tract
#
#
#

################################### Functions designed to accommodate custom atlases ############################
#' Get region ontology name from acronym
#' @description Get whole regional names from acronyms based on a custom lookup table imported to accomodate other ontologies.
#'
#' @param ontology (str, options: "unified")
#' @param x (str) Regional acronym or vector of regional acronyms
#' @return full ontology name
#' @noRd
name.from.acronym.custom <- function (x, ontology = "unified"){
  if (tolower(ontology) == "unified"){
    indices <- lapply(x, function(y){which(SMARTTR::ontology.unified$acronym %in% y)}) %>% unlist()
    return(SMARTTR::ontology.unified$name[indices])
    } else {
    stop("You did not enter a valid ontology name.")
  }
}


#' Get region ontology name from ID
#' @description Similar to wholebrain package's search functions to get whole regional names from a numerical ID lookup table. For custom ontologies.
#' @param ontology (str, options: "unified")
#' @param x (int) integer ID
#' @return full ontology name
#' @noRd

name.from.id.custom <- function (x, ontology = "unified"){
  if (tolower(ontology) == "unified"){
    unlist(lapply(x, function(y){SMARTTR::ontology.unified$name[SMARTTR::ontology.unified$id %in% y] %>% return()})) %>% return()
  }else {
    stop("You did not enter a valid ontology name.")
  }

}


#' Get region ontology  ID from acronym
#' @description Similar to wholebrain package's search functions to get whole regional names from a numerical ID lookup table. For custom ontologies.
#' @param ontology (str, options: "unified")
#' @param x (int) integer ID
#' @return return acronyms
#' @noRd
id.from.acronym.custom <- function (x, ontology = "unified"){
  # x <- na.omit(x)
  if (tolower(ontology) == "unified"){
    unlist(lapply(x, function(y){SMARTTR::ontology.unified$id[SMARTTR::ontology.unified$acronym %in% y] %>% return()})) %>% return()
  } else {
    stop("You did not enter a valid ontology name.")
  }
}


#' Get acronyms of child structures
#' @description Function to get the acronym of parent regions in the ontology
#' @param ontology (str, options: "unified")
#' @param x (str) Regional acronym
#' @return acronyms str
#' @noRd

get.acronym.child.custom <-  function (x, ontology = "unified")
{
  if (tolower(ontology) == "unified"){
    acronyms <- unlist(lapply(x, function(y) {
      if (length(which(SMARTTR::ontology.unified$parent_acronym == y) != 0)) {
        return(SMARTTR::ontology.unified$acronym[which(SMARTTR::ontology.unified$parent_acronym == y)])
      }
      else {
        return(NA)
      }
    }))
  } else {
    stop("You did not enter a valid ontology name.")
  }
  return(acronyms)
}

#' Get parent region acronyms
#' @description Function to get the acronym of parent regions
#' @param ontology (str, options: "unified")
#' @param x (str) Regional acronym
#' @return acronyms of parents
#' @noRd
get.acronym.parent.custom <- function (x, ontology = "unified")
{
  if (tolower(ontology) == "unified"){
    parents <- unlist(lapply(x, function(y) {
      if (length(which(SMARTTR::ontology.unified$acronym == y)) != 0) {
        if (y == "root") {
          return("")
        }
        else {
          return(SMARTTR::ontology.unified$parent_acronym[which(SMARTTR::ontology.unified$acronym == y)])
        }
      }
      else {
        return(NA)
      }
    }))
  } else {
    stop("You did not enter a valid ontology name.")
  }
  return(parents)
}

#' get acronym of region from its id
#'
#' @param x id of region
#' @param ontology which ontology to use
#' @return acronym of the region
#' @noRd

acronym.from.id.custom <- function (x, ontology = "unified"){
  if (tolower(ontology) == "unified"){
    acronyms <- unlist(lapply(x, function(y) {
      if (length(which(SMARTTR::ontology.unified$id == y)) != 0) {
        return(as.character(SMARTTR::ontology.unified$acronym[which(SMARTTR::ontology.unified$id == y)]))
      }
      else {
        return(NA)
      }
    }))
  } else {
    stop("You did not enter a valid ontology name.")
  }
  return(acronyms)
}



#' Get parent id from id
#' @param x id of region
#' @param ontology which ontology to use
#' @return id of parent region
#' @noRd

parentid.from.id.custom <- function (x, ontology = "unified"){
  if (tolower(ontology) == "unified"){
    acronyms <- unlist(lapply(x, function(y) {
      if (length(which(SMARTTR::ontology.unified$id == y)) != 0) {
        return(SMARTTR::ontology.unified$parent[which(SMARTTR::ontology.unified$id == y)])
      }
      else {
        return(NA)
      }
    }))
  } else {
    stop("You did not enter a valid ontology name.")
  }
  return(acronyms)
}

######### These functions are the recursive ones that should be used ###############
#' Get subregion acronyms
#' @description Search function to get ALL substructure acronyms from parent acronym. Note this functions is recursive.
#' @param ontology (str, options: "unified")
#' @param x (str vector) Regional acronyms in a vector
#' @return acronyms of substructures
#' @noRd
get.sub.structure.custom <- function(x, ontology = "unified"){
  if (tolower(ontology) == "unified"){
    tmp <- get.acronym.child.custom(x, ontology = ontology)
    if (sum(is.na(tmp)/length(tmp)) != 0) {
      tmp <- x
      return(tmp)
    }
    tmp2 <- tmp
    for (i in tmp) {
      # print(paste0("getting sub for ", i))
      # print(get.sub.structure.custom(i, ontology = ontology))
      tmp2 <- append(tmp2, get.sub.structure.custom(i, ontology = ontology))
    }
    tmp2 <- intersect(tmp2, SMARTTR::ontology.unified$acronym)
  } else {
    stop("You did not enter a valid ontology name.")
  }
  return(tmp2)
}

#' Get super parent region acronyms
#' @description Function to get the acronym of super parent regions in the ontology
#' @param x (str) Regional acronym
#' @param ontology (str, options: "unified")
#' @param matching.string (str vector, default = c("CTX", "CNU", "IB", "MB", "HB", "grey", "root", "VS", "fiber tracts")) Vector of the basest parent levels to stop at.
#' @return acronyms of parents
#' @noRd

get.sup.structure.custom <- function (x,
                                      ontology = "unified",
                                      matching.string = c("root", "grey", "CH", "VS", "CTX", "CNU", "MB", "HB"))
  {
  if (x %in% matching.string) {
    return(x)
  }
  tmp <- get.acronym.parent.custom(x, ontology = ontology)
  while (!(tmp %in% matching.string)) {
    tmp0 <- tmp
    tmp <- get.acronym.parent.custom(tmp0, ontology = ontology)
  }
  return(tmp)
}



#' Standard error function
#' @param x (vec)
#' @return numeric
#' @export
#' @examples sem(c(3,4,5))
sem <- function(x){
  sd(x) / sqrt(length(x))
}

