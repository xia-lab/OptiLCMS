#' Extract isotope patterns from an annotated empirical compounds JSON
#'
#' Parse an annotated JSON file (e.g. Annotated_empiricalCompounds.json) and
#' return a list of data.frames where each data.frame contains two columns:
#' the first column is `mz` and the second column is `peak_area`.
#'
#' The function uses a set of heuristics to locate isotope information inside
#' common JSON layouts produced by annotation pipelines. It will ignore entries
#' where no isotope information can be found.
#'
#' @param json_path Path to the Annotated_empiricalCompounds.json file.
#' @return A named list of data.frames (columns: `mz`, `peak_area`).
#' @export
ExtractIsotopePatternsFromAnnotatedJSON <- function(json_path) {
  if (!file.exists(json_path)) stop("json_path does not exist: ", json_path)

  # read JSON (prefer jsonlite, fallback to RJSONIO)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    obj <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  } else if (requireNamespace("RJSONIO", quietly = TRUE)) {
    obj <- RJSONIO::fromJSON(json_path)
  } else {
    stop("Please install 'jsonlite' or 'RJSONIO' to parse JSON files")
  }

  nms <- names(obj)
  is_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

  # try to extract mz/intensity pairs and optional isotope/modification from a candidate object
  extract_pairs <- function(x) {
    # x may be: list of lists each with mz/intensity fields, or two parallel vectors,
    # or an array of [mz,int]
    if (is.null(x)) return(NULL)
    # case: matrix-like (numeric vectors)
    if (is.matrix(x) && ncol(x) >= 2) {
      df <- data.frame(mz = as.numeric(x[,1]), peak_area = as.numeric(x[,2]), isotope = NA_character_, modification = NA_character_, stringsAsFactors = FALSE)
      return(df)
    }
    # case: list of numeric pairs or named entries
    if (is.list(x) && length(x) > 0) {
      # if it's a list of length-2 numeric vectors
      if (all(vapply(x, function(el) is.numeric(el) && length(el) >= 2, logical(1L)))) {
        mz <- vapply(x, function(el) as.numeric(el[[1]]), numeric(1L))
        pa <- vapply(x, function(el) as.numeric(el[[2]]), numeric(1L))
        return(data.frame(mz = mz, peak_area = pa, isotope = NA_character_, modification = NA_character_, stringsAsFactors = FALSE))
      }

      # if it's a list of objects with fields like mz & intensity/area
      mz_vec <- numeric(0); pa_vec <- numeric(0); iso_vec <- character(0); mod_vec <- character(0)
      for (el in x) {
        if (is.null(el)) next
        if (is.list(el) && length(el) > 0) {
          # try common field names
          mz_val <- NULL
          pa_val <- NULL
          iso_val <- NA_character_
          mod_val <- NA_character_
          nm <- names(el)
          if (!is.null(nm)) {
            if ("mz" %in% nm) mz_val <- el[[which(nm=="mz")[1]]]
            if ("mass" %in% nm) mz_val <- mz_val %||% el[[which(nm=="mass")[1]]]
            # intensity synonyms
            for (cand in c("peak_area","area","intensity","abundance","height","int")) {
              if (cand %in% nm && is.numeric(el[[which(nm==cand)[1]]])) { pa_val <- el[[which(nm==cand)[1]]]; break }
            }
            # isotope label candidates
            for (ic in c("isotope","label","iso","isotope_label","isotopePattern","isotope_id")) {
              if (ic %in% nm) { iso_val <- as.character(el[[which(nm==ic)[1]]]); break }
            }
            # modification/adduct candidates
            for (mc in c("modification","adduct","adducts","ion","annotation","mod")) {
              if (mc %in% nm) { mod_val <- as.character(el[[which(nm==mc)[1]]]); break }
            }
          }
          # fallback: if el is numeric vector of length>=2
          if (is.null(mz_val) && is.numeric(el) && length(el) >= 2) mz_val <- el[[1]]
          if (is.null(pa_val) && is.numeric(el) && length(el) >= 2) pa_val <- el[[2]]
          if (!is.null(mz_val) && !is.null(pa_val) && is.finite(as.numeric(mz_val)) && is.finite(as.numeric(pa_val))) {
            mz_vec <- c(mz_vec, as.numeric(mz_val))
            pa_vec <- c(pa_vec, as.numeric(pa_val))
            iso_vec <- c(iso_vec, ifelse(is.na(iso_val), NA_character_, as.character(iso_val)))
            mod_vec <- c(mod_vec, ifelse(is.na(mod_val), NA_character_, as.character(mod_val)))
          }
        }
      }
      if (length(mz_vec) > 0) return(data.frame(mz = mz_vec, peak_area = pa_vec, isotope = iso_vec, modification = mod_vec, stringsAsFactors = FALSE))
    }

    # case: two parallel numeric vectors in a named list
    if (is.list(x) && !is.null(names(x))) {
      nm <- names(x)
      # find mz-like and area-like fields
      mz_field <- nm[grepl("mz|mass", nm, ignore.case = TRUE)][1]
      pa_field <- nm[grepl("peak|area|intensity|abundance|height", nm, ignore.case = TRUE)][1]
      iso_field <- nm[grepl("isotope|label|iso", nm, ignore.case = TRUE)][1]
      mod_field <- nm[grepl("mod|adduct|ion|annotation|modification", nm, ignore.case = TRUE)][1]
      if (!is.null(mz_field) && !is.null(pa_field)) {
        mzv <- x[[mz_field]]; pav <- x[[pa_field]]
        iso_v <- if (!is.null(iso_field)) as.character(x[[iso_field]]) else rep(NA_character_, length(mzv))
        mod_v <- if (!is.null(mod_field)) as.character(x[[mod_field]]) else rep(NA_character_, length(mzv))
        if (is.numeric(mzv) && is.numeric(pav) && length(mzv) == length(pav)) {
          return(data.frame(mz = as.numeric(mzv), peak_area = as.numeric(pav), isotope = as.character(iso_v), modification = as.character(mod_v), stringsAsFactors = FALSE))
        }
      }
    }

    # nothing matched
    NULL
  }

  # null-coalescing helper
  `%||%` <- function(a,b) if (!is.null(a)) a else b

  results <- list()
  add_result <- function(df, name = NULL) {
    if (is.null(df) || nrow(df) == 0) return()
    # normalize to four columns: mz, peak_area, isotope, modification
    if (is.data.frame(df)) {
      if (all(c("mz","peak_area") %in% names(df))) {
        mzv <- as.numeric(df[["mz"]])
        pav <- as.numeric(df[["peak_area"]])
        iso_v <- if ("isotope" %in% names(df)) as.character(df[["isotope"]]) else rep(NA_character_, length(mzv))
        mod_v <- if ("modification" %in% names(df)) as.character(df[["modification"]]) else rep(NA_character_, length(mzv))
      } else {
        mzv <- as.numeric(df[[1]])
        pav <- as.numeric(df[[2]])
        iso_v <- if (ncol(df) >= 3) as.character(df[[3]]) else rep(NA_character_, length(mzv))
        mod_v <- if (ncol(df) >= 4) as.character(df[[4]]) else rep(NA_character_, length(mzv))
      }
    } else {
      mzv <- as.numeric(df[[1]])
      pav <- as.numeric(df[[2]])
      iso_v <- rep(NA_character_, length(mzv))
      mod_v <- rep(NA_character_, length(mzv))
    }
    out_df <- data.frame(mz = mzv, peak_area = pav, isotope = iso_v, modification = mod_v, stringsAsFactors = FALSE)
    if (!is.null(name) && nzchar(name)) results[[name]] <<- out_df else results[[length(results) + 1]] <<- out_df
  }

  # locate candidate containers inside the JSON
  candidates <- list()
  if (is.list(obj) && length(obj) > 0) {
    # if top-level contains a named container of compounds
    top_names <- names(obj)
    if (!is.null(top_names)) {
      for (nm in c("empiricalCompounds","Annotated_empiricalCompounds","compounds","results","items")) {
        if (nm %in% top_names) {
          candidates <- obj[[nm]]
          break
        }
      }
    }
    # fallback: if not found, and top-level appears to be a list of compounds
    if (length(candidates) == 0) {
      # detect if obj is a list of objects with isotope-like fields
      maybe <- obj
      # if elements look like compound entries (each is a list)
      if (is.list(maybe) && length(maybe) > 0 && all(vapply(maybe, is.list, logical(1L)))) candidates <- maybe
    }
  }

  # iterate candidates and attempt extraction
  if (length(candidates) > 0) {
    for (idx in seq_along(candidates)) {
      item <- candidates[[idx]]
      name <- NULL
      if (is.list(item) && !is.null(names(item))) {
        if ("name" %in% names(item)) name <- as.character(item[["name"]])
        if ("id" %in% names(item) && is.null(name)) name <- as.character(item[["id"]])
      }

      # common locations for isotope info
      df <- NULL
      # 1) explicit 'isotopes' field
      if (is.list(item) && "isotopes" %in% names(item)) df <- extract_pairs(item[["isotopes"]])
      # 2) explicit 'isotopePattern' or 'isotope_pattern'
      if (is.null(df) && is.list(item) && ("isotopePattern" %in% names(item) || "isotope_pattern" %in% names(item))) {
        key <- ifelse("isotopePattern" %in% names(item), "isotopePattern", "isotope_pattern")
        df <- extract_pairs(item[[key]])
      }
      # 3) peaks field
      if (is.null(df) && is.list(item) && "peaks" %in% names(item)) df <- extract_pairs(item[["peaks"]])
      # 4) if item itself looks like mz/area vectors
      if (is.null(df)) df <- extract_pairs(item)

      # 5) try nested structures: some formats store 'empiricalCompound' with nested 'features'
      if (is.null(df) && is.list(item)) {
        # search sublists for isotope-like arrays
        for (sub in item) {
          if (is.list(sub)) {
            df <- extract_pairs(sub)
            if (!is.null(df)) break
          }
        }
      }

      add_result(df, name)
    }
  }

  names(results)<- nms
  return(results)
}
