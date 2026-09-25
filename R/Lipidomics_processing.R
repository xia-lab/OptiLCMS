#' PerformLipidsAnalysis
#'
#' Perform lipidomics database searching after MS2 consensus generation.
#' This function queries a large sqlite lipidomics database by precursor m/z
#' (within ppm tolerance) using Rcpp, then scores candidate spectra with
#' \code{calculate_score}.
#'
#' Expected upstream steps:
#' \code{PerformMSnImport()} -> \code{PerformDDADeconvolution()}/\code{PerformDIADeconvolution()} -> \code{PerformSpectrumConsenus()}.
#'
#' @param mSet mSet object containing consensus spectra in \code{mSet@MSnResults$Concensus_spec}.
#' @param database_path Path to lipidomics sqlite database.
#' @param precursor_ppm Numeric ppm tolerance used to extract precursor-matched candidates.
#' @param max_candidates Integer, max candidates fetched per feature across all tables.
#' @param max_per_table Integer, max candidates fetched from each table.
#' @param top_n Integer, number of top-scoring hits kept for each feature.
#' @param parameter Optional search parameter list produced by \code{new_ms_ref_search_parameter()}.
#' @param source_type Character source tag passed to \code{calculate_score}.
#'
#' @return Updated \code{mSet} object with:
#' \itemize{
#'   \item \code{mSet@MSnResults$LipidMatchRes}: per-feature ranked scoring results.
#'   \item \code{mSet@MSnResults$LipidTopHits}: summary table of the top-\code{top_n} hits per feature.
#' }
#' @export
PerformLipidsAnalysis <- function(
  mSet = NULL,
  database_path = "/Volumes/ExtremeSSD/Lipidomics_projects/database_curation/lipidomics_complete.sqlite",
  precursor_ppm = 10,
  max_candidates = 300L,
  max_per_table = 100L,
  top_n = 5L,
  parameter = NULL,
  source_type = "MspDB",
  use_isotopes = FALSE) {

  if (is.null(mSet)) {
    stop("mSet is required. Please run PerformSpectrumConsenus() first.")
  }
  if (is.null(mSet@MSnResults[["Concensus_spec"]])) {
    stop("Consensus spectra missing. Please run PerformSpectrumConsenus() first.")
  }
  if (!file.exists(database_path)) {
    stop("database_path does not exist: ", database_path)
  }
  if (!is.numeric(precursor_ppm) || precursor_ppm <= 0) {
    stop("precursor_ppm must be a positive numeric value.")
  }
  if (!is.integer(max_candidates)) {
    max_candidates <- as.integer(max_candidates)
  }
  if (!is.integer(max_per_table)) {
    max_per_table <- as.integer(max_per_table)
  }
  if (!is.integer(top_n)) {
    top_n <- as.integer(top_n)
  }
  if (max_candidates <= 0L || max_per_table <= 0L || top_n <= 0L) {
    stop("max_candidates, max_per_table and top_n must all be positive integers.")
  }
  if (!exists("calculate_score", mode = "function", inherits = TRUE)) {
    stop("calculate_score() was not found. Please ensure LipidsScoring.R is loaded.")
  }
  if (!exists("lipidomics_db_search_precursor_cpp", mode = "function", inherits = TRUE)) {
    stop("lipidomics_db_search_precursor_cpp() was not found. Reinstall/rebuild package to register Rcpp exports.")
  }

  if (is.null(parameter)) {
    parameter <- new_ms_ref_search_parameter()
  }

  .assemble_candidates_from_list <- function(candidates_list) {
    if (length(candidates_list) == 0) return(data.frame())
    ncan <- length(candidates_list)
    Table <- character(ncan)
    ID <- integer(ncan)
    PrecursorMZ <- numeric(ncan)
    Adduct <- character(ncan)
    ExactMass <- numeric(ncan)
    MolecularWeight <- numeric(ncan)
    Name <- character(ncan)
    InChIKey <- character(ncan)
    CompoundClass <- character(ncan)
    Comment <- character(ncan)
    RetentionTime <- numeric(ncan)
    CollisionCrossSection <- numeric(ncan)
    MS2Peaks_list <- vector("list", ncan)
    AbbreName <- character(ncan)
    for (k in seq_len(ncan)) {
      cnd <- candidates_list[[k]]
      Table[k] <- if (!is.null(cnd$Table)) as.character(cnd$Table) else ""
      ID[k] <- if (!is.null(cnd$RecordIndex)) as.integer(cnd$RecordIndex) else NA_integer_
      PrecursorMZ[k] <- if (!is.null(cnd$PrecursorMZ)) as.numeric(cnd$PrecursorMZ) else NA_real_
      Adduct[k] <- if (!is.null(cnd$Adduct)) as.character(cnd$Adduct) else ""
      ExactMass[k] <- if (!is.null(cnd$ExactMass)) as.numeric(cnd$ExactMass) else NA_real_
      MolecularWeight[k] <- if (!is.null(cnd$MolecularWeight)) as.numeric(cnd$MolecularWeight) else NA_real_
      Name[k] <- if (!is.null(cnd$IUPACName) && nzchar(as.character(cnd$IUPACName))) as.character(cnd$IUPACName) else if (!is.null(cnd$AbbreName)) as.character(cnd$AbbreName) else ""
      AbbreName[k] <- if (!is.null(cnd$AbbreName)) as.character(cnd$AbbreName) else ""
      InChIKey[k] <- if (!is.null(cnd$InChIKey)) as.character(cnd$InChIKey) else ""
      CompoundClass[k] <- if (!is.null(cnd$SpecType)) as.character(cnd$SpecType) else ""
      Comment[k] <- ""
      RetentionTime[k] <- if (!is.null(cnd$RetentionTime)) as.numeric(cnd$RetentionTime) else NA_real_
      CollisionCrossSection[k] <- if (!is.null(cnd$CollisionCrossSection)) as.numeric(cnd$CollisionCrossSection) else NA_real_
      MS2Peaks_list[[k]] <- if (!is.null(cnd$Spectrum)) as.matrix(cnd$Spectrum) else matrix(numeric(0), ncol = 2)
    }
    df <- data.frame(
      Table = Table,
      ID = ID,
      Adduct = Adduct,
      ExactMass = ExactMass,
      MolecularWeight = MolecularWeight,
      AbbreName = AbbreName,
      PrecursorMZ = PrecursorMZ,
      Name = Name,
      InChIKey = InChIKey,
      CompoundClass = CompoundClass,
      Comment = Comment,
      RetentionTime = RetentionTime,
      CollisionCrossSection = CollisionCrossSection,
      stringsAsFactors = FALSE
    )
    df$MS2Peaks <- I(MS2Peaks_list)
    df
  }

  # if LipidDBMatchRes already computed, use it and skip DB searching
  precomputed_matches <- NULL
  if (!is.null(mSet@MSnResults[["LipidDBMatchRes"]])) {
    precomputed_matches <- mSet@MSnResults[["LipidDBMatchRes"]]
  }

  .matrix_to_peak_list <- function(spec_mtx) {
    if (is.null(spec_mtx) || length(spec_mtx) == 0) return(list())
    spec_mtx <- as.matrix(spec_mtx)
    if (nrow(spec_mtx) == 0 || ncol(spec_mtx) < 2) return(list())
    lapply(seq_len(nrow(spec_mtx)), function(i) {
      list(Mass = as.numeric(spec_mtx[i, 1]), Intensity = as.numeric(spec_mtx[i, 2]))
    })
  }

  .generate_realtime_oad_reference <- function(reference, candidate) {
    if (!identical(source_type, "GeneratedLipid") ||
        !identical(parameter$CollisionType, "OAD") ||
        !exists("generate_oad_lipid_spectrum", mode = "function", inherits = TRUE)) {
      return(reference)
    }

    lipid_name <- if (!is.null(candidate$AbbreName) &&
                      nzchar(as.character(candidate$AbbreName))) {
      as.character(candidate$AbbreName)
    } else {
      as.character(reference$Name)
    }
    adduct <- if (!is.null(candidate$Adduct) &&
                  nzchar(as.character(candidate$Adduct))) {
      as.character(candidate$Adduct)
    } else if (grepl("-$", ion_mode)) {
      "[M-H]-"
    } else {
      "[M+H]+"
    }

    neutral_mass <- suppressWarnings(as.numeric(candidate$ExactMass))
    if (!is.finite(neutral_mass)) {
      neutral_mass <- suppressWarnings(as.numeric(candidate$MolecularWeight))
    }
    if (!is.finite(neutral_mass)) {
      precursor <- suppressWarnings(as.numeric(reference$PrecursorMz))
      neutral_mass <- switch(adduct,
        "[M+H]+" = precursor - .oad_sg_proton,
        "[M+Na]+" = precursor - .oad_sg_na,
        "[M+NH4]+" = precursor - .oad_sg_nh4,
        "[M+H-H2O]+" = precursor - .oad_sg_proton + .oad_sg_water,
        "[M-H2O+H]+" = precursor - .oad_sg_proton + .oad_sg_water,
        "[M-H]-" = precursor + .oad_sg_proton,
        "[M+HCOO]-" = precursor - .oad_sg_formate + .oad_sg_proton,
        "[M+CH3COO]-" = precursor - .oad_sg_acetate + .oad_sg_proton,
        "[M+HCO3]-" = precursor - .oad_sg_bicarbonate + .oad_sg_proton,
        NA_real_
      )
    }
    if (!is.finite(neutral_mass)) return(reference)

    lipid <- list(
      Name = lipid_name,
      Mass = neutral_mass,
      InChIKey = reference$InChIKey
    )
    if (!isTRUE(can_generate_oad_lipid_spectrum(lipid, adduct))) return(reference)

    generated <- tryCatch(
      generate_oad_lipid_spectrum(lipid, adduct),
      error = function(e) NULL
    )
    if (is.null(generated) || is.null(generated$Spectrum) ||
        nrow(generated$Spectrum) == 0L) {
      return(reference)
    }

    reference$Spectrum <- .matrix_to_peak_list(generated$Spectrum[, c("Mass", "Intensity"), drop = FALSE])
    reference$PrecursorMz <- as.numeric(generated$PrecursorMz)
    reference$Name <- generated$Name
    reference$CompoundClass <- generated$CompoundClass
    reference$Adduct <- generated$AdductIonName
    reference
  }

  .parse_ms2peaks_local <- function(ms2peaks_str) {
    if (is.null(ms2peaks_str) || is.na(ms2peaks_str) || !nzchar(ms2peaks_str)) {
      return(matrix(numeric(0), ncol = 2))
    }
    if (exists("parse_ms2peaks", mode = "function", inherits = TRUE)) {
      return(parse_ms2peaks(ms2peaks_str))
    }

    lines <- strsplit(ms2peaks_str, "\n", fixed = TRUE)[[1]]
    lines <- lines[nzchar(lines)]
    if (length(lines) == 0) return(matrix(numeric(0), ncol = 2))

    vals <- lapply(lines, function(line) {
      sp <- strsplit(trimws(line), "[[:space:]]+", perl = TRUE)[[1]]
      nums <- suppressWarnings(as.numeric(sp))
      nums[!is.na(nums)]
    })
    vals <- vals[vapply(vals, function(x) length(x) >= 2, logical(1L))]
    if (length(vals) == 0) return(matrix(numeric(0), ncol = 2))

    mat <- do.call(rbind, lapply(vals, function(x) c(x[1], x[2])))
    matrix(as.numeric(mat), ncol = 2)
  }

  .get_feature_precursor_rt <- function(MSn_data, consensus_idx_zero_based) {
    row_idx <- consensus_idx_zero_based + 1L
    if (MSn_data$acquisitionMode == "DIA") {
      pk <- MSn_data$peak_mtx[[row_idx]]
      list(precursor = as.numeric(pk[1]), rt = as.numeric(pk[4]))
    } else {
      pk <- MSn_data$peak_mtx[row_idx, ]
      list(
        precursor = mean(as.numeric(pk[1:2]), na.rm = TRUE),
        rt = mean(as.numeric(pk[3:4]), na.rm = TRUE)
      )
    }
  }

  # isotope fetching is handled in C++: lipidomics_get_isotope_pattern_cpp(db_path, lipid_index)

  # Extract simple scan isotope envelope from consensus spec around precursor
  .extract_scan_isotopes <- function(spec_mtx, precursor_mz, ppm = 15, max_iso = 3) {
    if (is.null(spec_mtx) || length(spec_mtx) == 0 || !is.finite(precursor_mz)) return(matrix(numeric(0), ncol = 2))
    spec_mtx <- as.matrix(spec_mtx)
    if (nrow(spec_mtx) == 0 || ncol(spec_mtx) < 2) return(matrix(numeric(0), ncol = 2))
    iso_mass <- 1.003355
    out <- list()
    for (k in 0:max_iso) {
      target <- precursor_mz + k * iso_mass
      tol <- target * ppm / 1e6
      rows <- which(abs(spec_mtx[, 1] - target) <= tol)
      if (length(rows) == 0) next
      # sum intensities for this isotope bin
      mass <- mean(spec_mtx[rows, 1])
      inten <- sum(as.numeric(spec_mtx[rows, 2]), na.rm = TRUE)
      out[[length(out) + 1]] <- c(mass, inten)
    }
    if (length(out) == 0) return(matrix(numeric(0), ncol = 2))
    mat <- do.call(rbind, out)
    matrix(as.numeric(mat), ncol = 2)
  }

  consensus <- mSet@MSnResults[["Concensus_spec"]]
  feature_ids <- consensus[[1]]
  consensus_specs <- consensus[[2]]
  MSn_data <- mSet@MSnData
  ion_mode <- if (isTRUE(MSn_data$ion_mode == 1L)) "Positive" else "Negative"

  lipid_match_res <- vector("list", length(feature_ids))
  top_hit_rows <- vector("list", length(feature_ids))

  for (i in seq_along(feature_ids)) {
    
    if (i %% 100L == 0L) {
      message("PerformLipidsAnalysis processed ", i, " / ", length(feature_ids), " consensus spectra")
    }
    
    idx0 <- as.integer(feature_ids[[i]])
    spec_mtx <- consensus_specs[[i]][[1]]
    spec_list <- .matrix_to_peak_list(spec_mtx)

    if (length(spec_list) == 0) {
      lipid_match_res[[i]] <- list(
        FeatureIndex = idx0,
        CandidateCount = 0L,
        Scores = list(),
        CandidateMeta = data.frame()
      )
      top_hit_rows[[i]] <- NULL
      next
    }

    feature_info <- .get_feature_precursor_rt(MSn_data, idx0)
    precursor_mz <- feature_info$precursor
    feature_rt <- feature_info$rt

    if (!is.finite(precursor_mz) || precursor_mz <= 0) {
      lipid_match_res[[i]] <- list(
        FeatureIndex = idx0,
        CandidateCount = 0L,
        Scores = list(),
        CandidateMeta = data.frame()
      )
      top_hit_rows[[i]] <- NULL
      next
    }

    # If precomputed matches exist, use them; otherwise perform DB search
    if (!is.null(precomputed_matches)) {
      raw_cand <- precomputed_matches[[i]]
      if (is.data.frame(raw_cand)) {
        candidates <- raw_cand
        if (!"AbbreName" %in% names(candidates)) candidates$AbbreName <- ""
      } else if (is.list(raw_cand) && length(raw_cand) > 0) {
        candidates <- .assemble_candidates_from_list(raw_cand)
      } else {
        candidates <- data.frame()
      }
    } else {
      # Prefer blob-returning search (structured spectra) when available
      cpp_search_blob <- NULL
      if (exists("lipidomics_db_search_precursors_blob_cpp", mode = "function", inherits = TRUE)) {
        cpp_search_blob <- get("lipidomics_db_search_precursors_blob_cpp", mode = "function", inherits = TRUE)
      }

      if (!is.null(cpp_search_blob)) {
        spectra_table <- if (ion_mode == "Positive") "spectra_pos" else "spectra_neg"
        res_blob <- cpp_search_blob(
          database_path = database_path,
          spectra_table = spectra_table,
          precursor_mzs = as.numeric(c(precursor_mz)),
          ppm = precursor_ppm,
          max_candidates = max_candidates,
          max_per_table = max_per_table
        )
        # res_blob is a list per precursor; take first element
        candidates_list <- if (length(res_blob) >= 1) res_blob[[1]] else list()
        candidates <- .assemble_candidates_from_list(candidates_list)
      } else {
        cpp_search <- get("lipidomics_db_search_precursor_cpp", mode = "function", inherits = TRUE)
        candidates <- cpp_search(
          database_path = database_path,
          precursor_mz = precursor_mz,
          ppm = precursor_ppm,
          max_candidates = max_candidates,
          max_per_table = max_per_table
        )
        if (is.data.frame(candidates) && !"AbbreName" %in% names(candidates)) candidates$AbbreName <- ""
      }
    }

    if (nrow(candidates) == 0) {
      lipid_match_res[[i]] <- list(
        FeatureIndex = idx0,
        CandidateCount = 0L,
        Scores = list(),
        CandidateMeta = candidates
      )
      top_hit_rows[[i]] <- NULL
      next
    }

    scan <- list(
      Spectrum = spec_list,
      PrecursorMz = precursor_mz,
      IonMode = ion_mode
    )

    property <- list(
      PrecursorMz = precursor_mz,
      ChromXs_RT = feature_rt,
      CollisionCrossSection = NA_real_
    )

    # prepare isotope info if requested
    if (use_isotopes) {
      scan_isotopes_pre <- .extract_scan_isotopes(spec_mtx, precursor_mz, ppm = precursor_ppm, max_iso = 3)
      scan_isotopes_arg <- if (is.matrix(scan_isotopes_pre) && nrow(scan_isotopes_pre) > 0) list(scan_isotopes_pre) else list()
    } else {
      scan_isotopes_arg <- list()
    }

    scored <- vector("list", nrow(candidates))
    for (j in seq_len(nrow(candidates))) {
      ms2p_obj <- candidates$MS2Peaks[[j]]
      if (is.matrix(ms2p_obj)) {
        ref_mtx <- as.matrix(ms2p_obj)
      } else {
        ref_mtx <- .parse_ms2peaks_local(ms2p_obj)
      }
      ref_spec <- .matrix_to_peak_list(ref_mtx)

      ref_prec <- as.numeric(candidates$PrecursorMZ[[j]])
      if (!is.finite(ref_prec)) ref_prec <- precursor_mz

      reference <- list(
        Spectrum = ref_spec,
        PrecursorMz = ref_prec,
        Name = as.character(candidates$Name[[j]]),
        ScanID = as.integer(candidates$ID[[j]]),
        InChIKey = as.character(candidates$InChIKey[[j]]),
        CompoundClass = as.character(candidates$CompoundClass[[j]]),
        Comment = as.character(candidates$Comment[[j]]),
        ChromXs_RT = as.numeric(candidates$RetentionTime[[j]]),
        CollisionCrossSection = as.numeric(candidates$CollisionCrossSection[[j]])
      )

      reference <- .generate_realtime_oad_reference(reference, candidates[j, , drop = FALSE])
      if (length(reference$Spectrum) == 0) {
        scored[[j]] <- NULL
        next
      }

      if (use_isotopes) {
        ref_iso_pre <- NULL
        if (exists("lipidomics_get_isotope_pattern_cpp", mode = "function", inherits = TRUE)) {
          ref_iso_pre <- lipidomics_get_isotope_pattern_cpp(database_path, as.integer(candidates$ID[[j]]))
        }
        ref_isotopes_arg <- if (!is.null(ref_iso_pre) && is.matrix(ref_iso_pre) && nrow(ref_iso_pre) > 0) list(ref_iso_pre) else list()
      } else {
        ref_isotopes_arg <- list()
      }

      scored[[j]] <- calculate_score(
        property = property,
        scan = scan,
        scan_isotopes = scan_isotopes_arg,
        reference = reference,
        reference_isotopes = ref_isotopes_arg,
        parameter = parameter,
        annotator_id = as.character(candidates$Table[[j]]),
        source_type = source_type,
        priority = 0L,
        use_ms2 = TRUE
      )
    }

    keep <- !vapply(scored, is.null, logical(1L))
    scored <- scored[keep]
    cand_kept <- candidates[keep, , drop = FALSE]

    if (length(scored) > 0) {
      ord <- order(vapply(scored, function(x) x$TotalScore, numeric(1L)), decreasing = TRUE)
      ord <- ord[seq_len(min(length(ord), top_n))]
      scored <- scored[ord]
      cand_kept <- cand_kept[ord, , drop = FALSE]

      # keep up to top_n rows for this feature (one per retained scored candidate)
      nr <- length(scored)
      rows <- vector("list", nr)
      for (k in seq_len(nr)) {
        s <- scored[[k]]
        rows[[k]] <- data.frame(
          FeatureIndex = idx0,
          PrecursorMz = precursor_mz,
          RT = feature_rt,
          LibraryID = s$LibraryID,
          Name = s$Name,
          AbbreName = if ("AbbreName" %in% names(cand_kept)) as.character(cand_kept$AbbreName[k]) else "",
          InChIKey = s$InChIKey,
          TotalScore = s$TotalScore,
          WeightedDotProduct = s$WeightedDotProduct,
          ReverseDotProduct = s$ReverseDotProduct,
          MatchedPeaksPercentage = s$MatchedPeaksPercentage,
          SourceTable = cand_kept$Table[[k]],
          stringsAsFactors = FALSE
        )
      }
      top_hit_rows[[i]] <- do.call(rbind, rows)
    } else {
      top_hit_rows[[i]] <- NULL
    }

    lipid_match_res[[i]] <- list(
      FeatureIndex = idx0,
      CandidateCount = as.integer(nrow(candidates)),
      Scores = scored,
      CandidateMeta = cand_kept
    )

  }

  top_hit_rows <- top_hit_rows[!vapply(top_hit_rows, is.null, logical(1L))]
  lipid_top_hits <- if (length(top_hit_rows) > 0) do.call(rbind, top_hit_rows) else data.frame()

  mSet@MSnResults[["LipidMatchRes"]] <- lipid_match_res
  mSet@MSnResults[["LipidTopHits"]] <- lipid_top_hits
  mSet@MSnData[["lipidomics_database_path"]] <- database_path

  mSet
}

# Lipid-specific spectra searching wrapper (merged from R/lipid_db_search.R)
LipidSpectraSearching <- function(ConsensusRes, idxs, peak_matrix, ppm_ms1, ppm_ms2,
                                    rt_tol, rt_ms1, scan_ms1, ion_mode = 1,
                                    database_path = "/Volumes/ExtremeSSD/Lipidomics_projects/database_curation/lipidomics_complete.sqlite",
                                    database_table = "spectra_pos",
                                    use_rt = FALSE,
                                    enableNL = FALSE,
                                    useEntropy = FALSE) {
    if (!exists("lipidomics_db_search_precursors_blob_cpp", mode = "function", inherits = TRUE)) {
      stop("lipidomics_db_search_precursors_blob_cpp() not available. Rebuild package to register Rcpp exports.")
    }

    ft_idxs <- ConsensusRes[[1]]
    results <- vector("list", length(idxs))
    for (i in seq_along(idxs)) {
      this_idx <- ft_idxs[[idxs[[i]] + 1]]
      mz_min <- peak_matrix[this_idx + 1, 1]
      mz_max <- peak_matrix[this_idx + 1, 2]
      mz_med <- (mz_min + mz_max) / 2

      spectra_table <- database_table
      if (ion_mode == 1) spectra_table <- "spectra_pos" else spectra_table <- "spectra_neg"

      res_blob <- lipidomics_db_search_precursors_blob_cpp(database_path, spectra_table, as.numeric(c(mz_med)), ppm_ms1, max_candidates = 200, max_per_table = 100)
      if (length(res_blob) >= 1) {
        candidates <- res_blob[[1]]
      } else {
        candidates <- list()
      }
      results[[i]] <- candidates
    }
    return(results)
  }


  #' PerformLipidsDBSearchingBatch
  #'
  #' Batch database searching for lipidomics DB (sqlite). Stores results in mSet@MSnResults$LipidDBMatchRes
  PerformLipidsDBSearchingBatch <- function(mSet = NULL,
                                           ppm1 = 5,
                                           ppm2 = 15,
                                           rt_tol = 0,
                                           database_path = "/Volumes/ExtremeSSD/Lipidomics_projects/database_curation/lipidomics_complete.sqlite",
                                           use_rt = FALSE,
                                           ncores = 1L,
                                           useEntropy = FALSE) {
    require(parallel)
    if (is.null(mSet)) stop("mSet is required. Run PerformSpectrumConsenus() first.")
    if (!file.exists(database_path)) stop("Lipid database not found: ", database_path)

    ConsensusRes <- mSet@MSnResults[["Concensus_spec"]]
    if (is.null(ConsensusRes)) stop("Missing ConsensusRes. Run PerformSpectrumConsenus() first.")

    peak_matrix <- mSet@MSnData[["peak_mtx"]]
    num_fts <- length(ConsensusRes[[1]])
    seq_idx <- seq(0, num_fts - 1)

    if (ncores == 1) {
      res_list <- LipidSpectraSearching(ConsensusRes, seq_idx, peak_matrix, ppm1, ppm2, rt_tol,
                                        mSet@MSnData[["scanrts_ms1"]], mSet@MSnData[["scan_ms1"]],
                                        ion_mode = mSet@MSnData[["ion_mode"]],
                                        database_path = database_path,
                                        use_rt = use_rt,
                                        useEntropy = useEntropy)
      mSet@MSnResults[["LipidDBMatchRes"]] <- res_list
    } else {
      mem_num <- ceiling(num_fts / ncores)
      ft_idx_grps <- lapply(1:ncores, function(x) {
        rs <- seq_idx[(1 + (x - 1) * mem_num):(x * mem_num)]
        rs[!is.na(rs)]
      })
      cl <- makeCluster(getOption("cl.cores", ncores))
      clusterExport(cl, c("LipidSpectraSearching", "ConsensusRes", "peak_matrix", "ppm1", "ppm2", "rt_tol", "database_path", "use_rt", "useEntropy"), envir = environment())
      parts <- parLapply(cl, 1:ncores, function(x, ft_idx_grps, ConsensusRes, peak_matrix, ppm1, ppm2, rt_tol, db_path, use_rt, useEntropy, scanrts_ms1, scan_ms1, ion_mode) {
        LipidSpectraSearching(ConsensusRes, ft_idx_grps[[x]], peak_matrix, ppm1, ppm2, rt_tol,
                              scanrts_ms1, scan_ms1,
                              ion_mode = ion_mode,
                              database_path = db_path,
                              use_rt = use_rt,
                              useEntropy = useEntropy)
      }, ft_idx_grps = ft_idx_grps, ConsensusRes = ConsensusRes, peak_matrix = peak_matrix, ppm1 = ppm1, ppm2 = ppm2, rt_tol = rt_tol, db_path = database_path, use_rt = use_rt, useEntropy = useEntropy, scanrts_ms1 = mSet@MSnData[["scanrts_ms1"]], scan_ms1 = mSet@MSnData[["scan_ms1"]], ion_mode = mSet@MSnData[["ion_mode"]])
      stopCluster(cl)
      res <- do.call(c, parts)
      mSet@MSnResults[["LipidDBMatchRes"]] <- res
    }

    mSet@MSnData[["lipidomics_database_path"]] <- database_path
    return(mSet)
  }

#' PerformLipidSpecDeconvolution
#'
#' Deconvolve lipid DDA spectra with the native C++ deconvolution engine.
#' This function shares the same core deconvolution logic as
#' \code{PerformDDADeconvolution()}, but uses lipidomics sqlite schema
#' (\code{lipidomics_complete.sqlite}) for reference matching.
#'
#' @param mSet mSet object containing imported MSn data from \code{PerformMSnImport()}.
#' @param database_path Path to the lipidomics sqlite database.
#' @param ppm1 Numeric, MS1 precursor tolerance in ppm.
#' @param ppm2 Numeric, MS2 fragment tolerance in ppm.
#' @param sn Numeric, signal-to-noise threshold for deconvolution.
#' @param filtering Numeric, minimum intensity threshold used during preprocessing.
#' @param window_size Numeric, precursor isolation window size.
#' @param intensity_thresh Numeric, threshold for precursor intensity.
#' @param ncores Integer, number of CPU cores.
#' @param decoOn Logical, whether deconvolution is enabled.
#' @param useEntropy Logical, whether entropy similarity is used in candidate selection.
#'
#' @return mSet object with \code{mSet@MSnResults$LipidSpecDecRes} filled.
#' @export
PerformLipidSpecDeconvolution <- function(mSet = NULL,
                                          database_path = "/Volumes/ExtremeSSD/Lipidomics_projects/database_curation/lipidomics_complete.sqlite",
                                          ppm1 = 5,
                                          ppm2 = 15,
                                          sn = 12,
                                          filtering = 2000,
                                          window_size = 1,
                                          intensity_thresh = 1e3,
                                          ncores = 1L,
                                          decoOn = TRUE,
                                          useEntropy = FALSE) {
  require(parallel)
  if (is.null(mSet)) {
    stop("mSet is missing! Please import your MSn data first with function \"PerformMSnImport\" !")
  }
  if (ncores <= 0) {
    warning("At least 1 cpu core required for processing!")
  }
  if (!is.integer(ncores)) {
    stop("\"ncores\" must be an integer!")
  }
  if (!file.exists(database_path)) {
    stop("Your database doesnot exist in this path: ", database_path, ", please check!")
  }
  if (!exists("PerformLipidDDADeco", mode = "function", inherits = TRUE)) {
    stop("PerformLipidDDADeco() was not found. Reinstall/rebuild package to register Rcpp exports.")
  }
  cpp_deco <- get("PerformLipidDDADeco", mode = "function", inherits = TRUE)

  MSn_data <- mSet@MSnData
  idxVec <- seq_along(MSn_data[["scanrts_ms2"]])
  ion_mode <- MSn_data$ion_mode

  nPeaks <- nrow(MSn_data[["peak_mtx"]])
  if (ncores <= nPeaks) {
    Npbatch <- ceiling(nPeaks / ncores)
    batches <- split(c(1:nPeaks), ceiling(c(1:nPeaks) / Npbatch))
    if (length(batches) < ncores) {
      batches <- split(c(1:nPeaks), rep(c(1:ncores), Npbatch)[1:nPeaks])
    }
    col1 <- as.integer(vapply(idxVec, function(x) {rep(x, ncores)}, FUN.VALUE = integer(length = ncores)))
    col2 <- rep(1:ncores, length(idxVec))
  } else {
    batches <- list(seq(nPeaks))
    col1 <- idxVec
    col2 <- rep(1, length(idxVec))
  }

  ls_orches <- list(col1, col2, batches)

  if (ncores > 1) {
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl, c("ls_orches", "MSn_data",
                        "window_size", "ppm1", "ppm2",
                        "sn", "filtering", "intensity_thresh",
                        "ion_mode", "database_path", "decoOn", "useEntropy",
                        "cpp_deco"), envir = environment())

    DecRes <- parLapply(cl,
                        seq(col1),
                        function(x, ls_orches, winSize, ppm1, ppm2, sn, filt,
                                 intensity_thresh, ionmode, db_path, decoOn, useEntropy){
                          x1 <- ls_orches[[1]][x]
                          x2 <- ls_orches[[3]][[ls_orches[[2]][x]]]
                          if (any(x2 == 1)) {showOutput <- TRUE} else {showOutput <- FALSE}
                          if (length(x2) == 1) {
                            peak_matrix <- matrix(MSn_data[["peak_mtx"]][x2, ],
                                                  ncol = length(MSn_data[["peak_mtx"]][x2, ]))
                          } else {
                            peak_matrix <- MSn_data[["peak_mtx"]][x2, ]
                          }

                          prec_mzs <- MSn_data[["precursors"]][[x1]]
                          scant1 <- MSn_data[["scanrts_ms1"]][[x1]]
                          scant2 <- MSn_data[["scanrts_ms2"]][[x1]]
                          scanms1 <- MSn_data[["scan_ms1"]][[x1]]
                          scanms2 <- MSn_data[["scan_ms2"]][[x1]]
                          fileName <- MSn_data[["fileNames"]][x1]
                          thread_int <- x

                          cpp_deco(peak_matrix,
                                   scant1, scant2,
                                   scanms1, scanms2,
                                   prec_mzs, winSize,
                                   ppm1, ppm2,
                                   sn, filt, intensity_thresh,
                                   ionmode, db_path,
                                   decoOn, useEntropy,
                                   showOutput, thread_int, fileName)
                        },
                        ls_orches = ls_orches,
                        winSize = window_size,
                        ppm1 = ppm1,
                        ppm2 = ppm2,
                        sn = sn,
                        filt = filtering,
                        intensity_thresh = intensity_thresh,
                        ionmode = ion_mode,
                        db_path = database_path,
                        decoOn = decoOn,
                        useEntropy = useEntropy)
    stopCluster(cl)
  } else {
    DecRes <- lapply(seq(col1),
                     function(x, ls_orches, winSize, ppm1, ppm2, sn, filt,
                              intensity_thresh, ionmode, db_path, decoOn, useEntropy){
                       x1 <- ls_orches[[1]][x]
                       x2 <- ls_orches[[3]][[ls_orches[[2]][x]]]
                       if (any(x2 == 1)) {showOutput <- TRUE} else {showOutput <- FALSE}
                       if (length(x2) == 1) {
                         peak_matrix <- matrix(MSn_data[["peak_mtx"]][x2, ],
                                               ncol = length(MSn_data[["peak_mtx"]][x2, ]))
                       } else {
                         peak_matrix <- MSn_data[["peak_mtx"]][x2, ]
                       }
                       prec_mzs <- MSn_data[["precursors"]][[x1]]
                       scant1 <- MSn_data[["scanrts_ms1"]][[x1]]
                       scant2 <- MSn_data[["scanrts_ms2"]][[x1]]
                       scanms1 <- MSn_data[["scan_ms1"]][[x1]]
                       scanms2 <- MSn_data[["scan_ms2"]][[x1]]
                       fileName <- MSn_data[["fileNames"]][x1]
                       thread_int <- x

                       cpp_deco(peak_matrix,
                                scant1, scant2,
                                scanms1, scanms2,
                                prec_mzs, winSize,
                                ppm1, ppm2,
                                sn, filt, intensity_thresh,
                                ionmode, db_path,
                                decoOn, useEntropy,
                                showOutput, thread_int, fileName)
                     },
                     ls_orches = ls_orches,
                     winSize = window_size,
                     ppm1 = ppm1,
                     ppm2 = ppm2,
                     sn = sn,
                     filt = filtering,
                     intensity_thresh = intensity_thresh,
                     ionmode = ion_mode,
                     db_path = database_path,
                     decoOn = decoOn,
                     useEntropy = useEntropy)
  }

  cat("\nLipid deconvolution executed successfully!")

  DecResx <- lapply(idxVec, function(x){
    idx0 <- c(ls_orches[[1]] == x)
    DecRes0 <- DecRes[which(idx0)]
    Spectra <- lapply(DecRes0, function(u){u[[1]]})
    Spectra <- do.call(c, Spectra)
    Indicator <- lapply(DecRes0, function(u){u[[2]]})
    Indicator <- do.call(c, Indicator)
    FeatureIdx <- lapply(seq_along(batches), function(u){
      d0 <- batches[[u]]
      d0[DecRes0[[u]][["FeatureIdx"]] + 1] - 1
    })
    FeatureIdx <- do.call(c, FeatureIdx)

    list(Spectra = Spectra,
         Indicator = Indicator,
         FeatureIdx = FeatureIdx)
  })

  mSet@MSnResults[["LipidSpecDecRes"]] <- DecResx
  mSet@MSnData[["lipidomics_database_path"]] <- database_path
  mSet
}
