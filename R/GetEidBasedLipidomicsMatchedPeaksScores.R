# GetEidBasedLipidomicsMatchedPeaksScores.R
# Translation of MSDIAL EID scoring through EieioMsCharacterizationUtility.

.eid_value <- function(x, names, default = NULL) {
  if (is.null(x)) return(default)
  for (name in names) {
    value <- x[[name]]
    if (!is.null(value) && length(value) > 0L && !is.na(value[1L])) return(value[1L])
  }
  default
}

.eid_peak_rows <- function(spectrum) {
  if (is.null(spectrum) || length(spectrum) == 0L) return(data.frame())
  if (is.data.frame(spectrum)) {
    if (is.null(spectrum$Intensity)) spectrum$Intensity <- 0
    if (is.null(spectrum$SpectrumComment)) spectrum$SpectrumComment <- ""
    if (is.null(spectrum$Comment)) spectrum$Comment <- ""
    if (is.null(spectrum$AbsolutelyRequired)) spectrum$AbsolutelyRequired <- FALSE
    return(spectrum)
  }
  if (is.matrix(spectrum) && ncol(spectrum) >= 2L) {
    return(data.frame(Mass = as.numeric(spectrum[, 1L]), Intensity = as.numeric(spectrum[, 2L]), stringsAsFactors = FALSE))
  }
  if (is.list(spectrum) && !is.null(spectrum$Mass)) {
    return(data.frame(Mass = as.numeric(spectrum$Mass), Intensity = as.numeric(spectrum$Intensity), stringsAsFactors = FALSE))
  }
  if (is.list(spectrum)) {
    rows <- lapply(spectrum, function(peak) {
      if (!is.list(peak)) return(NULL)
      data.frame(
        Mass = as.numeric(.eid_value(peak, c("Mass", "mass", "Mz", "mz"), NA_real_)),
        Intensity = as.numeric(.eid_value(peak, c("Intensity", "intensity"), 0)),
        Comment = as.character(.eid_value(peak, c("Comment", "comment"), "")),
        SpectrumComment = as.character(.eid_value(peak, c("SpectrumComment", "spectrum_comment"), "")),
        AbsolutelyRequired = as.logical(.eid_value(peak, c("IsAbsolutelyRequiredFragmentForAnnotation", "AbsolutelyRequired"), FALSE)),
        stringsAsFactors = FALSE)
    })
    rows <- rows[!vapply(rows, is.null, logical(1L))]
    if (length(rows) > 0L) return(do.call(rbind, rows))
  }
  data.frame()
}

.eid_comment_has <- function(comment, flag) {
  if (is.null(comment) || length(comment) == 0L || is.na(comment)) return(FALSE)
  flag_value <- switch(flag,
    metaboliteclass = 0x4000L,
    acylchain = 0x8000L,
    doublebond = 0x10000L,
    snposition = 0x20000L,
    doublebond_high = 0x40000L,
    doublebond_low = 0x80000L,
    NA_integer_)
  numeric_comment <- suppressWarnings(as.numeric(comment[1L]))
  if (!is.na(numeric_comment) && !is.na(flag_value)) {
    return(bitwAnd(as.integer(numeric_comment), flag_value) != 0L)
  }
  grepl(paste0("(^|[,;|[:space:]])", flag, "([,;|[:space:]]|$)"), as.character(comment), perl = TRUE, ignore.case = TRUE)
}

.eid_get_matched_peaks <- function(scan, reference, tolerance, mz_begin, mz_end) {
  exp_peaks <- .eid_peak_rows(scan$Spectrum)
  ref_peaks <- .eid_peak_rows(reference$Spectrum)
  if (nrow(exp_peaks) == 0L || nrow(ref_peaks) == 0L) return(data.frame())
  exp_peaks <- exp_peaks[order(exp_peaks$Mass), , drop = FALSE]
  ref_peaks <- ref_peaks[order(ref_peaks$Mass), , drop = FALSE]
  ref_peaks <- ref_peaks[ref_peaks$Mass >= mz_begin & ref_peaks$Mass <= mz_end, , drop = FALSE]
  if (nrow(ref_peaks) == 0L) return(ref_peaks)
  out <- ref_peaks
  out$Resolution <- 0
  out$OriginalIntensity <- 0
  out$IsMatched <- FALSE
  for (i in seq_len(nrow(out))) {
    index <- which(abs(exp_peaks$Mass - out$Mass[i]) < tolerance)
    if (length(index) > 0L) {
      out$Resolution[i] <- sum(exp_peaks$Intensity[index], na.rm = TRUE)
      out$OriginalIntensity[i] <- out$Resolution[i]
      out$IsMatched[i] <- TRUE
    }
  }
  out
}

.eid_matched_coefficient <- function(peaks) {
  peaks <- peaks[peaks$IsMatched, , drop = FALSE]
  if (nrow(peaks) < 2L) return(0)
  resolution <- as.numeric(peaks$Resolution)
  intensity <- as.numeric(peaks$Intensity)
  if (sd(resolution) == 0 || sd(intensity) == 0) return(0)
  sum((resolution - mean(resolution)) * (intensity - mean(intensity))) / sqrt(sum((resolution - mean(resolution))^2) * sum((intensity - mean(intensity))^2))
}

.eid_default_score <- function(scan, reference, tolerance, mz_begin, mz_end, class_cutoff, chain_cutoff, position_cutoff, doublebond_cutoff) {
  matched <- .eid_get_matched_peaks(scan, reference, tolerance, mz_begin, mz_end)
  empty <- list(ClassIonsDetected = 0L, ChainIonsDetected = 0L, PositionIonsDetected = 0L, DoubleBondIonsDetected = 0L, IsClassIonsExisted = FALSE, IsChainIonsExisted = FALSE, IsPositionIonsExisted = FALSE, IsDoubleBondIonsExisted = FALSE, TotalScore = 0, TotalMatchedIonCount = 0)
  if (nrow(matched) == 0L) return(empty)

  classify <- function(flag) vapply(matched$SpectrumComment, .eid_comment_has, logical(1L), flag = flag)
  classions <- matched[classify("metaboliteclass"), , drop = FALSE]
  chainions <- matched[classify("acylchain"), , drop = FALSE]
  positionions <- matched[classify("snposition"), , drop = FALSE]
  dbions <- matched[classify("doublebond"), , drop = FALSE]
  dbhigh <- dbions[vapply(dbions$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_high"), , drop = FALSE]

  class_detected <- sum(classions$IsMatched)
  chain_detected <- sum(chainions$IsMatched)
  position_detected <- sum(positionions$IsMatched)
  db_detected <- sum(dbions$IsMatched)
  class_exists <- all(!classions$AbsolutelyRequired | classions$IsMatched) && class_detected >= class_cutoff
  chain_exists <- all(!chainions$AbsolutelyRequired | chainions$IsMatched) && chain_detected >= chain_cutoff
  position_exists <- position_cutoff <= 0 || (all(!positionions$AbsolutelyRequired | positionions$IsMatched) && position_detected >= position_cutoff)

  db_matched <- dbions[dbions$IsMatched, , drop = FALSE]
  db_percent <- db_detected / (nrow(dbions) + 1e-10)
  db_high_exists <- nrow(dbhigh) == sum(dbhigh$IsMatched)
  is_db <- db_high_exists && db_percent > doublebond_cutoff
  high <- db_matched[vapply(db_matched$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_high"), , drop = FALSE]
  low <- db_matched[vapply(db_matched$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_low"), , drop = FALSE]
  high_average <- if (nrow(high) > 0L) mean(high$Resolution) else 0
  low_average <- if (nrow(low) > 0L) mean(low$Resolution) else 0
  high_low_bonus <- if (high_average > low_average * 1.5) 0.5 else 0
  db_score <- if (is_db) .eid_matched_coefficient(dbions) + db_percent else 0
  if (is_db) db_score <- db_score + high_low_bonus

  class_score <- if (class_exists && nrow(classions) > 0L) class_detected / nrow(classions) else 0
  chain_score <- if (chain_exists && nrow(chainions) > 0L) chain_detected / nrow(chainions) else 0
  position_score <- if (position_exists && nrow(positionions) > 0L) position_detected / nrow(positionions) else 0
  list(ClassIonsDetected = class_detected, ChainIonsDetected = chain_detected, PositionIonsDetected = position_detected, DoubleBondIonsDetected = db_detected, IsClassIonsExisted = class_exists, IsChainIonsExisted = chain_exists, IsPositionIonsExisted = position_exists, IsDoubleBondIonsExisted = is_db, DoubleBondMatchedPercent = db_percent, ClassIonScore = class_score, ChainIonScore = chain_score, PositionIonScore = position_score, DoubleBondIonScore = db_score, TotalScore = class_score + chain_score + position_score + db_score, TotalMatchedIonCount = class_detected + chain_detected + position_detected + db_detected)
}

.eid_characterize <- function(scan, reference, tolerance, mz_begin, mz_end, group) {
  cutoffs <- switch(group, diacyl = c(1, 0, 1, 0.5), alkyl = c(1, 0, 1, 0.5), triacyl = c(1, 0, 1, 0.5), ceramide = c(1, 0, 0, 0.5), monoacyl = c(1, 0, 0, 0.5), single = c(1, 0, 0, 0.5))
  score <- do.call(.eid_default_score, c(list(scan, reference, tolerance, mz_begin, mz_end), as.list(cutoffs)))
  list(reference$Name, c(score$TotalScore, score$TotalMatchedIonCount), score)
}

extract_lipid_class <- function(lipid_name) sub("[ /].*$", "", trimws(as.character(lipid_name)[1L]))

get_eid_based_lipid_molecule_annotation_result <- function(scan, reference, tolerance, mz_begin, mz_end) {
  lipid_name <- as.character(if (is.null(reference$Name)) "" else reference$Name)
  lipid_class <- extract_lipid_class(lipid_name)
  diacyl <- c("PC", "PE", "PS", "PG", "PI", "PA", "DG", "BMP", "LPC", "LPS", "LPE", "LPG", "LPI", "DGTA", "DGTS", "LDGTA", "LDGTS", "DMEDFAHFA", "PC_d5", "PE_d5", "PS_d5", "PG_d5", "PI_d5", "LPC_d5", "LPS_d5", "LPE_d5", "LPG_d5", "LPI_d5", "DG_d5")
  ceramide <- c("SM", "Cer_NS", "Cer_NDS", "Cer_AS", "Cer_ADS", "Cer_BS", "Cer_BDS", "Cer_NP", "Cer_AP", "Cer_ABP", "HexCer_NS", "Hex2Cer", "SHexCer", "GM3", "SM_d9", "Cer_NS_d7")
  single <- c("CAR", "DMEDFA", "DMEDOxFA", "CE", "CE_d7")
  group <- if (lipid_class %in% diacyl) "diacyl" else if (lipid_class == "MG") "monoacyl" else if (lipid_class %in% single) "single" else if (lipid_class %in% ceramide) "ceramide" else if (lipid_class %in% c("HBMP", "TG", "TG_d5")) "triacyl" else if (lipid_class %in% c("EtherPC", "EtherPE")) "alkyl" else NULL
  if (is.null(group)) return(list(lipid_name, c(0, 0)))
  .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, group)[1:2]
}

get_eid_based_lipidomics_matched_peaks_scores <- function(scan, reference, tolerance, mz_begin, mz_end) {
  get_eid_based_lipid_molecule_annotation_result(scan, reference, tolerance, mz_begin, mz_end)[[2L]]
}

get_eid_default_score <- .eid_default_score
get_eid_matched_peaks <- .eid_get_matched_peaks
characterize_eid_diacylglycerols <- function(scan, reference, tolerance, mz_begin, mz_end) .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, "diacyl")[1:2]
characterize_eid_monoacylglycerols <- function(scan, reference, tolerance, mz_begin, mz_end) .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, "monoacyl")[1:2]
characterize_eid_single_acyl_chain <- function(scan, reference, tolerance, mz_begin, mz_end) .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, "single")[1:2]
characterize_eid_ceramides <- function(scan, reference, tolerance, mz_begin, mz_end) .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, "ceramide")[1:2]
characterize_eid_triacylglycerols <- function(scan, reference, tolerance, mz_begin, mz_end) .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, "triacyl")[1:2]
characterize_eid_alkyl_acyl_glycerols <- function(scan, reference, tolerance, mz_begin, mz_end) .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, "alkyl")[1:2]

# Compatibility names retained from the initialized class-specific framework.
.eid_alias <- function(group) function(scan, reference, tolerance, mz_begin, mz_end) {
  .eid_characterize(scan, reference, tolerance, mz_begin, mz_end, group)[1:2]
}
characterize_pc_eid <- .eid_alias("diacyl")
characterize_pe_eid <- .eid_alias("diacyl")
characterize_ps_eid <- .eid_alias("diacyl")
characterize_pg_eid <- .eid_alias("diacyl")
characterize_pi_eid <- .eid_alias("diacyl")
characterize_pa_eid <- .eid_alias("diacyl")
characterize_dg_eid <- .eid_alias("diacyl")
characterize_bmp_eid <- .eid_alias("diacyl")
characterize_lpc_eid <- .eid_alias("diacyl")
characterize_lps_eid <- .eid_alias("diacyl")
characterize_lpe_eid <- .eid_alias("diacyl")
characterize_lpg_eid <- .eid_alias("diacyl")
characterize_lpi_eid <- .eid_alias("diacyl")
characterize_dgta_eid <- .eid_alias("diacyl")
characterize_dgts_eid <- .eid_alias("diacyl")
characterize_ldgta_eid <- .eid_alias("diacyl")
characterize_ldgts_eid <- .eid_alias("diacyl")
characterize_sm_eid <- .eid_alias("ceramide")
characterize_ceramide_eid <- .eid_alias("ceramide")
characterize_hexcer_eid <- .eid_alias("ceramide")
characterize_hex2cer_eid <- .eid_alias("ceramide")
characterize_hbmp_eid <- .eid_alias("triacyl")
characterize_tg_eid <- .eid_alias("triacyl")
characterize_etherpc_eid <- .eid_alias("alkyl")
characterize_etherpe_eid <- .eid_alias("alkyl")
characterize_shexcer_eid <- .eid_alias("ceramide")
characterize_gm3_eid <- .eid_alias("ceramide")
characterize_ce_eid <- .eid_alias("single")
characterize_mg_eid <- .eid_alias("monoacyl")
characterize_car_eid <- .eid_alias("single")
characterize_dmedfahfa_eid <- .eid_alias("diacyl")
characterize_dmedfa_eid <- .eid_alias("single")
characterize_dmedoxfa_eid <- .eid_alias("single")
