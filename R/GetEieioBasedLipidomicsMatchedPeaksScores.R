# ─────────────────────────────────────────────────────────────────────────────
# GetEieioBasedLipidomicsMatchedPeaksScores.R
# 
# Translation of MsScanMatching.GetEieioBasedLipidomicsMatchedPeaksScores
# and GetEieioBasedLipidMoleculeAnnotationResult from C#
#
# Purpose: Analyze MS/MS spectra using EIEIO (Electron-Impact Excitation of 
#          Ions from Organics) collision energy for lipid characterization
#
# These functions parse the reference lipid name, identify the lipid class,
# and dispatch to class-specific characterization routines that match diagnostic
# fragments from the experimental scan against theoretical expectations.
#
# Returns: double[2] with [1] = MatchedPeaksPercentage, [2] = MatchedPeaksCount
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# get_eieio_based_lipid_molecule_annotation_result
# 
# Dispatcher function that routes to lipid-class-specific characterizers.
# This mimics the large switch statement in C#'s GetEieioBasedLipidMoleculeAnnotationResult.
#
# Arguments:
#   scan      : MS/MS scan with $Spectrum (list of peaks)
#   reference : Reference lipid with $Name (lipid name string), $Spectrum
#   tolerance : MS/MS mass tolerance (ppm or Da)
#   mz_begin  : Lower mass range boundary
#   mz_end    : Upper mass range boundary
#
# Returns: list(lipid_name, scores) where scores[1] = percentage, scores[2] = count
# ─────────────────────────────────────────────────────────────────────────────
get_eieio_based_lipid_molecule_annotation_result <- function(scan, reference, 
                                                             tolerance, mz_begin, mz_end) {
  # Parse lipid name to extract class and chain information
  lipid_name <- reference$Name
  lipid_class <- extract_lipid_class(lipid_name)
  
  # Dispatch to class-specific characterizer
  result <- switch(lipid_class,
    "PC"        = characterize_pc_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PE"        = characterize_pe_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PS"        = characterize_ps_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PG"        = characterize_pg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PI"        = characterize_pi_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PA"        = characterize_pa_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DG"        = characterize_dg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "BMP"       = characterize_bmp_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPC"       = characterize_lpc_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPS"       = characterize_lps_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPE"       = characterize_lpe_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPG"       = characterize_lpg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPI"       = characterize_lpi_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DGTA"      = characterize_dgta_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DGTS"      = characterize_dgts_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LDGTA"     = characterize_ldgta_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LDGTS"     = characterize_ldgts_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "SM"        = characterize_sm_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NS"    = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NDS"   = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_AS"    = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_ADS"   = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_BS"    = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_BDS"   = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NP"    = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_AP"    = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_ABP"   = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "HexCer_NS" = characterize_hexcer_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Hex2Cer"   = characterize_hex2cer_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "HBMP"      = characterize_hbmp_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "TG"        = characterize_tg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "EtherPC"   = characterize_etherpc_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "EtherPE"   = characterize_etherpe_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "SHexCer"   = characterize_shexcer_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "GM3"       = characterize_gm3_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "CE"        = characterize_ce_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "MG"        = characterize_mg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "CAR"       = characterize_car_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDFAHFA" = characterize_dmedfahfa_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDFA"    = characterize_dmedfa_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDOxFA"  = characterize_dmedoxfa_eieio(scan, reference, tolerance, mz_begin, mz_end),
    # Isotope-labeled variants (D5, D7, D9, etc.)
    "PC_d5"     = characterize_pc_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PE_d5"     = characterize_pe_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PS_d5"     = characterize_ps_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PG_d5"     = characterize_pg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "PI_d5"     = characterize_pi_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NS_d7" = characterize_ceramide_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "SM_d9"     = characterize_sm_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "TG_d5"     = characterize_tg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "CE_d7"     = characterize_ce_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "DG_d5"     = characterize_dg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPC_d5"    = characterize_lpc_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPS_d5"    = characterize_lps_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPE_d5"    = characterize_lpe_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPG_d5"    = characterize_lpg_eieio(scan, reference, tolerance, mz_begin, mz_end),
    "LPI_d5"    = characterize_lpi_eieio(scan, reference, tolerance, mz_begin, mz_end),
    # Default fallback
    list(lipid_name, c(0, 0))
  )
  
  if (is.null(result)) {
    result <- list(lipid_name, c(0, 0))
  }
  return(result)
}

# ─────────────────────────────────────────────────────────────────────────────
# get_eieio_based_lipidomics_matched_peaks_scores
#
# Main entry point: returns only the scores from the full annotation result
#
# Returns: double[2] with [1] = MatchedPeaksPercentage, [2] = MatchedPeaksCount
# ─────────────────────────────────────────────────────────────────────────────
get_eieio_based_lipidomics_matched_peaks_scores <- function(scan, reference, 
                                                            tolerance, mz_begin, mz_end) {
  if (nrow(.eid_peak_rows(scan$Spectrum)) == 0L ||
      nrow(.eid_peak_rows(reference$Spectrum)) == 0L) {
    return(c(-1, -1))
  }
  result <- get_eieio_based_lipid_molecule_annotation_result(
    scan, reference, tolerance, mz_begin, mz_end)
  # result is list(lipid_name, scores_vector)
  return(result[[2]])
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: extract_lipid_class
# Parses lipid name to extract the lipid class (e.g., "PC" from "PC 36:2")
# ─────────────────────────────────────────────────────────────────────────────
extract_lipid_class <- function(lipid_name) {
  # Remove leading/trailing whitespace
  lipid_name <- trimws(lipid_name)
  # Extract class (word before first space or slash)
  class_part <- sub(" .*", "", lipid_name)
  class_part <- sub("/.*", "", class_part)
  return(class_part)
}

# Shared EIEIO scoring path. The C# characterizers differ primarily in these
# cutoffs; fragment matching and double-bond handling are shared.
.eieio_default_score <- function(scan, reference, tolerance, mz_begin, mz_end,
                                 class_cutoff, chain_cutoff, position_cutoff,
                                 doublebond_cutoff, position_diagnostic = NULL) {
  matched <- .eid_get_matched_peaks(scan, reference, tolerance, mz_begin, mz_end)
  if (nrow(matched) == 0L) return(list(TotalScore = 0, TotalMatchedIonCount = 0))
  classify <- function(flag) vapply(matched$SpectrumComment, .eid_comment_has, logical(1L), flag = flag)
  classions <- matched[classify("metaboliteclass"), , drop = FALSE]
  chainions <- matched[classify("acylchain"), , drop = FALSE]
  positionions <- matched[classify("snposition"), , drop = FALSE]
  dbions <- matched[classify("doublebond"), , drop = FALSE]
  dbhigh <- dbions[vapply(dbions$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_high"), , drop = FALSE]
  detected <- function(peaks) sum(peaks$IsMatched)
  class_detected <- detected(classions)
  chain_detected <- detected(chainions)
  position_detected <- if (position_cutoff > 0) detected(positionions) else 0L
  db_detected <- detected(dbions)
  class_exists <- class_detected >= class_cutoff && all(!classions$AbsolutelyRequired | classions$IsMatched)
  chain_exists <- chain_detected >= chain_cutoff && all(!chainions$AbsolutelyRequired | chainions$IsMatched)
  minimum_intensity <- min(.eid_peak_rows(scan$Spectrum)$Intensity, na.rm = TRUE)
  if (!is.finite(minimum_intensity)) minimum_intensity <- 0
  position_advanced <- .eieio_diagnostic_ions_exist(positionions, position_diagnostic, minimum_intensity)
  position_exists <- position_cutoff <= 0 || (position_advanced && position_detected >= position_cutoff && all(!positionions$AbsolutelyRequired | positionions$IsMatched))
  db_percent <- db_detected / (nrow(dbions) + 1e-10)
  reference_high <- .eid_peak_rows(reference$Spectrum)
  reference_high <- reference_high[vapply(reference_high$SpectrumComment,
                                          .eid_comment_has, logical(1L),
                                          flag = "doublebond_high"), , drop = FALSE]
  db_high_exists <- .eieio_diagnostic_fragments_exist(
    dbhigh[dbhigh$IsMatched, , drop = FALSE], reference_high, tolerance,
    minimum_intensity)
  db_exists <- db_high_exists && db_percent > doublebond_cutoff
  db_score <- if (db_exists) db_percent + .eid_matched_coefficient(dbions) else 0
  hgrain_exists <- any(grepl("+H_p3", dbions$Comment, fixed = TRUE) &
                       vapply(dbions$SpectrumComment, .eid_comment_has,
                              logical(1L), flag = "doublebond_high"))
  if (db_exists) {
    high_resolution <- dbions$Resolution[vapply(dbions$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_high") & dbions$IsMatched]
    low_resolution <- dbions$Resolution[vapply(dbions$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_low") & dbions$IsMatched]
    high_mean <- if (length(high_resolution) > 0L) mean(high_resolution) else 0
    low_mean <- if (length(low_resolution) > 0L) mean(low_resolution) else 0
    if (!hgrain_exists && high_mean > low_mean * 1.5) db_score <- db_score + 0.5
    if (.eieio_hgrain_bonus(dbions) && db_exists) db_score <- db_score + 0.5
  }
  class_score <- if (class_exists && nrow(classions) > 0L) class_detected / nrow(classions) else 0
  chain_score <- if (chain_exists && nrow(chainions) > 0L) chain_detected / nrow(chainions) else 0
  position_score <- if (position_cutoff > 0 && position_exists && nrow(positionions) > 0L) position_detected / nrow(positionions) else 0
  list(TotalScore = class_score + chain_score + position_score + db_score,
       TotalMatchedIonCount = class_detected + chain_detected + position_detected + db_detected)
}

.eieio_diagnostic_fragments_exist <- function(peaks, expected, tolerance,
                                              minimum_intensity,
                                              minimum_intensity_factor = 2) {
  if (nrow(expected) == 0L) return(TRUE)
  if (nrow(peaks) == 0L) return(FALSE)
  threshold <- minimum_intensity * minimum_intensity_factor
  misses <- 0L
  for (index in seq_len(nrow(expected))) {
    present <- any(abs(peaks$Mass - expected$Mass[index]) < tolerance &
                   peaks$Resolution > 1e-7 &
                   peaks$Resolution > threshold &
                   peaks$OriginalIntensity > 12)
    if (!present) {
      misses <- misses + 1L
      if ((nrow(expected) < 5L && misses == 1L) ||
          (nrow(expected) >= 5L && misses == 2L)) return(FALSE)
    }
  }
  TRUE
}

    .eieio_diagnostic_ions_exist <- function(peaks, diagnostic_ions, minimum_intensity,
                     minimum_intensity_factor = 2) {
      if (is.null(diagnostic_ions) || nrow(diagnostic_ions) == 0L) return(TRUE)
      if (nrow(peaks) == 0L) return(FALSE)
      all(vapply(seq_len(nrow(diagnostic_ions)), function(index) {
        tolerance <- diagnostic_ions$MzTolerance[index]
        cutoff <- diagnostic_ions$IonAbundanceCutOff[index]
        any(abs(peaks$Mass - diagnostic_ions$Mz[index]) < tolerance &
        peaks$Resolution > cutoff &
        peaks$Resolution > minimum_intensity * minimum_intensity_factor &
        peaks$OriginalIntensity > 12)
      }, logical(1L)))
    }

.eieio_hgrain_bonus <- function(doublebond_ions) {
  hgrain <- doublebond_ions[grepl("+H_p3", doublebond_ions$Comment, fixed = TRUE), , drop = FALSE]
  if (nrow(hgrain) == 0L) return(FALSE)
  hgrain <- hgrain[order(hgrain$Mass, decreasing = TRUE), , drop = FALSE]
  count <- 0L
  peak_average <- 0
  for (index in rev(seq_len(nrow(hgrain)))) {
    count <- count + 1L
    peak <- hgrain[index, , drop = FALSE]
    if (.eid_comment_has(peak$SpectrumComment, "doublebond_high")) {
      peak_average <- peak_average / count
      return(peak$Resolution > peak_average * 3)
    }
    peak_average <- peak_average + peak$Resolution
  }
  FALSE
}

.eieio_characterize <- function(scan, reference, tolerance, mz_begin, mz_end, cutoffs) {
  score <- do.call(.eieio_default_score, c(list(scan, reference, tolerance, mz_begin, mz_end), as.list(cutoffs)))
  list(reference$Name, c(score$TotalScore, score$TotalMatchedIonCount))
}

.eieio_alias <- function(cutoffs) function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_characterize(scan, reference, tolerance, mz_begin, mz_end, cutoffs)
}

.eieio_chain_summary <- function(reference) {
  name <- as.character(reference$Name %||% "")
  matches <- regmatches(name, gregexpr("[0-9]+:[0-9]+", name, perl = TRUE))[[1L]]
  if (length(matches) == 0L || identical(matches, character(0))) return(matrix(numeric(), ncol = 2L))
  matrix(as.numeric(unlist(strsplit(matches, ":", fixed = TRUE))), ncol = 2L, byrow = TRUE)
}

.eieio_adduct_name <- function(reference) {
  adduct <- reference$AdductType
  if (is.list(adduct)) adduct <- adduct$AdductIonName %||% adduct$Name
  as.character(adduct %||% "")
}

.eieio_adduct_mass <- function(reference) {
  adduct <- reference$AdductType
  if (is.list(adduct)) {
    adduct <- adduct$AdductIonAccurateMass %||% adduct$AccurateMass %||% adduct$Mass
  }
  value <- suppressWarnings(as.numeric(adduct))
  if (length(value) == 0L || is.na(value[1L])) 0 else value[1L]
}

.eieio_fragment_intensity <- function(spectrum, mass, tolerance) {
  peaks <- .eid_peak_rows(spectrum)
  if (nrow(peaks) == 0L) return(0)
  matched <- abs(peaks$Mass - mass) < tolerance
  if (!any(matched)) 0 else max(peaks$Intensity[matched], na.rm = TRUE)
}

.eieio_fragment1_greater <- function(spectrum, fragment1, fragment2, tolerance) {
  .eieio_fragment_intensity(spectrum, fragment1, tolerance) >
    .eieio_fragment_intensity(spectrum, fragment2, tolerance)
}

.eieio_dgta_fragment_masses <- function(reference) {
  carbon <- 12.00000000000
  hydrogen <- 1.00782503207
  nitrogen <- 14.00307400480
  oxygen <- 15.99491461956
  adduct <- .eieio_adduct_mass(reference)
  c(dgta = 7 * carbon + 13 * hydrogen + nitrogen + 2 * oxygen + adduct,
    dgts = 6 * carbon + 11 * hydrogen + nitrogen + 2 * oxygen + adduct)
}

.eieio_dynamic_characterize <- function(scan, reference, tolerance, mz_begin, mz_end,
                                        cutoffs, equal_chain_cutoff = FALSE,
                                        position_zero_adducts = character(),
                                        triacyl_chain_cutoff = FALSE) {
  chains <- .eieio_chain_summary(reference)
  if (equal_chain_cutoff && nrow(chains) == 2L && all(chains[1L, ] == chains[2L, ])) {
    cutoffs[2L] <- 1
  }
  if (triacyl_chain_cutoff && nrow(chains) == 3L) {
    equal12 <- all(chains[1L, ] == chains[2L, ])
    equal23 <- all(chains[2L, ] == chains[3L, ])
    equal13 <- all(chains[1L, ] == chains[3L, ])
    cutoffs[2L] <- if (equal12 && equal23) 1 else if (equal12 || equal23 || equal13) 2 else 3
  }
  if (isTRUE(position_zero_adducts) ||
      identical(.eieio_adduct_name(reference), position_zero_adducts)) {
    cutoffs[3L] <- 0
  }
  .eieio_characterize(scan, reference, tolerance, mz_begin, mz_end, cutoffs)
}

.eieio_diacyl <- c(2, 2, 1, 0.5)
.eieio_diacyl_neutral <- c(0, 2, 1, 0.5)
.eieio_hbmp <- c(1, 2, 1, 0.5)
.eieio_lysophospholipid <- c(2, 1, 1, 0.5)
.eieio_single <- c(1, 1, 0, 2)
.eieio_mono <- c(2, 1, 1, 2)
.eieio_ceramide <- c(1, 1, 0, 0.5)
.eieio_hexceramide <- c(2, 2, 1, 0.5)
.eieio_triacyl <- c(0, 2, 1, 0.5)

characterize_pc_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_diacyl, equal_chain_cutoff = TRUE)
}
characterize_pe_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_diacyl, equal_chain_cutoff = TRUE,
                              position_zero_adducts = "[M+H]+")
}
characterize_ps_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              c(3, 2, 1, 0.5), equal_chain_cutoff = TRUE,
                              position_zero_adducts = "[M+H]+")
}
characterize_pg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              c(3, 2, 1, 0.5), equal_chain_cutoff = TRUE,
                              position_zero_adducts = "[M+NH4]+")
}
characterize_pi_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_diacyl, equal_chain_cutoff = TRUE,
                              position_zero_adducts = "[M+NH4]+")
}
characterize_pa_eieio <- .eieio_alias(.eieio_diacyl)
characterize_dg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_diacyl_neutral, equal_chain_cutoff = TRUE,
                              position_zero_adducts = "[M+NH4]+")
}
characterize_bmp_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_diacyl, equal_chain_cutoff = TRUE,
                              position_zero_adducts = "[M+H]+")
}
characterize_lpc_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lps_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lpe_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lpg_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lpi_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_dgta_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  masses <- .eieio_dgta_fragment_masses(reference)
  cutoff <- if (.eieio_fragment1_greater(scan$Spectrum, masses[["dgts"]], masses[["dgta"]], tolerance)) 10 else 2
  .eieio_characterize(scan, reference, tolerance, mz_begin, mz_end, c(cutoff, 2, 1, 0.5))
}
characterize_dgts_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  masses <- .eieio_dgta_fragment_masses(reference)
  cutoff <- if (.eieio_fragment1_greater(scan$Spectrum, masses[["dgta"]], masses[["dgts"]], tolerance)) 10 else 2
  .eieio_characterize(scan, reference, tolerance, mz_begin, mz_end, c(cutoff, 2, 1, 0.5))
}
characterize_ldgta_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  masses <- .eieio_dgta_fragment_masses(reference)
  cutoff <- if (.eieio_fragment1_greater(scan$Spectrum, masses[["dgts"]], masses[["dgta"]], tolerance)) 10 else 1
  .eieio_characterize(scan, reference, tolerance, mz_begin, mz_end, c(cutoff, 2, 1, 0.5))
}
characterize_ldgts_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  masses <- .eieio_dgta_fragment_masses(reference)
  cutoff <- if (.eieio_fragment1_greater(scan$Spectrum, masses[["dgta"]], masses[["dgts"]], tolerance)) 10 else 1
  .eieio_characterize(scan, reference, tolerance, mz_begin, mz_end, c(cutoff, 2, 1, 0.5))
}
characterize_sm_eieio <- .eieio_alias(.eieio_ceramide)
characterize_ceramide_eieio <- .eieio_alias(.eieio_ceramide)
characterize_hexcer_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_hex2cer_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_hbmp_eieio <- .eieio_alias(.eieio_hbmp)
characterize_tg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_triacyl, position_zero_adducts = "[M+NH4]+",
                              triacyl_chain_cutoff = TRUE)
}
characterize_etherpc_eieio <- .eieio_alias(.eieio_diacyl)
characterize_etherpe_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .eieio_dynamic_characterize(scan, reference, tolerance, mz_begin, mz_end,
                              .eieio_diacyl, position_zero_adducts = "[M+H]+")
}
characterize_shexcer_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_gm3_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_ce_eieio <- .eieio_alias(.eieio_single)
characterize_mg_eieio <- .eieio_alias(.eieio_mono)
characterize_car_eieio <- .eieio_alias(c(2, 1, 0, 0.5))
characterize_dmedfahfa_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  matched <- .eid_get_matched_peaks(scan, reference, tolerance, mz_begin, mz_end)
  position <- matched[vapply(matched$SpectrumComment, .eid_comment_has, logical(1L), flag = "snposition"), , drop = FALSE]
  diagnostic <- if (nrow(position) > 0L) {
    data.frame(Mz = position$Mass[1L], MzTolerance = 0.05,
               IonAbundanceCutOff = 10, stringsAsFactors = FALSE)
  } else NULL
  score <- .eieio_default_score(scan, reference, tolerance, mz_begin, mz_end,
                                1, 1, 1, 0.5, diagnostic)
  list(reference$Name, c(score$TotalScore, score$TotalMatchedIonCount))
}
characterize_dmedfa_eieio <- .eieio_alias(c(1, 1, 0, 0.5))
characterize_dmedoxfa_eieio <- .eieio_alias(c(1, 1, 0, 0.5))
