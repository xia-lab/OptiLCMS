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
                                 doublebond_cutoff) {
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
  position_detected <- detected(positionions)
  db_detected <- detected(dbions)
  class_exists <- class_detected >= class_cutoff && all(!classions$AbsolutelyRequired | classions$IsMatched)
  chain_exists <- chain_detected >= chain_cutoff && all(!chainions$AbsolutelyRequired | chainions$IsMatched)
  position_exists <- position_cutoff <= 0 || (position_detected >= position_cutoff && all(!positionions$AbsolutelyRequired | positionions$IsMatched))
  db_percent <- db_detected / (nrow(dbions) + 1e-10)
  db_high_exists <- nrow(dbhigh) == sum(dbhigh$IsMatched)
  db_exists <- db_high_exists && db_percent > doublebond_cutoff
  db_score <- if (db_exists) db_percent + .eid_matched_coefficient(dbions) else 0
  if (db_exists) {
    high_resolution <- dbions$Resolution[vapply(dbions$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_high") & dbions$IsMatched]
    low_resolution <- dbions$Resolution[vapply(dbions$SpectrumComment, .eid_comment_has, logical(1L), flag = "doublebond_low") & dbions$IsMatched]
    high_mean <- if (length(high_resolution) > 0L) mean(high_resolution) else 0
    low_mean <- if (length(low_resolution) > 0L) mean(low_resolution) else 0
    if (high_mean > low_mean * 1.5) db_score <- db_score + 0.5
  }
  class_score <- if (class_exists && nrow(classions) > 0L) class_detected / nrow(classions) else 0
  chain_score <- if (chain_exists && nrow(chainions) > 0L) chain_detected / nrow(chainions) else 0
  position_score <- if (position_cutoff > 0 && position_exists && nrow(positionions) > 0L) position_detected / nrow(positionions) else 0
  list(TotalScore = class_score + chain_score + position_score + db_score,
       TotalMatchedIonCount = class_detected + chain_detected + position_detected + db_detected)
}

.eieio_characterize <- function(scan, reference, tolerance, mz_begin, mz_end, cutoffs) {
  score <- do.call(.eieio_default_score, c(list(scan, reference, tolerance, mz_begin, mz_end), as.list(cutoffs)))
  list(reference$Name, c(score$TotalScore, score$TotalMatchedIonCount))
}

.eieio_alias <- function(cutoffs) function(scan, reference, tolerance, mz_begin, mz_end) {
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

characterize_pc_eieio <- .eieio_alias(.eieio_diacyl)
characterize_pe_eieio <- .eieio_alias(.eieio_diacyl)
characterize_ps_eieio <- .eieio_alias(.eieio_diacyl)
characterize_pg_eieio <- .eieio_alias(.eieio_diacyl)
characterize_pi_eieio <- .eieio_alias(.eieio_diacyl)
characterize_pa_eieio <- .eieio_alias(.eieio_diacyl)
characterize_dg_eieio <- .eieio_alias(.eieio_diacyl_neutral)
characterize_bmp_eieio <- .eieio_alias(.eieio_diacyl)
characterize_lpc_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lps_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lpe_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lpg_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_lpi_eieio <- .eieio_alias(.eieio_lysophospholipid)
characterize_dgta_eieio <- .eieio_alias(.eieio_diacyl)
characterize_dgts_eieio <- .eieio_alias(.eieio_diacyl)
characterize_ldgta_eieio <- .eieio_alias(.eieio_diacyl)
characterize_ldgts_eieio <- .eieio_alias(.eieio_diacyl)
characterize_sm_eieio <- .eieio_alias(.eieio_ceramide)
characterize_ceramide_eieio <- .eieio_alias(.eieio_ceramide)
characterize_hexcer_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_hex2cer_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_hbmp_eieio <- .eieio_alias(.eieio_hbmp)
characterize_tg_eieio <- .eieio_alias(.eieio_triacyl)
characterize_etherpc_eieio <- .eieio_alias(.eieio_diacyl)
characterize_etherpe_eieio <- .eieio_alias(.eieio_diacyl)
characterize_shexcer_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_gm3_eieio <- .eieio_alias(.eieio_hexceramide)
characterize_ce_eieio <- .eieio_alias(.eieio_single)
characterize_mg_eieio <- .eieio_alias(.eieio_mono)
characterize_car_eieio <- .eieio_alias(c(2, 1, 0, 0.5))
characterize_dmedfahfa_eieio <- .eieio_alias(c(1, 1, 1, 0.5))
characterize_dmedfa_eieio <- .eieio_alias(c(1, 1, 0, 0.5))
characterize_dmedoxfa_eieio <- .eieio_alias(c(1, 1, 0, 0.5))
