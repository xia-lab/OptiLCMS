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

# ─────────────────────────────────────────────────────────────────────────────
# Lipid-class-specific characterizers (EIEIO collision type)
# 
# Each of these functions should:
#   1. Parse the reference structure (chains, double bonds, etc.)
#   2. Calculate diagnostic m/z values for EIEIO fragmentation
#   3. Search the experimental scan for matches
#   4. Return list(lipid_name_annotated, c(percentage, count))
#
# TODO: Implement each characterizer based on EIEIO fragmentation patterns
# ─────────────────────────────────────────────────────────────────────────────

characterize_pc_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific PC fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pe_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific PE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ps_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific PS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific PG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pi_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific PI fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pa_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific PA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific DG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_bmp_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific BMP fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpc_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LPC fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lps_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LPS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpe_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LPE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LPG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpi_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LPI fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dgta_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific DGTA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dgts_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific DGTS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ldgta_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LDGTA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ldgts_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific LDGTS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_sm_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific SM fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ceramide_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific ceramide fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_hexcer_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific HexCer fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_hex2cer_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific Hex2Cer fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_hbmp_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific HBMP fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_tg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific TG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_etherpc_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific EtherPC fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_etherpe_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific EtherPE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_shexcer_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific SHexCer fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_gm3_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific GM3 fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ce_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific CE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_mg_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific MG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_car_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific CAR fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dmedfahfa_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific DMEDFAHFA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dmedfa_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific DMEDFA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dmedoxfa_eieio <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EIEIO-specific DMEDOxFA fragmentation logic
  list(reference$Name, c(0, 0))
}
