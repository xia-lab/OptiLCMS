# ─────────────────────────────────────────────────────────────────────────────
# GetEidBasedLipidomicsMatchedPeaksScores.R
# 
# Translation of MsScanMatching.GetEidBasedLipidomicsMatchedPeaksScores
# and GetEidBasedLipidMoleculeAnnotationResult from C#
#
# Purpose: Analyze MS/MS spectra using EID (Electron Ionization Dissociation)
#          collision energy for lipid characterization
#
# These functions parse the reference lipid name, identify the lipid class,
# and dispatch to class-specific characterization routines that match diagnostic
# fragments from the experimental scan against theoretical EID expectations.
#
# Returns: double[2] with [1] = MatchedPeaksPercentage, [2] = MatchedPeaksCount
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# get_eid_based_lipid_molecule_annotation_result
# 
# Dispatcher function that routes to lipid-class-specific characterizers.
# This mimics the large switch statement in C#'s GetEidBasedLipidMoleculeAnnotationResult.
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
get_eid_based_lipid_molecule_annotation_result <- function(scan, reference, 
                                                           tolerance, mz_begin, mz_end) {
  # Parse lipid name to extract class and chain information
  lipid_name <- reference$Name
  lipid_class <- extract_lipid_class(lipid_name)
  
  # Dispatch to class-specific characterizer
  # Note: EID has a different routing pattern than EIEIO in some lipid classes
  result <- switch(lipid_class,
    "PC"        = characterize_pc_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PE"        = characterize_pe_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PS"        = characterize_ps_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PG"        = characterize_pg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PI"        = characterize_pi_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PA"        = characterize_pa_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DG"        = characterize_dg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "BMP"       = characterize_bmp_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPC"       = characterize_lpc_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPS"       = characterize_lps_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPE"       = characterize_lpe_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPG"       = characterize_lpg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPI"       = characterize_lpi_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DGTA"      = characterize_dgta_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DGTS"      = characterize_dgts_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LDGTA"     = characterize_ldgta_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LDGTS"     = characterize_ldgts_eid(scan, reference, tolerance, mz_begin, mz_end),
    "SM"        = characterize_sm_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NS"    = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NDS"   = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_AS"    = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_ADS"   = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_BS"    = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_BDS"   = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NP"    = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_AP"    = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_ABP"   = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "HexCer_NS" = characterize_hexcer_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Hex2Cer"   = characterize_hex2cer_eid(scan, reference, tolerance, mz_begin, mz_end),
    "HBMP"      = characterize_hbmp_eid(scan, reference, tolerance, mz_begin, mz_end),
    "TG"        = characterize_tg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "EtherPC"   = characterize_etherpc_eid(scan, reference, tolerance, mz_begin, mz_end),
    "EtherPE"   = characterize_etherpe_eid(scan, reference, tolerance, mz_begin, mz_end),
    "SHexCer"   = characterize_shexcer_eid(scan, reference, tolerance, mz_begin, mz_end),
    "GM3"       = characterize_gm3_eid(scan, reference, tolerance, mz_begin, mz_end),
    "CE"        = characterize_ce_eid(scan, reference, tolerance, mz_begin, mz_end),
    "MG"        = characterize_mg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "CAR"       = characterize_car_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDFAHFA" = characterize_dmedfahfa_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDFA"    = characterize_dmedfa_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDOxFA"  = characterize_dmedoxfa_eid(scan, reference, tolerance, mz_begin, mz_end),
    # Isotope-labeled variants
    "PC_d5"     = characterize_pc_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PE_d5"     = characterize_pe_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PS_d5"     = characterize_ps_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PG_d5"     = characterize_pg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "PI_d5"     = characterize_pi_eid(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NS_d7" = characterize_ceramide_eid(scan, reference, tolerance, mz_begin, mz_end),
    "SM_d9"     = characterize_sm_eid(scan, reference, tolerance, mz_begin, mz_end),
    "TG_d5"     = characterize_tg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "CE_d7"     = characterize_ce_eid(scan, reference, tolerance, mz_begin, mz_end),
    "DG_d5"     = characterize_dg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPC_d5"    = characterize_lpc_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPS_d5"    = characterize_lps_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPE_d5"    = characterize_lpe_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPG_d5"    = characterize_lpg_eid(scan, reference, tolerance, mz_begin, mz_end),
    "LPI_d5"    = characterize_lpi_eid(scan, reference, tolerance, mz_begin, mz_end),
    # Default fallback
    list(lipid_name, c(0, 0))
  )
  
  if (is.null(result)) {
    result <- list(lipid_name, c(0, 0))
  }
  return(result)
}

# ─────────────────────────────────────────────────────────────────────────────
# get_eid_based_lipidomics_matched_peaks_scores
#
# Main entry point: returns only the scores from the full annotation result
#
# Returns: double[2] with [1] = MatchedPeaksPercentage, [2] = MatchedPeaksCount
# ─────────────────────────────────────────────────────────────────────────────
get_eid_based_lipidomics_matched_peaks_scores <- function(scan, reference, 
                                                          tolerance, mz_begin, mz_end) {
  result <- get_eid_based_lipid_molecule_annotation_result(
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
# Lipid-class-specific characterizers (EID collision type)
# 
# Each of these functions should:
#   1. Parse the reference structure (chains, double bonds, etc.)
#   2. Calculate diagnostic m/z values for EID fragmentation
#   3. Search the experimental scan for matches
#   4. Return list(lipid_name_annotated, c(percentage, count))
#
# TODO: Implement each characterizer based on EID fragmentation patterns
# ─────────────────────────────────────────────────────────────────────────────

characterize_pc_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific PC fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pe_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific PE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ps_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific PS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pg_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific PG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pi_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific PI fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_pa_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific PA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dg_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific DG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_bmp_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific BMP fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpc_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LPC fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lps_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LPS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpe_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LPE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpg_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LPG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_lpi_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LPI fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dgta_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific DGTA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dgts_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific DGTS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ldgta_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LDGTA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ldgts_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific LDGTS fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_sm_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific SM fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ceramide_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific ceramide fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_hexcer_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific HexCer fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_hex2cer_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific Hex2Cer fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_hbmp_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific HBMP fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_tg_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific TG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_etherpc_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific EtherPC fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_etherpe_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific EtherPE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_shexcer_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific SHexCer fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_gm3_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific GM3 fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_ce_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific CE fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_mg_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific MG fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_car_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific CAR fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dmedfahfa_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific DMEDFAHFA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dmedfa_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific DMEDFA fragmentation logic
  list(reference$Name, c(0, 0))
}

characterize_dmedoxfa_eid <- function(scan, reference, tolerance, mz_begin, mz_end) {
  # TODO: Implement EID-specific DMEDOxFA fragmentation logic
  list(reference$Name, c(0, 0))
}
