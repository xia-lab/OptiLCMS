# ─────────────────────────────────────────────────────────────────────────────
# GetOadBasedLipidomicsMatchedPeaksScores.R
# 
# Translation of MsScanMatching.GetOadBasedLipidomicsMatchedPeaksScores
# and GetOadBasedLipidMoleculeAnnotationResult from C#
#
# Purpose: Analyze MS/MS spectra using OAD (Ozone-Induced Dissociation / O-PTAD)
#          collision energy for lipid characterization with double-bond position
#          localization capability
#
# These functions parse the reference lipid name, identify the lipid class,
# and dispatch to class-specific characterization routines that match diagnostic
# fragments from the experimental scan against theoretical OAD expectations.
# OAD is particularly useful for localizing double-bond positions in fatty acid chains.
#
# Returns: double[2] with [1] = MatchedPeaksPercentage, [2] = MatchedPeaksCount
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# get_oad_based_lipid_molecule_annotation_result
# 
# Dispatcher function that routes to lipid-class-specific characterizers.
# This mimics the large switch statement in C#'s GetOadBasedLipidMoleculeAnnotationResult.
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
get_oad_based_lipid_molecule_annotation_result <- function(scan, reference, 
                                                           tolerance, mz_begin, mz_end) {
  # Parse lipid name to extract class and chain information
  lipid_name <- reference$Name
  lipid_class <- extract_lipid_class(lipid_name)
  
  # Dispatch to class-specific characterizer
  result <- switch(lipid_class,
    "PC"        = characterize_pc_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PE"        = characterize_pe_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PS"        = characterize_ps_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PG"        = characterize_pg_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PI"        = characterize_pi_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PA"        = characterize_pa_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DG"        = characterize_dg_oad(scan, reference, tolerance, mz_begin, mz_end),
    "BMP"       = characterize_bmp_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPC"       = characterize_lpc_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPS"       = characterize_lps_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPE"       = characterize_lpe_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPG"       = characterize_lpg_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPI"       = characterize_lpi_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DGTA"      = characterize_dgta_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DGTS"      = characterize_dgts_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LDGTA"     = characterize_ldgta_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LDGTS"     = characterize_ldgts_oad(scan, reference, tolerance, mz_begin, mz_end),
    "SM"        = characterize_sm_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Sph"       = characterize_sphingosine_oad(scan, reference, tolerance, mz_begin, mz_end),
    "SPB"       = characterize_sphingosine_oad(scan, reference, tolerance, mz_begin, mz_end),
    "ASM"       = characterize_acylsm_oad(scan, reference, tolerance, mz_begin, mz_end),
    "NAGly"     = characterize_nacylglyoxfa_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NS"    = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NDS"   = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_AS"    = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_ADS"   = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_BS"    = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_BDS"   = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NP"    = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_AP"    = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_ABP"   = characterize_ceramide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "HexCer_NS" = characterize_hexcer_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Hex2Cer"   = characterize_hex2cer_oad(scan, reference, tolerance, mz_begin, mz_end),
    "HBMP"      = characterize_hbmp_oad(scan, reference, tolerance, mz_begin, mz_end),
    "TG"        = characterize_tg_oad(scan, reference, tolerance, mz_begin, mz_end),
    "EtherLPC"  = characterize_etherlysopc_oad(scan, reference, tolerance, mz_begin, mz_end),
    "EtherLPE"  = characterize_etherlysope_oad(scan, reference, tolerance, mz_begin, mz_end),
    "EtherPC"   = characterize_etherpc_oad(scan, reference, tolerance, mz_begin, mz_end),
    "EtherPE"   = characterize_etherpe_oad(scan, reference, tolerance, mz_begin, mz_end),
    "SHexCer"   = characterize_shexcer_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PI-Cer"    = characterize_picermide_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PhytoSph"  = characterize_sphingomyelin_phyto_oad(scan, reference, tolerance, mz_begin, mz_end),
    "GM3"       = characterize_gm3_oad(scan, reference, tolerance, mz_begin, mz_end),
    "CE"        = characterize_ce_oad(scan, reference, tolerance, mz_begin, mz_end),
    "MG"        = characterize_mg_oad(scan, reference, tolerance, mz_begin, mz_end),
    "CAR"       = characterize_car_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDFAHFA" = characterize_dmedfahfa_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDFA"    = characterize_dmedfa_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DMEDOxFA"  = characterize_dmedoxfa_oad(scan, reference, tolerance, mz_begin, mz_end),
    # Isotope-labeled variants
    "PC_d5"     = characterize_pc_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PE_d5"     = characterize_pe_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PS_d5"     = characterize_ps_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PG_d5"     = characterize_pg_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "PI_d5"     = characterize_pi_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "Cer_NS_d7" = characterize_cer_ns_d7_oad(scan, reference, tolerance, mz_begin, mz_end),
    "SM_d9"     = characterize_sm_d9_oad(scan, reference, tolerance, mz_begin, mz_end),
    "TG_d5"     = characterize_tg_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "CE_d7"     = characterize_ce_d7_oad(scan, reference, tolerance, mz_begin, mz_end),
    "DG_d5"     = characterize_dg_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPC_d5"    = characterize_lpc_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPS_d5"    = characterize_lps_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPE_d5"    = characterize_lpe_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPG_d5"    = characterize_lpg_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    "LPI_d5"    = characterize_lpi_d5_oad(scan, reference, tolerance, mz_begin, mz_end),
    # Default fallback
    list(lipid_name, c(0, 0), 0L, NULL)
  )
  
  if (is.null(result)) {
    result <- list(lipid_name, c(0, 0), 0L, NULL)
  }
  return(result)
}

# ─────────────────────────────────────────────────────────────────────────────
# C#-equivalent helpers used by the score-only entry point
# ─────────────────────────────────────────────────────────────────────────────
.oad_is_compared_available <- function(scan, reference) {
  !is.null(scan) && !is.null(reference) &&
    !is.null(scan$Spectrum) && !is.null(reference$Spectrum) &&
    length(scan$Spectrum) > 0L && length(reference$Spectrum) > 0L
}

.oad_spectrum_columns <- function(spectrum) {
  if (is.data.frame(spectrum)) {
    mass_name <- intersect(c("Mass", "mass", "mz", "Mz"), names(spectrum))[1L]
    intensity_name <- intersect(c("Intensity", "intensity", "abundance"), names(spectrum))[1L]
    if (length(mass_name) == 0L || length(intensity_name) == 0L ||
      is.na(mass_name) || is.na(intensity_name)) {
      stop("Spectrum data frame must contain mass/mz and intensity columns.")
    }
    return(list(mass = spectrum[[mass_name]], intensity = spectrum[[intensity_name]]))
  }

  if (is.list(spectrum) && !is.null(spectrum$Mass) && !is.null(spectrum$Intensity)) {
    return(list(mass = spectrum$Mass, intensity = spectrum$Intensity))
  }

  if (is.matrix(spectrum) && ncol(spectrum) >= 2L) {
    return(list(mass = spectrum[, 1L], intensity = spectrum[, 2L]))
  }

  stop("Spectrum must be a data frame, a list with Mass/Intensity, or a two-column matrix.")
}

.oad_get_matched_peaks_scores <- function(scan, reference, bin, mass_begin, mass_end) {
  if (!.oad_is_compared_available(scan, reference)) return(c(-1, -1))

  measured <- .oad_spectrum_columns(scan$Spectrum)
  library <- .oad_spectrum_columns(reference$Spectrum)
  if (length(measured$mass) == 0L || length(library$mass) == 0L) return(c(-1, -1))

  min_mz <- library$mass[1L]
  max_mz <- library$mass[length(library$mass)]
  if (mass_begin > min_mz) min_mz <- mass_begin
  if (max_mz > mass_end) max_mz <- mass_end

  focused_mz <- min_mz
  max_library_intensity <- max(library$intensity)
  measured_index <- 1L
  library_index <- 1L
  matched_count <- 0L
  library_count <- 0L

  while (focused_mz <= max_mz) {
    library_sum <- 0
    for (index in seq.int(library_index, length(library$mass))) {
      if (library$mass[index] < focused_mz - bin) next
      if (library$mass[index] < focused_mz + bin) {
        library_sum <- library_sum + library$intensity[index]
      } else {
        library_index <- index
        break
      }
    }
    if (library_sum >= 0.01 * max_library_intensity) library_count <- library_count + 1L

    measured_sum <- 0
    for (index in seq.int(measured_index, length(measured$mass))) {
      if (measured$mass[index] < focused_mz - bin) next
      if (measured$mass[index] < focused_mz + bin) {
        measured_sum <- measured_sum + measured$intensity[index]
      } else {
        measured_index <- index
        break
      }
    }
    if (measured_sum > 0 && library_sum >= 0.01 * max_library_intensity) {
      matched_count <- matched_count + 1L
    }

    if (focused_mz + bin > library$mass[length(library$mass)]) break
    if (library_index > length(library$mass)) break
    focused_mz <- library$mass[library_index]
  }

  if (library_count == 0L) return(c(0, 0))
  c(matched_count / library_count, matched_count)
}

# ─────────────────────────────────────────────────────────────────────────────
# get_oad_based_lipidomics_matched_peaks_scores
#
# Main entry point corresponding to C#'s
# GetLipidomicsMoleculerSpeciesLevelAnnotationPeaksScoresForOAD.
# ─────────────────────────────────────────────────────────────────────────────
get_oad_based_lipidomics_matched_peaks_scores <- function(scan, reference, 
                                                          tolerance, mz_begin, mz_end) {
  if (!.oad_is_compared_available(scan, reference)) return(c(-1, -1))

  result <- .oad_get_matched_peaks_scores(scan, reference, tolerance, mz_begin, mz_end)
  compound_class <- reference$CompoundClass %||% extract_lipid_class(reference$Name)
  comment <- reference$Comment %||% ""

  if (comment != "SPLASH" && !compound_class %in% c("Unknown", "Others")) {
    annotation <- get_oad_based_lipid_molecule_annotation_result(
      scan, reference, tolerance, mz_begin, mz_end)
    if (!is.null(annotation) && length(annotation) >= 3L) {
      annotation_level <- annotation[[3L]] %||% 0
      if (annotation_level == 1L) {
        result[1L] <- if (compound_class == "SM" &&
                           grepl("3O|O3", reference$Name)) 2 else 1
      } else if (annotation_level >= 2L) {
        result[1L] <- 2
      }
    }
  } else if (comment == "SPLASH" && compound_class == "CE") {
    return(c(-1, -1))
  }
  result
}

`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L || is.na(value[1L])) fallback else value[1L]
}

.oad_field <- function(object, names, default = NULL) {
  for (name in names) {
    value <- object[[name]]
    if (!is.null(value) && length(value) > 0L && !is.na(value[1L])) return(value[1L])
  }
  default
}

.oad_reference_molecule <- function(reference) {
  molecule <- reference$Molecule %||% reference
  lipid_name <- .oad_field(molecule, c("LipidName", "Name"), reference$Name %||% "")
  lipid_class <- .oad_field(molecule, c("LipidClass", "CompoundClass"), extract_lipid_class(lipid_name))
  chain_matches <- gregexpr("[0-9]+:[0-9]+", lipid_name, perl = TRUE)
  chain_text <- regmatches(lipid_name, chain_matches)[[1L]]
  chain_values <- do.call(rbind, lapply(chain_text, function(value) {
    as.integer(strsplit(value, ":", fixed = TRUE)[[1L]])
  }))
  if (length(chain_values) == 0L) chain_values <- matrix(c(0L, 0L), ncol = 2L)

  total_carbon <- as.integer(.oad_field(molecule, "TotalCarbonCount", sum(chain_values[, 1L])))
  total_double <- as.integer(.oad_field(molecule, "TotalDoubleBondCount", sum(chain_values[, 2L])))
  total_oxidized <- as.integer(.oad_field(molecule, c("TotalOxidizedCount", "TotalOxidized"), 0L))
  if (nrow(chain_values) > 1L) {
    first_chain <- if (nrow(chain_values) > 2L) 1L else nrow(chain_values) - 1L
    sn1_carbon <- as.integer(.oad_field(molecule, "Sn1CarbonCount", chain_values[first_chain, 1L]))
    sn1_double <- as.integer(.oad_field(molecule, "Sn1DoubleBondCount", chain_values[first_chain, 2L]))
    sn2_index <- if (nrow(chain_values) > 2L) 2L else nrow(chain_values)
    sn2_carbon <- as.integer(.oad_field(molecule, "Sn2CarbonCount", chain_values[sn2_index, 1L]))
    sn2_double <- as.integer(.oad_field(molecule, "Sn2DoubleBondCount", chain_values[sn2_index, 2L]))
    sn3_carbon <- as.integer(.oad_field(molecule, c("Sn3CarbonCount", "Sn3Carbon"), if (nrow(chain_values) > 2L) chain_values[3L, 1L] else 0L))
    sn3_double <- as.integer(.oad_field(molecule, c("Sn3DoubleBondCount", "Sn3DoubleBond"), if (nrow(chain_values) > 2L) chain_values[3L, 2L] else 0L))
  } else {
    sn1_carbon <- as.integer(.oad_field(molecule, "Sn1CarbonCount", floor(total_carbon / 2L)))
    sn1_double <- as.integer(.oad_field(molecule, "Sn1DoubleBondCount", floor(total_double / 2L)))
    sn2_carbon <- total_carbon - sn1_carbon
    sn2_double <- total_double - sn1_double
    sn3_carbon <- 0L
    sn3_double <- 0L
  }

  list(
    name = lipid_name,
    class = lipid_class,
    mz = as.numeric(.oad_field(molecule, c("Mz", "PrecursorMz"), reference$Mz %||% reference$PrecursorMz)),
    total_carbon = total_carbon,
    total_double = total_double,
    total_oxidized = total_oxidized,
    sn1_carbon = sn1_carbon,
    sn1_double = sn1_double,
    sn2_carbon = sn2_carbon,
    sn2_double = sn2_double,
    sn3_carbon = sn3_carbon,
    sn3_double = sn3_double,
    adduct = if (is.list(reference$Adduct)) reference$Adduct else list(
      IonMode = .oad_field(molecule, "IonMode", reference$IonMode %||% "Positive"),
      AdductIonName = .oad_field(molecule, "AdductIonName", reference$AdductIonName %||% reference$Adduct %||% "[M+H]+")
    )
  )
}

.oad_annotation_result <- function(reference, result) {
  if (is.null(result)) return(NULL)
  list(
    result$LipidName %||% reference$Name,
    c(0, 0),
    as.integer(result$AnnotationLevel %||% 0L),
    result
  )
}

.oad_default_characterization_result <- function(scan, reference, tolerance, mz_begin,
                     mz_end, level = 1L) {
  score <- .oad_get_matched_peaks_scores(scan, reference, tolerance, mz_begin, mz_end)
  molecule <- .oad_reference_molecule(reference)
  list(reference$Name %||% "", score, as.integer(level),
    LipidName = molecule$name, AnnotationLevel = as.integer(level),
    TotalCarbonCount = molecule$total_carbon,
    TotalDoubleBondCount = molecule$total_double,
    Score = score[1L], MatchedPeakCount = score[2L])
}

.oad_mass <- function(carbon = 0, hydrogen = 0, nitrogen = 0, oxygen = 0, phosphorus = 0) {
  carbon * MASS_DIFF$carbon + hydrogen * MASS_DIFF$hydrogen +
    nitrogen * MASS_DIFF$nitrogen + oxygen * MASS_DIFF$oxygen +
    phosphorus * MASS_DIFF$phosphorus
}

.oad_fatty_acid_product_ion <- function(carbon, double_bond) {
  acyl_chain_mass(carbon, double_bond) + MASS_DIFF$oxygen - PROTON
}

.oad_spectrum_vectors <- function(scan) {
  spectrum <- scan$Spectrum
  if (is.data.frame(spectrum)) {
    mass_name <- intersect(c("Mass", "mass", "mz", "Mz"), names(spectrum))[1L]
    intensity_name <- intersect(c("Intensity", "intensity", "abundance"), names(spectrum))[1L]
    return(list(mass = as.numeric(spectrum[[mass_name]]),
                intensity = as.numeric(spectrum[[intensity_name]])))
  }
  if (is.list(spectrum) && !is.null(spectrum$Mass)) {
    return(list(mass = as.numeric(spectrum$Mass), intensity = as.numeric(spectrum$Intensity)))
  }
  list(mass = vapply(spectrum, `[[`, numeric(1), "Mass"),
       intensity = vapply(spectrum, `[[`, numeric(1), "Intensity"))
}

.oad_peak <- function(scan, tolerance, diagnostic_mz, threshold = 0) {
  spectrum <- .oad_spectrum_vectors(scan)
  any(spectrum$intensity > threshold & abs(spectrum$mass - diagnostic_mz) < tolerance)
}

.oad_fragments <- function(scan, tolerance, fragments, threshold = 0) {
  spectrum <- .oad_spectrum_vectors(scan)
  hits <- vapply(fragments, function(diagnostic_mz) {
    any(spectrum$intensity > threshold & abs(spectrum$mass - diagnostic_mz) < tolerance)
  }, logical(1))
  list(count = sum(hits), average_intensity = if (any(hits)) {
    mean(vapply(fragments[hits], function(diagnostic_mz) {
      max(spectrum$intensity[abs(spectrum$mass - diagnostic_mz) < tolerance])
    }, numeric(1)))
  } else 0)
}

.oad_call_oad_characterizer <- function(scan, reference, tolerance, fn, mode = c("two_chain", "lyso", "ether_lyso", "single_chain", "single_chain_oxidized", "phyto", "three_chain", "two_chain_oxidized", "range_two_chain")) {
  mode <- match.arg(mode)
  molecule <- .oad_reference_molecule(reference)
  if (length(molecule$mz) != 1L || !is.finite(molecule$mz) || !exists(fn, mode = "function")) return(NULL)
  if (mode == "two_chain") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$sn1_carbon, molecule$total_carbon - molecule$sn1_carbon,
      molecule$sn1_double, molecule$total_double - molecule$sn1_double, molecule$adduct
    ))
  } else if (mode == "range_two_chain") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      as.integer(.oad_field(molecule, "min_sn_carbon", molecule$sn1_carbon)),
      as.integer(.oad_field(molecule, "max_sn_carbon", molecule$sn1_carbon)),
      as.integer(.oad_field(molecule, "min_sn_double", molecule$sn1_double)),
      as.integer(.oad_field(molecule, "max_sn_double", molecule$sn1_double)), molecule$adduct
    ))
  } else if (mode == "lyso") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$sn1_carbon, molecule$sn1_double, molecule$adduct
    ))
  } else if (mode == "ether_lyso") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$sn1_carbon, molecule$sn1_carbon, molecule$sn1_double,
      molecule$sn1_double, molecule$adduct
    ))
  } else if (mode == "three_chain") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$sn1_carbon, molecule$sn2_carbon %||% 0L,
      molecule$sn3_carbon,
      molecule$sn1_double, molecule$sn2_double %||% 0L,
      molecule$sn3_double,
      molecule$adduct
    ))
  } else if (mode == "two_chain_oxidized") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$sn1_carbon, molecule$total_carbon - molecule$sn1_carbon,
      molecule$sn1_double, molecule$total_double - molecule$sn1_double,
      molecule$adduct, molecule$total_oxidized
    ))
  } else if (mode == "single_chain_oxidized") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$total_oxidized, molecule$adduct
    ))
  } else if (mode == "phyto") {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$sn1_carbon, molecule$sn1_carbon, molecule$sn1_double,
      molecule$sn1_double, molecule$adduct
    ))
  } else {
    result <- do.call(fn, list(
      scan, tolerance, molecule$mz, molecule$total_carbon, molecule$total_double,
      molecule$adduct
    ))
  }
  .oad_annotation_result(reference, result)
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
# Lipid-class-specific characterizers (OAD collision type)
#
# The existing R MS/MS characterizers are used for switch branches whose C#
# implementation delegates to LipidMsmsCharacterization. OAD-only branches
# remain explicit in the map below and return no annotation until their
# LipidOadMsmsCharacterization method has been translated.
# ─────────────────────────────────────────────────────────────────────────────

characterize_pc_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylcholine_oad")
}

characterize_pe_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylethanolamine_oad")
}

characterize_ps_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylserine_oad")
}

characterize_pg_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylglycerol_oad")
}

characterize_pi_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylinositol_oad")
}

characterize_pa_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_dg_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_dag_oad")
}

characterize_bmp_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_bismonoacylglycerophosphate_oad")
}

characterize_lpc_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysopc_oad", "lyso")
}

characterize_lps_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysops_oad", "lyso")
}

characterize_lpe_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysope_oad", "lyso")
}

characterize_lpg_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysopg_oad", "lyso")
}

characterize_lpi_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysopi_oad", "lyso")
}

characterize_etherlysopc_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_etherlysopc_oad", "ether_lyso")
}

characterize_etherlysope_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_dgta_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  characterize_dgts_oad(scan, reference, tolerance, mz_begin, mz_end)
}

characterize_dgts_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_dgts_oad", "range_two_chain")
}

characterize_ldgta_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  characterize_ldgts_oad(scan, reference, tolerance, mz_begin, mz_end)
}

characterize_ldgts_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_ldgts_oad", "range_two_chain")
}

characterize_sm_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  if (grepl("O3|3O", reference$Name %||% "")) {
    return(characterize_sphingomyelin_phyto_oad(scan, reference, tolerance, mz_begin, mz_end))
  }
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_sphingomyelin_oad")
}

characterize_sphingosine_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_sphingosine_oad", "single_chain")
}

characterize_acylsm_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_acylsm_oad")
}

characterize_nacylglyoxfa_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_nacylglyoxfa_oad", "single_chain_oxidized")
}

characterize_ceramide_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_hexcer_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_hex2cer_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_hbmp_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_bismonoacylglycerophosphate_oad")
}

characterize_tg_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_etherpc_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_etherpc_oad")
}

characterize_etherpe_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_etherpe_oad")
}

characterize_shexcer_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_shexcer_oad", "two_chain_oxidized")
}

characterize_sphingomyelin_phyto_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_sphingomyelin_phyto_oad", "phyto")
}

characterize_picermide_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_picermide_oad", "two_chain_oxidized")
}

characterize_gm3_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_ce_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_cholesteryl_ester_oad", "single_chain")
}

characterize_mg_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_car_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_acylcarnitine_oad", "single_chain")
}

characterize_dmedfahfa_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_dmedfa_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

characterize_dmedoxfa_oad <- function(scan, reference, tolerance, mz_begin, mz_end) {
  .oad_default_characterization_result(scan, reference, tolerance, mz_begin, mz_end, 1L)
}

.oad_deuterium <- 2.01410177812
.oad_cd2 <- 12 + 2 * .oad_deuterium
.oad_mass_d <- function(carbon = 0, hydrogen = 0, deuterium = 0, nitrogen = 0,
                        oxygen = 0, phosphorus = 0) {
  .oad_mass(carbon, hydrogen, nitrogen, oxygen, phosphorus) + deuterium * .oad_deuterium
}

.oad_isotope_result <- function(label, level, theoretical_mz, adduct, candidates = list(),
                                suffix = "") {
  value <- .oad_direct_result(label, level, theoretical_mz, adduct, candidates)
  value$Candidates <- candidates
  value$Suffix <- suffix
  value
}

.oad_isotope_two_chain <- function(scan, tolerance, theoretical_mz, total_carbon,
                                   total_double_bond, sn1_carbon, sn2_carbon,
                                   sn1_double, sn2_double, adduct, label,
                                   fragments, candidate_label = label,
                                   threshold = 0.1, required = 2L) {
  hits <- .oad_fragments(scan, tolerance, fragments, threshold)
  candidates <- if (hits$count >= required) list(list(
    LipidName = candidate_label, Sn1Carbon = sn1_carbon,
    Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
    Sn2DoubleBond = sn2_double, AverageIntensity = hits$average_intensity
  )) else list()
  .oad_isotope_result(label, 2L, theoretical_mz, adduct, candidates)
}

judge_if_phosphatidylcholine_d5_oad <- function(scan, tolerance, theoretical_mz,
                                                total_carbon, total_double_bond,
                                                sn1_carbon, sn2_carbon, sn1_double,
                                                sn2_double, adduct) {
  if (!.oad_is_compared_available(scan, list(Spectrum = scan$Spectrum))) return(NULL)
  c5h14no4p <- .oad_mass_d(5, 14, 0, 1, 4, 1)
  c3h9n <- .oad_mass(3, 9, 1)
  gly_c <- .oad_mass_d(8, 13, 5, 1, 4, 1)
  gly_o <- .oad_mass_d(7, 13, 3, 1, 5, 1)
  name <- adduct$AdductIonName %||% ""
  result <- function(candidates) .oad_isotope_result("PC_d5", 2L, theoretical_mz, adduct, candidates)
  if (name == "[M+H]+") {
    if (sn1_carbon < 10L || sn2_carbon < 10L || sn1_double > 6L || sn2_double > 6L) return(NULL)
    if (!.oad_peak(scan, tolerance, c5h14no4p + PROTON, 3) ||
        .oad_peak(scan, tolerance, theoretical_mz - c3h9n, 10) ||
        .oad_peak(scan, tolerance, theoretical_mz - 141.019094261, 5)) return(NULL)
    f <- c(theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen,
           theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen - H2O,
           theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen,
           theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen - H2O)
    hits <- .oad_fragments(scan, tolerance, f, 0.1)
    ch2 <- c(f[2L] - .oad_cd2, f[4L] - .oad_cd2)
    candidates <- if (hits$count >= 2L && .oad_fragments(scan, tolerance, ch2, 1)$count > 0L) list(list(
      Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
      Sn2DoubleBond = sn2_double, AverageIntensity = hits$average_intensity)) else list()
    return(result(candidates))
  }
  if (name == "[M+Na]+") {
    if (!.oad_peak(scan, tolerance, theoretical_mz - c5h14no4p, 3) ||
        !.oad_peak(scan, tolerance, theoretical_mz - c3h9n, 3) ||
        !.oad_peak(scan, tolerance, c5h14no4p + NA_MASS, 3) ||
        !.oad_peak(scan, tolerance, gly_c + NA_MASS, 3) || !.oad_peak(scan, tolerance, gly_o + NA_MASS, 3)) return(NULL)
    candidates <- list()
    for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
      carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1
      if (double2 >= 7L) next
      f <- c(theoretical_mz - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen,
             theoretical_mz - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen - H2O,
             theoretical_mz - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen,
             theoretical_mz - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen - H2O)
      hits <- .oad_fragments(scan, tolerance, f, 0.1)
      ch2 <- c(f[2L] + MASS_DIFF$hydrogen - .oad_cd2, f[4L] + MASS_DIFF$hydrogen - .oad_cd2)
      if (hits$count >= 2L && .oad_fragments(scan, tolerance, ch2, 1)$count > 0L) candidates <- c(candidates, list(list(
        Sn1Carbon = carbon1, Sn1DoubleBond = double1, Sn2Carbon = carbon2,
        Sn2DoubleBond = double2, AverageIntensity = hits$average_intensity)))
    }
    return(result(candidates))
  }
  if (name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(scan, tolerance, diagnostic, 10) ||
        ((name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) && .oad_peak(scan, tolerance, theoretical_mz - 60.021129369, 10))) return(NULL)
    f <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    return(result(if (.oad_fragments(scan, tolerance, f, 0)$count == 2L) list(list(
      Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
      Sn2DoubleBond = sn2_double)) else list()))
  }
  NULL
}

judge_if_phosphatidylethanolamine_d5_oad <- function(scan, tolerance, theoretical_mz,
                                                     total_carbon, total_double_bond,
                                                     sn1_carbon, sn2_carbon, sn1_double,
                                                     sn2_double, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL)
  name <- adduct$AdductIonName %||% ""; c2h8no4p <- .oad_mass(2, 8, 1, 4, 1)
  gly_c <- .oad_mass_d(5, 7, 5, 1, 4, 1); gly_o <- .oad_mass_d(4, 7, 3, 1, 5, 1)
  result <- function(candidates) .oad_isotope_result("PE_d5", 2L, theoretical_mz, adduct, candidates)
  if (name == "[M+H]+") {
    if (!.oad_peak(scan, tolerance, theoretical_mz - c2h8no4p, 2.5)) return(NULL)
    f <- c(acyl_chain_mass(sn1_carbon, sn1_double) - ELECTRON, acyl_chain_mass(sn2_carbon, sn2_double) - ELECTRON)
    hits <- .oad_fragments(scan, tolerance, f, 0.1)
    ch2 <- c(theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen - H2O - .oad_mass(1, 2),
             theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen - H2O - .oad_mass(1, 2))
    return(result(if (hits$count == 2L && .oad_fragments(scan, tolerance, ch2, 1)$count > 0L) list(list(
      Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
      Sn2DoubleBond = sn2_double, AverageIntensity = hits$average_intensity)) else list()))
  }
  if (name == "[M+Na]+") {
    if (!.oad_peak(scan, tolerance, theoretical_mz - c2h8no4p, 3) || !.oad_peak(scan, tolerance, gly_c + NA_MASS, 3) || !.oad_peak(scan, tolerance, gly_o + NA_MASS, 3)) return(NULL)
    candidates <- list()
    for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
      carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1; if (double2 >= 7L) next
      f <- c(theoretical_mz - acyl_chain_mass(carbon1, double1) - MASS_DIFF$hydrogen,
             theoretical_mz - acyl_chain_mass(carbon1, double1) - MASS_DIFF$hydrogen - MASS_DIFF$oxygen,
             theoretical_mz - acyl_chain_mass(carbon2, double2) - MASS_DIFF$hydrogen,
             theoretical_mz - acyl_chain_mass(carbon2, double2) - MASS_DIFF$hydrogen - MASS_DIFF$oxygen)
      hits <- .oad_fragments(scan, tolerance, f, 0.1)
      ch2 <- c(f[2L] - .oad_cd2, f[4L] - .oad_cd2)
      if (hits$count >= 2L && .oad_fragments(scan, tolerance, ch2, 1)$count > 0L) candidates <- c(candidates, list(list(
        Sn1Carbon = carbon1, Sn1DoubleBond = double1, Sn2Carbon = carbon2,
        Sn2DoubleBond = double2, AverageIntensity = hits$average_intensity)))
    }
    return(result(candidates))
  }
  if (name == "[M-H]-") {
    if (!.oad_peak(scan, tolerance, 196.03803, 5) || .oad_peak(scan, tolerance, 152.995833871, 5)) return(NULL)
    f <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    if (.oad_fragments(scan, tolerance, f, 1)$count != 2L) return(NULL)
    return(.oad_isotope_result("PE_d5", 2L, theoretical_mz, adduct))
  }
  NULL
}

.oad_phospholipid_d5_judge <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond,
                                       sn1_carbon, sn2_carbon, sn1_double, sn2_double, adduct, label) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL)
  name <- adduct$AdductIonName %||% ""; candidates <- list()
  if (label == "PG_d5") {
    gly_c <- .oad_mass_d(6, 8, 5, 0, 6, 1); gly_o <- .oad_mass_d(5, 8, 3, 0, 7, 1); header <- .oad_mass(3, 9, 0, 6, 1)
    if (name == "[M+NH4]+") {
      diagnostic <- theoretical_mz - header - .oad_mass(0, 3, 1)
      if (!.oad_peak(scan, tolerance, diagnostic, 10)) return(NULL)
      f <- diagnostic - c(acyl_chain_mass(sn1_carbon, sn1_double), acyl_chain_mass(sn2_carbon, sn2_double)) + MASS_DIFF$hydrogen
      hits <- .oad_fragments(scan, tolerance, c(f, f - H2O), 0.1)
      if (hits$count >= 2L && hits$average_intensity < 30) {
        ch2 <- c(f[1L] - .oad_cd2, f[1L] - .oad_cd2)
        if (.oad_fragments(scan, tolerance, ch2, 1)$count > 0L) candidates <- list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double, AverageIntensity = hits$average_intensity))
      }
    } else if (name == "[M+Na]+" && .oad_peak(scan, tolerance, theoretical_mz - header, 10) && .oad_peak(scan, tolerance, gly_c + NA_MASS, 5) && .oad_peak(scan, tolerance, gly_o + NA_MASS, 5)) {
      return(.oad_isotope_result(label, 2L, theoretical_mz, adduct))
    } else if (name == "[M-H]-" && .oad_peak(scan, tolerance, 152.995833871, 0)) {
      f <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
      if (.oad_fragments(scan, tolerance, f, 1)$count == 2L) candidates <- list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double))
    } else return(NULL)
  } else if (label == "PS_d5") {
    c3h8no6p <- .oad_mass(3, 8, 1, 6, 1); gly_c <- .oad_mass_d(6, 7, 5, 1, 6, 1); gly_o <- .oad_mass_d(5, 7, 3, 1, 7, 1)
    if (name == "[M+H]+" && !.oad_peak(scan, tolerance, theoretical_mz - c3h8no6p, 10)) return(NULL)
    if (name == "[M+Na]+" && (!.oad_peak(scan, tolerance, theoretical_mz - c3h8no6p, 3) || !.oad_peak(scan, tolerance, gly_c + NA_MASS, 3) || !.oad_peak(scan, tolerance, gly_o + NA_MASS, 3))) return(NULL)
    if (name == "[M-H]-" && (!.oad_peak(scan, tolerance, theoretical_mz - 87.032029, 10) || .oad_peak(scan, tolerance, theoretical_mz - 63.008491, 30))) return(NULL)
    if (!name %in% c("[M+H]+", "[M+Na]+", "[M-H]-")) return(NULL)
    diagnostic <- if (name == "[M+H]+") theoretical_mz - c3h8no6p else theoretical_mz
    if (name == "[M-H]-") f <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    else if (name == "[M+H]+") f <- diagnostic - c(acyl_chain_mass(sn1_carbon, sn1_double), acyl_chain_mass(sn2_carbon, sn2_double)) + MASS_DIFF$hydrogen
    else { candidates <- list(); for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) { carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1; if (double2 >= 7L) next; f <- theoretical_mz - c(acyl_chain_mass(carbon1, double1), acyl_chain_mass(carbon2, double2)) + MASS_DIFF$hydrogen; if (.oad_fragments(scan, tolerance, c(f, f - MASS_DIFF$oxygen - MASS_DIFF$hydrogen), 0.1)$count >= 2L) candidates <- c(candidates, list(list(Sn1Carbon = carbon1, Sn1DoubleBond = double1, Sn2Carbon = carbon2, Sn2DoubleBond = double2, AverageIntensity = .oad_fragments(scan, tolerance, c(f, f - MASS_DIFF$oxygen - MASS_DIFF$hydrogen), 0)$average_intensity))) }; return(.oad_isotope_result(label, 2L, theoretical_mz, adduct, candidates)) }
    hits <- .oad_fragments(scan, tolerance, c(f, f - H2O), if (name == "[M-H]-") 0 else 0.1)
    if (name == "[M-H]-") hits <- .oad_fragments(scan, tolerance, f, 0.01)
    if (hits$count >= 2L || (name == "[M-H]-" && hits$count == 2L)) candidates <- list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double, AverageIntensity = hits$average_intensity))
  } else if (label == "PI_d5") {
    c6h13o9p <- .oad_mass(6, 13, 0, 9, 1); gly_c <- .oad_mass_d(9, 12, 5, 0, 9, 1); gly_o <- .oad_mass_d(8, 12, 3, 0, 10, 1)
    if (name == "[M-H]-") { if (!.oad_peak(scan, tolerance, 241.01188 + ELECTRON, 0) && !.oad_peak(scan, tolerance, 297.037548 + ELECTRON, 0)) return(NULL); f <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double)); if (.oad_fragments(scan, tolerance, f, 0)$count == 2L) candidates <- list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double))
    } else if (name %in% c("[M+NH4]+", "[M+Na]+")) { if (name == "[M+NH4]+" && !.oad_peak(scan, tolerance, theoretical_mz - c6h13o9p - .oad_mass(0, 3, 1), 10)) return(NULL); if (name == "[M+Na]+" && (!.oad_peak(scan, tolerance, theoretical_mz - .oad_mass(6, 10, 0, 5), 5) || !.oad_peak(scan, tolerance, gly_c + NA_MASS, 5) || !.oad_peak(scan, tolerance, gly_o + NA_MASS, 5))) return(NULL); diagnostic <- if (name == "[M+NH4]+") theoretical_mz - c6h13o9p - .oad_mass(0, 3, 1) else theoretical_mz; for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) { carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1; f <- diagnostic - c(acyl_chain_mass(carbon1, double1), acyl_chain_mass(carbon2, double2)) + MASS_DIFF$hydrogen; hits <- .oad_fragments(scan, tolerance, c(f, f - H2O), 0.1); if (hits$count >= 2L) candidates <- c(candidates, list(list(Sn1Carbon = carbon1, Sn1DoubleBond = double1, Sn2Carbon = carbon2, Sn2DoubleBond = double2, AverageIntensity = hits$average_intensity))) } } else return(NULL)
  }
  .oad_isotope_result(label, 2L, theoretical_mz, adduct, candidates)
}

judge_if_phosphatidylglycerol_d5_oad <- function(...) .oad_phospholipid_d5_judge(..., label = "PG_d5")
judge_if_phosphatidylserine_d5_oad <- function(...) .oad_phospholipid_d5_judge(..., label = "PS_d5")
judge_if_phosphatidylinositol_d5_oad <- function(...) .oad_phospholipid_d5_judge(..., label = "PI_d5")

judge_if_lysopc_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L || total_carbon > 28L) return(NULL)
  name <- adduct$AdductIonName %||% ""; result <- function(suffix = "") .oad_isotope_result("LPC_d5", 1L, theoretical_mz, adduct, list(), suffix)
  if (name == "[M+H]+" && .oad_peak(scan, tolerance, 184.07332, 2) && !.oad_peak(scan, tolerance, theoretical_mz - 59.0735, 10)) {
    vectors <- .oad_spectrum_vectors(scan)
    diagnostic_intensity <- max(vectors$intensity[abs(vectors$mass - 184.07332) < tolerance & vectors$intensity > 2], 0)
    diagnostic2_intensity <- max(vectors$intensity[abs(vectors$mass - 104.106990) < tolerance & vectors$intensity > 2], 0)
    return(result(if (diagnostic2_intensity / diagnostic_intensity > 0.3) "/0:0" else ""))
  }
  if (name == "[M+Na]+" && .oad_peak(scan, tolerance, theoretical_mz - 59.072951, 10)) return(result())
  if (name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) { d <- theoretical_mz - if (name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369; if (.oad_peak(scan, tolerance, d, 1) && .oad_peak(scan, tolerance, .oad_fatty_acid_product_ion(total_carbon, total_double_bond), 1)) return(result()) }
  NULL
}

judge_if_lysope_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L || total_carbon > 28L) return(NULL)
  name <- adduct$AdductIonName %||% ""; c2h8no4p <- .oad_mass(2, 8, 1, 4, 1); result <- function() .oad_isotope_result("LPE_d5", 1L, theoretical_mz, adduct)
  if (name == "[M+H]+" && .oad_peak(scan, tolerance, theoretical_mz - c2h8no4p, 0.5)) { alkyl <- sn_carbon * MASS_DIFF$carbon + MASS_DIFF$hydrogen * (2 * sn_carbon - 2 * sn_double_bond + 1); d <- theoretical_mz - c2h8no4p; if (!.oad_peak(scan, tolerance, d - alkyl + PROTON, 0.5) && !.oad_peak(scan, tolerance, alkyl + 2 * MASS_DIFF$hydrogen + 139.00290, 0.5)) return(result()) }
  if (name == "[M+Na]+" && .oad_peak(scan, tolerance, theoretical_mz - 141.019094, 10)) { alkyl <- sn_carbon * MASS_DIFF$carbon + MASS_DIFF$hydrogen * (2 * sn_carbon - 2 * sn_double_bond + 1); if (!.oad_peak(scan, tolerance, theoretical_mz - 141.019094 - alkyl + PROTON, 10) && !.oad_peak(scan, tolerance, alkyl + 139.00290 + 2 * MASS_DIFF$hydrogen, 10)) return(result()) }
  if (name == "[M-H]-" && .oad_peak(scan, tolerance, theoretical_mz - 159.0691, 10)) return(result())
  NULL
}

judge_if_lysopg_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""; d <- acyl_chain_mass(sn_carbon, sn_double_bond) + (12 * 3 + .oad_deuterium * 5 + MASS_DIFF$oxygen * 2) + PROTON
  if ((name == "[M-H]-" && .oad_peak(scan, tolerance, 152.99583, 1) && .oad_peak(scan, tolerance, .oad_fatty_acid_product_ion(total_carbon, total_double_bond), 1)) || (name %in% c("[M+H]+", "[M+NH4]+") && .oad_peak(scan, tolerance, d, 1))) return(.oad_isotope_result("LPG_d5", 1L, theoretical_mz, adduct)); NULL
}

judge_if_lysopi_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""; d <- acyl_chain_mass(sn_carbon, sn_double_bond) + (12 * 3 + .oad_deuterium * 5 + MASS_DIFF$oxygen * 2) + PROTON
  if ((name == "[M-H]-" && .oad_peak(scan, tolerance, 241.0118806 + ELECTRON, 0) && .oad_peak(scan, tolerance, 320.087237, 0)) || (name %in% c("[M+H]+", "[M+NH4]+") && .oad_peak(scan, tolerance, d, 1))) return(.oad_isotope_result("LPI_d5", 1L, theoretical_mz, adduct)); NULL
}

judge_if_lysops_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""; d <- acyl_chain_mass(sn_carbon, sn_double_bond) + (12 * 3 + .oad_deuterium * 5 + MASS_DIFF$oxygen * 2) + PROTON
  if ((name == "[M-H]-" && .oad_peak(scan, tolerance, 158.0249, 5) && .oad_peak(scan, tolerance, theoretical_mz - 87.032029, 1)) || (name == "[M+H]+" && .oad_peak(scan, tolerance, d, 1))) return(.oad_isotope_result("LPS_d5", 1L, theoretical_mz, adduct)); NULL
}

judge_if_dag_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn1_carbon, sn2_carbon, sn1_double, sn2_double, adduct) {
  if (total_carbon > 52L || is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""
  if (name == "[M+NH4]+") { if (sn2_double >= 7L || !.oad_peak(scan, tolerance, theoretical_mz - 17.026549, 1)) return(NULL); f <- theoretical_mz - 17.026549 - c(acyl_chain_mass(sn1_carbon, sn1_double), acyl_chain_mass(sn2_carbon, sn2_double)) - H2O + MASS_DIFF$hydrogen; if (.oad_fragments(scan, tolerance, f, 0)$count != 2L) return(NULL); return(.oad_isotope_result("DG_d5", 2L, theoretical_mz, adduct)) }
  if (name == "[M+Na]+") { candidates <- list(); for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) { carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1; if (double2 >= 7L) next; f <- theoretical_mz - c(acyl_chain_mass(carbon1, double1), acyl_chain_mass(carbon2, double2)) - H2O + MASS_DIFF$hydrogen; hits <- .oad_fragments(scan, tolerance, f, 0); if (hits$count == 2L) candidates <- c(candidates, list(list(Sn1Carbon = carbon1, Sn1DoubleBond = double1, Sn2Carbon = carbon2, Sn2DoubleBond = double2, AverageIntensity = hits$average_intensity))) }; if (length(candidates) == 0L) return(NULL); return(.oad_isotope_result("DG_d5", 2L, theoretical_mz, adduct, candidates)) }
  NULL
}

judge_if_triacylglycerol_d5_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sn1_carbon, sn2_carbon, sn3_carbon, sn1_double, sn2_double, sn3_double, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""; d <- if (name == "[M+NH4]+") theoretical_mz - 17.026549 else theoretical_mz; f <- d - c(acyl_chain_mass(sn1_carbon, sn1_double), acyl_chain_mass(sn2_carbon, sn2_double), acyl_chain_mass(sn3_carbon, sn3_double)) - H2O + MASS_DIFF$hydrogen
  if (name == "[M+NH4]+" && (sn1_carbon == 18L && sn1_double == 5L || sn2_carbon == 18L && sn2_double == 5L || sn3_carbon == 18L && sn3_double == 5L)) return(NULL)
  if (name == "[M+Na]+") { hits <- .oad_fragments(scan, tolerance, f, 0.1); if (hits$count < 3L) { d <- theoretical_mz - 22.9892207 + MASS_DIFF$hydrogen; f <- d - c(acyl_chain_mass(sn1_carbon, sn1_double), acyl_chain_mass(sn2_carbon, sn2_double), acyl_chain_mass(sn3_carbon, sn3_double)) - H2O + MASS_DIFF$hydrogen; hits <- .oad_fragments(scan, tolerance, f, 0.1) } } else if (name == "[M+NH4]+") hits <- .oad_fragments(scan, tolerance, f, 3) else return(NULL)
  if (hits$count != 3L) return(NULL); .oad_isotope_result("TG_d5", 3L, theoretical_mz, adduct, list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double, Sn3Carbon = sn3_carbon, Sn3DoubleBond = sn3_double, AverageIntensity = hits$average_intensity)))
}

judge_if_sphingomyelin_d9_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sph_carbon, acyl_carbon, sph_double, acyl_double, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""; c5h5d9no4p <- .oad_mass_d(5, 5, 9, 1, 4, 1); c2h2n <- .oad_mass(2, 2, 1)
  if (name == "[M+H]+" && .oad_peak(scan, tolerance, c5h5d9no4p + PROTON, 1) && sph_carbon > 13L && !(sph_carbon == 16L && sph_double >= 3L) && acyl_carbon >= 8L) { f <- c(acyl_chain_mass(acyl_carbon, acyl_double) + c2h2n + MASS_DIFF$hydrogen + PROTON, acyl_chain_mass(acyl_carbon, acyl_double) + c2h2n + MASS_DIFF$hydrogen + PROTON + c5h5d9no4p - MASS_DIFF$hydrogen); hits <- .oad_fragments(scan, tolerance, f, 0); if (hits$count == 2L) return(.oad_isotope_result("SM_d9", 2L, theoretical_mz, adduct, list(list(SphCarbon = sph_carbon, SphDoubleBond = sph_double, AcylCarbon = acyl_carbon, AcylDoubleBond = acyl_double, AverageIntensity = hits$average_intensity)), "d")) }
  if (name == "[M+Na]+" && .oad_peak(scan, tolerance, theoretical_mz - 59.0735, 20) && .oad_peak(scan, tolerance, c5h5d9no4p + NA_MASS, 30)) return(.oad_isotope_result("SM_d9", 2L, theoretical_mz, adduct, list(), "d"))
  if (name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) { d <- theoretical_mz - if (name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 77.060255 else 60.021129369; if (.oad_peak(scan, tolerance, d, 50) && .oad_peak(scan, tolerance, 174.089522, 0.01)) { hits <- .oad_fragments(scan, tolerance, d - acyl_chain_mass(acyl_carbon, acyl_double) + MASS_DIFF$hydrogen, 0); return(.oad_isotope_result("SM_d9", 2L, theoretical_mz, adduct, if (hits$count == 1L) list(list(SphCarbon = sph_carbon, SphDoubleBond = sph_double, AcylCarbon = acyl_carbon, AcylDoubleBond = acyl_double, AverageIntensity = hits$average_intensity)) else list(), "d")) } }
  NULL
}

judge_if_ceramidens_d7_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, sph_carbon, acyl_carbon, sph_double, acyl_double, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); name <- adduct$AdductIonName %||% ""
  if (name %in% c("[M+H]+", "[M+H-H2O]+")) { d <- if (name == "[M+H]+") theoretical_mz - H2O else theoretical_mz; if (!.oad_peak(scan, tolerance, d, 5) || acyl_double >= 7L || !.oad_peak(scan, tolerance, d - acyl_chain_mass(acyl_carbon, acyl_double) + MASS_DIFF$hydrogen - H2O, 5)) return(NULL); f <- c(d - acyl_chain_mass(acyl_carbon, acyl_double) + MASS_DIFF$hydrogen, d - acyl_chain_mass(acyl_carbon, acyl_double) + MASS_DIFF$hydrogen - H2O - 12); hits <- .oad_fragments(scan, tolerance, f, 1); if (hits$count < if (acyl_carbon < 12L) 2L else 1L) return(NULL); return(.oad_isotope_result("Cer_d7", 2L, theoretical_mz, adduct, list(list(SphCarbon = sph_carbon, SphDoubleBond = sph_double, AcylCarbon = acyl_carbon, AcylDoubleBond = acyl_double, AverageIntensity = hits$average_intensity)), "d")) }
  if (name == "[M+Na]+" && !.oad_peak(scan, tolerance, theoretical_mz - 162.052833 - H2O, 1)) { f <- .oad_mass(sph_carbon, 2 * sph_carbon - 2 * sph_double + 1) + MASS_DIFF$hydrogen - MASS_DIFF$oxygen + 7 * MASS_DIFF$hydrogen - H2O + PROTON; if (.oad_peak(scan, tolerance, f, 1)) return(.oad_isotope_result("Cer_d7", 2L, theoretical_mz, adduct, list(), "d")) }
  if (name %in% c("[M-H]-", "[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) { d <- if (name == "[M-H]-") theoretical_mz else theoretical_mz - MASS_DIFF$hydrogen - if (name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 59.013864 else 44.998214; f <- c(d - 12 - H2O, d - 12 - 2 * H2O, ((sph_carbon - 2) * 12) + MASS_DIFF$hydrogen * ((sph_carbon - 2) * 2 - 2 * sph_double) - 1 + MASS_DIFF$oxygen + 7 * MASS_DIFF$hydrogen, .oad_fatty_acid_product_ion(acyl_carbon, acyl_double) - MASS_DIFF$oxygen - 2 * MASS_DIFF$hydrogen); if ((.oad_peak(scan, tolerance, f[1], 0) || .oad_peak(scan, tolerance, f[2], 0)) && !(name %in% c("[M-H]-") && .oad_peak(scan, tolerance, theoretical_mz - MASS_DIFF$hydrogen - 44.998214, 50)) && .oad_fragments(scan, tolerance, f[3:4], 0)$count >= 1L) return(.oad_isotope_result("Cer_d7", 2L, theoretical_mz, adduct, list(), "d")) }
  NULL
}

judge_if_cholesteryl_ester_d7_oad <- function(scan, tolerance, theoretical_mz, total_carbon, total_double_bond, adduct) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL); skeleton <- .oad_mass_d(27, 39, 7, 0, 1); name <- adduct$AdductIonName %||% ""
  if (name == "[M+NH4]+" && total_carbon < 41L && total_double_bond < 4L && .oad_peak(scan, tolerance, skeleton - H2O + PROTON, 1)) return(.oad_isotope_result("CE_d7", 1L, theoretical_mz, adduct))
  if (name == "[M+Na]+" && .oad_peak(scan, tolerance, skeleton - H2O, 10)) return(.oad_isotope_result("CE_d7", 1L, theoretical_mz, adduct))
  NULL
}

characterize_pc_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylcholine_d5_oad", "two_chain") }
characterize_pe_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylethanolamine_d5_oad", "two_chain") }
characterize_pg_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylglycerol_d5_oad", "two_chain") }
characterize_ps_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylserine_d5_oad", "two_chain") }
characterize_pi_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_phosphatidylinositol_d5_oad", "two_chain") }
characterize_lpc_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysopc_d5_oad", "lyso") }
characterize_lpe_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysope_d5_oad", "lyso") }
characterize_lpg_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysopg_d5_oad", "lyso") }
characterize_lpi_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysopi_d5_oad", "lyso") }
characterize_lps_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_lysops_d5_oad", "lyso") }
characterize_dg_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_dag_d5_oad", "two_chain") }
characterize_tg_d5_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_triacylglycerol_d5_oad", "three_chain") }
characterize_sm_d9_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_sphingomyelin_d9_oad", "two_chain") }
characterize_cer_ns_d7_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_ceramidens_d7_oad", "two_chain") }
characterize_ce_d7_oad <- function(scan, reference, tolerance, mz_begin, mz_end) { .oad_call_oad_characterizer(scan, reference, tolerance, "judge_if_cholesteryl_ester_d7_oad", "single_chain") }

# OAD-local JudgeIf translations. These use the C# argument shape and only the
# OAD spectrum primitives above; regular CID judge functions are deliberately
# not used here.
.oad_phospholipid_molecule_level2 <- function(lipid_class, sn1_carbon, sn1_double,
                                              sn2_carbon, sn2_double, score) {
  chains <- rbind(c(sn1_carbon, sn1_double), c(sn2_carbon, sn2_double))
  chains <- chains[order(chains[, 2L], chains[, 1L]), , drop = FALSE]
  sn1_chain <- paste0(chains[1L, 1L], ":", chains[1L, 2L])
  sn2_chain <- paste0(chains[2L, 1L], ":", chains[2L, 2L])
  total_carbon <- sn1_carbon + sn2_carbon
  total_double <- sn1_double + sn2_double
  total_string <- paste0(total_carbon, ":", total_double)
  list(LipidClass = lipid_class, AnnotationLevel = 2L,
       SublevelLipidName = paste(lipid_class, total_string),
       LipidName = paste(lipid_class, paste0(sn1_chain, "_", sn2_chain)),
       TotalCarbonCount = total_carbon, TotalDoubleBondCount = total_double,
       TotalChainString = total_string, Score = score,
       Sn1CarbonCount = chains[1L, 1L], Sn1DoubleBondCount = chains[1L, 2L],
       Sn1AcylChainString = sn1_chain, Sn2CarbonCount = chains[2L, 1L],
       Sn2DoubleBondCount = chains[2L, 2L], Sn2AcylChainString = sn2_chain)
}

.oad_phospholipid_molecule_level3 <- function(lipid_class, sn1_carbon, sn1_double,
                                              sn2_carbon, sn2_double, score) {
  candidate <- .oad_phospholipid_molecule_level2(
    lipid_class, sn1_carbon, sn1_double, sn2_carbon, sn2_double, score)
  candidate$AnnotationLevel <- 3L
  candidate$LipidName <- paste(lipid_class, paste0(sn1_carbon, ":", sn1_double,
                                                   "/", sn2_carbon, ":", sn2_double))
  candidate$Sn1CarbonCount <- sn1_carbon
  candidate$Sn1DoubleBondCount <- sn1_double
  candidate$Sn1AcylChainString <- paste0(sn1_carbon, ":", sn1_double)
  candidate$Sn2CarbonCount <- sn2_carbon
  candidate$Sn2DoubleBondCount <- sn2_double
  candidate$Sn2AcylChainString <- paste0(sn2_carbon, ":", sn2_double)
  candidate
}

.oad_single_chain_molecule_level2 <- function(lipid_class, carbon, double_bond,
                                              score, suffix = "") {
  chain <- paste0(carbon, ":", double_bond, suffix)
  total_string <- paste0(carbon, ":", double_bond)
  list(LipidClass = lipid_class, AnnotationLevel = 2L,
       SublevelLipidName = paste(lipid_class, total_string),
       LipidName = paste(lipid_class, chain), TotalCarbonCount = carbon,
       TotalDoubleBondCount = double_bond, TotalChainString = total_string,
       Score = score, Sn1CarbonCount = carbon, Sn1DoubleBondCount = double_bond,
       Sn1AcylChainString = chain)
}

.oad_direct_result <- function(label, level, theoretical_mz, adduct, candidates = list()) {
  result <- if (length(candidates) > 0L && is.list(candidates[[1L]])) {
    candidates[[which.max(vapply(candidates, function(candidate) {
      candidate$Score %||% 0
    }, numeric(1)))]]
  } else {
    list(LipidName = label, AnnotationLevel = as.integer(level))
  }
  result$LipidName <- result$LipidName %||% label
  result$AnnotationLevel <- as.integer(result$AnnotationLevel %||% level)
  result$Mz <- theoretical_mz
  result$Adduct <- adduct
  result
}

.oad_two_chain_result <- function(scan, tolerance, theoretical_mz, total_carbon,
                                  total_double_bond, sn1_carbon, sn2_carbon,
                                  sn1_double, sn2_double, adduct, label,
                                  head_mz = numeric(), head_threshold = 0,
                                  mode = "positive", diagnostic_mz = NULL,
                                  diagnostic_threshold = 0,
                                  require_fragments = 2L) {
  if (is.null(scan$Spectrum) || length(scan$Spectrum) == 0L) return(NULL)
  if (length(head_mz) > 0L && !all(vapply(head_mz, function(mz) {
    .oad_peak(scan, tolerance, mz, head_threshold)
  }, logical(1)))) return(NULL)
  if (!is.null(diagnostic_mz) && !.oad_peak(scan, tolerance, diagnostic_mz, diagnostic_threshold)) return(NULL)
  fragments <- if (mode == "positive") {
    c(theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen,
      theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen - H2O,
      theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen,
      theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen - H2O)
  } else {
    c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double),
      .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
  }
  hits <- .oad_fragments(scan, tolerance, fragments, 0)
  if (hits$count < require_fragments) return(NULL)
  candidates <- if (label %in% c("PC", "PE", "PG", "PS", "PI", "BMP")) {
    list(.oad_phospholipid_molecule_level2(
      label, sn1_carbon, sn1_double, sn2_carbon, sn2_double,
      hits$average_intensity))
  } else list()
  .oad_direct_result(label, 2L, theoretical_mz, adduct, candidates)
}

judge_if_phosphatidylcholine_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                             total_carbon, total_double_bond, sn1_carbon,
                                             sn2_carbon, sn1_double, sn2_double, adduct) {
  adduct_name <- adduct$AdductIonName %||% "[M+H]+"
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  result <- function(candidates, level = if (length(candidates) > 0L) 2L else 1L) {
    .oad_direct_result("PC", level, theoretical_mz, adduct, candidates)
  }
  count_hits <- function(fragments, threshold) sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, threshold), logical(1)))
  add_candidate <- function(candidates, fragments, threshold, required = 2L,
                            carbon1 = sn1_carbon, double1 = sn1_double,
                            carbon2 = sn2_carbon, double2 = sn2_double,
                            level = 2L) {
    if (count_hits(fragments, threshold) >= required) {
      candidate <- if (level == 3L) {
        .oad_phospholipid_molecule_level3("PC", carbon1, double1, carbon2, double2,
                                          .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)
      } else {
        .oad_phospholipid_molecule_level2("PC", carbon1, double1, carbon2, double2,
                                          .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)
      }
      c(candidates, list(candidate))
    } else candidates
  }
  c5h14no4p <- .oad_mass(5, 14, 1, 4, 1)
  c3h9n <- .oad_mass(3, 9, 1, 0, 0)
  gly_c <- .oad_mass(8, 18, 1, 4, 1)
  gly_o <- .oad_mass(7, 16, 1, 5, 1)
  if (adduct_name == "[M+H]+") {
    if (sn1_carbon < 10L || sn2_carbon < 10L || sn1_double > 6L || sn2_double > 6L) return(NULL)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p + PROTON, 3) ||
        .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c3h9n, 10) ||
        .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 141.019094261, 5)) return(NULL)
    fragments <- c(theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen,
      theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen - H2O,
      theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen,
      theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen - H2O)
    return(result(add_candidate(list(), fragments, 0.1)))
  }
  if (adduct_name == "[M+Na]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c5h14no4p, 3) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c3h9n, 3) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p + NA_MASS, 3) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, gly_c + NA_MASS, 3) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o + NA_MASS, 3)) return(NULL)
    candidates <- list()
    for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
      carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1
      if (double2 >= 7L) next
      fragments <- c(theoretical_mz - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen,
        theoretical_mz - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen - H2O,
        theoretical_mz - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen,
        theoretical_mz - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen - H2O)
      candidates <- add_candidate(candidates, fragments, 0.1, carbon1 = carbon1,
                  double1 = double1, carbon2 = carbon2, double2 = double2)
    }
    return(result(candidates))
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 10)) return(NULL)
    if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-") && .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 60.021129369, 10)) return(NULL)
    return(result(add_candidate(list(), c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double)), 0.1)))
  }
  if (adduct_name == "[M+HCO3]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - PROTON - .oad_mass(4, 10, 1, 3), 10)) return(NULL)
    fragments <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.1), logical(1))) < 1L) return(result(list()))
    flag <- .oad_fatty_acid_product_ion(sn1_carbon, sn1_double) + .oad_mass(5, 9, 0, 1, 1)
    level <- if (.oad_peak(ms_scan_prop, ms2_tolerance, flag, 0.1)) 3L else 2L
    return(result(add_candidate(list(), fragments, 0.1, required = 1L, level = level), level))
  }
  NULL
}

judge_if_phosphatidylethanolamine_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                                  total_carbon, total_double_bond, sn1_carbon,
                                                  sn2_carbon, sn1_double, sn2_double, adduct) {
  adduct_name <- adduct$AdductIonName %||% "[M+H]+"
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  result <- function(candidates, level = if (length(candidates) > 0L) 2L else 1L) {
    .oad_direct_result("PE", level, theoretical_mz, adduct, candidates)
  }
  c2h8no4p <- .oad_mass(2, 8, 1, 4, 1)
  gly_c <- .oad_mass(5, 12, 1, 4, 1)
  gly_o <- .oad_mass(4, 10, 1, 5, 1)
  if (adduct_name == "[M+H]+") {
    if (sn1_carbon < 10L || sn2_carbon < 10L || sn1_double > 6L || sn2_double > 6L) return(NULL)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c2h8no4p, 2.5)) return(NULL)
    fragments <- c(acyl_chain_mass(sn1_carbon, sn1_double) - ELECTRON, acyl_chain_mass(sn2_carbon, sn2_double) - ELECTRON)
    hits <- sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.1), logical(1)))
    return(result(if (hits == 2L) list(.oad_phospholipid_molecule_level2(
      "PE", sn1_carbon, sn1_double, sn2_carbon, sn2_double,
      .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)) else list()))
  }
  if (adduct_name == "[M+Na]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c2h8no4p, 3) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_c + NA_MASS, 3) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o + NA_MASS, 3)) return(NULL)
    candidates <- list()
    for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
      carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1
      if (double2 >= 7L) next
      fragments <- c(theoretical_mz - acyl_chain_mass(carbon1, double1) - 2 * MASS_DIFF$hydrogen,
        theoretical_mz - acyl_chain_mass(carbon1, double1) - 2 * MASS_DIFF$hydrogen - MASS_DIFF$oxygen,
        theoretical_mz - acyl_chain_mass(carbon2, double2) - 2 * MASS_DIFF$hydrogen,
        theoretical_mz - acyl_chain_mass(carbon2, double2) - 2 * MASS_DIFF$hydrogen - MASS_DIFF$oxygen)
      if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.1), logical(1))) >= 2L) candidates <- c(candidates, list(.oad_phospholipid_molecule_level2(
        "PE", carbon1, double1, carbon2, double2,
        .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)))
    }
    if (length(candidates) == 0L) return(NULL)
    return(result(candidates))
  }
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 196.03803, 5) || .oad_peak(ms_scan_prop, ms2_tolerance, 152.995833871, 5)) return(NULL)
    fragments <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 1), logical(1))) != 2L) return(NULL)
    return(.oad_direct_result("PE", 2L, theoretical_mz, adduct,
      list(.oad_phospholipid_molecule_level2("PE", sn1_carbon, sn1_double,
        sn2_carbon, sn2_double, .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity))))
  }
  NULL
}

judge_if_phosphatidylglycerol_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                              total_carbon, total_double_bond, sn1_carbon,
                                              sn2_carbon, sn1_double, sn2_double, adduct) {
  adduct_name <- adduct$AdductIonName %||% "[M+H]+"
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  result <- function(candidates) .oad_direct_result("PG", if (length(candidates) > 0L) 2L else 1L, theoretical_mz, adduct, candidates)
  c3h9o6p <- .oad_mass(3, 9, 0, 6, 1)
  gly_c <- .oad_mass(6, 13, 0, 6, 1)
  gly_o <- .oad_mass(5, 11, 0, 7, 1)
  if (adduct_name == "[M+NH4]+") {
    diagnostic <- theoretical_mz - c3h9o6p - .oad_mass(0, 3, 1, 0)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 10)) return(NULL)
    nl1 <- diagnostic - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen
    nl2 <- diagnostic - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen
    v <- .oad_spectrum_vectors(ms_scan_prop)
    nl1_intensity <- max(v$intensity[abs(v$mass - nl1) < ms2_tolerance], 0)
    nl2_intensity <- max(v$intensity[abs(v$mass - nl2) < ms2_tolerance], 0)
    diagnostic_intensity <- max(v$intensity[abs(v$mass - diagnostic) < ms2_tolerance], 0)
    if (nl1_intensity > diagnostic_intensity && nl2_intensity > diagnostic_intensity) return(NULL)
    fragments <- c(nl1, nl1 - H2O, nl2, nl2 - H2O)
    candidates <- if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 1), logical(1))) >= 2L) list(.oad_phospholipid_molecule_level2(
      "PG", sn1_carbon, sn1_double, sn2_carbon, sn2_double,
      .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)) else list()
    ch2 <- c(nl1 - H2O - .oad_mass(1, 2, 0, 0) + MASS_DIFF$hydrogen,
         nl1 - H2O - .oad_mass(1, 2, 0, 0) + MASS_DIFF$hydrogen)
    if (length(candidates) > 0L && sum(vapply(ch2, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 1), logical(1))) == 0L) candidates <- list()
    return(result(candidates))
  }
  if (adduct_name == "[M+Na]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c3h9o6p, 10) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, gly_c + NA_MASS, 5) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o + NA_MASS, 5)) return(NULL)
    return(result(list()))
  }
  if (adduct_name == "[M-H]-") {
    fragments <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    candidates <- if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.01), logical(1))) == 2L) list(.oad_phospholipid_molecule_level2(
      "PG", sn1_carbon, sn1_double, sn2_carbon, sn2_double,
      .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)) else list()
    return(result(candidates))
  }
  NULL
}

judge_if_bismonoacylglycerophosphate_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                                    total_carbon, total_double_bond, sn1_carbon,
                                                    sn2_carbon, sn1_double, sn2_double, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  if ((adduct$AdductIonName %||% "") != "[M+NH4]+") return(NULL)
  diagnostic <- theoretical_mz - .oad_mass(3, 9, 0, 6, 1) - .oad_mass(0, 3, 1, 0)
  if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 1)) return(NULL)
  fragments <- c(diagnostic - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen,
                 diagnostic - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen)
  v <- .oad_spectrum_vectors(ms_scan_prop)
  if (all(vapply(fragments, function(mz) max(v$intensity[abs(v$mass - diagnostic) < ms2_tolerance], 0) > max(v$intensity[abs(v$mass - mz) < ms2_tolerance], 0), logical(1)))) return(NULL)
  if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 10), logical(1))) != 2L) return(NULL)
  .oad_direct_result("BMP", 2L, theoretical_mz, adduct,
    list(.oad_phospholipid_molecule_level2("BMP", sn1_carbon, sn1_double,
      sn2_carbon, sn2_double, .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)))
}

judge_if_phosphatidylserine_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                            total_carbon, total_double_bond, sn1_carbon,
                                            sn2_carbon, sn1_double, sn2_double, adduct) {
  adduct_name <- adduct$AdductIonName %||% "[M+H]+"
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  result <- function(candidates) .oad_direct_result("PS", if (length(candidates) > 0L) 2L else 1L, theoretical_mz, adduct, candidates)
  c3h8no6p <- .oad_mass(3, 8, 1, 6, 1)
  gly_c <- .oad_mass(6, 12, 1, 6, 1)
  gly_o <- .oad_mass(5, 10, 1, 7, 1)
  if (adduct_name == "[M+H]+" && !.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c3h8no6p, 10)) return(NULL)
  if (adduct_name == "[M+Na]+" && (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c3h8no6p, 3) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_c + NA_MASS, 3) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o + NA_MASS, 3))) return(NULL)
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 87.032029, 3) || .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 63.008491, 30)) return(NULL)
    fragments <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    return(result(if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 1), logical(1))) == 2L) list(.oad_phospholipid_molecule_level2(
      "PS", sn1_carbon, sn1_double, sn2_carbon, sn2_double,
      .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)) else list()))
  }
  if (adduct_name == "[M+H]+") fragments <- c(theoretical_mz - c3h8no6p - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen, theoretical_mz - c3h8no6p - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen - H2O, theoretical_mz - c3h8no6p - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen, theoretical_mz - c3h8no6p - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen - H2O)
  else if (adduct_name == "[M+Na]+") {
    candidates <- list()
    for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
      carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1
      if (double2 >= 7L) next
      f <- c(theoretical_mz - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen, theoretical_mz - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen - MASS_DIFF$oxygen - MASS_DIFF$hydrogen, theoretical_mz - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen, theoretical_mz - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen - MASS_DIFF$oxygen - MASS_DIFF$hydrogen)
      if (sum(vapply(f, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.1), logical(1))) >= 2L) candidates <- c(candidates, list(.oad_phospholipid_molecule_level2(
        "PS", carbon1, double1, carbon2, double2,
        .oad_fragments(ms_scan_prop, ms2_tolerance, f, 0)$average_intensity)))
    }
    return(result(candidates))
  } else return(NULL)
  result(if (sum(vapply(fragments, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 1), logical(1))) >= 2L) list(.oad_phospholipid_molecule_level2(
    "PS", sn1_carbon, sn1_double, sn2_carbon, sn2_double,
    .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$average_intensity)) else list())
}

judge_if_phosphatidylinositol_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                              total_carbon, total_double_bond, sn1_carbon,
                                              sn2_carbon, sn1_double, sn2_double, adduct) {
  adduct_name <- adduct$AdductIonName %||% "[M+H]+"
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  c6h13o9p <- .oad_mass(6, 13, 0, 9, 1)
  gly_c <- .oad_mass(9, 17, 0, 9, 1)
  gly_o <- .oad_mass(8, 15, 0, 10, 1)
  result <- function(candidates, empty_level = 1L) .oad_direct_result("PI", if (length(candidates) > 0L) 2L else empty_level, theoretical_mz, adduct, candidates)
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 241.01188 + ELECTRON, 0.01) && !.oad_peak(ms_scan_prop, ms2_tolerance, 297.037548 + ELECTRON, 0.01)) return(NULL)
    f <- c(.oad_fatty_acid_product_ion(sn1_carbon, sn1_double), .oad_fatty_acid_product_ion(sn2_carbon, sn2_double))
    return(result(if (sum(vapply(f, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.01), logical(1))) == 2L) list(.oad_phospholipid_molecule_level2(
      "PI", sn1_carbon, sn1_double, sn2_carbon, sn2_double,
      .oad_fragments(ms_scan_prop, ms2_tolerance, f, 0)$average_intensity)) else list()))
  }
  if (adduct_name == "[M+NH4]+" && !.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c6h13o9p - .oad_mass(0, 3, 1, 0), 10)) return(NULL)
  if (adduct_name == "[M+Na]+" && (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - .oad_mass(6, 10, 0, 5), 5) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_c + NA_MASS, 5) || !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o + NA_MASS, 5))) return(NULL)
  if (!adduct_name %in% c("[M+NH4]+", "[M+Na]+")) return(NULL)
  diagnostic <- if (adduct_name == "[M+NH4]+") theoretical_mz - c6h13o9p - .oad_mass(0, 3, 1, 0) else theoretical_mz
  candidates <- list()
  for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
    carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1
    f <- c(diagnostic - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen, diagnostic - acyl_chain_mass(carbon1, double1) + MASS_DIFF$hydrogen - H2O, diagnostic - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen, diagnostic - acyl_chain_mass(carbon2, double2) + MASS_DIFF$hydrogen - H2O)
    if (sum(vapply(f, function(mz) .oad_peak(ms_scan_prop, ms2_tolerance, mz, 0.1), logical(1))) >= 2L) candidates <- c(candidates, list(.oad_phospholipid_molecule_level2(
      "PI", carbon1, double1, carbon2, double2,
      .oad_fragments(ms_scan_prop, ms2_tolerance, f, 0)$average_intensity)))
  }
  result(candidates, 1L)
}

judge_if_lysopc_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  sn_carbon <- min(sn_carbon, total_carbon)
  sn_double_bond <- min(sn_double_bond, total_double_bond)
  adduct_name <- adduct$AdductIonName %||% ""
  positive_result <- function(candidates, label = "LPC", suffix = "", level = 1L) {
    if (length(candidates) > 0L && !is.null(candidates[[1L]]$Carbon)) {
      candidates <- lapply(candidates, function(candidate) {
        .oad_single_chain_molecule_level2(label, candidate$Carbon, candidate$DoubleBond,
                                          candidate$Score %||% 0, candidate$Suffix %||% suffix)
      })
      level <- 2L
    }
    result <- .oad_direct_result(label, level, theoretical_mz, adduct, candidates)
    result$Candidates <- candidates
    result$Suffix <- suffix
    result
  }
  if (adduct_name == "[M+H]+" && total_carbon <= 28L) {
    diagnostic <- 184.07332
    diagnostic2 <- 104.106990
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 2) ||
        .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 59.0735, 10)) return(NULL)
    pe_loss <- theoretical_mz - 141.019094261 + PROTON
    if (.oad_peak(ms_scan_prop, ms2_tolerance, pe_loss, 3)) {
      vectors <- .oad_spectrum_vectors(ms_scan_prop)
      pe_intensity <- max(vectors$intensity[abs(vectors$mass - pe_loss) < ms2_tolerance], 0)
      diagnostic_intensity <- max(vectors$intensity[abs(vectors$mass - diagnostic) < ms2_tolerance], 0)
      if (pe_intensity > diagnostic_intensity) return(NULL)
    }
    vectors <- .oad_spectrum_vectors(ms_scan_prop)
    diagnostic_intensity <- max(vectors$intensity[abs(vectors$mass - diagnostic) < ms2_tolerance & vectors$intensity > 2], 0)
    diagnostic2_intensity <- max(vectors$intensity[abs(vectors$mass - diagnostic2) < ms2_tolerance & vectors$intensity > 2], 0)
    suffix <- if (diagnostic2_intensity / diagnostic_intensity > 0.3) "/0:0" else ""
    score <- if (total_carbon < 30L) 1 else 0
    return(positive_result(list(list(Carbon = total_carbon, DoubleBond = total_double_bond,
                     Score = score, Suffix = suffix)), "LPC", suffix, 1L))
  }
  if (adduct_name == "[M+Na]+" && total_carbon <= 28L) {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 59.072951, 3)) return(NULL)
    score <- if (total_carbon < 30L) 1 else 0
    return(positive_result(list(list(Carbon = total_carbon, DoubleBond = total_double_bond,
                     Score = score)), "LPC", "", 1L))
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "M+CH3COO]-") && total_carbon <= 28L) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 3) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, .oad_fatty_acid_product_ion(total_carbon, total_double_bond), 3)) return(NULL)
    return(positive_result(list(), "LPC", "", 1L))
  }
  if (adduct_name == "[M+HCO3]-" && total_carbon <= 28L) {
    diagnostic <- theoretical_mz - (12 + 3 * MASS_DIFF$oxygen + MASS_DIFF$hydrogen) - PROTON - .oad_mass(3, 9, 1)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 5)) return(NULL)
    return(positive_result(list(), "LPC", "", 1L))
  }
  NULL
}

judge_if_etherlysopc_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                     total_double_bond, min_sn_carbon, max_sn_carbon,
                                     min_sn_double_bond, max_sn_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  max_sn_carbon <- min(max_sn_carbon, total_carbon)
  max_sn_double_bond <- min(max_sn_double_bond, total_double_bond)
  adduct_name <- adduct$AdductIonName %||% ""
  result <- function() {
    value <- .oad_direct_result("LPC", 1L, theoretical_mz, adduct)
    value$Suffix <- "e"
    value
  }
  if (adduct_name == "[M+H]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 104.106990, 20)) return(NULL)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 184.07332, 1) &&
        !.oad_peak(ms_scan_prop, ms2_tolerance, 124.99982, 1)) return(NULL)
    return(result())
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 10) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic - 89.08461258, 10)) return(NULL)
    return(result())
  }
  if (adduct_name == "[M+HCO3]-") {
    diagnostic <- theoretical_mz - (12 + 3 * MASS_DIFF$oxygen + MASS_DIFF$hydrogen) - PROTON - .oad_mass(3, 9, 1)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 5)) return(NULL)
    return(result())
  }
  NULL
}

judge_if_lysope_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                total_double_bond, sn1_carbon, sn1_double, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  result <- function(label, level) {
    value <- .oad_direct_result(label, level, theoretical_mz, adduct)
    if (label %in% c("PE", "LPC")) value$Suffix <- "e"
    value
  }
  c2h8no4p <- .oad_mass(2, 8, 1, 4, 1)
  ether_reject <- function(diagnostic, threshold) {
    alkyl <- sn1_carbon * MASS_DIFF$carbon + MASS_DIFF$hydrogen * (2 * sn1_carbon - 2 * sn1_double + 1)
    nl <- diagnostic - alkyl + PROTON
    rearrange <- alkyl + 139.00290 + 2 * MASS_DIFF$hydrogen
    .oad_peak(ms_scan_prop, ms2_tolerance, nl, threshold) ||
      .oad_peak(ms_scan_prop, ms2_tolerance, rearrange, threshold)
  }
  if (adduct_name == "[M+H]+") {
    if (total_carbon > 28L || !.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c2h8no4p, 0.5)) return(NULL)
    if (ether_reject(theoretical_mz - c2h8no4p, 0.5)) return(NULL)
    return(if (total_carbon > 30L) result("PE", 1L) else result("LPE", 1L))
  }
  if (adduct_name == "[M+Na]+") {
    diagnostic <- theoretical_mz - 141.019094
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 3) || ether_reject(diagnostic, 3)) return(NULL)
    return(if (total_carbon > 30L) result("PE", 2L) else result("LPE", 1L))
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 10) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic - 89.08461258, 10)) return(NULL)
    return(.oad_direct_result("LPC", 1L, theoretical_mz, adduct))
  }
  if (adduct_name == "[M+HCO3]-") {
    diagnostic <- theoretical_mz - (12 + 3 * MASS_DIFF$oxygen + MASS_DIFF$hydrogen) - PROTON - .oad_mass(3, 9, 1)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 5)) return(NULL)
    return(.oad_direct_result("LPC", 1L, theoretical_mz, adduct))
  }
  NULL
}

judge_if_lysopg_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  diagnostic <- if (adduct_name == "[M-H]-") 152.99583 else
    acyl_chain_mass(sn_carbon, sn_double_bond) + .oad_mass(3, 5, 0, 2) + PROTON
  threshold <- if (adduct_name == "[M-H]-") 1 else 1
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, threshold) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, .oad_fatty_acid_product_ion(total_carbon, total_double_bond), threshold)) return(NULL)
  } else if (adduct_name %in% c("[M+H]+", "[M+NH4]+")) {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, threshold)) return(NULL)
  } else return(NULL)
  .oad_direct_result("LPG", 1L, theoretical_mz, adduct)
}

judge_if_lysopi_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 241.0118806 + ELECTRON, 1) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, 315.048656, 1) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, .oad_fatty_acid_product_ion(total_carbon, total_double_bond), 10)) return(NULL)
  } else if (adduct_name %in% c("[M+H]+", "[M+NH4]+")) {
    diagnostic <- acyl_chain_mass(sn_carbon, sn_double_bond) + .oad_mass(3, 5, 0, 2) + PROTON
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 1)) return(NULL)
  } else return(NULL)
  .oad_direct_result("LPI", 1L, theoretical_mz, adduct)
}

judge_if_lysops_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                total_double_bond, sn_carbon, sn_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 152.99583, 5) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 87.032029, 1)) return(NULL)
  } else if (adduct_name == "[M+H]+") {
    diagnostic <- acyl_chain_mass(sn_carbon, sn_double_bond) + .oad_mass(3, 5, 0, 2) + PROTON
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 1)) return(NULL)
  } else return(NULL)
  .oad_direct_result("LPS", 1L, theoretical_mz, adduct)
}

judge_if_etherpe_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                 total_double_bond, sn1_carbon, sn2_carbon,
                                 sn1_double, sn2_double, adduct) {
  if (total_carbon <= 28L || is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  c2h8no4p <- .oad_mass(2, 8, 1, 4, 1)
  result <- function(candidates) {
    value <- .oad_direct_result("PE", 2L, theoretical_mz, adduct, candidates)
    value$Candidates <- candidates
    value$Suffix <- "e"
    value
  }
  if (adduct_name == "[M+H]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - c2h8no4p, 1)) return(NULL)
    if (sn1_double >= 5L) return(NULL)
    alkyl <- sn1_carbon * MASS_DIFF$carbon + MASS_DIFF$hydrogen * (2 * sn1_carbon - 2 * sn1_double + 1)
    nl1 <- theoretical_mz - alkyl - MASS_DIFF$oxygen
    nl2 <- theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) - H2O
    rearrange <- alkyl + c2h8no4p
    fragments <- c(nl1, nl2, rearrange)
    if (.oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)$count >= 2L) {
      suffix <- if (.oad_peak(ms_scan_prop, ms2_tolerance, rearrange, 5)) "p" else "e"
      ch2 <- .oad_mass(1, 2)
      nl1_ch2 <- nl1 - ch2
      nl2_ch2 <- nl2 - ch2 + 2 * MASS_DIFF$hydrogen
      if (.oad_fragments(ms_scan_prop, ms2_tolerance, c(nl1_ch2, nl2_ch2), 0)$count > 0L) {
        return(result(list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = if (suffix == "p") sn1_double - 1L else sn1_double,
                                Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double,
                                Suffix = suffix))))
      }
    }
    return(result(list()))
  }
  if (adduct_name == "[M-H]-") {
    if (sn1_carbon >= 24L && sn1_double >= 5L) return(NULL)
    sn2 <- .oad_fatty_acid_product_ion(sn2_carbon, sn2_double)
    nl2 <- theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, c(sn2, nl2), 0)
    if (hits$count != 2L) return(NULL)
    return(result(list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double,
                            Sn2Carbon = sn2_carbon, Sn2DoubleBond = sn2_double,
                            Suffix = "e", AverageIntensity = hits$average_intensity))))
  }
  NULL
}

judge_if_triacylglycerol_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                         total_carbon, total_double_bond, sn1_carbon,
                                         sn2_carbon, sn1_double, sn2_double, adduct) {
  judge_if_dag_oad(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
    total_double_bond, sn1_carbon, sn2_carbon, sn1_double, sn2_double, adduct)
}

.oad_fragment1_greater_than_fragment2 <- function(scan, tolerance, fragment1, fragment2) {
  vectors <- .oad_spectrum_vectors(scan)
  fragment1_intensity <- max(vectors$intensity[abs(vectors$mass - fragment1) < tolerance], 0)
  fragment2_intensity <- max(vectors$intensity[abs(vectors$mass - fragment2) < tolerance], 0)
  fragment1_intensity > fragment2_intensity
}

judge_if_dgts_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                              total_carbon, total_double_bond, min_sn_carbon,
                              max_sn_carbon, min_sn_double_bond, max_sn_double_bond,
                              adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  max_sn_carbon <- min(max_sn_carbon, total_carbon)
  max_sn_double_bond <- min(max_sn_double_bond, total_double_bond)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name != "[M+H]+" || max_sn_carbon < min_sn_carbon ||
      max_sn_double_bond < min_sn_double_bond) return(NULL)

  gly_c <- .oad_mass(10, 19, 1, 3) + PROTON
  gly_o <- .oad_mass(9, 17, 1, 4) + PROTON
  if (!.oad_peak(ms_scan_prop, ms2_tolerance, gly_c, 1) &&
      !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o, 1)) return(NULL)

  candidates <- list()
  for (sn1_carbon in min_sn_carbon:max_sn_carbon) {
    for (sn1_double in min_sn_double_bond:max_sn_double_bond) {
      sn2_carbon <- total_carbon - sn1_carbon
      sn2_double <- total_double_bond - sn1_double
      fragments <- c(
        theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen,
        theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) + MASS_DIFF$hydrogen - H2O,
        theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen,
        theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) + MASS_DIFF$hydrogen - H2O
      )
      hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0.01)
      if (hits$count >= 2L) {
        label <- if (.oad_fragment1_greater_than_fragment2(ms_scan_prop, ms2_tolerance, 130.0862, 144.10191)) "DGTS" else "DGTA"
        candidates <- c(candidates, list(list(
          LipidName = label, Carbon = sn1_carbon + sn2_carbon,
          DoubleBond = sn1_double + sn2_double, Sn1Carbon = sn1_carbon,
          Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
          Sn2DoubleBond = sn2_double, AverageIntensity = hits$average_intensity
        )))
      }
    }
  }
  value <- .oad_direct_result("DGTS", 2L, theoretical_mz, adduct, candidates)
  value$Candidates <- candidates
  value$Suffix <- ""
  value
}

judge_if_ldgts_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                               total_carbon, total_double_bond, min_sn_carbon,
                               max_sn_carbon, min_sn_double_bond, max_sn_double_bond,
                               adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  max_sn_carbon <- min(max_sn_carbon, total_carbon)
  max_sn_double_bond <- min(max_sn_double_bond, total_double_bond)
  adduct_name <- adduct$AdductIonName %||% ""
  gly_c <- .oad_mass(10, 19, 1, 3) + PROTON
  gly_o <- .oad_mass(9, 17, 1, 4) + PROTON
  if (adduct_name == "[M+H]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, gly_c, 1) &&
        !.oad_peak(ms_scan_prop, ms2_tolerance, gly_o, 1)) return(NULL)
    label <- if (.oad_fragment1_greater_than_fragment2(ms_scan_prop, ms2_tolerance, 130.0862, 144.10191)) "LDGTS" else "LDGTA"
    value <- .oad_direct_result(label, 1L, theoretical_mz, adduct)
    value$Candidates <- list()
    value$Suffix <- ""
    return(value)
  }
  if (!adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-") ||
      max_sn_carbon < min_sn_carbon || max_sn_double_bond < min_sn_double_bond) return(NULL)
  diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 60.02167792 else 46.00602785
  if (!.oad_peak(ms_scan_prop, ms2_tolerance,
                 diagnostic - .oad_mass(3, 5), 50)) return(NULL)
  for (sn1_carbon in min_sn_carbon:max_sn_carbon) {
    for (sn1_double in min_sn_double_bond:max_sn_double_bond) {
      sn1 <- .oad_fatty_acid_product_ion(sn1_carbon, sn1_double)
      if (.oad_peak(ms_scan_prop, ms2_tolerance, sn1, 0.01)) {
        value <- .oad_direct_result("LDGTS", 1L, theoretical_mz, adduct)
        value$Candidates <- list()
        value$Suffix <- ""
        return(value)
      }
    }
  }
  NULL
}

# Exact OAD translations for the next C# group.  These definitions intentionally
# remain separate from the regular CID characterizers.
judge_if_sphingomyelin_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                       total_carbon, total_double_bond, sph_carbon,
                                       acyl_carbon, sph_double, acyl_double, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  c5h14no4p <- .oad_mass(5, 14, 1, 4, 1)
  adduct_name <- adduct$AdductIonName %||% ""
  candidate <- function(sph, sph_db, acyl, acyl_db, score) list(
    SphCarbon = sph, SphDoubleBond = sph_db, AcylCarbon = acyl,
    AcylDoubleBond = acyl_db, Suffix = "d", AverageIntensity = score
  )
  result <- function(candidates) {
    value <- .oad_direct_result("SM", 2L, theoretical_mz, adduct, candidates)
    value$Candidates <- candidates
    value$Suffix <- "d"
    value
  }
  if (adduct_name == "[M+H]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p + PROTON, 5) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p + PROTON + .oad_mass(2, 3, 1), 1) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p + PROTON + .oad_mass(3, 3, 1, 1), 0.5) ||
        sph_carbon <= 13L || (sph_carbon == 16L && sph_double >= 3L) || acyl_carbon < 8L) return(NULL)
    chain1 <- acyl_chain_mass(acyl_carbon, acyl_double) + .oad_mass(2, 2, 1) + PROTON
    chain2 <- chain1 + c5h14no4p + PROTON - MASS_DIFF$hydrogen
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, c(chain1, chain2), 0)
    if (hits$count < 1L) return(NULL)
    return(result(list(candidate(sph_carbon, sph_double, acyl_carbon, acyl_double, hits$average_intensity))))
  }
  if (adduct_name == "[M+Na]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 59.0735, 20) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p + NA_MASS, 30)) return(NULL)
    candidates <- list()
    for (sph in 6L:total_carbon) for (sph_db in 0L:total_double_bond) {
      acyl <- total_carbon - sph; acyl_db <- total_double_bond - sph_db
      if (acyl_db >= 7L) next
      chain1 <- acyl_chain_mass(acyl, acyl_db) + .oad_mass(2, 2, 1) + MASS_DIFF$hydrogen + NA_MASS
      chain2 <- chain1 + c5h14no4p - MASS_DIFF$hydrogen
      hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, c(chain1, chain2), 20)
      if (hits$count >= 1L) candidates <- c(candidates, list(candidate(sph, sph_db, acyl, acyl_db, hits$average_intensity)))
    }
    if (length(candidates) == 0L) return(NULL)
    return(result(candidates))
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 50) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, 168.042572 + ELECTRON, 0.01) || acyl_carbon < 8L) return(NULL)
    sph_fragment <- diagnostic - acyl_chain_mass(acyl_carbon, acyl_double) + MASS_DIFF$hydrogen
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, sph_fragment, 0)
    candidates <- if (hits$count == 1L) list(candidate(sph_carbon, sph_double, acyl_carbon, acyl_double, hits$average_intensity)) else list()
    value <- result(candidates)
    if (length(candidates) == 0L) return(value)
    return(value)
  }
  if (adduct_name == "[M+HCO3]-") {
    diagnostic <- theoretical_mz - (12 + 3 * MASS_DIFF$oxygen + MASS_DIFF$hydrogen) - PROTON - .oad_mass(3, 9, 1)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 5) ||
        !.oad_peak(ms_scan_prop, ms2_tolerance, 182.057671, 0.5) || acyl_carbon < 8L) return(NULL)
    sph_fragment <- diagnostic - sphingo_chain_mass(sph_carbon, sph_double) + .oad_mass(2, 4, 1, 1)
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, sph_fragment, 0)
    candidates <- if (hits$count == 1L) list(candidate(sph_carbon, sph_double, acyl_carbon, acyl_double, hits$average_intensity)) else list()
    value <- result(candidates)
    return(value)
  }
  NULL
}

judge_if_acylsm_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon,
                                total_double_bond, sph_carbon, ext_carbon, sph_double,
                                ext_double, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  acyl_carbon <- total_carbon - sph_carbon - ext_carbon
  acyl_double <- total_double_bond - sph_double - ext_double
  result <- function(candidates) {
    value <- .oad_direct_result("SM", 3L, theoretical_mz, adduct, candidates)
    value$Candidates <- candidates; value$Suffix <- "d"; value
  }
  if (adduct_name == "[M+H]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 184.073, 1)) return(NULL)
    loss <- theoretical_mz - .oad_fatty_acid_product_ion(ext_carbon, ext_double) - MASS_DIFF$hydrogen + ELECTRON
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, loss, 0)
    candidates <- if (hits$count == 1L) list(list(SphCarbon = sph_carbon + acyl_carbon, SphDoubleBond = sph_double + acyl_double,
      EsterCarbon = ext_carbon, EsterDoubleBond = ext_double, Suffix = "d", AverageIntensity = hits$average_intensity)) else list()
    return(result(candidates))
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 50)) return(NULL)
    if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-") && .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 60.021129369, 50)) return(NULL)
    class_ion3 <- .oad_peak(ms_scan_prop, ms2_tolerance, 168.0431, 0.01)
    ext_loss <- diagnostic - .oad_fatty_acid_product_ion(ext_carbon, ext_double) - MASS_DIFF$hydrogen
    ext_fa <- .oad_fatty_acid_product_ion(ext_carbon, ext_double)
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, c(ext_loss, ext_fa), 0)
    candidates <- list()
    if (hits$count == 2L) {
      amide <- .oad_fatty_acid_product_ion(acyl_carbon, acyl_double) + MASS_DIFF$nitrogen + MASS_DIFF$hydrogen + ELECTRON + MASS_DIFF$oxygen
      if (.oad_peak(ms_scan_prop, ms2_tolerance, amide, 0.01)) {
        candidates <- list(list(SphCarbon = sph_carbon, SphDoubleBond = sph_double, AcylCarbon = acyl_carbon,
          AcylDoubleBond = acyl_double, EsterCarbon = ext_carbon, EsterDoubleBond = ext_double, Suffix = "d", AverageIntensity = hits$average_intensity))
      } else {
        candidates <- list(list(SphCarbon = sph_carbon + acyl_carbon, SphDoubleBond = sph_double + acyl_double,
          EsterCarbon = ext_carbon, EsterDoubleBond = ext_double, Suffix = "d", AverageIntensity = hits$average_intensity))
      }
    }
    if (length(candidates) == 0L && !class_ion3) return(NULL)
    return(result(candidates))
  }
  NULL
}

judge_if_sphingosine_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                     total_carbon, total_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  if ((adduct$AdductIonName %||% "") != "[M+H]+") return(NULL)
  f <- theoretical_mz - c(H2O, 2 * H2O, 2 * H2O + 12)
  if (!.oad_peak(ms_scan_prop, ms2_tolerance, f[1L], 10) || !.oad_peak(ms_scan_prop, ms2_tolerance, f[2L], 1) || !.oad_peak(ms_scan_prop, ms2_tolerance, f[3L], 1)) return(NULL)
  value <- .oad_direct_result("SPB", 1L, theoretical_mz, adduct)
  value$Oxidized <- 2L
  value
}

judge_if_nacylglyoxfa_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz,
                                      total_carbon, total_double_bond, total_oxidized, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name == "[M+H]+" && .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - .oad_mass(1, 1, 0, 2), 1)) {
    value <- .oad_direct_result("NAGly", 1L, theoretical_mz, adduct); value$Oxidized <- total_oxidized; return(value)
  }
  if (adduct_name == "[M-H]-" && .oad_peak(ms_scan_prop, ms2_tolerance, 74.024752, 10)) {
    value <- .oad_direct_result("NAGly", 1L, theoretical_mz, adduct); value$Oxidized <- total_oxidized; return(value)
  }
  NULL
}

judge_if_dag_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond,
                             sn1_carbon, sn2_carbon, sn1_double, sn2_double, adduct) {
  if (total_carbon > 52L || is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name == "[M+NH4]+") {
    diagnostic <- theoretical_mz - 17.026549
    if (sn2_double <= 7L) return(NULL)
    fragments <- diagnostic - c(acyl_chain_mass(sn1_carbon, sn1_double), acyl_chain_mass(sn2_carbon, sn2_double)) - H2O + MASS_DIFF$hydrogen
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)
    if (hits$count != 2L) return(NULL)
    return(.oad_direct_result("DG", 2L, theoretical_mz, adduct))
  }
  if (adduct_name == "[M+Na]+") {
    if (sn2_double >= 7L) return(NULL)
    candidates <- list()
    for (carbon1 in 6L:total_carbon) for (double1 in 0L:total_double_bond) {
      carbon2 <- total_carbon - carbon1; double2 <- total_double_bond - double1
      if (double2 >= 7L) next
      diagnostic <- theoretical_mz
      fragments <- diagnostic - c(acyl_chain_mass(carbon1, double1), acyl_chain_mass(carbon2, double2)) - H2O + 2 * MASS_DIFF$hydrogen
      hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, fragments, 0)
      if (hits$count == 2L) candidates <- c(candidates, list(list(Sn1Carbon = carbon1, Sn1DoubleBond = double1, Sn2Carbon = carbon2, Sn2DoubleBond = double2, AverageIntensity = hits$average_intensity)))
    }
    if (length(candidates) == 0L) return(NULL)
    value <- .oad_direct_result("DG", 2L, theoretical_mz, adduct, candidates); value$Candidates <- candidates; return(value)
  }
  NULL
}

judge_if_etherpc_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond,
                                 sn1_carbon, sn2_carbon, sn1_double, sn2_double, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  result <- function(candidates = list()) {
    value <- .oad_direct_result("PC", 2L, theoretical_mz, adduct, candidates); value$Candidates <- candidates; value$Suffix <- "e"; value
  }
  c5h14no4p_h <- .oad_mass(5, 14, 1, 4, 1) + PROTON
  if (adduct_name == "[M+H]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, c5h14no4p_h, 1) || .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 59.0735, 10)) return(NULL)
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, c(theoretical_mz - acyl_chain_mass(sn2_carbon, sn2_double) - H2O,
      theoretical_mz - acyl_chain_mass(sn1_carbon, sn1_double) - 2 * MASS_DIFF$hydrogen), 0)
    if (hits$count < 1L) return(result())
    return(result(list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
      Sn2DoubleBond = sn2_double, Suffix = "e", AverageIntensity = hits$average_intensity))))
  }
  if (adduct_name == "[M+Na]+") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 59.072951, 10)) return(NULL)
    return(result())
  }
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 10)) return(NULL)
    if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-") && .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 60.021129369, 30)) return(NULL)
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, .oad_fatty_acid_product_ion(sn2_carbon, sn2_double), 0)
    if (hits$count != 1L) return(NULL)
    return(result(list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
      Sn2DoubleBond = sn2_double, Suffix = "e", AverageIntensity = hits$average_intensity))))
  }
  if (adduct_name == "[M+HCO3]-") {
    diagnostic <- theoretical_mz - PROTON - .oad_mass(4, 10, 1, 3)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 10)) return(NULL)
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, .oad_fatty_acid_product_ion(sn2_carbon, sn2_double), 0)
    if (hits$count != 1L) return(NULL)
    return(result(list(list(Sn1Carbon = sn1_carbon, Sn1DoubleBond = sn1_double, Sn2Carbon = sn2_carbon,
      Sn2DoubleBond = sn2_double, Suffix = "e", AverageIntensity = hits$average_intensity))))
  }
  NULL
}

judge_if_cholesteryl_ester_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  skeleton <- .oad_mass(27, 46, 0, 1)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name == "[M+NH4]+" &&
      .oad_peak(ms_scan_prop, ms2_tolerance, skeleton - H2O + PROTON, 20) &&
      !(total_carbon >= 41L && total_double_bond >= 4L)) return(.oad_direct_result("CE", 1L, theoretical_mz, adduct))
  if (adduct_name == "[M+Na]+" && .oad_peak(ms_scan_prop, ms2_tolerance, skeleton - H2O, 10)) return(.oad_direct_result("CE", 1L, theoretical_mz, adduct))
  NULL
}

judge_if_acylcarnitine_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  if (!((adduct$AdductIonName %||% "") %in% c("[M+H]+", "[M]+"))) return(NULL)
  acyl <- acyl_chain_mass(total_carbon, total_double_bond)
  loss <- theoretical_mz - acyl + PROTON
  if (!.oad_peak(ms_scan_prop, ms2_tolerance, loss, 0.01) && !.oad_peak(ms_scan_prop, ms2_tolerance, loss - H2O, 0.01)) return(NULL)
  .oad_direct_result("CAR", 1L, theoretical_mz, adduct)
}

judge_if_sphingomyelin_phyto_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond,
                                             min_sph_carbon, max_sph_carbon, min_sph_double, max_sph_double, adduct) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  max_sph_carbon <- min(max_sph_carbon, total_carbon); max_sph_double <- min(max_sph_double, total_double_bond)
  adduct_name <- adduct$AdductIonName %||% ""
  if (adduct_name == "[M+H]+" && !.oad_peak(ms_scan_prop, ms2_tolerance, 184.07332, 30)) return(NULL)
  if (adduct_name %in% c("[M+FA-H]-", "[M+Hac-H]-", "[M+HCOO]-", "[M+CH3COO]-")) {
    diagnostic <- theoretical_mz - if (adduct_name %in% c("[M+CH3COO]-", "[M+Hac-H]-")) 74.036779433 else 60.021129369
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic, 50) || !.oad_peak(ms_scan_prop, ms2_tolerance, 168.042572 + ELECTRON, 0.01)) return(NULL)
  } else if (adduct_name != "[M+H]+") return(NULL)
  value <- .oad_direct_result("SM", 2L, theoretical_mz, adduct); value$Suffix <- "t"; value
}

judge_if_shexcer_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond,
                                 sph_carbon, acyl_carbon, sph_double, acyl_double, adduct, total_oxidized) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""
  diagnostic <- NULL
  if (adduct_name %in% c("[M+H]+", "[M+NH4]+")) {
    base <- if (adduct_name == "[M+NH4]+") theoretical_mz - .oad_mass(0, 3, 1) else theoretical_mz
    diagnostic <- c(base - MASS_DIFF$sulfur - 3 * MASS_DIFF$oxygen - H2O - ELECTRON, base - MASS_DIFF$sulfur - 3 * MASS_DIFF$oxygen - H2O - ELECTRON - 162.052833)
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic[1L], 1) && !.oad_peak(ms_scan_prop, ms2_tolerance, diagnostic[2L], 1)) return(NULL)
    acyl_oxidized <- total_oxidized - 2L
    sph1 <- diagnostic[2L] - acyl_chain_mass(acyl_carbon, acyl_double) + MASS_DIFF$hydrogen - acyl_oxidized * MASS_DIFF$oxygen
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, c(sph1, sph1 - H2O, sph1 - H2O - 12), 0)
    candidates <- if (hits$count >= 1L) list(list(SphCarbon = sph_carbon, SphDoubleBond = sph_double, AcylCarbon = acyl_carbon, AcylDoubleBond = acyl_double, Oxidized = acyl_oxidized, Suffix = "d", AverageIntensity = hits$average_intensity)) else list()
    value <- .oad_direct_result("SHexCer", 2L, theoretical_mz, adduct, candidates); value$Candidates <- candidates; value$Suffix <- "d"; value$Oxidized <- acyl_oxidized; return(value)
  }
  if (adduct_name == "[M-H]-" && .oad_peak(ms_scan_prop, ms2_tolerance, 96.960103266, 0.1)) {
    value <- .oad_direct_result("SHexCer", 2L, theoretical_mz, adduct); value$Candidates <- list(); value$Suffix <- "d"; value$Oxidized <- total_oxidized - 2L; return(value)
  }
  NULL
}

judge_if_picermide_oad <- function(ms_scan_prop, ms2_tolerance, theoretical_mz, total_carbon, total_double_bond,
                                   sph_carbon, acyl_carbon, sph_double, acyl_double, adduct, total_oxidized) {
  if (is.null(ms_scan_prop$Spectrum) || length(ms_scan_prop$Spectrum) == 0L) return(NULL)
  adduct_name <- adduct$AdductIonName %||% ""; acyl_oxidized <- total_oxidized - 2L
  if (adduct_name == "[M+H]+" ) {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 260.029722, 1) || .oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - MASS_DIFF$sulfur - 3 * MASS_DIFF$oxygen - H2O - ELECTRON, 1)) return(NULL)
    value <- .oad_direct_result("PI-Cer", 2L, theoretical_mz, adduct); value$Candidates <- list(); value$Suffix <- "d"; value$Oxidized <- acyl_oxidized; return(value)
  }
  if (adduct_name == "[M-H]-") {
    if (!.oad_peak(ms_scan_prop, ms2_tolerance, 241.01188, 5) || !.oad_peak(ms_scan_prop, ms2_tolerance, theoretical_mz - 162.05282, 1)) return(NULL)
    loss <- theoretical_mz - acyl_chain_mass(acyl_carbon, acyl_double) - acyl_oxidized * MASS_DIFF$oxygen + PROTON
    hits <- .oad_fragments(ms_scan_prop, ms2_tolerance, loss, 0)
    candidates <- if (hits$count == 1L) list(list(SphCarbon = sph_carbon, SphDoubleBond = sph_double, AcylCarbon = acyl_carbon, AcylDoubleBond = acyl_double, Oxidized = acyl_oxidized, Suffix = "d", AverageIntensity = hits$average_intensity)) else list()
    value <- .oad_direct_result("PI-Cer", 2L, theoretical_mz, adduct, candidates); value$Candidates <- candidates; value$Suffix <- "d"; value$Oxidized <- acyl_oxidized; return(value)
  }
  NULL
}
