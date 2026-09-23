# =============================================================================
# GetLipidomicsMatchedPeaksScores.R
#
# R translation of the following C# functions from MSDIAL v5.5:
#   MsScanMatching.cs        — IsComparedAvailable, GetMatchedPeaksScores,
#                               GetLipidMoleculeAnnotationResult,
#                               GetLipidomicsMatchedPeaksScores
#   LipidomicsConverter.cs   — ConvertMsdialLipidnameToLipidMoleculeObjectVS2,
#                               acylChainStringSeparatorVS2, setChainPropertiesVS2,
#                               getTotalChainString, GetLipidHeaderString,
#                               setMono/Di/Tri/TetraAcylChainProperty,
#                               SetSingleLipidStructure, SetCoqMolecule,
#                               SetLipidAcylChainProperties,
#                               SetBasicMoleculerProperties,
#                               ConvertMsdialClassDefinitionToSuperClassVS2,
#                               ConvertMsdialClassDefinitionToLbmClassEnumVS2
#
# Functions NOT repeated here (already in LipidMsmsCharacterization.R):
#   - All judge_if_* functions
#   - return_annotation_result()
#   - count_fragment_existence() and related helpers
#   - Mass constants (ELECTRON, PROTON, H2O, etc.)
#
# This script assumes LipidMsmsCharacterization.R has been sourced.
#
# R data structure conventions:
#   LipidMolecule  — named list; fields described in new_lipid_molecule()
#   MoleculeMsReference (mol_ms_ref) — named list with:
#       $Name, $CompoundClass, $PrecursorMz, $SMILES, $InChIKey, $Formula (string),
#       $AdductType (AdductIon list OR plain string "[M+H]+" OR NULL),
#       $adduct_name (string, snake_case fallback from Python-extracted lbm2),
#       $adduct_ion_mode (string, snake_case fallback), $IonMode (string "Positive"/"Negative"),
#       $Spectrum (list of lists, each with $Mass and $Intensity), $Comment
#   IMSScanProperty (ms_scan_prop) — named list with:
#       $Spectrum (same format), $IonMode (string), $PrecursorMz
# =============================================================================


# -----------------------------------------------------------------------------
# 1.  LipidMolecule constructor — equivalent to: new LipidMolecule()
# -----------------------------------------------------------------------------
new_lipid_molecule <- function() {
  list(
    LipidName          = "",
    SublevelLipidName  = "",
    LipidClass         = "Undefined",   # string, e.g. "PC", "Cer-NS"
    AnnotationLevel    = 0L,
    Adduct             = NULL,           # AdductIon list: list(AdductIonName, IonMode, FormatCheck, ...)
    Mz                 = 0.0,
    Rt                 = 0.0,
    Smiles             = "",
    Formula            = "",
    InChIKey           = "",
    IonMode            = "Positive",
    LipidSubclass      = "",
    LipidCategory      = "",
    IsValidatedFormat  = FALSE,
    Score              = 0.0,
    TotalChainString      = "",
    Sn1AcylChainString    = "",
    Sn2AcylChainString    = "",
    Sn3AcylChainString    = "",
    Sn4AcylChainString    = "",
    TotalCarbonCount      = 0L,
    TotalDoubleBondCount  = 0L,
    TotalOxidizedCount    = 0L,
    Sn1CarbonCount        = 0L,
    Sn1DoubleBondCount    = 0L,
    Sn1Oxidizedount       = 0L,   # note: original C# field name has typo "ount"
    Sn2CarbonCount        = 0L,
    Sn2DoubleBondCount    = 0L,
    Sn2Oxidizedount       = 0L,
    Sn3CarbonCount        = 0L,
    Sn3DoubleBondCount    = 0L,
    Sn3Oxidizedount       = 0L,
    Sn4CarbonCount        = 0L,
    Sn4DoubleBondCount    = 0L,
    Sn4Oxidizedount       = 0L
  )
}


# -----------------------------------------------------------------------------
# 2.  Helper: extract named regex capture groups (base-R, perl=TRUE)
# -----------------------------------------------------------------------------
.named_captures <- function(pattern, text) {
  m <- regexpr(pattern, text, perl = TRUE)
  if (m[1] == -1L) return(NULL)
  starts  <- attr(m, "capture.start")[1, ]
  lengths <- attr(m, "capture.length")[1, ]
  nms     <- attr(m, "capture.names")
  result  <- vector("list", length(nms))
  names(result) <- nms
  for (i in seq_along(nms)) {
    if (nchar(nms[i]) > 0 && starts[i] > 0)
      result[[nms[i]]] <- substr(text, starts[i], starts[i] + lengths[i] - 1L)
    else
      result[[nms[i]]] <- NA_character_
  }
  result
}


# -----------------------------------------------------------------------------
# 3.  setChainPropertiesVS2
#     Parses a chain string (e.g. "18:1;2O", "O-16:0", "P-18:1", "d18:1")
#     Returns a named list with carbon, double_bond, oxidized.
# -----------------------------------------------------------------------------
set_chain_properties_vs2 <- function(chain_string) {
  carbon     <- 0L
  double_bond <- 0L
  oxidized   <- 0L

  is_plasmenyl <- grepl("P-", chain_string, fixed = TRUE)

  # Strip prefixes / letter modifiers
  s <- chain_string
  s <- gsub("O-", "",  s, fixed = TRUE)
  s <- gsub("P-", "",  s, fixed = TRUE)
  s <- gsub("N-", "",  s, fixed = TRUE)
  s <- gsub("n-", "",  s, fixed = TRUE)
  s <- gsub("e",  "",  s, fixed = TRUE)
  s <- gsub("p",  "",  s, fixed = TRUE)
  s <- gsub("m",  "",  s, fixed = TRUE)
  s <- gsub("d",  "",  s, fixed = TRUE)
  s <- gsub("t",  "",  s, fixed = TRUE)

  # Remove double-bond position annotations, e.g. (9E,11Z) or (9Z)
  s <- gsub("\\([0-9]+(?:[EZ],)*[0-9]*[EZ]*\\)", "", s, perl = TRUE)

  # Parse oxidized moiety
  if (grepl(";", s, fixed = TRUE)) {
    # e.g. 18:2;2O  or  18:2;(2OH)  or  18:2;O
    parts <- strsplit(s, ";", fixed = TRUE)[[1]]
    s     <- parts[1]
    ox    <- gsub("(", "", parts[2], fixed = TRUE)
    ox    <- gsub(")", "", ox, fixed = TRUE)
    if (grepl("OH", ox, fixed = TRUE)) {
      oxidized <- length(strsplit(ox, ",", fixed = TRUE)[[1]])
    } else if (grepl("O", ox, fixed = TRUE)) {
      ox_num <- suppressWarnings(as.integer(gsub("O", "", ox, fixed = TRUE)))
      oxidized <- if (is.na(ox_num) || ox_num == 0L) 1L else ox_num
    }
  } else if (grepl("+", s, fixed = TRUE)) {
    # e.g. 20:3+3O
    parts    <- strsplit(s, "+", fixed = TRUE)[[1]]
    s        <- parts[1]
    ox_str   <- gsub("O", "", parts[2], fixed = TRUE)
    if (nchar(ox_str) == 0) ox_str <- "1"
    oxidized <- suppressWarnings(as.integer(ox_str))
    if (is.na(oxidized)) oxidized <- 0L
  } else if (grepl("(", s, fixed = TRUE)) {
    # e.g. 18:1(1OH,3OH)  or  18:2(2OH)
    parts    <- strsplit(s, "(", fixed = TRUE)[[1]]
    s        <- parts[1]
    ox_moiety <- gsub(")", "", parts[2], fixed = TRUE)
    oxidized <- length(strsplit(ox_moiety, ",", fixed = TRUE)[[1]])
  }

  # Strip SN position markers
  s <- gsub("-SN1", "", s, fixed = TRUE)
  s <- gsub("-SN2", "", s, fixed = TRUE)

  # Now s should be "carbon:doublebond"
  colon_parts <- strsplit(s, ":", fixed = TRUE)[[1]]
  if (length(colon_parts) >= 2) {
    carbon      <- suppressWarnings(as.integer(colon_parts[1]))
    double_bond <- suppressWarnings(as.integer(colon_parts[2]))
    if (is.na(carbon))      carbon      <- 0L
    if (is.na(double_bond)) double_bond <- 0L
  }

  if (is_plasmenyl) double_bond <- double_bond + 1L

  list(carbon = carbon, double_bond = double_bond, oxidized = oxidized)
}


# -----------------------------------------------------------------------------
# 4.  retrieve_sterol_header_chain_strings
#     For sterol lipids (SE, ST, SG, BA, ASG), extracts the chain string.
#     Returns list(header_string="", chain_string=...).
# -----------------------------------------------------------------------------
retrieve_sterol_header_chain_strings <- function(molecule_string) {
  header_string <- ""
  chain_string  <- ""
  if (grepl("/", molecule_string, fixed = TRUE)) {
    parts        <- strsplit(molecule_string, "/", fixed = TRUE)[[1]]
    chain_string <- trimws(parts[length(parts)])
  } else {
    parts        <- strsplit(molecule_string, " ", fixed = TRUE)[[1]]
    chain_string <- trimws(parts[length(parts)])
  }
  list(header_string = header_string, chain_string = chain_string)
}


# -----------------------------------------------------------------------------
# 5.  acyl_chain_string_separator_vs2
#     Parses a full lipid name string and returns a character vector of
#     individual acyl chain strings.  Returns NULL if parsing fails.
# -----------------------------------------------------------------------------
acyl_chain_string_separator_vs2 <- function(molecule_string) {
  if (length(strsplit(molecule_string, " ", fixed = TRUE)[[1]]) == 1L)
    return(NULL)

  header_string <- trimws(strsplit(molecule_string, " ", fixed = TRUE)[[1]][1])
  chain_string  <- ""

  if (header_string %in% c("SE", "ST", "SG", "BA", "ASG")) {
    res          <- retrieve_sterol_header_chain_strings(molecule_string)
    header_string <- res$header_string
    chain_string  <- res$chain_string
  } else {
    chain_string <- substring(molecule_string, nchar(header_string) + 2L)
  }

  # d-substituted compound: pipe character → unsupported
  if (grepl("|", chain_string, fixed = TRUE)) return(NULL)

  # Remove deuterium label annotations, e.g. (d9)
  chain_string <- gsub("\\(d[0-9]*\\)", "", chain_string, perl = TRUE)

  # Named-group regex patterns (PCRE)
  pattern4   <- "(?P<chain1>.+?)/(?P<chain2>.+?)\\(FA (?P<chain3>.+?)\\)"
  pattern3   <- "(?P<chain1>.+?)\\(FA (?P<chain2>.+?)\\)"
  pattern2   <- "\\((?P<chain1>.+?)\\)(?P<chain2>.+?)[/_](?P<chain3>.+?)$"
  pattern2_1 <- "\\((?P<chain1>.+?)\\)(?P<chain2>.+?)$"
  pattern12  <- "\\(FA (?P<chain2>.+?)\\)(?P<chain1>.+?)$"

  chains <- NULL

  if (grepl("/", chain_string, fixed = TRUE) && grepl("(FA", chain_string, fixed = TRUE)) {
    # pattern 4: e.g. "Cer 14:0;2O/12:0;(3OH)(FA 12:0)"
    g <- .named_captures(pattern4, chain_string)
    if (!is.null(g))
      chains <- c(g[["chain1"]], g[["chain2"]], g[["chain3"]])

  } else if (grepl("(FA", chain_string, fixed = TRUE)) {
    # pattern 3 or 12: e.g. "SM 30:1;2O(FA 14:0)"
    g <- .named_captures(pattern3, chain_string)
    if (!is.null(g)) {
      chain1strings <- g[["chain1"]]
      if (!is.na(chain1strings) && grepl("_", chain1strings, fixed = TRUE)) {
        # pattern 11 variant: chain1 is itself multi-chain, e.g. "TG 16:0_16:1_18:0;O(FA 16:0)"
        base_chains <- strsplit(chain1strings, "_", fixed = TRUE)[[1]]
        chains      <- c(base_chains, g[["chain2"]])
      } else if (!is.na(chain1strings) && nchar(chain1strings) == 0) {
        # pattern 12: e.g. "LPE-N (FA 16:0)18:1"
        g12 <- .named_captures(pattern12, chain_string)
        if (!is.null(g12))
          chains <- c(g12[["chain1"]], g12[["chain2"]])
      } else {
        chains <- c(g[["chain1"]], g[["chain2"]])
      }
    }

  } else if (grepl("(O-", chain_string, fixed = TRUE)) {
    # pattern 2 / 2_1: e.g. "AHexCer (O-14:0)16:1;2O/14:0;O"
    if (grepl("/", chain_string, fixed = TRUE) || grepl("_", chain_string, fixed = TRUE)) {
      g <- .named_captures(pattern2, chain_string)
      if (!is.null(g))
        chains <- c(g[["chain1"]], g[["chain2"]], g[["chain3"]])
    } else {
      g <- .named_captures(pattern2_1, chain_string)
      if (!is.null(g))
        chains <- c(g[["chain1"]], g[["chain2"]])
    }

  } else {
    # Default: split by _ (or / converted to _)
    s       <- gsub("/", "_", chain_string, fixed = TRUE)
    parts   <- strsplit(s, "_", fixed = TRUE)[[1]]
    chains  <- parts[nchar(parts) > 0][seq_len(min(4L, sum(nchar(parts) > 0)))]
  }

  # Remove any NA/empty strings
  chains <- chains[!is.na(chains) & nchar(chains) > 0]
  if (length(chains) == 0L) return(NULL)
  chains
}


# -----------------------------------------------------------------------------
# 6.  get_lipid_header_string
#     Returns the class header from a lipid name string.
#     For sterols (SE/ST/SG/BA/ASG) returns "" (matching C# behaviour).
# -----------------------------------------------------------------------------
get_lipid_header_string <- function(lipidname) {
  lipid_header <- strsplit(lipidname, " ", fixed = TRUE)[[1]][1]
  if (lipid_header %in% c("SE", "ST", "SG", "BA", "ASG")) {
    res          <- retrieve_sterol_header_chain_strings(lipidname)
    lipid_header <- res$header_string   # "" for sterols
  }
  lipid_header
}


# -----------------------------------------------------------------------------
# 7.  get_total_chain_string
#     Computes the total acyl chain string (e.g. "36:2;O").
#     Mirrors getTotalChainString() in LipidomicsConverter.cs.
# -----------------------------------------------------------------------------
get_total_chain_string <- function(carbon, rdb, oxidized, lipid_class,
                                   chain_prefix = "", acyl_chain_count = 2L) {
  rdb_val <- rdb

  # Classes that add +1 to rdb and oxidized (due to ester bond)
  ester_adjusted_classes <- c(
    "Cer-EODS", "Cer-EBDS", "ASM",
    "FAHFA", "NAGly", "NAGlySer", "NAOrn",
    "TG_EST", "DMEDFAHFA", "NATryA"
  )
  if (lipid_class %in% ester_adjusted_classes) {
    rdb_val  <- rdb + 1L
    oxidized <- oxidized + 1L
  }

  # Cer_EOS / HexCer_EOS: only adjust when not exactly 2 chains
  if (lipid_class %in% c("Cer-EOS", "HexCer-EOS") && acyl_chain_count != 2L) {
    rdb_val  <- rdb + 1L
    oxidized <- oxidized + 1L
  }

  # Plasmenyl: subtract 1 from rdb (P- prefix)
  if (chain_prefix == "P-") {
    rdb_val <- rdb - 1L
  }

  ox_string <- if (oxidized == 0L) {
    ""
  } else if (oxidized == 1L) {
    ";O"
  } else {
    paste0(";O", oxidized)
  }

  paste0(chain_prefix, carbon, ":", rdb_val, ox_string)
}


# -----------------------------------------------------------------------------
# 8.  set_basic_moleculer_properties
#     Copies Mz, Smiles, InChIKey, Formula, Adduct, IonMode from query.
#     Modifies molecule in place (returned); mirrors SetBasicMoleculerProperties.
# -----------------------------------------------------------------------------
set_basic_moleculer_properties <- function(molecule, query) {
  molecule$Mz      <- query$PrecursorMz
  molecule$Smiles  <- if (!is.null(query$SMILES))    query$SMILES    else ""
  molecule$InChIKey <- if (!is.null(query$InChIKey)) query$InChIKey  else ""
  # Formula may be a string or a list with $FormulaString
  molecule$Formula <- if (is.character(query$Formula)) {
    query$Formula
  } else if (is.list(query$Formula) && !is.null(query$Formula$FormulaString)) {
    query$Formula$FormulaString
  } else ""

  # Build an AdductIon list (with $AdductIonName and $IonMode) from query.
  # query$AdductType may be:
  #   (a) already an AdductIon list  → use as-is
  #   (b) a plain adduct string      → wrap in list
  #   (c) NULL                       → try flat snake_case fields (adduct_name /
  #                                    adduct_ion_mode) from Python-extracted lbm2
  adduct_val <- query$AdductType
  if (is.null(adduct_val)) {
    adduct_name_str <- query$adduct_name
    if (!is.null(adduct_name_str) && nzchar(adduct_name_str)) {
      ion_mode_str <- if (!is.null(query$adduct_ion_mode)) query$adduct_ion_mode
                      else if (!is.null(query$IonMode))    query$IonMode
                      else if (!is.null(query$ion_mode))   query$ion_mode
                      else NA_character_
      adduct_val <- list(AdductIonName = adduct_name_str,
                         IonMode       = ion_mode_str,
                         FormatCheck   = TRUE)
    }
  } else if (is.character(adduct_val)) {
    ion_mode_str <- if (!is.null(query$IonMode))    query$IonMode
                    else if (!is.null(query$ion_mode)) query$ion_mode
                    else NA_character_
    adduct_val <- list(AdductIonName = adduct_val,
                       IonMode       = ion_mode_str,
                       FormatCheck   = TRUE)
  }
  # else: already an AdductIon list — use as-is

  molecule$Adduct  <- adduct_val
  molecule$IonMode <- query$IonMode
  molecule
}


# -----------------------------------------------------------------------------
# 9.  set_single_lipid_structure
#     For sterols / vitamins / bile acids with no acyl chain detail.
# -----------------------------------------------------------------------------
set_single_lipid_structure <- function(molecule, query) {
  molecule               <- set_basic_moleculer_properties(molecule, query)
  lipidinfo              <- query$Name
  molecule$LipidName        <- lipidinfo
  molecule$SublevelLipidName <- lipidinfo
  molecule$TotalChainString  <- "0:0"
  molecule$TotalCarbonCount  <- 0L
  molecule$TotalDoubleBondCount <- 0L
  molecule$Sn1AcylChainString   <- "0:0"
  molecule$Sn1CarbonCount       <- 0L
  molecule$Sn1DoubleBondCount   <- 0L
  molecule
}


# -----------------------------------------------------------------------------
# 10. set_coq_molecule
#     For CoQ (Coenzyme Q) molecules — CoQ10, CoQ9, etc.
# -----------------------------------------------------------------------------
set_coq_molecule <- function(molecule, query) {
  molecule              <- set_basic_moleculer_properties(molecule, query)
  lipidinfo             <- query$Name
  carbon_count_str      <- substring(lipidinfo, 4)   # "CoQ10" -> "10"
  carbon_count          <- suppressWarnings(as.integer(carbon_count_str))
  if (is.na(carbon_count)) carbon_count <- 0L
  molecule$LipidName           <- lipidinfo
  molecule$SublevelLipidName   <- lipidinfo
  molecule$TotalChainString    <- carbon_count_str
  molecule$TotalCarbonCount    <- carbon_count
  molecule$TotalDoubleBondCount <- 0L
  molecule$Sn1AcylChainString  <- carbon_count_str
  molecule$Sn1CarbonCount      <- carbon_count
  molecule$Sn1DoubleBondCount  <- 0L
  molecule
}


# -----------------------------------------------------------------------------
# 11. set_mono_acyl_chain_property
# -----------------------------------------------------------------------------
set_mono_acyl_chain_property <- function(molecule, lipidname, ontology, chain_strings) {
  if (length(chain_strings) != 1L) return(molecule)
  cp   <- set_chain_properties_vs2(chain_strings[1])
  tot_c  <- cp$carbon
  tot_db <- cp$double_bond
  tot_ox <- cp$oxidized
  total_chain <- chain_strings[1]   # use the original string as-is
  molecule$LipidName           <- lipidname
  molecule$SublevelLipidName   <- lipidname
  molecule$TotalChainString    <- total_chain
  molecule$TotalCarbonCount    <- tot_c
  molecule$TotalDoubleBondCount <- tot_db
  molecule$TotalOxidizedCount  <- tot_ox
  molecule$Sn1AcylChainString  <- chain_strings[1]
  molecule$Sn1CarbonCount      <- tot_c
  molecule$Sn1DoubleBondCount  <- tot_db
  molecule$Sn1Oxidizedount     <- tot_ox
  molecule
}


# -----------------------------------------------------------------------------
# 12. set_di_acyl_chain_property
# -----------------------------------------------------------------------------
set_di_acyl_chain_property <- function(molecule, lipidname, ontology, chain_strings) {
  if (length(chain_strings) != 2L) return(molecule)
  lbm_class   <- convert_msdial_class_to_lbm_class_enum_vs2(ontology)
  lipid_header <- get_lipid_header_string(lipidname)
  cp1 <- set_chain_properties_vs2(chain_strings[1])
  cp2 <- set_chain_properties_vs2(chain_strings[2])
  tot_c  <- cp1$carbon  + cp2$carbon
  tot_db <- cp1$double_bond + cp2$double_bond
  tot_ox <- cp1$oxidized + cp2$oxidized
  chain_prefix <- if (grepl("^O-", chain_strings[1])) "O-" else
                  if (grepl("^P-", chain_strings[1])) "P-" else ""
  total_chain    <- get_total_chain_string(tot_c, tot_db, tot_ox, lbm_class, chain_prefix, 2L)
  sublevel_name  <- paste0(lipid_header, " ", total_chain)
  molecule$SublevelLipidName   <- sublevel_name
  molecule$LipidName           <- lipidname
  molecule$TotalChainString    <- total_chain
  molecule$TotalCarbonCount    <- tot_c
  molecule$TotalDoubleBondCount <- tot_db
  molecule$TotalOxidizedCount  <- tot_ox
  molecule$Sn1AcylChainString  <- chain_strings[1]
  molecule$Sn1CarbonCount      <- cp1$carbon
  molecule$Sn1DoubleBondCount  <- cp1$double_bond
  molecule$Sn1Oxidizedount     <- cp1$oxidized
  molecule$Sn2AcylChainString  <- chain_strings[2]
  molecule$Sn2CarbonCount      <- cp2$carbon
  molecule$Sn2DoubleBondCount  <- cp2$double_bond
  molecule$Sn2Oxidizedount     <- cp2$oxidized
  molecule
}


# -----------------------------------------------------------------------------
# 13. set_tri_acyl_chain_property
# -----------------------------------------------------------------------------
set_tri_acyl_chain_property <- function(molecule, lipidname, ontology, chain_strings) {
  if (length(chain_strings) != 3L) return(molecule)
  lbm_class    <- convert_msdial_class_to_lbm_class_enum_vs2(ontology)
  lipid_header <- get_lipid_header_string(lipidname)
  cp1 <- set_chain_properties_vs2(chain_strings[1])
  cp2 <- set_chain_properties_vs2(chain_strings[2])
  cp3 <- set_chain_properties_vs2(chain_strings[3])
  tot_c  <- cp1$carbon  + cp2$carbon  + cp3$carbon
  tot_db <- cp1$double_bond + cp2$double_bond + cp3$double_bond
  tot_ox <- cp1$oxidized + cp2$oxidized + cp3$oxidized
  chain_prefix <- if (grepl("Ether", lbm_class, fixed = TRUE)) {
    if (grepl("^O-", chain_strings[1])) "O-" else
    if (grepl("^P-", chain_strings[1])) "P-" else ""
  } else ""
  total_chain   <- get_total_chain_string(tot_c, tot_db, tot_ox, lbm_class, chain_prefix, 3L)
  sublevel_name <- paste0(lipid_header, " ", total_chain)
  molecule$SublevelLipidName   <- sublevel_name
  molecule$LipidName           <- lipidname
  molecule$TotalChainString    <- total_chain
  molecule$TotalCarbonCount    <- tot_c
  molecule$TotalDoubleBondCount <- tot_db
  molecule$TotalOxidizedCount  <- tot_ox
  molecule$Sn1AcylChainString  <- chain_strings[1]
  molecule$Sn1CarbonCount      <- cp1$carbon
  molecule$Sn1DoubleBondCount  <- cp1$double_bond
  molecule$Sn1Oxidizedount     <- cp1$oxidized
  molecule$Sn2AcylChainString  <- chain_strings[2]
  molecule$Sn2CarbonCount      <- cp2$carbon
  molecule$Sn2DoubleBondCount  <- cp2$double_bond
  molecule$Sn2Oxidizedount     <- cp2$oxidized
  molecule$Sn3AcylChainString  <- chain_strings[3]
  molecule$Sn3CarbonCount      <- cp3$carbon
  molecule$Sn3DoubleBondCount  <- cp3$double_bond
  molecule$Sn3Oxidizedount     <- cp3$oxidized
  molecule
}


# -----------------------------------------------------------------------------
# 14. set_tetra_acyl_chain_property
# -----------------------------------------------------------------------------
set_tetra_acyl_chain_property <- function(molecule, lipidname, ontology, chain_strings) {
  if (length(chain_strings) != 4L) return(molecule)
  lbm_class    <- convert_msdial_class_to_lbm_class_enum_vs2(ontology)
  lipid_header <- get_lipid_header_string(lipidname)
  cp1 <- set_chain_properties_vs2(chain_strings[1])
  cp2 <- set_chain_properties_vs2(chain_strings[2])
  cp3 <- set_chain_properties_vs2(chain_strings[3])
  cp4 <- set_chain_properties_vs2(chain_strings[4])
  tot_c  <- cp1$carbon  + cp2$carbon  + cp3$carbon  + cp4$carbon
  tot_db <- cp1$double_bond + cp2$double_bond + cp3$double_bond + cp4$double_bond
  tot_ox <- cp1$oxidized + cp2$oxidized + cp3$oxidized + cp4$oxidized
  chain_prefix <- if (grepl("^O-", chain_strings[1])) "O-" else
                  if (grepl("^P-", chain_strings[1])) "P-" else ""
  total_chain   <- get_total_chain_string(tot_c, tot_db, tot_ox, lbm_class, chain_prefix, 4L)
  sublevel_name <- paste0(lipid_header, " ", total_chain)
  molecule$SublevelLipidName   <- sublevel_name
  molecule$LipidName           <- lipidname
  molecule$TotalChainString    <- total_chain
  molecule$TotalCarbonCount    <- tot_c
  molecule$TotalDoubleBondCount <- tot_db
  molecule$TotalOxidizedCount  <- tot_ox
  molecule$Sn1AcylChainString  <- chain_strings[1]
  molecule$Sn1CarbonCount      <- cp1$carbon
  molecule$Sn1DoubleBondCount  <- cp1$double_bond
  molecule$Sn1Oxidizedount     <- cp1$oxidized
  molecule$Sn2AcylChainString  <- chain_strings[2]
  molecule$Sn2CarbonCount      <- cp2$carbon
  molecule$Sn2DoubleBondCount  <- cp2$double_bond
  molecule$Sn2Oxidizedount     <- cp2$oxidized
  molecule$Sn3AcylChainString  <- chain_strings[3]
  molecule$Sn3CarbonCount      <- cp3$carbon
  molecule$Sn3DoubleBondCount  <- cp3$double_bond
  molecule$Sn3Oxidizedount     <- cp3$oxidized
  molecule$Sn4AcylChainString  <- chain_strings[4]
  molecule$Sn4CarbonCount      <- cp4$carbon
  molecule$Sn4DoubleBondCount  <- cp4$double_bond
  molecule$Sn4Oxidizedount     <- cp4$oxidized
  molecule
}


# -----------------------------------------------------------------------------
# 15. set_lipid_acyl_chain_properties
#     Dispatcher: parses chain strings and routes to mono/di/tri/tetra setter.
# -----------------------------------------------------------------------------
set_lipid_acyl_chain_properties <- function(molecule, query) {
  lipidname    <- trimws(query$Name)
  chain_strings <- acyl_chain_string_separator_vs2(lipidname)
  ontology     <- query$CompoundClass
  molecule     <- set_basic_moleculer_properties(molecule, query)
  if (is.null(chain_strings)) return(molecule)
  n <- length(chain_strings)
  if (n == 1L) molecule <- set_mono_acyl_chain_property(molecule, lipidname, ontology, chain_strings)
  else if (n == 2L) molecule <- set_di_acyl_chain_property(molecule, lipidname, ontology, chain_strings)
  else if (n == 3L) molecule <- set_tri_acyl_chain_property(molecule, lipidname, ontology, chain_strings)
  else if (n == 4L) molecule <- set_tetra_acyl_chain_property(molecule, lipidname, ontology, chain_strings)
  molecule
}


# -----------------------------------------------------------------------------
# 16. convert_msdial_class_to_super_class_vs2
#     Maps lipid class string to super-category string.
#     Mirrors ConvertMsdialClassDefinitionToSuperClassVS2.
# -----------------------------------------------------------------------------
convert_msdial_class_to_super_class_vs2 <- function(lipidclass) {
  fatty_acyls <- c(
    "NAE","NAGly","NAGlySer","NAOrn","NAPhe","NATau","NATryA","NA5HT",
    "NAAla","NAGln","NALeu","NAVal","NASer","NAAnt","NAGABA","WE",
    "CAR","FA","OxFA","FAHFA","DMEDFAHFA","DMEDFA","DMEDOxFA"
  )
  glycerolipids <- c(
    "MG","DG","TG","OxTG","TG_EST","EtherDG","EtherTG",
    "LDGTS","LDGTA","LDGCC",
    "DGDG","MGDG","SQDG","DGTS","DGTA","DGCC","DGGA","ADGGA",
    "EtherMGDG","EtherDGDG","EtherSDGDG","EtherSMGDG","SMGDG",
    "MGMG","DGMG"
  )
  glycerophospholipids <- c(
    "LPC","LPA","LPE","LPG","LPI","LPS","BisMeLPA",
    "PC","PA","PE","PG","PI","PS","PT","BMP","HBMP","CL","DLCL","MLCL",
    "EtherPC","EtherPE","EtherPE_O","EtherPE_P","EtherPS","EtherPG","EtherPI",
    "EtherLPC","EtherLPE","EtherLPS","EtherLPG","EtherLPI",
    "OxPC","OxPE","OxPG","OxPI","OxPS",
    "EtherOxPC","EtherOxPE",
    "PMeOH","PEtOH","PBtOH","MMPE","DMPE",
    "LNAPE","LNAPS",
    "Ac2PIM1","Ac2PIM2","Ac3PIM2","Ac4PIM2",
    "GPNAE",
    "LPC_d5","LPE_d5","LPG_d5","LPI_d5","LPS_d5",
    "PC_d5","PE_d5","PG_d5","PI_d5","PS_d5","bmPC"
  )
  sterol_lipids <- c(
    "CE","Cholesterol","CholesterolSulfate",
    "SHex","SSulfate","BAHex","BASulfate",
    "BRSE","CASE","SISE","STSE","EGSE","DEGSE","DSMSE",
    "AHexCS","AHexBRS","AHexCAS","AHexSIS","AHexSTS",
    "BRSLPHex","BRSPHex","CASLPHex","CASPHex","CSLPHex","CSPHex",
    "SISLPHex","SISPHex","STSLPHex","STSPHex","SPE","SPEHex","SPGHex",
    "DCAE","GDCAE","GLCAE","TDCAE","TLCAE","LCAE","KLCAE","KDCAE",
    "Vitamin_D","Vitamin D","BileAcid","ST",
    "CE_d7"
  )
  prenol_lipids <- c("CoQ","Vitamin_E","Vitamin E","VAE")
  sphingolipids <- c(
    "PhytoSph","DHSph","Sph",
    "Cer_ADS","Cer_AS","Cer_BDS","Cer_BS","Cer_NDS","Cer_NS","Cer_NP","Cer_AP",
    "Cer_ABP","Cer_NH","Cer_AH","Cer_NH_d9","Cer_AH_d9",
    "Cer_OS","Cer_HS","Cer_HDS","Cer_NDOS","Cer_EODS","Cer_EOS",
    "Cer-ADS","Cer-AS","Cer-BDS","Cer-BS","Cer-NDS","Cer-NS","Cer-NP","Cer-AP",
    "Cer-ABP","Cer-NH","Cer-AH","Cer-NH_d9","Cer-AH_d9",
    "Cer-OS","Cer-HS","Cer-HDS","Cer-NDOS","CerP","Cer-EODS","Cer-EOS",
    "HexCer-NS","HexCer-NDS","HexCer-AP","HexCer-HS","HexCer-HDS","HexCer-EOS",
    "HexCer_NS","HexCer_NDS","HexCer_AP","HexCer_HS","HexCer_HDS","HexCer_EOS",
    "Hex2Cer","Hex3Cer",
    "Cer_EBDS","Cer-EBDS","AHexCer","ASHexCer","ASM",
    "PI-Cer","PE-Cer","PI_Cer","PE_Cer",
    "PI-Cer+O","PE-Cer+O","PI_Cer+O","PE_Cer+O","MIPC",
    "SM","SHexCer","SL","SM+O","SHexCer+O","SL+O","GM3","GM3[NeuAc]",
    "GD1a","GD1b","GD2","GD3","GM1","GT1b","GQ1b","NGcGM3",
    "Cer_NS_d7","SM_d9"
  )

  if (lipidclass %in% fatty_acyls)        return("FattyAcyls")
  if (lipidclass %in% glycerolipids)      return("Glycerolipids")
  if (lipidclass %in% glycerophospholipids) return("Glycerophospholipids")
  if (lipidclass %in% sterol_lipids)      return("SterolLipids")
  if (lipidclass %in% prenol_lipids)      return("PrenolLipids")
  if (lipidclass %in% sphingolipids)      return("Sphingolipids")
  "Unassigned lipid"
}


# -----------------------------------------------------------------------------
# 17. convert_msdial_class_to_lbm_class_enum_vs2
#     Maps lipid class string to canonical class string (normalises variants).
#     In R we use strings rather than enums; this function returns the
#     canonical form used throughout the dispatcher.
#     Mirrors ConvertMsdialClassDefinitionToLbmClassEnumVS2 → then
#     ConvertLbmClassEnumToMsdialClassDefinitionVS2.
# -----------------------------------------------------------------------------
convert_msdial_class_to_lbm_class_enum_vs2 <- function(lipidclass) {
  # Normalise dash ↔ underscore variants for ceramides
  canonical <- switch(lipidclass,
    # Ceramides (underscore input → canonical dash form)
    "Cer_ADS"  = "Cer-ADS",  "Cer_AS"   = "Cer-AS",
    "Cer_BDS"  = "Cer-BDS",  "Cer_BS"   = "Cer-BS",
    "Cer_NDS"  = "Cer-NDS",  "Cer_NS"   = "Cer-NS",
    "Cer_NP"   = "Cer-NP",   "Cer_AP"   = "Cer-AP",
    "Cer_ABP"  = "Cer-ABP",  "Cer_NH"   = "Cer-NH",
    "Cer_AH"   = "Cer-AH",   "Cer_NH_d9" = "Cer-NH_d9",
    "Cer_AH_d9" = "Cer-AH_d9",
    "Cer_OS"   = "Cer-OS",   "Cer_HS"   = "Cer-HS",
    "Cer_HDS"  = "Cer-HDS",  "Cer_NDOS" = "Cer-NDOS",
    "Cer_EODS" = "Cer-EODS", "Cer_EOS"  = "Cer-EOS",
    "Cer_EBDS" = "Cer-EBDS",
    "HexCer_NS"  = "HexCer-NS",  "HexCer_NDS" = "HexCer-NDS",
    "HexCer_AP"  = "HexCer-AP",  "HexCer_HS"  = "HexCer-HS",
    "HexCer_HDS" = "HexCer-HDS", "HexCer_EOS" = "HexCer-EOS",
    "PE_Cer"   = "PE-Cer",   "PI_Cer"   = "PI-Cer",
    "PE_Cer+O" = "PE-Cer",   "PI_Cer+O" = "PI-Cer",
    "PE-Cer+O" = "PE-Cer",   "PI-Cer+O" = "PI-Cer",
    "SM+O"     = "SM",
    "SHexCer+O" = "SHexCer", "SL+O"    = "SL",
    # EtherPE alias forms
    "EtherPE_O" = "EtherPE", "EtherPE_P" = "EtherPE",
    # Vitamin aliases
    "Vitamin E" = "Vitamin_E", "Vitamin D" = "Vitamin_D",
    # MLCL old name
    "LCL" = "MLCL",
    # Default: return as-is (most classes are already canonical)
    lipidclass
  )
  if (is.null(canonical)) lipidclass else canonical
}


# -----------------------------------------------------------------------------
# 18. convert_msdial_lipidname_to_lipid_molecule_object_vs2
#     Creates and populates a LipidMolecule list from a MoleculeMsReference.
#     Mirrors ConvertMsdialLipidnameToLipidMoleculeObjectVS2(MoleculeMsReference).
# -----------------------------------------------------------------------------
convert_msdial_lipidname_to_lipid_molecule_object_vs2 <- function(mol_ms_ref) {
  molecule   <- new_lipid_molecule()
  lipidclass <- mol_ms_ref$CompoundClass

  lipid_category <- convert_msdial_class_to_super_class_vs2(lipidclass)
  lbm_class      <- convert_msdial_class_to_lbm_class_enum_vs2(lipidclass)

  if (lipid_category %in% c("FattyAcyls", "Glycerolipids",
                              "Glycerophospholipids", "Sphingolipids") ||
      lbm_class == "VAE") {
    molecule <- set_lipid_acyl_chain_properties(molecule, mol_ms_ref)

  } else if (lipid_category %in% c("SterolLipids", "PrenolLipids")) {
    is_simple_sterol <- lbm_class %in% c("Vitamin_D", "Vitamin_E",
                                          "SHex", "SSulfate", "BAHex",
                                          "BASulfate", "BileAcid") ||
                        mol_ms_ref$Name %in% c("Cholesterol", "CholesterolSulfate")
    if (is_simple_sterol) {
      molecule <- set_single_lipid_structure(molecule, mol_ms_ref)
    } else if (lbm_class == "CoQ") {
      molecule <- set_coq_molecule(molecule, mol_ms_ref)
    } else {
      molecule <- set_lipid_acyl_chain_properties(molecule, mol_ms_ref)
    }
  }

  molecule$LipidClass    <- lbm_class
  molecule$LipidCategory <- lipid_category
  molecule$LipidSubclass <- lipidclass

  if (is.null(molecule$LipidName) || nchar(molecule$LipidName) == 0 ||
      is.null(molecule$Adduct)) {
    molecule$IsValidatedFormat <- FALSE
  } else {
    molecule$IsValidatedFormat <- TRUE
  }

  molecule
}


# -----------------------------------------------------------------------------
# 19. is_compared_available
#     Returns TRUE iff both objects have non-NULL, non-empty Spectrum lists.
#     Mirrors IsComparedAvailable(IMSScanProperty, IMSScanProperty).
# -----------------------------------------------------------------------------
is_compared_available <- function(obj1, obj2) {
  !is.null(obj1$Spectrum) && !is.null(obj2$Spectrum) &&
    length(obj1$Spectrum) > 0L  && length(obj2$Spectrum) > 0L
}


# -----------------------------------------------------------------------------
# 20. get_matched_peaks_scores
#     Sliding-window fragment-matching score.
#     Returns numeric[2]: c(matched_ratio, matched_count), or c(-1,-1) / c(0,0).
#     Mirrors GetMatchedPeaksScores(List<SpectrumPeak>, List<SpectrumPeak>, ...).
#
#     scan     — list with $Spectrum (list of lists, each with $Mass/$Intensity)
#     ref_spec — same format
# -----------------------------------------------------------------------------
get_matched_peaks_scores <- function(scan, ref_spec, bin, mass_begin, mass_end) {
  if (!is_compared_available(scan, ref_spec)) return(c(-1, -1))

  # Extract mass/intensity vectors
  .extract <- function(spec) {
    if (is.data.frame(spec)) {
      list(mass = spec$Mass, intensity = spec$Intensity)
    } else {
      list(mass      = vapply(spec, `[[`, numeric(1), "Mass"),
           intensity = vapply(spec, `[[`, numeric(1), "Intensity"))
    }
  }
  p1 <- .extract(scan$Spectrum)
  p2 <- .extract(ref_spec$Spectrum)
  m1 <- p1$mass;  i1 <- p1$intensity
  m2 <- p2$mass;  i2 <- p2$intensity
  n1 <- length(m1); n2 <- length(m2)

  min_mz <- m2[1L];   max_mz  <- m2[n2]
  if (mass_begin > min_mz) min_mz <- mass_begin
  if (max_mz > mass_end)   max_mz <- mass_end

  max_lib_int <- max(i2)
  focused_mz  <- min_mz
  ri_m <- 1L;  ri_l <- 1L
  counter <- 0L;  lib_counter <- 0L

  while (focused_mz <= max_mz) {
    # ---- library window sum ----
    sum_l    <- 0.0
    new_ri_l <- ri_l
    for (i in ri_l:n2) {
      mi <- m2[i]
      if (mi < focused_mz - bin) next
      else if (mi < focused_mz + bin) { sum_l <- sum_l + i2[i] }
      else { new_ri_l <- i; break }
    }
    ri_l <- new_ri_l
    if (sum_l >= 0.01 * max_lib_int) lib_counter <- lib_counter + 1L

    # ---- experimental window sum ----
    sum_m    <- 0.0
    new_ri_m <- ri_m
    for (i in ri_m:n1) {
      mi <- m1[i]
      if (mi < focused_mz - bin) next
      else if (mi < focused_mz + bin) { sum_m <- sum_m + i1[i] }
      else { new_ri_m <- i; break }
    }
    ri_m <- new_ri_m
    if (sum_m > 0 && sum_l >= 0.01 * max_lib_int) counter <- counter + 1L

    if (focused_mz + bin > m2[n2]) break
    focused_mz <- m2[ri_l]
  }

  if (lib_counter == 0L) c(0, 0) else c(counter / lib_counter, counter)
}


# -----------------------------------------------------------------------------
# 21. get_lipid_molecule_annotation_result
#     Dispatches to the appropriate judge_if_* function based on lipid class.
#     Returns a LipidMolecule result list (from judge_if_*) or NULL.
#     Requires LipidMsmsCharacterization.R to be sourced.
#     Mirrors GetLipidMoleculeAnnotationResult().
# -----------------------------------------------------------------------------
get_lipid_molecule_annotation_result <- function(ms_scan_prop, molecule, ms2_tol) {

  lc  <- molecule$LipidClass      # string, e.g. "PC", "Cer-NS"
  mz  <- molecule$Mz
  add <- molecule$Adduct           # string, e.g. "[M+H]+"

  tot_c  <- molecule$TotalCarbonCount
  tot_db <- molecule$TotalDoubleBondCount
  tot_ox <- molecule$TotalOxidizedCount

  sn1_c  <- molecule$Sn1CarbonCount
  sn1_db <- molecule$Sn1DoubleBondCount
  sn1_ox <- molecule$Sn1Oxidizedount
  sn2_ox <- molecule$Sn2Oxidizedount

  # Common 2-sn call (sn1 min==sn1 max, same for db)
  .call2 <- function(fn)
    fn(ms_scan_prop, ms2_tol, mz,
       tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add)

  if (lc == "PC")
    return(judge_if_phosphatidylcholine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PE")
    return(judge_if_phosphatidylethanolamine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PS")
    return(judge_if_phosphatidylserine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PG")
    return(judge_if_phosphatidylglycerol(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "BMP")
    return(judge_if_bismonoacylglycerophosphate(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PI")
    return(judge_if_phosphatidylinositol(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SM") {
    if (grepl("O3", molecule$TotalChainString, fixed = TRUE))
      return(judge_if_sphingomyelin_phyto(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
    else
      return(judge_if_sphingomyelin(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  }
  if (lc == "LNAPE")
    return(judge_if_nacylphosphatidylethanolamine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LNAPS")
    return(judge_if_nacylphosphatidylserine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "CE")
    return(judge_if_cholesteryl_ester(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "CAR")
    return(judge_if_acylcarnitine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "DG")
    return(judge_if_dag(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "MG")
    return(judge_if_mag(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "MGDG")
    return(judge_if_mgdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DGDG")
    return(judge_if_dgdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PMeOH")
    return(judge_if_pmeoh(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PEtOH")
    return(judge_if_petoh(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PBtOH")
    return(judge_if_pbtoh(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPC")
    return(judge_if_lysopc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPE")
    return(judge_if_lysope(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PA")
    return(judge_if_phosphatidic_acid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPA")
    return(judge_if_lysopa(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPG")
    return(judge_if_lysopg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPI")
    return(judge_if_lysopi(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPS")
    return(judge_if_lysops(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherPC")
    return(judge_if_etherpc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherPE")
    return(judge_if_etherpe(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherLPC")
    return(judge_if_etherlysopc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherLPE")
    return(judge_if_etherlysope(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "OxPC")
    return(judge_if_oxpc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "OxPE")
    return(judge_if_oxpe(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "OxPG")
    return(judge_if_oxpg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "OxPI")
    return(judge_if_oxpi(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "OxPS")
    return(judge_if_oxps(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "EtherMGDG")
    return(judge_if_ethermgdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherDGDG")
    return(judge_if_etherdgdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DGTS")
    return(judge_if_dgts(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LDGTS")
    return(judge_if_ldgts(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DGCC")
    return(judge_if_dgcc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LDGCC")
    return(judge_if_ldgcc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DGGA")
    return(judge_if_glcadg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SQDG")
    return(judge_if_sqdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DLCL")
    return(judge_if_dilysocardiolipin(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "FA")
    return(judge_if_fattyacid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "OxFA")
    return(judge_if_oxfattyacid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add, tot_ox))
  if (lc == "FAHFA")
    return(judge_if_fahfa(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DMEDFAHFA")
    return(judge_if_fahfa_dmed(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DMEDFA")
    return(judge_if_dmed_fattyacid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DMEDOxFA")
    return(judge_if_dmed_oxfattyacid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add, tot_ox))
  if (lc == "EtherOxPC")
    return(judge_if_etheroxpc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "EtherOxPE")
    return(judge_if_etheroxpe(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add,
             tot_ox, sn1_ox, sn2_ox))
  if (lc == "Cer-NS")
    return(judge_if_ceramidens(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-NDS")
    return(judge_if_ceramidends(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "HexCer-NS")
    return(judge_if_hexceramidens(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "HexCer-NDS")
    return(judge_if_hexceramidends(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Hex2Cer")
    return(judge_if_hexhexceramidens(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Hex3Cer")
    return(judge_if_hexhexhexceramidens(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-AP")
    return(judge_if_ceramideap(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-ABP")
    return(judge_if_ceramideabp(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-AH")
    return(judge_if_ceramideah(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-AH_d9")
    return(judge_if_ceramideahd9(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-NH_d9")
    return(judge_if_ceramidenhd9(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-NH")
    return(judge_if_ceramidenh(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "HexCer-AP")
    return(judge_if_hexceramideap(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SHexCer")
    return(judge_if_shexcer(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add, tot_ox))
  if (lc == "GM3")
    return(judge_if_gm3(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DHSph")
    return(judge_if_sphinganine(ms_scan_prop, ms2_tol, mz,
             molecule$TotalCarbonCount, molecule$TotalDoubleBondCount, add))
  if (lc == "Sph")
    return(judge_if_sphingosine(ms_scan_prop, ms2_tol, mz,
             molecule$TotalCarbonCount, molecule$TotalDoubleBondCount, add))
  if (lc == "PhytoSph")
    return(judge_if_phytosphingosine(ms_scan_prop, ms2_tol, mz,
             molecule$TotalCarbonCount, molecule$TotalDoubleBondCount, add))

  # 3-chain dispatches (need sn2)
  if (lc %in% c("TG","ADGGA","HBMP","EtherTG","MLCL",
                 "Cer-EOS","Cer-EODS","HexCer-EOS","ASM","Cer-EBDS",
                 "AHexCer","ASHexCer","OxTG","TG_EST","TG_d5")) {
    sn2_c  <- molecule$Sn2CarbonCount
    sn2_db <- molecule$Sn2DoubleBondCount
    if (lc == "TG")
      return(judge_if_triacylglycerol(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "ADGGA")
      return(judge_if_acylglcadg(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "HBMP")
      return(judge_if_hemiismonoacylglycerophosphate(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "EtherTG")
      return(judge_if_ethertag(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "MLCL")
      return(judge_if_lysocardiolipin(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "Cer-EOS")
      return(judge_if_ceramideeos(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "Cer-EODS")
      return(judge_if_ceramideeods(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "HexCer-EOS")
      return(judge_if_hexceramideeos(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "ASM")
      return(judge_if_acylsm(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "Cer-EBDS")
      return(judge_if_acylcerbds(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "AHexCer")
      return(judge_if_acylhexcer(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, tot_ox, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "ASHexCer")
      return(judge_if_ashexcer(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, tot_ox, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
    if (lc == "OxTG")
      return(judge_if_ox_triacylglycerol(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, tot_ox, add))
    if (lc == "TG_EST") {
      sn3_c  <- molecule$Sn3CarbonCount
      sn3_db <- molecule$Sn3DoubleBondCount
      return(judge_if_fahfa_triacylglycerol(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db,
               sn3_c, sn3_c, sn3_db, sn3_db, add))
    }
    if (lc == "TG_d5")
      return(judge_if_tg_d5(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db, add))
  }

  # CL (cardiolipin)
  if (lc == "CL") {
    sn2_c  <- molecule$Sn2CarbonCount
    sn2_db <- molecule$Sn2DoubleBondCount
    sn3_c  <- molecule$Sn3CarbonCount
    sn3_db <- molecule$Sn3DoubleBondCount
    if (sn3_c < 1L)
      return(judge_if_cardiolipin_simple(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
    else
      return(judge_if_cardiolipin(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db,
               sn2_c, sn2_c, sn2_db, sn2_db,
               sn3_c, sn3_c, sn3_db, sn3_db, add))
  }

  if (lc == "EtherPI")
    return(judge_if_etherpi(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherPS")
    return(judge_if_etherps(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherDG")
    return(judge_if_ether_dag(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PI-Cer")
    return(judge_if_picermide(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add, tot_ox))
  if (lc == "PE-Cer")
    return(judge_if_pecermide(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add, tot_ox))
  if (lc == "DCAE")
    return(judge_if_dcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "GDCAE")
    return(judge_if_gdcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "GLCAE")
    return(judge_if_glcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "TDCAE")
    return(judge_if_tdcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "TLCAE")
    return(judge_if_tlcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "NAE")
    return(judge_if_anandamide(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "NAGly") {
    if (tot_c == sn1_c)
      return(judge_if_nacyl_gly_ox_fa(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, tot_ox, add))
    else
      return(judge_if_fahfamidegly(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  }
  if (lc == "NAGlySer") {
    if (tot_c == sn1_c)
      return(judge_if_nacyl_glyser_ox_fa(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, tot_ox, add))
    else
      return(judge_if_fahfamideglyser(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  }
  if (lc == "SL")
    return(judge_if_sulfonolipid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add, tot_ox))
  if (lc == "EtherPG")
    return(judge_if_etherpg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherLPG")
    return(judge_if_etherlysopg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "CoQ")
    return(judge_if_coenzymeq(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "Vitamin_E")
    return(judge_if_vitamin_e_molecules(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "Vitamin_D")
    return(judge_if_vitamin_d_molecules(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "VAE")
    return(judge_if_vitaminaestermolecules(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "NAOrn") {
    if (tot_c == sn1_c)
      return(judge_if_nacyl_orn_ox_fa(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, tot_ox, add))
    else
      return(judge_if_fahfamideorn(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  }
  if (lc == "BRSE")
    return(judge_if_brse_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "CASE")
    return(judge_if_case_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "SISE")
    return(judge_if_sise_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "STSE")
    return(judge_if_stse_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "AHexBRS")
    return(judge_if_ahexbrs_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "AHexCAS")
    return(judge_if_ahexcase_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "AHexCS")
    return(judge_if_ahexce_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "AHexSIS")
    return(judge_if_ahexsise_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "AHexSTS")
    return(judge_if_ahexstse_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "Cer-AS")
    return(judge_if_ceramideas(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-ADS")
    return(judge_if_ceramideads(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-BS")
    return(judge_if_ceramidebs(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-BDS")
    return(judge_if_ceramidebds(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-NP")
    return(judge_if_ceramidenp(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-OS")
    return(judge_if_ceramideos(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc %in% c("Cer-HS", "Cer-HDS"))
    return(judge_if_ceramideo(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "Cer-NDOS")
    return(judge_if_ceramidedos(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc %in% c("HexCer-HS", "HexCer-HDS"))
    return(judge_if_hexceramideo(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SHex")
    return(judge_if_sterol_hexoside(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc == "BAHex")
    return(judge_if_sterol_hexoside(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc == "SSulfate")
    return(judge_if_sterol_sulfate(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc == "BASulfate")
    return(judge_if_sterol_sulfate(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc == "CerP")
    return(judge_if_ceramide_phosphate(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SMGDG")
    return(judge_if_smgdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EtherSMGDG")
    return(judge_if_ethersmgdg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LCAE")
    return(judge_if_lcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "KLCAE")
    return(judge_if_klcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "KDCAE")
    return(judge_if_kdcae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add, tot_ox))
  if (lc == "DMPE")
    return(judge_if_di_methyl_pe(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "MMPE")
    return(judge_if_mono_methyl_pe(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "MIPC")
    return(judge_if_mipc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "EGSE")
    return(judge_if_egse_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "DEGSE")
    return(judge_if_degse_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "DSMSE")
    return(judge_if_dsmse_species(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "GPNAE")
    return(judge_if_gpnae(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "MGMG")
    return(judge_if_mgmg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "DGMG")
    return(judge_if_dgmg(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "GD1a")
    return(judge_if_gd1a(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "GD1b")
    return(judge_if_gd1b(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "GD2")
    return(judge_if_gd2(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "GD3")
    return(judge_if_gd3(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "GM1")
    return(judge_if_gm1(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "GQ1b")
    return(judge_if_gq1b(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "GT1b")
    return(judge_if_gt1b(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "NGcGM3")
    return(judge_if_ngc_gm3(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "ST")
    return(judge_if_nochain_sterol(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc %in% c("CSLPHex","BRSLPHex","CASLPHex","SISLPHex","STSLPHex"))
    return(judge_if_steroid_with_lpa(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc %in% c("CSPHex","BRSPHex","CASPHex","SISPHex","STSPHex"))
    return(judge_if_steroid_with_pa(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SPE")
    return(judge_if_spe_species(molecule$LipidName, molecule$LipidClass,
             ms_scan_prop, ms2_tol, mz, tot_c, tot_db, add))
  if (lc == "NAPhe")
    return(judge_if_nacyl_phe_fa(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NATau")
    return(judge_if_nacyl_tau_fa(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "PT")
    return(judge_if_phosphatidylthreonine(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))

  # Stable-isotope labelled standards
  if (lc == "PC_d5")
    return(judge_if_pc_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PE_d5")
    return(judge_if_pe_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PS_d5")
    return(judge_if_ps_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PG_d5")
    return(judge_if_pg_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "PI_d5")
    return(judge_if_pi_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPC_d5")
    return(judge_if_lysopc_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPE_d5")
    return(judge_if_lysope_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPG_d5")
    return(judge_if_lysopg_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPI_d5")
    return(judge_if_lysopi_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "LPS_d5")
    return(judge_if_lysops_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "DG_d5")
    return(judge_if_dag_d5(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "SM_d9")
    return(judge_if_sm_d9(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "CE_d7")
    return(judge_if_ce_d7(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, add))
  if (lc == "Cer_NS_d7")
    return(judge_if_cer_ns_d7(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "bmPC")
    return(judge_if_beta_methyl_pc(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "NA5HT")
    return(judge_if_nacyl_5ht(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "WE")
    return(judge_if_wax_ester(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, sn1_c, sn1_c, sn1_db, sn1_db, add))
  if (lc == "NAAla")
    return(judge_if_nacyl_ala(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NAGln")
    return(judge_if_nacyl_gln(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NALeu")
    return(judge_if_nacyl_leu(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NAVal")
    return(judge_if_nacyl_val(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NASer")
    return(judge_if_nacyl_ser(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NAAnt")
    return(judge_if_nacyl_anthranilicacid(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "BisMeLPA")
    return(judge_if_bismelpa(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "NATryA") {
    if (tot_c < 29L)
      return(judge_if_nacyl_trya(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, tot_ox, sn1_c, sn1_c, sn1_db, sn1_db, add))
    else
      return(judge_if_fahfamide_trya(ms_scan_prop, ms2_tol, mz,
               tot_c, tot_db, sn1_c, sn1_c, sn1_db, sn1_db, add))
  }
  if (lc == "NAGABA")
    return(judge_if_nacyl_gaba(ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "SPEHex")
    return(judge_if_spehex(molecule$LipidName, ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))
  if (lc == "SPGHex")
    return(judge_if_spghex(molecule$LipidName, ms_scan_prop, ms2_tol, mz,
             tot_c, tot_db, tot_ox, add))

  NULL   # default
}


# -----------------------------------------------------------------------------
# 22. get_lipidomics_matched_peaks_scores  — MAIN FUNCTION
#     Returns numeric[2]: c(score, matched_count).
#     score may be -1 (unavailable), 0 (no match), 1 (class-level bonus),
#     or 2 (molecular-species-level bonus).
#     Requires LipidMsmsCharacterization.R to be sourced.
#
#     ms_scan_prop — IMSScanProperty: list($Spectrum, $IonMode, $PrecursorMz)
#     mol_ms_ref   — MoleculeMsReference: list($Spectrum, $CompoundClass,
#                       $Comment, $Name, $AdductType, $IonMode, $PrecursorMz,
#                       $SMILES, $InChIKey, $Formula)
#     bin          — mass bin width (Da) for fragment matching
#     mass_begin   — lower m/z limit for fragment window
#     mass_end     — upper m/z limit for fragment window
# -----------------------------------------------------------------------------
get_lipidomics_matched_peaks_scores <- function(ms_scan_prop, mol_ms_ref,
                                                bin, mass_begin, mass_end) {
  if (!is_compared_available(ms_scan_prop, mol_ms_ref))
    return(c(-1, -1))

  result_array <- get_matched_peaks_scores(ms_scan_prop, mol_ms_ref,
                                           bin, mass_begin, mass_end)
  comp_class <- mol_ms_ref$CompoundClass
  comment    <- mol_ms_ref$Comment

  if (!identical(comment, "SPLASH") &&
      !identical(comp_class, "Unknown") &&
      !identical(comp_class, "Others")) {

    molecule <- convert_msdial_lipidname_to_lipid_molecule_object_vs2(mol_ms_ref)
    if (is.null(molecule) || is.null(molecule$Adduct))
      return(result_array)

    # Special case: EtherPE with only 3 library peaks in positive mode
    # (too few fragments for reliable substructure scoring)
    ion_mode <- if (is.numeric(ms_scan_prop$IonMode))
                  ifelse(ms_scan_prop$IonMode == 0, "Negative", "Positive")
                else ms_scan_prop$IonMode
    if (molecule$LipidClass == "EtherPE" &&
        length(mol_ms_ref$Spectrum) == 3L &&
        ion_mode == "Positive")
      return(result_array)

    result <- tryCatch(
      get_lipid_molecule_annotation_result(ms_scan_prop, molecule, bin),
      error = function(e) NULL
    )

    if (is.list(result) && length(result) > 0L && !is.null(result$AnnotationLevel)) {
      annotation_level <- suppressWarnings(as.integer(result$AnnotationLevel[[1]]))
      if (is.na(annotation_level)) annotation_level <- 0L

      if (annotation_level == 1L) {
        # SM with phyto-sphingoid base (3O) gets extra bonus
        if (comp_class == "SM" &&
            (grepl("3O", molecule$LipidName, fixed = TRUE) ||
             grepl("O3", molecule$LipidName, fixed = TRUE))) {
          result_array[1] <- 2.0
        } else {
          result_array[1] <- 1.0
        }
      } else if (annotation_level >= 2L) {
        result_array[1] <- 2.0
      }
    }
  }
  # else: SPLASH / Unknown / Others — return raw matched-peaks score as-is
  result_array
}
