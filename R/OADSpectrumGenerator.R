# OAD theoretical spectrum generation translated from MSDIAL.

.oad_sg_mass_diff <- list(
  carbon = 12,
  hydrogen = 1.00782503223,
  nitrogen = 14.00307400443,
  oxygen = 15.99491461957,
  phosphorus = 30.97376199842
)
.oad_sg_proton <- 1.007276466621
.oad_sg_electron <- 0.00054858026
.oad_sg_deuterium <- 2.01410177812
.oad_sg_water <- 2 * .oad_sg_mass_diff$hydrogen + .oad_sg_mass_diff$oxygen
.oad_sg_na <- 22.989218
.oad_sg_nh4 <- 18.033823
.oad_sg_nh3 <- 3 * .oad_sg_mass_diff$hydrogen + .oad_sg_mass_diff$nitrogen
.oad_sg_formate <- 44.998201
.oad_sg_acetate <- 59.013851
.oad_sg_bicarbonate <- 60.021129

.oad_sg_or <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L || is.na(value[1L])) fallback else value[1L]
}

.oad_sg_mass <- function(carbon = 0, hydrogen = 0, nitrogen = 0,
                         oxygen = 0, phosphorus = 0) {
  d <- .oad_sg_mass_diff
  carbon * d$carbon + hydrogen * d$hydrogen + nitrogen * d$nitrogen +
    oxygen * d$oxygen + phosphorus * d$phosphorus
}

.oad_sg_field <- function(x, names, default = NULL) {
  if (is.null(x)) return(default)
  for (name in names) {
    value <- x[[name]]
    if (!is.null(value) && length(value) > 0L && !is.na(value[1L])) {
      return(value[1L])
    }
  }
  default
}

.oad_sg_adduct_name <- function(adduct, lipid = NULL) {
  value <- if (is.character(adduct)) adduct else .oad_sg_field(adduct, c("AdductIonName", "name"))
  .oad_sg_or(value, .oad_sg_field(lipid, c("AdductIonName", "PrecursorType"), "[M+H]+"))
}

.oad_sg_ion_mz <- function(mass, adduct) {
  name <- .oad_sg_adduct_name(adduct)
  shift <- switch(name,
    "[M+H]+" = .oad_sg_proton,
    "[M+Na]+" = .oad_sg_na,
    "[M+NH4]+" = .oad_sg_nh4,
    "[M+H-H2O]+" = .oad_sg_proton - .oad_sg_water,
    "[M-H2O+H]+" = .oad_sg_proton - .oad_sg_water,
    "[M-H]-" = -.oad_sg_proton,
    "[M+HCOO]-" = .oad_sg_formate - .oad_sg_proton,
    "[M+CH3COO]-" = .oad_sg_acetate - .oad_sg_proton,
    "[M+HCO3]-" = .oad_sg_bicarbonate - .oad_sg_proton,
    .oad_sg_proton
  )
  mass + shift
}

.oad_sg_extract_class <- function(name) {
  sub("[[:space:](].*$", "", as.character(name)[1L])
}

.oad_sg_parse_chains <- function(name, lipid = NULL) {
  matches <- gregexpr("[0-9]+:[0-9]+(?:\\([^)]*\\))?", name, perl = TRUE)
  values <- regmatches(name, matches)[[1L]]
  if (length(values) == 0L || identical(values, character(0))) {
    carbon <- as.integer(.oad_sg_field(lipid, c("TotalCarbonCount", "TotalCarbon"), 0L))
    db <- as.integer(.oad_sg_field(lipid, c("TotalDoubleBondCount", "TotalDoubleBond"), 0L))
    return(data.frame(carbon = carbon, double_bond = db, positions = I(list(integer())),
                      stringsAsFactors = FALSE))
  }
  rows <- lapply(values, function(value) {
    parts <- regexec("^([0-9]+):([0-9]+)(?:\\(([^)]*)\\))?$", value, perl = TRUE)
    fields <- regmatches(value, parts)[[1L]]
    positions <- integer()
    if (length(fields) >= 4L && nzchar(fields[4L])) {
      positions <- as.integer(strsplit(gsub("[^0-9,]", "", fields[4L]), ",", fixed = TRUE)[[1L]])
      positions <- positions[!is.na(positions)]
    }
    list(carbon = as.integer(fields[2L]), double_bond = as.integer(fields[3L]), positions = positions)
  })
  data.frame(carbon = vapply(rows, `[[`, integer(1), "carbon"),
             double_bond = vapply(rows, `[[`, integer(1), "double_bond"),
             positions = I(lapply(rows, `[[`, "positions")),
             stringsAsFactors = FALSE)
}

.oad_sg_chain_mass <- function(carbon, double_bond) {
  .oad_sg_mass(carbon, 2 * carbon - 2 * double_bond, oxygen = 2)
}

.oad_sg_peak <- function(mass, intensity, comment, oad_id = NA_character_, type = "doublebond") {
  data.frame(Mass = as.numeric(mass), Intensity = as.numeric(intensity),
             Comment = as.character(comment), OadId = as.character(oad_id),
             SpectrumComment = type, stringsAsFactors = FALSE)
}

.oad_sg_peak_mz <- function(mz, intensity, comment, oad_id = NA_character_, type = "doublebond") {
  data.frame(Mass = as.numeric(mz), Intensity = as.numeric(intensity),
             Comment = as.character(comment), OadId = as.character(oad_id),
             SpectrumComment = type, stringsAsFactors = FALSE)
}

.oad_sg_merge <- function(peaks) {
  if (length(peaks) == 0L) {
    return(data.frame(Mass = numeric(), Intensity = numeric(), Comment = character(),
                      OadId = character(), SpectrumComment = character(), stringsAsFactors = FALSE))
  }
  peaks <- do.call(rbind, peaks)
  key <- formatC(peaks$Mass, format = "f", digits = 10)
  groups <- split(seq_len(nrow(peaks)), key)
  result <- lapply(groups, function(index) {
    first <- peaks[index[1L], , drop = FALSE]
    first$Intensity <- sum(peaks$Intensity[index])
    first$Comment <- paste(unique(peaks$Comment[index]), collapse = ", ")
    first$OadId <- paste(unique(peaks$OadId[index][!is.na(peaks$OadId[index])]), collapse = ",")
    if (!nzchar(first$OadId)) first$OadId <- NA_character_
    first
  })
  result <- do.call(rbind, result)
  result[order(result$Mass), , drop = FALSE]
}

.oad_sg_class_peaks <- function(mass, class, adduct, chains = NULL) {
  peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", NA, "precursor"))
  name <- .oad_sg_adduct_name(adduct)
  header <- switch(class,
    PC = .oad_sg_mass(5, 14, 1, 4, 1),
    LPC = .oad_sg_mass(5, 14, 1, 4, 1),
    PE = .oad_sg_mass(2, 8, 1, 4, 1),
    LPE = .oad_sg_mass(2, 8, 1, 4, 1),
    PG = .oad_sg_mass(3, 9, 0, 6, 1),
    LPG = .oad_sg_mass(3, 9, 0, 6, 1),
    PS = .oad_sg_mass(3, 8, 1, 6, 1),
    LPS = .oad_sg_mass(3, 8, 1, 6, 1),
    PI = .oad_sg_mass(6, 13, 0, 9, 1),
    LPI = .oad_sg_mass(6, 13, 0, 9, 1),
    PA = .oad_sg_mass(3, 6, 0, 5, 1),
    LPA = .oad_sg_mass(3, 6, 0, 5, 1),
    LPS = 0,
    0
  )
  if (header > 0 && class != "LPS") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 100,
                                                  "Header", NA, "metaboliteclass")
  }
  if (name %in% c("[M+H]+", "[M+Na]+", "[M+NH4]+") && class != "LPS") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass - .oad_sg_water, adduct), 100,
                                                  "Precursor -H2O", NA, "metaboliteclass")
  }
  if (class == "TG" && !is.null(chains) && nrow(chains) > 0L) {
    for (i in seq_len(nrow(chains))) {
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak(
        .oad_sg_ion_mz(mass - .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]) -
                         .oad_sg_mass(oxygen = 1), adduct),
        200, paste0(chains$carbon[i], ":", chains$double_bond[i], " loss"),
        NA_character_, "acylchain")
    }
  }
  if (class == "DG" && !is.null(chains) && nrow(chains) > 0L) {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(mass + .oad_sg_proton, 200,
                                                "[M+H]+", NA, "metaboliteclass")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(
      mass - .oad_sg_water + .oad_sg_proton, 999, "[M+H]+ -H2O",
      NA, "metaboliteclass")
    for (i in seq_len(nrow(chains))) {
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak(
        .oad_sg_ion_mz(mass - .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]) -
                         .oad_sg_mass(oxygen = 1) - .oad_sg_electron, adduct),
        50, paste0(chains$carbon[i], ":", chains$double_bond[i], " loss"),
        NA_character_, "acylchain")
    }
  }
  peaks
}

.oad_sg_dg_class_peaks <- function(mass, adduct, chains, labelled = FALSE) {
  peaks <- list()
  name <- .oad_sg_adduct_name(adduct)
  if (name == "[M+NH4]+") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999,
                                                  "Precursor", NA, "precursor")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(mass + .oad_sg_proton, 50,
                                                   "[M+H]+", "", "metaboliteclass")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      mass - .oad_sg_water + .oad_sg_proton, 200, "[M+H]+ -H2O", "", "metaboliteclass")
    if (!is.null(chains) && nrow(chains) > 0L) {
      for (i in seq_len(nrow(chains))) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
          mass - chain_mass - .oad_sg_mass_diff$oxygen - .oad_sg_electron,
          if (labelled) 500 else 50,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "", "acylchain")
        if (labelled) {
          peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
            mass - chain_mass - .oad_sg_mass_diff$oxygen - .oad_sg_electron +
              .oad_sg_mass_diff$hydrogen,
            200, paste0("-", chains$carbon[i], ":", chains$double_bond[i], "+H"), "", "acylchain")
        }
      }
    }
  } else {
    peaks[[1L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999,
                                 "Precursor", NA, "precursor")
  }
  peaks
}

.oad_sg_tgd5_class_peaks <- function(mass, adduct, chains) {
  peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 500,
                             "Precursor", NA, "precursor"))
  if (.oad_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass - .oad_sg_nh3, adduct),
                                                 250, "[M+H]+", NA, "precursor")
    for (i in seq_len(nrow(chains))) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
        mass - chain_mass - .oad_sg_mass_diff$oxygen, 999,
        paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), NA, "acylchain")
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
        chain_mass, 250, paste0(chains$carbon[i], ":", chains$double_bond[i]), NA, "acylchain")
    }
  }
  peaks
}

.oad_sg_ced7_class_peaks <- function(mass, adduct) {
  skeleton <- .oad_sg_mass(27, 38, 0, 0) + 7 * .oad_sg_deuterium
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999,
                             "Precursor", NA, "precursor"))
  if (name %in% c("[M+NH4]+", "[M-H]-")) {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      skeleton + .oad_sg_proton, 500, "skelton", NA, "metaboliteclass")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(
      .oad_sg_ion_mz(mass - .oad_sg_water, adduct), 50,
      "Precursor -H2O", NA, "metaboliteclass")
  } else if (name == "[M+Na]+") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      skeleton, 500, "skelton", NA, "metaboliteclass")
  }
  peaks
}

.oad_sg_cerns_d7_class_peaks <- function(mass, adduct, chains) {
  d7_balance <- 7 * (.oad_sg_deuterium - .oad_sg_mass_diff$hydrogen)
  h2o <- .oad_sg_water
  ch4o2 <- .oad_sg_mass(1, 4, 0, 2)
  ch5o <- .oad_sg_mass(1, 5, 0, 1)
  c2h3n <- .oad_sg_mass(2, 3, 1)
  acetate <- .oad_sg_acetate
  effective_mass <- mass + d7_balance
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(effective_mass, adduct), 999,
                             "Precursor", NA, "precursor"))
  if (name == "[M+H]+") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      .oad_sg_ion_mz(effective_mass, adduct) - 2 * h2o, 100,
      "Precursor-2H2O", NA, "metaboliteclass")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      .oad_sg_ion_mz(effective_mass, adduct) - h2o, 200,
      "Precursor-H2O", NA, "metaboliteclass")
    if (nrow(chains) >= 1L) {
      sphingo_mass <- .oad_sg_chain_mass(chains$carbon[1L], chains$double_bond[1L])
      for (loss in c(ch4o2, h2o, 2 * h2o)) {
        intensity <- if (identical(loss, 2 * h2o)) 500 else 200
        peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
          sphingo_mass + d7_balance - loss + 2 * .oad_sg_mass_diff$hydrogen,
          intensity, paste0(chains$carbon[1L], ":", chains$double_bond[1L], " loss"), NA, "acylchain")
      }
    }
  } else if (name == "[M+CH3COO]-") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      .oad_sg_ion_mz(effective_mass, adduct) - acetate - .oad_sg_mass_diff$hydrogen,
      200, "Precursor-CH3COO", NA, "metaboliteclass")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      .oad_sg_ion_mz(effective_mass, adduct) - acetate - ch5o,
      200, "Precursor-CH3COO-CH5O", NA, "metaboliteclass")
    if (nrow(chains) >= 2L) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
        chain_mass - 2 * .oad_sg_mass_diff$hydrogen, 50,
        paste0(chains$carbon[2L], ":", chains$double_bond[2L], "-CH4O2"), NA, "acylchain")
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
        chain_mass + c2h3n, 100,
        paste0(chains$carbon[2L], ":", chains$double_bond[2L], "-C2H3N"), NA, "acylchain")
    }
  } else if (name == "[M-H]-") {
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      .oad_sg_ion_mz(effective_mass, adduct) - h2o - ch5o - 2 * .oad_sg_mass_diff$hydrogen,
      200, "Precursor-H2O-CH5O", NA, "metaboliteclass")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
      .oad_sg_ion_mz(effective_mass, adduct) - h2o - 2 * .oad_sg_mass_diff$hydrogen,
      200, "Precursor-H2O", NA, "metaboliteclass")
    if (nrow(chains) >= 2L) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
        chain_mass - 2 * .oad_sg_mass_diff$hydrogen, 50,
        paste0(chains$carbon[2L], ":", chains$double_bond[2L], "-CH4O2"), NA, "acylchain")
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(
        chain_mass + c2h3n, 100,
        paste0(chains$carbon[2L], ":", chains$double_bond[2L], "-C2H3N"), NA, "acylchain")
    }
  }
  list(peaks = peaks, mass = effective_mass, d7_balance = d7_balance)
}

.oad_sg_cer_class_peaks <- function(mass, adduct, chains) {
  h2o <- .oad_sg_water
  ch4o2 <- .oad_sg_mass(1, 4, 0, 2)
  ch2o <- .oad_sg_mass(1, 2, 0, 1)
  ch5o <- .oad_sg_mass(1, 5, 0, 1)
  c2h3n <- .oad_sg_mass(2, 3, 1)
  c2h7no <- .oad_sg_mass(2, 7, 1, 1)
  peaks <- list()
  name <- .oad_sg_adduct_name(adduct)
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass, adduct) - h2o, 500, "Precursor -H2O")
    if (nrow(chains) >= 1L) {
      sphingo_mass <- .oad_sg_chain_mass(chains$carbon[1L], chains$double_bond[1L])
      add(sphingo_mass - ch4o2 + 2 * .oad_sg_mass_diff$hydrogen, 200, "Sphingo-CH4O2", "acylchain")
      add(sphingo_mass - h2o + 2 * .oad_sg_mass_diff$hydrogen, 200, "Sphingo-H2O", "acylchain")
      add(sphingo_mass - 2 * h2o + 2 * .oad_sg_mass_diff$hydrogen, 500, "Sphingo-2H2O", "acylchain")
    }
  } else if (name == "[M+CH3COO]-") {
    add(.oad_sg_ion_mz(mass, adduct) - .oad_sg_acetate - .oad_sg_mass_diff$hydrogen,
        200, "Precursor-CH3COO")
    add(.oad_sg_ion_mz(mass, adduct) - .oad_sg_acetate - ch5o, 200, "Precursor-CH3COO-CH5O")
    if (nrow(chains) >= 2L) {
      acyl_mass <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      add(acyl_mass - 2 * .oad_sg_mass_diff$hydrogen, 50, "Acyl-2H", "acylchain")
      add(acyl_mass + c2h3n, 100, "Acyl+C2H3N", "acylchain")
    }
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct) - ch2o, 300, "Precursor-CH2O")
    add(.oad_sg_ion_mz(mass, adduct) - ch2o - 2 * .oad_sg_mass_diff$hydrogen,
        200, "Precursor-CH2O-2H")
    add(.oad_sg_ion_mz(mass, adduct) - ch2o - h2o, 200, "Precursor-CH2O-H2O")
    if (nrow(chains) >= 2L) {
      acyl_mass <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      add(acyl_mass - 2 * .oad_sg_mass_diff$hydrogen, 700, "Acyl-2H", "acylchain")
      add(acyl_mass + c2h3n, 500, "Acyl+C2H3N", "acylchain")
    }
    if (nrow(chains) >= 1L) {
      sphingo_mass <- .oad_sg_chain_mass(chains$carbon[1L], chains$double_bond[1L])
      add(sphingo_mass - c2h7no, 700, "Sphingo-C2H7NO", "acylchain")
    }
  } else {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass, adduct) - h2o, 200, "Precursor-H2O")
    add(.oad_sg_ion_mz(mass, adduct) - h2o + .oad_sg_mass_diff$hydrogen,
        70, "Precursor-H2O+H")
  }
  peaks
}

.oad_sg_ps_class_peaks <- function(mass, adduct, chains) {
  c3h8no6p <- .oad_sg_mass(3, 8, 1, 6, 1)
  c3h5no2 <- .oad_sg_mass(3, 5, 1, 2)
  h2o <- .oad_sg_water
  electron <- .oad_sg_electron
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 100, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - c3h8no6p, adduct), 500, "Precursor -C3H8NO6P")
    for (i in seq_len(nrow(chains))) {
      cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(.oad_sg_ion_mz(mass - cm + .oad_sg_mass_diff$hydrogen, adduct), 30, "-Acyl", "acylchain")
      add(.oad_sg_ion_mz(mass - cm + .oad_sg_mass_diff$hydrogen - h2o, adduct), 15, "-Acyl-H2O", "acylchain")
      add(.oad_sg_ion_mz(mass - cm + 2 * .oad_sg_mass_diff$hydrogen, adduct), 10, "-Acyl+H", "acylchain")
      add(.oad_sg_ion_mz(mass - cm + 2 * .oad_sg_mass_diff$hydrogen - h2o, adduct), 10, "-Acyl-H2O+H", "acylchain")
    }
  } else if (name == "[M+Na]+") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(mass - c3h8no6p, 100, "Precursor -C3H8NO6P")
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(mass - c3h5no2, 100, "Precursor -C3H5NO2")
    for (i in seq_len(nrow(chains))) {
      cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(cm + .oad_sg_mass_diff$oxygen + electron, 30, "Acyl FA", "acylchain")
      add(cm + .oad_sg_mass_diff$oxygen + electron + .oad_sg_mass_diff$hydrogen,
          10, "Acyl FA +H", "acylchain")
    }
  }
  peaks
}

.oad_sg_lpa_class_peaks <- function(mass, adduct) {
  c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1)
  list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 300, "Precursor", NA, "precursor"),
       .oad_sg_peak_mz(c3h6o5p, 999, "Phosphoglycerol - H2O", NA, "metaboliteclass"))
}

.oad_sg_etherpc_class_peaks <- function(mass, adduct) {
  header <- .oad_sg_mass(5, 14, 1, 4, 1)
  name <- .oad_sg_adduct_name(adduct)
  if (name == "[M+H]+") {
    return(list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", NA, "precursor"),
                .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 999, "Header", NA, "metaboliteclass")))
  }
  list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), if (name == "[M+HCO3]-") 500 else 999,
                    "Precursor", NA, "precursor"))
}

.oad_sg_etherpe_class_peaks <- function(mass, adduct, chains) {
  head <- .oad_sg_mass(2, 8, 1, 4, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(head, adduct), 100, "Header")
  } else if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 100, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - head, adduct), 999, "-C2H8NO4P", "acylchain")
    add(.oad_sg_ion_mz(head, adduct), 100, "Header")
  } else {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_dmedfahfa_class_peaks <- function(mass, adduct) {
  list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999,
                    "Precursor(DMED derv.)", NA, "precursor"))
}

.oad_sg_lps_class_peaks <- function(mass, adduct, chains) {
  d <- .oad_sg_mass_diff
  c3h8no6p <- .oad_sg_mass(3, 8, 1, 6, 1)
  c3h5no2 <- .oad_sg_mass(3, 5, 1, 2)
  c3h9o6p <- .oad_sg_mass(3, 9, 0, 6, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 999, "Precursor", "precursor")
    add(c3h8no6p, 100, "Header")
    add(mass - c3h8no6p, 500, "Precursor -C3H8NO6P")
    add(mass - .oad_sg_water, 100, "Precursor -H2O")
    add(c3h9o6p, 100, "Phosphoglycerol")
    add(c3h9o6p - .oad_sg_water, 100, "Phosphoglycerol -H2O")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        add(mass - .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]) + d$hydrogen,
            30, paste0("-", chains$carbon[i], ":", chains$double_bond[i]))
      }
    }
  } else if (name == "[M-H]-") {
    add(mass, 500, "Precursor", "precursor")
    add(mass - c3h5no2, 999, "Precursor -C3H5NO2")
    add(c3h9o6p, 100, "Phosphoglycerol")
    add(c3h9o6p - .oad_sg_water, 500, "Phosphoglycerol -H2O")
  }
  peaks
}

.oad_sg_lpsd5_class_peaks <- function(mass, adduct, chains) {
  c3h8no6p <- .oad_sg_mass(3, 8, 1, 6, 1)
  c3h5no2 <- .oad_sg_mass(3, 5, 1, 2)
  c3h4d5o6p <- .oad_sg_mass(3, 4, 0, 6, 1) + 5 * .oad_sg_deuterium
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 999, "Precursor", "precursor")
    add(c3h8no6p, 100, "Header")
    add(mass - c3h8no6p, 500, "Precursor -C3H8NO6P")
    add(mass - .oad_sg_water, 100, "Precursor -H2O")
    add(c3h4d5o6p, 100, "C3H4D5O6P")
    add(c3h4d5o6p - .oad_sg_water, 100, "C3H4D5O6P - H2O")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        add(mass - .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]) +
              .oad_sg_mass_diff$hydrogen, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]))
      }
    }
  } else if (name == "[M-H]-") {
    add(mass, 500, "Precursor", "precursor")
    add(mass - c3h5no2, 999, "Precursor -C3H5NO2")
    add(c3h4d5o6p, 100, "C3H4D5O6P")
    add(c3h4d5o6p - .oad_sg_water, 500, "C3H4D5O6P - H2O")
  }
  peaks
}

.oad_sg_lped5_class_peaks <- function(mass, adduct, chains) {
  c2h8no4p <- .oad_sg_mass(2, 8, 1, 4, 1)
  c5h6d5no5p <- .oad_sg_mass(5, 6, 1, 5, 1) + 5 * .oad_sg_deuterium
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - .oad_sg_water, adduct), 100, "Precursor -H2O")
    add(.oad_sg_ion_mz(mass - c2h8no4p, adduct), 100, "Precursor -C2H8NO4P")
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - c5h6d5no5p, adduct), 100,
        "Precursor -C5H6D5NO5P")
  }
  peaks
}

.oad_sg_lpgd5_class_peaks <- function(mass, adduct, chains) {
  c3h9o6p <- .oad_sg_mass(3, 9, 0, 6, 1)
  c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - .oad_sg_water, adduct), 100, "Precursor -H2O")
    add(.oad_sg_ion_mz(c3h9o6p, adduct), 100, "Header")
  } else if (name == "[M+NH4]+") {
    add(.oad_sg_ion_mz(mass, adduct), 100, "Precursor", "precursor")
    add(mass, 100, "[M+H]+", "precursor")
    add(mass - c3h9o6p + .oad_sg_proton, 999, "[M+H]+ -C3H9O6P")
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(c3h6o5p + .oad_sg_electron, 30, "Header")
  } else if (name == "[M+Na]+") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_lpid5_class_peaks <- function(mass, adduct, chains) {
  c6h13o9p <- .oad_sg_mass(6, 13, 0, 9, 1)
  c3h2d5o5p <- .oad_sg_mass(3, 2, 0, 5, 1) + 5 * .oad_sg_deuterium
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 50, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - .oad_sg_water, adduct), 50, "Precursor -H2O")
    add(mass - c6h13o9p, 999, "[M+H]+ -Header")
  } else if (name == "[M+NH4]+") {
    add(.oad_sg_ion_mz(mass, adduct), 50, "Precursor", "precursor")
    add(mass, 50, "[M+H]+", "precursor")
    add(mass - c6h13o9p + .oad_sg_mass_diff$hydrogen, 999, "[M+H]+ -Header")
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(c6h13o9p - .oad_sg_water, adduct), 500, "Phosphoinositol -H2O")
    add(.oad_sg_ion_mz(c3h2d5o5p, adduct), 100, "Characteristic fragment")
  }
  peaks
}

.oad_sg_pgd5_class_peaks <- function(mass, adduct, chains) {
  c3h9o6p <- .oad_sg_mass(3, 9, 0, 6, 1)
  c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1)
  nl_mass <- c3h9o6p + .oad_sg_nh3
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+NH4]+") {
    add(.oad_sg_ion_mz(mass, adduct), 100, "Precursor", "precursor")
    add(mass - c3h9o6p + .oad_sg_proton, 999, "Precursor -C3H9O6P")
    for (i in seq_len(nrow(chains))) {
      cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(.oad_sg_ion_mz(mass - nl_mass - cm + .oad_sg_mass_diff$hydrogen, adduct),
          50, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      add(.oad_sg_ion_mz(.oad_sg_mass_diff$hydrogen * 2 - nl_mass, adduct),
          20, paste0(chains$carbon[i], ":", chains$double_bond[i], "-D"), "acylchain")
    }
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(c3h6o5p + .oad_sg_electron, 30, "Header-")
    for (i in seq_len(nrow(chains))) {
      cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(cm + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 50,
          paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain")
      add(mass - cm + .oad_sg_electron, 20,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
    }
  }
  peaks
}

.oad_sg_pid5_class_peaks <- function(mass, adduct, chains) {
  c6h13o9p <- .oad_sg_mass(6, 13, 0, 9, 1)
  c6h10o5 <- .oad_sg_mass(6, 10, 0, 5)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+NH4]+") {
    nl_mass <- c6h13o9p + .oad_sg_nh3
    add(.oad_sg_ion_mz(mass, adduct), 100, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - .oad_sg_nh3, adduct), 100, "[M+H]+")
    add(mass - c6h13o9p + .oad_sg_proton, 999, "Precursor -C6H13O9P")
    for (i in seq_len(nrow(chains))) {
      cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(.oad_sg_ion_mz(mass - nl_mass - cm + .oad_sg_mass_diff$hydrogen, adduct),
          50, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
    }
  } else if (name == "[M-H]-") {
    add(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    add(.oad_sg_ion_mz(c6h13o9p - .oad_sg_water, adduct), 30, "Header-")
    for (i in seq_len(nrow(chains))) {
      cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(cm + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 50,
          paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain")
      add(mass - c6h10o5 - .oad_sg_water - cm + .oad_sg_electron, 20,
          paste0("-Header-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      add(mass - cm - .oad_sg_water + .oad_sg_electron, 20,
          paste0("-Header-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
    }
  }
  peaks
}

.oad_sg_smd9_class_peaks <- function(mass, adduct, chains) {
  c5h5d9no4p <- .oad_sg_mass(5, 5, 1, 4, 1) + 9 * .oad_sg_deuterium
  cd3 <- .oad_sg_mass(1, 0, 0, 0) + 3 * .oad_sg_deuterium
  c3d9n <- .oad_sg_mass(3, 0, 1) + 9 * .oad_sg_deuterium
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type)
  }
  if (name == "[M+Na]+") {
    add(.oad_sg_ion_mz(mass, adduct) - c3d9n, 150, "Precursor-C3D9N")
    add(.oad_sg_ion_mz(mass, adduct) - c5h5d9no4p, 150, "NL of header")
  } else if (name == "[M+H]+") {
    add(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", "precursor")
    add(.oad_sg_ion_mz(c5h5d9no4p, adduct), 999, "Header")
    if (nrow(chains) >= 1L) {
      cm <- .oad_sg_chain_mass(chains$carbon[1L], chains$double_bond[1L])
      add(cm - 2 * .oad_sg_water + 2 * .oad_sg_mass_diff$hydrogen, 50,
          paste0(chains$carbon[1L], ":", chains$double_bond[1L], "-CH4O2"), "acylchain")
    }
  } else if (name == "[M+CH3COO]-") {
    add(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", "precursor")
    add(.oad_sg_ion_mz(mass - .oad_sg_acetate - cd3, adduct), 999, "Precursor-CH3COO-CD3")
    add(c5h5d9no4p - cd3, 300, "Characteristic fragment")
    if (nrow(chains) >= 2L) {
      cm <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      add(.oad_sg_ion_mz(mass - .oad_sg_acetate - cd3 - cm + .oad_sg_mass_diff$hydrogen, adduct),
          10, paste0("NL of CD3 and ", chains$carbon[2L], ":", chains$double_bond[2L]), "acylchain")
    }
  }
  peaks
}

.oad_sg_lpe_class_peaks <- function(mass, adduct, chains) {
  c2h8no4p <- .oad_sg_mass(2, 8, 1, 4, 1)
  c5h11no5p <- .oad_sg_mass(5, 11, 1, 5, 1)
  c2h2o <- .oad_sg_mass(2, 2, 0, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 999, "Precursor", "precursor")
    add(mass - .oad_sg_water, 100, "Precursor -H2O")
    add(mass - c2h8no4p, 100, "Precursor -C2H8NO4P")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(mass - chain_mass + .oad_sg_mass_diff$hydrogen, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
        add(chain_mass + c2h2o, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], " +C2H2O"), "acylchain")
      }
    }
  } else if (name == "[M-H]-") {
    add(mass, 999, "Precursor", "precursor")
    add(mass - c5h11no5p, 100, "Precursor -C5H11NO5P")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(mass - chain_mass + .oad_sg_mass_diff$hydrogen, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
        add(mass - chain_mass + .oad_sg_mass_diff$hydrogen - .oad_sg_water, 100,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), "acylchain")
      }
    }
  } else {
    add(mass, 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_lpc_class_peaks <- function(mass, adduct, chains) {
  header <- .oad_sg_mass(5, 14, 1, 4, 1)
  c5h14no <- .oad_sg_mass(5, 13, 1, 1) + .oad_sg_proton
  ch3 <- .oad_sg_mass(1, 3)
  acetate <- .oad_sg_acetate
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 999, "Precursor", "precursor")
    add(header, 500, "Header")
    add(mass - .oad_sg_water, 100, "Precursor - H2O")
    add(c5h14no, 50, "C5H14NO")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(mass - chain_mass + .oad_sg_mass_diff$hydrogen, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
        add(mass - chain_mass + 2 * .oad_sg_mass_diff$hydrogen, 10,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], "+H"), "acylchain")
      }
    }
  } else if (name == "[M+CH3COO]-") {
    add(mass, 999, "Precursor", "precursor")
    add(mass - acetate, 300, "NL of CH3COO", "precursor")
    add(mass - ch3 - acetate, 100, "NL of CH3+CH3COO")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(chain_mass + .oad_sg_mass_diff$hydrogen, 30,
            paste0(chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      }
    }
  } else {
    add(mass, 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_pc_class_peaks <- function(mass, adduct, chains) {
  header <- .oad_sg_mass(5, 14, 1, 4, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 500, "Precursor", "precursor")
    add(header, 999, "Header")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(mass - chain_mass + .oad_sg_mass_diff$hydrogen, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
        add(mass - chain_mass + .oad_sg_mass_diff$hydrogen - .oad_sg_water, 15,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), "acylchain")
        add(mass - chain_mass + 2 * .oad_sg_mass_diff$hydrogen, 10,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], "+H"), "acylchain")
        add(mass - chain_mass + 2 * .oad_sg_mass_diff$hydrogen - .oad_sg_water, 10,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O +H"), "acylchain")
      }
    }
  } else {
    add(mass, 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_pe_class_peaks <- function(mass, adduct, chains) {
  header <- .oad_sg_mass(2, 8, 1, 4, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 500, "Precursor", "precursor")
    add(mass - header, 999, "Precursor -C2H8NO4P")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(chain_mass, 50, paste0(chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      }
    }
  } else if (name == "[M-H]-") {
    add(mass, 999, "Precursor", "precursor")
    add(header, 30, "Header-")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(chain_mass + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 30,
            paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain")
        add(mass - chain_mass + .oad_sg_electron, 30,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
        add(mass - chain_mass + .oad_sg_electron - .oad_sg_water, 15,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), "acylchain")
      }
    }
  } else {
    add(mass, 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_pg_class_peaks <- function(mass, adduct, chains) {
  c3h9o6p <- .oad_sg_mass(3, 9, 0, 6, 1)
  c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1)
  nh3 <- .oad_sg_nh3
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  add_mz <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA_character_, type)
  }
  if (name == "[M+NH4]+") {
    add(mass, 100, "Precursor", "precursor")
    add_mz(mass - c3h9o6p + .oad_sg_proton, 999, "Precursor -C3H9O6P")
    nl_mass <- c3h9o6p + nh3
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(mass - nl_mass - chain_mass + .oad_sg_mass_diff$hydrogen, 50,
            paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      }
    }
  } else if (name == "[M-H]-") {
    add(mass, 999, "Precursor", "precursor")
    add_mz(c3h6o5p + .oad_sg_electron, 30, "Header-")
    for (i in seq_len(nrow(chains))) {
      if (chains$carbon[i] != 0L) {
        chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add_mz(chain_mass + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 50,
               paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain")
        add_mz(mass - chain_mass + .oad_sg_electron, 20,
               paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      }
    }
  } else {
    add(mass, 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_lpg_class_peaks <- function(mass, adduct, chains) {
  c3h9o6p <- .oad_sg_mass(3, 9, 0, 6, 1)
  c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1)
  c2h2o <- .oad_sg_mass(2, 2, 0, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  add_mz <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+") {
    add(mass, 999, "Precursor", "precursor")
    add(mass - .oad_sg_water, 100, "Precursor -H2O")
    add(c3h9o6p, 100, "Header")
    for (i in seq_len(nrow(chains))) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(mass - chain_mass + .oad_sg_mass_diff$hydrogen, 30,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      add(chain_mass + c2h2o, 30,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i], " +C2H2O"), "acylchain")
    }
  } else if (name == "[M+NH4]+") {
    add(mass, 100, "Precursor", "precursor")
    add_mz(mass, 100, "[M+H]+", "precursor")
    add_mz(mass - c3h9o6p + .oad_sg_proton, 999, "[M+H]+ -C3H9O6P")
    for (i in seq_len(nrow(chains))) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(mass - chain_mass + .oad_sg_mass_diff$hydrogen, 50,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
    }
  } else if (name == "[M-H]-") {
    add(mass, 999, "Precursor", "precursor")
    add_mz(c3h6o5p + .oad_sg_electron, 30, "Header")
    for (i in seq_len(nrow(chains))) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add_mz(chain_mass + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 50,
             paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain")
      add_mz(mass - chain_mass + .oad_sg_electron, 20,
             paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
    }
  } else {
    add(mass, 999, "Precursor", "precursor")
  }
  peaks
}

.oad_sg_lpi_class_peaks <- function(mass, adduct, chains) {
  c6h13o9p <- .oad_sg_mass(6, 13, 0, 9, 1)
  c3h7o5p <- .oad_sg_mass(3, 7, 0, 5, 1)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(mz, adduct), intensity, comment, NA_character_, type)
  }
  add_mz <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA_character_, type)
  }
  if (name == "[M+H]+" || name == "[M+NH4]+") {
    add(mass, 50, "Precursor", "precursor")
    if (name == "[M+NH4]+") add_mz(mass, 50, "[M+H]+", "precursor")
    add_mz(mass - c6h13o9p + if (name == "[M+NH4]+") .oad_sg_proton else 0, 999, "[M+H]+ -Header")
    if (name == "[M+H]+") add(mass - .oad_sg_water, 50, "Precursor -H2O")
  } else if (name == "[M-H]-") {
    add(mass, 999, "Precursor", "precursor")
    add(c6h13o9p - .oad_sg_water, 500, "Phosphoinositol -H2O")
    add(c3h7o5p, 100, "Characteristic fragment")
    for (i in seq_len(nrow(chains))) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      add(mass - chain_mass + .oad_sg_mass_diff$hydrogen - .oad_sg_water, 90,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), "acylchain")
      add(mass - chain_mass + 2 * .oad_sg_mass_diff$hydrogen - .oad_sg_water, 30,
          paste0("-", chains$carbon[i], ":", chains$double_bond[i], "+H-H2O"), "acylchain")
    }
  } else add(mass, 999, "Precursor", "precursor")
  peaks
}

.oad_sg_pa_class_peaks <- function(mass, adduct, chains) {
  c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1)
  peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", NA, "precursor"))
  peaks[[2L]] <- .oad_sg_peak_mz(c3h6o5p, 200, "Phosphoglycerol - H2O")
  for (i in seq_len(nrow(chains))) {
    chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(chain_mass + .oad_sg_mass_diff$oxygen, 300,
                                                   paste0(chains$carbon[i], ":", chains$double_bond[i]), NA, "acylchain")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass - chain_mass + .oad_sg_mass_diff$hydrogen, adduct), 100,
                                                 paste0("-", chains$carbon[i], ":", chains$double_bond[i]), NA, "acylchain")
    peaks[[length(peaks) + 1L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass - chain_mass - .oad_sg_water + .oad_sg_mass_diff$hydrogen, adduct), 100,
                                                 paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), NA, "acylchain")
  }
  peaks
}

.oad_sg_sm_class_peaks <- function(mass, adduct, chains) {
  header <- .oad_sg_mass(5, 14, 1, 4, 1)
  c3h9n <- .oad_sg_mass(3, 9, 1)
  ch3 <- .oad_sg_mass(1, 3)
  c2h2 <- .oad_sg_mass(2, 2)
  hco3 <- .oad_sg_mass(1, 1, 0, 3)
  name <- .oad_sg_adduct_name(adduct)
  peaks <- list()
  add <- function(mz, intensity, comment, type = "metaboliteclass") {
    peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(.oad_sg_ion_mz(mz, adduct), intensity, comment, NA, type)
  }
  if (name == "[M+H]+") {
    add(mass, 999, "Precursor", "precursor"); add(header, 100, "Header")
    if (nrow(chains) >= 1L) {
      sphingo_mass <- .oad_sg_chain_mass(chains$carbon[1L], chains$double_bond[1L])
      add(sphingo_mass - 2 * .oad_sg_water + 2 * .oad_sg_mass_diff$hydrogen, 50,
          paste0(chains$carbon[1L], ":", chains$double_bond[1L], "-CH4O2"), "acylchain")
    }
  } else if (name == "[M+Na]+") {
    add(mass, 999, "Precursor", "precursor"); add(mass - header, 100, "NL of header"); add(mass - c3h9n, 100, "NL of C3H9N")
  } else if (name == "[M+CH3COO]-") {
    add(mass, 500, "Precursor", "precursor"); add(mass - .oad_sg_acetate - ch3, 999, "Precursor-CH3COO-CH3")
    add(header - ch3, 300, "Characteristic fragment")
    if (nrow(chains) >= 2L) {
      cm <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      add(mass - .oad_sg_acetate - ch3 - cm + .oad_sg_mass_diff$hydrogen, 10,
          paste0("NL of CH3 and ", chains$carbon[2L], ":", chains$double_bond[2L]), "acylchain")
    }
  } else if (name == "[M+HCO3]-") {
    add(mass, 50, "Precursor", "precursor"); add(mass - hco3 - .oad_sg_mass_diff$hydrogen, 200, "NL of HCO3")
    add(mass - hco3 - .oad_sg_mass_diff$hydrogen - c3h9n, 500, "NL of HCO3 and C3H9N")
    add(mass - hco3 - .oad_sg_mass_diff$hydrogen - c3h9n - c2h2, 999, "NL of HCO3 and C3H9N and C2H2")
    add(header - .oad_sg_mass_diff$hydrogen, 500, "C3H9N")
    if (nrow(chains) >= 2L) {
      cm <- .oad_sg_chain_mass(chains$carbon[2L], chains$double_bond[2L])
      add(mass - hco3 - cm - c2h2 - 6 * .oad_sg_mass_diff$hydrogen, 10,
          paste0("-", chains$carbon[2L], ":", chains$double_bond[2L], " - C2H8"), "acylchain")
    }
  } else add(mass, 999, "Precursor", "precursor")
  peaks
}

.oad_sg_tg_class_peaks <- function(mass, adduct, chains) {
  peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", NA, "precursor"))
  if (.oad_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks[[2L]] <- .oad_sg_peak(.oad_sg_ion_mz(mass - .oad_sg_nh3, adduct), 250, "[M+H]+", NA, "precursor")
    for (i in seq_len(nrow(chains))) {
      chain_mass <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(mass - chain_mass - .oad_sg_mass_diff$oxygen, 999,
                                                     paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), NA, "acylchain")
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak_mz(chain_mass, 250,
                                                     paste0(chains$carbon[i], ":", chains$double_bond[i]), NA, "acylchain")
    }
  }
  peaks
}

.oad_sg_labelled_class_peaks <- function(mass, class, adduct, chains) {
  base_class <- sub("_d[579]$", "", class)
  name <- .oad_sg_adduct_name(adduct)
  if (class == "PC_d5") {
    header <- .oad_sg_mass(5, 9, 1, 4, 1) + 5 * .oad_sg_deuterium
    peaks <- list(); add <- function(mz, intensity, comment, type = "metaboliteclass") { peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(.oad_sg_ion_mz(mz, adduct), intensity, comment, NA, type) }; add_mz <- function(mz, intensity, comment, type = "metaboliteclass") { peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type) }
    if (name == "[M+H]+") { add(mass, 500, "Precursor", "precursor"); add(header, 999, "Header"); for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add(mass - cm + .oad_sg_mass_diff$hydrogen, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain"); add(mass - cm + .oad_sg_mass_diff$hydrogen - .oad_sg_water, 15, paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), "acylchain"); add(mass - cm + 2 * .oad_sg_mass_diff$hydrogen, 10, paste0("-", chains$carbon[i], ":", chains$double_bond[i], "+H"), "acylchain"); add(mass - cm + 2 * .oad_sg_mass_diff$hydrogen - .oad_sg_water, 10, paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O +H"), "acylchain") } }
    else if (name %in% c("[M+HCOO]-", "[M+CH3COO]-")) { add(mass, 999, "Precursor", "precursor"); add_mz(mass - .oad_sg_mass(1, 3), 100, "[M-CH3]-"); for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add_mz(cm + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 30, paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain"); add_mz(cm + .oad_sg_mass_diff$oxygen + .oad_sg_electron + .oad_sg_mass_diff$hydrogen, 10, paste0(chains$carbon[i], ":", chains$double_bond[i], " FA +H"), "acylchain") } }
    else add(mass, 999, "Precursor", "precursor")
    return(peaks)
  }
  if (class == "PE_d5" && name == "[M-H]-") {
    header <- .oad_sg_mass(5, 11, 1, 5, 1) + 5 * .oad_sg_deuterium
    return(list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", NA, "precursor"),
                .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 30, "Header-", NA, "metaboliteclass")))
  }
  if (class == "LPC_d5") {
    header <- .oad_sg_mass(5, 14, 1, 4, 1)
    c5h14no <- .oad_sg_mass(5, 13, 1, 1) + .oad_sg_proton
    peaks <- list()
    add <- function(mz, intensity, comment, type = "metaboliteclass") {
      peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(.oad_sg_ion_mz(mz, adduct), intensity, comment, NA, type)
    }
    if (name == "[M+H]+") {
      add(mass, 999, "Precursor", "precursor"); add(header, 500, "Header"); add(mass - .oad_sg_water, 100, "Precursor - H2O"); add(c5h14no, 50, "C5H14NO")
      for (i in seq_len(nrow(chains))) {
        cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(mass - cm + .oad_sg_mass_diff$hydrogen, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
        add(mass - cm + 2 * .oad_sg_mass_diff$hydrogen, 10, paste0("-", chains$carbon[i], ":", chains$double_bond[i], "+H"), "acylchain")
      }
    } else if (name == "[M+CH3COO]-") {
      add(mass, 999, "Precursor", "precursor"); add(mass - .oad_sg_acetate, 300, "NL of CH3COO", "precursor")
      add(mass - .oad_sg_mass(1, 3) - .oad_sg_acetate, 100, "NL of CH3+CH3COO")
      for (i in seq_len(nrow(chains))) {
        cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i])
        add(cm + .oad_sg_mass_diff$hydrogen, 30, paste0(chains$carbon[i], ":", chains$double_bond[i]), "acylchain")
      }
    } else add(mass, 999, "Precursor", "precursor")
    return(peaks)
  }
  if (class == "LPE_d5") {
    c2h8no4p <- .oad_sg_mass(2, 8, 1, 4, 1)
    c5h6d5no5p <- .oad_sg_mass(5, 6, 1, 5, 1) + 5 * .oad_sg_deuterium
    c2h2o <- .oad_sg_mass(2, 2, 0, 1)
    peaks <- list(); add <- function(mz, intensity, comment, type = "metaboliteclass") {
      peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(.oad_sg_ion_mz(mz, adduct), intensity, comment, NA, type)
    }
    if (name == "[M+H]+") {
      add(mass, 999, "Precursor", "precursor"); add(mass - .oad_sg_water, 100, "Precursor -H2O"); add(mass - c2h8no4p, 100, "Precursor -C2H8NO4P")
      for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add(mass - cm + .oad_sg_mass_diff$hydrogen, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain"); add(cm + c2h2o, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i], " +C2H2O"), "acylchain") }
    } else if (name == "[M-H]-") {
      add(mass, 999, "Precursor", "precursor"); add(mass - c5h6d5no5p, 100, "Precursor -C5H6D5NO5P")
      for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add(mass - cm + .oad_sg_mass_diff$hydrogen, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain"); add(mass - cm + .oad_sg_mass_diff$hydrogen - .oad_sg_water, 100, paste0("-", chains$carbon[i], ":", chains$double_bond[i], "-H2O"), "acylchain") }
    } else add(mass, 999, "Precursor", "precursor")
    return(peaks)
  }
  if (class == "LPG_d5") {
    c3h9o6p <- .oad_sg_mass(3, 9, 0, 6, 1); c3h6o5p <- .oad_sg_mass(3, 6, 0, 5, 1); c2h2o <- .oad_sg_mass(2, 2, 0, 1)
    peaks <- list(); add <- function(mz, intensity, comment, type = "metaboliteclass") { peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(.oad_sg_ion_mz(mz, adduct), intensity, comment, NA, type) }; add_mz <- function(mz, intensity, comment, type = "metaboliteclass") { peaks[[length(peaks) + 1L]] <<- .oad_sg_peak_mz(mz, intensity, comment, NA, type) }
    if (name == "[M+H]+") { add(mass, 999, "Precursor", "precursor"); add(mass - .oad_sg_water, 100, "Precursor -H2O"); add(c3h9o6p, 100, "Header"); for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add(mass - cm + .oad_sg_mass_diff$hydrogen, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain"); add(cm + c2h2o, 30, paste0("-", chains$carbon[i], ":", chains$double_bond[i], " +C2H2O"), "acylchain") } }
    else if (name == "[M+NH4]+") { add(mass, 100, "Precursor", "precursor"); add_mz(mass, 100, "[M+H]+", "precursor"); add_mz(mass - c3h9o6p + .oad_sg_proton, 999, "[M+H]+ -C3H9O6P"); for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add(mass - cm + .oad_sg_mass_diff$hydrogen, 50, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain") } }
    else if (name == "[M-H]-") { add(mass, 999, "Precursor", "precursor"); add_mz(c3h6o5p + .oad_sg_electron, 30, "Header-"); for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add_mz(cm + .oad_sg_mass_diff$oxygen + .oad_sg_electron, 50, paste0(chains$carbon[i], ":", chains$double_bond[i], " FA"), "acylchain"); add_mz(mass - cm + .oad_sg_electron, 20, paste0("-", chains$carbon[i], ":", chains$double_bond[i]), "acylchain") } }
    else add(mass, 999, "Precursor", "precursor")
    return(peaks)
  }
  if (class == "LPI_d5") {
    peaks <- .oad_sg_lpi_class_peaks(mass, adduct, chains)
    if (name == "[M-H]-") { peaks[[3L]]$Mass <- .oad_sg_ion_mz(.oad_sg_mass(3, 2, 0, 5, 1) + 5 * .oad_sg_deuterium, adduct); peaks[[3L]]$Comment <- "Characteristic fragment" }
    return(peaks)
  }
  if (class == "SM_d9" && name == "[M+H]+") {
    header <- .oad_sg_mass(5, 5, 1, 4, 1) + 9 * .oad_sg_deuterium
    return(list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", NA, "precursor"),
                .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 999, "Header", NA, "metaboliteclass")))
  }
  if (class == "SM_d9") {
    header <- .oad_sg_mass(5, 5, 1, 4, 1) + 9 * .oad_sg_deuterium
    cd3 <- .oad_sg_mass(1, 0, 0, 0) + 3 * .oad_sg_deuterium
    c3d9n <- .oad_sg_mass(3, 0, 1) + 9 * .oad_sg_deuterium
    peaks <- list(); add <- function(mz, intensity, comment, type = "metaboliteclass") { peaks[[length(peaks) + 1L]] <<- .oad_sg_peak(.oad_sg_ion_mz(mz, adduct), intensity, comment, NA, type) }
    if (name == "[M+Na]+") { add(mass, 150, "Precursor", "precursor"); add(mass - c3d9n, 150, "Precursor-C3D9N"); add(mass - header, 150, "NL of header") }
    else if (name == "[M+CH3COO]-") { add(mass, 500, "Precursor", "precursor"); add(mass - .oad_sg_acetate - cd3, 999, "Precursor-CH3COO-CD3"); add(header - cd3, 300, "Characteristic fragment"); for (i in seq_len(nrow(chains))) { cm <- .oad_sg_chain_mass(chains$carbon[i], chains$double_bond[i]); add(mass - .oad_sg_acetate - cd3 - cm + .oad_sg_mass_diff$hydrogen, 10, paste0("NL of CD3 and ", chains$carbon[i], ":", chains$double_bond[i]), "acylchain") } }
    return(peaks)
  }
  .oad_sg_class_peaks(mass, base_class, adduct, chains)
}

.oad_sg_double_bond_peaks <- function(lipid_mass, chain, adduct, nl_mass = 0,
                                      abundance = 40, oad_ids = NULL) {
  positions <- chain$positions
  if (is.list(positions)) positions <- positions[[1L]]
  if (chain$carbon <= 0L || length(positions) == 0L) return(list())
  d <- .oad_sg_mass_diff
  ch2 <- d$carbon + 2 * d$hydrogen
  diffs <- rep(ch2, chain$carbon)
  for (position in positions) {
    if (position >= 1L && position < chain$carbon) {
      diffs[position] <- diffs[position] - d$hydrogen
      diffs[position + 1L] <- diffs[position + 1L] - d$hydrogen
    }
  }
  diffs <- cumsum(diffs)
  chain_loss <- lipid_mass - .oad_sg_chain_mass(chain$carbon, chain$double_bond) - nl_mass
  ids <- c("OAD01", "OAD02", "OAD02+O", "OAD03", "OAD04", "OAD05", "OAD06",
           "OAD07", "OAD08", "OAD09", "OAD10", "OAD11", "OAD12", "OAD13",
           "OAD14", "OAD15", "OAD15+O", "OAD16", "OAD17", "OAD12+O",
           "OAD12+O+H", "OAD12+O+2H", "OAD01+H")
  if (!is.null(oad_ids)) ids <- intersect(ids, oad_ids)
  result <- list()
  add <- function(id, neutral_mass, factor, suffix) {
    if (id %in% ids) result[[length(result) + 1L]] <<- .oad_sg_peak(
      .oad_sg_ion_mz(neutral_mass, adduct), abundance * factor,
      paste0(chain$carbon, ":", chain$double_bond, " C", bond, " ", suffix, " ", id), id)
  }
  for (bond in positions) {
    if (bond <= 1L || bond > chain$carbon - 1L) next
    higher <- chain_loss + diffs[bond + 1L] + d$hydrogen + d$oxygen
    middle <- chain_loss + diffs[bond]
    lower <- chain_loss + diffs[bond - 1L] + d$hydrogen
    add("OAD01", higher + d$hydrogen, 0.3, "+C +O +H")
    add("OAD02", higher, 1, "+O")
    add("OAD02+O", higher + d$oxygen, 0.3, "+O")
    add("OAD03", higher - d$hydrogen, 0.5, "+C +O -H")
    add("OAD04", higher - 2 * d$hydrogen, 0.2, "+C +O -2H")
    add("OAD05", higher - .oad_sg_water + d$hydrogen, 0.15, "+C +O +H -H2O")
    add("OAD06", higher - .oad_sg_water, 0.3, "+C +O -H2O")
    add("OAD07", higher - .oad_sg_water - d$hydrogen, 0.1, "+C +O -H -H2O")
    add("OAD08", middle + d$oxygen, 0.2, "+O")
    add("OAD09", middle + 2 * d$hydrogen, 0.2, "")
    add("OAD10", middle + d$hydrogen, 0.3, "-H")
    add("OAD11", middle, 0.1, "-2H")
    add("OAD12", lower + d$oxygen - d$hydrogen, 0.3, "+O -H")
    add("OAD13", lower + d$oxygen - 2 * d$hydrogen, 0.3, "+O -2H")
    add("OAD14", lower + d$hydrogen, 0.2, "-C +H")
    add("OAD15", lower, 0.8, "-C")
    add("OAD15+O", lower + d$oxygen, 0.2, "-C")
    add("OAD16", lower - d$hydrogen, 0.2, "-C -H")
    add("OAD17", lower - 2 * d$hydrogen, 0.15, "-C -2H")
    add("OAD12+O", lower + 2 * d$oxygen - d$hydrogen, 0.3, "+O -H")
    add("OAD12+O+H", lower + 2 * d$oxygen, 0.3, "+O -H")
    add("OAD12+O+2H", lower + 2 * d$oxygen + d$hydrogen, 0.2, "+O -H")
    add("OAD01+H", higher + 2 * d$hydrogen, 0.2, "+C +O +H")
  }
  result
}

.oad_sg_generate_class <- function(lipid, class, adduct, oad_ids = NULL,
                                    abundance = 40, nl_mass = 0) {
  name <- .oad_sg_field(lipid, "Name", class)
  mass <- as.numeric(.oad_sg_field(lipid, c("Mass", "MolecularMass", "NeutralMass"), NA_real_))
  if (is.na(mass)) stop("lipid must contain a neutral Mass for OAD generation")
  chains <- .oad_sg_parse_chains(name, lipid)
  cer_ns_d7 <- NULL
  peaks <- if (class == "Cer_NS_d7") {
    cer_ns_d7 <- .oad_sg_cerns_d7_class_peaks(mass, adduct, chains)
    mass <- cer_ns_d7$mass
    cer_ns_d7$peaks
  } else if (class == "CE_d7") {
    .oad_sg_ced7_class_peaks(mass, adduct)
  } else if (class %in% c("Cer_NS", "Cer_NDS")) {
    .oad_sg_cer_class_peaks(mass, adduct, chains)
  } else if (class %in% c("PS", "PS_d5")) {
    .oad_sg_ps_class_peaks(mass, adduct, chains)
  } else if (class == "LPA") {
    .oad_sg_lpa_class_peaks(mass, adduct)
  } else if (class == "EtherPC") {
    .oad_sg_etherpc_class_peaks(mass, adduct)
  } else if (class == "EtherPE") {
    .oad_sg_etherpe_class_peaks(mass, adduct, chains)
  } else if (class == "DMEDFAHFA") {
    .oad_sg_dmedfahfa_class_peaks(mass, adduct)
  } else if (class == "DG") {
    .oad_sg_dg_class_peaks(mass, adduct, chains)
  } else if (class == "DG_d5") {
    .oad_sg_dg_class_peaks(mass, adduct, chains, labelled = TRUE)
  } else if (class == "TG_d5") {
    .oad_sg_tgd5_class_peaks(mass, adduct, chains)
  } else if (class == "LPS") {
    .oad_sg_lps_class_peaks(mass, adduct, chains)
  } else if (class %in% c("PI", "PI_d5")) {
    .oad_sg_pid5_class_peaks(mass, adduct, chains)
  } else if (class == "PG_d5") {
    .oad_sg_pgd5_class_peaks(mass, adduct, chains)
  } else if (class == "LPS_d5") {
    .oad_sg_lpsd5_class_peaks(mass, adduct, chains)
  } else if (class == "LPE_d5") {
    .oad_sg_lped5_class_peaks(mass, adduct, chains)
  } else if (class == "LPG_d5") {
    .oad_sg_lpgd5_class_peaks(mass, adduct, chains)
  } else if (class == "LPI_d5") {
    .oad_sg_lpid5_class_peaks(mass, adduct, chains)
  } else if (class == "LPC") {
    .oad_sg_lpc_class_peaks(mass, adduct, chains)
  } else if (class == "PC") {
    .oad_sg_pc_class_peaks(mass, adduct, chains)
  } else if (class == "PE") {
    .oad_sg_pe_class_peaks(mass, adduct, chains)
  } else if (class == "PG") {
    .oad_sg_pg_class_peaks(mass, adduct, chains)
  } else if (class == "LPG") {
    .oad_sg_lpg_class_peaks(mass, adduct, chains)
  } else if (class == "LPI") {
    .oad_sg_lpi_class_peaks(mass, adduct, chains)
  } else if (class == "PA") {
    .oad_sg_pa_class_peaks(mass, adduct, chains)
  } else if (class == "SM") {
    .oad_sg_sm_class_peaks(mass, adduct, chains)
  } else if (class == "SM_d9") {
    .oad_sg_smd9_class_peaks(mass, adduct, chains)
  } else if (class == "TG") {
    .oad_sg_tg_class_peaks(mass, adduct, chains)
  } else if (class == "LPE") {
    .oad_sg_lpe_class_peaks(mass, adduct, chains)
  } else if (grepl("_d[579]$", class)) {
    .oad_sg_labelled_class_peaks(mass, class, adduct, chains)
  } else {
    .oad_sg_class_peaks(mass, class, adduct, chains)
  }
  for (i in seq_len(nrow(chains))) {
    chain_nl_mass <- nl_mass
    if (class %in% c("DG", "DG_d5", "TG_d5") && nrow(chains) > 1L) {
      next_index <- if (i == nrow(chains)) 1L else i + 1L
      chain_nl_mass <- .oad_sg_chain_mass(chains$carbon[next_index], chains$double_bond[next_index])
      if (class == "TG_d5") {
        chain_nl_mass <- chain_nl_mass + .oad_sg_mass(1, 0, 0, 1) + 3 * .oad_sg_deuterium
      }
    }
    if (class == "Cer_NS_d7") {
      chain_nl_mass <- chain_nl_mass + cer_ns_d7$d7_balance
    }
    peaks <- c(peaks, .oad_sg_double_bond_peaks(mass, chains[i, , drop = FALSE], adduct,
                                                 chain_nl_mass, abundance, oad_ids))
  }
    adduct_name <- .oad_sg_adduct_name(adduct)
    ion_mode <- if (grepl("-$", adduct_name)) "Negative" else "Positive"
    charge <- if (grepl("2[+-]$", adduct_name)) 2L else 1L
    list(Name = name, PrecursorMz = .oad_sg_ion_mz(mass, adduct),
       Mz = .oad_sg_ion_mz(mass, adduct), AdductIonName = .oad_sg_adduct_name(adduct),
      AdductType = adduct, IonMode = ion_mode, Charge = charge,
      CompoundClass = class, Ontology = .oad_sg_field(lipid, "Ontology"),
      Spectrum = .oad_sg_merge(peaks), Formula = .oad_sg_field(lipid, "Formula"),
      InChIKey = .oad_sg_field(lipid, "InChIKey"), SMILES = .oad_sg_field(lipid, "SMILES"))
}

# C#-style class generator entry points.
generate_pc_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PC", adduct, c("OAD01", "OAD02", "OAD03", "OAD14", "OAD15", "OAD16", "OAD12+O", "OAD12+O+H", "OAD12+O+2H"), 40)
generate_lpc_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPC", adduct, c("OAD01", "OAD02", "OAD02+O", "OAD03", "OAD04", "OAD08", "OAD13", "OAD15", "OAD16", "OAD17", "OAD12+O", "OAD12+O+H", "OAD01+H"), 40)
generate_pe_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PE", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD09", "OAD10", "OAD11", "OAD12", "OAD13", "OAD14", "OAD15", "OAD16", "OAD17", "OAD01+H"), 30)
generate_lpe_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPE", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD06", "OAD08", "OAD13", "OAD14", "OAD15", "OAD16", "OAD17", "OAD12+O"), 30)
generate_pg_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PG", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD08", "OAD10", "OAD11", "OAD12", "OAD13", "OAD14", "OAD15", "OAD16", "OAD17", "OAD01+H"), 30, ifelse(.oad_sg_adduct_name(adduct) == "[M+NH4]+", .oad_sg_mass(3, 9, oxygen = 6, phosphorus = 1) + 18.033823, 0))
generate_lpg_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPG", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16", "OAD17", "OAD12+O", "OAD12+O+H", "OAD01+H"), 30)
generate_pi_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PI", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16", "OAD12+O", "OAD12+O+H", "OAD12+O+2H"), 30)
generate_lpi_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPI", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD08", "OAD13", "OAD15", "OAD16", "OAD17", "OAD12+O", "OAD01+H"), 40)
generate_ps_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PS", adduct, c("OAD01", "OAD02", "OAD02+O", "OAD03", "OAD04", "OAD08", "OAD12", "OAD14", "OAD15", "OAD15+O", "OAD16", "OAD17", "OAD12+O", "OAD12+O+H", "OAD12+O+2H", "OAD01+H"), 40)
generate_lps_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPS", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16", "OAD17", "OAD12+O", "OAD12+O+H", "OAD01+H"), 30)
generate_lpsd5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPS_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16", "OAD17", "OAD12+O", "OAD12+O+H", "OAD01+H"), 30)
generate_pc_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PC_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD14", "OAD15", "OAD16", "OAD12+O", "OAD12+O+H", "OAD12+O+2H"), 40)
generate_pe_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PE_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD09", "OAD10", "OAD11", "OAD12", "OAD13", "OAD14", "OAD15", "OAD16", "OAD17", "OAD01+H"), 30)
generate_pg_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PG_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD08", "OAD10", "OAD11", "OAD12", "OAD13", "OAD14", "OAD15", "OAD16", "OAD17", "OAD01+H"), 30)
generate_pi_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PI_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16", "OAD12+O", "OAD12+O+H", "OAD12+O+2H"), 30)
generate_ps_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PS_d5", adduct, c("OAD01", "OAD02", "OAD02+O", "OAD03", "OAD04", "OAD08", "OAD12", "OAD14", "OAD15", "OAD15+O", "OAD16", "OAD17", "OAD12+O", "OAD12+O+H", "OAD12+O+2H", "OAD01+H"), 40)
generate_lpc_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPC_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD14", "OAD15", "OAD16"), 40)
generate_lpe_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPE_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16"), 30)
generate_lpg_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPG_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16"), 30)
generate_lpi_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPI_d5", adduct, c("OAD01", "OAD02", "OAD03", "OAD04", "OAD14", "OAD15", "OAD16"), 30)
generate_dg_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "DG_d5", adduct, NULL, 40)
generate_tg_d5_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "TG_d5", adduct, NULL, 40)
generate_sm_d9_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "SM_d9", adduct, c("OAD02", "OAD03", "OAD05", "OAD06", "OAD07", "OAD08", "OAD15", "OAD16"), 40)
generate_cer_ns_d7_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "Cer_NS_d7", adduct, c("OAD02", "OAD03", "OAD05", "OAD06", "OAD07", "OAD15", "OAD16"), 40)
generate_ce_d7_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "CE_d7", adduct, NULL, 40)
generate_pa_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "PA", adduct, c("OAD01", "OAD02", "OAD03", "OAD14", "OAD15", "OAD16", "OAD12+O", "OAD12+O+H", "OAD12+O+2H"), 30)
generate_dg_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "DG", adduct, NULL, 40)
generate_tg_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "TG", adduct, NULL, 40)
generate_sm_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "SM", adduct, c("OAD02", "OAD03", "OAD05", "OAD06", "OAD07", "OAD08", "OAD15", "OAD16"), 40)
generate_ceramide_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, .oad_sg_extract_class(.oad_sg_field(lipid, "Name", "Cer_NS")), adduct, c("OAD02", "OAD03", "OAD05", "OAD06", "OAD07", "OAD15", "OAD16"), 40)
generate_dmedfahfa_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "DMEDFAHFA", adduct, NULL, 40)
generate_oad_default_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, .oad_sg_extract_class(.oad_sg_field(lipid, "Name", "Unknown")), adduct, NULL, 40)

generate_etherpc_oad_spectrum <- generate_pc_oad_spectrum
generate_etherpe_oad_spectrum <- generate_pe_oad_spectrum
generate_etherlyso_oad_spectrum <- generate_oad_default_spectrum
generate_lpa_oad_spectrum <- function(lipid, adduct, molecule = NULL) .oad_sg_generate_class(lipid, "LPA", adduct, c("OAD01", "OAD02", "OAD03", "OAD14", "OAD15", "OAD16"), 30)

generate_oad_lipid_spectrum <- function(lipid, adduct = NULL, molecule = NULL) {
  if (is.null(lipid) || is.null(lipid$Name)) stop("lipid must contain Name")
  class <- .oad_sg_extract_class(lipid$Name)
  adduct <- .oad_sg_or(adduct, .oad_sg_adduct_name(NULL, lipid))
  generator <- switch(class,
    PC = generate_pc_oad_spectrum, LPC = generate_lpc_oad_spectrum,
    EtherPC = generate_etherpc_oad_spectrum, EtherLPC = generate_etherlyso_oad_spectrum,
    PE = generate_pe_oad_spectrum, LPE = generate_lpe_oad_spectrum,
    EtherPE = generate_etherpe_oad_spectrum, EtherLPE = generate_etherlyso_oad_spectrum,
    PG = generate_pg_oad_spectrum, LPG = generate_lpg_oad_spectrum,
    PI = generate_pi_oad_spectrum, LPI = generate_lpi_oad_spectrum,
    PS = generate_ps_oad_spectrum, LPS = generate_lps_oad_spectrum,
    PA = generate_pa_oad_spectrum, LPA = generate_lpa_oad_spectrum,
    DG = generate_dg_oad_spectrum, TG = generate_tg_oad_spectrum,
    SM = generate_sm_oad_spectrum,
    Cer_NS = generate_ceramide_oad_spectrum, Cer_NDS = generate_ceramide_oad_spectrum,
    DMEDFAHFA = generate_dmedfahfa_oad_spectrum,
    PC_d5 = generate_pc_d5_oad_spectrum, PE_d5 = generate_pe_d5_oad_spectrum,
    PG_d5 = generate_pg_d5_oad_spectrum, PI_d5 = generate_pi_d5_oad_spectrum,
    PS_d5 = generate_ps_d5_oad_spectrum, LPC_d5 = generate_lpc_d5_oad_spectrum,
    LPE_d5 = generate_lpe_d5_oad_spectrum, LPG_d5 = generate_lpg_d5_oad_spectrum,
    LPI_d5 = generate_lpi_d5_oad_spectrum, LPS_d5 = generate_lpsd5_oad_spectrum,
    DG_d5 = generate_dg_d5_oad_spectrum, TG_d5 = generate_tg_d5_oad_spectrum,
    SM_d9 = generate_sm_d9_oad_spectrum, Cer_NS_d7 = generate_cer_ns_d7_oad_spectrum,
    CE_d7 = generate_ce_d7_oad_spectrum,
    generate_oad_default_spectrum
  )
  generator(lipid, adduct, molecule)
}

get_oad_lipid_spectrum <- generate_oad_lipid_spectrum
get_oad_based_lipid_spectrum <- generate_oad_lipid_spectrum
can_generate_oad_lipid_spectrum <- function(lipid, adduct = NULL) {
  if (is.null(lipid) || is.null(lipid$Name)) return(FALSE)
  class <- .oad_sg_extract_class(lipid$Name)
  name <- .oad_sg_adduct_name(adduct, lipid)
  allowed <- switch(class,
    PC = c("[M+H]+", "[M+Na]+", "[M+HCOO]-", "[M+CH3COO]-", "[M+HCO3]-"),
    LPC = c("[M+H]+", "[M+CH3COO]-"), LPE = c("[M+H]+", "[M-H]-"),
    PG = c("[M+NH4]+", "[M-H]-"), PI = c("[M+NH4]+", "[M-H]-"),
    PE = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    LPS = c("[M+H]+", "[M-H]-"),
    EtherPC = c("[M+H]+", "[M+HCOO]-", "[M+CH3COO]-", "[M+HCO3]-"),
    EtherPE = c("[M+H]+", "[M-H]-"), LPA = c("[M-H]-"),
    DMEDFAHFA = c("[M+H]+"),
    LPG = c("[M+NH4]+", "[M-H]-"), LPI = c("[M+NH4]+", "[M-H]-"),
    PS = c("[M+H]+", "[M+Na]+", "[M-H]-"), LPS = c("[M+H]+", "[M-H]-"),
    PA = c("[M+NH4]+", "[M-H]-"), SM = c("[M+H]+", "[M+Na]+", "[M+CH3COO]-"),
    DG = c("[M+NH4]+"), TG = c("[M+NH4]+"),
    PC_d5 = c("[M+H]+", "[M+Na]+", "[M+HCOO]-", "[M+CH3COO]-"),
    PE_d5 = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    PG_d5 = c("[M+NH4]+", "[M-H]-"), PI_d5 = c("[M+NH4]+", "[M-H]-"),
    PS_d5 = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    LPC_d5 = c("[M+H]+", "[M+CH3COO]-"), LPE_d5 = c("[M+H]+", "[M-H]-"),
    LPG_d5 = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M-H]-"),
    LPI_d5 = c("[M+H]+", "[M+NH4]+", "[M-H]-"),
    LPS_d5 = c("[M+H]+", "[M-H]-"), DG_d5 = c("[M+NH4]+"),
    TG_d5 = c("[M+NH4]+"), SM_d9 = c("[M+H]+", "[M+Na]+", "[M+CH3COO]-"),
    Cer_NS = c("[M+H]+", "[M+Na]+", "[M+CH3COO]-", "[M-H]-"),
    Cer_NDS = c("[M+H]+", "[M+Na]+", "[M+CH3COO]-", "[M-H]-"),
    Cer_NS_d7 = c("[M+H]+", "[M-H]-", "[M+CH3COO]-"), CE_d7 = c("[M+NH4]+", "[M+Na]+", "[M-H]-"),
    c("[M+H]+", "[M+Na]+", "[M-H]-"))
  name %in% allowed
}

# C#-style public aliases retained for callers ported from MSDIAL.
OadLipidSpectrumGenerator <- generate_oad_lipid_spectrum
CanGenerateOadSpectrum <- can_generate_oad_lipid_spectrum
OadDefaultSpectrumGenerator <- generate_oad_default_spectrum
EtherPCOadSpectrumGenerator <- generate_etherpc_oad_spectrum
EtherPEOadSpectrumGenerator <- generate_etherpe_oad_spectrum
EtherLPCOadSpectrumGenerator <- generate_etherlyso_oad_spectrum
EtherLPEOadSpectrumGenerator <- generate_etherlyso_oad_spectrum
PCOadSpectrumGenerator <- generate_pc_oad_spectrum
LPCOadSpectrumGenerator <- generate_lpc_oad_spectrum
PEOadSpectrumGenerator <- generate_pe_oad_spectrum
LPEOadSpectrumGenerator <- generate_lpe_oad_spectrum
PGOadSpectrumGenerator <- generate_pg_oad_spectrum
LPGOadSpectrumGenerator <- generate_lpg_oad_spectrum
PIOadSpectrumGenerator <- generate_pi_oad_spectrum
LPIOadSpectrumGenerator <- generate_lpi_oad_spectrum
PSOadSpectrumGenerator <- generate_ps_oad_spectrum
LPSOadSpectrumGenerator <- generate_lps_oad_spectrum
LPSd5OadSpectrumGenerator <- generate_lpsd5_oad_spectrum
PCd5OadSpectrumGenerator <- generate_pc_d5_oad_spectrum
PEd5OadSpectrumGenerator <- generate_pe_d5_oad_spectrum
PGd5OadSpectrumGenerator <- generate_pg_d5_oad_spectrum
PId5OadSpectrumGenerator <- generate_pi_d5_oad_spectrum
PSd5OadSpectrumGenerator <- generate_ps_d5_oad_spectrum
LPCd5OadSpectrumGenerator <- generate_lpc_d5_oad_spectrum
LPEd5OadSpectrumGenerator <- generate_lpe_d5_oad_spectrum
LPGd5OadSpectrumGenerator <- generate_lpg_d5_oad_spectrum
LPId5OadSpectrumGenerator <- generate_lpi_d5_oad_spectrum
DGd5OadSpectrumGenerator <- generate_dg_d5_oad_spectrum
TGd5OadSpectrumGenerator <- generate_tg_d5_oad_spectrum
SMd9OadSpectrumGenerator <- generate_sm_d9_oad_spectrum
CerNSd7OadSpectrumGenerator <- generate_cer_ns_d7_oad_spectrum
CEd7OadSpectrumGenerator <- generate_ce_d7_oad_spectrum
PAOadSpectrumGenerator <- generate_pa_oad_spectrum
LPAOadSpectrumGenerator <- generate_lpa_oad_spectrum
TGOadSpectrumGenerator <- generate_tg_oad_spectrum
DGOadSpectrumGenerator <- generate_dg_oad_spectrum
SMOadSpectrumGenerator <- generate_sm_oad_spectrum
CeramideOadSpectrumGenerator <- generate_ceramide_oad_spectrum
DMEDFAHFAOadSpectrumGenerator <- generate_dmedfahfa_oad_spectrum

OadSpectrumPeakGenerator <- function(lipid_mass, chain, adduct, nl_mass = 0,
                                     abundance = 40, oad_id = NULL) {
  .oad_sg_merge(.oad_sg_double_bond_peaks(lipid_mass, chain, adduct, nl_mass,
                                           abundance, oad_id))
}
GetAcylDoubleBondSpectrum <- OadSpectrumPeakGenerator
GetAlkylDoubleBondSpectrum <- OadSpectrumPeakGenerator
GetSphingoDoubleBondSpectrum <- OadSpectrumPeakGenerator
