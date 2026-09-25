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
        800, paste0(chains$carbon[i], ":", chains$double_bond[i], " loss"),
        NA_character_, "acylchain")
    }
  }
  peaks
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
  nh3 <- .oad_sg_nh4
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

.oad_sg_labelled_class_peaks <- function(mass, class, adduct, chains) {
  base_class <- sub("_d[579]$", "", class)
  name <- .oad_sg_adduct_name(adduct)
  if (class == "PC_d5") {
    header <- .oad_sg_mass(5, 9, 1, 4, 1) + 5 * .oad_sg_deuterium
    peaks <- list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", NA, "precursor"),
                  .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 999, "Header", NA, "metaboliteclass"))
    if (name %in% c("[M+H]+", "[M+Na]+")) {
      peaks[[length(peaks) + 1L]] <- .oad_sg_peak(
        .oad_sg_ion_mz(mass - .oad_sg_water, adduct), 100,
        "Precursor -H2O", NA, "metaboliteclass")
    }
    return(peaks)
  }
  if (class == "PE_d5" && name == "[M-H]-") {
    header <- .oad_sg_mass(5, 11, 1, 5, 1) + 5 * .oad_sg_deuterium
    return(list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 999, "Precursor", NA, "precursor"),
                .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 30, "Header-", NA, "metaboliteclass")))
  }
  if (class == "SM_d9" && name == "[M+H]+") {
    header <- .oad_sg_mass(5, 5, 1, 4, 1) + 9 * .oad_sg_deuterium
    return(list(.oad_sg_peak(.oad_sg_ion_mz(mass, adduct), 500, "Precursor", NA, "precursor"),
                .oad_sg_peak(.oad_sg_ion_mz(header, adduct), 999, "Header", NA, "metaboliteclass")))
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
  peaks <- if (class == "LPS") {
    .oad_sg_lps_class_peaks(mass, adduct, chains)
  } else if (class == "LPS_d5") {
    .oad_sg_lpsd5_class_peaks(mass, adduct, chains)
  } else if (class == "LPC") {
    .oad_sg_lpc_class_peaks(mass, adduct, chains)
  } else if (class == "PC") {
    .oad_sg_pc_class_peaks(mass, adduct, chains)
  } else if (class == "PE") {
    .oad_sg_pe_class_peaks(mass, adduct, chains)
  } else if (class == "PG") {
    .oad_sg_pg_class_peaks(mass, adduct, chains)
  } else if (class == "LPE") {
    .oad_sg_lpe_class_peaks(mass, adduct, chains)
  } else if (grepl("_d[579]$", class)) {
    .oad_sg_labelled_class_peaks(mass, class, adduct, chains)
  } else {
    .oad_sg_class_peaks(mass, class, adduct, chains)
  }
  for (i in seq_len(nrow(chains))) {
    peaks <- c(peaks, .oad_sg_double_bond_peaks(mass, chains[i, , drop = FALSE], adduct,
                                                 nl_mass, abundance, oad_ids))
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
    PE = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    LPS = c("[M+H]+", "[M-H]-"),
    DG = c("[M+NH4]+"), TG = c("[M+NH4]+"),
    PC_d5 = c("[M+H]+", "[M+Na]+", "[M+HCOO]-", "[M+CH3COO]-"),
    PE_d5 = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    PG_d5 = c("[M+NH4]+", "[M-H]-"), PI_d5 = c("[M+NH4]+", "[M-H]-"),
    PS_d5 = c("[M+H]+", "[M+Na]+", "[M-H]-"),
    LPC_d5 = c("[M+H]+", "[M+CH3COO]-"), LPE_d5 = c("[M+H]+", "[M-H]-"),
    LPG_d5 = c("[M+NH4]+", "[M-H]-"), LPI_d5 = c("[M+NH4]+", "[M-H]-"),
    LPS_d5 = c("[M+H]+", "[M-H]-"), DG_d5 = c("[M+NH4]+"),
    TG_d5 = c("[M+NH4]+"), SM_d9 = c("[M+H]+", "[M+Na]+", "[M+CH3COO]-"),
    Cer_NS_d7 = c("[M+H]+", "[M-H]-"), CE_d7 = c("[M+NH4]+", "[M-H]-"),
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
