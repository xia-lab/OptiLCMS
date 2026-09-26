# EID lipid spectrum generation translated from MSDIAL's EidSpectrumGenerator.

.eid_sg_mass <- function(carbon = 0, hydrogen = 0, nitrogen = 0, oxygen = 0,
                         phosphorus = 0, sulfur = 0) {
  carbon * 12 + hydrogen * 1.00782503223 + nitrogen * 14.00307400443 +
    oxygen * 15.99491461957 + phosphorus * 30.97376199842 +
    sulfur * 31.9720711744
}
.eid_sg_proton <- 1.007276466621
.eid_sg_electron <- 0.00054858026
.eid_sg_h2o <- .eid_sg_mass(hydrogen = 2, oxygen = 1)
.eid_sg_nh3 <- .eid_sg_mass(hydrogen = 3, nitrogen = 1)
.eid_sg_ch2 <- .eid_sg_mass(carbon = 1, hydrogen = 2)
.eid_sg_o <- .eid_sg_mass(oxygen = 1)
.eid_sg_h <- .eid_sg_mass(hydrogen = 1)
.eid_sg_header <- list(
  PC = .eid_sg_mass(5, 14, 1, 4, 1),
  PC_FRAGMENT = .eid_sg_mass(5, 13, 1, 1) + .eid_sg_proton,
  PE = .eid_sg_mass(2, 8, 1, 4, 1),
  PG = .eid_sg_mass(3, 9, 0, 6, 1),
  PS = .eid_sg_mass(3, 8, 1, 6, 1),
  PI = .eid_sg_mass(6, 13, 0, 9, 1),
  SM = .eid_sg_mass(5, 14, 1, 4, 1),
  PC_GLY_C = .eid_sg_mass(8, 18, 1, 4, 1),
  PC_GLY_O = .eid_sg_mass(7, 16, 1, 5, 1),
  PE_GLY_C = .eid_sg_mass(5, 12, 1, 4, 1),
  PE_GLY_O = .eid_sg_mass(4, 10, 1, 5, 1),
  PG_GLY_C = .eid_sg_mass(6, 13, 0, 6, 1),
  PG_GLY_O = .eid_sg_mass(5, 11, 0, 7, 1),
  PS_GLY_C = .eid_sg_mass(6, 12, 1, 6, 1),
  PS_GLY_O = .eid_sg_mass(5, 10, 1, 7, 1),
  PI_GLY_C = .eid_sg_mass(8, 18, 0, 9, 1),
  PI_GLY_O = .eid_sg_mass(7, 16, 0, 10, 1),
  PA_H3PO4 = .eid_sg_mass(0, 3, 0, 4, 1),
  PA_GLY_C = .eid_sg_mass(3, 7, 0, 4, 1),
  PA_GLY_O = .eid_sg_mass(2, 5, 0, 5, 1),
  PG_NEUTRAL = .eid_sg_mass(3, 6, 0, 2),
  PS_CHO2 = .eid_sg_mass(1, 1, 0, 2),
  PI_C6H10O5 = .eid_sg_mass(6, 10, 0, 5)
)
.eid_sg_sm_c3h9n <- .eid_sg_mass(3, 9, 1)
.eid_sg_sm_c2h3no <- .eid_sg_mass(2, 3, 1, 1)
.eid_sg_sm_c2h2n <- .eid_sg_mass(2, 2, 1)
.eid_sg_cho2 <- .eid_sg_mass(1, 1, 0, 2)
.eid_sg_c3h6o2 <- .eid_sg_mass(3, 6, 0, 2)
.eid_sg_c3h3o2 <- .eid_sg_mass(3, 3, 0, 2)

.eid_sg_field <- function(x, names, default = NULL) {
  if (is.null(x)) return(default)
  if (!is.list(x) && !is.data.frame(x)) return(default)
  for (name in names) if (!is.null(x[[name]])) return(x[[name]])
  default
}
.eid_sg_adduct_name <- function(adduct, lipid = NULL) {
  if (is.character(adduct) && length(adduct) > 0L) return(adduct[1L])
  value <- .eid_sg_field(adduct, c("AdductIonName", "name", "Name"), adduct)
  if (is.list(value)) value <- .eid_sg_field(value, c("AdductIonName", "name", "Name"), NULL)
  if (!is.null(value) && length(value) == 1L && is.character(value)) return(value)
  .eid_sg_field(lipid, c("AdductIonName", "Adduct", "adduct"), "[M+H]+")
}
.eid_sg_ion_mz <- function(mass, adduct) {
  name <- .eid_sg_adduct_name(adduct)
  if (name == "[M+H]+") return(mass + .eid_sg_proton)
  if (name == "[M+Na]+") return(mass + .eid_sg_mass(hydrogen = -1) + 22.989218)
  if (name == "[M+NH4]+") return(mass + .eid_sg_mass(hydrogen = 3, nitrogen = 1) + .eid_sg_proton)
  if (name %in% c("[M+H-H2O]+", "[M-H2O+H]+")) return(mass - .eid_sg_h2o + .eid_sg_proton)
  if (name == "[M-H]-") return(mass - .eid_sg_mass(hydrogen = 1) + .eid_sg_electron)
  if (name == "[M+CH3COO]-") return(mass + .eid_sg_mass(carbon = 2, hydrogen = 3, oxygen = 2) - .eid_sg_electron)
  if (name == "[M+HCOO]-") return(mass + .eid_sg_mass(carbon = 1, hydrogen = 1, oxygen = 2) - .eid_sg_electron)
  mass
}
.eid_sg_peak <- function(mz, intensity, comment, type = "none") {
  data.frame(Mass = as.numeric(mz), Intensity = as.numeric(intensity),
             Comment = as.character(comment), SpectrumComment = as.character(type),
             stringsAsFactors = FALSE)
}
.eid_sg_add <- function(peaks, mz, intensity, comment, type = "none") {
  rbind(peaks, .eid_sg_peak(mz, intensity, comment, type))
}
.eid_sg_parse_chains <- function(lipid) {
  name <- as.character(.eid_sg_field(lipid, c("Name", "name"), ""))[1L]
  explicit <- .eid_sg_field(lipid, c("Chains", "chains"), NULL)
  if (is.data.frame(explicit)) {
    names(explicit) <- tolower(names(explicit))
    if (!"carbon" %in% names(explicit)) names(explicit)[1L] <- "carbon"
    if (!"double_bond" %in% names(explicit)) {
      db_name <- intersect(c("doublebond", "double_bonds", "db"), names(explicit))[1L]
      if (!is.na(db_name)) names(explicit)[names(explicit) == db_name] <- "double_bond"
    }
    if (!"positions" %in% names(explicit)) explicit$positions <- I(replicate(nrow(explicit), integer(), simplify = FALSE))
    explicit$carbon <- as.integer(explicit$carbon)
    explicit$double_bond <- as.integer(explicit$double_bond)
    return(explicit[, c("carbon", "double_bond", "positions"), drop = FALSE])
  }
  tokens <- regmatches(name, gregexpr("[0-9]+:[0-9]+(?:\\([0-9,]+\\))?", name, perl = TRUE))[[1L]]
  if (!length(tokens)) return(data.frame(carbon = integer(), double_bond = integer(), positions = I(list())))
  rows <- lapply(tokens, function(token) {
    base <- sub("\\(.*", "", token)
    parts <- strsplit(base, ":", fixed = TRUE)[[1L]]
    positions <- integer()
    if (grepl("\\(", token)) positions <- as.integer(strsplit(sub(".*\\((.*)\\).*", "\\1", token), ",", fixed = TRUE)[[1L]])
    data.frame(carbon = as.integer(parts[1L]), double_bond = as.integer(parts[2L]), positions = I(list(positions)))
  })
  do.call(rbind, rows)
}
.eid_sg_chain_mass <- function(chain) {
  .eid_sg_mass(carbon = chain$carbon, hydrogen = 2 * chain$carbon - 2 * chain$double_bond)
}
.eid_sg_class <- function(name) {
  value <- sub("[ _].*", "", as.character(name)[1L])
  sub("/.*", "", value)
}

.eid_sg_specific <- function(mass, chain, adduct, nl_mass, intensity) {
  positions <- unlist(chain$positions, use.names = FALSE)
  if (chain$double_bond <= 0L || !length(positions)) return(data.frame())
  if (chain$carbon <= 0L) return(data.frame())
  adjusted_nl <- nl_mass
  adjusted_nl <- adjusted_nl - .eid_sg_o + 2 * .eid_sg_h
  chain_loss <- mass - .eid_sg_chain_mass(chain) - adjusted_nl
  diffs <- rep(.eid_sg_ch2, chain$carbon)
  for (bond in positions) {
    if (bond >= 1L && bond < chain$carbon) {
      diffs[bond] <- diffs[bond] - .eid_sg_h
      diffs[bond + 1L] <- diffs[bond + 1L] - .eid_sg_h
    }
  }
  if (length(diffs) > 1L) diffs <- cumsum(diffs)
  peaks <- data.frame()
  add <- function(index, value, factor, bond) {
    if (index <= 0L || index > length(diffs)) return()
    peaks <<- .eid_sg_add(peaks, .eid_sg_ion_mz(chain_loss + value, adduct),
      intensity * factor, paste0(chain$carbon, ":", chain$double_bond, " db", bond, " EID specific(c", index, ")"), "doublebond")
  }
  if (chain$double_bond < 3L) {
    for (bond in positions) {
      if (bond == 1L) next
      factors <- c(.5, .5, .5, .7, 1, .7, .5)
      for (index in (bond - 2L):min(length(diffs), bond + 4L)) {
        offset <- index - (bond - 2L) + 1L
        value <- diffs[index]
        if (index != bond) value <- value + if (index >= bond) .eid_sg_h else -.eid_sg_h
        add(index, value, factors[offset], bond)
      }
    }
  } else if (length(positions) >= 3L && all(c(max(positions) - 3L, max(positions) - 6L) %in% positions)) {
    bond <- max(positions) - 6L
    if (length(positions) == 4L) {
      add(bond - 2L, diffs[bond - 3L] + .eid_sg_h, .5, bond)
      add(bond - 1L, diffs[bond - 2L], .5, bond)
      add(bond, diffs[bond - 1L], .75, bond)
      add(bond + 1L, diffs[bond] + .eid_sg_h, 1, bond)
      add(bond + 2L, diffs[bond + 1L] + .eid_sg_h, .5, bond)
    } else {
      add(bond, diffs[bond - 1L], .25, bond)
      add(bond + 1L, diffs[bond] + .eid_sg_h, 1, bond)
    }
  }
  if (all(c(5L, 8L, 11L) %in% positions)) {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(chain_loss + diffs[3L] - .eid_sg_h, adduct),
      intensity * .5, paste0(chain$carbon, ":", chain$double_bond, " C3 specific"), "doublebond")
  }
  peaks
}

.eid_sg_chain_peaks <- function(peaks, mass, chain, adduct, nl_mass, intensity = 40) {
  chain_mass <- .eid_sg_chain_mass(chain)
  label <- paste0(chain$carbon, ":", chain$double_bond)
  peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass - chain_mass + .eid_sg_h, adduct), intensity, paste0("-", label), "acylchain")
  peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass - chain_mass + 2 * .eid_sg_h, adduct), intensity, paste0("-", label, " +H"), "acylchain")
  peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass - chain_mass - .eid_sg_h2o + .eid_sg_h, adduct), intensity, paste0("-", label, " -H2O"), "acylchain")
  specific <- .eid_sg_specific(mass, chain, adduct, nl_mass, intensity)
  if (nrow(specific)) peaks <- rbind(peaks, specific)
  peaks
}
.eid_sg_adduct_mass <- function(adduct) {
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") return(.eid_sg_proton)
  .eid_sg_ion_mz(0, adduct)
}
.eid_sg_add_ion <- function(peaks, mass, adduct, intensity, comment, type = "metaboliteclass") {
  .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct), intensity, comment, type)
}
.eid_sg_get_pc_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC, adduct, 200, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC_GLY_C, adduct, 50, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC_GLY_O, adduct, 50, "Gly-O")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PC - 2 * .eid_sg_h, adduct, 25, "Precursor-Header")
  .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) / 2, 100, "[Precursor]2+", "metaboliteclass")
}
.eid_sg_get_pe_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PE, adduct, 500, "Precursor -C2H8NO4P")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE, adduct, 100, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_C, adduct, 100, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_O, adduct, 100, "Gly-O")
  .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) / 2, 100, "[Precursor]2+", "metaboliteclass")
}
.eid_sg_get_pg_spectrum <- function(mass, adduct) {
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_h2o, adduct, 100, "Precursor -H2O")
  peaks <- .eid_sg_add(peaks, mass - .eid_sg_header$PG + adduct_mass, 999, "Precursor -C3H9O6P", "metaboliteclass")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PG_NEUTRAL, adduct, 100, "Precursor -C3H6O2")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PG + adduct_mass, 100, "Header", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, (mass - .eid_sg_header$PG + .eid_sg_proton + adduct_mass) / 2, 150, "[Precursor -C3H9O6P]2+", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PG_GLY_C + adduct_mass, 100, "Gly-C", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PG_GLY_O + adduct_mass, 100, "Gly-O", "metaboliteclass")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 200, "[M+H]+")
  peaks
}
.eid_sg_get_pa_spectrum <- function(mass, adduct) {
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  header_loss <- .eid_sg_header$PA_H3PO4
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add(peaks, mass - header_loss + adduct_mass, 500, "- Header", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PG + adduct_mass, 150, "C3H9O6P", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PA_GLY_C + adduct_mass, 100, "Gly-C", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PA_GLY_O + adduct_mass, 100, "Gly-O", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, (mass - header_loss + adduct_mass) / 2, 200, "[Precursor-Header]2+", "metaboliteclass")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 200, "[M+H]+")
  peaks
}
.eid_sg_get_ps_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_h2o, adduct, 100, "Precursor -H2O")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PS_CHO2, adduct, 200, "Precursor -CHO2")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PS, adduct, 500, "Precursor -C3H8NO6P")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PS, adduct, 100, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PS_GLY_C, adduct, 100, "Gly-C")
  .eid_sg_add_ion(peaks, .eid_sg_header$PS_GLY_O, adduct, 100, "Gly-O")
}
.eid_sg_get_pi_spectrum <- function(mass, adduct) {
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PI + adduct_mass, 100, "Header", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PI_GLY_C + adduct_mass, 100, "Gly-C", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, .eid_sg_header$PI_GLY_O + adduct_mass, 100, "Gly-O", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, mass - .eid_sg_header$PI_C6H10O5 + adduct_mass, 100, "Precursor -C6H10O5", "metaboliteclass")
  if (.eid_sg_adduct_name(adduct) %in% c("[M+H]+", "[M+NH4]+")) peaks <- .eid_sg_add(peaks, mass - .eid_sg_header$PI + .eid_sg_proton, 700, "Precursor -Header", "metaboliteclass")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 200, "[M+H]+")
  if (.eid_sg_adduct_name(adduct) == "[M+Na]+") peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PI_C6H10O5, adduct, 500, "Precursor -C6H10O5")
  peaks
}
.eid_sg_get_lpc_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC, adduct, 250, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC_GLY_C, adduct, 100, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC_GLY_O, adduct, 100, "Gly-O")
  if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_h2o, adduct, 200, "C5H14NO+")
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_mass(5, 13, 1, 1), adduct, 200, "C5H14NO+")
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC - .eid_sg_h2o, adduct, 50, "Header - H2O")
  }
  peaks
}
.eid_sg_get_lpe_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_h2o, adduct, 200, "Precursor - H2O")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PE, adduct, 250, "Precursor -C2H8NO4P")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE, adduct, 200, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_C, adduct, 100, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_O, adduct, 200, "Gly-O")
  if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG, adduct, 100, "C3H9O6P")
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG - .eid_sg_h2o, adduct, 100, "C3H9O6P - H2O")
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_O - .eid_sg_h2o, adduct, 100, "Gly-O - H2O")
  }
  peaks
}
.eid_sg_get_tg_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) / 2, 150, "[Precursor]2+", "precursor")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_h2o, 150, "Precursor-H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 150, "[M+H]+", "metaboliteclass")
  } else if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_h2o, 100, "Precursor-H2O", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_etherpc_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC_GLY_C, adduct, 400, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC_GLY_O, adduct, 400, "Gly-O")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PC, adduct, 500, "Header")
  .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) / 2, 200, "[Precursor]2+", "metaboliteclass")
}
.eid_sg_get_dg_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_h2o, 150, "Precursor-H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 150, "[M+H]+", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton - .eid_sg_h2o, 150, "[M+H]+ -H2O", "metaboliteclass")
  } else if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_h2o, 100, "Precursor-H2O", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_sm_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$SM, adduct, 500,
    "C5H14NO4P (Header)")
  peaks <- .eid_sg_add_ion(peaks,
    .eid_sg_header$SM + .eid_sg_sm_c2h3no - .eid_sg_o, adduct, 150,
    "C7H18N2O4P (Header+C2H3N)")
  if (.eid_sg_adduct_name(adduct) == "[M+Na]+") {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_sm_c3h9n,
      150, "Precursor-C3H9N")
    peaks <- .eid_sg_add(peaks,
      .eid_sg_ion_mz(mass, adduct) - .eid_sg_header$SM, 150,
      "Precursor-Header")
  } else if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add_ion(peaks,
      .eid_sg_header$SM + .eid_sg_sm_c2h3no + .eid_sg_mass(1), adduct, 150,
      "C7H18N2O4P (Header+C3H3NO)")
  }
  peaks
}
.eid_sg_sm_chain_peaks <- function(peaks, mass, chain, adduct, position) {
  chain_mass <- .eid_sg_chain_mass(chain)
  if (position == 1L && .eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add(peaks,
      chain_mass + .eid_sg_h + .eid_sg_proton - 2 * .eid_sg_h2o, 100,
      "[sph+H]+ -Header -H2O", "acylchain")
  } else if (position == 2L) {
    charged_chain <- .eid_sg_ion_mz(chain_mass + .eid_sg_h, adduct)
    peaks <- .eid_sg_add(peaks, charged_chain + .eid_sg_sm_c2h2n, 200,
      "[FAA+C2H+adduct]+", "acylchain")
    peaks <- .eid_sg_add(peaks,
      charged_chain + .eid_sg_header$SM + .eid_sg_sm_c2h2n - .eid_sg_h,
      200, "[FAA+C2H+Header+adduct]+", "acylchain")
  }
  peaks
}
.eid_sg_get_sphingo_spectrum <- function(mass, chain, adduct) {
  .eid_sg_sm_chain_peaks(data.frame(), mass, chain, adduct, 1L)
}
.eid_sg_get_acyl_spectrum <- function(mass, chain, adduct) {
  .eid_sg_sm_chain_peaks(data.frame(), mass, chain, adduct, 2L)
}
.eid_sg_get_bmp_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG, adduct, 200, "C3H9O6P")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PG, adduct, 100, "Precursor -C3H9O6P")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG - .eid_sg_h2o, adduct, 200, "C3H9O6P - H2O")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 200, "[M+H]+", "metaboliteclass")
  }
  if (.eid_sg_adduct_name(adduct) %in% c("[M+H]+", "[M+NH4]+")) {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_h2o, 100, "Precursor -H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - 2 * .eid_sg_h2o, 100, "Precursor -2H2O", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_hbmp_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  precursor <- .eid_sg_ion_mz(mass, adduct)
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add(peaks, precursor - .eid_sg_h2o, 200, "Precursor -H2O", "metaboliteclass")
  peaks <- .eid_sg_add(peaks, precursor / 2, 200, "[Precursor]2+", "precursor")
  peaks <- .eid_sg_add(peaks, (precursor - .eid_sg_h2o) / 2, 200, "[Precursor -H2O]2+", "metaboliteclass")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG, adduct, 50, "C3H9O6P")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG - .eid_sg_h2o, adduct, 50, "C3H9O6P - H2O")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 200, "[M+H]+", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, (mass + .eid_sg_proton) / 2, 200, "[M+2H]2+", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_cl_spectrum <- function(mass, adduct) {
  precursor <- .eid_sg_ion_mz(mass, adduct)
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_h2o, adduct, 999, "Precursor -H2O")
  peaks <- .eid_sg_add(peaks, precursor / 2, 150, "[Precursor]2+", "precursor")
  peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass - .eid_sg_h2o, adduct) / 2, 150, "[Precursor -H2O]2+", "metaboliteclass")
  peaks
}
.eid_sg_get_mg_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton, 750, "[M+H]+", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + .eid_sg_proton - .eid_sg_h2o, 250, "[M+H]+ -H2O", "metaboliteclass")
  } else if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) - .eid_sg_h2o, 150, "Precursor-H2O", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_etherpe_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE, adduct, 100, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_C, adduct, 100, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PE_GLY_O, adduct, 100, "Gly-O")
  peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct) / 2, 100, "[Precursor]2+", "metaboliteclass")
  .eid_sg_add_ion(peaks, mass - .eid_sg_header$PE, adduct, 500, "Precursor -C2H8NO4P")
}
.eid_sg_get_lpg_spectrum <- function(mass, adduct) {
  name <- .eid_sg_adduct_name(adduct)
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  peaks <- data.frame()
  if (name == "[M+H]+") peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  if (name == "[M+NH4]+") {
    peaks <- .eid_sg_add_ion(peaks, mass, adduct, 500, "Precursor", "precursor")
    peaks <- .eid_sg_add(peaks, mass + adduct_mass, 999, "[M+H]+", "metaboliteclass")
  }
  if (name %in% c("[M+H]+", "[M+NH4]+")) {
    peaks <- .eid_sg_add(peaks, mass + adduct_mass - .eid_sg_c3h6o2, 100, "Precursor -C3H6O2", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + adduct_mass - .eid_sg_h2o, 500, "[M+H]+ -H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + adduct_mass - 2 * .eid_sg_h2o, 200, "[M+H]+ -2H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass + adduct_mass - .eid_sg_c3h6o2 - .eid_sg_h2o, 100, "[M+H]+ -C3H6O2 -H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, .eid_sg_header$PG + .eid_sg_proton, 300, "Header", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, .eid_sg_header$PG - .eid_sg_h2o + .eid_sg_proton, 200, "Header - H2O", "metaboliteclass")
    peaks <- .eid_sg_add(peaks, mass - .eid_sg_header$PG + .eid_sg_proton, 500, "Precursor -C3H9O6P", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_lpi_spectrum <- function(mass, adduct) {
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  peaks <- data.frame()
  add_mass <- function(value, intensity, comment, type = "metaboliteclass") {
    peaks <<- .eid_sg_add(peaks, value + adduct_mass, intensity, comment, type)
  }
  add_mass(.eid_sg_header$PI, 300, "Header")
  add_mass(.eid_sg_header$PI_GLY_C, 100, "Gly-C")
  add_mass(.eid_sg_header$PI_GLY_O, 200, "Gly-O")
  add_mass(.eid_sg_header$PG, 300, "C3H9O6P")
  add_mass(.eid_sg_header$PG - .eid_sg_h2o, 300, "C3H9O6P - H2O")
  add_mass(mass - .eid_sg_header$PI, 500, "[M+H]+ -Header")
  add_mass(mass - .eid_sg_mass(6, 10, 0, 5), 100, "[M+H]+ -C6H10O5")
  add_mass(mass - .eid_sg_mass(6, 10, 0, 5) - .eid_sg_h2o, 150, "[M+H]+ -C6H12O6")
  add_mass(mass - .eid_sg_h2o, 800, "[M+H]+ -H2O")
  if (.eid_sg_adduct_name(adduct) == "[M+H]+") peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  if (.eid_sg_adduct_name(adduct) == "[M+NH4]+") {
    peaks <- .eid_sg_add_ion(peaks, mass, adduct, 800, "Precursor", "precursor")
    peaks <- .eid_sg_add(peaks, mass + adduct_mass, 999, "[M+H]+", "metaboliteclass")
  }
  peaks
}
.eid_sg_get_lps_spectrum <- function(mass, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass, adduct, 999, "Precursor", "precursor")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PS, adduct, 100, "Header")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PS_GLY_C, adduct, 100, "Gly-C")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PS_GLY_O, adduct, 100, "Gly-O")
  if (.eid_sg_adduct_name(adduct) == "[M+H]+") {
    peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_header$PS, adduct, 500, "Precursor -C3H8NO6P")
    peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_cho2, adduct, 200, "Precursor -CHO2")
    peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_h2o, adduct, 100, "Precursor -H2O")
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG, adduct, 100, "C3H9O6P")
    peaks <- .eid_sg_add_ion(peaks, .eid_sg_header$PG - .eid_sg_h2o, adduct, 100, "C3H9O6P - H2O")
  }
  peaks
}
.eid_sg_lyso_chain_peaks <- function(peaks, mass, chain, adduct, class) {
  chain_mass <- .eid_sg_chain_mass(chain) - .eid_sg_h
  if (!is.finite(chain_mass) || chain_mass == 0) return(peaks)
  label <- paste0(chain$carbon, ":", chain$double_bond)
  add_ion <- function(value, intensity, comment) {
    peaks <<- .eid_sg_add_ion(peaks, value, adduct, intensity, comment, "acylchain")
  }
  if (class == "LPC") {
    add_ion(chain_mass + .eid_sg_proton, 100, paste0(label, " acyl"))
    add_ion(mass - chain_mass, 100, paste0("-", label))
    add_ion(mass - chain_mass - .eid_sg_h2o, 100, paste0("-", label, "-O"))
  } else if (class == "LPE") {
    add_ion(chain_mass + .eid_sg_proton, 100, paste0(label, " acyl"))
    add_ion(mass - chain_mass, 50, paste0("-", label))
  } else if (class == "LPG") {
    add_ion(mass - chain_mass, 200, paste0("-", label))
    add_ion(chain_mass + .eid_sg_proton, 100, paste0(label, " acyl"))
  } else if (class == "LPI") {
    add_ion(mass - chain_mass, 150, paste0("-", label))
    add_ion(mass - chain_mass - .eid_sg_h2o, 100, paste0("-", label, " -H2O"))
    add_ion(chain_mass + .eid_sg_proton, 100, paste0(label, " acyl+"))
  } else if (class == "LPS") {
    add_ion(chain_mass + .eid_sg_proton, 100, paste0(label, " acyl"))
    add_ion(mass - chain_mass, 100, paste0("-", label))
  }
  peaks
}
.eid_sg_lyso_position_peaks <- function(peaks, mass, chain, adduct, class) {
  if (class == "LPE") {
    return(.eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(chain) - .eid_sg_o - .eid_sg_ch2 - .eid_sg_h,
      adduct, 100, "-CH2(Sn1)", "snposition"))
  }
  if (class == "LPG") {
    peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(chain) - .eid_sg_o - .eid_sg_ch2,
      adduct, 100, "-CH2(Sn1)", "snposition")
    return(.eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(chain) - .eid_sg_h2o - .eid_sg_h - .eid_sg_o - .eid_sg_ch2,
      adduct, 100, "-H2O -CH2(Sn1)", "snposition"))
  }
  if (class == "LPI") {
    return(.eid_sg_add_ion(peaks, mass - .eid_sg_h - .eid_sg_chain_mass(chain) - .eid_sg_o - .eid_sg_ch2,
      adduct, 100, "-CH2(Sn1)", "snposition"))
  }
  .eid_sg_add_ion(peaks, mass - (.eid_sg_chain_mass(chain) - .eid_sg_h) - .eid_sg_h2o - .eid_sg_ch2,
    adduct, 100, "-CH2(Sn1)", "snposition")
}
.eid_sg_etherpc_chain_peaks <- function(peaks, mass, alkyl, acyl, adduct) {
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o, adduct, 100, paste0("-", alkyl$carbon, ":", alkyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_o - .eid_sg_h, adduct, 200, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl), adduct, 100, paste0("-", alkyl$carbon, ":", alkyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) + .eid_sg_h, adduct, 100, paste0("-", acyl$carbon, ":", acyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_h2o, adduct, 100, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
  .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o - .eid_sg_ch2, adduct, 150, "-CH2(Sn1)", "snposition")
}
.eid_sg_get_etherpc_p_spectrum <- function(mass, alkyl, acyl, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o, adduct, 100, paste0("-", alkyl$carbon, ":", alkyl$double_bond, "-O"), "acylchain")
  .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_o - .eid_sg_h, adduct, 200, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
}
.eid_sg_get_etherpc_o_spectrum <- function(mass, alkyl, acyl, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl), adduct, 100, paste0("-", alkyl$carbon, ":", alkyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) + .eid_sg_h, adduct, 100, paste0("-", acyl$carbon, ":", acyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o, adduct, 100, paste0("-", alkyl$carbon, ":", alkyl$double_bond, "-O"), "acylchain")
  .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_h2o, adduct, 100, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
}
.eid_sg_etherpe_chain_peaks <- function(peaks, mass, alkyl, acyl, adduct) {
  pe <- .eid_sg_header$PE
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) + .eid_sg_h, adduct, 30, paste0("-", alkyl$carbon, ":", alkyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl), adduct, 30, paste0("-", acyl$carbon, ":", acyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o, adduct, 150, paste0("-", alkyl$carbon, ":", alkyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_h2o, adduct, 30, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_chain_mass(alkyl) + pe - .eid_sg_h, adduct, 250, "Sn1Ether+C2H8NO3P", "acylchain")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_chain_mass(alkyl) + pe - .eid_sg_mass(3, 3, 0, 4, 1) - .eid_sg_h, adduct, 150, "Sn1Ether+C2H8NO3P-H3PO4", "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - pe - .eid_sg_chain_mass(alkyl) + .eid_sg_h, adduct, 300, "NL of C2H8NO4P+Sn1", "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o - .eid_sg_ch2, adduct, 50, "-CH2(Sn1)", "snposition")
  .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - pe - .eid_sg_o - .eid_sg_ch2, adduct, 200, "- Header -CH2(Sn1)", "snposition")
}
.eid_sg_get_etherpe_p_spectrum <- function(mass, alkyl, acyl, adduct) {
  pe <- .eid_sg_header$PE
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) + .eid_sg_h, adduct, 30, paste0("-", alkyl$carbon, ":", alkyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl), adduct, 30, paste0("-", acyl$carbon, ":", acyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o, adduct, 150, paste0("-", alkyl$carbon, ":", alkyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_h2o, adduct, 30, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_chain_mass(alkyl) + pe - .eid_sg_h, adduct, 250, "Sn1Ether+C2H8NO3P", "acylchain")
  peaks <- .eid_sg_add_ion(peaks, .eid_sg_chain_mass(alkyl) + pe - .eid_sg_mass(3, 3, 0, 4) - .eid_sg_h, adduct, 150, "Sn1Ether+C2H8NO3P-H3PO4", "acylchain")
  .eid_sg_add_ion(peaks, mass - pe - .eid_sg_chain_mass(alkyl) + .eid_sg_h, adduct, 300, "NL of C2H8NO4P+Sn1", "acylchain")
}
.eid_sg_get_etherpe_o_spectrum <- function(mass, alkyl, acyl, adduct) {
  peaks <- data.frame()
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) + .eid_sg_h, adduct, 50, paste0("-", alkyl$carbon, ":", alkyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) + .eid_sg_h, adduct, 50, paste0("-", acyl$carbon, ":", acyl$double_bond), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_o, adduct, 200, paste0("-", alkyl$carbon, ":", alkyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_h2o, adduct, 200, paste0("-", acyl$carbon, ":", acyl$double_bond, "-O"), "acylchain")
  peaks <- .eid_sg_add(peaks, .eid_sg_chain_mass(acyl), 50, paste0(acyl$carbon, ":", acyl$double_bond, " acyl+"), "acylchain")
  peaks <- .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(alkyl) - .eid_sg_header$PE + .eid_sg_h, adduct, 50, paste0("- Header -", alkyl$carbon, ":", alkyl$double_bond), "acylchain")
  .eid_sg_add_ion(peaks, mass - .eid_sg_chain_mass(acyl) - .eid_sg_header$PE + .eid_sg_h, adduct, 50, paste0("- Header -", acyl$carbon, ":", acyl$double_bond), "acylchain")
}
.eid_sg_get_lyso_acyl_level_spectrum <- function(mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain) - .eid_sg_h
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  label <- paste0(chain$carbon, ":", chain$double_bond)
  peaks <- data.frame()
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_c3h6o2 + adduct_mass, 450, paste0(label, "+C3H4O2+H"), "acylchain")
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_c3h6o2 - .eid_sg_h2o + adduct_mass, 100, paste0(label, " + C3H4O2 -H2O +H"), "acylchain")
  peaks <- .eid_sg_add(peaks, mass - chain_mass - .eid_sg_h + adduct_mass, 50, paste0("-", label), "acylchain")
  .eid_sg_add(peaks, mass - chain_mass - .eid_sg_h2o + adduct_mass, 50, paste0("-", label, "-H2O"), "acylchain")
}
.eid_sg_get_sn1_position_spectrum <- function(mass, chain, adduct, neutral_loss = 0) {
  .eid_sg_add(data.frame(), .eid_sg_ion_mz(mass - .eid_sg_chain_mass(chain) - .eid_sg_o - .eid_sg_ch2 - neutral_loss, adduct), 150, "-CH2(Sn1)", "snposition")
}
.eid_sg_get_acyl_level_spectrum <- function(mass, chain, adduct, class = NULL) {
  peaks <- data.frame()
  if (!is.null(class) && class == "TG") return(.eid_sg_tg_chain_peaks(peaks, mass, chain, adduct))
  if (!is.null(class) && class == "DG") return(.eid_sg_dg_chain_peaks(peaks, mass, chain, adduct))
  .eid_sg_chain_peaks(peaks, mass, chain, adduct, 0, 40)
}
.eid_sg_get_acyl_position_spectrum <- function(mass, chain, adduct, neutral_loss = 0) {
  .eid_sg_get_sn1_position_spectrum(mass, chain, adduct, neutral_loss)
}
.eid_sg_get_hbmp_acyl_position_spectrum <- function(mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain) - .eid_sg_h
  .eid_sg_add(data.frame(), mass + .eid_sg_adduct_mass(adduct) - chain_mass - .eid_sg_h2o - .eid_sg_ch2,
    100, "-CH2(Sn1)", "snposition")
}
.eid_sg_get_acyl_double_bond_spectrum <- function(mass, chain, adduct, neutral_loss = 0, intensity = 10) {
  .eid_sg_specific(mass, chain, adduct, neutral_loss, intensity)
}
.eid_sg_dg_chain_peaks <- function(peaks, mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain)
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  lipid_mass <- mass + adduct_mass
  if (.eid_sg_adduct_name(adduct) == "[M+Na]+") {
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h, 50, paste0("-", chain$carbon, ":", chain$double_bond), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_h - .eid_sg_o, 200, paste0("-", chain$carbon, ":", chain$double_bond, "-O"), "acylchain")
  } else {
    peaks <- .eid_sg_add(peaks, chain_mass - .eid_sg_h + .eid_sg_proton, 100, paste0(chain$carbon, ":", chain$double_bond, " acyl"), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h, 50, paste0("-", chain$carbon, ":", chain$double_bond), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h - .eid_sg_h2o, 400, paste0("-", chain$carbon, ":", chain$double_bond, "-O"), "acylchain")
  }
  peaks
}
.eid_sg_tg_chain_peaks <- function(peaks, mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain)
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  lipid_mass <- mass + adduct_mass
  label <- paste0(chain$carbon, ":", chain$double_bond)
  if (.eid_sg_adduct_name(adduct) == "[M+Na]+") {
    peaks <- .eid_sg_add(peaks, chain_mass + adduct_mass + .eid_sg_mass(2, 3, 0, 1) + .eid_sg_o, 100, paste0(label, "+C2H3O2"), "acylchain")
    peaks <- .eid_sg_add(peaks, chain_mass + adduct_mass + .eid_sg_mass(2, 3, 0, 1) + .eid_sg_ch2, 100, paste0(label, "+C3H5O"), "acylchain")
    peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_electron, 100, paste0(label, "+"), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h, 50, paste0("-", label), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h - .eid_sg_h2o, 200, paste0("-", label, "-O"), "acylchain")
  } else {
    peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_electron, 100, paste0(label, "+"), "acylchain")
    peaks <- .eid_sg_add(peaks, chain_mass + adduct_mass + .eid_sg_mass(3, 5, 0, 2), 100, paste0(label, "+C3H5O2"), "acylchain")
    peaks <- .eid_sg_add(peaks, chain_mass + adduct_mass + .eid_sg_mass(3, 5, 0, 2) - .eid_sg_h2o, 100, paste0(label, "+C3H3O"), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h, 50, paste0("-", label), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h - .eid_sg_h2o, 200, paste0("-", label, "-O"), "acylchain")
  }
  peaks
}
.eid_sg_tg_position_peaks <- function(peaks, mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain) + .eid_sg_adduct_mass(adduct)
  lipid_mass <- mass + .eid_sg_adduct_mass(adduct)
  if (.eid_sg_adduct_name(adduct) == "[M+Na]+") {
    peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_mass(2, 3, 0, 1), 100, "Sn2 diagnostics", "snposition")
    peaks <- .eid_sg_add(peaks, .eid_sg_chain_mass(chain) + .eid_sg_electron, 100, paste0(chain$carbon, ":", chain$double_bond, "+ Sn2"), "acylchain")
    return(.eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_o, 200, "-chain-O Sn2", "acylchain"))
  }
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_mass(2, 3, 0, 1), 100, "Sn2 diagnostics", "snposition")
  peaks <- .eid_sg_add(peaks, .eid_sg_chain_mass(chain) + .eid_sg_electron, 100, paste0(chain$carbon, ":", chain$double_bond, "+ Sn2"), "acylchain")
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_mass(3, 5, 0, 2), 100, paste0(chain$carbon, "+C3H5O2 Sn2"), "acylchain")
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_mass(3, 5, 0, 2) - .eid_sg_h2o, 100, paste0(chain$carbon, "+C3H3O Sn2"), "acylchain")
  peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + 2 * .eid_sg_h, 50, paste0("-", chain$carbon, ":", chain$double_bond, " Sn2"), "acylchain")
  .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_o, 200, paste0("-", chain$carbon, ":", chain$double_bond, "-O Sn2"), "acylchain")
}
.eid_sg_bmp_chain_peaks <- function(peaks, mass, chain, adduct, class) {
  chain_mass <- .eid_sg_chain_mass(chain)
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  lipid_mass <- mass + adduct_mass
  label <- paste0(chain$carbon, ":", chain$double_bond)
  if (class == "BMP") {
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h, 50, paste0("-", label), "acylchain")
    if (.eid_sg_adduct_name(adduct) == "[M+Na]+") {
      peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_h - .eid_sg_o, 50, paste0("-", label, "-OH"), "acylchain")
      peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_c3h6o2, 400, paste0("-C3H6O2 -", label), "acylchain")
      return(.eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_c3h6o2 - .eid_sg_h2o, 50, paste0("-C3H6O2 -", label, "-H2O"), "acylchain"))
    }
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h - .eid_sg_h2o, 50, paste0("-", label, "-H2O"), "acylchain")
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h - .eid_sg_header$PG, 400, paste0("-C3H9O6P -", label), "acylchain")
    return(.eid_sg_add(peaks, lipid_mass - chain_mass + .eid_sg_h - .eid_sg_header$PG - .eid_sg_h2o, 50, paste0("-C3H9O6P -", label, "-H2O"), "acylchain"))
  }
  peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass, 50, paste0("-", label), "acylchain")
  peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_h2o, 50, paste0("-", label, "-H2O"), "acylchain")
  .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_header$PG, 450, paste0("-C3H9O6P -", label), "acylchain")
}
.eid_sg_hbmp_lyso_chain_peaks <- function(peaks, mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain) - .eid_sg_h
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  label <- paste0(chain$carbon, ":", chain$double_bond)
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_c3h6o2 + adduct_mass, 450, paste0(label, "+C3H4O2+H"), "acylchain")
  peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_c3h6o2 - .eid_sg_h2o + adduct_mass, 100, paste0(label, " + C3H4O2 -H2O +H"), "acylchain")
  peaks <- .eid_sg_add(peaks, mass - chain_mass - .eid_sg_h + adduct_mass, 50, paste0("-", label), "acylchain")
  .eid_sg_add(peaks, mass - chain_mass - .eid_sg_h2o + adduct_mass, 50, paste0("-", label, "-H2O"), "acylchain")
}
.eid_sg_bmp_position_peaks <- function(peaks, mass, chain, adduct) {
  chain_mass <- .eid_sg_chain_mass(chain) - .eid_sg_h
  lipid_mass <- mass + .eid_sg_adduct_mass(adduct)
  peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass - .eid_sg_h2o - .eid_sg_ch2, 100, "-CH2(Sn1)", "snposition")
  if (.eid_sg_adduct_name(adduct) != "[M+Na]+") {
    peaks <- .eid_sg_add(peaks, lipid_mass - chain_mass - 2 * .eid_sg_h2o - .eid_sg_ch2, 50, "-H2O -CH2(Sn1)", "snposition")
  }
  peaks
}
.eid_sg_cl_chain_peaks <- function(peaks, mass, chains, adduct) {
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  if (nrow(chains) < 4L) return(peaks)
  pair1 <- .eid_sg_chain_mass(chains[1L, , drop = FALSE]) + .eid_sg_chain_mass(chains[2L, , drop = FALSE]) + .eid_sg_c3h3o2 + .eid_sg_h
  pair2 <- .eid_sg_chain_mass(chains[3L, , drop = FALSE]) + .eid_sg_chain_mass(chains[4L, , drop = FALSE]) + .eid_sg_c3h3o2 + .eid_sg_h
  peaks <- .eid_sg_add(peaks, mass - pair1 + .eid_sg_proton, 30, "[M-Sn1-Sn2-C3H3O2+H]+", "acylchain")
  peaks <- .eid_sg_add(peaks, mass - pair2 + .eid_sg_proton, 30, "[M-Sn3-Sn4-C3H3O2+H]+", "acylchain")
  peaks <- .eid_sg_add(peaks, mass - pair1 + .eid_sg_proton - .eid_sg_h2o, 30, "[M-Sn1-Sn2-C3H3O2-H2O+H]+", "acylchain")
  peaks <- .eid_sg_add(peaks, mass - pair2 + .eid_sg_proton - .eid_sg_h2o, 30, "[M-Sn3-Sn4-C3H3O2-H2O+H]+", "acylchain")
  peaks <- .eid_sg_add(peaks, pair1 + .eid_sg_proton, 500, "[Sn1+Sn2+C3H3O2+H]+", "acylchain")
  peaks <- .eid_sg_add(peaks, pair2 + .eid_sg_proton, 500, "[Sn3+Sn4+C3H3O2+H]+", "acylchain")
  for (i in seq_len(nrow(chains))) {
    chain_mass <- .eid_sg_chain_mass(chains[i, , drop = FALSE]) - .eid_sg_h
    label <- paste0(chains$carbon[i], ":", chains$double_bond[i])
    peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_proton, 50, paste0(label, " acyl+"), "acylchain")
    peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_c3h6o2 + .eid_sg_proton, 100, paste0(label, " +C3H6O2"), "acylchain")
    peaks <- .eid_sg_add(peaks, chain_mass + .eid_sg_c3h6o2 - .eid_sg_o + .eid_sg_proton, 50, paste0(label, " +C3H6O"), "acylchain")
    peaks <- .eid_sg_add(peaks, mass - chain_mass + .eid_sg_proton - adduct_mass, 20, paste0("-", label), "acylchain")
    peaks <- .eid_sg_add(peaks, mass - chain_mass + .eid_sg_proton - adduct_mass - .eid_sg_h2o, 20, paste0("-", label, "- H2O"), "acylchain")
  }
  peaks
}
.eid_sg_class_spectrum <- function(class, mass, adduct) {
  switch(class,
    PC = .eid_sg_get_pc_spectrum(mass, adduct),
    PE = .eid_sg_get_pe_spectrum(mass, adduct),
    PG = .eid_sg_get_pg_spectrum(mass, adduct),
    PA = .eid_sg_get_pa_spectrum(mass, adduct),
    PS = .eid_sg_get_ps_spectrum(mass, adduct),
    PI = .eid_sg_get_pi_spectrum(mass, adduct),
    LPC = .eid_sg_get_lpc_spectrum(mass, adduct),
    LPE = .eid_sg_get_lpe_spectrum(mass, adduct),
    TG = .eid_sg_get_tg_spectrum(mass, adduct),
    EtherPC = .eid_sg_get_etherpc_spectrum(mass, adduct),
    DG = .eid_sg_get_dg_spectrum(mass, adduct),
    SM = .eid_sg_get_sm_spectrum(mass, adduct),
    BMP = .eid_sg_get_bmp_spectrum(mass, adduct),
    HBMP = .eid_sg_get_hbmp_spectrum(mass, adduct),
    CL = .eid_sg_get_cl_spectrum(mass, adduct),
    MG = .eid_sg_get_mg_spectrum(mass, adduct),
    EtherPE = .eid_sg_get_etherpe_spectrum(mass, adduct),
    LPG = .eid_sg_get_lpg_spectrum(mass, adduct),
    LPI = .eid_sg_get_lpi_spectrum(mass, adduct),
    LPS = .eid_sg_get_lps_spectrum(mass, adduct),
    NULL)
}
.eid_sg_class_chain_peaks <- function(peaks, mass, chain, adduct, class) {
  chain_mass <- .eid_sg_chain_mass(chain)
  acyl_mass <- chain_mass - .eid_sg_h
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  add <- function(value, intensity, comment, type = "acylchain") {
    peaks <<- .eid_sg_add(peaks, value, intensity, comment, type)
  }
  if (class == "PC") {
    add(.eid_sg_ion_mz(mass - acyl_mass, adduct), 50, paste0("-", chain$carbon, ":", chain$double_bond))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_h2o, adduct), 50, paste0("-", chain$carbon, ":", chain$double_bond, "-O"))
  } else if (class == "PE") {
    add(acyl_mass + .eid_sg_proton, 50, paste0(chain$carbon, ":", chain$double_bond, " acyl+"))
    add(.eid_sg_ion_mz(mass - acyl_mass, adduct), 100, paste0("-", chain$carbon, ":", chain$double_bond))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_header$PE, adduct), 100, paste0("-Header -", chain$carbon, ":", chain$double_bond))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_o - .eid_sg_h, adduct), 100, paste0("-", chain$carbon, ":", chain$double_bond, "-O"))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_header$PE - .eid_sg_h2o, adduct), 100, paste0("-Header -", chain$carbon, ":", chain$double_bond, "-O"))
  } else if (class == "PG") {
    add(acyl_mass + .eid_sg_proton, 100, paste0(chain$carbon, ":", chain$double_bond, " acyl"))
    add(mass - acyl_mass - .eid_sg_h2o + adduct_mass, 100, paste0("-", chain$carbon, ":", chain$double_bond, "-H2O"))
    add(mass - acyl_mass - .eid_sg_header$PG + adduct_mass, 400, paste0("-Header -", chain$carbon, ":", chain$double_bond))
    add(mass - acyl_mass - .eid_sg_header$PG - .eid_sg_h2o + adduct_mass, 100, paste0("-Header -", chain$carbon, ":", chain$double_bond, "-H2O"))
  } else if (class == "PA") {
    add(acyl_mass + .eid_sg_proton, 100, paste0(chain$carbon, ":", chain$double_bond, " acyl"))
    add(mass - acyl_mass + adduct_mass, 50, paste0("-", chain$carbon, ":", chain$double_bond))
    add(mass - acyl_mass - .eid_sg_h2o + adduct_mass, 100, paste0("-", chain$carbon, ":", chain$double_bond, "-O"))
    add(mass - acyl_mass - .eid_sg_header$PA_H3PO4 + adduct_mass, 200, paste0("- Header -", chain$carbon, ":", chain$double_bond))
    add(mass - acyl_mass - .eid_sg_header$PA_H3PO4 - .eid_sg_h2o + adduct_mass, 50, paste0("- Header -", chain$carbon, ":", chain$double_bond, "-O"))
  } else if (class == "PS") {
    add(acyl_mass + .eid_sg_proton, 100, paste0(chain$carbon, ":", chain$double_bond, " acyl"))
    add(.eid_sg_ion_mz(mass - acyl_mass, adduct), 50, paste0("-", chain$carbon, ":", chain$double_bond))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_h - .eid_sg_o, adduct), 50, paste0("-", chain$carbon, ":", chain$double_bond, "-O"))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_header$PS, adduct), 150, paste0("-Header -", chain$carbon, ":", chain$double_bond))
    add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_header$PS - .eid_sg_h - .eid_sg_o, adduct), 50, paste0("-Header -", chain$carbon, ":", chain$double_bond, "-O"))
  } else if (class == "PI") {
    if (.eid_sg_adduct_name(adduct) %in% c("[M+H]+", "[M+NH4]+")) {
      add(.eid_sg_ion_mz(mass - acyl_mass, adduct), 100, paste0("-", chain$carbon, ":", chain$double_bond))
      add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_header$PI, adduct), 300, paste0("-C6H13O9P -", chain$carbon, ":", chain$double_bond))
    } else if (.eid_sg_adduct_name(adduct) == "[M+Na]+") {
      add(.eid_sg_ion_mz(mass - acyl_mass, adduct), 100, paste0("-", chain$carbon, ":", chain$double_bond))
      add(.eid_sg_ion_mz(mass - acyl_mass - .eid_sg_header$PI_C6H10O5, adduct), 100, paste0("-C6H10O5 -", chain$carbon, ":", chain$double_bond))
    }
  } else {
    return(.eid_sg_chain_peaks(peaks, mass, chain, adduct, 0, 40))
  }
  peaks
}
.eid_sg_class_position_peaks <- function(peaks, mass, chain, adduct, class) {
  chain_mass <- .eid_sg_chain_mass(chain) - .eid_sg_h
  adduct_mass <- .eid_sg_adduct_mass(adduct)
  add <- function(value, intensity, comment) {
    peaks <<- .eid_sg_add(peaks, value, intensity, comment, "snposition")
  }
  if (class == "PC") {
    add(.eid_sg_ion_mz(mass - .eid_sg_chain_mass(chain) - .eid_sg_o - .eid_sg_ch2, adduct), 50, "-CH2(Sn1)")
  } else if (class == "PE") {
    add(.eid_sg_ion_mz(mass - .eid_sg_chain_mass(chain) - .eid_sg_o - .eid_sg_ch2, adduct), 100, "-CH2(Sn1)")
    add(.eid_sg_ion_mz(mass - .eid_sg_chain_mass(chain) - .eid_sg_header$PE - .eid_sg_o - .eid_sg_ch2, adduct), 100, "-Header -CH2(Sn1)")
  } else if (class == "PG") {
    add(mass + adduct_mass - chain_mass - .eid_sg_h2o - .eid_sg_ch2, 100, "-CH2(Sn1)")
    add(mass + adduct_mass - chain_mass - 2 * .eid_sg_h2o - .eid_sg_ch2, 100, "-H2O -CH2(Sn1)")
  } else if (class == "PA") {
    add(mass - chain_mass - .eid_sg_h2o - .eid_sg_ch2 + adduct_mass, 40, "-CH2(Sn1)")
    add(mass - chain_mass - .eid_sg_header$PA_H3PO4 - .eid_sg_h2o - .eid_sg_ch2 + adduct_mass, 40, "- Header -CH2(Sn1)")
  } else if (class == "PS") {
    add(mass + adduct_mass - chain_mass - .eid_sg_h2o - .eid_sg_ch2, 100, "-CH2(Sn1)")
  } else if (class == "PI") {
    add(mass + adduct_mass - chain_mass - .eid_sg_h2o - .eid_sg_ch2, 100, "-CH2(Sn1)")
  }
  peaks
}
.eid_sg_all_chains <- function(peaks, mass, chains, adduct, nl_mass = 0, intensity = 40) {
  if (!nrow(chains)) return(peaks)
  for (i in seq_len(nrow(chains))) {
    next_nl <- nl_mass
    if (nrow(chains) > 1L && i < nrow(chains)) next_nl <- .eid_sg_chain_mass(chains[i + 1L, , drop = FALSE])
    peaks <- .eid_sg_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct, next_nl, intensity)
  }
  peaks
}
.eid_sg_header_peaks <- function(class, mass, adduct, chains) {
  name <- .eid_sg_adduct_name(adduct)
  peaks <- data.frame()
  add <- function(mz, intensity, comment, type = "metaboliteclass") peaks <<- .eid_sg_add(peaks, .eid_sg_ion_mz(mz, adduct), intensity, comment, type)
  if (class %in% c("PC", "SM", "EtherPC")) {
    if (name %in% c("[M+H]+", "[M+NH4]+")) { add(mass, 999, "Precursor", "precursor"); add(.eid_sg_header$PC_FRAGMENT, 500, "C5H14NO") }
  } else if (class %in% c("PE", "EtherPE")) {
    add(mass, 999, "Precursor", "precursor"); add(.eid_sg_header$PE, 500, "C2H8NO4P")
  } else if (class %in% c("PG", "LPG")) {
    add(mass, 999, "Precursor", "precursor"); add(.eid_sg_header$PG, 500, "C3H9O6P")
  } else if (class == "PS") {
    add(mass, 999, "Precursor", "precursor"); add(.eid_sg_header$PS, 500, "C3H8NO6P")
  } else if (class %in% c("PI", "LPI")) {
    add(mass, 999, "Precursor", "precursor"); add(.eid_sg_header$PI, 500, "C6H13O9P")
  } else add(mass, 999, "Precursor", "precursor")
  peaks
}
.eid_sg_class_peaks <- function(mass, class, adduct, chains) {
  specific <- .eid_sg_class_spectrum(class, mass, adduct)
  if (!is.null(specific)) return(specific)
  peaks <- .eid_sg_header_peaks(class, mass, adduct, chains)
  name <- .eid_sg_adduct_name(adduct)
  if (class == "CE") return(.eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor"))
  if (class %in% c("DG", "TG", "MG", "BMP", "HBMP", "CL")) {
    peaks <- data.frame(); peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    if (name == "[M+NH4]+") peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass - .eid_sg_nh3, adduct), 300, "[M+H]+", "precursor")
  }
  if (class %in% c("Cer_NS", "HexCer_NS", "SM")) {
    peaks <- data.frame(); peaks <- .eid_sg_add(peaks, .eid_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
    if (class == "SM") peaks <- .eid_sg_add(peaks, .eid_sg_header$PC_FRAGMENT, 50, "C5H14NO")
  }
  peaks
}
.eid_sg_generate_class <- function(lipid, class, adduct, intensity = 40, molecule = NULL) {
  mass <- as.numeric(.eid_sg_field(lipid, c("Mass", "MolecularMass", "NeutralMass"), NA_real_))
  if (is.na(mass)) stop("lipid must contain a neutral Mass for EID generation")
  chains <- .eid_sg_parse_chains(lipid)
  peaks <- .eid_sg_class_peaks(mass, class, adduct, chains)
  if (class == "SM") {
    if (nrow(chains) >= 1L) peaks <- rbind(peaks, .eid_sg_get_sphingo_spectrum(mass, chains[1L, , drop = FALSE], adduct))
    if (nrow(chains) >= 2L) peaks <- rbind(peaks, .eid_sg_get_acyl_spectrum(mass, chains[2L, , drop = FALSE], adduct))
    return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  }
  if (class %in% c("EtherPC", "EtherPE")) {
    if (nrow(chains) >= 2L) {
      alkyl <- chains[1L, , drop = FALSE]
      acyl <- chains[2L, , drop = FALSE]
      has_position_one <- 1 %in% unlist(alkyl$positions, use.names = FALSE)
      if (class == "EtherPC") {
        peaks <- rbind(peaks, .eid_sg_get_sn1_position_spectrum(mass, alkyl, adduct))
        if (has_position_one) peaks <- rbind(peaks, .eid_sg_get_etherpc_p_spectrum(mass, alkyl, acyl, adduct))
        else peaks <- rbind(peaks, .eid_sg_get_etherpc_o_spectrum(mass, alkyl, acyl, adduct))
      }
      if (class == "EtherPE") {
        peaks <- rbind(peaks, .eid_sg_get_sn1_position_spectrum(mass, alkyl, adduct))
        if (has_position_one) peaks <- rbind(peaks, .eid_sg_get_etherpe_p_spectrum(mass, alkyl, acyl, adduct))
        else peaks <- rbind(peaks, .eid_sg_get_etherpe_o_spectrum(mass, alkyl, acyl, adduct))
      }
    }
    return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  }
  if (class %in% c("DG", "TG")) {
    if (nrow(chains)) {
      position <- if (class == "TG" && nrow(chains) >= 2L) 2L else 1L
      for (i in seq_len(nrow(chains))) {
        if (class == "TG" && i == position) next
        if (class == "DG") peaks <- .eid_sg_dg_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct)
        if (class == "TG") peaks <- .eid_sg_tg_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct)
      }
      if (class == "TG" && nrow(chains) >= 2L) peaks <- .eid_sg_tg_position_peaks(peaks, mass, chains[2L, , drop = FALSE], adduct)
      if (class == "DG") peaks <- .eid_sg_add(peaks, mass + .eid_sg_adduct_mass(adduct) - .eid_sg_chain_mass(chains[1L, , drop = FALSE]) + .eid_sg_h - .eid_sg_h2o - .eid_sg_ch2, 100, "-CH2(Sn1)", "snposition")
    }
    return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  }
  if (class %in% c("BMP", "HBMP")) {
    if (nrow(chains)) {
      for (i in seq_len(nrow(chains))) {
        if (class == "HBMP" && i == 1L) peaks <- .eid_sg_hbmp_lyso_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct)
        else peaks <- .eid_sg_bmp_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct, class = if (class == "HBMP") "BMP" else class)
      }
      if (class == "HBMP") {
        if (nrow(chains) >= 1L) peaks <- rbind(peaks, .eid_sg_get_hbmp_acyl_position_spectrum(mass, chains[1L, , drop = FALSE], adduct))
        if (nrow(chains) >= 2L) peaks <- rbind(peaks, .eid_sg_get_hbmp_acyl_position_spectrum(mass, chains[2L, , drop = FALSE], adduct))
      } else if (nrow(chains) >= 1L) {
        peaks <- .eid_sg_bmp_position_peaks(peaks, mass, chains[1L, , drop = FALSE], adduct)
      }
    }
    return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  }
  if (class == "CL") {
    peaks <- .eid_sg_cl_chain_peaks(peaks, mass, chains, adduct)
    return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  }
  if (class %in% c("Cer_NS", "HexCer_NS")) return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  if (class %in% c("LPC", "LPE", "LPG", "LPI", "LPS")) {
    if (nrow(chains)) {
      for (i in seq_len(nrow(chains))) {
        peaks <- .eid_sg_lyso_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct, class)
        if (i == 1L) peaks <- .eid_sg_lyso_position_peaks(peaks, mass, chains[i, , drop = FALSE], adduct, class)
      }
    }
    specific_intensity <- switch(class, LPC = 250, LPE = 300, LPS = 150, LPG = 100, LPI = 100, 100)
    if (nrow(chains)) {
      for (i in seq_len(nrow(chains))) {
        specific <- .eid_sg_specific(mass, chains[i, , drop = FALSE], adduct, 0, specific_intensity)
        if (nrow(specific)) peaks <- rbind(peaks, specific)
      }
    }
    return(.eid_sg_reference(lipid, mass, class, adduct, peaks, molecule))
  }
  if (class %in% c("PC", "PE", "PG", "PA", "PS", "PI")) {
    if (nrow(chains)) {
      for (i in seq_len(nrow(chains))) {
        peaks <- .eid_sg_class_chain_peaks(peaks, mass, chains[i, , drop = FALSE], adduct, class)
        if (i == 1L) peaks <- .eid_sg_class_position_peaks(peaks, mass, chains[i, , drop = FALSE], adduct, class)
        specific_intensity <- switch(class, PC = 200, PE = 100, PG = 200, PA = 50, PS = 30, PI = 30, 0)
        if (specific_intensity > 0) {
          specific <- .eid_sg_specific(mass, chains[i, , drop = FALSE], adduct, 0, specific_intensity)
          if (nrow(specific)) peaks <- rbind(peaks, specific)
        }
      }
    }
  } else {
    peaks <- .eid_sg_all_chains(peaks, mass, chains, adduct, 0, intensity)
  }
  .eid_sg_reference(lipid, mass, class, adduct, peaks, molecule)
}
.eid_sg_reference <- function(lipid, mass, class, adduct, peaks, molecule = NULL) {
  if (!nrow(peaks)) peaks <- .eid_sg_peak(.eid_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
  key <- paste(round(peaks$Mass, 8), peaks$Comment, peaks$SpectrumComment, sep = "\r")
  peaks <- peaks[!duplicated(key), , drop = FALSE]
  peaks <- peaks[order(peaks$Mass), , drop = FALSE]
  name <- as.character(.eid_sg_field(lipid, c("Name", "name"), ""))[1L]
  list(Name = name, PrecursorMz = .eid_sg_ion_mz(mass, adduct), Mz = .eid_sg_ion_mz(mass, adduct),
       AdductIonName = .eid_sg_adduct_name(adduct, lipid), AdductType = adduct,
       IonMode = if (grepl("-$", .eid_sg_adduct_name(adduct, lipid))) "Negative" else "Positive",
       Charge = if (grepl("2[+-]$", .eid_sg_adduct_name(adduct, lipid))) 2L else 1L,
       CompoundClass = class, Ontology = .eid_sg_field(molecule %||% lipid, "Ontology"),
       Formula = .eid_sg_field(molecule %||% lipid, "Formula"), InChIKey = .eid_sg_field(molecule %||% lipid, "InChIKey"),
       SMILES = .eid_sg_field(molecule %||% lipid, "SMILES"), Spectrum = peaks)
}
`%||%` <- function(x, y) if (is.null(x)) y else x

generate_pc_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "PC", adduct, molecule = molecule)
generate_pe_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "PE", adduct, molecule = molecule)
generate_ps_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "PS", adduct, molecule = molecule)
generate_pg_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "PG", adduct, molecule = molecule)
generate_pi_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "PI", adduct, molecule = molecule)
generate_pa_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "PA", adduct, molecule = molecule)
generate_dg_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "DG", adduct, molecule = molecule)
generate_tg_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "TG", adduct, molecule = molecule)
generate_mg_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "MG", adduct, molecule = molecule)
generate_bmp_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "BMP", adduct, molecule = molecule)
generate_hbmp_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "HBMP", adduct, molecule = molecule)
generate_cl_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "CL", adduct, molecule = molecule)
generate_ceramide_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "Cer_NS", adduct, molecule = molecule)
generate_hexcer_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "HexCer_NS", adduct, molecule = molecule)
generate_sm_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "SM", adduct, molecule = molecule)
generate_etherpc_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "EtherPC", adduct, molecule = molecule)
generate_etherpe_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "EtherPE", adduct, molecule = molecule)
generate_lpc_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "LPC", adduct, molecule = molecule)
generate_lpe_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "LPE", adduct, molecule = molecule)
generate_lpg_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "LPG", adduct, molecule = molecule)
generate_lpi_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "LPI", adduct, molecule = molecule)
generate_lps_eid_spectrum <- function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, "LPS", adduct, molecule = molecule)

get_eid_lipid_spectrum <- function(lipid, adduct = NULL, molecule = NULL) {
  if (is.null(lipid) || is.null(lipid$Name)) stop("lipid must contain Name")
  if (is.null(adduct)) adduct <- .eid_sg_field(lipid, c("AdductIonName", "Adduct", "adduct"), "[M+H]+")
  class <- .eid_sg_class(lipid$Name)
  generator <- switch(class,
    PC = generate_pc_eid_spectrum, PE = generate_pe_eid_spectrum,
    PS = generate_ps_eid_spectrum, PG = generate_pg_eid_spectrum,
    PI = generate_pi_eid_spectrum, PA = generate_pa_eid_spectrum,
    DG = generate_dg_eid_spectrum, TG = generate_tg_eid_spectrum,
    MG = generate_mg_eid_spectrum, BMP = generate_bmp_eid_spectrum,
    HBMP = generate_hbmp_eid_spectrum, CL = generate_cl_eid_spectrum,
    Cer_NS = generate_ceramide_eid_spectrum, HexCer_NS = generate_hexcer_eid_spectrum,
    SM = generate_sm_eid_spectrum, EtherPC = generate_etherpc_eid_spectrum,
    EtherPE = generate_etherpe_eid_spectrum, LPC = generate_lpc_eid_spectrum,
    LPE = generate_lpe_eid_spectrum, LPG = generate_lpg_eid_spectrum,
    LPI = generate_lpi_eid_spectrum, LPS = generate_lps_eid_spectrum,
    function(lipid, adduct, molecule = NULL) .eid_sg_generate_class(lipid, class, adduct, molecule = molecule))
  generator(lipid, adduct, molecule)
}
generate_eid_lipid_spectrum <- get_eid_lipid_spectrum
can_generate_eid_lipid_spectrum <- function(lipid, adduct = NULL) {
  if (is.null(lipid) || is.null(lipid$Name)) return(FALSE)
  name <- .eid_sg_adduct_name(adduct, lipid)
  name %in% c("[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+H-H2O]+", "[M-H2O+H]+")
}

EidLipidSpectrumGenerator <- generate_eid_lipid_spectrum
EidDefaultSpectrumGenerator <- generate_eid_lipid_spectrum
EidSpecificSpectrumGenerator <- .eid_sg_specific
CanGenerateEidSpectrum <- can_generate_eid_lipid_spectrum
PC_EidSpectrumGenerator <- generate_pc_eid_spectrum
PE_EidSpectrumGenerator <- generate_pe_eid_spectrum
PGEidSpectrumGenerator <- generate_pg_eid_spectrum
PIEidSpectrumGenerator <- generate_pi_eid_spectrum
DGEidSpectrumGenerator <- generate_dg_eid_spectrum
TGEidSpectrumGenerator <- generate_tg_eid_spectrum
MGEidSpectrumGenerator <- generate_mg_eid_spectrum
