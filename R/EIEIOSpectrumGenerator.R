 # EIEIO/default lipid reference-spectrum generation translated from
 # MSDIAL's FacadeLipidSpectrumGenerator.Default and SpectrumPeakGenerator.

.eieio_sg_mass <- function(carbon = 0, hydrogen = 0, nitrogen = 0,
													 oxygen = 0, phosphorus = 0, sulfur = 0) {
	carbon * 12 + hydrogen * 1.00782503223 + nitrogen * 14.00307400443 +
		oxygen * 15.99491461957 + phosphorus * 30.97376199842 +
		sulfur * 31.9720711744
}
.eieio_sg_h <- 1.00782503223
.eieio_sg_proton <- 1.00727641974
.eieio_sg_electron <- 0.00054858026
.eieio_sg_o <- .eieio_sg_mass(oxygen = 1)
.eieio_sg_ch2 <- .eieio_sg_mass(carbon = 1, hydrogen = 2)
.eieio_sg_h2o <- .eieio_sg_mass(hydrogen = 2, oxygen = 1)
.eieio_sg_cho2 <- .eieio_sg_mass(carbon = 1, hydrogen = 1, oxygen = 2)
.eieio_sg_c3h8no6p <- .eieio_sg_mass(carbon = 3, hydrogen = 8, nitrogen = 1, oxygen = 6, phosphorus = 1)
.eieio_sg_c3h5no2 <- .eieio_sg_mass(carbon = 3, hydrogen = 5, nitrogen = 1, oxygen = 2)
.eieio_sg_c3h9o6p <- .eieio_sg_mass(carbon = 3, hydrogen = 9, oxygen = 6, phosphorus = 1)
.eieio_sg_lps_gly_c <- .eieio_sg_mass(carbon = 6, hydrogen = 12, nitrogen = 1, oxygen = 6, phosphorus = 1)
.eieio_sg_lps_gly_o <- .eieio_sg_mass(carbon = 5, hydrogen = 10, nitrogen = 1, oxygen = 7, phosphorus = 1)
.eieio_sg_c5h14no4p <- .eieio_sg_mass(carbon = 5, hydrogen = 14, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_c3h9n <- .eieio_sg_mass(carbon = 3, hydrogen = 9, nitrogen = 1)
.eieio_sg_c5h13no <- .eieio_sg_mass(carbon = 5, hydrogen = 13, nitrogen = 1, oxygen = 1)
.eieio_sg_c5h11n <- .eieio_sg_mass(carbon = 5, hydrogen = 11, nitrogen = 1)
.eieio_sg_lpc_gly_c <- .eieio_sg_mass(carbon = 8, hydrogen = 18, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_lpc_gly_o <- .eieio_sg_mass(carbon = 7, hydrogen = 16, nitrogen = 1, oxygen = 5, phosphorus = 1)
.eieio_sg_c2h8no4p <- .eieio_sg_mass(carbon = 2, hydrogen = 8, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_c2h5n <- .eieio_sg_mass(carbon = 2, hydrogen = 5, nitrogen = 1)
.eieio_sg_lpe_gly_c <- .eieio_sg_mass(carbon = 5, hydrogen = 12, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_lpe_gly_o <- .eieio_sg_mass(carbon = 4, hydrogen = 10, nitrogen = 1, oxygen = 5, phosphorus = 1)
.eieio_sg_lpg_gly_c <- .eieio_sg_mass(carbon = 6, hydrogen = 13, oxygen = 6, phosphorus = 1)
.eieio_sg_lpg_gly_o <- .eieio_sg_mass(carbon = 5, hydrogen = 11, oxygen = 7, phosphorus = 1)
.eieio_sg_c3h6o2 <- .eieio_sg_mass(carbon = 3, hydrogen = 6, oxygen = 2)
.eieio_sg_c6h13o9p <- .eieio_sg_mass(carbon = 6, hydrogen = 13, oxygen = 9, phosphorus = 1)
.eieio_sg_c6h10o5 <- .eieio_sg_mass(carbon = 6, hydrogen = 10, oxygen = 5)
.eieio_sg_lpi_gly_c <- .eieio_sg_mass(carbon = 9, hydrogen = 17, oxygen = 9, phosphorus = 1)
.eieio_sg_lpi_gly_o <- .eieio_sg_mass(carbon = 8, hydrogen = 15, oxygen = 10, phosphorus = 1)
.eieio_sg_c3h8no6p <- .eieio_sg_mass(carbon = 3, hydrogen = 8, nitrogen = 1, oxygen = 6, phosphorus = 1)
.eieio_sg_c3h5no2 <- .eieio_sg_mass(carbon = 3, hydrogen = 5, nitrogen = 1, oxygen = 2)
.eieio_sg_pg_gly_c <- .eieio_sg_mass(carbon = 6, hydrogen = 13, oxygen = 6, phosphorus = 1)
.eieio_sg_pg_gly_o <- .eieio_sg_mass(carbon = 5, hydrogen = 11, oxygen = 7, phosphorus = 1)
.eieio_sg_c3h3o2 <- .eieio_sg_mass(carbon = 3, hydrogen = 3, oxygen = 2)
.eieio_sg_h3po4 <- .eieio_sg_mass(hydrogen = 3, oxygen = 4, phosphorus = 1)
.eieio_sg_ps_gly_c <- .eieio_sg_mass(carbon = 6, hydrogen = 12, nitrogen = 1, oxygen = 6, phosphorus = 1)
.eieio_sg_ps_gly_o <- .eieio_sg_mass(carbon = 5, hydrogen = 10, nitrogen = 1, oxygen = 7, phosphorus = 1)
.eieio_sg_c7h13no2 <- .eieio_sg_mass(carbon = 7, hydrogen = 13, nitrogen = 1, oxygen = 2)
.eieio_sg_c8h16no3 <- .eieio_sg_mass(carbon = 8, hydrogen = 16, nitrogen = 1, oxygen = 3)
.eieio_sg_c6h11no2 <- .eieio_sg_mass(carbon = 6, hydrogen = 11, nitrogen = 1, oxygen = 2)
.eieio_sg_c2h5no <- .eieio_sg_mass(carbon = 2, hydrogen = 5, nitrogen = 1, oxygen = 1)
.eieio_sg_c6h10o5 <- .eieio_sg_mass(carbon = 6, hydrogen = 10, oxygen = 5)
.eieio_sg_c5h10o4 <- .eieio_sg_mass(carbon = 5, hydrogen = 10, oxygen = 4)
.eieio_sg_ch3o <- .eieio_sg_mass(carbon = 1, hydrogen = 3, oxygen = 1)
.eieio_sg_ch6o2 <- .eieio_sg_mass(carbon = 1, hydrogen = 6, oxygen = 2)
.eieio_sg_c2h8no4p <- .eieio_sg_mass(carbon = 2, hydrogen = 8, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_c5h14no4p <- .eieio_sg_mass(carbon = 5, hydrogen = 14, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_c3h9n <- .eieio_sg_mass(carbon = 3, hydrogen = 9, nitrogen = 1)
.eieio_sg_c5h11n <- .eieio_sg_mass(carbon = 5, hydrogen = 11, nitrogen = 1)
.eieio_sg_c5h14no4p <- .eieio_sg_mass(carbon = 5, hydrogen = 14, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_c7h15no3 <- .eieio_sg_mass(carbon = 7, hydrogen = 15, nitrogen = 1, oxygen = 3)
.eieio_sg_c11h15no7 <- .eieio_sg_mass(carbon = 11, hydrogen = 15, nitrogen = 1, oxygen = 7)
.eieio_sg_c8h16no6 <- .eieio_sg_mass(carbon = 8, hydrogen = 16, nitrogen = 1, oxygen = 6)
.eieio_sg_c2h8no4p <- .eieio_sg_mass(carbon = 2, hydrogen = 8, nitrogen = 1, oxygen = 4, phosphorus = 1)
.eieio_sg_so3 <- .eieio_sg_mass(sulfur = 1, oxygen = 3)
.eieio_sg_c3h5o2 <- .eieio_sg_mass(carbon = 3, hydrogen = 5, oxygen = 2)
.eieio_sg_c2h3o <- .eieio_sg_mass(carbon = 2, hydrogen = 3, oxygen = 1)
.eieio_sg_c2h2n <- .eieio_sg_mass(carbon = 2, hydrogen = 2, nitrogen = 1)
.eieio_sg_c2h3n <- .eieio_sg_mass(carbon = 2, hydrogen = 3, nitrogen = 1)
.eieio_sg_c2h3no <- .eieio_sg_mass(carbon = 2, hydrogen = 3, nitrogen = 1, oxygen = 1)
.eieio_sg_c2nh7 <- .eieio_sg_mass(carbon = 2, nitrogen = 1, hydrogen = 7)
.eieio_sg_c3h4no <- .eieio_sg_mass(carbon = 3, hydrogen = 4, nitrogen = 1, oxygen = 1)

.eieio_sg_field <- function(x, names, default = NULL) {
	if (is.null(x)) return(default)
	for (name in names) {
		value <- x[[name]]
		if (is.data.frame(value)) return(value)
		if (!is.null(value) && length(value) > 0L && !is.na(value[1L])) return(value[1L])
	}
	default
}

.eieio_sg_adduct_name <- function(adduct, lipid = NULL) {
	if (is.character(adduct) && length(adduct)) return(adduct[1L])
	value <- .eieio_sg_field(adduct, c("AdductIonName", "name", "Name"), NULL)
	if (!is.null(value)) return(as.character(value[1L]))
	as.character(.eieio_sg_field(lipid, c("AdductIonName", "Adduct", "adduct"), "[M+H]+"))[1L]
}

.eieio_sg_ion_mz <- function(mass, adduct) {
	switch(.eieio_sg_adduct_name(adduct),
		"[M+H]+" = mass + .eieio_sg_proton,
		"[M+Na]+" = mass + 22.989218,
		"[M+NH4]+" = mass + 18.033823,
		"[M+H-H2O]+" = mass - .eieio_sg_h2o + .eieio_sg_proton,
		"[M-H2O+H]+" = mass - .eieio_sg_h2o + .eieio_sg_proton,
		"[M-H]-" = mass - .eieio_sg_proton + .eieio_sg_electron,
		"[M+HCOO]-" = mass + 44.998201 - .eieio_sg_proton,
		"[M+CH3COO]-" = mass + 59.013851 - .eieio_sg_proton,
		"[M+HCO3]-" = mass + 60.021129 - .eieio_sg_proton,
		mass
	)
}

.eieio_sg_adduct_mass <- function(adduct) {
	switch(.eieio_sg_adduct_name(adduct),
		"[M+H]+" = .eieio_sg_proton,
		"[M+Na]+" = 22.989218,
		"[M+NH4]+" = 18.033823,
		"[M+H-H2O]+" = .eieio_sg_proton - .eieio_sg_h2o,
		0)
}

.eieio_sg_isotope_supported_adduct <- function(class, adduct) {
	if (!class %in% c("PC_d5", "PE_d5", "PG_d5", "PI_d5", "PS_d5", "LPC_d5", "LPE_d5", "LPG_d5", "LPI_d5", "LPS_d5", "DG_d5", "TG_d5", "CE_d7", "SM_d9", "Cer_NS_d7", "Cer_NP", "Cer_AP", "Cer_ABP")) return(TRUE)
	name <- .eieio_sg_adduct_name(adduct)
	allowed <- switch(class,
		PC_d5 = c("[M+H]+", "[M+Na]+"),
		PE_d5 = c("[M+H]+", "[M+Na]+"),
		PG_d5 = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		PI_d5 = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		PS_d5 = c("[M+H]+", "[M+Na]+"),
		LPC_d5 = c("[M+H]+", "[M+Na]+"),
		LPE_d5 = c("[M+H]+", "[M+Na]+"),
		LPG_d5 = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		LPI_d5 = c("[M+H]+", "[M+NH4]+"),
		LPS_d5 = c("[M+H]+", "[M+Na]+"),
		DG_d5 = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		TG_d5 = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		CE_d7 = c("[M+NH4]+", "[M+Na]+"),
		SM_d9 = c("[M+H]+", "[M+Na]+"),
		Cer_NS_d7 = c("[M+H]+", "[M+H-H2O]+", "[M+Na]+"),
		if (class %in% c("Cer_NP", "Cer_AP", "Cer_ABP")) c("[M+H]+", "[M+H-H2O]+") else character()
	)
	name %in% allowed
}

.eieio_sg_supported_adduct <- function(class, adduct) {
	name <- .eieio_sg_adduct_name(adduct)
	allowed <- switch(class,
		EtherPC = c("[M+H]+", "[M+Na]+"), EtherPE = c("[M+H]+", "[M+Na]+"),
		LPC = c("[M+H]+", "[M+Na]+"), LPE = c("[M+H]+", "[M+Na]+"),
		LPG = c("[M+H]+", "[M+NH4]+", "[M+Na]+"), LPI = c("[M+H]+", "[M+NH4]+"),
		LPS = c("[M+H]+", "[M+Na]+"), PA = c("[M+H]+", "[M+NH4]+"),
		PC = c("[M+H]+", "[M+Na]+"), PE = c("[M+H]+", "[M+Na]+"),
		PG = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		PI = c("[M+H]+", "[M+NH4]+", "[M+Na]+"), PS = c("[M+H]+", "[M+Na]+"),
		CL = c("[M+H]+", "[M+NH4]+"), BMP = c("[M+H]+", "[M+Na]+", "[M+NH4]+"),
		HBMP = c("[M+H]+", "[M+NH4]+"), DGTA = "[M+H]+", DGTS = "[M+H]+",
		LDGTA = "[M+H]+", LDGTS = "[M+H]+", CAR = "[M+H]+",
		DG = c("[M+H]+", "[M+NH4]+", "[M+Na]+"), TG = c("[M+H]+", "[M+NH4]+", "[M+Na]+"),
				Cer_AS = c("[M+H]+", "[M+H-H2O]+"), Cer_ADS = c("[M+H]+", "[M+H-H2O]+"),
				Cer_BS = c("[M+H]+", "[M+H-H2O]+"), Cer_BDS = c("[M+H]+", "[M+H-H2O]+"),
				Cer_NS = c("[M+H]+", "[M+H-H2O]+", "[M+Na]+"), Cer_NDS = c("[M+H]+", "[M+H-H2O]+", "[M+Na]+"),
				HexCer_NS = c("[M+H]+", "[M+Na]+", "[M+H-H2O]+"), HexCer_NDS = c("[M+H]+", "[M+Na]+", "[M+H-H2O]+"),
				GM3 = c("[M+H]+", "[M+NH4]+", "[M+Na]+"), SHexCer = c("[M+H]+", "[M+Na]+"),
				SM = c("[M+H]+", "[M+Na]+"), DMEDFA = "[M+H]+", DMEDFAHFA = "[M+H]+", DMEDOxFA = "[M+H]+",
		character())
	if (!length(allowed)) return(TRUE)
	name %in% allowed
}

.eieio_sg_empty <- function() {
	data.frame(Mass = numeric(), Intensity = numeric(), Comment = character(),
						 SpectrumComment = character(), IsAbsolutelyRequiredFragmentForAnnotation = logical(),
						 stringsAsFactors = FALSE)
}
.eieio_sg_peak <- function(mass, intensity, comment, type = "none", required = FALSE) {
	data.frame(Mass = as.numeric(mass), Intensity = as.numeric(intensity),
						 Comment = as.character(comment), SpectrumComment = as.character(type),
						 IsAbsolutelyRequiredFragmentForAnnotation = as.logical(required),
						 stringsAsFactors = FALSE)
}
.eieio_sg_rbind <- function(peaks, peak) {
	if (is.null(peaks) || !nrow(peaks)) return(peak)
	rbind(peaks, peak)
}
.eieio_sg_merge <- function(peaks) {
	if (is.null(peaks) || !nrow(peaks)) return(.eieio_sg_empty())
	key <- formatC(peaks$Mass, format = "f", digits = 8)
	groups <- split(seq_len(nrow(peaks)), key)
	merged <- lapply(groups, function(index) {
		first <- peaks[index[1L], , drop = FALSE]
		first$Intensity <- sum(peaks$Intensity[index])
		first$Comment <- paste(unique(peaks$Comment[index]), collapse = ", ")
		first$SpectrumComment <- paste(unique(peaks$SpectrumComment[index]), collapse = ",")
		first$IsAbsolutelyRequiredFragmentForAnnotation <- any(peaks$IsAbsolutelyRequiredFragmentForAnnotation[index])
		first
	})
	result <- do.call(rbind, merged)
	result[order(result$Mass), , drop = FALSE]
}

.eieio_sg_parse_chains <- function(lipid) {
	supplied <- .eieio_sg_field(lipid, c("Chains", "chains"), NULL)
	if (is.data.frame(supplied)) {
		names(supplied) <- tolower(names(supplied))
		carbon <- intersect(c("carbon", "carbon_count"), names(supplied))[1L]
		db <- intersect(c("double_bond", "doublebond", "db", "double_bond_count"), names(supplied))[1L]
		if (!is.na(carbon) && !is.na(db)) {
			positions <- intersect(c("positions", "bond_positions", "double_bond_positions"), names(supplied))[1L]
			if (is.na(positions)) supplied$positions <- I(replicate(nrow(supplied), integer(), simplify = FALSE)) else names(supplied)[names(supplied) == positions] <- "positions"
			supplied$carbon <- as.integer(supplied[[carbon]])
			supplied$double_bond <- as.integer(supplied[[db]])
			supplied$type <- if ("type" %in% names(supplied)) as.character(supplied$type) else "acyl"
			return(supplied[, intersect(c("carbon", "double_bond", "positions", "type", "mass", "oxidized"), names(supplied)), drop = FALSE])
		}
	}
	name <- as.character(.eieio_sg_field(lipid, c("Name", "name", "LipidName"), ""))[1L]
	tokens <- regmatches(name, gregexpr("[0-9]+:[0-9]+(?:\\([^)]*\\))?", name, perl = TRUE))[[1L]]
	if (!length(tokens)) return(data.frame(carbon = integer(), double_bond = integer(), positions = I(list()), type = character(), stringsAsFactors = FALSE))
	rows <- lapply(tokens, function(token) {
		base <- sub("\\(.*", "", token)
		parts <- strsplit(base, ":", fixed = TRUE)[[1L]]
		positions <- integer()
		if (grepl("\\(", token)) positions <- as.integer(strsplit(sub(".*\\((.*)\\).*", "\\1", token), ",", fixed = TRUE)[[1L]])
		ccarbon <- as.integer(parts[1L]); cdb <- as.integer(parts[2L])
		data.frame(carbon = ccarbon, double_bond = cdb, positions = I(list(positions)),
							 type = if (grepl("(^|[ _/])(?:O-|P-)", name, perl = TRUE) && length(rows) == 0L) "alkyl" else "acyl", stringsAsFactors = FALSE)
	})
	do.call(rbind, rows)
}

.eieio_sg_chain_mass <- function(chain) {
	carbon <- chain$carbon[1L]; db <- chain$double_bond[1L]
	type <- if ("type" %in% names(chain)) tolower(chain$type[1L]) else "acyl"
	if (type == "sphingo") .eieio_sg_mass(carbon, 2 * carbon - 2 * db, nitrogen = 1, oxygen = 2)
	else if (type == "alkyl") .eieio_sg_mass(carbon, 2 * carbon + 2 - 2 * db, oxygen = 1)
	else .eieio_sg_mass(carbon, 2 * carbon - 2 * db, oxygen = 2)
}
.eieio_sg_chain_label <- function(chain) paste0(chain$carbon[1L], ":", chain$double_bond[1L])
.eieio_sg_positions <- function(chain) unlist(chain$positions[[1L]], use.names = FALSE)

.eieio_sg_double_bond <- function(lipid_mass, chain, adduct, nl_mass = 0, abundance = 10, sphingo = FALSE) {
	positions <- .eieio_sg_positions(chain)
	carbon <- chain$carbon[1L]
	if (!length(positions) || carbon <= 0L) return(.eieio_sg_empty())
	chain_loss <- lipid_mass - .eieio_sg_chain_mass(chain) - nl_mass
	if (sphingo) chain_loss <- chain_loss + .eieio_sg_mass(2, 5, 1, 1)
	diffs <- rep(.eieio_sg_ch2, carbon)
	for (bond in positions) if (bond >= 1L && bond < carbon) {
		diffs[bond] <- diffs[bond] - .eieio_sg_h
		diffs[bond + 1L] <- diffs[bond + 1L] - .eieio_sg_h
	}
	diffs <- cumsum(diffs)
	start <- if (sphingo) 3L else 1L
	result <- .eieio_sg_empty()
	for (i in start:(carbon - 1L)) {
		factor <- 1; factor_hloss <- .5; factor_hgain <- if (sphingo) .2 else .05
		comment <- "doublebond"
		if ((i - 1L) %in% positions) { factor <- 4; factor_hloss <- 2; factor_hgain <- .05; comment <- "doublebond,doublebond_high" }
		else if (!sphingo && (i + 1L) %in% positions) { factor <- .5; factor_hloss <- .25; comment <- "doublebond,doublebond_low" }
		else if (sphingo && (i + 1L) %in% positions) { factor <- .25; factor_hloss <- .5; factor_hgain <- .05; comment <- "doublebond,doublebond_low" }
		else if ((i + 3L) %in% positions) { factor <- if (sphingo) 1 else 2; factor_hloss <- 4; comment <- "doublebond,doublebond_high" }
		if (!sphingo && all(c(i, i + 3L, i + 6L) %in% positions)) { factor_hgain <- 4; comment <- "doublebond,doublebond_high" }
		if (!sphingo && i == 3L) {
			if (1L %in% positions) { factor <- 2.5; factor_hloss <- .5; factor_hgain <- 0 } else { factor <- .75; factor_hloss <- 2; factor_hgain <- .5 }
		}
		if (!sphingo && i == 2L) { factor <- 1.5; factor_hloss <- .5; factor_hgain <- 2 }
		label <- paste0(.eieio_sg_chain_label(chain), " C", i, if (carbon >= 1L && length(positions) >= 3L) "_p3" else "")
		result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(chain_loss + diffs[i] - .eieio_sg_h, adduct), factor_hloss * abundance, paste0(label, "-H"), comment))
		result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(chain_loss + diffs[i], adduct), factor * abundance, label, comment))
		result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(chain_loss + diffs[i] + .eieio_sg_h, adduct), factor_hgain * abundance, paste0(label, "+H"), comment))
	}
	result
}

.eieio_sg_sphingolipid_double_bond <- function(lipid_mass, chain, adduct, class, abundance = 30) {
	balance <- if (class == "Cer_NS_d7") 7 * (.eieio_sg_d - .eieio_sg_h) else 0
	.eieio_sg_double_bond(lipid_mass + balance, chain, adduct,
		nl_mass = if (class == "Cer_NS_d7") -balance else 0,
		abundance = abundance, sphingo = tolower(chain$type[1L]) == "sphingo")
}

.eieio_sg_phyto_double_bond <- function(lipid_mass, chain, adduct, abundance = 30) {
	.eieio_sg_double_bond(lipid_mass, chain, adduct, nl_mass = 0,
		abundance = abundance, sphingo = tolower(chain$type[1L]) == "sphingo")
}

.eieio_sg_chain_peaks <- function(peaks, lipid_mass, chain, adduct, abundance = 50) {
	label <- .eieio_sg_chain_label(chain); chain_mass <- .eieio_sg_chain_mass(chain)
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass, adduct), abundance, paste0("-", label), "acylchain"))
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_h2o, adduct), abundance, paste0("-", label, "-O"), "acylchain"))
	peaks
}

.eieio_sg_header <- function(class, mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); result <- .eieio_sg_peak(.eieio_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
	formula <- switch(class,
		PC = c(5, 14, 1, 4, 1), PE = c(2, 8, 1, 4, 1), PG = c(3, 9, 0, 6, 1),
		PI = c(6, 13, 0, 9, 1), PS = c(3, 8, 1, 6, 1), PA = c(3, 7, 0, 6, 1),
		c(0, 0, 0, 0, 0))
	if (any(formula != 0)) {
		header <- .eieio_sg_mass(formula[1], formula[2], formula[3], formula[4], formula[5])
		result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(header, adduct), 300, "Header", "metaboliteclass", class %in% c("PC", "PE")))
		if (class == "PC") result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(mass - header - 2 * .eieio_sg_h, adduct), 25, "Precursor-Header", "metaboliteclass"))
		if (class == "PE") result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(mass - header, adduct), 500, "Precursor -C2H8NO4P", "metaboliteclass", TRUE))
	}
	if (name == "[M+Na]+" && class == "PC") result <- .eieio_sg_rbind(result, .eieio_sg_peak(.eieio_sg_ion_mz(mass - .eieio_sg_mass(3, 9, 1), adduct), 200, "Precursor-C3H9N", "metaboliteclass"))
	result
}

.eieio_sg_get_lps_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mass, intensity, comment, required = FALSE) {
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(mass, adduct), intensity,
			comment, "metaboliteclass", required))
	}
	add(.eieio_sg_c3h8no6p, 300, "Header", TRUE)
	add(.eieio_sg_lps_gly_c, 100, "Gly-C")
	add(.eieio_sg_lps_gly_o, 100, "Gly-O")
	if (name == "[M+H]+") {
		add(lipid_mass - .eieio_sg_c3h8no6p, 500, "Precursor -C3H8NO6P")
		add(lipid_mass - .eieio_sg_cho2, 200, "Precursor -CHO2")
		add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O")
		add(.eieio_sg_c3h9o6p, 300, "C3H9O6P")
		add(.eieio_sg_c3h9o6p - .eieio_sg_h2o, 400, "C3H9O6P - H2O")
	} else if (name == "[M+Na]+") {
		add(lipid_mass - .eieio_sg_c3h5no2, 200, "Precursor -C3H5NO2")
	}
	peaks
}

.eieio_sg_get_lpc_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE) {
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(value, adduct), intensity, comment, "metaboliteclass", required))
	}
	add(.eieio_sg_c5h14no4p, 250, "Header", TRUE); add(.eieio_sg_lpc_gly_c, 100, "Gly-C"); add(.eieio_sg_lpc_gly_o, 100, "Gly-O")
	if (name == "[M+H]+") {
		add(lipid_mass - .eieio_sg_h2o, 200, "C5H14NO+"); add(.eieio_sg_c5h13no, 200, "C5H14NO+"); add(.eieio_sg_c5h14no4p - .eieio_sg_h2o, 50, "Header - H2O")
	} else if (name == "[M+Na]+") {
		add(lipid_mass - .eieio_sg_c3h9n, 800, "Precursor -C3H9N"); add(lipid_mass - .eieio_sg_c5h11n, 200, "Precursor -C5H11N")
	}
	peaks
}

.eieio_sg_get_lpe_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE, direct = FALSE) {
		mz <- if (direct) value else .eieio_sg_ion_mz(value, adduct)
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	}
	add(lipid_mass - .eieio_sg_h2o, 200, "Precursor - H2O")
	add(lipid_mass - .eieio_sg_c2h8no4p + .eieio_sg_proton, 250, "Precursor -C2H8NO4P", TRUE, TRUE)
	add(.eieio_sg_c2h8no4p, 200, "Header"); add(.eieio_sg_lpe_gly_c, 100, "Gly-C"); add(.eieio_sg_lpe_gly_o, 200, "Gly-O")
	if (name == "[M+H]+") { add(.eieio_sg_c3h9o6p, 100, "[C3H9O6P+H]+"); add(.eieio_sg_c3h9o6p - .eieio_sg_h2o, 100, "[C3H9O6P+H]+ -H2O"); add(.eieio_sg_lpe_gly_o - .eieio_sg_h2o, 100, "Gly-O - H2O") }
	else if (name == "[M+Na]+") add(lipid_mass - .eieio_sg_c2h5n, 250, "Precursor -C2H5N")
	peaks
}

.eieio_sg_get_lpg_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift_mass <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct)
	peaks <- .eieio_sg_empty(); add <- function(value, intensity, comment, required = FALSE) {
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(value + shift_mass, intensity, comment, "metaboliteclass", required))
	}
	if (name == "[M+H]+") add(lipid_mass, 999, "Precursor")
	if (name == "[M+NH4]+") { add(lipid_mass, 500, "Precursor"); add(lipid_mass + .eieio_sg_proton, 999, "[M+H]+") }
	if (name == "[M+Na]+") { add(.eieio_sg_c3h9o6p, 300, "Header", TRUE); add(.eieio_sg_c3h9o6p - .eieio_sg_h2o, 200, "Header - H2O", TRUE); add(.eieio_sg_lpg_gly_c, 100, "Gly-C"); add(.eieio_sg_lpg_gly_o, 100, "Gly-O"); add(lipid_mass, 999, "Precursor") }
	if (name %in% c("[M+H]+", "[M+NH4]+")) { add(lipid_mass - .eieio_sg_c3h6o2, 100, "Precursor -C3H6O2"); add(lipid_mass - .eieio_sg_h2o, 500, "[M+H]+ -H2O"); add(lipid_mass - 2 * .eieio_sg_h2o, 200, "[M+H]+ -2H2O"); add(lipid_mass - .eieio_sg_c3h6o2 - .eieio_sg_h2o, 100, "[M+H]+ -C3H6O2 -H2O"); add(.eieio_sg_lpg_gly_c, 100, "Gly-C"); add(.eieio_sg_lpg_gly_o, 100, "Gly-O"); add(.eieio_sg_c3h9o6p, 300, "Header", TRUE); add(.eieio_sg_c3h9o6p - .eieio_sg_h2o, 200, "Header - H2O", TRUE); add(lipid_mass - .eieio_sg_c3h9o6p, 500, "Precursor -C3H9O6P") }
	peaks
}

.eieio_sg_lysophospholipid_chain_peaks <- function(peaks, lipid_mass, chain, adduct, class) {
	chain_mass <- .eieio_sg_chain_mass(chain) - .eieio_sg_h
	if (isTRUE(all.equal(chain_mass, 0))) return(peaks)
	label <- .eieio_sg_chain_label(chain)
	intensity <- switch(class, LPE = 50, LPG = 100, LPI = 100, 100)
	loss_intensity <- switch(class, LPG = 200, LPI = 150, 100)
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(chain_mass + .eieio_sg_proton, intensity,
		paste0(label, if (class == "LPI") " acyl+" else " acyl"), "acylchain"))
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass, adduct), loss_intensity,
		paste0("-", label), "acylchain"))
	.eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_h2o, adduct),
		if (class == "LPE") 50 else 100, paste0("-", label, "-O"), "acylchain"))
}

.eieio_sg_lysophospholipid_position <- function(lipid_mass, chain, adduct, class) {
	chain_mass <- .eieio_sg_chain_mass(chain)
	shift <- if (.eieio_sg_adduct_name(adduct) == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_adduct_mass(adduct)
	base <- lipid_mass - chain_mass - .eieio_sg_h - .eieio_sg_o - .eieio_sg_ch2 + shift
	if (class == "LPG") {
		return(.eieio_sg_rbind(
			.eieio_sg_peak(base, 100, "-CH2(Sn1)", "snposition"),
			.eieio_sg_peak(base - .eieio_sg_h2o, 100, "-H2O -CH2(Sn1)", "snposition")
		))
	}
	.eieio_sg_peak(base, 100, "-CH2(Sn1)", "snposition")
}

.eieio_sg_pe_chain_peaks <- function(peaks, lipid_mass, chain, adduct) {
	chain_mass <- .eieio_sg_chain_mass(chain) - .eieio_sg_h
	if (isTRUE(all.equal(chain_mass, 0))) return(peaks)
	label <- .eieio_sg_chain_label(chain)
	name <- .eieio_sg_adduct_name(adduct)
	if (name == "[M+H]+") {
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(chain_mass + .eieio_sg_proton, 100, paste0(label, " acyl+"), "acylchain"))
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass, adduct), 100, paste0("-", label), "acylchain"))
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_c2h8no4p, adduct), 100, paste0("-Header -", label), "acylchain"))
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_o - .eieio_sg_h, adduct), 100, paste0("-", label, "-O"), "acylchain"))
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_c2h8no4p - .eieio_sg_h2o, adduct), 100, paste0("-Header -", label, "-O"), "acylchain"))
	} else if (name == "[M+Na]+") {
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - 2 * .eieio_sg_h, adduct), 100, paste0("-", label), "acylchain"))
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_h2o, adduct), 100, paste0("-", label, "-O"), "acylchain"))
	}
	peaks
}

.eieio_sg_pg_chain_peaks <- function(peaks, lipid_mass, chain, adduct) {
	chain_mass <- .eieio_sg_chain_mass(chain) - .eieio_sg_h
	if (isTRUE(all.equal(chain_mass, 0))) return(peaks)
	label <- .eieio_sg_chain_label(chain)
	shift <- if (.eieio_sg_adduct_name(adduct) == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_adduct_mass(adduct)
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(chain_mass + .eieio_sg_proton, 100, paste0(label, " acyl"), "acylchain"))
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass - chain_mass + shift, 100, paste0("-", label), "acylchain"))
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass - chain_mass - .eieio_sg_h2o + shift, 100, paste0("-", label, "-H2O"), "acylchain"))
	if (.eieio_sg_adduct_name(adduct) == "[M+NH4]+") {
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass - chain_mass - .eieio_sg_c3h9o6p + shift, 200, paste0("-Header -", label), "acylchain"))
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass - chain_mass - .eieio_sg_c3h9o6p - .eieio_sg_h2o + shift, 100, paste0("-Header -", label, "-H2O"), "acylchain"))
	}
	peaks
}

.eieio_sg_get_lpi_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift_mass <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct); ion_mass <- lipid_mass + shift_mass
	peaks <- .eieio_sg_empty(); add <- function(value, intensity, comment, required = FALSE, direct = FALSE) {
		mz <- if (direct) value else value + shift_mass
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	}
	add(.eieio_sg_c6h13o9p, 300, "Header", TRUE); add(.eieio_sg_c6h13o9p - .eieio_sg_h2o, 300, "Header -H2O"); add(.eieio_sg_lpi_gly_c, 100, "Gly-C"); add(.eieio_sg_lpi_gly_o, 200, "Gly-O"); add(.eieio_sg_c3h9o6p, 300, "C3H9O6P"); add(.eieio_sg_c3h9o6p - .eieio_sg_h2o, 300, "C3H9O6P - H2O"); add(ion_mass - .eieio_sg_c6h13o9p, 500, "[M+H]+ -Header", direct = TRUE); add(ion_mass - .eieio_sg_c6h10o5, 100, "[M+H]+ -C6H10O5", direct = TRUE); add(ion_mass - .eieio_sg_c6h10o5 - .eieio_sg_h2o, 150, "[M+H]+ -C6H12O6", direct = TRUE); add(ion_mass - .eieio_sg_h2o, 800, "[M+H]+ -H2O", direct = TRUE)
	if (name == "[M+H]+") add(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", direct = TRUE)
	else if (name == "[M+NH4]+") { add(.eieio_sg_ion_mz(lipid_mass, adduct), 800, "Precursor", direct = TRUE); add(ion_mass, 999, "[M+H]+", direct = TRUE) }
	peaks
}

.eieio_sg_get_pg_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE, direct = FALSE) { mz <- if (direct) value else value + shift; peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required)) }
	add(lipid_mass - .eieio_sg_c3h9o6p, 500, "Precursor -C3H9O6P", TRUE); add(lipid_mass - .eieio_sg_c3h6o2, 100, "Precursor -C3H6O2"); add(.eieio_sg_c3h9o6p, 100, "Header"); add(.eieio_sg_pg_gly_c, 400, "Gly-C"); add(.eieio_sg_pg_gly_o, 350, "Gly-O")
	if (name == "[M+NH4]+") { add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O", FALSE, FALSE); add(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", FALSE, TRUE) }
	else if (name == "[M+Na]+") add(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 50, "Precursor2+", FALSE, TRUE)
	peaks
}

.eieio_sg_get_pi_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE, direct = FALSE) { mz <- if (direct) value else value + shift; peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required)) }
	add(.eieio_sg_c6h13o9p, 500, "Header"); add(.eieio_sg_lpi_gly_c, 250, "Gly-C"); add(.eieio_sg_lpi_gly_o, 200, "Gly-O")
	if (name == "[M+Na]+") add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c6h10o5, adduct), 400, "Precursor -C6H10O5", TRUE, TRUE)
	else { add(lipid_mass - .eieio_sg_c6h13o9p + .eieio_sg_proton, 500, "Precursor -Header", FALSE, TRUE); add(lipid_mass - .eieio_sg_c6h10o5, 100, "Precursor -C6H10O5", TRUE) }
	if (name == "[M+NH4]+") add(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", FALSE, TRUE)
	peaks
}

.eieio_sg_get_ps_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(value, intensity, comment, required = FALSE) { peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(value, adduct), intensity, comment, "metaboliteclass", required)) }
	add(lipid_mass - .eieio_sg_c3h8no6p, 500, "Precursor -C3H8NO6P", TRUE); add(.eieio_sg_ps_gly_c, 150, "Gly-C"); add(.eieio_sg_ps_gly_o, 100, "Gly-O")
	if (name == "[M+Na]+") { add(lipid_mass - .eieio_sg_c3h8no6p + .eieio_sg_proton, 200, "Precursor -C3H8NO6P -Na+"); add(lipid_mass - .eieio_sg_c3h5no2, 200, "Precursor -C3H5NO2"); add(.eieio_sg_c3h8no6p, 500, "Header"); add(lipid_mass - .eieio_sg_cho2, 50, "Precursor -CHO2") }
	else { add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O"); add(lipid_mass - .eieio_sg_cho2, 200, "Precursor -CHO2"); add(.eieio_sg_c3h8no6p, 100, "Header") }
	peaks
}

.eieio_sg_get_pa_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(value, intensity, comment, required = FALSE, direct = FALSE) { mz <- if (direct) value else value + shift; peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required)) }
	add(lipid_mass - .eieio_sg_h3po4, 500, "- Header"); add(.eieio_sg_c3h9o6p, 150, "C3H9O6P", TRUE); add(.eieio_sg_pg_gly_c, 100, "Gly-C"); add(.eieio_sg_pg_gly_o, 100, "Gly-O")
	if (name == "[M+NH4]+") add(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", FALSE, TRUE)
	peaks
}

.eieio_sg_get_cl_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_h2o, adduct), 100, "Precursor -H2O", "metaboliteclass", TRUE))
	if (name == "[M+NH4]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", "metaboliteclass"))
	peaks
}

.eieio_sg_get_bmp_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE, direct = FALSE) { mz <- if (direct) value else value + shift; peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required)) }
	add(lipid_mass - .eieio_sg_c3h9o6p, 150, "Precursor -C3H9O6P", direct = TRUE)
	add(.eieio_sg_c3h9o6p, 200, "C3H9O6P"); add(.eieio_sg_c3h9o6p - .eieio_sg_h2o, 200, "C3H9O6P - H2O")
	if (name == "[M+NH4]+") add(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", direct = TRUE)
	if (name %in% c("[M+H]+", "[M+NH4]+")) { add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O", direct = TRUE); add(lipid_mass - 2 * .eieio_sg_h2o, 100, "Precursor -2H2O", direct = TRUE) }
	peaks
}

.eieio_sg_get_hbmp_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_h2o, adduct), 100, "Precursor -H2O", "metaboliteclass", TRUE))
	if (name == "[M+NH4]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton, 100, "[M+H]+", "metaboliteclass"))
	peaks
}

.eieio_sg_get_dgta_spectrum <- function(lipid_mass, adduct, class = "DGTA") {
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(value, adduct), intensity, comment, "metaboliteclass", required))
	add(lipid_mass - .eieio_sg_cho2, 200, "Precursor - CO2")
	if (class == "DGTA") { add(.eieio_sg_c7h13no2, 100, "Header", TRUE); add(.eieio_sg_c7h13no2 + .eieio_sg_o, 100, "Header +O"); add(.eieio_sg_c7h13no2 + .eieio_sg_h2o, 100, "Header +H2O") }
	else { add(.eieio_sg_c8h16no3, 200, "C8H16NO3", TRUE); add(.eieio_sg_c6h11no2, 200, "Header - CH2", TRUE) }
	add(.eieio_sg_lps_gly_c, if (class == "DGTA") 200 else 150, "Gly-C", TRUE); add(.eieio_sg_lps_gly_o, if (class == "DGTA") 200 else 150, "Gly-O", TRUE)
	peaks
}

.eieio_sg_get_ceramide_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	if (name == "[M+Na]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_ch3o, 150, "Precursor-CH3O", "metaboliteclass", TRUE))
	else { peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_h2o, 200, "Precursor-H2O", "metaboliteclass", TRUE)); peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton - .eieio_sg_ch6o2, 100, "[M+H]+ -CH6O2", "metaboliteclass", TRUE)); if (name == "[M+H]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - 2 * .eieio_sg_h2o, 100, "Precursor-2H2O", "metaboliteclass")) }
	peaks
}

.eieio_sg_get_hexcer_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(.eieio_sg_c6h10o5 + .eieio_sg_c2h5no, adduct), 400, "Hex + C2H5NO", TRUE)
	if (name == "[M+Na]+") { add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c5h10o4, 150, "Precursor-C5H10O4"); add(.eieio_sg_ion_mz(.eieio_sg_c5h10o4, adduct) - .eieio_sg_h, 100, "[C5H9O4+adduct]+"); add(.eieio_sg_ion_mz(.eieio_sg_c6h10o5 + .eieio_sg_h2o, adduct), 100, "[Hex+adduct+H]+"); add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c6h10o5 - .eieio_sg_h2o, 100, "Precursor-Hex-H2O", TRUE) }
	else { add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_h2o, 200, "Precursor-H2O"); add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c6h10o5, 300, "Precursor-Hex"); add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c6h10o5 - .eieio_sg_h2o, 400, "Precursor-Hex-H2O", TRUE) }
	peaks
}

.eieio_sg_get_etherpc_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 40, "Precursor2+"); add(.eieio_sg_ion_mz(.eieio_sg_lpc_gly_c, adduct), 400, "Gly-C", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_lpc_gly_o, adduct), 400, "Gly-O", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_c5h14no4p, adduct), 500, "Header", TRUE)
	if (name == "[M+Na]+") { add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c3h9n, adduct), 500, "Precursor -C3H9N"); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c5h11n, adduct), 200, "Precursor -C5H11N") }
	peaks
}

.eieio_sg_get_etherpe_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(.eieio_sg_c2h8no4p, adduct), 100, "Header", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_lpe_gly_c, adduct), 300, "Gly-C", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_lpe_gly_o, adduct), 300, "Gly-O", TRUE); add(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 100, "[Precursor]2+")
	if (name == "[M+Na]+") add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c2h5n, adduct), 150, "Precursor -C2H5N") else add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c2h8no4p, adduct), 500, "Precursor -C2H8NO4P")
	peaks
}

.eieio_sg_get_gm3_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); shift <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct); h2o_loss <- lipid_mass + shift - .eieio_sg_h2o
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(h2o_loss - .eieio_sg_h, 100, "[M-H2O-H+H]+"); add(h2o_loss - .eieio_sg_c11h15no7 - .eieio_sg_c6h10o5 - .eieio_sg_h2o, 150, "[M-H2O-C17H27NO13+adduct]+"); add(.eieio_sg_c11h15no7 + 2 * .eieio_sg_c6h10o5 + 2 * .eieio_sg_h2o + .eieio_sg_c2h3n + shift, 300, "[C23H35NO17+C2H3N+2H2O+adduct]+")
	if (name %in% c("[M+H]+", "[M+NH4]+")) { add(h2o_loss - .eieio_sg_c11h15no7, 100, "[M-H2O-C11H15NO7+H]+"); add(h2o_loss - .eieio_sg_c11h15no7 - .eieio_sg_c6h10o5, 150, "[M-H2O-C17H25NO12+H]+"); add(h2o_loss - .eieio_sg_c11h15no7 - 2 * .eieio_sg_c6h10o5 - .eieio_sg_h2o, 300, "[M-H2O-C23H37NO18+H]+", TRUE); add(.eieio_sg_c11h15no7 + .eieio_sg_c6h10o5 + .eieio_sg_h2o + shift, 300, "[C17H27NO13+H]+", TRUE); add(.eieio_sg_c11h15no7 + shift, 300, "[C11H15NO7+H]+") }
	if (name == "[M+NH4]+") add(lipid_mass + shift, 500, "[M+H]+")
	peaks
}

.eieio_sg_get_shexcer_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_so3 - .eieio_sg_h2o, adduct), 200, "[M-H2SO4+H]+", TRUE)
	if (name == "[M+H]+") { add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_so3, adduct), 100, "[M-SO3+H]+"); add(.eieio_sg_ion_mz(.eieio_sg_c8h16no6, adduct), 150, "[C8H16NO6+H]+", TRUE); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_so3 - .eieio_sg_c6h10o5 - .eieio_sg_h2o, adduct), 250, "[M-C6H12O9S+H]+", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_c8h16no6 + .eieio_sg_so3, adduct), 150, "[C8H16NO9S+H]+") }
	else if (name == "[M+Na]+") { add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_so3, adduct), 700, "[M-SO3+Na]+"); add(.eieio_sg_ion_mz(.eieio_sg_c6h10o5 + .eieio_sg_h2o, adduct), 400, "[Hex+Na]+"); add(.eieio_sg_ion_mz(.eieio_sg_c8h16no6, adduct), 500, "[C8H16NO6+Na]+", TRUE); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_so3 - .eieio_sg_c6h10o5 - .eieio_sg_h2o, adduct), 150, "[M-C6H12O9S+Na]+", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_c8h16no6 + .eieio_sg_so3, adduct), 100, "[C8H16NO9S+Na]+") }
	peaks
}

.eieio_sg_get_sm_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(.eieio_sg_c5h14no4p, adduct), 500, "C5H14NO4P (Header)", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_c5h14no4p + .eieio_sg_c2h3no - .eieio_sg_o, adduct), 400, "C5H14NO4P(Header) + C2H3N", TRUE)
	if (name == "[M+Na]+") { add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c3h9n, 300, "Precursor-C3H9N", TRUE); add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c5h11n, 150, "Precursor-C5H11N", TRUE); add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_c5h14no4p, 150, "Precursor-Header") }
	peaks
}

.eieio_sg_get_car_spectrum <- function(lipid_mass, adduct) { peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required)); add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_cho2, 150, "Precursor-CHO2"); add(.eieio_sg_ion_mz(.eieio_sg_c7h15no3, adduct), 150, "Header", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_c7h15no3 - .eieio_sg_h2o, adduct), 150, "Header-H2O", TRUE); peaks }

.eieio_sg_get_dg_spectrum <- function(lipid_mass, adduct) { name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); add <- function(mz, intensity, comment) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass")); if (name %in% c("[M+H]+", "[M+NH4]+")) { add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_h2o, 100, "Precursor-H2O"); if (name == "[M+NH4]+") { add(lipid_mass + .eieio_sg_proton, 150, "[M+H]+"); add(lipid_mass + .eieio_sg_proton - .eieio_sg_h2o, 150, "[M+H]+ -H2O") } }; peaks }

.eieio_sg_get_tg_spectrum <- function(lipid_mass, adduct) { name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor"); peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 50, "[Precursor]2+", "precursor")); if (name %in% c("[M+H]+", "[M+NH4]+")) { peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_h2o, 100, "Precursor-H2O", "metaboliteclass")); if (name == "[M+NH4]+") { peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton, 150, "[M+H]+", "metaboliteclass")); peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton - .eieio_sg_h2o, 150, "[M+H]+ -H2O", "metaboliteclass")) } }; peaks }

.eieio_sg_get_dmedfa_spectrum <- function(lipid_mass, adduct, oxidized = FALSE, hfa = FALSE) {
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor(DMED derv.)", "precursor")
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c2nh7, adduct), if (hfa) 50 else 200, "Precursor - C2NH7", "metaboliteclass"))
	if (oxidized && !hfa) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c2nh7 - .eieio_sg_h2o, adduct), 50, "Precursor - C2NH7 - H2O", "metaboliteclass"))
	peaks
}

.eieio_sg_get_pc_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(.eieio_sg_c5h14no4p, adduct), 300, "Header", TRUE); add(.eieio_sg_ion_mz(.eieio_sg_lpc_gly_c, adduct), 150, "Gly-C"); add(.eieio_sg_ion_mz(.eieio_sg_lpc_gly_o, adduct), 100, "Gly-O"); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c5h14no4p - 2 * .eieio_sg_h, adduct), 25, "Precursor-Header"); add(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 150, "[Precursor]2+")
	if (name == "[M+Na]+") { add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c5h14no4p, adduct), 300, "Precursor-Header"); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c3h9n, adduct), 200, "Precursor-C3H9N") }
	peaks
}

.eieio_sg_get_phyto_cer_spectrum <- function(lipid_mass, adduct, abp = FALSE) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	if (name == "[M+Na]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_ch3o, 150, "Precursor-CH3O", "metaboliteclass", TRUE))
	else { peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_h2o, adduct), 200, "Precursor-H2O", "metaboliteclass", TRUE)); if (!abp) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton - .eieio_sg_ch6o2, 100, "[M+H]+ -CH6O2", "metaboliteclass", TRUE)); if (name == "[M+H]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - 2 * .eieio_sg_h2o, 100, "Precursor-2H2O", "metaboliteclass")) }
	peaks
}

.eieio_sg_get_pe_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); header <- .eieio_sg_c2h8no4p
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(lipid_mass - header, adduct), 500, "Precursor -C2H8NO4P", TRUE); add(.eieio_sg_ion_mz(header, adduct), 100, "Header"); add(.eieio_sg_ion_mz(.eieio_sg_lpe_gly_c, adduct), 400, "Gly-C"); add(.eieio_sg_ion_mz(.eieio_sg_lpe_gly_o, adduct), 300, "Gly-O"); add(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 100, "[Precursor]2+")
	if (name == "[M+Na]+") { add(lipid_mass - header + .eieio_sg_proton, 250, "Precursor -C2H8NO4P"); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c2h5n + .eieio_sg_ch2, adduct), 100, "Precursor -CH4N"); add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_c2h5n, adduct), 300, "Precursor -C2H5N") }
	peaks
}

.eieio_sg_d <- 2.01410177812
.eieio_sg_mass_d <- function(carbon = 0, hydrogen = 0, nitrogen = 0,
								 oxygen = 0, phosphorus = 0, deuterium = 0) {
	.eieio_sg_mass(carbon, hydrogen, nitrogen, oxygen, phosphorus) +
		deuterium * (.eieio_sg_d - .eieio_sg_h)
}

.eieio_sg_isotope_pl_spectrum <- function(class, lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct)
	shift <- if (name == "[M+NH4]+") .eieio_sg_proton else .eieio_sg_ion_mz(0, adduct)
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(value, intensity, comment, required = FALSE, direct = FALSE) {
		mz <- if (direct) value else value + shift
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	}
	if (class == "PC_d5") {
		header <- .eieio_sg_mass_d(5, 14, 1, 4, 1, 5)
		glyc <- .eieio_sg_mass_d(8, 13, 1, 4, 1, 5); glyo <- .eieio_sg_mass_d(7, 13, 1, 5, 1, 3)
		add(header, 300, "Header", TRUE); add(glyc, 150, "Gly-C"); add(glyo, 100, "Gly-O"); add(lipid_mass - header - 2 * .eieio_sg_h, 25, "Precursor-Header")
		if (name == "[M+Na]+") { add(lipid_mass - header, 300, "Precursor-Header"); add(lipid_mass - .eieio_sg_mass(3, 9, 1), 200, "Precursor-C3H9N") }
	} else if (class == "PE_d5") {
		header <- .eieio_sg_c2h8no4p; glyc <- .eieio_sg_mass_d(5, 7, 1, 4, 1, 5); glyo <- .eieio_sg_mass_d(4, 7, 1, 5, 1, 3)
		add(lipid_mass - header + shift, 500, "Precursor -C2H8NO4P", TRUE, TRUE); add(header, 100, "Header"); add(glyc, 400, "Gly-C"); add(glyo, 300, "Gly-O"); add(lipid_mass, 100, "[Precursor]2+", FALSE, TRUE)
		if (name == "[M+Na]+") { add(lipid_mass - header + .eieio_sg_proton, 250, "Precursor -C2H8NO4P", FALSE, TRUE); add(lipid_mass - .eieio_sg_mass(1, 5, 1) + .eieio_sg_h, 100, "Precursor -CH4N"); add(lipid_mass - .eieio_sg_mass(1, 5, 1) - .eieio_sg_mass(1, 0), 300, "Precursor -C2H5N") }
	} else if (class == "PG_d5") {
		header <- .eieio_sg_c3h9o6p; glyc <- .eieio_sg_mass_d(6, 8, 0, 6, 1, 5); glyo <- .eieio_sg_mass_d(5, 8, 0, 7, 1, 3)
		add(lipid_mass - header + shift, 500, "Precursor -C3H9O6P", TRUE, TRUE); add(lipid_mass - .eieio_sg_mass(3, 6, 0, 2) + shift, 100, "Precursor -C3H6O2", FALSE, TRUE); add(header, 100, "Header"); add(glyc, 400, "Gly-C"); add(glyo, 350, "Gly-O")
		if (name == "[M+NH4]+") { add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O", FALSE, TRUE); add(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", FALSE, TRUE) } else if (name == "[M+Na]+") add(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 50, "Precursor2+", FALSE, TRUE)
	} else if (class == "PI_d5") {
		header <- .eieio_sg_mass_d(6, 13, 0, 9, 1, 5); glyc <- .eieio_sg_mass_d(9, 12, 0, 9, 1, 5); glyo <- .eieio_sg_mass_d(8, 12, 0, 10, 1, 3)
		add(header, 500, "Header"); add(glyc, 250, "Gly-C"); add(glyo, 200, "Gly-O")
		if (name == "[M+Na]+") add(.eieio_sg_ion_mz(lipid_mass - .eieio_sg_mass(6, 10, 0, 5), adduct), 400, "Precursor -C6H10O5", TRUE, TRUE)
		else { add(lipid_mass - header + .eieio_sg_proton, 500, "Precursor -Header", FALSE, TRUE); add(lipid_mass - .eieio_sg_mass(6, 10, 0, 5), 100, "Precursor -C6H10O5", TRUE) }
		if (name == "[M+NH4]+") add(lipid_mass + .eieio_sg_proton, 200, "[M+H]+", FALSE, TRUE)
	} else if (class == "PS_d5") {
		header <- .eieio_sg_mass_d(3, 8, 1, 6, 1, 5); glyc <- .eieio_sg_mass_d(6, 7, 1, 6, 1, 5); glyo <- .eieio_sg_mass_d(5, 7, 1, 7, 1, 3)
		add(lipid_mass - header + shift, 500, "Precursor -C3H8NO6P", TRUE, TRUE); add(glyc, 150, "Gly-C"); add(glyo, 100, "Gly-O")
		if (name == "[M+Na]+") { add(lipid_mass - header + .eieio_sg_proton, 200, "Precursor -C3H8NO6P -Na+", FALSE, TRUE); add(lipid_mass - .eieio_sg_mass(3, 5, 1, 2), 200, "Precursor -C3H5NO2"); add(header, 500, "Header"); add(lipid_mass - .eieio_sg_cho2, 50, "Precursor -CHO2") } else { add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O"); add(lipid_mass - .eieio_sg_cho2, 200, "Precursor -CHO2"); add(header, 100, "Header") }
	} else if (class %in% c("LPC_d5", "LPE_d5", "LPG_d5", "LPI_d5", "LPS_d5")) {
		base <- sub("_d5$", "", class)
		if (base == "LPC") { header <- .eieio_sg_mass_d(5, 14, 1, 4, 1, 5); glyc <- .eieio_sg_mass_d(8, 13, 1, 4, 1, 5); glyo <- .eieio_sg_mass_d(7, 13, 1, 5, 1, 3); add(header, 250, "Header", TRUE); add(glyc, 100, "Gly-C"); add(glyo, 100, "Gly-O"); if (name == "[M+H]+") { add(lipid_mass - .eieio_sg_h2o, 200, "Precursor -H2O"); add(.eieio_sg_mass_d(5, 13, 1, 1, 0, 5), 200, "C5H14NO+"); add(header - .eieio_sg_h2o, 50, "Header - H2O") } else if (name == "[M+Na]+") { add(lipid_mass - .eieio_sg_c3h9n, 800, "Precursor -C3H9N"); add(lipid_mass - .eieio_sg_c5h11n, 200, "Precursor -C5H11N") }
		} else if (base == "LPE") { header <- .eieio_sg_c2h8no4p; glyc <- .eieio_sg_mass_d(5, 7, 1, 4, 1, 5); glyo <- .eieio_sg_mass_d(4, 7, 1, 5, 1, 3); add(lipid_mass - .eieio_sg_h2o, 200, "Precursor - H2O"); add(lipid_mass - header + .eieio_sg_proton, 250, "Precursor -C2H8NO4P", TRUE, TRUE); add(header, 200, "Header"); add(glyc, 100, "Gly-C"); add(glyo, 200, "Gly-O"); if (name == "[M+H]+") { add(.eieio_sg_mass_d(3, 4, 0, 6, 1, 5), 100, "C3H4D5O6P"); add(.eieio_sg_mass_d(3, 4, 0, 6, 1, 5) - .eieio_sg_h2o, 100, "C3H4D5O6P - H2O"); add(glyo - .eieio_sg_h2o, 100, "Gly-O - H2O") } else if (name == "[M+Na]+") add(lipid_mass - .eieio_sg_c2h5n, 250, "Precursor -C2H5N")
		} else if (base == "LPG") { header <- .eieio_sg_c3h9o6p; glyc <- .eieio_sg_mass_d(6, 8, 0, 6, 1, 5); glyo <- .eieio_sg_mass_d(5, 8, 0, 7, 1, 3); if (name == "[M+Na]+") { add(header, 300, "Header", TRUE); add(header - .eieio_sg_h2o, 200, "Header - H2O", TRUE); add(glyc, 100, "Gly-C"); add(glyo, 100, "Gly-O") } else { add(lipid_mass - .eieio_sg_c3h6o2, 100, "Precursor -C3H6O2", FALSE, TRUE); add(lipid_mass - .eieio_sg_h2o, 500, "[M+H]+ -H2O", FALSE, TRUE); add(lipid_mass - 2*.eieio_sg_h2o, 200, "[M+H]+ -2H2O", FALSE, TRUE); add(lipid_mass - .eieio_sg_c3h6o2 - .eieio_sg_h2o, 100, "[M+H]+ -C3H6O2", FALSE, TRUE); add(glyc, 100, "Gly-C"); add(glyo, 100, "Gly-O"); add(header, 300, "Header", TRUE); add(header - .eieio_sg_h2o, 200, "Header - H2O", TRUE); add(lipid_mass - header, 500, "Precursor -C3H9O6P", FALSE, TRUE) } }
		else if (base == "LPI") { header <- .eieio_sg_mass_d(6, 13, 0, 9, 1, 5); glyc <- .eieio_sg_mass_d(9, 12, 0, 9, 1, 5); glyo <- .eieio_sg_mass_d(8, 12, 0, 10, 1, 3); add(header, 300, "Header", TRUE); add(header - .eieio_sg_h2o, 300, "Header -H2O"); add(glyc, 100, "Gly-C"); add(glyo, 200, "Gly-O"); add(.eieio_sg_mass_d(3, 4, 0, 6, 1, 5), 300, "C3H9O6P"); add(.eieio_sg_mass_d(3, 4, 0, 6, 1, 5) - .eieio_sg_h2o, 300, "C3H9O6P - H2O"); add(lipid_mass + shift - .eieio_sg_c6h13o9p, 500, "[M+H]+ -Header", FALSE, TRUE); add(lipid_mass + shift - .eieio_sg_c6h10o5, 100, "[M+H]+ -C6H10O5", FALSE, TRUE); add(lipid_mass + shift - .eieio_sg_c6h10o5 - .eieio_sg_h2o, 150, "[M+H]+ -C6H12O6", FALSE, TRUE); add(lipid_mass + shift - .eieio_sg_h2o, 800, "[M+H]+ -H2O", FALSE, TRUE) }
			else { header <- .eieio_sg_mass_d(3, 8, 1, 6, 1, 5); glyc <- .eieio_sg_mass_d(6, 7, 1, 6, 1, 5); glyo <- .eieio_sg_mass_d(5, 7, 1, 7, 1, 3); add(header, 100, "Header", TRUE); add(glyc, 100, "Gly-C"); add(glyo, 100, "Gly-O"); if (name == "[M+H]+") { add(lipid_mass - header + shift, 500, "Precursor -C3H8NO6P", FALSE, TRUE); add(lipid_mass - .eieio_sg_cho2, 200, "Precursor -CHO2"); add(lipid_mass - .eieio_sg_h2o, 100, "Precursor -H2O"); add(.eieio_sg_mass_d(3, 4, 0, 6, 1, 5), 100, "C3H4D5O6P"); add(.eieio_sg_mass_d(3, 4, 0, 6, 1, 5) - .eieio_sg_h2o, 100, "C3H4D5O6P - H2O") } else if (name == "[M+Na]+") add(lipid_mass - .eieio_sg_mass(3, 5, 1, 2) + shift, 200, "Precursor -C3H5NO2", FALSE, TRUE) }
	}
	peaks
}

.eieio_sg_get_dg_d5_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	if (name == "[M+NH4]+") { peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton, 150, "[M+H]+", "metaboliteclass")); peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton - .eieio_sg_h2o, 150, "[M+H]+ -H2O", "metaboliteclass")) }
	else if (name == "[M+H]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_mass(0, 0, 0, 1) - 2 * (.eieio_sg_d - .eieio_sg_h), 100, "Precursor-D2O", "metaboliteclass"))
	peaks
}

.eieio_sg_get_tg_d5_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) / 2, 50, "[Precursor]2+", "precursor"))
	if (name == "[M+NH4]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(lipid_mass + .eieio_sg_proton, 150, "[M+H]+", "metaboliteclass"))
	peaks
}

.eieio_sg_get_ce_d7_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); skeleton <- .eieio_sg_mass(27, 46, 0, 1) + 7 * (.eieio_sg_d - .eieio_sg_h) - .eieio_sg_h2o
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	if (name %in% c("[M+NH4]+", "[M+Na]+")) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(skeleton + if (name == "[M+NH4]+") .eieio_sg_proton else 0, 500, "skelton", "metaboliteclass", TRUE))
	if (name == "[M+NH4]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_h2o, 50, "Precursor -H2O", "metaboliteclass"))
	peaks
}

.eieio_sg_get_sm_d9_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); header <- .eieio_sg_mass_d(5, 5, 1, 4, 1, 9); peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass, adduct), 999, "Precursor", "precursor")
	add <- function(mz, intensity, comment, required = FALSE) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mz, intensity, comment, "metaboliteclass", required))
	add(.eieio_sg_ion_mz(header, adduct), 500, "Header", TRUE); add(.eieio_sg_ion_mz(header + .eieio_sg_mass(2, 3, 1, 1) - .eieio_sg_o, adduct), 400, "Header+C2H3N", TRUE)
	if (name == "[M+Na]+") { add(.eieio_sg_ion_mz(lipid_mass, adduct) - .eieio_sg_mass_d(3, 0, 1, 0, 0, 9), 300, "Precursor-C3D9N", TRUE); add(.eieio_sg_ion_mz(lipid_mass, adduct) - header, 150, "Precursor-Header") }
	else if (name == "[M+H]+") add(.eieio_sg_ion_mz(header + .eieio_sg_mass(2, 3, 1, 1) + .eieio_sg_mass(1, 0), adduct), 150, "Header+C3H3NO")
	peaks
}

.eieio_sg_get_cer_ns_d7_spectrum <- function(lipid_mass, adduct) {
	name <- .eieio_sg_adduct_name(adduct); balance <- 7 * (.eieio_sg_d - .eieio_sg_h); labeled <- lipid_mass + balance
	peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(labeled, adduct), 999, "Precursor", "precursor")
	if (name == "[M+Na]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(labeled, adduct) - .eieio_sg_ch3o, 150, "Precursor-CH3O", "metaboliteclass", TRUE))
	else { peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(labeled, adduct) - .eieio_sg_h2o, 200, "Precursor-H2O", "metaboliteclass", TRUE)); peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(labeled + .eieio_sg_proton - .eieio_sg_ch6o2, 100, "[M+H]+ -CH6O2", "metaboliteclass")); if (name == "[M+H]+") peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(labeled, adduct) - 2 * .eieio_sg_h2o, 100, "Precursor-2H2O", "metaboliteclass")) }
	peaks
}

.eieio_sg_isotope_chain_peaks <- function(peaks, lipid_mass, chain, adduct, class, abundance = 50, position = 1L) {
	label <- .eieio_sg_chain_label(chain)
	chain_mass <- .eieio_sg_chain_mass(chain)
	chain_loss <- chain_mass - .eieio_sg_h
	cd2 <- .eieio_sg_mass_d(1, 0, 0, 0, 0, 2)
	add <- function(mass, intensity, comment, type = "acylchain", required = FALSE) {
		peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mass, intensity, comment, type, required))
	}
	if (class == "TG_d5") {
		if (position == 2L) {
			if (.eieio_sg_adduct_name(adduct) == "[M+Na]+") add(.eieio_sg_ion_mz(chain_mass + .eieio_sg_mass_d(2, 0, 0, 1, 0, 3), adduct), 100, "Sn2 diagnostics", "snposition", TRUE)
			else {
				add(chain_mass + .eieio_sg_electron, 100, paste0(label, "+ Sn2"))
				add(chain_mass + .eieio_sg_adduct_mass(adduct) + .eieio_sg_mass_d(3, 0, 0, 2, 0, 5), 100, paste0(label, "+C3D5O2 Sn2"))
				add(chain_mass + .eieio_sg_adduct_mass(adduct) + .eieio_sg_mass_d(3, 0, 0, 2, 0, 5) - .eieio_sg_h2o, 100, paste0(label, "+C3D3O Sn2"))
				add(.eieio_sg_ion_mz(lipid_mass - chain_loss + 2 * .eieio_sg_h, adduct), 50, paste0("-", label, " Sn2"))
				add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_o, adduct), 200, paste0("-", label, "-O Sn2"))
			}
		} else add(.eieio_sg_ion_mz(chain_mass + .eieio_sg_mass(2, 3, 0, 1), adduct), 100, paste0(label, "+C2D3O"))
		return(peaks)
	}
	if (class == "CE_d7") {
		if (position == 1L && .eieio_sg_adduct_name(adduct) == "[M+Na]+") add(.eieio_sg_ion_mz(chain_loss + .eieio_sg_h2o, adduct), 200, paste0(label, "+O"))
		return(peaks)
	}
	if (class == "PC_d5") {
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss, adduct), 50, paste0("-", label))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_h2o, adduct), 50, paste0("-", label, "-O"))
	} else if (class == "PE_d5" && .eieio_sg_adduct_name(adduct) == "[M+H]+") {
		add(chain_loss + .eieio_sg_proton, 100, paste0(label, " acyl+"))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss, adduct), 100, paste0("-", label))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_c2h8no4p, adduct), 100, paste0("-Header -", label))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_o - .eieio_sg_h, adduct), 100, paste0("-", label, "-O"))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_c2h8no4p - .eieio_sg_h2o, adduct), 100, paste0("-Header -", label, "-O"))
	} else if (class == "PE_d5") {
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - 2 * .eieio_sg_h, adduct), 100, paste0("-", label))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_h2o, adduct), 100, paste0("-", label, "-O"))
	} else if (.eieio_sg_adduct_name(adduct) == "[M+H]+") {
		add(chain_loss + .eieio_sg_proton, 100, paste0(label, " acyl+"))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss, adduct), 100, paste0("-", label))
		add(.eieio_sg_ion_mz(lipid_mass - chain_loss - .eieio_sg_h2o, adduct), 100, paste0("-", label, "-O"))
	}
	if (class != "LPC_d5" || position == 1L) add(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_o - cd2 - if (class == "LPE_d5" && position == 1L) 0 else 0, adduct), 100, "-CD2(Sn1)", "snposition")
	peaks
}

.eieio_sg_isotope_sphingolipid_chain_peaks <- function(peaks, lipid_mass, chain, adduct, class) {
	chain_mass <- .eieio_sg_chain_mass(chain); balance <- if (class == "Cer_NS_d7") 7 * (.eieio_sg_d - .eieio_sg_h) else 0
	if (tolower(chain$type[1L]) == "sphingo") {
		loss <- chain_mass + .eieio_sg_h + balance
		if (.eieio_sg_adduct_name(adduct) != "[M+Na]+" && class == "Cer_NS_d7") {
			peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(loss + .eieio_sg_proton - .eieio_sg_h2o, 200, paste0("[sph+H]+ -H2O"), "acylchain"))
			peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(loss + .eieio_sg_proton - 2 * .eieio_sg_h2o, 500, paste0("[sph+H]+ -2H2O"), "acylchain"))
		} else if (class == "SM_d9") {
			peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(loss + .eieio_sg_proton - 2 * .eieio_sg_h2o, 100, "[sph+H]+ -Header -H2O", "acylchain"))
		}
	} else if (class == "Cer_NS_d7") {
		peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(chain_mass + .eieio_sg_h + balance, adduct) + .eieio_sg_mass(2, 3, 1, 1), 150, "[FAA+C2H2O+H]+", "acylchain"))
	}
	peaks
}

.eieio_sg_phyto_chain_peaks <- function(peaks, chain, adduct) {
	chain_mass <- .eieio_sg_chain_mass(chain)
	add <- function(mass, intensity, comment) peaks <<- .eieio_sg_rbind(peaks, .eieio_sg_peak(mass, intensity, comment, "acylchain"))
	if (tolower(chain$type[1L]) == "sphingo") {
		add(chain_mass + .eieio_sg_proton - .eieio_sg_h2o, 200, "[sph+H]+ -H2O")
		add(chain_mass + .eieio_sg_proton - 2 * .eieio_sg_h2o, 500, "[sph+H]+ -2H2O")
		add(chain_mass + .eieio_sg_proton - .eieio_sg_ch6o2, 150, "[sph+H]+ -CH6O2")
		add(chain_mass + .eieio_sg_proton, 150, "[sph+H]+ ")
	} else if (.eieio_sg_adduct_name(adduct) == "[M+H]+") {
		add(.eieio_sg_ion_mz(chain_mass - .eieio_sg_h, adduct) + .eieio_sg_c2h3no + .eieio_sg_h, 150, "[FAA+C2H2O+H]+")
		add(chain_mass + .eieio_sg_proton + .eieio_sg_c2h2n, 200, "[FAA+H]+")
		add(chain_mass + .eieio_sg_proton + .eieio_sg_c2h3no + .eieio_sg_h2o + 12, 150, "[FAA+C3H4O2+H]+")
		add(chain_mass + .eieio_sg_proton + .eieio_sg_c3h4no, 150, "[FAA+C3H3O+H]+")
		add(chain_mass + .eieio_sg_proton + .eieio_sg_c2h3no, 150, "[FAA+C2H2O+H]+")
		add(chain_mass + .eieio_sg_proton + .eieio_sg_c2h3no - .eieio_sg_h - .eieio_sg_o, 200, "[FAA+C2H+H]+")
	} else {
		add(chain_mass + .eieio_sg_proton + .eieio_sg_c2h2n, 200, "[FAA+H]+")
	}
	peaks
}

GetPGSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_pg_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetPISpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_pi_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetPSSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_ps_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetPASpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_pa_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)

GetLPCSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_lpc_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetLPESpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_lpe_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetLPGSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_lpg_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetLPISpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_lpi_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetCLSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_cl_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetBMPSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_bmp_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetHBMPSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_hbmp_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetDGTASpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_dgta_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct, "DGTA")
GetDGTSSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_dgta_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct, "DGTS")
GetCeramideSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_ceramide_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetHexCerSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_hexcer_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetEtherPCSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_etherpc_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetEtherPESpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_etherpe_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetGM3Spectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_gm3_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetSHexCerSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_shexcer_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetSMSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_sm_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetCARSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_car_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetDGSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_dg_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)
GetTGSpectrum <- function(lipid, adduct = "[M+H]+") .eieio_sg_get_tg_spectrum(as.numeric(.eieio_sg_field(lipid, c("Mass", "mass"))), adduct)

GetLPSSpectrum <- function(lipid, adduct = "[M+H]+") {
	mass <- as.numeric(.eieio_sg_field(lipid, c("Mass", "mass", "PrecursorMass"), NA_real_))[1L]
	if (is.na(mass)) stop("lipid must provide Mass")
	.eieio_sg_get_lps_spectrum(mass, adduct)
}

.eieio_sg_lps_chain_level <- function(peaks, lipid_mass, chain, adduct) {
	chain_mass <- .eieio_sg_chain_mass(chain) - .eieio_sg_h
	if (isTRUE(all.equal(chain_mass, 0))) return(peaks)
	label <- .eieio_sg_chain_label(chain)
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(chain_mass + .eieio_sg_proton, 100,
		paste0(label, " acyl+"), "acylchain"))
	peaks <- .eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass, adduct), 100,
		paste0("-", label), "acylchain"))
	.eieio_sg_rbind(peaks, .eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_h2o, adduct), 100,
		paste0("-", label, "-O"), "acylchain"))
}

.eieio_sg_lps_position <- function(lipid_mass, chain, adduct) {
	chain_mass <- .eieio_sg_chain_mass(chain) - .eieio_sg_h
	.eieio_sg_peak(.eieio_sg_ion_mz(lipid_mass - chain_mass - .eieio_sg_h2o - .eieio_sg_ch2, adduct),
		100, "-CH2(Sn1)", "snposition")
}

.eieio_sg_registered_classes <- c("EtherPC", "EtherPE", "LPC", "LPE", "LPG", "LPI", "LPS", "PA", "PC", "PE", "PG", "PI", "PS", "CL", "BMP", "HBMP", "DGTS", "LDGTA", "LDGTS", "DGTA", "DMEDFAHFA", "DMEDFA", "DMEDOxFA", "CAR", "DG", "TG", "SM", "Cer_NS", "Cer_NDS", "Cer_NP", "Cer_AS", "Cer_ADS", "Cer_AP", "Cer_BS", "Cer_BDS", "Cer_ABP", "Cer_AH", "Cer_NH", "HexCer_NS", "HexCer_NDS", "GM3", "SHexCer", "PC_d5", "PE_d5", "PG_d5", "PI_d5", "PS_d5", "LPC_d5", "LPE_d5", "LPG_d5", "LPI_d5", "LPS_d5", "DG_d5", "TG_d5", "CE_d7", "SM_d9", "Cer_NS_d7")

# Generate the active FacadeLipidSpectrumGenerator.Default spectrum. The R
# input accepts a lipid list/data frame with Mass, Name and Chains fields.
generate_eieio_spectrum <- function(lipid, adduct = "[M+H]+", abundance = NULL) {
	mass <- as.numeric(.eieio_sg_field(lipid, c("Mass", "mass", "PrecursorMass"), NA_real_))[1L]
	if (is.na(mass)) stop("lipid must provide Mass")
	class <- as.character(.eieio_sg_field(lipid, c("LipidClass", "class", "Class"), NA_character_))[1L]
	if (is.na(class) || !nzchar(class)) class <- sub("[ _].*$", "", as.character(.eieio_sg_field(lipid, c("Name", "name"), ""))[1L])
	if (!class %in% .eieio_sg_registered_classes) return(.eieio_sg_empty())
	if (!.eieio_sg_isotope_supported_adduct(class, adduct)) return(.eieio_sg_empty())
	if (!.eieio_sg_supported_adduct(class, adduct)) return(.eieio_sg_empty())
	base_class <- sub("_d[579]$", "", class)
	if (class %in% c("PC_d5", "PE_d5", "PG_d5", "PI_d5", "PS_d5", "LPC_d5", "LPE_d5", "LPG_d5", "LPI_d5", "LPS_d5")) peaks <- .eieio_sg_isotope_pl_spectrum(class, mass, adduct)
	else if (class == "DG_d5") peaks <- .eieio_sg_get_dg_d5_spectrum(mass, adduct)
	else if (class == "TG_d5") peaks <- .eieio_sg_get_tg_d5_spectrum(mass, adduct)
	else if (class == "CE_d7") peaks <- .eieio_sg_get_ce_d7_spectrum(mass, adduct)
	else if (class == "SM_d9") peaks <- .eieio_sg_get_sm_d9_spectrum(mass, adduct)
	else if (class == "Cer_NS_d7") peaks <- .eieio_sg_get_cer_ns_d7_spectrum(mass, adduct)
	else if (base_class == "LPS") peaks <- .eieio_sg_get_lps_spectrum(mass, adduct)
	else if (base_class == "LPC") peaks <- .eieio_sg_get_lpc_spectrum(mass, adduct)
	else if (base_class == "LPE") peaks <- .eieio_sg_get_lpe_spectrum(mass, adduct)
	else if (base_class == "LPG") peaks <- .eieio_sg_get_lpg_spectrum(mass, adduct)
	else if (base_class == "LPI") peaks <- .eieio_sg_get_lpi_spectrum(mass, adduct)
	else if (base_class == "PG") peaks <- .eieio_sg_get_pg_spectrum(mass, adduct)
	else if (base_class == "PI") peaks <- .eieio_sg_get_pi_spectrum(mass, adduct)
	else if (base_class == "PS") peaks <- .eieio_sg_get_ps_spectrum(mass, adduct)
	else if (base_class == "PA") peaks <- .eieio_sg_get_pa_spectrum(mass, adduct)
	else if (base_class == "CL") peaks <- .eieio_sg_get_cl_spectrum(mass, adduct)
	else if (base_class == "BMP") peaks <- .eieio_sg_get_bmp_spectrum(mass, adduct)
	else if (base_class == "HBMP") peaks <- .eieio_sg_get_hbmp_spectrum(mass, adduct)
	else if (base_class %in% c("DGTA", "LDGTA")) peaks <- .eieio_sg_get_dgta_spectrum(mass, adduct, "DGTA")
	else if (base_class %in% c("DGTS", "LDGTS")) peaks <- .eieio_sg_get_dgta_spectrum(mass, adduct, "DGTS")
	else if (base_class %in% c("Cer_NP", "Cer_AP", "Cer_ABP")) peaks <- .eieio_sg_get_phyto_cer_spectrum(mass, adduct, base_class == "Cer_ABP")
	else if (grepl("^Cer_", base_class)) peaks <- .eieio_sg_get_ceramide_spectrum(mass, adduct)
	else if (base_class == "HexCer") peaks <- .eieio_sg_get_hexcer_spectrum(mass, adduct)
	else if (base_class == "EtherPC") peaks <- .eieio_sg_get_etherpc_spectrum(mass, adduct)
	else if (base_class == "EtherPE") peaks <- .eieio_sg_get_etherpe_spectrum(mass, adduct)
	else if (base_class == "GM3") peaks <- .eieio_sg_get_gm3_spectrum(mass, adduct)
	else if (base_class == "SHexCer") peaks <- .eieio_sg_get_shexcer_spectrum(mass, adduct)
	else if (base_class == "SM") peaks <- .eieio_sg_get_sm_spectrum(mass, adduct)
	else if (base_class == "CAR") peaks <- .eieio_sg_get_car_spectrum(mass, adduct)
	else if (base_class == "DG") peaks <- .eieio_sg_get_dg_spectrum(mass, adduct)
	else if (base_class == "TG") peaks <- .eieio_sg_get_tg_spectrum(mass, adduct)
	else if (base_class == "PC") peaks <- .eieio_sg_get_pc_spectrum(mass, adduct)
	else if (base_class == "PE") peaks <- .eieio_sg_get_pe_spectrum(mass, adduct)
	else if (base_class == "DMEDFA") peaks <- .eieio_sg_get_dmedfa_spectrum(mass, adduct)
	else if (base_class == "DMEDFAHFA") peaks <- .eieio_sg_get_dmedfa_spectrum(mass, adduct, hfa = TRUE)
	else if (base_class == "DMEDOxFA") peaks <- .eieio_sg_get_dmedfa_spectrum(mass, adduct, oxidized = TRUE)
	else if (base_class %in% c("PC", "PE", "PG", "PI", "PS", "PA")) peaks <- .eieio_sg_header(base_class, mass, adduct)
	else peaks <- .eieio_sg_peak(.eieio_sg_ion_mz(mass, adduct), 999, "Precursor", "precursor")
	chains <- .eieio_sg_parse_chains(lipid)
	if (nrow(chains)) {
		for (i in seq_len(nrow(chains))) {
			chain <- chains[i, , drop = FALSE]
						if (class %in% c("PC_d5", "PE_d5", "PG_d5", "PI_d5", "PS_d5", "LPC_d5", "LPE_d5", "LPG_d5", "LPI_d5", "LPS_d5", "DG_d5", "TG_d5", "CE_d7")) peaks <- .eieio_sg_isotope_chain_peaks(peaks, mass, chain, adduct, class, if (is.null(abundance)) 50 else abundance, i)
						else if (class %in% c("SM_d9", "Cer_NS_d7")) peaks <- .eieio_sg_isotope_sphingolipid_chain_peaks(peaks, mass, chain, adduct, class)
						else if (base_class %in% c("Cer_NP", "Cer_AP", "Cer_ABP")) peaks <- .eieio_sg_phyto_chain_peaks(peaks, chain, adduct)
						else if (base_class == "PE") peaks <- .eieio_sg_pe_chain_peaks(peaks, mass, chain, adduct)
						else if (base_class == "PG") peaks <- .eieio_sg_pg_chain_peaks(peaks, mass, chain, adduct)
						else if (base_class %in% c("LPC", "LPE", "LPG", "LPI")) peaks <- .eieio_sg_lysophospholipid_chain_peaks(peaks, mass, chain, adduct, base_class)
						else if (base_class == "LPS") peaks <- .eieio_sg_lps_chain_level(peaks, mass, chain, adduct)
			else peaks <- .eieio_sg_chain_peaks(peaks, mass, chain, adduct, if (is.null(abundance)) 50 else abundance)
			if (base_class %in% c("LPC", "LPE", "LPG", "LPI") && i == 1L) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_lysophospholipid_position(mass, chain, adduct, base_class))
			if (base_class == "PG" && i == 1L) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_rbind(
				.eieio_sg_peak(.eieio_sg_ion_mz(mass - .eieio_sg_chain_mass(chain) + 2 * .eieio_sg_h - .eieio_sg_h2o - .eieio_sg_ch2, adduct), 100, "-CH2(Sn1)", "snposition"),
				.eieio_sg_rbind(
					.eieio_sg_peak(.eieio_sg_ion_mz(mass - .eieio_sg_chain_mass(chain) + 2 * .eieio_sg_h - 2 * .eieio_sg_h2o - .eieio_sg_ch2, adduct), 100, "-H2O -CH2(Sn1)", "snposition"),
					.eieio_sg_peak(.eieio_sg_ion_mz(mass - .eieio_sg_chain_mass(chain) + 2 * .eieio_sg_h - .eieio_sg_c3h9o6p - .eieio_sg_h2o - .eieio_sg_ch2, adduct), 100, "-Header -CH2(Sn1)", "snposition"))))
			type <- if ("type" %in% names(chain)) tolower(chain$type[1L]) else "acyl"
						db_abundance <- if (base_class == "LPS") 30 else if (is.null(abundance)) 10 else abundance
						nl_mass <- 0
						if (class == "PG_d5") nl_mass <- if (.eieio_sg_adduct_name(adduct) == "[M+Na]+") 0 else .eieio_sg_ion_mz(.eieio_sg_c3h9o6p, adduct) - .eieio_sg_proton
						if (class == "PI_d5") nl_mass <- if (.eieio_sg_adduct_name(adduct) == "[M+Na]+") 0 else .eieio_sg_ion_mz(.eieio_sg_c6h13o9p, adduct) - .eieio_sg_proton
						if (class == "PS_d5") nl_mass <- if (.eieio_sg_adduct_name(adduct) == "[M+Na]+") 0 else .eieio_sg_c3h8no6p
						if (class == "DG_d5") nl_mass <- if (.eieio_sg_adduct_name(adduct) == "[M+Na]+") 0 else .eieio_sg_adduct_mass(adduct) + .eieio_sg_h2o - .eieio_sg_proton
						if (class == "CE_d7") nl_mass <- NA_real_
						if (class %in% c("SM_d9", "Cer_NS_d7")) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_sphingolipid_double_bond(mass, chain, adduct, class, db_abundance))
						else if (base_class %in% c("Cer_NP", "Cer_AP", "Cer_ABP")) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_phyto_double_bond(mass, chain, adduct, db_abundance))
												else if (base_class == "PE") {
							peaks <- .eieio_sg_rbind(peaks, .eieio_sg_double_bond(mass, chain, adduct, nl_mass = 0, abundance = db_abundance, sphingo = type == "sphingo"))
							peaks <- .eieio_sg_rbind(peaks, .eieio_sg_double_bond(mass, chain, adduct, nl_mass = .eieio_sg_c2h8no4p, abundance = db_abundance, sphingo = type == "sphingo"))
						}
						else if (!is.na(nl_mass)) peaks <- .eieio_sg_rbind(peaks, .eieio_sg_double_bond(mass, chain, adduct, nl_mass = nl_mass, abundance = db_abundance, sphingo = type == "sphingo"))
		}
	}
	.eieio_sg_merge(peaks)
}

EIEIO_generate_spectrum <- generate_eieio_spectrum
get_eieio_spectrum <- generate_eieio_spectrum

for (.eieio_class in .eieio_sg_registered_classes) {
	local({
		class_name <- .eieio_class
		generator <- function(lipid, adduct = "[M+H]+", abundance = NULL) {
			if (is.list(lipid) || is.data.frame(lipid)) lipid$LipidClass <- class_name
			generate_eieio_spectrum(lipid, adduct, abundance)
		}
		assign(paste0("generate_", class_name, "_spectrum"), generator, envir = .GlobalEnv)
	})
}
