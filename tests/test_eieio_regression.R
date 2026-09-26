if (requireNamespace("OptiLCMS", quietly = TRUE) &&
    exists("get_eieio_based_lipidomics_matched_peaks_scores", envir = asNamespace("OptiLCMS"))) {
  ns <- asNamespace("OptiLCMS")
  get_eieio_based_lipidomics_matched_peaks_scores <- get("get_eieio_based_lipidomics_matched_peaks_scores", envir = ns)
  .eieio_dgta_fragment_masses <- get(".eieio_dgta_fragment_masses", envir = ns)
  new_ms_ref_search_parameter <- get("new_ms_ref_search_parameter", envir = ns)
  calculate_score <- get("calculate_score", envir = ns)
} else {
  find_r_dir <- function() {
    for (cand in c(".", "..", "../..")) {
      p <- file.path(cand, "R")
      if (dir.exists(p) && file.exists(file.path(p, "GetEieioBasedLipidomicsMatchedPeaksScores.R"))) return(p)
    }
    stop("Cannot find R source directory")
  }
  r_dir <- find_r_dir()
  source(file.path(r_dir, "GetEidBasedLipidomicsMatchedPeaksScores.R"))
  source(file.path(r_dir, "GetEieioBasedLipidomicsMatchedPeaksScores.R"))
  source(file.path(r_dir, "LipidsScoring.R"))
}

make_peak <- function(mass, intensity, comment, spectrum_comment,
                      resolution = intensity, charge = intensity) {
  data.frame(
    Mass = mass,
    Intensity = intensity,
    Resolution = resolution,
    Charge = charge,
    Comment = comment,
    SpectrumComment = spectrum_comment,
    AbsolutelyRequired = FALSE,
    stringsAsFactors = FALSE
  )
}

reference <- list(
  Name = "PC 36:2",
  AdductType = list(AdductIonName = "[M+H]+", AdductIonAccurateMass = 1.007276),
  Spectrum = do.call(rbind, list(
    make_peak(100, 100, "", "metaboliteclass"),
    make_peak(110, 100, "", "metaboliteclass"),
    make_peak(120, 100, "", "acylchain"),
    make_peak(130, 100, "", "acylchain"),
    make_peak(140, 100, "", "snposition"),
    make_peak(150, 1000, "", "doublebond,doublebond_high"),
    make_peak(160, 100, "", "doublebond,doublebond_low")
  ))
)

scan <- reference
scan$Spectrum$Charge <- c(100, 100, 100, 100, 100, 1000, 100)

stopifnot(identical(
  get_eieio_based_lipidomics_matched_peaks_scores(
    list(Spectrum = data.frame()), reference, 0.01, 0, 2000
  ),
  c(-1, -1)
))

pc_score <- get_eieio_based_lipidomics_matched_peaks_scores(
  scan, reference, 0.01, 0, 2000
)
stopifnot(length(pc_score) == 2L, pc_score[1] > 0, pc_score[2] == 7)

pe_reference <- reference
pe_reference$Name <- "PE 36:2"
pe_scan <- scan
pe_score <- get_eieio_based_lipidomics_matched_peaks_scores(
  pe_scan, pe_reference, 0.01, 0, 2000
)
stopifnot(length(pe_score) == 2L, pe_score[2] == 6)

dgta_reference <- reference
dgta_reference$Name <- "DGTA 36:2"
dgta_masses <- .eieio_dgta_fragment_masses(dgta_reference)
dgta_scan <- scan
dgta_scan$Spectrum <- rbind(
  dgta_scan$Spectrum,
  make_peak(dgta_masses[["dgts"]], 200, "", "")
)
dgta_score <- get_eieio_based_lipidomics_matched_peaks_scores(
  dgta_scan, dgta_reference, 0.01, 0, 2000
)
stopifnot(length(dgta_score) == 2L, dgta_score[2] == 7)

as_peak_list <- function(peaks) {
  lapply(seq_len(nrow(peaks)), function(index) as.list(peaks[index, , drop = FALSE]))
}

end_to_end_scan <- list(
  Spectrum = as_peak_list(scan$Spectrum),
  PrecursorMz = 800,
  IonMode = "Positive"
)
end_to_end_reference <- reference
end_to_end_reference$Spectrum <- as_peak_list(reference$Spectrum)
end_to_end_reference$ScanID <- 1L
end_to_end_reference$PrecursorMz <- 800
end_to_end_reference$ChromXs_RT <- 1
end_to_end_reference$CollisionCrossSection <- 1
end_to_end_reference$InChIKey <- "X"
end_to_end_parameter <- new_ms_ref_search_parameter()
end_to_end_parameter$CollisionType <- "EIEIO"
end_to_end_result <- calculate_score(
  list(PrecursorMz = 800, ChromXs_RT = 1, CollisionCrossSection = 1),
  end_to_end_scan,
  list(),
  end_to_end_reference,
  list(),
  end_to_end_parameter
)
stopifnot(
  end_to_end_result$MatchedPeaksCount == 7,
  end_to_end_result$MatchedPeaksPercentage > 0,
  is.finite(end_to_end_result$TotalScore)
)

cat("EIEIO regression checks passed\n")