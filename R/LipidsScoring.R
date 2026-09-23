# CalculateScore.R
# Translation of MsReferenceScorer.CalculateScore and all sub-functions from C#.
# Only the lipidomics CID/HCD path (line 59: GetLipidomicsMatchedPeaksScores) is implemented.
# Functions in GetLipidomicsMatchedPeaksScores.R are NOT repeated here.
#
# Data conventions:
#   spectrum  : list of lists, each with $Mass and $Intensity (sorted by Mass)
#   scan/prop : named list with $Spectrum (as above) and $PrecursorMz
#   property  : named list with $PrecursorMz, $ChromXs_RT (RT value), $CollisionCrossSection
#   reference : named list with $Spectrum, $PrecursorMz, $Name, $ScanID, $InChIKey,
#               $CompoundClass, $Comment, $ChromXs_RT, $CollisionCrossSection
#   isotopes  : list of lists, each with $RelativeAbundance
#   parameter : named list — see new_ms_ref_search_parameter() for all fields

# ─────────────────────────────────────────────────────────────────────────────
# 1. fix_mass_tolerance
#    MolecularFormulaUtility.FixMassTolerance
# ─────────────────────────────────────────────────────────────────────────────
fix_mass_tolerance <- function(tolerance, mass) {
  if (mass <= 500) return(tolerance)
  # ppm = abs(PpmCalculator(500, 500 + tolerance)) = tolerance / 500 * 1e6
  ppm <- abs(tolerance / 500 * 1e6)
  # ConvertPpmToMassAccuracy(mass, ppm) = ppm * mass / 1e6
  return(ppm * mass / 1e6)
}

# ─────────────────────────────────────────────────────────────────────────────
# 2. get_gaussian_similarity
#    MsScanMatching.GetGaussianSimilarity(double actual, double reference, double tolerance)
# ─────────────────────────────────────────────────────────────────────────────
get_gaussian_similarity <- function(actual, reference, tolerance) {
  exp(-0.5 * ((actual - reference) / tolerance)^2)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. get_isotope_ratio_similarity
#    MsScanMatching.GetIsotopeRatioSimilarity
# ─────────────────────────────────────────────────────────────────────────────
get_isotope_ratio_similarity <- function(peaks1, peaks2, targeted_mz, tolerance) {
  # peaks1, peaks2: list of lists, each with $RelativeAbundance
  is_available1 <- !is.null(peaks1) && length(peaks1) > 0
  is_available2 <- !is.null(peaks2) && length(peaks2) > 0
  if (!is_available1 || !is_available2) return(-1)
  if (peaks1[[1]]$RelativeAbundance <= 0 || peaks2[[1]]$RelativeAbundance <= 0) return(-1)

  minimum <- min(length(peaks1), length(peaks2))
  if (minimum < 2) return(1)
  similarity <- 0
  for (i in 2:minimum) {
    ratio1 <- peaks1[[i]]$RelativeAbundance / peaks1[[1]]$RelativeAbundance
    ratio2 <- peaks2[[i]]$RelativeAbundance / peaks2[[1]]$RelativeAbundance

    if (ratio1 <= 1 && ratio2 <= 1) {
      similarity <- similarity + abs(ratio1 - ratio2)
    } else {
      if (ratio1 > ratio2) {
        similarity <- similarity + (1 - ratio2 / ratio1)
      } else if (ratio2 > ratio1) {
        similarity <- similarity + (1 - ratio1 / ratio2)
      }
    }
  }
  return(1 - similarity)
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. get_binned_spectrum
#    SpectrumHandler.GetBinnedSpectrum(spectrum, bin)   [2-parameter overload]
#    Used internally by get_spectral_entropy_similarity.
#    Returns list of lists with $Mass and $Intensity.
# ─────────────────────────────────────────────────────────────────────────────
get_binned_spectrum <- function(spectrum, bin) {
  if (length(spectrum) == 0) return(list())

  # group peaks by floor(mass / bin)
  bins <- list()
  for (pk in spectrum) {
    massframe <- as.character(floor(pk$Mass / bin))
    if (is.null(bins[[massframe]])) {
      bins[[massframe]] <- list()
    }
    bins[[massframe]] <- c(bins[[massframe]], list(pk))
  }

  result <- list()
  for (group in bins) {
    intensities <- sapply(group, function(p) p$Intensity)
    masses      <- sapply(group, function(p) p$Mass)
    max_mass    <- masses[which.max(intensities)]
    sum_int     <- sum(intensities)
    result <- c(result, list(list(Mass = max_mass, Intensity = sum_int)))
  }
  return(result)
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. get_normalized_by_total_intensity_peaks
#    SpectrumHandler.GetNormalizedByTotalIntensityPeaks
# ─────────────────────────────────────────────────────────────────────────────
get_normalized_by_total_intensity_peaks <- function(spectrum) {
  if (length(spectrum) == 0) return(list())
  sum_int <- sum(sapply(spectrum, function(p) p$Intensity))
  if (sum_int == 0) return(spectrum)
  lapply(spectrum, function(p) list(Mass = p$Mass, Intensity = p$Intensity / sum_int))
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. get_combined_spectrum
#    SpectrumHandler.GetCombinedSpectrum
#    Both input spectra should already be normalized.
# ─────────────────────────────────────────────────────────────────────────────
get_combined_spectrum <- function(peaks1, peaks2, bin) {
  bins <- list()

  add_to_bins <- function(pk) {
    massframe <- as.character(floor(pk$Mass / bin))
    if (is.null(bins[[massframe]])) {
      bins[[massframe]] <<- list()
    }
    bins[[massframe]] <<- c(bins[[massframe]], list(pk))
  }

  for (pk in peaks1) add_to_bins(pk)
  for (pk in peaks2) add_to_bins(pk)

  result <- list()
  for (group in bins) {
    intensities <- sapply(group, function(p) p$Intensity)
    masses      <- sapply(group, function(p) p$Mass)
    max_mass    <- masses[which.max(intensities)]
    sum_int     <- sum(intensities) * 0.5
    result <- c(result, list(list(Mass = max_mass, Intensity = sum_int)))
  }
  return(result)
}

# ─────────────────────────────────────────────────────────────────────────────
# 7. get_spectral_entropy
#    MsScanMatching.GetSpectralEntropy
# ─────────────────────────────────────────────────────────────────────────────
get_spectral_entropy <- function(peaks) {
  if (length(peaks) == 0) return(0)
  intensities <- sapply(peaks, function(p) p$Intensity)
  sum_int <- sum(intensities)
  if (sum_int == 0) return(0)
  p <- intensities / sum_int
  p <- p[p > 0]
  -sum(p * log2(p))
}

# ─────────────────────────────────────────────────────────────────────────────
# 8. get_spectral_entropy_similarity
#    MsScanMatching.GetSpectralEntropySimilarity
# ─────────────────────────────────────────────────────────────────────────────
get_spectral_entropy_similarity <- function(peaks1, peaks2, bin) {
  is_available1 <- !is.null(peaks1) && length(peaks1) > 0
  is_available2 <- !is.null(peaks2) && length(peaks2) > 0
  if (!is_available1 || !is_available2) return(-1)

  norm1 <- get_normalized_by_total_intensity_peaks(peaks1)
  norm2 <- get_normalized_by_total_intensity_peaks(peaks2)
  combined <- get_combined_spectrum(norm1, norm2, bin)

  entropy12 <- get_spectral_entropy(combined)
  entropy1  <- get_spectral_entropy(get_binned_spectrum(peaks1, bin))
  entropy2  <- get_spectral_entropy(get_binned_spectrum(peaks2, bin))

  1 - (2 * entropy12 - entropy1 - entropy2) * 0.5
}

# Convert list-of-peaks spectrum to numeric vectors once.
.spectrum_to_vectors <- function(spectrum) {
  n <- length(spectrum)
  if (n == 0) {
    return(list(mass = numeric(0), intensity = numeric(0)))
  }
  mass <- numeric(n)
  intensity <- numeric(n)
  for (i in seq_len(n)) {
    pk <- spectrum[[i]]
    mass[i] <- pk$Mass
    intensity[i] <- pk$Intensity
  }
  list(mass = mass, intensity = intensity)
}

# Fast builder from numeric vectors; calls Rcpp when available.
.build_dot_product_buffer_from_vectors <- function(mass1, intensity1, mass2, intensity2,
                                                   bin, min_mz, max_mz, drive = "both") {
  n1 <- length(mass1)
  n2 <- length(mass2)
  if (n1 == 0 || n2 == 0) {
    return(list(mz = numeric(0), m = numeric(0), r = numeric(0),
                baseM = 0, baseR = 0, size = 0L))
  }

  if (exists("lipids_build_dot_buffer_cpp", mode = "function", inherits = TRUE)) {
    drive_reference <- identical(drive, "reference")
    cpp_builder <- get("lipids_build_dot_buffer_cpp", mode = "function", inherits = TRUE)
    return(cpp_builder(
      mass1, intensity1, mass2, intensity2,
      as.double(bin), as.double(min_mz), as.double(max_mz), drive_reference
    ))
  }

  focused_mz <- min_mz
  remain_m <- 1L
  remain_l <- 1L

  buf_mz <- numeric(n1 + n2)
  buf_m <- numeric(n1 + n2)
  buf_r <- numeric(n1 + n2)
  size <- 0L

  base_m <- -.Machine$double.xmax
  base_r <- -.Machine$double.xmax

  while (focused_mz <= max_mz) {
    sum_l <- 0
    for (i in seq(remain_l, n2)) {
      if (mass2[i] < focused_mz - bin) next
      if (mass2[i] < focused_mz + bin) {
        sum_l <- sum_l + intensity2[i]
      } else {
        remain_l <- i
        break
      }
      if (i == n2) remain_l <- n2 + 1L
    }

    sum_m <- 0
    for (i in seq(remain_m, n1)) {
      if (mass1[i] < focused_mz - bin) next
      if (mass1[i] < focused_mz + bin) {
        sum_m <- sum_m + intensity1[i]
      } else {
        remain_m <- i
        break
      }
      if (i == n1) remain_m <- n1 + 1L
    }

    size <- size + 1L
    buf_mz[size] <- focused_mz
    buf_m[size] <- sum_m
    buf_r[size] <- sum_l
    if (sum_m > base_m) base_m <- sum_m
    if (sum_l > base_r) base_r <- sum_l

    if (drive == "reference") {
      if (focused_mz + bin > mass2[n2]) break
      if (remain_l <= n2) {
        focused_mz <- mass2[remain_l]
      } else {
        break
      }
    } else {
      last_mass <- max(mass1[n1], mass2[n2])
      if (focused_mz + bin > last_mass) break

      next_m <- if (remain_m <= n1) mass1[remain_m] else Inf
      next_l <- if (remain_l <= n2) mass2[remain_l] else Inf

      if (next_l > focused_mz + bin && next_m <= focused_mz + bin) {
        focused_mz <- next_m
      } else if (next_l <= focused_mz + bin && next_m > focused_mz + bin) {
        focused_mz <- next_l
      } else {
        focused_mz <- min(next_m, next_l)
      }
    }
  }

  list(
    mz = buf_mz[seq_len(size)],
    m = buf_m[seq_len(size)],
    r = buf_r[seq_len(size)],
    baseM = base_m,
    baseR = base_r,
    size = size
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build sliding-window buffer for dot-product functions.
# Returns a list with $measured (data.frame: FocusedMz, Intensity) and
# $reference (same), $baseM, $baseR.
#
# This implements the shared pattern used by Weighted, Reverse, Enhanced, and
# Simple dot products: advance focusedMz through both spectra simultaneously.
# ─────────────────────────────────────────────────────────────────────────────
.build_dot_product_buffer <- function(peaks1, peaks2, bin, min_mz, max_mz, drive = "both") {
  vec1 <- .spectrum_to_vectors(peaks1)
  vec2 <- .spectrum_to_vectors(peaks2)
  .build_dot_product_buffer_from_vectors(
    vec1$mass, vec1$intensity,
    vec2$mass, vec2$intensity,
    bin, min_mz, max_mz, drive
  )
}

get_reverse_dot_product_from_vectors <- function(mass1, intensity1, mass2, intensity2,
                                                 bin, mass_begin, mass_end) {
  n1 <- length(mass1)
  n2 <- length(mass2)
  if (n1 == 0 || n2 == 0) return(-1)

  min_mz <- mass2[1]
  max_mz <- mass2[n2]
  if (mass_begin > min_mz) min_mz <- mass_begin
  if (max_mz > mass_end) max_mz <- mass_end

  buf <- .build_dot_product_buffer_from_vectors(
    mass1, intensity1, mass2, intensity2, bin, min_mz, max_mz, drive = "reference"
  )
  if (buf$baseM == 0 || buf$baseR == 0) return(0)

  norm_m <- buf$m / buf$baseM
  norm_r <- buf$r / buf$baseR
  l_counter <- sum(norm_r > 0.1)
  penalty <- .peak_count_penalty(l_counter)

  cutoff <- 0.01
  scalar_m <- 0
  scalar_r <- 0
  covariance <- 0

  for (i in seq_len(buf$size)) {
    if (norm_r[i] < cutoff) next
    scalar_m <- scalar_m + norm_m[i] * buf$mz[i]
    scalar_r <- scalar_r + norm_r[i] * buf$mz[i]
    covariance <- covariance + sqrt(norm_m[i] * norm_r[i]) * buf$mz[i]
  }

  if (scalar_m == 0 || scalar_r == 0) return(0)
  (covariance^2) / scalar_m / scalar_r * penalty
}

get_enhanced_dot_product_from_vectors <- function(mass1, intensity1, mass2, intensity2,
                                                  bin, mass_begin, mass_end, penalty = 0.6) {
  n1 <- length(mass1)
  n2 <- length(mass2)
  if (n1 == 0 || n2 == 0) return(-1)

  min_mz <- max(mass_begin, min(mass1[1], mass2[1]))
  max_mz <- min(mass_end, max(mass1[n1], mass2[n2]))

  buf <- .build_dot_product_buffer_from_vectors(
    mass1, intensity1, mass2, intensity2, bin, min_mz, max_mz, drive = "both"
  )
  if (buf$baseM == 0 || buf$baseR == 0) return(0)

  norm_m <- buf$m / buf$baseM
  norm_r <- buf$r / buf$baseR
  l_counter <- sum(norm_r > 0.1)
  pcp <- .peak_count_penalty(l_counter)

  cutoff <- 0.01
  scalar_m <- 0
  scalar_r <- 0
  covariance <- 0

  for (i in seq_len(buf$size)) {
    if (norm_r[i] < cutoff) next
    if (norm_m[i] == 0) {
      scalar_m <- scalar_m + norm_m[i] * (1 - penalty) * buf$mz[i]
    } else {
      scalar_m <- scalar_m + norm_m[i] * buf$mz[i]
    }
    scalar_r <- scalar_r + norm_r[i] * buf$mz[i]
    covariance <- covariance + sqrt(norm_m[i] * norm_r[i]) * buf$mz[i]
  }

  if (scalar_m == 0 || scalar_r == 0) return(0)
  (covariance^2) / scalar_m / scalar_r * pcp
}

get_weighted_dot_product_from_vectors <- function(mass1, intensity1, mass2, intensity2,
                                                  bin, mass_begin, mass_end) {
  n1 <- length(mass1)
  n2 <- length(mass2)
  if (n1 == 0 || n2 == 0) return(-1)

  min_mz <- min(mass1[1], mass2[1])
  max_mz <- max(mass1[n1], mass2[n2])
  if (mass_begin > min_mz) min_mz <- mass_begin
  if (max_mz > mass_end) max_mz <- mass_end

  buf <- .build_dot_product_buffer_from_vectors(
    mass1, intensity1, mass2, intensity2, bin, min_mz, max_mz, drive = "both"
  )
  if (buf$baseM == 0 || buf$baseR == 0) return(0)

  norm_m <- buf$m / buf$baseM
  norm_r <- buf$r / buf$baseR
  l_counter <- sum(norm_r > 0.1)
  penalty <- .peak_count_penalty(l_counter)

  cutoff <- 0.01
  scalar_m <- 0
  scalar_r <- 0
  covariance <- 0

  for (i in seq_len(buf$size)) {
    if (norm_m[i] < cutoff) next
    scalar_m <- scalar_m + norm_m[i] * buf$mz[i]
    scalar_r <- scalar_r + norm_r[i] * buf$mz[i]
    covariance <- covariance + sqrt(norm_m[i] * norm_r[i]) * buf$mz[i]
  }

  if (scalar_m == 0 || scalar_r == 0) return(0)
  (covariance^2) / scalar_m / scalar_r * penalty
}

get_simple_dot_product_from_vectors <- function(mass1, intensity1, mass2, intensity2,
                                                bin, mass_begin, mass_end) {
  n1 <- length(mass1)
  n2 <- length(mass2)
  if (n1 == 0 || n2 == 0) return(-1)

  max_mz <- max(mass1[n1], mass2[n2])
  if (max_mz > mass_end) max_mz <- mass_end

  remain_m <- 1L
  while (remain_m <= n1 && mass1[remain_m] < mass_begin - bin) remain_m <- remain_m + 1L
  remain_l <- 1L
  while (remain_l <= n2 && mass2[remain_l] < mass_begin - bin) remain_l <- remain_l + 1L

  next_m_mz <- if (remain_m <= n1) mass1[remain_m] else Inf
  next_l_mz <- if (remain_l <= n2) mass2[remain_l] else Inf
  focused_mz <- min(next_m_mz, next_l_mz)

  buf_m <- numeric(n1 + n2)
  buf_r <- numeric(n1 + n2)
  size <- 0L
  base_m <- -.Machine$double.xmax
  base_r <- -.Machine$double.xmax

  while (focused_mz <= max_mz) {
    sum_m <- 0
    i <- remain_m
    while (i <= n1) {
      if (mass1[i] < focused_mz - bin) { remain_m <- i + 1L; i <- i + 1L; next }
      if (mass1[i] < focused_mz + bin) { sum_m <- sum_m + intensity1[i]; remain_m <- i + 1L }
      else break
      i <- i + 1L
    }

    sum_r <- 0
    i <- remain_l
    while (i <= n2) {
      if (mass2[i] < focused_mz - bin) { remain_l <- i + 1L; i <- i + 1L; next }
      if (mass2[i] < focused_mz + bin) { sum_r <- sum_r + intensity2[i]; remain_l <- i + 1L }
      else break
      i <- i + 1L
    }

    size <- size + 1L
    buf_m[size] <- sum_m
    buf_r[size] <- sum_r
    if (sum_m > base_m) base_m <- sum_m
    if (sum_r > base_r) base_r <- sum_r

    if (focused_mz + bin > max(mass1[n1], mass2[n2])) break

    if (remain_m > n1 || remain_l > n2) {
      focused_mz <- if (remain_l > n2) mass1[remain_m] else mass2[remain_l]
      next
    }
    next_m <- mass1[remain_m]
    next_l <- mass2[remain_l]
    if (next_l > focused_mz + bin && next_m <= focused_mz + bin) {
      focused_mz <- next_m
    } else if (next_l <= focused_mz + bin && next_m > focused_mz + bin) {
      focused_mz <- next_l
    } else {
      focused_mz <- min(next_m, next_l)
    }
  }

  if (base_m == 0 || base_r == 0) return(0)
  norm_m <- buf_m[seq_len(size)] / base_m * 999
  norm_r <- buf_r[seq_len(size)] / base_r * 999

  scalar_m <- sum(norm_m)
  scalar_r <- sum(norm_r)
  covariance <- sum(sqrt(norm_m * norm_r))
  if (scalar_m == 0 || scalar_r == 0) return(0)
  (covariance^2) / scalar_m / scalar_r
}

# ─────────────────────────────────────────────────────────────────────────────
# Peak-count penalty (used by Weighted, Reverse, Enhanced)
# lSpectrumCounter = number of reference bins with normalized intensity > 0.1
# ─────────────────────────────────────────────────────────────────────────────
.peak_count_penalty <- function(l_spectrum_counter) {
  if      (l_spectrum_counter == 1) 0.75
  else if (l_spectrum_counter == 2) 0.88
  else if (l_spectrum_counter == 3) 0.94
  else if (l_spectrum_counter == 4) 0.97
  else                               1.0
}

# ─────────────────────────────────────────────────────────────────────────────
# 9. get_reverse_dot_product
#    MsScanMatching.GetReverseDotProduct
#    prop1 = measured scan, prop2 = reference scan
# ─────────────────────────────────────────────────────────────────────────────
get_reverse_dot_product <- function(prop1, prop2, bin, mass_begin, mass_end) {
  peaks1 <- prop1$Spectrum
  peaks2 <- prop2$Spectrum
  is_available1 <- !is.null(peaks1) && length(peaks1) > 0
  is_available2 <- !is.null(peaks2) && length(peaks2) > 0
  if (!is_available1 || !is_available2) return(-1)

  vec1 <- .spectrum_to_vectors(peaks1)
  vec2 <- .spectrum_to_vectors(peaks2)
  get_reverse_dot_product_from_vectors(
    vec1$mass, vec1$intensity,
    vec2$mass, vec2$intensity,
    bin, mass_begin, mass_end
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 10. get_enhanced_dot_product
#     MsScanMatching.GetEnhancedDotProduct
# ─────────────────────────────────────────────────────────────────────────────
get_enhanced_dot_product <- function(prop1, prop2, bin, mass_begin, mass_end, penalty = 0.6) {
  peaks1 <- prop1$Spectrum
  peaks2 <- prop2$Spectrum
  is_available1 <- !is.null(peaks1) && length(peaks1) > 0
  is_available2 <- !is.null(peaks2) && length(peaks2) > 0
  if (!is_available1 || !is_available2) return(-1)

  vec1 <- .spectrum_to_vectors(peaks1)
  vec2 <- .spectrum_to_vectors(peaks2)
  get_enhanced_dot_product_from_vectors(
    vec1$mass, vec1$intensity,
    vec2$mass, vec2$intensity,
    bin, mass_begin, mass_end, penalty
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 11. get_weighted_dot_product
#     MsScanMatching.GetWeightedDotProduct
# ─────────────────────────────────────────────────────────────────────────────
get_weighted_dot_product <- function(prop1, prop2, bin, mass_begin, mass_end) {
  peaks1 <- prop1$Spectrum
  peaks2 <- prop2$Spectrum
  is_available1 <- !is.null(peaks1) && length(peaks1) > 0
  is_available2 <- !is.null(peaks2) && length(peaks2) > 0
  if (!is_available1 || !is_available2) return(-1)

  vec1 <- .spectrum_to_vectors(peaks1)
  vec2 <- .spectrum_to_vectors(peaks2)
  get_weighted_dot_product_from_vectors(
    vec1$mass, vec1$intensity,
    vec2$mass, vec2$intensity,
    bin, mass_begin, mass_end
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 12. get_simple_dot_product
#     MsScanMatching.GetSimpleDotProduct
#     No mz-weighting; normalizes to 999; no peak-count penalty.
# ─────────────────────────────────────────────────────────────────────────────
get_simple_dot_product <- function(prop1, prop2, bin, mass_begin, mass_end) {
  peaks1 <- prop1$Spectrum
  peaks2 <- prop2$Spectrum
  is_available1 <- !is.null(peaks1) && length(peaks1) > 0
  is_available2 <- !is.null(peaks2) && length(peaks2) > 0
  if (!is_available1 || !is_available2) return(-1)

  vec1 <- .spectrum_to_vectors(peaks1)
  vec2 <- .spectrum_to_vectors(peaks2)
  get_simple_dot_product_from_vectors(
    vec1$mass, vec1$intensity,
    vec2$mass, vec2$intensity,
    bin, mass_begin, mass_end
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 13. get_refined_lipid_annotation_level
#     MsScanMatching.GetRefinedLipidAnnotationLevel
#     Returns a named list:
#       $name, $is_lipid_class_matched, $is_lipid_chain_matched,
#       $is_lipid_position_matched, $is_others
#     Depends on convert_msdial_lipidname_to_lipid_molecule_object_vs2 and
#     get_lipid_molecule_annotation_result from GetLipidomicsMatchedPeaksScores.R
# ─────────────────────────────────────────────────────────────────────────────
get_refined_lipid_annotation_level <- function(ms_scan_prop, mol_ms_ref, bin) {
  result_empty <- list(
    name                    = "",
    is_lipid_class_matched  = FALSE,
    is_lipid_chain_matched  = FALSE,
    is_lipid_position_matched = FALSE,
    is_others               = FALSE
  )

  is_available_scan <- !is.null(ms_scan_prop) &&
    !is.null(ms_scan_prop$Spectrum) && length(ms_scan_prop$Spectrum) > 0
  is_available_ref  <- !is.null(mol_ms_ref) &&
    !is.null(mol_ms_ref$Spectrum) && length(mol_ms_ref$Spectrum) > 0
  if (!is_available_scan || !is_available_ref) return(result_empty)

  comp_class <- if (!is.null(mol_ms_ref$CompoundClass)) mol_ms_ref$CompoundClass else ""
  comment    <- if (!is.null(mol_ms_ref$Comment))       mol_ms_ref$Comment       else ""

  if (comment != "SPLASH" && comp_class != "Unknown" && comp_class != "Others") {
    # special classes — return reference name as-is
    special_classes <- c("Cholesterol", "CholesterolSulfate", "Undefined", "BileAcid",
                         "Ac2PIM1", "Ac2PIM2", "Ac3PIM2", "Ac4PIM2", "LipidA")
    if (comp_class %in% special_classes) {
      return(list(
        name                    = mol_ms_ref$Name,
        is_lipid_class_matched  = FALSE,
        is_lipid_chain_matched  = FALSE,
        is_lipid_position_matched = FALSE,
        is_others               = TRUE
      ))
    }

    molecule <- tryCatch(
      convert_msdial_lipidname_to_lipid_molecule_object_vs2(mol_ms_ref),
      error = function(e) NULL
    )
    if (is.null(molecule) || is.null(molecule$Adduct)) {
      return(list(
        name                    = mol_ms_ref$Name,
        is_lipid_class_matched  = FALSE,
        is_lipid_chain_matched  = FALSE,
        is_lipid_position_matched = FALSE,
        is_others               = TRUE
      ))
    }

    annot_result <- tryCatch(
      get_lipid_molecule_annotation_result(ms_scan_prop, molecule, bin),
      error = function(e) NULL
    )
    if (is.list(annot_result) && length(annot_result) > 0L && !is.null(annot_result$AnnotationLevel)) {
      annotation_level <- suppressWarnings(as.integer(annot_result$AnnotationLevel[[1]]))
      if (is.na(annotation_level)) annotation_level <- 0L

      if (annotation_level == 1L) {
        return(list(
          name                    = annot_result$SublevelLipidName,
          is_lipid_class_matched  = TRUE,
          is_lipid_chain_matched  = FALSE,
          is_lipid_position_matched = FALSE,
          is_others               = FALSE
        ))
      } else if (annotation_level >= 2L) {
        refined_name <- if (annot_result$SublevelLipidName == annot_result$LipidName) {
          annot_result$SublevelLipidName
        } else {
          paste0(annot_result$SublevelLipidName, "|", annot_result$LipidName)
        }
        return(list(
          name                    = refined_name,
          is_lipid_class_matched  = TRUE,
          is_lipid_chain_matched  = TRUE,
          is_lipid_position_matched = FALSE,
          is_others               = FALSE
        ))
      } else {
        return(result_empty)
      }
    } else {
      return(result_empty)
    }
  } else {
    # SPLASH / Unknown / Others — return reference name
    return(list(
      name                    = mol_ms_ref$Name,
      is_lipid_class_matched  = FALSE,
      is_lipid_chain_matched  = FALSE,
      is_lipid_position_matched = FALSE,
      is_others               = TRUE
    ))
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# 14. new_ms_scan_match_result
#     Constructor for MsScanMatchResult (named list with default values)
# ─────────────────────────────────────────────────────────────────────────────
new_ms_scan_match_result <- function() {
  list(
    Name                        = "",
    LibraryID                   = 0L,
    InChIKey                    = "",
    Source                      = "",
    AnnotatorID                 = "",
    Priority                    = 0L,

    # squared scores (stored directly)
    SquaredWeightedDotProduct   = 0,
    SquaredSimpleDotProduct     = 0,
    SquaredReverseDotProduct    = 0,

    # derived display scores: sqrt(max(squared, 0))
    WeightedDotProduct          = 0,
    SimpleDotProduct            = 0,
    ReverseDotProduct           = 0,

    EnhancedDotProduct          = 0,
    SpectralEntropy             = 0,

    MatchedPeaksPercentage      = 0,
    MatchedPeaksCount           = 0,
    AcurateMassSimilarity       = 0,
    IsotopeSimilarity           = 0,
    RtSimilarity                = 0,
    CcsSimilarity               = 0,
    TotalScore                  = 0,

    # Boolean flags
    IsPrecursorMzMatch          = FALSE,
    IsSpectrumMatch             = FALSE,
    IsRtMatch                   = FALSE,
    IsCcsMatch                  = FALSE,
    IsLipidClassMatch           = FALSE,
    IsLipidChainsMatch          = FALSE,
    IsLipidPositionMatch        = FALSE,
    IsOtherLipidMatch           = FALSE,
    IsReferenceMatched          = FALSE,
    IsAnnotationSuggested       = FALSE
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 15. validate_base
#     MsReferenceScorer.ValidateBase (lipidomics path only)
# ─────────────────────────────────────────────────────────────────────────────
validate_base <- function(result, property, reference, parameter) {
  # IsSpectrumMatch (lipidomics)
  result$IsSpectrumMatch <-
    result$SquaredWeightedDotProduct >= parameter$SquaredWeightedDotProductCutOff ||
    result$SquaredSimpleDotProduct   >= parameter$SquaredSimpleDotProductCutOff   ||
    result$SquaredReverseDotProduct  >= parameter$SquaredReverseDotProductCutOff

  comp_class <- if (!is.null(reference$CompoundClass)) reference$CompoundClass else ""
  if ((comp_class == "EtherTG" || comp_class == "EtherDG") &&
      result$SquaredSimpleDotProduct < parameter$SquaredSimpleDotProductCutOff) {
    result$IsSpectrumMatch <- FALSE
  }

  ms1_tol <- fix_mass_tolerance(parameter$Ms1Tolerance, property$PrecursorMz)
  result$IsPrecursorMzMatch <- abs(property$PrecursorMz - reference$PrecursorMz) <= ms1_tol

  if (!is.null(parameter$IsUseTimeForAnnotationScoring) && parameter$IsUseTimeForAnnotationScoring) {
    result$IsRtMatch <- abs(property$ChromXs_RT - reference$ChromXs_RT) <= parameter$RtTolerance
  }

  if (!is.null(parameter$IsUseCcsForAnnotationScoring) && parameter$IsUseCcsForAnnotationScoring) {
    result$IsCcsMatch <- abs(property$CollisionCrossSection - reference$CollisionCrossSection) <= parameter$CcsTolerance
  }

  result
}

# ─────────────────────────────────────────────────────────────────────────────
# 16. validate_on_lipidomics
#     MsReferenceScorer.ValidateOnLipidomics (default / CID path only)
# ─────────────────────────────────────────────────────────────────────────────
validate_on_lipidomics <- function(result, scan, reference, parameter) {
  annot <- get_refined_lipid_annotation_level(scan, reference, parameter$Ms2Tolerance)

  result$IsLipidChainsMatch   <- annot$is_lipid_chain_matched
  result$IsLipidClassMatch    <- annot$is_lipid_class_matched
  result$IsLipidPositionMatch <- annot$is_lipid_position_matched
  result$IsOtherLipidMatch    <- annot$is_others

  result$IsSpectrumMatch <- result$IsSpectrumMatch &&
    (annot$is_lipid_chain_matched |
     annot$is_lipid_class_matched |
     annot$is_lipid_position_matched |
     annot$is_others)

  if (result$IsOtherLipidMatch) {
    result$Name <- if (nchar(annot$name) > 0) annot$name else reference$Name
    return(result)
  }

  if (!result$IsSpectrumMatch) {
    lipid <- convert_msdial_lipidname_to_lipid_molecule_object_vs2(reference)
    if (!is.null(lipid) &&
        !is.null(lipid$SublevelLipidName) && nchar(lipid$SublevelLipidName) > 0 &&
        lipid$LipidName != lipid$SublevelLipidName) {
      annot$name <- paste0(lipid$SublevelLipidName, "|", lipid$LipidName)
    }
  }

  result$Name <- if (nchar(annot$name) > 0) annot$name else reference$Name
  result
}

# ─────────────────────────────────────────────────────────────────────────────
# 17. validate
#     MsReferenceScorer.Validate (simplified to lipidomics default path)
# ─────────────────────────────────────────────────────────────────────────────
validate <- function(result, property, scan, reference, parameter, use_ms2 = TRUE) {
  result <- validate_base(result, property, reference, parameter)
  result <- validate_on_lipidomics(result, scan, reference, parameter)

  use_rt  <- !is.null(parameter$IsUseTimeForAnnotationScoring) && parameter$IsUseTimeForAnnotationScoring
  use_ccs <- !is.null(parameter$IsUseCcsForAnnotationScoring) && parameter$IsUseCcsForAnnotationScoring

  result$IsReferenceMatched <- result$IsPrecursorMzMatch &&
    (!use_rt  || result$IsRtMatch)  &&
    (!use_ccs || result$IsCcsMatch) &&
    (!use_ms2 || result$IsSpectrumMatch)

  result$IsAnnotationSuggested <- result$IsPrecursorMzMatch &&
    (!use_rt  || result$IsRtMatch)  &&
    (!use_ccs || result$IsCcsMatch) &&
    !result$IsReferenceMatched

  result
}

# ─────────────────────────────────────────────────────────────────────────────
# 18. calculate_score
#     MsReferenceScorer.CalculateScore
#     Only the lipidomics CID/HCD path (GetLipidomicsMatchedPeaksScores) is
#     implemented. get_lipidomics_matched_peaks_scores is defined in
#     GetLipidomicsMatchedPeaksScores.R and must be sourced before calling this.
#
#  Arguments:
#   property          : named list with $PrecursorMz, $ChromXs_RT, $CollisionCrossSection
#   scan              : named list with $Spectrum (list of lists $Mass/$Intensity),
#                       $PrecursorMz, $IonMode ("Positive" or "Negative")
#   scan_isotopes     : list of lists, each with $RelativeAbundance
#   reference         : named list with $Spectrum, $PrecursorMz, $Name, $ScanID,
#                       $InChIKey, $CompoundClass, $Comment, $ChromXs_RT,
#                       $CollisionCrossSection
#   reference_isotopes: list of lists, each with $RelativeAbundance
#   parameter         : named list (see new_ms_ref_search_parameter())
#   annotator_id      : character string, e.g. the library ID string
#   source_type       : character, e.g. "GeneratedLipid" or "MspDB"
#   priority          : integer
#   use_ms2           : logical (default TRUE)
# ─────────────────────────────────────────────────────────────────────────────
calculate_score <- function(property, scan, scan_isotopes, reference, reference_isotopes,
                            parameter, annotator_id = "", source_type = "GeneratedLipid",
                            priority = 0L, use_ms2 = TRUE) {

  ms2_tol    <- parameter$Ms2Tolerance
  mass_begin <- parameter$MassRangeBegin
  mass_end   <- parameter$MassRangeEnd

  scan_vec <- .spectrum_to_vectors(scan$Spectrum)
  ref_vec  <- .spectrum_to_vectors(reference$Spectrum)

  # ── dot-product scores ───────────────────────────────────────────────────
  sq_wdp <- get_weighted_dot_product_from_vectors(
    scan_vec$mass, scan_vec$intensity,
    ref_vec$mass, ref_vec$intensity,
    ms2_tol, mass_begin, mass_end
  )
  sq_sdp <- get_simple_dot_product_from_vectors(
    scan_vec$mass, scan_vec$intensity,
    ref_vec$mass, ref_vec$intensity,
    ms2_tol, mass_begin, mass_end
  )
  sq_rdp <- get_reverse_dot_product_from_vectors(
    scan_vec$mass, scan_vec$intensity,
    ref_vec$mass, ref_vec$intensity,
    ms2_tol, mass_begin, mass_end
  )
  sq_edp <- get_enhanced_dot_product_from_vectors(
    scan_vec$mass, scan_vec$intensity,
    ref_vec$mass, ref_vec$intensity,
    ms2_tol, mass_begin, mass_end, 0.6
  )

  spectrum_entropy <- get_spectral_entropy_similarity(
    scan$Spectrum, reference$Spectrum, ms2_tol)

  # ── lipidomics matched-peaks scores ────────────────────────────────────────
  # Route to EAD-specific or standard CID/HCD path based on source_type and collision type
  # 
  # For GeneratedLipid source:
  #   - EIEIO: get_eieio_based_lipidomics_matched_peaks_scores
  #   - EID: get_eid_based_lipidomics_matched_peaks_scores
  #   - OAD: get_oad_based_lipidomics_matched_peaks_scores
  #   - other: standard get_lipidomics_matched_peaks_scores
  # For MspDB or other sources: always standard path
  #
  # Function definitions:
  # - get_lipidomics_matched_peaks_scores in GetLipidomicsMatchedPeaksScores.R
  # - get_eieio_based_lipidomics_matched_peaks_scores in GetEieioBasedLipidomicsMatchedPeaksScores.R
  # - get_eid_based_lipidomics_matched_peaks_scores in GetEidBasedLipidomicsMatchedPeaksScores.R
  # - get_oad_based_lipidomics_matched_peaks_scores in GetOadBasedLipidomicsMatchedPeaksScores.R
  #
  matched_peaks_scores <- NULL
  if (source_type == "GeneratedLipid" && !is.null(parameter$CollisionType)) {
    collision_type <- parameter$CollisionType
    if (collision_type == "EIEIO") {
      matched_peaks_scores <- get_eieio_based_lipidomics_matched_peaks_scores(
        scan, reference, ms2_tol, mass_begin, mass_end)
    }
    else if (collision_type == "EID") {
      matched_peaks_scores <- get_eid_based_lipidomics_matched_peaks_scores(
        scan, reference, ms2_tol, mass_begin, mass_end)
    }
    else if (collision_type == "OAD") {
      matched_peaks_scores <- get_oad_based_lipidomics_matched_peaks_scores(
        scan, reference, ms2_tol, mass_begin, mass_end)
    }
    else {
      # Default to standard CID/HCD path (CID, HCD, or unrecognized)
      matched_peaks_scores <- get_lipidomics_matched_peaks_scores(
        scan, reference, ms2_tol, mass_begin, mass_end)
    }
  }
  else {
    # MspDB or other source, or no collision type specified: standard path
    matched_peaks_scores <- get_lipidomics_matched_peaks_scores(
      scan, reference, ms2_tol, mass_begin, mass_end)
  }
  # matched_peaks_scores[1] = MatchedPeaksPercentage
  # matched_peaks_scores[2] = MatchedPeaksCount

  # ── precursor m/z similarity ─────────────────────────────────────────────
  ms1_tol     <- fix_mass_tolerance(parameter$Ms1Tolerance, property$PrecursorMz)
  ms1_sim     <- get_gaussian_similarity(property$PrecursorMz, reference$PrecursorMz, ms1_tol)
  isotope_sim <- get_isotope_ratio_similarity(
    scan_isotopes, reference_isotopes, property$PrecursorMz, ms1_tol)

  # ── build result ─────────────────────────────────────────────────────────
  result <- new_ms_scan_match_result()
  result$Name                      <- reference$Name
  result$LibraryID                 <- reference$ScanID
  result$InChIKey                  <- if (!is.null(reference$InChIKey)) reference$InChIKey else ""
  result$Source                    <- source_type
  result$AnnotatorID               <- annotator_id
  result$Priority                  <- priority

  result$SquaredWeightedDotProduct <- as.double(sq_wdp)
  result$SquaredSimpleDotProduct   <- as.double(sq_sdp)
  result$SquaredReverseDotProduct  <- as.double(sq_rdp)

  result$WeightedDotProduct        <- sqrt(max(sq_wdp, 0))
  result$SimpleDotProduct          <- sqrt(max(sq_sdp, 0))
  result$ReverseDotProduct         <- sqrt(max(sq_rdp, 0))

  result$EnhancedDotProduct        <- sqrt(max(sq_edp, 0))
  result$SpectralEntropy           <- as.double(spectrum_entropy)

  result$MatchedPeaksPercentage    <- as.double(matched_peaks_scores[1])
  result$MatchedPeaksCount         <- as.double(matched_peaks_scores[2])
  result$AcurateMassSimilarity     <- as.double(ms1_sim)
  result$IsotopeSimilarity         <- as.double(isotope_sim)

  # ── optional RT / CCS similarity ─────────────────────────────────────────
  if (!is.null(parameter$IsUseTimeForAnnotationScoring) && parameter$IsUseTimeForAnnotationScoring) {
    rt_sim <- get_gaussian_similarity(
      property$ChromXs_RT, reference$ChromXs_RT, parameter$RtTolerance)
    result$RtSimilarity <- as.double(rt_sim)
  }

  if (!is.null(parameter$IsUseCcsForAnnotationScoring) && parameter$IsUseCcsForAnnotationScoring) {
    ccs_sim <- get_gaussian_similarity(
      property$CollisionCrossSection, reference$CollisionCrossSection, parameter$CcsTolerance)
    result$CcsSimilarity <- as.double(ccs_sim)
  }

  # ── scoring factors (lipidomics CID/HCD) ─────────────────────────────────
  dot_product_factor         <- 1.0
  reverse_dot_prod_factor    <- 2.0
  presence_percentage_factor <- 3.0
  msms_factor                <- 3.0
  rt_factor                  <- 0.5
  ccs_factor                 <- 0.5
  mass_factor                <- 1.0
  # Normalize matched-peaks component from [0, 2] to [0, 1].
  matched_peaks_norm <- min(max(result$MatchedPeaksPercentage, 0), 2) / 2

  # msms composite score (normalized to [0, 1])
  msms_score <- (
    result$WeightedDotProduct     * dot_product_factor         +
    result$SimpleDotProduct       * dot_product_factor         +
    result$ReverseDotProduct      * reverse_dot_prod_factor    +
    matched_peaks_norm            * presence_percentage_factor
  ) / (dot_product_factor + dot_product_factor + reverse_dot_prod_factor + presence_percentage_factor)

  # accumulate factor-weighted scores for TotalScore using weighted average
  weighted_sum <- 0
  total_weight <- 0

  if (result$AcurateMassSimilarity >= 0 && mass_factor > 0) {
    weighted_sum <- weighted_sum + result$AcurateMassSimilarity * mass_factor
    total_weight <- total_weight + mass_factor
  }
  if (result$WeightedDotProduct >= 0 && result$SimpleDotProduct >= 0 && result$ReverseDotProduct >= 0) {
    weighted_sum <- weighted_sum + msms_score * msms_factor
    total_weight <- total_weight + msms_factor
  }
  if (!is.null(parameter$IsUseTimeForAnnotationScoring) && parameter$IsUseTimeForAnnotationScoring &&
      result$RtSimilarity >= 0 && rt_factor > 0) {
    weighted_sum <- weighted_sum + result$RtSimilarity * rt_factor
    total_weight <- total_weight + rt_factor
  }
  if (!is.null(parameter$IsUseCcsForAnnotationScoring) && parameter$IsUseCcsForAnnotationScoring &&
      result$CcsSimilarity >= 0 && ccs_factor > 0) {
    weighted_sum <- weighted_sum + result$CcsSimilarity * ccs_factor
    total_weight <- total_weight + ccs_factor
  }
  # isotope_factor == 0, so isotopeSimilarity is not included

  result$TotalScore <- if (total_weight == 0) 0 else weighted_sum / total_weight

  # penalise if no InChIKey
  if (is.null(result$InChIKey) || nchar(result$InChIKey) == 0) {
    result$TotalScore <- result$TotalScore * 0.9
  }

  # Keep score strictly in [0, 1]
  result$TotalScore <- min(max(result$TotalScore, 0), 1)

  # ── validate ─────────────────────────────────────────────────────────────
  result <- validate(result, property, scan, reference, parameter, use_ms2)

  result
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: default search parameter constructor
# ─────────────────────────────────────────────────────────────────────────────
new_ms_ref_search_parameter <- function(
    ms1_tolerance                       = 0.01,
    ms2_tolerance                       = 0.025,
    mass_range_begin                    = 0,
    mass_range_end                      = 2000,
    rt_tolerance                        = 0.5,
    ccs_tolerance                       = 5.0,
    is_use_time_for_annotation_scoring  = FALSE,
    is_use_ccs_for_annotation_scoring   = FALSE,
    squared_weighted_dot_product_cutoff = 0.36,
    squared_simple_dot_product_cutoff   = 0.36,
    squared_reverse_dot_product_cutoff  = 0.36,
    matched_peaks_percentage_cutoff     = 0.0,
    minimum_spectrum_match              = 0
) {
  list(
    Ms1Tolerance                       = ms1_tolerance,
    Ms2Tolerance                       = ms2_tolerance,
    MassRangeBegin                     = mass_range_begin,
    MassRangeEnd                       = mass_range_end,
    RtTolerance                        = rt_tolerance,
    CcsTolerance                       = ccs_tolerance,
    IsUseTimeForAnnotationScoring      = is_use_time_for_annotation_scoring,
    IsUseCcsForAnnotationScoring       = is_use_ccs_for_annotation_scoring,
    SquaredWeightedDotProductCutOff    = squared_weighted_dot_product_cutoff,
    SquaredSimpleDotProductCutOff      = squared_simple_dot_product_cutoff,
    SquaredReverseDotProductCutOff     = squared_reverse_dot_product_cutoff,
    MatchedPeaksPercentageCutOff       = matched_peaks_percentage_cutoff,
    MinimumSpectrumMatch               = minimum_spectrum_match
  )
}
