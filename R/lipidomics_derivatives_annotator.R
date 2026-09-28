

#' Centralized MS2 Diagnostic Rule Dictionary for Lipid Classes
#' Uses Regular Expressions to capture families of lipids (e.g., "PC", "EtherPC", "OxPC")
#' 
#' @return A list of diagnostic rules
#' Comprehensive MS2 Diagnostic Rule Dictionary for 200+ Lipid Classes
#' Uses Regular Expressions to capture specific headgroups, fragments, and neutral losses.
#' 
#' @return A list of diagnostic rules (exact masses for high-resolution MS)
get_lipid_rules_old <- function() {
  list(
    # ================= 1. N-ACYL AMIDES & AMINO ACIDS =================
    # Characteristic cleavage of the amide bond yielding [Amine + H]+ or specific negative fragments
    list(pattern = "NAGABA", mode = "POS", type = "fragment", val = 104.0706),  # GABA + H
    list(pattern = "NAGlySer", mode = "POS", type = "fragment", val = 163.0713),# Gly-Ser + H
    list(pattern = "^NAGly", mode = "POS", type = "fragment", val = 76.0393),   # Glycine + H
    list(pattern = "NATau", mode = "NEG", type = "fragment", val = 124.0074),   # Taurine - H
    list(pattern = "NASer", mode = "POS", type = "fragment", val = 106.0499),   # Serine + H
    list(pattern = "NAAla", mode = "POS", type = "fragment", val = 90.0550),    # Alanine + H
    list(pattern = "NALeu", mode = "POS", type = "fragment", val = 132.1019),   # Leucine + H
    list(pattern = "NAVal", mode = "POS", type = "fragment", val = 118.0863),   # Valine + H
    list(pattern = "NAPhe", mode = "POS", type = "fragment", val = 166.0863),   # Phenylalanine + H
    list(pattern = "NAOrn", mode = "POS", type = "fragment", val = 133.0972),   # Ornithine + H
    list(pattern = "NA5HT", mode = "POS", type = "fragment", val = 177.0917),   # Serotonin + H
    list(pattern = "NATryA", mode = "POS", type = "fragment", val = 161.0968),  # Tryptamine + H
    list(pattern = "NAGln", mode = "POS", type = "fragment", val = 147.0764),   # Glutamine + H
    list(pattern = "NAAnt", mode = "POS", type = "fragment", val = 138.0550),   # Anthranilate + H
    list(pattern = "NAE|GPNAE", mode = "POS", type = "fragment", val = 62.0600),# Ethanolamine + H
    list(pattern = "CAR", mode = "POS", type = "fragment", val = 85.0289),      # Carnitine characteristic backbone
    
    # ================= 2. PHOSPHOLIPIDS & SPHINGOMYELINS =================
    list(pattern = "PC|SM|LPC|bmPC", mode = "POS", type = "fragment", val = 184.0739), # Phosphocholine
    list(pattern = "PE|LPE|MMPE|DMPE|LNAPE", mode = "POS", type = "nl", val = 141.0191),        # Phosphoethanolamine NL
    list(pattern = "PE|LPE|MMPE|DMPE|LNAPE", mode = "NEG", type = "fragment", val = 140.0110),  # Phosphoethanolamine frag
    list(pattern = "PS|LPS|LNAPS", mode = "POS", type = "nl", val = 185.0089),        # Phosphoserine NL
    list(pattern = "PS|LPS|LNAPS", mode = "NEG", type = "nl", val = 87.0320),         # Serine headgroup NL
    list(pattern = "PI|LPI|PIM", mode = "NEG", type = "fragment", val = 241.0115),  # Phosphoinositol
    list(pattern = "PT", mode = "POS", type = "nl", val = 199.0240),            # Phosphothreonine NL
    
    # Glycerophosphate backbone lipids (Cardiolipins, BMP, PG, PA)
    list(pattern = "PG|LPG|PA|LPA|CL|DLCL|MLCL|BMP|HBMP|PMeOH|PEtOH|PBtOH", mode = "NEG", type = "fragment", val = 152.9953), 
    list(pattern = "PG|LPG", mode = "POS", type = "nl", val = 189.0402),
    
    # Phosphatidyl Alcohols
    list(pattern = "PMeOH", mode = "POS", type = "nl", val = 112.9998),         # Phosphomethanol NL
    list(pattern = "PEtOH", mode = "POS", type = "nl", val = 127.0155),         # Phosphoethanol NL
    list(pattern = "PBtOH", mode = "POS", type = "nl", val = 155.0468),         # Phosphobutanol NL
    
    # Methylated PEs
    list(pattern = "DMPE", mode = "POS", type = "nl", val = 169.0504),
    list(pattern = "MMPE", mode = "POS", type = "nl", val = 155.0347),
    
    # ================= 3. GLYCOSPHINGOLIPIDS & GANGLIOSIDES =================
    # Gangliosides (Sialic acid / NeuAc fragments)
    list(pattern = "GM1|GM3|GD1|GD2|GD3|GT1|GQ1", mode = "NEG", type = "fragment", val = 290.0878), # NeuAc - H
    list(pattern = "NGcGM3", mode = "NEG", type = "fragment", val = 306.0827),  # N-glycolyl NeuAc - H
    
    # Ceramides with Hexose
    list(pattern = "HexCer|AHexCer|ASHexCer", mode = "POS", type = "nl", val = 162.0528),
    list(pattern = "Hex2Cer", mode = "POS", type = "nl", val = 324.1056),
    list(pattern = "Hex3Cer", mode = "POS", type = "nl", val = 486.1584),
    
    # ================= 4. GLYCEROLIPIDS =================
    list(pattern = "MGDG|EtherMGDG|MGMG", mode = "POS", type = "nl", val = 179.0555),
    list(pattern = "DGDG|EtherDGDG|DGMG", mode = "POS", type = "nl", val = 341.1084),
    list(pattern = "DGTS|DGTA|LDGTS|LDGTA", mode = "POS", type = "fragment", val = 236.1498), # Betaine headgroup
    list(pattern = "DGCC|LDGCC", mode = "POS", type = "fragment", val = 104.1075),            # Choline-like headgroup
    list(pattern = "DGGA|ADGGA", mode = "NEG", type = "nl", val = 176.0321),                  # Glucuronide NL
    list(pattern = "SQDG|SMGDG", mode = "NEG", type = "fragment", val = 225.0069),
    
    # ================= 5. SULFATED LIPIDS & BILE ACIDS =================
    # Universal sulfate group fragment [HSO4]-
    list(pattern = "Sulfate|SHex|SL|ASM|SQDG|SMGDG", mode = "NEG", type = "fragment", val = 96.9596),
    list(pattern = "SHexCer", mode = "NEG", type = "fragment", val = 241.0018),               # Galactose-sulfate fragment
    
    # ================= 6. STEROLS =================
    # Steryl esters typically lose the fatty acid, but CE consistently shows the cholesterol backbone
    list(pattern = "^CE|CE_d7|Cholesterol", mode = "POS", type = "fragment", val = 369.3516),             # Cholesterol - H2O + H
    list(pattern = "BRSE", mode = "POS", type = "fragment", val = 381.3516),                  # Brassicasterol backbone
    list(pattern = "CASE", mode = "POS", type = "fragment", val = 383.3672),                  # Campesterol backbone
    list(pattern = "SISE", mode = "POS", type = "fragment", val = 397.3829),                  # Sitosterol backbone
    list(pattern = "STSE", mode = "POS", type = "fragment", val = 395.3672),                  # Stigmasterol backbone
    # ================= 7. PRENOLS & VITAMINS =================
    # Coenzyme Q (Ubiquinone) characteristic tropylium ion fragment
    list(pattern = "CoQ", mode = "POS", type = "fragment", val = 197.0813),
    
    # Vitamin E (alpha-Tocopherol) characteristic chromanol ring fragment
    list(pattern = "Vitamin_E", mode = "POS", type = "fragment", val = 165.0915),
    
    # Vitamin D (D2/D3) often relies on side-chain cleavage + water loss
    # E.g., m/z 259.2056 is a highly conserved steroid core fragment for Vit D
    list(pattern = "Vitamin_D", mode = "POS", type = "fragment", val = 259.2056),
    
    # Vitamin A (Retinol esters) highly conserved conjugated triene fragment
    list(pattern = "VAE", mode = "POS", type = "fragment", val = 269.2269),
    
    # ================= 8. BILE ACIDS & ESTERS =================
    # Tauro-conjugated bile acids (TDCAE, TLCAE)
    list(pattern = "^T[A-Z]{3,4}E", mode = "NEG", type = "fragment", val = 124.0074), # Taurine fragment
    list(pattern = "^T[A-Z]{3,4}E", mode = "POS", type = "nl", val = 125.0147),       # Taurine NL
    
    # Glyco-conjugated bile acids (GDCAE, GLCAE)
    list(pattern = "^G[A-Z]{3,4}E", mode = "NEG", type = "fragment", val = 74.0249),  # Glycine fragment
    list(pattern = "^G[A-Z]{3,4}E", mode = "POS", type = "nl", val = 75.0320),        # Glycine NL
    
    # Unconjugated / Keto Bile Acid Esters (DCAE, LCAE, KLCAE, KDCAE)
    # These lack amino acid headgroups and mainly undergo sequential water losses.
    # While they can be assigned generic water loss NLs (-18.0105), it is mathematically 
    # much safer to let the Modified Cosine algorithm handle their structural confirmation.
    list(pattern = "KLCAE|KDCAE|DCAE|LCAE", mode = "POS", type = "nl", val = 18.0105),
    
    # ================= 9. COMPLEX STEROL CONJUGATES =================
    # Steryl phosphohexosides (BRSPHex, CASPHex, CSPHex, SISPHex, STSPHex)
    list(pattern = "PHex", mode = "POS", type = "nl", val = 242.0191),                # Phosphohexose NL
    
    # Acylated hexosyl sterols (AHexCS, AHexBRS, etc.) & Steryl hexosides
    list(pattern = "AHex|SHex$", mode = "POS", type = "nl", val = 162.0528),          # Hexose NL
    
    # Steryl phosphoethanolamines (SPE)
    list(pattern = "SPE", mode = "POS", type = "nl", val = 141.0191)                  # PE headgroup NL
  )
}


# ==============================================================================
# MODULE 1: Exhaustive Lipid Rules Dictionary (Supporting 150+ Classes)
# ==============================================================================
get_lipid_rules <- function() {
  list(
    # ================= 1. N-ACYL AMIDES & AMINO ACIDS =================
    list(pattern = "NAGABA", mode = "POS", type = "fragment", val = 104.0706),
    list(pattern = "NAGlySer", mode = "POS", type = "fragment", val = 163.0713),
    list(pattern = "^NAGly", mode = "POS", type = "fragment", val = 76.0393),
    list(pattern = "NATau", mode = "NEG", type = "fragment", val = 124.0074),
    list(pattern = "NASer", mode = "POS", type = "fragment", val = 106.0499),
    list(pattern = "NAAla", mode = "POS", type = "fragment", val = 90.0550),
    list(pattern = "NALeu", mode = "POS", type = "fragment", val = 132.1019),
    list(pattern = "NAVal", mode = "POS", type = "fragment", val = 118.0863),
    list(pattern = "NAPhe", mode = "POS", type = "fragment", val = 166.0863),
    list(pattern = "NAOrn", mode = "POS", type = "fragment", val = 133.0972),
    list(pattern = "NA5HT", mode = "POS", type = "fragment", val = 177.0917),
    list(pattern = "NATryA", mode = "POS", type = "fragment", val = 161.0968),
    list(pattern = "NAGln", mode = "POS", type = "fragment", val = 147.0764),
    list(pattern = "NAAnt", mode = "POS", type = "fragment", val = 138.0550),
    list(pattern = "NAE|GPNAE", mode = "POS", type = "fragment", val = 62.0600),
    list(pattern = "CAR", mode = "POS", type = "fragment", val = 85.0289),
    
    # ================= 2. PHOSPHOLIPIDS & SPHINGOMYELINS =================
    list(pattern = "PC|SM|LPC|bmPC|ASM", mode = "POS", type = "fragment", val = 184.0739), 
    list(pattern = "PE|LPE|MMPE|DMPE|LNAPE|PE-Cer", mode = "POS", type = "nl", val = 141.0191),
    list(pattern = "PE|LPE|MMPE|DMPE|LNAPE|PE-Cer", mode = "NEG", type = "fragment", val = 140.0110),
    list(pattern = "PS|LPS|LNAPS", mode = "POS", type = "nl", val = 185.0089),
    list(pattern = "PS|LPS|LNAPS", mode = "NEG", type = "nl", val = 87.0320),
    list(pattern = "PI|LPI|PIM|PI-Cer|MIPC", mode = "NEG", type = "fragment", val = 241.0115),
    list(pattern = "PT", mode = "POS", type = "nl", val = 199.0240),
    
    # Glycerophosphate backbone lipids
    list(pattern = "PG|LPG|PA|LPA|CL|DLCL|MLCL|BMP|HBMP|PMeOH|PEtOH|PBtOH", mode = "NEG", type = "fragment", val = 152.9953), 
    list(pattern = "PG|LPG", mode = "POS", type = "nl", val = 189.0402),
    list(pattern = "PMeOH", mode = "POS", type = "nl", val = 112.9998),
    list(pattern = "PEtOH", mode = "POS", type = "nl", val = 127.0155),
    list(pattern = "PBtOH", mode = "POS", type = "nl", val = 155.0468),
    list(pattern = "DMPE", mode = "POS", type = "nl", val = 169.0504),
    list(pattern = "MMPE", mode = "POS", type = "nl", val = 155.0347),
    list(pattern = "CerP", mode = "NEG", type = "fragment", val = 96.9698), # H2PO4-
    
    # ================= 3. GLYCOSPHINGOLIPIDS & GANGLIOSIDES =================
    list(pattern = "GM1|GM3|GD1|GD2|GD3|GT1|GQ1", mode = "NEG", type = "fragment", val = 290.0878),
    list(pattern = "NGcGM3", mode = "NEG", type = "fragment", val = 306.0827),
    list(pattern = "HexCer|AHexCer|ASHexCer", mode = "POS", type = "nl", val = 162.0528),
    list(pattern = "Hex2Cer", mode = "POS", type = "nl", val = 324.1056),
    list(pattern = "Hex3Cer", mode = "POS", type = "nl", val = 486.1584),
    
    # ================= 4. GLYCEROLIPIDS =================
    list(pattern = "MGDG|EtherMGDG|MGMG", mode = "POS", type = "nl", val = 179.0555),
    list(pattern = "DGDG|EtherDGDG|DGMG", mode = "POS", type = "nl", val = 341.1084),
    list(pattern = "DGTS|DGTA|LDGTS|LDGTA", mode = "POS", type = "fragment", val = 236.1498),
    list(pattern = "DGCC|LDGCC", mode = "POS", type = "fragment", val = 104.1075),
    list(pattern = "DGGA|ADGGA", mode = "NEG", type = "nl", val = 176.0321),
    
    # ================= 5. SULFATED LIPIDS & BILE ACIDS =================
    list(pattern = "Sulfate|SHex|SL|ASM|SQDG|SMGDG|EtherSDGDG|EtherSMGDG", mode = "NEG", type = "fragment", val = 96.9596),
    list(pattern = "SHexCer", mode = "NEG", type = "fragment", val = 241.0018),
    
    # ================= 6. STEROLS =================
    list(pattern = "^CE|CE_d7|Cholesterol", mode = "POS", type = "fragment", val = 369.3516),
    list(pattern = "BRSE", mode = "POS", type = "fragment", val = 381.3521),
    list(pattern = "CASE", mode = "POS", type = "fragment", val = 383.3678),
    list(pattern = "SISE", mode = "POS", type = "fragment", val = 397.3834),
    list(pattern = "STSE", mode = "POS", type = "fragment", val = 395.3678),
    list(pattern = "EGSE", mode = "POS", type = "fragment", val = 379.3365), # Ergosteryl
    list(pattern = "DEGSE", mode = "POS", type = "fragment", val = 377.3208),# Dehydroergosteryl
    list(pattern = "DSMSE", mode = "POS", type = "fragment", val = 367.3365),# Desmosteryl
    
    # ================= 7. PRENOLS & VITAMINS =================
    list(pattern = "CoQ", mode = "POS", type = "fragment", val = 197.0813),
    list(pattern = "Vitamin_E", mode = "POS", type = "fragment", val = 165.0915),
    list(pattern = "Vitamin_D", mode = "POS", type = "fragment", val = 259.2056),
    list(pattern = "VAE", mode = "POS", type = "fragment", val = 269.2269),
    
    # ================= 8. BILE ACIDS & ESTERS =================
    list(pattern = "^T[A-Z]{3,4}E", mode = "NEG", type = "fragment", val = 124.0074),
    list(pattern = "^T[A-Z]{3,4}E", mode = "POS", type = "nl", val = 125.0147),
    list(pattern = "^G[A-Z]{3,4}E", mode = "NEG", type = "fragment", val = 74.0249),
    list(pattern = "^G[A-Z]{3,4}E", mode = "POS", type = "nl", val = 75.0320),
    list(pattern = "KLCAE|KDCAE|DCAE|LCAE", mode = "POS", type = "nl", val = 18.0105),
    
    # ================= 9. COMPLEX STEROL CONJUGATES =================
    list(pattern = "PHex", mode = "POS", type = "nl", val = 242.0191),
    list(pattern = "AHex|SHex$|BAHex", mode = "POS", type = "nl", val = 162.0528),
    list(pattern = "SPE", mode = "POS", type = "nl", val = 141.0191)
  )
}

# ==============================================================================
# MODULE 2: Universal Class-Specific Backbone Generator (For In-Silico Matching)
# ==============================================================================
#' Generate lipid class-specific backbone fragments dynamically using Regex
#' 
#' @param precursor_mz Numeric, accurate derivatized precursor m/z
#' @param lipid_class Character, exact lipid class abbreviation
#' @param ion_mode Character, "POS" or "NEG"
#' @param chains_dt data.table containing chain info
#' @param method Character, derivatization method to adjust FA loss mass
#' @return list of data.tables containing theoretical backbone peaks
generate_class_backbone_peaks <- function(precursor_mz, lipid_class, ion_mode, chains_dt, method) {
  peaks <- list()
  
  # Basic masses for calculation
  c_mass <- 12.00000000000
  h_mass <- 1.00782503223
  o_mass <- 15.99491461957
  
  add_peak <- function(mz, int, comment) {
    peaks[[length(peaks) + 1]] <<- data.table(mz = mz, theo_int = int, comment = comment)
  }
  
  # 1. Choline headgroup (PC, LPC, SM, etc.)
  if (grepl("PC|LPC|SM|bmPC|ASM", lipid_class, ignore.case = TRUE) && ion_mode == "POS") {
    add_peak(184.0739, 100, "Phosphocholine Headgroup")
  } 
  
  # 2. Ethanolamine headgroup (PE, LPE, etc.)
  else if (grepl("PE|LPE|MMPE|DMPE|LNAPE|PE-Cer", lipid_class, ignore.case = TRUE)) {
    if (ion_mode == "POS") {
      add_peak(precursor_mz - 141.0191, 100, "PE Neutral Loss (-C2H8NO4P)")
    } else if (ion_mode == "NEG") {
      add_peak(140.0110, 100, "Phosphoethanolamine Headgroup")
    }
  } 
  
  # 3. Serine headgroup (PS, LPS)
  else if (grepl("PS|LPS|LNAPS", lipid_class, ignore.case = TRUE)) {
    if (ion_mode == "POS") {
      add_peak(precursor_mz - 185.0089, 100, "PS Neutral Loss")
    } else if (ion_mode == "NEG") {
      add_peak(precursor_mz - 87.0320, 100, "Serine Headgroup NL")
    }
  }
  
  # 4. Inositol & Glycerophosphate (PI, PG, PA, CL, BMP)
  else if (grepl("PI|LPI|PIM|PI-Cer|MIPC", lipid_class, ignore.case = TRUE) && ion_mode == "NEG") {
    add_peak(241.0115, 100, "Phosphoinositol Fragment")
  }
  else if (grepl("PG|LPG|PA|LPA|CL|MLCL|DLCL|BMP|HBMP", lipid_class, ignore.case = TRUE) && ion_mode == "NEG") {
    add_peak(152.9953, 100, "Glycerophosphate Fragment")
  }
  
  # 5. Glycolipids & Gangliosides
  else if (grepl("HexCer|Hex2Cer|Hex3Cer|AHex|SHex|PHex", lipid_class, ignore.case = TRUE) && ion_mode == "POS") {
    add_peak(precursor_mz - 162.0528, 80, "Hexose Neutral Loss")
  }
  else if (grepl("GM|GD|GT|GQ", lipid_class, ignore.case = TRUE) && ion_mode == "NEG") {
    add_peak(290.0878, 100, "NeuAc Fragment")
  }
  
  # 6. Neutral Glycerolipids (TG, DG, MG - excluding glycolipids)
  else if (grepl("^TG|^DG|^MG|OxTG", lipid_class, ignore.case = TRUE) && !grepl("MGDG|DGDG|SQDG|DGTS|DGTA|DGCC", lipid_class, ignore.case = TRUE)) {
    if (ion_mode == "POS") {
      add_peak(precursor_mz - 18.0338, 20, "[M+H]+ (Ammonia loss)")
      for (i in 1:nrow(chains_dt)) {
        if (chains_dt$carbons[i] == 0) next
        fa_nl_mass <- (chains_dt$carbons[i] * c_mass) + 
          ((2 * chains_dt$carbons[i] - 2 * chains_dt$db_count[i]) * h_mass) + 
          (2 * o_mass)
        
        if (chains_dt$db_count[i] > 0) {
          reagent_mass <- switch(method, "PB" = 58.04186, "DMDS" = 93.9884, "mCPBA" = 15.9949)
          fa_nl_mass <- fa_nl_mass + (chains_dt$db_count[i] * reagent_mass)
        }
        add_peak(precursor_mz - 18.0338 - fa_nl_mass, 90, sprintf("NL of FA %s", chains_dt$chain_id[i]))
      }
    }
  } 
  
  # 7. Ceramides & Basic Sphingolipids
  else if (grepl("^Cer|Sph|DHSph|PhytoSph", lipid_class, ignore.case = TRUE) && !grepl("Hex|SM|CerP", lipid_class, ignore.case = TRUE)) {
    if (ion_mode == "POS") {
      add_peak(precursor_mz - 18.0105, 100, "Sphingoid Water Loss (-H2O)")
      add_peak(precursor_mz - 36.0210, 60, "Sphingoid Double Water Loss (-2H2O)")
    }
  }
  
  # 8. Sterols 
  else if (grepl("^CE|^Cholesterol", lipid_class, ignore.case = TRUE) && ion_mode == "POS") {
    add_peak(369.3516, 100, "Cholesterol Backbone Fragment")
  }
  else if (grepl("EGSE", lipid_class, ignore.case = TRUE) && ion_mode == "POS") add_peak(379.3365, 100, "Ergosteryl Backbone")
  else if (grepl("DSMSE", lipid_class, ignore.case = TRUE) && ion_mode == "POS") add_peak(367.3365, 100, "Desmosteryl Backbone")
  
  # 9. Sulfates
  if (grepl("Sulfate|SHex|SL|SQDG|SMGDG", lipid_class, ignore.case = TRUE) && ion_mode == "NEG") {
    add_peak(96.9596, 100, "Sulfate Fragment")
  }
  
  # 10. Universal FA Carboxylate Anions (NEG Mode for oxygenated acyl chains)
  if (ion_mode == "NEG" && !grepl("^CE|Cholesterol|BRSE|CASE|SISE|STSE|EGSE|DSMSE|CoQ|Vitamin", lipid_class, ignore.case = TRUE)) {
    for (i in 1:nrow(chains_dt)) {
      if (chains_dt$carbons[i] == 0) next
      chain_mass <- (chains_dt$carbons[i] * c_mass) + 
        ((2 * chains_dt$carbons[i] - 1 - 2 * chains_dt$db_count[i]) * h_mass) + 
        (2 * o_mass)
      
      if (chains_dt$db_count[i] > 0) {
        reagent_mass <- switch(method, "PB" = 58.04186, "DMDS" = 93.9884, "mCPBA" = 15.9949)
        chain_mass <- chain_mass + (chains_dt$db_count[i] * reagent_mass)
      }
      add_peak(chain_mass, 80, sprintf("Carboxylate Anion %s", chains_dt$chain_id[i]))
    }
  }
  
  return(peaks)
}


#' Check for class-specific polar headgroup in MS2 spectrum
#'
#' @param lipid_annotation Character, full annotation string (e.g., "EtherPE 36:2")
#' @param ms2_spec data.table/data.frame of derivatized MS2 (columns: mz, intensity)
#' @param precursor_mz Numeric, precursor m/z (needed for neutral loss calculation)
#' @param ion_mode Character, "POS" or "NEG" (Positive or Negative ionization mode)
#' @param mz_tol Numeric, absolute m/z tolerance, default 0.02
#' @return Logical, TRUE if characteristic headgroup rule is met, FALSE otherwise
check_headgroup_presence <- function(lipid_annotation, ms2_spec, precursor_mz, 
                                     ion_mode = "POS", mz_tol = 0.02) {
  
  if (is.null(ms2_spec) || nrow(ms2_spec) == 0 || is.na(lipid_annotation) || lipid_annotation == "") {
    return(FALSE)
  }
  #save(lipid_annotation, ms2_spec, precursor_mz, ion_mode, mz_tol, file = "debugging.rda")
  # Extract the main lipid class (everything before the first space or underscore if formatted differently)
  # e.g., "EtherPE 36:2" -> "EtherPE"
  main_class <- strsplit(trimws(lipid_annotation), " ")[[1]][1]
  #cat("main class is ==> ", main_class, "\n")
  # Load rules dictionary
  rules <- get_lipid_rules()
  
  # Find applicable rules for this lipid class AND current ionization mode
  match_found <- has_headgroup <- FALSE
  
  for (rule in rules) {
    if (rule$mode == ion_mode && grepl(rule$pattern, main_class, ignore.case = TRUE)) {
      
      if (rule$type == "fragment") {
        # Check for direct characteristic fragment
        match_found <- any(abs(ms2_spec[[1]] - rule$val) <= mz_tol)
      } else if (rule$type == "nl") {
        # Check for characteristic neutral loss
        match_found <- any(abs(ms2_spec[[1]] - (precursor_mz - rule$val)) <= mz_tol)
      }
      
      if(is.na(match_found)){match_found <- FALSE}
      
      # If any applicable rule matches, mark as found and exit loop
      if (match_found) {
        has_headgroup <- TRUE
        break
      }
    }
  }
  
  return(has_headgroup)
}

#' Calculate Modified Cosine Score
calculate_modified_cosine <- function(spec1, spec2, delta_m, mz_tol = 0.02) {
  if(is.null(spec1)){
    cat("spec1 is empty null..")
  }
  if(nrow(spec1) == 0){
    cat("spec1 is empty..")
  }
  if(nrow(spec2) == 0){
    cat("spec2 is empty..")
  }
  if (is.null(spec1) || is.null(spec2) || nrow(spec1) == 0 || nrow(spec2) == 0) return(0)
  #cat("Performing Cosine Score..")
  # if(nrow(spec1) >1 & nrow(spec2) >1){
  #   cat("spec1 \n")
  #   cat(spec1)
  #   cat("\nspec2 \n")
  #   cat(spec2)
  #   save(spec1, spec2,delta_m, mz_tol, file = "debugging.rda")
  # }

  s1 <- data.table(mz = spec1[,1], int = sqrt(spec1[,2]))
  s2 <- data.table(mz = spec2[,1], int = sqrt(spec2[,2]))
  
  grid <- as.data.table(expand.grid(i = 1:nrow(s1), j = 1:nrow(s2)))
  grid[, mz1 := s1$mz[i]]
  grid[, mz2 := s2$mz[j]]
  
  grid[, direct_match := abs(mz1 - mz2) <= mz_tol]
  grid[, shifted_match := abs((mz1 + delta_m) - mz2) <= mz_tol]
  
  aligned <- grid[direct_match | shifted_match]
  if (nrow(aligned) == 0) return(0)
  
  aligned[, error := pmin(ifelse(direct_match, abs(mz1 - mz2), Inf), 
                          ifelse(shifted_match, abs((mz1 + delta_m) - mz2), Inf))]
  aligned <- aligned[order(error)]
  aligned <- aligned[!duplicated(i) & !duplicated(j)]
  
  aligned[, int1 := s1$int[i]]
  aligned[, int2 := s2$int[j]]
  dot_product <- sum(aligned$int1 * aligned$int2)
  
  norm1 <- sum(s1$int ^ 2)
  norm2 <- sum(s2$int ^ 2)
  if (norm1 == 0 || norm2 == 0) return(0)
  
  return(dot_product / sqrt(norm1 * norm2))
}

#' MS1 Pairing & Hybrid MS2 Scoring Function (Headgroup + Cosine)
#'
#' @param input_data Data containing MS1 features and a 'lipid_annotation' column
find_derivatization_pairs_hybrid_ms2 <- function(input_data, 
                                                 ms2_list = NULL,
                                                 db_hits = NULL,
                                                 ion_mode = "POS",
                                                 method = c("PB", "mCPBA", "DMDS"), 
                                                 intensity_cols = NULL, 
                                                 mz_tol_ppm = 10) {
  
  method <- match.arg(method)
  data <- as.data.table(copy(input_data))
  
  # formating ms2_list
  if(is.null(ms2_list)){
    stop("ms2_list is required for the following analysis!")
  } else {
    ms2_list0 <- ms2_list
    ms2_list <- lapply(0:(nrow(input_data)-1), function(x){
      idx <- which(x == ms2_list0[[1]])
      if(length(idx)==0){
        return(matrix(NA, nrow = 1, ncol = 2))
      } else {
        return(ms2_list0[[2]][[idx]][[1]])
      }
    })
  }
  
  # formatting db_hits to lipids_annotation
  if(is.null(db_hits)){
    stop("db_hits is required for the following analysis!")
  } else {
    lipid_annotation <- vapply(0:(nrow(input_data)-1), function(x){
      idx <- which(x == db_hits$FeatureIndex)
      if(length(idx)==0){
        return("")
      } else {
        #return(paste0(db_hits$AbbreName[idx], collapse = " "))
        return(db_hits$AbbreName[idx][1])
      }
    }, FUN.VALUE = character(length = 1L), USE.NAMES = F)
  }
  data <- cbind(data, lipid_annotation = lipid_annotation)
  
  # Ensure feature_id exists
  if (!"feature_id" %in% names(data)) data[, feature_id := as.character(1:.N)]
  # Ensure lipid_annotation column exists for headgroup check, default to empty if missing
  if (!"lipid_annotation" %in% names(data)) data[, lipid_annotation := ""]
  
  # Auto-generate feature_id if missing
  if (!"feature_id" %in% names(data)) {
    data[, feature_id := 1:.N]
  }
  
  # Handle NA values in intensity matrix to prevent correlation calculation errors
  if (!is.null(intensity_cols)) {
    for (col in intensity_cols) {
      set(data, which(is.na(data[[col]])), col, 0)
    }
    int_mat <- as.matrix(data[, ..intensity_cols])
    rownames(int_mat) <- as.character(data$feature_id)
  }
  
  params <- switch(method,
                   "PB" = list(target_dmz = 58.04186, rt_shift_min = -20, rt_shift_max = 20),
                   "mCPBA" = list(target_dmz = 15.9949, rt_shift_min = -120, rt_shift_max = -10),
                   "DMDS" = list(target_dmz = 93.9884, rt_shift_min = 30, rt_shift_max = 180)
  )
  
  native_dt <- copy(data)
  deriv_dt <- copy(data)
  
  # native_search <- native_dt[, .(native_id = feature_id, mz_lower, mz_upper, rt_lower, rt_upper, 
  #                                native_mz = mz, native_rt = rt, native_adduct = adduct)]
  # deriv_target <- deriv_dt[, .(deriv_id = feature_id, mz, rt, deriv_adduct = adduct)]
  # 
  # # Fast pairing based on mz and rt tolerances
  # pairs <- deriv_target[native_search, 
  #                       .(native_id = i.native_id, deriv_id = x.deriv_id, 
  #                         native_mz = i.native_mz, deriv_mz = x.mz, 
  #                         native_rt = i.native_rt, deriv_rt = x.rt, 
  #                         native_adduct = i.native_adduct, deriv_adduct = x.deriv_adduct),
  #                       on = .(mz >= mz_lower, mz <= mz_upper, 
  #                              rt >= rt_lower, rt <= rt_upper),
  #                       nomatch = NULL]
  # 
  
  
  native_dt[, mz_target := mz + params$target_dmz]
  native_dt[, mz_lower := mz_target - (mz_target * mz_tol_ppm / 1e6)]
  native_dt[, mz_upper := mz_target + (mz_target * mz_tol_ppm / 1e6)]
  native_dt[, rt_lower := rt + params$rt_shift_min]
  native_dt[, rt_upper := rt + params$rt_shift_max]
  
  native_search <- native_dt[, .(native_id = feature_id, mz_lower, mz_upper, rt_lower, rt_upper, 
                                 native_mz = mz, native_rt = rt, native_adduct = adduct, native_class = lipid_annotation)]
  deriv_target <- deriv_dt[, .(deriv_id = feature_id, mz, rt, deriv_adduct = adduct)]
  
  pairs <- deriv_target[native_search, 
                        .(native_id = i.native_id, deriv_id = x.deriv_id, 
                          native_mz = i.native_mz, deriv_mz = x.mz, 
                          native_rt = i.native_rt, deriv_rt = x.rt,
                          native_adduct = i.native_adduct, deriv_adduct = x.deriv_adduct,
                          native_class = i.native_class),
                        on = .(mz >= mz_lower, mz <= mz_upper, 
                               rt >= rt_lower, rt <= rt_upper),
                        nomatch = NULL]
  
  if (nrow(pairs) == 0) return(data.table())
  
  # A. MS1 Mass Accuracy (Max 30 points)
  pairs[, mz_error_ppm := abs(deriv_mz - (native_mz + params$target_dmz)) / (native_mz + params$target_dmz) * 1e6]
  pairs[, score_mz := 30 * pmax(0, 1 - (mz_error_ppm / mz_tol_ppm))]
  
  # B. MS1 RT Shift (Max 15 points)
  # Calculate dynamic RT score based on relative lipid mass (LipidIN-inspired concept)
  pairs[, score_rt := calculate_dynamic_rt_score(native_mz, native_rt, deriv_mz, deriv_rt, method)/2]
  
  # pairs[, rt_shift := deriv_rt - native_rt] 
  # if (method == "PB") {
  #   pairs[, score_rt := 35 * exp(-abs(rt_shift) / 10)]
  # } else if (method == "mCPBA") {
  #   pairs[, score_rt := fcase(rt_shift >= -90 & rt_shift <= -12, 35,
  #                             rt_shift < -90 | (rt_shift > -12 & rt_shift <= 0), 15,
  #                             rt_shift > 0, -35, default = 0)]
  # } else if (method == "DMDS") {
  #   pairs[, score_rt := fcase(rt_shift >= 30 & rt_shift <= 120, 35,
  #                             rt_shift > 120 | (rt_shift >= 0 & rt_shift < 30), 15,
  #                             rt_shift < 0, -35, default = 0)]
  # }
  
  # C. Cohort correlation score (Max 15)
  if (!is.null(intensity_cols) && length(intensity_cols) > 5) {
    pairs[, correlation := sapply(1:nrow(pairs), function(k) {
      v1 <- int_mat[as.character(native_id[k]), ]
      v2 <- int_mat[as.character(deriv_id[k]), ]
      if (sd(v1) == 0 || sd(v2) == 0) return(0)
      return(cor(v1, v2))
    })]
    pairs[, score_cor := fcase(
      correlation >= 0.8, 15,
      correlation >= 0.5, 5,
      correlation < 0, -15,
      default = 0
    )]
  } else if (!is.null(intensity_cols) && length(intensity_cols) <= 5) {
    # C2. Small-sample Co-detection Check (Max 10)
    pairs[, score_cor := sapply(1:nrow(pairs), function(k) {
      v1 <- int_mat[as.character(native_id[k]), ]
      v2 <- int_mat[as.character(deriv_id[k]), ]
      
      both_detected <- sum(v1 > 0 & v2 > 0)
      native_detected <- sum(v1 > 0)
      
      # If native peak is completely undetected, return 0
      if (native_detected == 0) return(0)
      
      # Reward 10 points if co-detection rate is >= 75%
      if (both_detected / native_detected >= 0.75) {
        return(10)
      } else if (both_detected / native_detected >= 0.66) {
        return(7)
      } else if (both_detected / native_detected >= 0.5) {
        return(5)
      } else if (both_detected / native_detected >= 0.33) {
        return(3)
      } else {
        return(0)
      }
    })]
  } else {
    pairs[, score_cor := 0] 
  }
  
  # D. Adduct co-occurrence score (Max 10)
  pairs[, score_adduct := ifelse(!is.na(native_adduct) & native_adduct != "" & 
                                   native_adduct == deriv_adduct, 10, 0)]
  
  # E. Hybrid MS2 Scoring (10 points Headgroup + 20 points Cosine)
  pairs[, score_headgroup := sapply(1:nrow(pairs), function(k) {
    s2 <- ms2_list[[ as.integer(deriv_id[k]) ]]
    found <- check_headgroup_presence(native_class[k], s2, deriv_mz[k], ion_mode = ion_mode, mz_tol = 0.02)
    return(ifelse(found, 10, 0))
  })]
  
  pairs[, score_cosine := sapply(1:nrow(pairs), function(k) {
    s1 <- ms2_list[[ as.integer(native_id[k]) ]]
    s2 <- ms2_list[[ as.integer(deriv_id[k]) ]]
    sim <- calculate_modified_cosine(s1, s2, params$target_dmz, mz_tol = 0.02)
    return(20 * sim) # Max 20 points
  })]
  
  pairs[, score_ms2 := score_headgroup + score_cosine]
  
  # Total Score integration (Max 100)
  pairs[, total_score := pmax(0, pmin(100, score_mz + score_rt + score_cor + score_adduct + score_ms2))]
  
  pairs[, confidence_tier := fcase(
    total_score > 85, "Tier 1 (High Confidence)",
    total_score > 60, "Tier 2 (Moderate Confidence)",
    total_score > 40, "Tier 3 (Warning: Edge Case)",
    default = "Tier 4 (Discarded)"
  )]
  
  return(pairs[total_score > 40][order(-total_score)])
}

#' Predict Equivalent Carbon Number (ECN) based on basic RT vs MZ linear approximation
#' (Simplified integration of LipidIN's relative RT rule concepts)
#' 
#' @param mz Numeric vector of mass-to-charge ratios
#' @param rt Numeric vector of retention times
#' @return A numeric vector representing the predicted ECN score penalty
calculate_dynamic_rt_score <- function(native_mz, native_rt, deriv_mz, deriv_rt, method) {
  
  rt_shift <- deriv_rt - native_rt
  
  # Base RT score matrix
  score_rt <- numeric(length(native_mz))
  
  if (method == "mCPBA") {
    # mCPBA adds an oxygen. ECN drops, polarity increases, RT shifts earlier.
    # LipidIN logic: The shift is NOT static. Larger lipids (higher mz) have stronger 
    # hydrophobic retention, so the relative impact of one oxygen is smaller.
    # Thus, delta_RT should be inversely proportional to native_mz.
    
    expected_shift_base <- -80 # Base expected shift in seconds
    mass_correction_factor <- 800 / native_mz # Normalize around typical lipid mass (e.g., 800 Da)
    
    expected_rt_shift <- expected_shift_base * mass_correction_factor
    
    # Calculate penalty based on deviation from the dynamically predicted shift
    shift_deviation <- abs(rt_shift - expected_rt_shift)
    
    # Gaussian-like decay for scoring
    score_rt <- 30 * exp(-(shift_deviation^2) / (2 * 20^2))
    
    # Absolute veto: if it shifts later, it violates fundamental chromatography rules
    score_rt[rt_shift > 0] <- -30
    
  } else if (method == "DMDS") {
    # DMDS adds two -SCH3 groups. ECN increases significantly, RT shifts later.
    # Relative impact rule: smaller lipids experience a massive RT jump, 
    # larger lipids experience a slightly less dramatic relative jump.
    
    expected_shift_base <- 120 
    mass_correction_factor <- 800 / native_mz 
    
    expected_rt_shift <- expected_shift_base * mass_correction_factor
    shift_deviation <- abs(rt_shift - expected_rt_shift)
    
    score_rt <- 30 * exp(-(shift_deviation^2) / (2 * 30^2))
    score_rt[rt_shift < 0] <- -30
    
  } else if (method == "PB") {
    # PB adds a cyclobutane ring. ECN remains relatively stable.
    # Shift should be near zero regardless of mass.
    score_rt <- 30 * exp(-abs(rt_shift) / 10)
  }
  
  return(score_rt)
}



#' MS1 Derivatization Feature Pairing and Scoring Function
#'
#' @param input_data Data frame or data.table containing mz, rt, adduct, and intensity columns
#' @param method Derivatization method: "PB", "mCPBA", or "DMDS"
#' @param intensity_cols Vector of column names containing sample intensities (for correlation)
#' @param mz_tol_ppm Mass tolerance in ppm, default is 10
#' @return A data.table containing paired features, scores, and confidence tiers
find_derivatization_pairs_ms1_only <- function(input_data, 
                                               method = c("PB", "mCPBA", "DMDS"), 
                                               intensity_cols = NULL, 
                                               mz_tol_ppm = 10) {
  
  require(data.table)
  
  method <- match.arg(method)
  data <- as.data.table(copy(input_data))
  
  # Auto-generate feature_id if missing
  if (!"feature_id" %in% names(data)) {
    data[, feature_id := 1:.N]
  }
  
  # Handle NA values in intensity matrix to prevent correlation calculation errors
  if (!is.null(intensity_cols)) {
    for (col in intensity_cols) {
      set(data, which(is.na(data[[col]])), col, 0)
    }
    int_mat <- as.matrix(data[, ..intensity_cols])
    rownames(int_mat) <- as.character(data$feature_id)
  }
  
  # ==========================================
  # 1. Configure parameters based on method
  # ==========================================
  # Note: rt_shift_min and rt_shift_max are defined in seconds
  params <- switch(method,
                   "PB" = list(
                     target_dmz = 58.04186,
                     rt_shift_min = -20,   
                     rt_shift_max = 20
                   ),
                   "mCPBA" = list(
                     target_dmz = 15.9949,
                     rt_shift_min = -120,  
                     rt_shift_max = -10
                   ),
                   "DMDS" = list(
                     target_dmz = 93.9884,
                     rt_shift_min = 30,    
                     rt_shift_max = 180
                   )
  )
  
  # ==========================================
  # 2. Funnel Layer 1: Fast Non-equi join
  # ==========================================
  native_dt <- copy(data)
  deriv_dt <- copy(data)
  
  native_dt[, mz_target := mz + params$target_dmz]
  native_dt[, mz_lower := mz_target - (mz_target * mz_tol_ppm / 1e6)]
  native_dt[, mz_upper := mz_target + (mz_target * mz_tol_ppm / 1e6)]
  native_dt[, rt_lower := rt + params$rt_shift_min]
  native_dt[, rt_upper := rt + params$rt_shift_max]
  
  native_search <- native_dt[, .(native_id = feature_id, mz_lower, mz_upper, rt_lower, rt_upper, 
                                 native_mz = mz, native_rt = rt, native_adduct = adduct)]
  deriv_target <- deriv_dt[, .(deriv_id = feature_id, mz, rt, deriv_adduct = adduct)]
  
  # Fast pairing based on mz and rt tolerances
  pairs <- deriv_target[native_search, 
                        .(native_id = i.native_id, deriv_id = x.deriv_id, 
                          native_mz = i.native_mz, deriv_mz = x.mz, 
                          native_rt = i.native_rt, deriv_rt = x.rt, 
                          native_adduct = i.native_adduct, deriv_adduct = x.deriv_adduct),
                        on = .(mz >= mz_lower, mz <= mz_upper, 
                               rt >= rt_lower, rt <= rt_upper),
                        nomatch = NULL]
  
  if (nrow(pairs) == 0) {
    message("No matching pairs found.")
    return(data.table())
  }
  
  # ==========================================
  # 3. Orthogonal Validation & Scoring Engine
  # ==========================================
  pairs[, mz_error_ppm := abs(deriv_mz - (native_mz + params$target_dmz)) / (native_mz + params$target_dmz) * 1e6]
  pairs[, rt_shift := deriv_rt - native_rt] 
  
  # A. Mass accuracy score (Max 40, linear decay)
  pairs[, score_mz := 40 * pmax(0, 1 - (mz_error_ppm / mz_tol_ppm))]
  
  # B. Method-specific RT shift score (Max 30)
  # Calculate dynamic RT score based on relative lipid mass (LipidIN-inspired concept)
  pairs[, score_rt := calculate_dynamic_rt_score(native_mz, native_rt, deriv_mz, deriv_rt, method)]
  
  # if (method == "PB") {
  #   pairs[, score_rt := 30 * exp(-abs(rt_shift) / 10)]
  # } else if (method == "mCPBA") {
  #   pairs[, score_rt := fcase(
  #     rt_shift >= -90 & rt_shift <= -12, 30,
  #     rt_shift < -90 | (rt_shift > -12 & rt_shift <= 0), 15,
  #     rt_shift > 0, -30,
  #     default = 0
  #   )]
  # } else if (method == "DMDS") {
  #   pairs[, score_rt := fcase(
  #     rt_shift >= 30 & rt_shift <= 120, 30,
  #     rt_shift > 120 | (rt_shift >= 0 & rt_shift < 30), 15,
  #     rt_shift < 0, -30,
  #     default = 0
  #   )]
  # }
  
  # C. Cohort correlation score (Max 20)
  if (!is.null(intensity_cols) && length(intensity_cols) > 5) {
    pairs[, correlation := sapply(1:nrow(pairs), function(k) {
      v1 <- int_mat[as.character(native_id[k]), ]
      v2 <- int_mat[as.character(deriv_id[k]), ]
      if (sd(v1) == 0 || sd(v2) == 0) return(0)
      return(cor(v1, v2))
    })]
    pairs[, score_cor := fcase(
      correlation >= 0.8, 20,
      correlation >= 0.5, 10,
      correlation < 0, -20,
      default = 0
    )]
  } else if (!is.null(intensity_cols) && length(intensity_cols) <= 5) {
    # C2. Small-sample Co-detection Check (Max 10)
    pairs[, score_cor := sapply(1:nrow(pairs), function(k) {
      v1 <- int_mat[as.character(native_id[k]), ]
      v2 <- int_mat[as.character(deriv_id[k]), ]
      
      both_detected <- sum(v1 > 0 & v2 > 0)
      native_detected <- sum(v1 > 0)
      
      # If native peak is completely undetected, return 0
      if (native_detected == 0) return(0)
      
      # Reward 10 points if co-detection rate is >= 75%
      if (both_detected / native_detected >= 0.75) {
        return(10)
      } else if (both_detected / native_detected >= 0.66) {
        return(7)
      } else if (both_detected / native_detected >= 0.5) {
        return(5)
      } else if (both_detected / native_detected >= 0.33) {
        return(3)
      } else {
        return(0)
      } 
    })]
  } else {
    pairs[, score_cor := 0] 
  }
  
  # D. Adduct co-occurrence score (Max 10)
  pairs[, score_adduct := ifelse(!is.na(native_adduct) & native_adduct != "" & 
                                   native_adduct == deriv_adduct, 10, 0)]
  
  # ==========================================
  # 4. Tier Classification & Output
  # ==========================================
  # Limit total score between 0 and 100
  pairs[, total_score := pmax(0, pmin(100, score_mz + score_rt + score_cor + score_adduct))]
  
  # Assign confidence tiers
  pairs[, confidence_tier := fcase(
    total_score > 80, "Tier 1 (High Confidence)",
    total_score > 60, "Tier 2 (Moderate Confidence)",
    total_score > 40, "Tier 3 (Warning: Edge Case)",
    default = "Tier 4 (Discarded)"
  )]
  
  # Filter low confidence pairs and sort descending by total_score
  final_results <- pairs[total_score > 0][order(-total_score)]
  
  return(final_results)
}

{
  
  # ==============================================================================
  # Comprehensive C=C Localization Pipeline (Platform-Independent & Multi-Chain)
  # ==============================================================================
  
  # ------------------------------------------------------------------------------
  # 1. CONSTANTS & EXACT MASS DICTIONARY
  # ------------------------------------------------------------------------------
  .exact_mass <- list(
    carbon = 12.00000000000,
    hydrogen = 1.00782503223,
    nitrogen = 14.00307400443,
    oxygen = 15.99491461957,
    sulfur = 31.97207100000 
  )
  
  
  # ------------------------------------------------------------------------------
  # 2. HELPER FUNCTIONS (Parsing, Merging, & Scoring)
  # ------------------------------------------------------------------------------
  extract_carbons <- function(annotation) {
    if (is.na(annotation) || annotation == "") return(0)
    matches <- regmatches(annotation, gregexpr("[0-9]+(?=:)", annotation, perl = TRUE))
    if (length(matches[[1]]) == 0) return(0)
    return(sum(as.numeric(matches[[1]])))
  }
  
  extract_db_count <- function(annotation) {
    if (is.na(annotation) || annotation == "") return(0)
    matches <- regmatches(annotation, gregexpr("(?<=:)[0-9]+", annotation, perl = TRUE))
    if (length(matches[[1]]) == 0) return(0)
    return(sum(as.numeric(matches[[1]])))
  }
  
  merge_theoretical_peaks <- function(peaks_dt, mz_tol = 0.005) {
    if (nrow(peaks_dt) == 0) return(peaks_dt)
    peaks_dt <- peaks_dt[order(mz)]
    merged_peaks <- list()
    current_peak <- peaks_dt[1]
    
    if (nrow(peaks_dt) > 1) {
      for (i in 2:nrow(peaks_dt)) {
        if (abs(peaks_dt$mz[i] - current_peak$mz) <= mz_tol) {
          current_peak$theo_int <- current_peak$theo_int + peaks_dt$theo_int[i]
          current_peak$comment <- paste(current_peak$comment, peaks_dt$comment[i], sep = " & ")
        } else {
          merged_peaks[[length(merged_peaks) + 1]] <- current_peak
          current_peak <- peaks_dt[i]
        }
      }
    }
    merged_peaks[[length(merged_peaks) + 1]] <- current_peak
    return(rbindlist(merged_peaks))
  }
  
  calculate_spectral_match <- function(exp_spec, theo_spec, mz_tol = 0.02) {
    if (nrow(exp_spec) == 0 || nrow(theo_spec) == 0) return(0)
    setnames(exp_spec, c("mz", "intensity"))
    exp_spec[, norm_int := (intensity / max(intensity)) * 100]
    
    matched_intensity_sum <- 0
    for (i in 1:nrow(theo_spec)) {
      matches <- exp_spec[abs(mz - theo_spec$mz[i]) <= mz_tol]
      if (nrow(matches) > 0) {
        best_match <- matches[order(-norm_int)][1]
        matched_intensity_sum <- matched_intensity_sum + (theo_spec$theo_int[i] * best_match$norm_int)
      }
    }
    exp_intensity_sum <- sum(exp_spec$norm_int^2)
    if (exp_intensity_sum == 0) return(0)
    return(matched_intensity_sum / sqrt(sum(theo_spec$theo_int^2) * exp_intensity_sum))
  }
  
  # ------------------------------------------------------------------------------
  # 3. ROUTE A: DIFFERENTIAL MS2 SUBTRACTION
  # ------------------------------------------------------------------------------
  # ==============================================================================
  # ROUTE A: DIFFERENTIAL MS2 SUBTRACTION (Upgraded for PB, DMDS, mCPBA)
  # ==============================================================================
  #' Differential MS2 Spectral Subtraction for C=C Localization
  #' 
  #' @param spec_native data.table, native MS2 spectrum
  #' @param spec_deriv data.table, derivatized MS2 spectrum
  #' @param precursor_mz Numeric, accurate derivatized precursor m/z
  #' @param reagent_mass Numeric, exact mass of the derivatization tag
  #' @param mz_tol Numeric, m/z tolerance for peak matching
  #' @param method Character, derivatization method ("PB", "DMDS", "mCPBA")
  #' @return A list containing status and annotation result
  differential_ms2_localization <- function(spec_native, spec_deriv, precursor_mz, 
                                            reagent_mass, mz_tol = 0.02, method = "PB") {
    
    s_nat <- as.data.table(spec_native)
    s_der <- as.data.table(spec_deriv)
    setnames(s_nat, c("mz", "int"))
    setnames(s_der, c("mz", "int"))
    
    # Filter noise (1% threshold)
    s_nat <- s_nat[int > max(int) * 0.01]
    s_der <- s_der[int > max(int) * 0.01]
    s_der[, is_background := FALSE]
    
    # Subtract background (Native peaks + Intact shifted peaks)
    for (i in 1:nrow(s_der)) {
      target_mz <- s_der$mz[i]
      direct_match <- any(abs(s_nat$mz - target_mz) <= mz_tol)
      shifted_match <- any(abs((s_nat$mz + reagent_mass) - target_mz) <= mz_tol)
      if (direct_match || shifted_match) s_der$is_background[i] <- TRUE
    }
    
    # Isolate unique diagnostic fragments
    delta_spec <- s_der[is_background == FALSE]
    
    if (nrow(delta_spec) == 0) {
      return(list(status = "Unconfirmed", annotation = "No unique diagnostic fragments survived subtraction."))
    }
    
    status <- "Unconfirmed"
    annotation_result <- "Unique fragments found, but rules did not match."
    
    ch2_mass <- 14.01565
    
    if (method == "PB") {
      # PB produces a pair of fragments separated by ~25.9793 Da
      pb_gap <- 25.9793
      n_peaks <- nrow(delta_spec)
      
      if (n_peaks >= 2) {
        grid <- as.data.table(expand.grid(i = 1:(n_peaks-1), j = 2:n_peaks))
        grid <- grid[i < j]
        grid[, mz1 := delta_spec$mz[i]]
        grid[, mz2 := delta_spec$mz[j]]
        grid[, gap := mz2 - mz1]
        
        valid_pairs <- grid[abs(gap - pb_gap) <= mz_tol]
        if (nrow(valid_pairs) > 0) {
          best_pair <- valid_pairs[order(- (delta_spec$int[i] + delta_spec$int[j]))][1]
          status <- "Confirmed (Differential)"
          annotation_result <- sprintf("PB Pair isolated: m/z %.4f and %.4f", best_pair$mz1, best_pair$mz2)
        }
      }
      
    } else if (method == "DMDS") {
      # DMDS neutral loss corresponds to the distal tail + SCH3 group (47.0192 Da)
      dmds_offset <- 47.0192 
      delta_spec[, nl_mass := precursor_mz - mz]
      
      # Estimate the number of carbons in the distal tail (omega)
      delta_spec[, est_omega := round((nl_mass - dmds_offset) / ch2_mass)]
      
      # Filter for chemically meaningful omega values (e.g., 1 to 20)
      valid_frags <- delta_spec[est_omega > 0 & est_omega < 25]
      
      if (nrow(valid_frags) > 0) {
        best_frag <- valid_frags[order(-int)][1]
        status <- "Confirmed (Differential)"
        annotation_result <- sprintf("DMDS Cleavage: m/z %.4f (Est. omega-tail: w=%d)", 
                                     best_frag$mz, best_frag$est_omega)
      }
      
    } else if (method == "mCPBA") {
      # mCPBA neutral loss corresponds to the distal tail + Oxygen - Hydrogen (14.9871 Da)
      mcpba_offset <- 14.9871 
      delta_spec[, nl_mass := precursor_mz - mz]
      
      # Estimate the number of carbons in the distal tail (omega)
      delta_spec[, est_omega := round((nl_mass - mcpba_offset) / ch2_mass)]
      
      valid_frags <- delta_spec[est_omega > 0 & est_omega < 25]
      
      if (nrow(valid_frags) > 0) {
        best_frag <- valid_frags[order(-int)][1]
        status <- "Confirmed (Differential)"
        annotation_result <- sprintf("mCPBA Aldehyde: m/z %.4f (Est. omega-tail: w=%d)", 
                                     best_frag$mz, best_frag$est_omega)
      }
    }
    
    return(list(status = status, annotation = annotation_result))
  }
  # ------------------------------------------------------------------------------
  # 4. ROUTE B: IN-SILICO MULTI-CHAIN PREDICTION ENGINE
  # ------------------------------------------------------------------------------
  generate_biological_isomers <- function(chain_carbons, db_count) {
    isomers <- list()
    if (db_count == 0 || chain_carbons < 4) return(isomers)
    
    if (db_count == 1) {
      valid_starts <- 4:(chain_carbons - 2)
      for (pos in valid_starts) isomers[[paste0("Delta_", pos)]] <- c(pos)
      return(isomers)
    }
    
    span <- (db_count - 1) * 3
    max_start <- chain_carbons - span - 2
    
    if (max_start >= 4) {
      for (start_pos in 4:max_start) {
        pos_sequence <- seq(from = start_pos, by = 3, length.out = db_count)
        isomers[[paste0("Delta_", paste(pos_sequence, collapse = ","))]] <- pos_sequence
      }
    }
    return(isomers)
  }
  
  unroll_sum_composition <- function(total_c, total_db, min_chain = 10, max_chain = 32) {
    valid_combinations <- list()
    
    # 拦截判定：仅当总碳数大于等于两条最小链之和时，才尝试进行双链拆解
    # 对于 FA 10:1 (total_c = 10)，这里 10 < 24，直接跳过循环，完美避免报错
    if (total_c >= (min_chain * 2)) {
      for (c1 in seq(min_chain, min(max_chain, total_c - min_chain), by = 2)) {
        c2 <- total_c - c1
        
        if (c2 %% 2 != 0 || c2 < min_chain || c2 > max_chain || c1 > c2) next
        
        for (db1 in 0:total_db) {
          db2 <- total_db - db1
          
          if (db1 > 6 || db2 > 6) next
          if (db1 > 0 && c1 < (db1 * 3 + 3)) next 
          if (db2 > 0 && c2 < (db2 * 3 + 3)) next
          
          combo_dt <- data.table(
            chain_id = c(sprintf("FA%d:%d", c1, db1), sprintf("FA%d:%d", c2, db2)),
            carbons = c(c1, c2),
            db_count = c(db1, db2)
          )
          label <- sprintf("C%d:%d_C%d:%d", c1, db1, c2, db2)
          valid_combinations[[label]] <- combo_dt
        }
      }
    }
    
    # Fallback for mono-acyl lipids (FA, LPC, LPE, etc.) or if sum logic fails
    if (length(valid_combinations) == 0) {
      combo_dt <- data.table(
        chain_id = sprintf("FA%d:%d", total_c, total_db),
        carbons = total_c, 
        db_count = total_db
      )
      valid_combinations[[sprintf("Single_C%d:%d", total_c, total_db)]] <- combo_dt
    }
    
    return(valid_combinations)
  }
  generate_multichain_theoretical_spectrum <- function(precursor_mz, chains_dt, hypothesis_list, 
                                                       lipid_class = "Unknown",ion_mode = "POS",
                                                       headgroup_mz = NA_real_, method = "PB") {
    peaks <- list()
    d <- .exact_mass
    
    peaks[[1]] <- data.table(mz = precursor_mz, theo_int = 20, comment = "Precursor")
    # if (!is.na(headgroup_mz)) {
    #   peaks[[2]] <- data.table(mz = headgroup_mz, theo_int = 100, comment = "Polar Headgroup")
    # }
    # Step 2: Dynamically add lipid class-specific backbone fragments
    backbone_peaks <- generate_class_backbone_peaks(precursor_mz, lipid_class, ion_mode, chains_dt, method)
    if (length(backbone_peaks) > 0) {
      peaks <- c(peaks, backbone_peaks)
    }
    
    for (i in 1:nrow(chains_dt)) {
      c_len <- chains_dt$carbons[i]
      c_db <- chains_dt$db_count[i]
      c_id <- chains_dt$chain_id[i]
      if (c_db == 0) next 
      
      db_pos <- sort(unlist(hypothesis_list[[c_id]]))
      
      for (j in seq_along(db_pos)) {
        target_bond <- db_pos[j]
        omega <- c_len - target_bond
        eta <- length(db_pos) - j + 1 
        base_tail_mass <- (omega * d$carbon) + ((2 * omega + 1 - 2 * (eta - 1)) * d$hydrogen)
        
        if (method == "PB") {
          mass_NO <- (omega * d$carbon) + ((2 * omega + 2 - 2 * eta) * d$hydrogen) + d$oxygen
          mass_NA <- ((omega + 3) * d$carbon) + ((2 * omega - 2 * eta + 8) * d$hydrogen)
          peaks[[length(peaks) + 1]] <- data.table(mz = precursor_mz - mass_NO, theo_int = 100, 
                                                   comment = sprintf("[%s] PB FA(w=%d) C%d", c_id, omega, target_bond))
          peaks[[length(peaks) + 1]] <- data.table(mz = precursor_mz - mass_NA, theo_int = 80, 
                                                   comment = sprintf("[%s] PB FO(w=%d) C%d", c_id, omega, target_bond))
        } else if (method == "DMDS") {
          sch3_mass <- d$carbon + (3 * d$hydrogen) + d$sulfur
          peaks[[length(peaks) + 1]] <- data.table(mz = precursor_mz - (base_tail_mass + sch3_mass), theo_int = 100, 
                                                   comment = sprintf("[%s] DMDS Cleavage C%d", c_id, target_bond))
        } else if (method == "mCPBA") {
          aldehyde_shift <- d$oxygen - d$hydrogen
          peaks[[length(peaks) + 1]] <- data.table(mz = precursor_mz - (base_tail_mass + aldehyde_shift), theo_int = 100, 
                                                   comment = sprintf("[%s] mCPBA Aldehyde C%d", c_id, target_bond))
        }
      }
    }
    return(merge_theoretical_peaks(rbindlist(peaks), mz_tol = 0.005))
  }
  
  #' Predict and Score Theoretical C=C Isomers against Experimental MS2
  predict_and_score_cc_isomers <- function(exp_ms2, precursor_mz, total_c, total_db,
                                           lipid_class = "Unknown",ion_mode = "POS",
                                           headgroup_mz = NA_real_, method = "PB") {
    # Step 1: Unroll generic "PC 48:6" into potential diacyl topologies
    chain_hypotheses <- unroll_sum_composition(total_c, total_db)
    
    if (length(chain_hypotheses) == 0) {
      return("Prediction Failed: No valid biological chains generated")
    }
    
    results <- list()
    
    # Step 2: Iterate over all possible chain length combinations
    for (combo_name in names(chain_hypotheses)) {
      chains_dt <- chain_hypotheses[[combo_name]]
      db_hypotheses <- list()
      
      # Step 3: Generate positional isomers for each chain
      for (i in 1:nrow(chains_dt)) {
        c_id <- chains_dt$chain_id[i]
        if (chains_dt$db_count[i] > 0) {
          db_hypotheses[[c_id]] <- generate_biological_isomers(chains_dt$carbons[i], chains_dt$db_count[i])
        } else {
          db_hypotheses[[c_id]] <- list("Saturated" = numeric(0))
        }
      }
      
      # Step 4: Create a grid of all topological combinations
      isomer_grid <- expand.grid(db_hypotheses, stringsAsFactors = FALSE)
      
      # Step 5: Generate theoretical spectrum for each combination and match
      for (row_idx in 1:nrow(isomer_grid)) {
        current_test <- as.list(isomer_grid[row_idx, , drop = FALSE])
        names(current_test) <- names(db_hypotheses)
        
        theo_spec <- generate_multichain_theoretical_spectrum(precursor_mz, chains_dt, current_test, lipid_class, ion_mode, headgroup_mz, method)
        score <- calculate_spectral_match(as.data.table(exp_ms2), theo_spec)
        
        specific_label <- paste(sapply(names(current_test), function(x) {
          paste0(x, "(", paste(current_test[[x]], collapse = ","), ")")
        }), collapse = " | ")
        
        results[[length(results) + 1]] <- data.table(Topology = paste(combo_name, specific_label, sep = " -> "), Score = round(score * 100, 2))
      }
    }
    
    # Step 6: Rank candidates
    final_results <- rbindlist(results)
    final_results <- final_results[order(-Score)]
    
    if (nrow(final_results) == 0 || final_results$Score[1] < 10) return("Unresolved / Score Too Low")
    return(sprintf("Predicted: %s (Cosine: %.2f)", final_results$Topology[1], final_results$Score[1]))
  }
  
  # ------------------------------------------------------------------------------
  # 5. ROUTE C: NATIVE HEADGROUP VALIDATION
  # ------------------------------------------------------------------------------
  check_headgroup_presence2 <- function(lipid_annotation, ms2_spec, precursor_mz, ion_mode = "POS", mz_tol = 0.02) {
    if (is.null(ms2_spec) || nrow(ms2_spec) == 0 || is.na(lipid_annotation) || lipid_annotation == "") return(FALSE)
    main_class <- strsplit(trimws(lipid_annotation), " ")[[1]][1]
    rules <- get_lipid_rules()
    
    for (rule in rules) {
      if (rule$mode == ion_mode && grepl(rule$pattern, main_class, ignore.case = TRUE)) {
        if (rule$type == "fragment" && any(abs(ms2_spec[[1]] - rule$val) <= mz_tol)) return(TRUE)
        if (rule$type == "nl" && any(abs(ms2_spec[[1]] - (precursor_mz - rule$val)) <= mz_tol)) return(TRUE)
      }
    }
    return(FALSE)
  }
  
  # ------------------------------------------------------------------------------
  # 6. ADAPTIVE ROUTER & MAIN ENTRY PIPELINE
  # ------------------------------------------------------------------------------
  #' Adaptive Engine to Route C=C Localization based on Data Availability
  adaptive_cc_localization <- function(ms1_pair, spec_native, spec_deriv,lipid_class,ion_mode, method = "PB") {
    has_native <- !is.null(spec_native) && nrow(spec_native) >= 3
    has_deriv <- !is.null(spec_deriv) && nrow(spec_deriv) >= 3
    
    result <- list(engine_used = "None", cc_position = "Unknown", confidence = "Tier 4 (MS1 Only)")
    
    total_c <- extract_carbons(ms1_pair$native_class) 
    total_db <- extract_db_count(ms1_pair$native_class)
    
    # Pre-fetch theoretical headgroup m/z for prediction module (if class known)
    # For simplicity, assumed NA here, but can be linked to exact rules
    headgroup_mz_val <- NA_real_ 
    
    # ROUTE A: Differential Subtraction (Both MS2 present)
    if (has_native && has_deriv) {
      reagent_mass <- switch(method, "PB" = 58.04186, "DMDS" = 93.9884, "mCPBA" = 15.9949)
      
      # Notice we pass ms1_pair$deriv_mz to enable DMDS/mCPBA neutral loss calculations
      diff_result <- differential_ms2_localization(
        spec_native = spec_native, 
        spec_deriv = spec_deriv, 
        precursor_mz = ms1_pair$deriv_mz, 
        reagent_mass = reagent_mass, 
        method = method
      )
      
      if (grepl("Confirmed", diff_result$status)) {
        result$engine_used <- "Route A: Differential Subtraction"
        result$cc_position <- diff_result$annotation
        result$confidence <- "Tier 1a (Platinum)"
        return(result)
      }
    }
    
    # ROUTE B: In-Silico Prediction (Deriv MS2 only, or Route A fallback)
    if (has_deriv && total_c > 0 && total_db > 0) {
      pred_result <- predict_and_score_cc_isomers(
        exp_ms2 = spec_deriv, precursor_mz = ms1_pair$deriv_mz, 
        total_c = total_c, total_db = total_db, 
        lipid_class = lipid_class,ion_mode = ion_mode,
        headgroup_mz = headgroup_mz_val, method = method
      )
      
      if (grepl("Predicted", pred_result)) {
        result$engine_used <- "Route B: In-Silico Prediction"
        result$cc_position <- pred_result
        result$confidence <- "Tier 1b (Gold Prediction)"
        return(result)
      }
    }
    
    # ROUTE C: Native Headgroup Check Only (Missing Deriv MS2)
    if (has_native) {
      headgroup_found <- check_headgroup_presence2(ms1_pair$native_class, spec_native, ms1_pair$native_mz)
      result$engine_used <- "Route C: Headgroup Validator"
      result$cc_position <- "C=C Unknown (Missing Deriv MS2)"
      result$confidence <- ifelse(headgroup_found, "Tier 2 (Class Confirmed)", "Tier 3 (Unconfirmed)")
      return(result)
    }
    
    return(result)
  }
  
  #' Master Pipeline Entry Function
  #' @param matched_pairs data.table from MS1 matching
  #' @param ms2_list Named list containing all MS2 spectra
  #' @param method Character, "PB", "mCPBA", or "DMDS"
  run_cc_annotation_pipeline <- function(matched_pairs, ms2_list, ion_mode = "POS", method = c("PB", "mCPBA", "DMDS")) {
    method <- match.arg(method)
    results_dt <- copy(as.data.table(matched_pairs))
    
    results_dt[, c("cc_engine", "cc_position", "cc_confidence") := .("Pending", "Pending", "Pending")]
    
    for (i in 1:nrow(results_dt)) {
      row_data <- results_dt[i]
      
      s_native <- ms2_list[[as.integer(row_data$native_id)]]
      s_deriv <- ms2_list[[as.integer(row_data$deriv_id)]]
      lipid_class <- strsplit(row_data$native_class," ")[[1]][1]
      
      loc_res <- adaptive_cc_localization(row_data, s_native, s_deriv, lipid_class, ion_mode, method)
      
      set(results_dt, i, "cc_engine", loc_res$engine_used)
      set(results_dt, i, "cc_position", loc_res$cc_position)
      set(results_dt, i, "cc_confidence", loc_res$confidence)
    }
    
    return(results_dt[order(cc_confidence)])
  }
  
  
  
}

# ─────────────────────────────────────────────────────────────────────────────
# GetPbBasedLipidomicsMatchedPeaksScores.R
# 
# Adapted for Paternò-Büchi (PB) chemical derivatization double bond localization
# based on the OAD algorithm framework.
# Based on LipidOA rules: PB reaction generates an oxetane ring, and characteristic
# CID cleavage produces FA and FO ions[cite: 1, 2].
# ─────────────────────────────────────────────────────────────────────────────

# --- Constants definition (based on LipidOA and common exact masses) ---
.pb_mass_diff <- list(
  carbon = 12.000000,
  hydrogen = 1.007825,
  oxygen = 15.994915,
  nitrogen = 14.003074,
  phosphorus = 30.973762
)
.pb_proton <- 1.007276

# Default PB reagent is Acetone
.pb_acetone_mass <- .pb_mass_diff$carbon * 3 + .pb_mass_diff$hydrogen * 6 + .pb_mass_diff$oxygen

# ─────────────────────────────────────────────────────────────────────────────
# Core diagnostic peak generation engine: Generate theoretical spectra for PB
# characteristic cleavages[cite: 2].
# ─────────────────────────────────────────────────────────────────────────────
.pb_sg_double_bond_peaks <- function(lipid_mass, chain, adduct, nl_mass = 0, abundance = 40, reagent = "acetone") {
  positions <- chain$positions
  if (is.list(positions)) positions <- positions[[1L]]
  if (chain$carbon <= 0L || length(positions) == 0L) return(list())
  
  result <- list()
  
  # Helper function: Add theoretical peak
  add <- function(id, neutral_mass, factor, comment_suffix) {
    # Simulate ion generation logic from the OAD script [source: 4]
    mz <- neutral_mass + .pb_proton # Simply assume positive ion mode plus proton; adjust according to adduct in practice
    result[[length(result) + 1L]] <<- data.frame(
      Mass = as.numeric(mz), 
      Intensity = as.numeric(abundance * factor),
      Comment = as.character(paste0(chain$carbon, ":", chain$double_bond, " ", comment_suffix)), 
      PbId = as.character(id),
      SpectrumComment = "doublebond", 
      stringsAsFactors = FALSE
    )
  }
  
  # Iterate through each double bond position on the chain (Delta nomenclature, counted from the carboxyl end)
  for (i in seq_along(positions)) {
    bond <- positions[i]
    if (bond <= 1L || bond >= chain$carbon) next
    
    # Convert to Omega (w) position counted from the methyl end[cite: 2]
    omega <- chain$carbon - bond
    # Calculate cumulative unsaturation up to the current double bond (eta)
    eta <- length(positions) - i + 1 
    
    # Theoretical Neutral Loss calculation formula based on LipidOA[cite: 2]
    if (reagent == "acetone") {
      # N_O: C_{omega} H_{2*omega+2-2*eta} O[cite: 2]
      mass_NO <- omega * .pb_mass_diff$carbon + 
        (2 * omega + 2 - 2 * eta) * .pb_mass_diff$hydrogen + 
        .pb_mass_diff$oxygen
      
      # N_A: C_{omega+3} H_{2*omega-2*eta+8}[cite: 2]
      mass_NA <- (omega + 3) * .pb_mass_diff$carbon + 
        (2 * omega - 2 * eta + 8) * .pb_mass_diff$hydrogen
    }
    
    # Mass of complementary diagnostic ions = precursor mass after PB derivatization - neutral loss[cite: 1]
    pb_precursor_mass <- lipid_mass + .pb_acetone_mass
    mass_FA <- pb_precursor_mass - mass_NO
    mass_FO <- pb_precursor_mass - mass_NA
    
    # Generate paired characteristic diagnostic peaks[cite: 2]
    add("FA", mass_FA, 1.0, paste0("FA(w=", omega, ")"))
    add("FO", mass_FO, 0.8, paste0("FO(w=", omega, ")"))
  }
  
  return(result)
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper function: Merge theoretical peak list [source: 4]
# ─────────────────────────────────────────────────────────────────────────────
.pb_sg_merge <- function(peaks) {
  if (length(peaks) == 0L) {
    return(data.frame(Mass = numeric(), Intensity = numeric(), Comment = character(),
                      PbId = character(), SpectrumComment = character(), stringsAsFactors = FALSE))
  }
  peaks_df <- do.call(rbind, peaks)
  peaks_df <- peaks_df[order(peaks_df$Mass), , drop = FALSE]
  return(peaks_df)
}

# ─────────────────────────────────────────────────────────────────────────────
# Lipid class specific theoretical spectrum generator (e.g., PC) [source: 4]
# ─────────────────────────────────────────────────────────────────────────────
.pb_sg_generate_class <- function(lipid, class, adduct, abundance = 40, nl_mass = 0) {
  mass <- as.numeric(lipid$NeutralMass)
  chains <- lipid$Chains # Assume parsed into a data.frame containing carbon, double_bond, positions
  
  peaks <- list()
  
  # 1. Inherit basic CID characteristic peaks (polar head, etc.)[cite: 1, 2]
  if (class == "PC" && adduct == "[M+H]+") {
    peaks[[1]] <- data.frame(Mass = 184.0733, Intensity = 999, Comment = "PC Header", 
                             PbId = NA, SpectrumComment = "metaboliteclass", stringsAsFactors = FALSE)
  }
  
  # 2. Append double bond diagnostic peaks from PB reaction[cite: 2]
  if (!is.null(chains) && nrow(chains) > 0L) {
    for (i in seq_len(nrow(chains))) {
      # Pass the mass of the natural lipid without derivatization; acetone mass is added internally
      pb_peaks <- .pb_sg_double_bond_peaks(mass, chains[i, , drop = FALSE], adduct, nl_mass, abundance)
      peaks <- c(peaks, pb_peaks)
    }
  }
  
  list(
    Name = lipid$Name, 
    CompoundClass = class,
    Spectrum = .pb_sg_merge(peaks)
  )
}

# Wrapper
generate_pc_pb_spectrum <- function(lipid, adduct, molecule = NULL) {
  .pb_sg_generate_class(lipid, "PC", adduct, 40)
}

# ─────────────────────────────────────────────────────────────────────────────
# Spectrum matching and scoring engine (equivalent to GetOadBasedLipidomicsMatchedPeaksScores) [source: 3]
# ─────────────────────────────────────────────────────────────────────────────
get_pb_based_lipidomics_matched_peaks_scores <- function(scan, reference_pb, tolerance) {
  
  measured <- scan$Spectrum
  library <- reference_pb$Spectrum
  
  if (nrow(measured) == 0L || nrow(library) == 0L) return(c(0, 0))
  
  matched_count <- 0L
  library_count <- sum(!is.na(library$PbId)) # Only evaluate the hit rate of double bond diagnostic peaks[cite: 2]
  
  if (library_count == 0L) return(c(0, 0))
  
  # Look for paired FA and FO
  for (i in seq_len(nrow(library))) {
    if (is.na(library$PbId[i])) next
    
    target_mz <- library$Mass[i]
    
    # Find matching peaks in the experimental spectrum within tolerance [source: 3]
    match_idx <- which(abs(measured$Mass - target_mz) <= tolerance)
    
    # Simple abundance threshold filtering (similar to S/N > 3 filtering in LipidOA[cite: 1])
    if (length(match_idx) > 0 && any(measured$Intensity[match_idx] > 100)) {
      matched_count <- matched_count + 1L
    }
  }
  
  c(matched_count / library_count, matched_count)
}

# ─────────────────────────────────────────────────────────────────────────────
# Top-level Dispatcher [source: 3]
# ─────────────────────────────────────────────────────────────────────────────
get_pb_based_lipid_molecule_annotation_result <- function(scan, lipid_info, adduct, tolerance) {
  
  lipid_class <- sub(" .*$", "", trimws(lipid_info$Name))
  
  # 1. Generate theoretical spectrum with PB characteristic peaks according to class [source: 4]
  reference_pb <- switch(lipid_class,
                         "PC" = generate_pc_pb_spectrum(lipid_info, adduct),
                         # "PE" = generate_pe_pb_spectrum(lipid_info, adduct),
                         NULL
  )
  
  if (is.null(reference_pb)) {
    return(list(Name = lipid_info$Name, Scores = c(0, 0), Level = 0L))
  }
  
  # 2. Match and score the experimental spectrum against the generated theoretical PB spectrum [source: 3]
  scores <- get_pb_based_lipidomics_matched_peaks_scores(scan, reference_pb, tolerance)
  
  level <- if (scores[2] >= 2) 2L else 1L # Assign Level 2 if matching criteria are met[cite: 2]
  
  return(list(
    Name = lipid_info$Name,
    Class = lipid_class,
    Scores = scores,      # scores[1] = percentage, scores[2] = matched count
    Level = level,
    TheoreticalSpectrum = reference_pb$Spectrum
  ))
}