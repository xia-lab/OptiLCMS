\name{mSet}
\alias{mSet}
\docType{data}
\title{
Trimmed IBD mSet data.
}
\description{
A small example data set for testing purposes, containing 10 trimmed spectra (UPLC-Q/E-ESI-, C18 neg) organized into three groups (Healthy, Crohn's Disease and QC) from (Lloyd-Price et al.) 
}
\usage{data("mSet")}
\format{
  The format is:
Formal class 'mSet' [package "OptiLCMS"] with 14 slots
  ..@ rawfiles        : chr(0) 
  ..@ rawOnDisk       :Formal class 'OnDiskMSnExp' [package "MSnbase"] with 10 slots
  .. .. ..@ spectraProcessingQueue: list()
  .. .. ..@ backend               : chr(0) 
  .. .. ..@ assayData             :<environment: 0x55cab1c25500> 
  .. .. ..@ phenoData             :Formal class 'NAnnotatedDataFrame' [package "MSnbase"] with 6 slots
  .. .. .. .. ..@ multiplex        : num 1
  .. .. .. .. ..@ multiLabels      : chr "Single run"
  .. .. .. .. ..@ varMetadata      :'data.frame':	2 obs. of  1 variable:
  .. .. .. .. .. ..$ labelDescription: chr [1:2] NA NA
  .. .. .. .. ..@ data             :'data.frame':	10 obs. of  2 variables:
  .. .. .. .. .. ..$ sample_name : chr [1:10] "CD-6KUCT" "CD-77FXR" "CD-9OS5Y" "CD-9WOBP" ...
  .. .. .. .. .. ..$ sample_group: Factor w/ 3 levels "Control","Disease",..: 1 1 1 1 2 2 2 2 3 3
  .. .. .. .. ..@ dimLabels        : chr [1:2] "rowNames" "columnNames"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 0 0 3
  .. .. ..@ featureData           :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. ..@ varMetadata      :'data.frame':	35 obs. of  1 variable:
  .. .. .. .. .. ..$ labelDescription: chr [1:35] NA NA NA NA ...
  .. .. .. .. ..@ data             :'data.frame':	4461 obs. of  35 variables:
  .. .. .. .. .. ..$ fileIdx                   : int [1:4461] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. ..$ spIdx                     : int [1:4461] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. .. .. ..$ smoothed                  : logi [1:4461] NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ seqNum                    : int [1:4461] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. .. .. ..$ acquisitionNum            : int [1:4461] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. .. .. ..$ msLevel                   : int [1:4461] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. ..$ polarity                  : int [1:4461] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. ..$ originalPeaksCount        : int [1:4461] 540 587 559 551 516 582 535 505 521 506 ...
  .. .. .. .. .. ..$ totIonCurrent             : num [1:4461] 83909668 86605992 86937473 90614298 87715841 ...
  .. .. .. .. .. ..$ retentionTime             : num [1:4461] 480 480 481 481 481 ...
  .. .. .. .. .. ..$ basePeakMZ                : num [1:4461] 405 405 405 405 405 ...
  .. .. .. .. .. ..$ basePeakIntensity         : num [1:4461] 24036730 22788144 20415082 19514576 17573978 ...
  .. .. .. .. .. ..$ collisionEnergy           : num [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ ionisationEnergy          : num [1:4461] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. ..$ lowMZ                     : num [1:4461] 251 251 251 250 252 ...
  .. .. .. .. .. ..$ highMZ                    : num [1:4461] 600 600 600 600 600 ...
  .. .. .. .. .. ..$ precursorScanNum          : int [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ precursorMZ               : num [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ precursorCharge           : int [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ precursorIntensity        : num [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ mergedScan                : int [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ mergedResultScanNum       : int [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ mergedResultStartScanNum  : int [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ mergedResultEndScanNum    : int [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ injectionTime             : num [1:4461] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. .. .. ..$ filterString              : chr [1:4461] NA NA NA NA ...
  .. .. .. .. .. ..$ spectrumId                : chr [1:4461] "scan=1" "scan=2" "scan=3" "scan=4" ...
  .. .. .. .. .. ..$ centroided                : logi [1:4461] TRUE TRUE TRUE TRUE TRUE TRUE ...
  .. .. .. .. .. ..$ ionMobilityDriftTime      : num [1:4461] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
  .. .. .. .. .. ..$ isolationWindowTargetMZ   : num [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ isolationWindowLowerOffset: num [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ isolationWindowUpperOffset: num [1:4461] NA NA NA NA NA NA NA NA NA NA ...
  .. .. .. .. .. ..$ scanWindowLowerLimit      : num [1:4461] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
  .. .. .. .. .. ..$ scanWindowUpperLimit      : num [1:4461] NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ...
  .. .. .. .. .. ..$ spectrum                  : int [1:4461] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. .. ..@ dimLabels        : chr [1:2] "rowNames" "columnNames"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. ..@ experimentData        :Formal class 'MIAPE' [package "MSnbase"] with 30 slots
  .. .. .. .. ..@ title                   : chr ""
  .. .. .. .. ..@ url                     : chr ""
  .. .. .. .. ..@ abstract                : chr ""
  .. .. .. .. ..@ pubMedIds               : chr ""
  .. .. .. .. ..@ samples                 : list()
  .. .. .. .. ..@ preprocessing           : list()
  .. .. .. .. ..@ other                   : list()
  .. .. .. .. ..@ dateStamp               : chr(0) 
  .. .. .. .. ..@ name                    : chr ""
  .. .. .. .. ..@ lab                     : chr ""
  .. .. .. .. ..@ contact                 : chr ""
  .. .. .. .. ..@ email                   : chr ""
  .. .. .. .. ..@ instrumentModel         : chr [1:10] "instrument model" "instrument model" "instrument model" "instrument model" ...
  .. .. .. .. ..@ instrumentManufacturer  : chr [1:10] "Unknown" "Unknown" "Unknown" "Unknown" ...
  .. .. .. .. ..@ instrumentCustomisations: chr(0) 
  .. .. .. .. ..@ softwareName            : chr(0) 
  .. .. .. .. ..@ softwareVersion         : chr(0) 
  .. .. .. .. ..@ switchingCriteria       : chr(0) 
  .. .. .. .. ..@ isolationWidth          : num(0) 
  .. .. .. .. ..@ parameterFile           : chr(0) 
  .. .. .. .. ..@ ionSource               : chr [1:10] "Unknown" "Unknown" "Unknown" "Unknown" ...
  .. .. .. .. ..@ ionSourceDetails        : chr(0) 
  .. .. .. .. ..@ analyser                : chr [1:10] "Unknown" "Unknown" "Unknown" "Unknown" ...
  .. .. .. .. ..@ analyserDetails         : chr(0) 
  .. .. .. .. ..@ collisionGas            : chr(0) 
  .. .. .. .. ..@ collisionPressure       : num(0) 
  .. .. .. .. ..@ collisionEnergy         : chr(0) 
  .. .. .. .. ..@ detectorType            : chr [1:10] "Unknown" "Unknown" "Unknown" "Unknown" ...
  .. .. .. .. ..@ detectorSensitivity     : chr(0) 
  .. .. .. .. ..@ .__classVersion__       :Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. .. .. ..$ : int [1:3] 0 2 2
  .. .. ..@ protocolData          :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. .. .. ..$ labelDescription: chr(0) 
  .. .. .. .. ..@ data             :'data.frame':	0 obs. of  0 variables
  .. .. .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. ..@ processingData        :Formal class 'MSnProcess' [package "MSnbase"] with 10 slots
  .. .. .. .. ..@ files            : chr [1:10] "/tmp/RtmpkZ19sH/SpectraData/Control/CD-6KUCT.mzML" "/tmp/RtmpkZ19sH/SpectraData/Control/CD-77FXR.mzML" "/tmp/RtmpkZ19sH/SpectraData/Control/CD-9OS5Y.mzML" "/tmp/RtmpkZ19sH/SpectraData/Control/CD-9WOBP.mzML" ...
  .. .. .. .. ..@ processing       : chr [1:2] "Data loaded [Thu Oct 15 16:54:13 2020]" "Filter: select MS level(s) 1 [Thu Oct 15 16:54:13 2020]"
  .. .. .. .. ..@ merged           : logi(0) 
  .. .. .. .. ..@ cleaned          : logi(0) 
  .. .. .. .. ..@ removedPeaks     : chr(0) 
  .. .. .. .. ..@ smoothed         : logi NA
  .. .. .. .. ..@ trimmed          : num(0) 
  .. .. .. .. ..@ normalised       : logi(0) 
  .. .. .. .. ..@ MSnbaseVersion   : chr "2.15.3"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 0 1 3
  .. .. ..@ .cache                :<environment: 0x55caafe6c9a8> 
  .. .. ..@ .__classVersion__     :Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 5
  .. .. .. .. .. ..$ : int [1:3] 4 0 0
  .. .. .. .. .. ..$ : int [1:3] 2 49 0
  .. .. .. .. .. ..$ : int [1:3] 0 1 1
  .. .. .. .. .. ..$ : int [1:3] 0 3 1
  .. .. .. .. .. ..$ : int [1:3] 0 0 1
  ..@ rawInMemory     :Formal class 'MSnExp' [package "MSnbase"] with 8 slots
  .. .. ..@ assayData        :<environment: 0x55caafe6d930> 
  .. .. ..@ phenoData        :Formal class 'NAnnotatedDataFrame' [package "MSnbase"] with 6 slots
  .. .. .. .. ..@ multiplex        : num 1
  .. .. .. .. ..@ multiLabels      : chr "Single run"
  .. .. .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. .. .. ..$ labelDescription: chr NA
  .. .. .. .. ..@ data             :'data.frame':	2 obs. of  1 variable:
  .. .. .. .. .. ..$ sample_name: chr [1:2] "QC_PREFA02.mzML" "QC_PREFB02.mzML"
  .. .. .. .. ..@ dimLabels        : chr [1:2] "rowNames" "columnNames"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 0 0 3
  .. .. ..@ featureData      :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. ..@ varMetadata      :'data.frame':	1 obs. of  1 variable:
  .. .. .. .. .. ..$ labelDescription: chr NA
  .. .. .. .. ..@ data             :'data.frame':	446 obs. of  1 variable:
  .. .. .. .. .. ..$ spectrum: num [1:446] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. .. ..@ dimLabels        : chr [1:2] "rowNames" "columnNames"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. ..@ experimentData   :Formal class 'MIAPE' [package "MSnbase"] with 30 slots
  .. .. .. .. ..@ title                   : chr ""
  .. .. .. .. ..@ url                     : chr ""
  .. .. .. .. ..@ abstract                : chr ""
  .. .. .. .. ..@ pubMedIds               : chr ""
  .. .. .. .. ..@ samples                 : list()
  .. .. .. .. ..@ preprocessing           : list()
  .. .. .. .. ..@ other                   : list()
  .. .. .. .. ..@ dateStamp               : chr(0) 
  .. .. .. .. ..@ name                    : chr ""
  .. .. .. .. ..@ lab                     : chr ""
  .. .. .. .. ..@ contact                 : chr ""
  .. .. .. .. ..@ email                   : chr ""
  .. .. .. .. ..@ instrumentModel         : chr [1:2] "instrument model" "instrument model"
  .. .. .. .. ..@ instrumentManufacturer  : chr [1:2] "Unknown" "Unknown"
  .. .. .. .. ..@ instrumentCustomisations: chr(0) 
  .. .. .. .. ..@ softwareName            : chr(0) 
  .. .. .. .. ..@ softwareVersion         : chr(0) 
  .. .. .. .. ..@ switchingCriteria       : chr(0) 
  .. .. .. .. ..@ isolationWidth          : num(0) 
  .. .. .. .. ..@ parameterFile           : chr(0) 
  .. .. .. .. ..@ ionSource               : chr [1:2] "Unknown" "Unknown"
  .. .. .. .. ..@ ionSourceDetails        : chr(0) 
  .. .. .. .. ..@ analyser                : chr [1:2] "Unknown" "Unknown"
  .. .. .. .. ..@ analyserDetails         : chr(0) 
  .. .. .. .. ..@ collisionGas            : chr(0) 
  .. .. .. .. ..@ collisionPressure       : num(0) 
  .. .. .. .. ..@ collisionEnergy         : chr(0) 
  .. .. .. .. ..@ detectorType            : chr [1:2] "Unknown" "Unknown"
  .. .. .. .. ..@ detectorSensitivity     : chr(0) 
  .. .. .. .. ..@ .__classVersion__       :Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 2
  .. .. .. .. .. .. .. ..$ : int [1:3] 1 0 0
  .. .. .. .. .. .. .. ..$ : int [1:3] 0 2 2
  .. .. ..@ protocolData     :Formal class 'AnnotatedDataFrame' [package "Biobase"] with 4 slots
  .. .. .. .. ..@ varMetadata      :'data.frame':	0 obs. of  1 variable:
  .. .. .. .. .. ..$ labelDescription: chr(0) 
  .. .. .. .. ..@ data             :'data.frame':	0 obs. of  0 variables
  .. .. .. .. ..@ dimLabels        : chr [1:2] "sampleNames" "sampleColumns"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 1 1 0
  .. .. ..@ processingData   :Formal class 'MSnProcess' [package "MSnbase"] with 10 slots
  .. .. .. .. ..@ files            : chr [1:2] "/tmp/RtmpkZ19sH/SpectraData/QC/QC_PREFA02.mzML" "/tmp/RtmpkZ19sH/SpectraData/QC/QC_PREFB02.mzML"
  .. .. .. .. ..@ processing       : chr "Data loaded: Thu Oct 15 16:50:41 2020"
  .. .. .. .. ..@ merged           : logi(0) 
  .. .. .. .. ..@ cleaned          : logi(0) 
  .. .. .. .. ..@ removedPeaks     : chr(0) 
  .. .. .. .. ..@ smoothed         : logi NA
  .. .. .. .. ..@ trimmed          : num(0) 
  .. .. .. .. ..@ normalised       : logi(0) 
  .. .. .. .. ..@ MSnbaseVersion   : chr "2.15.3"
  .. .. .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. .. .. ..@ .Data:List of 1
  .. .. .. .. .. .. .. ..$ : int [1:3] 0 1 3
  .. .. ..@ .cache           :<environment: 0x55caafc3e3d8> 
  .. .. ..@ .__classVersion__:Formal class 'Versions' [package "Biobase"] with 1 slot
  .. .. .. .. ..@ .Data:List of 4
  .. .. .. .. .. ..$ : int [1:3] 4 0 0
  .. .. .. .. .. ..$ : int [1:3] 2 49 0
  .. .. .. .. .. ..$ : int [1:3] 0 1 1
  .. .. .. .. .. ..$ : int [1:3] 0 3 1
  ..@ params          :List of 25
  .. ..$ Peak_method       : chr "centWave"
  .. ..$ RT_method         : chr "loess"
  .. ..$ ppm               : num 2.24
  .. ..$ peakwidth         : num [1:2] 10 15
  .. ..$ snthresh          : num 7.5
  .. ..$ prefilter         : num [1:2] 2 21862
  .. ..$ roiList           : list()
  .. ..$ firstBaselineCheck: logi TRUE
  .. ..$ roiScales         : num(0) 
  .. ..$ mzCenterFun       : chr "wMean"
  .. ..$ integrate         : num 1
  .. ..$ mzdiff            : num 0.012
  .. ..$ fitgauss          : logi FALSE
  .. ..$ noise             : num 12576
  .. ..$ verboseColumns    : logi FALSE
  .. ..$ binSize           : num 0.25
  .. ..$ bw                : num 5
  .. ..$ minFraction       : num 0.8
  .. ..$ minSamples        : num 1
  .. ..$ maxFeatures       : num 100
  .. ..$ extraPeaks        : num 1
  .. ..$ smooth            : chr "loess"
  .. ..$ span              : num 0.4
  .. ..$ family            : chr "gaussian"
  .. ..$ subsetAdjust      : chr "average"
  ..@ peakpicking     :List of 2
  .. ..$ chromPeaks   : num [1:4692, 1:11] 430 540 350 554 567 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4692] "CP0001" "CP0002" "CP0003" "CP0004" ...
  .. .. .. ..$ : chr [1:11] "mz" "mzmin" "mzmax" "rt" ...
  .. ..$ chromPeakData:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  .. .. .. ..@ rownames       : chr [1:4692] "CP0001" "CP0002" "CP0003" "CP0004" ...
  .. .. .. ..@ nrows          : int 4692
  .. .. .. ..@ listData       :List of 2
  .. .. .. .. ..$ ms_level : int [1:4692] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. ..$ is_filled: logi [1:4692] FALSE FALSE FALSE FALSE FALSE FALSE ...
  .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  ..@ peakgrouping    :List of 2
  .. ..$ :Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  .. .. .. ..@ rownames       : chr [1:431] "FT001" "FT002" "FT003" "FT004" ...
  .. .. .. ..@ nrows          : int 431
  .. .. .. ..@ listData       :List of 11
  .. .. .. .. ..$ mzmed  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ mzmin  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ mzmax  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ rtmed  : num [1:431] 503 498 498 511 576 ...
  .. .. .. .. ..$ rtmin  : num [1:431] 502 498 498 511 569 ...
  .. .. .. .. ..$ rtmax  : num [1:431] 509 498 498 523 578 ...
  .. .. .. .. ..$ npeaks : num [1:431] 6 10 10 9 5 12 2 5 10 8 ...
  .. .. .. .. ..$ Control: num [1:431] 1 4 4 2 3 4 0 1 3 2 ...
  .. .. .. .. ..$ Disease: num [1:431] 3 4 4 3 0 4 0 2 4 4 ...
  .. .. .. .. ..$ QC     : num [1:431] 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. ..$ peakidx:List of 431
  .. .. .. .. .. ..$ : int [1:6] 1091 2419 2856 3304 3688 4216
  .. .. .. .. .. ..$ : int [1:10] 108 594 1120 1582 1959 2428 2839 3305 3716 4224
  .. .. .. .. .. ..$ : int [1:10] 74 561 1083 1555 1935 2414 2803 3285 3684 4212
  .. .. .. .. .. ..$ : int [1:9] 135 1133 2031 2442 2487 3314 3351 3750 4276
  .. .. .. .. .. ..$ : int [1:5] 949 1411 1814 4039 4575
  .. .. .. .. .. ..$ : int [1:12] 155 635 636 1068 1167 1541 2014 2440 2910 3358 ...
  .. .. .. .. .. ..$ : int [1:2] 3726 4254
  .. .. .. .. .. ..$ : int [1:5] 183 2012 2891 3767 4299
  .. .. .. .. .. ..$ : int [1:10] 72 1123 1524 1642 1971 2416 2861 3318 3727 4242
  .. .. .. .. .. ..$ : int [1:8] 841 1343 2202 2600 3101 3483 3846 4462
  .. .. .. .. .. ..$ : int [1:8] 248 730 1250 2105 2506 2987 3401 3819
  .. .. .. .. .. ..$ : int [1:9] 69 576 625 1085 1545 1924 2426 3294 3694
  .. .. .. .. .. ..$ : int [1:4] 2311 3580 4079 4627
  .. .. .. .. .. ..$ : int [1:11] 306 842 1340 1753 2092 2192 2588 3108 3478 3952 ...
  .. .. .. .. .. ..$ : int [1:3] 2064 3692 4226
  .. .. .. .. .. ..$ : int [1:11] 289 822 1319 1738 2173 2573 3095 3458 3935 4416 ...
  .. .. .. .. .. ..$ : int [1:12] 169 647 1163 1605 2068 2174 2511 2900 3019 3374 ...
  .. .. .. .. .. ..$ : int [1:10] 826 1865 2175 2561 3020 3071 3467 3880 3885 4405
  .. .. .. .. .. ..$ : int [1:7] 396 414 1410 1819 2729 4059 4616
  .. .. .. .. .. ..$ : int [1:3] 1437 4096 4661
  .. .. .. .. .. ..$ : int [1:7] 513 1045 1900 2393 2790 3655 4181
  .. .. .. .. .. ..$ : int [1:10] 465 1004 1462 1866 2350 2735 3208 3594 4105 4670
  .. .. .. .. .. ..$ : int [1:7] 299 834 1317 1749 2579 3940 4476
  .. .. .. .. .. ..$ : int [1:7] 452 974 1857 2721 3588 4071 4635
  .. .. .. .. .. ..$ : int [1:9] 260 784 1282 1714 2125 2558 3433 3876 4414
  .. .. .. .. .. ..$ : int [1:2] 3851 4377
  .. .. .. .. .. ..$ : int [1:9] 262 785 1286 1716 2131 2562 3434 3889 4420
  .. .. .. .. .. ..$ : int [1:9] 100 1538 1955 2406 2464 3308 3334 3742 4269
  .. .. .. .. .. ..$ : int [1:5] 850 1345 3428 3865 4427
  .. .. .. .. .. ..$ : int [1:9] 256 776 1273 1708 2121 2555 3424 3881 4402
  .. .. .. .. .. ..$ : int [1:2] 3728 4248
  .. .. .. .. .. ..$ : int [1:12] 400 934 988 1467 1830 2280 2301 2699 3211 3563 ...
  .. .. .. .. .. ..$ : int [1:9] 291 820 1311 1739 2171 3076 3454 3927 4467
  .. .. .. .. .. ..$ : int [1:11] 114 598 599 1119 1960 2400 2411 2823 3267 3640 ...
  .. .. .. .. .. ..$ : int [1:7] 109 530 2379 2783 3242 3710 4149
  .. .. .. .. .. ..$ : int [1:10] 292 845 1328 1779 2207 2583 3081 3085 3447 4451
  .. .. .. .. .. ..$ : int [1:5] 236 1274 2544 3866 4393
  .. .. .. .. .. ..$ : int [1:13] 228 246 761 1258 2103 2126 2551 2994 3414 3857 ...
  .. .. .. .. .. ..$ : int [1:11] 61 497 568 1055 1079 1526 1918 2418 3250 3685 ...
  .. .. .. .. .. ..$ : int [1:8] 551 1071 1530 2372 2402 3362 3664 4183
  .. .. .. .. .. ..$ : int [1:3] 1132 3797 4319
  .. .. .. .. .. ..$ : int [1:6] 1262 2164 2575 3439 3908 4437
  .. .. .. .. .. ..$ : int [1:10] 531 1098 1933 2403 2804 2815 3264 3627 3654 4177
  .. .. .. .. .. ..$ : int [1:8] 111 523 1036 1539 1907 2384 3651 4175
  .. .. .. .. .. ..$ : int [1:10] 300 843 1292 2205 2580 3087 3406 3476 3955 4497
  .. .. .. .. .. ..$ : int [1:9] 273 804 1290 2139 2565 3040 3440 3898 4425
  .. .. .. .. .. ..$ : int [1:6] 347 1364 2283 3139 3991 4529
  .. .. .. .. .. ..$ : int [1:3] 3407 3770 4308
  .. .. .. .. .. ..$ : int [1:12] 349 835 1349 1790 2195 2197 2352 2648 3096 3474 ...
  .. .. .. .. .. ..$ : int [1:11] 301 830 1330 1754 2176 2260 2589 3066 3468 3936 ...
  .. .. .. .. .. ..$ : int [1:10] 639 666 1171 2018 2472 2474 2928 3348 3771 4302
  .. .. .. .. .. ..$ : int [1:5] 2257 2307 2651 4006 4554
  .. .. .. .. .. ..$ : int [1:3] 2006 3757 4293
  .. .. .. .. .. ..$ : int [1:11] 533 1111 1508 1919 1949 2410 2766 2848 3253 3635 ...
  .. .. .. .. .. ..$ : int [1:5] 261 757 3848 4385 4469
  .. .. .. .. .. ..$ : int [1:3] 3561 4054 4594
  .. .. .. .. .. ..$ : int [1:6] 392 847 1313 1734 3944 4495
  .. .. .. .. .. ..$ : int [1:12] 305 1322 2180 2530 2581 3065 3072 3455 3915 3941 ...
  .. .. .. .. .. ..$ : int [1:12] 88 584 1080 1551 1905 1944 2412 2767 2834 3306 ...
  .. .. .. .. .. ..$ : int [1:5] 297 836 1304 1736 3928
  .. .. .. .. .. ..$ : int [1:7] 610 1129 1974 1978 2443 3735 4262
  .. .. .. .. .. ..$ : int [1:2] 3641 4161
  .. .. .. .. .. ..$ : int [1:2] 3913 4445
  .. .. .. .. .. ..$ : int [1:9] 233 831 1314 2136 2564 3054 3436 3891 4424
  .. .. .. .. .. ..$ : int [1:5] 720 1895 3252 3639 4159
  .. .. .. .. .. ..$ : int [1:10] 220 682 781 2112 2496 2553 3022 3430 3872 4403
  .. .. .. .. .. ..$ : int [1:3] 928 4037 4587
  .. .. .. .. .. ..$ : int [1:2] 3909 4442
  .. .. .. .. .. ..$ : int [1:8] 416 973 1436 2705 3162 3182 4072 4621
  .. .. .. .. .. ..$ : int [1:10] 402 951 1423 2696 3170 3556 4040 4106 4604 4664
  .. .. .. .. .. ..$ : int [1:17] 62 81 560 563 1086 1553 1937 1963 2407 2413 ...
  .. .. .. .. .. ..$ : int [1:8] 268 808 2618 3043 3902 3957 4452 4498
  .. .. .. .. .. ..$ : int [1:12] 42 58 548 549 1534 1926 2395 2398 2802 3282 ...
  .. .. .. .. .. ..$ : int [1:5] 332 872 1367 1780 2258
  .. .. .. .. .. ..$ : int [1:4] 40 517 3636 4160
  .. .. .. .. .. ..$ : int [1:8] 293 762 1288 2669 3052 3524 3948 4431
  .. .. .. .. .. ..$ : int [1:7] 112 608 1168 2286 3051 3779 4305
  .. .. .. .. .. ..$ : int [1:3] 1938 3701 4227
  .. .. .. .. .. ..$ : int [1:14] 468 758 797 1241 1263 2093 2140 2538 2584 3425 ...
  .. .. .. .. .. ..$ : int [1:4] 470 1872 3972 4645
  .. .. .. .. .. ..$ : int [1:7] 237 731 809 2156 2574 3027 3432
  .. .. .. .. .. ..$ : int [1:8] 460 1008 1483 2349 2734 3592 4114 4671
  .. .. .. .. .. ..$ : int [1:5] 143 2884 3315 3737 4257
  .. .. .. .. .. ..$ : int [1:4] 1892 2374 2750 3237
  .. .. .. .. .. ..$ : int [1:6] 914 2276 2682 3533 4025 4576
  .. .. .. .. .. ..$ : int [1:6] 1002 2346 2730 3589 4107 4641
  .. .. .. .. .. ..$ : int [1:5] 2203 3082 3475 3949 4492
  .. .. .. .. .. ..$ : int [1:5] 2256 2647 2666 4014 4562
  .. .. .. .. .. ..$ : int [1:7] 566 1034 1105 2420 2791 3618 4143
  .. .. .. .. .. ..$ : int [1:5] 2193 3057 3469 3917 4465
  .. .. .. .. .. ..$ : int [1:14] 279 296 814 1308 1735 2157 2159 2188 2236 3067 ...
  .. .. .. .. .. ..$ : int [1:11] 125 129 170 630 667 1169 1595 2466 2935 3776 ...
  .. .. .. .. .. ..$ : int [1:7] 391 411 952 1425 2690 4060 4610
  .. .. .. .. .. ..$ : int [1:2] 4035 4589
  .. .. .. .. .. ..$ : int [1:10] 102 141 159 631 1159 1635 2448 2911 3762 4295
  .. .. .. .. .. ..$ : int [1:7] 2065 2512 2885 3785 3882 4320 4407
  .. .. .. .. .. ..$ : int [1:2] 3847 4382
  .. .. .. .. .. ..$ : int [1:7] 759 2083 2513 2539 2998 3411 4389
  .. .. .. .. .. ..$ : int [1:4] 359 2329 4061 4617
  .. .. .. .. .. .. [list output truncated]
  .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  .. ..$ :Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  .. .. .. ..@ rownames       : chr [1:431] "FT001" "FT002" "FT003" "FT004" ...
  .. .. .. ..@ nrows          : int 431
  .. .. .. ..@ listData       :List of 11
  .. .. .. .. ..$ mzmed  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ mzmin  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ mzmax  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ rtmed  : num [1:431] 503 498 498 511 577 ...
  .. .. .. .. ..$ rtmin  : num [1:431] 502 497 497 511 570 ...
  .. .. .. .. ..$ rtmax  : num [1:431] 509 499 499 523 578 ...
  .. .. .. .. ..$ npeaks : num [1:431] 6 10 10 9 5 12 2 5 10 8 ...
  .. .. .. .. ..$ Control: num [1:431] 1 4 4 2 3 4 0 1 3 2 ...
  .. .. .. .. ..$ Disease: num [1:431] 3 4 4 3 0 4 0 2 4 4 ...
  .. .. .. .. ..$ QC     : num [1:431] 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. ..$ peakidx:List of 431
  .. .. .. .. .. ..$ : num [1:6] 1091 2419 2856 3304 3688 ...
  .. .. .. .. .. ..$ : num [1:10] 108 594 1120 1582 1959 ...
  .. .. .. .. .. ..$ : num [1:10] 74 561 1083 1555 1935 ...
  .. .. .. .. .. ..$ : num [1:9] 135 1133 2031 2442 2487 ...
  .. .. .. .. .. ..$ : num [1:5] 949 1411 1814 4039 4575
  .. .. .. .. .. ..$ : num [1:12] 155 635 636 1068 1167 ...
  .. .. .. .. .. ..$ : num [1:2] 3726 4254
  .. .. .. .. .. ..$ : num [1:5] 183 2012 2891 3767 4299
  .. .. .. .. .. ..$ : num [1:10] 72 1123 1524 1642 1971 ...
  .. .. .. .. .. ..$ : num [1:8] 841 1343 2202 2600 3101 ...
  .. .. .. .. .. ..$ : num [1:8] 248 730 1250 2105 2506 ...
  .. .. .. .. .. ..$ : num [1:9] 69 576 625 1085 1545 ...
  .. .. .. .. .. ..$ : num [1:4] 2311 3580 4079 4627
  .. .. .. .. .. ..$ : num [1:11] 306 842 1340 1753 2092 ...
  .. .. .. .. .. ..$ : num [1:3] 2064 3692 4226
  .. .. .. .. .. ..$ : num [1:11] 289 822 1319 1738 2173 ...
  .. .. .. .. .. ..$ : num [1:12] 169 647 1163 1605 2068 ...
  .. .. .. .. .. ..$ : num [1:10] 826 1865 2175 2561 3020 ...
  .. .. .. .. .. ..$ : num [1:7] 396 414 1410 1819 2729 ...
  .. .. .. .. .. ..$ : num [1:3] 1437 4096 4661
  .. .. .. .. .. ..$ : num [1:7] 513 1045 1900 2393 2790 ...
  .. .. .. .. .. ..$ : num [1:10] 465 1004 1462 1866 2350 ...
  .. .. .. .. .. ..$ : num [1:7] 299 834 1317 1749 2579 ...
  .. .. .. .. .. ..$ : num [1:7] 452 974 1857 2721 3588 ...
  .. .. .. .. .. ..$ : num [1:9] 260 784 1282 1714 2125 ...
  .. .. .. .. .. ..$ : num [1:2] 3851 4377
  .. .. .. .. .. ..$ : num [1:9] 262 785 1286 1716 2131 ...
  .. .. .. .. .. ..$ : num [1:8] 100 1538 1955 2406 3308 ...
  .. .. .. .. .. ..$ : num [1:6] 850 1345 2464 3428 3865 ...
  .. .. .. .. .. ..$ : num [1:9] 256 776 1273 1708 2121 ...
  .. .. .. .. .. ..$ : num [1:2] 3728 4248
  .. .. .. .. .. ..$ : num [1:12] 400 934 988 1467 1830 ...
  .. .. .. .. .. ..$ : num [1:9] 291 820 1311 1739 2171 ...
  .. .. .. .. .. ..$ : num [1:11] 114 598 599 1119 1960 ...
  .. .. .. .. .. ..$ : num [1:7] 109 530 2379 2783 3242 ...
  .. .. .. .. .. ..$ : num [1:10] 292 845 1328 1779 2207 ...
  .. .. .. .. .. ..$ : num [1:5] 236 1274 2544 3866 4393
  .. .. .. .. .. ..$ : num [1:13] 228 246 761 1258 2103 ...
  .. .. .. .. .. ..$ : num [1:11] 61 497 568 1055 1079 ...
  .. .. .. .. .. ..$ : num [1:8] 551 1071 1530 2372 2402 ...
  .. .. .. .. .. ..$ : num [1:3] 1132 3797 4319
  .. .. .. .. .. ..$ : num [1:6] 1262 2164 2575 3439 3908 ...
  .. .. .. .. .. ..$ : num [1:10] 531 1098 1933 2403 2804 ...
  .. .. .. .. .. ..$ : num [1:8] 111 523 1036 1539 1907 ...
  .. .. .. .. .. ..$ : num [1:10] 300 843 1292 2205 2580 ...
  .. .. .. .. .. ..$ : num [1:9] 273 804 1290 2139 2565 ...
  .. .. .. .. .. ..$ : num [1:6] 347 1364 2283 3139 3991 ...
  .. .. .. .. .. ..$ : num [1:3] 3407 3770 4308
  .. .. .. .. .. ..$ : num [1:12] 349 835 1349 1790 2195 ...
  .. .. .. .. .. ..$ : num [1:11] 301 830 1330 1754 2176 ...
  .. .. .. .. .. ..$ : num [1:10] 639 666 1171 2018 2472 ...
  .. .. .. .. .. ..$ : num [1:5] 2257 2307 2651 4006 4554
  .. .. .. .. .. ..$ : num [1:3] 2006 3757 4293
  .. .. .. .. .. ..$ : num [1:11] 533 1111 1508 1919 1949 ...
  .. .. .. .. .. ..$ : num [1:4] 261 757 3848 4385
  .. .. .. .. .. ..$ : num [1:3] 3561 4054 4594
  .. .. .. .. .. ..$ : num [1:6] 392 847 1313 1734 3944 ...
  .. .. .. .. .. ..$ : num [1:12] 305 1322 2180 2530 2581 ...
  .. .. .. .. .. ..$ : num [1:12] 88 584 1080 1551 1905 ...
  .. .. .. .. .. ..$ : num [1:5] 297 836 1304 1736 3928
  .. .. .. .. .. ..$ : num [1:7] 610 1129 1974 1978 2443 ...
  .. .. .. .. .. ..$ : num [1:2] 3913 4445
  .. .. .. .. .. ..$ : num [1:2] 3641 4161
  .. .. .. .. .. ..$ : num [1:9] 233 831 1314 2136 2564 ...
  .. .. .. .. .. ..$ : num [1:5] 720 1895 3252 3639 4159
  .. .. .. .. .. ..$ : num [1:10] 220 682 781 2112 2496 ...
  .. .. .. .. .. ..$ : num [1:3] 928 4037 4587
  .. .. .. .. .. ..$ : num [1:2] 3909 4442
  .. .. .. .. .. ..$ : num [1:8] 416 973 1436 2705 3162 ...
  .. .. .. .. .. ..$ : num [1:10] 402 951 1423 2696 3170 ...
  .. .. .. .. .. ..$ : num [1:17] 62 81 560 563 1086 ...
  .. .. .. .. .. ..$ : num [1:8] 268 808 2618 3043 3902 ...
  .. .. .. .. .. ..$ : num [1:12] 42 58 548 549 1534 ...
  .. .. .. .. .. ..$ : num [1:5] 332 872 1367 1780 2258
  .. .. .. .. .. ..$ : num [1:4] 40 517 3636 4160
  .. .. .. .. .. ..$ : num [1:8] 293 762 1288 2669 3052 ...
  .. .. .. .. .. ..$ : num [1:7] 112 608 1168 2286 3051 ...
  .. .. .. .. .. ..$ : num [1:3] 1938 3701 4227
  .. .. .. .. .. ..$ : num [1:14] 468 758 797 1241 1263 ...
  .. .. .. .. .. ..$ : num [1:4] 470 1872 3972 4645
  .. .. .. .. .. ..$ : num [1:7] 237 731 809 2156 2574 ...
  .. .. .. .. .. ..$ : num [1:8] 460 1008 1483 2349 2734 ...
  .. .. .. .. .. ..$ : num [1:5] 143 2884 3315 3737 4257
  .. .. .. .. .. ..$ : num [1:4] 1892 2374 2750 3237
  .. .. .. .. .. ..$ : num [1:6] 1002 2346 2730 3589 4107 ...
  .. .. .. .. .. ..$ : num [1:6] 914 2276 2682 3533 4025 ...
  .. .. .. .. .. ..$ : num [1:5] 2203 3082 3475 3949 4492
  .. .. .. .. .. ..$ : num [1:5] 2256 2647 2666 4014 4562
  .. .. .. .. .. ..$ : num [1:7] 566 1034 1105 2420 2791 ...
  .. .. .. .. .. ..$ : num [1:5] 2193 3057 3469 3917 4465
  .. .. .. .. .. ..$ : num [1:14] 279 296 814 1308 1735 ...
  .. .. .. .. .. ..$ : num [1:11] 125 129 170 630 667 ...
  .. .. .. .. .. ..$ : num [1:7] 391 411 952 1425 2690 ...
  .. .. .. .. .. ..$ : num [1:2] 4035 4589
  .. .. .. .. .. ..$ : num [1:10] 102 141 159 631 1159 ...
  .. .. .. .. .. ..$ : num [1:4] 2065 2512 3882 4407
  .. .. .. .. .. ..$ : num [1:3] 2885 3785 4320
  .. .. .. .. .. ..$ : num [1:2] 3847 4382
  .. .. .. .. .. ..$ : num [1:7] 759 2083 2513 2539 2998 ...
  .. .. .. .. .. .. [list output truncated]
  .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  ..@ peakRTcorrection:List of 3
  .. ..$ pkGrpMat_Raw: num [1:74, 1:10] NA 491 507 495 496 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:74] "FT104" "FT128" "FT044" "FT187" ...
  .. .. .. ..$ : chr [1:10] "CD-6KUCT" "CD-77FXR" "CD-9OS5Y" "CD-9WOBP" ...
  .. ..$ chromPeaks  : num [1:4692, 1:11] 430 540 350 554 567 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:4692] "CP0001" "CP0002" "CP0003" "CP0004" ...
  .. .. .. ..$ : chr [1:11] "mz" "mzmin" "mzmax" "rt" ...
  .. ..$ adjustedRT  :List of 10
  .. .. ..$ 1 : Named num [1:446] 481 481 481 481 482 ...
  .. .. .. ..- attr(*, "names")= chr [1:446] "F01.S001" "F01.S002" "F01.S003" "F01.S004" ...
  .. .. ..$ 2 : Named num [1:446] 481 482 482 482 482 ...
  .. .. .. ..- attr(*, "names")= chr [1:446] "F02.S001" "F02.S002" "F02.S003" "F02.S004" ...
  .. .. ..$ 3 : Named num [1:446] 481 481 481 482 482 ...
  .. .. .. ..- attr(*, "names")= chr [1:446] "F03.S001" "F03.S002" "F03.S003" "F03.S004" ...
  .. .. ..$ 4 : Named num [1:444] 480 480 481 481 481 ...
  .. .. .. ..- attr(*, "names")= chr [1:444] "F04.S001" "F04.S002" "F04.S003" "F04.S004" ...
  .. .. ..$ 5 : Named num [1:447] 481 481 482 482 482 ...
  .. .. .. ..- attr(*, "names")= chr [1:447] "F05.S001" "F05.S002" "F05.S003" "F05.S004" ...
  .. .. ..$ 6 : Named num [1:447] 480 480 480 481 481 ...
  .. .. .. ..- attr(*, "names")= chr [1:447] "F06.S001" "F06.S002" "F06.S003" "F06.S004" ...
  .. .. ..$ 7 : Named num [1:446] 481 481 482 482 482 ...
  .. .. .. ..- attr(*, "names")= chr [1:446] "F07.S001" "F07.S002" "F07.S003" "F07.S004" ...
  .. .. ..$ 8 : Named num [1:446] 478 479 479 479 480 ...
  .. .. .. ..- attr(*, "names")= chr [1:446] "F08.S001" "F08.S002" "F08.S003" "F08.S004" ...
  .. .. ..$ 9 : Named num [1:447] 480 480 481 481 481 ...
  .. .. .. ..- attr(*, "names")= chr [1:447] "F09.S001" "F09.S002" "F09.S003" "F09.S004" ...
  .. .. ..$ 10: Named num [1:446] 480 480 481 481 481 ...
  .. .. .. ..- attr(*, "names")= chr [1:446] "F10.S001" "F10.S002" "F10.S003" "F10.S004" ...
  ..@ peakfilling     :List of 2
  .. ..$ msFeatureData    :List of 4
  .. .. ..$ pkGrpMat_Raw : num [1:74, 1:10] NA 491 507 495 496 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:74] "FT104" "FT128" "FT044" "FT187" ...
  .. .. .. .. ..$ : chr [1:10] "CD-6KUCT" "CD-77FXR" "CD-9OS5Y" "CD-9WOBP" ...
  .. .. ..$ chromPeaks   : num [1:5711, 1:11] 430 540 350 554 567 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:5711] "CP0001" "CP0002" "CP0003" "CP0004" ...
  .. .. .. .. ..$ : chr [1:11] "mz" "mzmin" "mzmax" "rt" ...
  .. .. ..$ adjustedRT   :List of 10
  .. .. .. ..$ 1 : Named num [1:446] 481 481 481 481 482 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:446] "F01.S001" "F01.S002" "F01.S003" "F01.S004" ...
  .. .. .. ..$ 2 : Named num [1:446] 481 482 482 482 482 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:446] "F02.S001" "F02.S002" "F02.S003" "F02.S004" ...
  .. .. .. ..$ 3 : Named num [1:446] 481 481 481 482 482 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:446] "F03.S001" "F03.S002" "F03.S003" "F03.S004" ...
  .. .. .. ..$ 4 : Named num [1:444] 480 480 481 481 481 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:444] "F04.S001" "F04.S002" "F04.S003" "F04.S004" ...
  .. .. .. ..$ 5 : Named num [1:447] 481 481 482 482 482 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:447] "F05.S001" "F05.S002" "F05.S003" "F05.S004" ...
  .. .. .. ..$ 6 : Named num [1:447] 480 480 480 481 481 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:447] "F06.S001" "F06.S002" "F06.S003" "F06.S004" ...
  .. .. .. ..$ 7 : Named num [1:446] 481 481 482 482 482 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:446] "F07.S001" "F07.S002" "F07.S003" "F07.S004" ...
  .. .. .. ..$ 8 : Named num [1:446] 478 479 479 479 480 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:446] "F08.S001" "F08.S002" "F08.S003" "F08.S004" ...
  .. .. .. ..$ 9 : Named num [1:447] 480 480 481 481 481 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:447] "F09.S001" "F09.S002" "F09.S003" "F09.S004" ...
  .. .. .. ..$ 10: Named num [1:446] 480 480 481 481 481 ...
  .. .. .. .. ..- attr(*, "names")= chr [1:446] "F10.S001" "F10.S002" "F10.S003" "F10.S004" ...
  .. .. ..$ chromPeakData:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  .. .. .. .. ..@ rownames       : chr [1:5711] "CP0001" "CP0002" "CP0003" "CP0004" ...
  .. .. .. .. ..@ nrows          : int 5711
  .. .. .. .. ..@ listData       :List of 2
  .. .. .. .. .. ..$ ms_level : int [1:5711] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. .. .. ..$ is_filled: logi [1:5711] FALSE FALSE FALSE FALSE FALSE FALSE ...
  .. .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. .. ..@ elementMetadata: NULL
  .. .. .. .. ..@ metadata       : list()
  .. ..$ FeatureGroupTable:Formal class 'DFrame' [package "S4Vectors"] with 6 slots
  .. .. .. ..@ rownames       : chr [1:431] "FT001" "FT002" "FT003" "FT004" ...
  .. .. .. ..@ nrows          : int 431
  .. .. .. ..@ listData       :List of 11
  .. .. .. .. ..$ mzmed  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ mzmin  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ mzmax  : num [1:431] 253 253 254 255 255 ...
  .. .. .. .. ..$ rtmed  : num [1:431] 503 498 498 511 577 ...
  .. .. .. .. ..$ rtmin  : num [1:431] 502 497 497 511 570 ...
  .. .. .. .. ..$ rtmax  : num [1:431] 509 499 499 523 578 ...
  .. .. .. .. ..$ npeaks : num [1:431] 6 10 10 9 5 12 2 5 10 8 ...
  .. .. .. .. ..$ Control: num [1:431] 1 4 4 2 3 4 0 1 3 2 ...
  .. .. .. .. ..$ Disease: num [1:431] 3 4 4 3 0 4 0 2 4 4 ...
  .. .. .. .. ..$ QC     : num [1:431] 2 2 2 2 2 2 2 2 2 2 ...
  .. .. .. .. ..$ peakidx:List of 431
  .. .. .. .. .. ..$ : num [1:8] 1091 2419 2856 3304 3688 ...
  .. .. .. .. .. ..$ : num [1:10] 108 594 1120 1582 1959 ...
  .. .. .. .. .. ..$ : num [1:10] 74 561 1083 1555 1935 ...
  .. .. .. .. .. ..$ : num [1:11] 135 1133 2031 2442 2487 ...
  .. .. .. .. .. ..$ : num [1:7] 949 1411 1814 4039 4575 ...
  .. .. .. .. .. ..$ : num [1:12] 155 635 636 1068 1167 ...
  .. .. .. .. .. ..$ : num [1:8] 3726 4254 4695 4817 4911 ...
  .. .. .. .. .. ..$ : num [1:10] 183 2012 2891 3767 4299 ...
  .. .. .. .. .. ..$ : num [1:11] 72 1123 1524 1642 1971 ...
  .. .. .. .. .. ..$ : num [1:10] 841 1343 2202 2600 3101 ...
  .. .. .. .. .. ..$ : num [1:9] 248 730 1250 2105 2506 ...
  .. .. .. .. .. ..$ : num [1:11] 69 576 625 1085 1545 ...
  .. .. .. .. .. ..$ : num [1:10] 2311 3580 4079 4627 4697 ...
  .. .. .. .. .. ..$ : num [1:11] 306 842 1340 1753 2092 ...
  .. .. .. .. .. ..$ : num [1:8] 2064 3692 4226 4698 4821 ...
  .. .. .. .. .. ..$ : num [1:11] 289 822 1319 1738 2173 ...
  .. .. .. .. .. ..$ : num [1:12] 169 647 1163 1605 2068 ...
  .. .. .. .. .. ..$ : num [1:12] 826 1865 2175 2561 3020 ...
  .. .. .. .. .. ..$ : num [1:11] 396 414 1410 1819 2729 ...
  .. .. .. .. .. ..$ : num [1:3] 1437 4096 4661
  .. .. .. .. .. ..$ : num [1:10] 513 1045 1900 2393 2790 ...
  .. .. .. .. .. ..$ : num [1:10] 465 1004 1462 1866 2350 ...
  .. .. .. .. .. ..$ : num [1:10] 299 834 1317 1749 2579 ...
  .. .. .. .. .. ..$ : num [1:10] 452 974 1857 2721 3588 ...
  .. .. .. .. .. ..$ : num [1:10] 260 784 1282 1714 2125 ...
  .. .. .. .. .. ..$ : num [1:3] 3851 4377 4823
  .. .. .. .. .. ..$ : num [1:9] 262 785 1286 1716 2131 ...
  .. .. .. .. .. ..$ : num [1:11] 100 1538 1955 2406 3308 ...
  .. .. .. .. .. ..$ : num [1:10] 850 1345 2464 3428 3865 ...
  .. .. .. .. .. ..$ : num [1:9] 256 776 1273 1708 2121 ...
  .. .. .. .. .. ..$ : num [1:6] 3728 4248 4702 5157 5423 ...
  .. .. .. .. .. ..$ : num [1:13] 400 934 988 1467 1830 ...
  .. .. .. .. .. ..$ : num [1:10] 291 820 1311 1739 2171 ...
  .. .. .. .. .. ..$ : num [1:12] 114 598 599 1119 1960 ...
  .. .. .. .. .. ..$ : num [1:8] 109 530 2379 2783 3242 ...
  .. .. .. .. .. ..$ : num [1:11] 292 845 1328 1779 2207 ...
  .. .. .. .. .. ..$ : num [1:10] 236 1274 2544 3866 4393 ...
  .. .. .. .. .. ..$ : num [1:13] 228 246 761 1258 2103 ...
  .. .. .. .. .. ..$ : num [1:12] 61 497 568 1055 1079 ...
  .. .. .. .. .. ..$ : num [1:11] 551 1071 1530 2372 2402 ...
  .. .. .. .. .. ..$ : num [1:10] 1132 3797 4319 4704 4826 ...
  .. .. .. .. .. ..$ : num [1:7] 1262 2164 2575 3439 3908 ...
  .. .. .. .. .. ..$ : num [1:12] 531 1098 1933 2403 2804 ...
  .. .. .. .. .. ..$ : num [1:9] 111 523 1036 1539 1907 ...
  .. .. .. .. .. ..$ : num [1:11] 300 843 1292 2205 2580 ...
  .. .. .. .. .. ..$ : num [1:9] 273 804 1290 2139 2565 ...
  .. .. .. .. .. ..$ : num [1:8] 347 1364 2283 3139 3991 ...
  .. .. .. .. .. ..$ : num [1:9] 3407 3770 4308 4706 4828 ...
  .. .. .. .. .. ..$ : num [1:12] 349 835 1349 1790 2195 ...
  .. .. .. .. .. ..$ : num [1:11] 301 830 1330 1754 2176 ...
  .. .. .. .. .. ..$ : num [1:12] 639 666 1171 2018 2472 ...
  .. .. .. .. .. ..$ : num [1:11] 2257 2307 2651 4006 4554 ...
  .. .. .. .. .. ..$ : num [1:9] 2006 3757 4293 4709 4830 ...
  .. .. .. .. .. ..$ : num [1:12] 533 1111 1508 1919 1949 ...
  .. .. .. .. .. ..$ : num [1:10] 261 757 3848 4385 4921 ...
  .. .. .. .. .. ..$ : num [1:6] 3561 4054 4594 4831 4922 ...
  .. .. .. .. .. ..$ : num [1:7] 392 847 1313 1734 3944 ...
  .. .. .. .. .. ..$ : num [1:12] 305 1322 2180 2530 2581 ...
  .. .. .. .. .. ..$ : num [1:12] 88 584 1080 1551 1905 ...
  .. .. .. .. .. ..$ : num [1:6] 297 836 1304 1736 3928 ...
  .. .. .. .. .. ..$ : num [1:8] 610 1129 1974 1978 2443 ...
  .. .. .. .. .. ..$ : num [1:3] 3913 4445 5279
  .. .. .. .. .. ..$ : num [1:2] 3641 4161
  .. .. .. .. .. ..$ : num [1:9] 233 831 1314 2136 2564 ...
  .. .. .. .. .. ..$ : num [1:8] 720 1895 3252 3639 4159 ...
  .. .. .. .. .. ..$ : num [1:10] 220 682 781 2112 2496 ...
  .. .. .. .. .. ..$ : num [1:7] 928 4037 4587 4712 4924 ...
  .. .. .. .. .. ..$ : num [1:2] 3909 4442
  .. .. .. .. .. ..$ : num [1:11] 416 973 1436 2705 3162 ...
  .. .. .. .. .. ..$ : num [1:11] 402 951 1423 2696 3170 ...
  .. .. .. .. .. ..$ : num [1:17] 62 81 560 563 1086 ...
  .. .. .. .. .. ..$ : num [1:12] 268 808 2618 3043 3902 ...
  .. .. .. .. .. ..$ : num [1:12] 42 58 548 549 1534 ...
  .. .. .. .. .. ..$ : num [1:8] 332 872 1367 1780 2258 ...
  .. .. .. .. .. ..$ : num [1:4] 40 517 3636 4160
  .. .. .. .. .. ..$ : num [1:10] 293 762 1288 2669 3052 ...
  .. .. .. .. .. ..$ : num [1:10] 112 608 1168 2286 3051 ...
  .. .. .. .. .. ..$ : num [1:3] 1938 3701 4227
  .. .. .. .. .. ..$ : num [1:14] 468 758 797 1241 1263 ...
  .. .. .. .. .. ..$ : num [1:9] 470 1872 3972 4645 4832 ...
  .. .. .. .. .. ..$ : num [1:10] 237 731 809 2156 2574 ...
  .. .. .. .. .. ..$ : num [1:10] 460 1008 1483 2349 2734 ...
  .. .. .. .. .. ..$ : num [1:10] 143 2884 3315 3737 4257 ...
  .. .. .. .. .. ..$ : num [1:4] 1892 2374 2750 3237
  .. .. .. .. .. ..$ : num [1:9] 1002 2346 2730 3589 4107 ...
  .. .. .. .. .. ..$ : num [1:7] 914 2276 2682 3533 4025 ...
  .. .. .. .. .. ..$ : num [1:7] 2203 3082 3475 3949 4492 ...
  .. .. .. .. .. ..$ : num [1:7] 2256 2647 2666 4014 4562 ...
  .. .. .. .. .. ..$ : num [1:11] 566 1034 1105 2420 2791 ...
  .. .. .. .. .. ..$ : num [1:6] 2193 3057 3469 3917 4465 ...
  .. .. .. .. .. ..$ : num [1:14] 279 296 814 1308 1735 ...
  .. .. .. .. .. ..$ : num [1:13] 125 129 170 630 667 ...
  .. .. .. .. .. ..$ : num [1:10] 391 411 952 1425 2690 ...
  .. .. .. .. .. ..$ : num [1:8] 4035 4589 4716 4836 4930 ...
  .. .. .. .. .. ..$ : num [1:12] 102 141 159 631 1159 ...
  .. .. .. .. .. ..$ : num [1:5] 2065 2512 3882 4407 4837
  .. .. .. .. .. ..$ : num [1:7] 2885 3785 4320 4717 4838 ...
  .. .. .. .. .. ..$ : num [1:5] 3847 4382 4932 5176 5290
  .. .. .. .. .. ..$ : num [1:9] 759 2083 2513 2539 2998 ...
  .. .. .. .. .. .. [list output truncated]
  .. .. .. ..@ elementType    : chr "ANY"
  .. .. .. ..@ elementMetadata: NULL
  .. .. .. ..@ metadata       : list()
  ..@ peakAnnotation  :List of 6
  .. ..$ AnnotateObject     :List of 12
  .. .. ..$ sample        : num NA
  .. .. ..$ groupInfo     : num [1:431, 1:20] 253 253 254 255 255 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:20] "mz" "mzmin" "mzmax" "rt" ...
  .. .. ..$ polarity      : chr "negative"
  .. .. ..$ runParallel   :List of 1
  .. .. .. ..$ enable: num 0
  .. .. ..$ annoID        : num [1:67, 1:4] 173 223 245 241 236 394 389 234 262 312 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:4] "id" "grpID" "ruleID" "parentID"
  .. .. ..$ annoGrp       : num [1:32, 1:4] 1 2 3 4 5 6 7 8 9 10 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:4] "id" "mass" "ips" "psgrp"
  .. .. ..$ isoID         : num [1:76, 1:4] 2 4 6 10 14 22 28 36 45 49 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : chr [1:4] "mpeak" "isopeak" "iso" "charge"
  .. .. ..$ pspectra      :List of 274
  .. .. .. ..$ : int [1:2] 36 37
  .. .. .. ..$ : int [1:14] 173 174 175 223 224 225 226 245 250 254 ...
  .. .. .. ..$ : int [1:4] 77 117 284 287
  .. .. .. ..$ : int [1:4] 236 239 241 410
  .. .. .. ..$ : int [1:3] 48 211 399
  .. .. .. ..$ : int [1:6] 10 11 386 387 389 394
  .. .. .. ..$ : int [1:5] 229 234 262 264 288
  .. .. .. ..$ : int [1:3] 258 260 350
  .. .. .. ..$ : int [1:3] 243 247 248
  .. .. .. ..$ : int [1:5] 51 107 169 345 405
  .. .. .. ..$ : int [1:3] 2 127 312
  .. .. .. ..$ : int [1:3] 111 129 132
  .. .. .. ..$ : int [1:2] 341 342
  .. .. .. ..$ : int [1:2] 170 242
  .. .. .. ..$ : int [1:4] 109 110 272 276
  .. .. .. ..$ : int [1:3] 80 141 316
  .. .. .. ..$ : int [1:5] 14 16 33 91 116
  .. .. .. ..$ : int [1:7] 119 120 122 145 269 303 305
  .. .. .. ..$ : int [1:4] 275 314 371 414
  .. .. .. ..$ : int [1:2] 156 157
  .. .. .. ..$ : int [1:2] 270 274
  .. .. .. ..$ : int [1:4] 47 146 401 404
  .. .. .. ..$ : int 136
  .. .. .. ..$ : int [1:7] 62 68 128 131 161 165 210
  .. .. .. ..$ : int [1:8] 67 151 152 206 207 302 354 355
  .. .. .. ..$ : int [1:6] 13 100 296 297 300 304
  .. .. .. ..$ : int [1:4] 6 9 246 377
  .. .. .. ..$ : int 325
  .. .. .. ..$ : int [1:3] 56 244 344
  .. .. .. ..$ : int [1:3] 64 66 99
  .. .. .. ..$ : int [1:6] 82 85 278 279 408 409
  .. .. .. ..$ : int 324
  .. .. .. ..$ : int [1:3] 25 27 30
  .. .. .. ..$ : int [1:2] 397 398
  .. .. .. ..$ : int 61
  .. .. .. ..$ : int [1:2] 28 31
  .. .. .. ..$ : int [1:2] 313 315
  .. .. .. ..$ : int 381
  .. .. .. ..$ : int 144
  .. .. .. ..$ : int [1:2] 121 158
  .. .. .. ..$ : int 26
  .. .. .. ..$ : int 38
  .. .. .. ..$ : int 55
  .. .. .. ..$ : int 76
  .. .. .. ..$ : int 79
  .. .. .. ..$ : int 126
  .. .. .. ..$ : int 149
  .. .. .. ..$ : int 150
  .. .. .. ..$ : int 186
  .. .. .. ..$ : int 227
  .. .. .. ..$ : int 233
  .. .. .. ..$ : int 307
  .. .. .. ..$ : int [1:2] 334 335
  .. .. .. ..$ : int 178
  .. .. .. ..$ : int 290
  .. .. .. ..$ : int 295
  .. .. .. ..$ : int 327
  .. .. .. ..$ : int 425
  .. .. .. ..$ : int 426
  .. .. .. ..$ : int [1:2] 171 281
  .. .. .. ..$ : int 114
  .. .. .. ..$ : int 125
  .. .. .. ..$ : int 18
  .. .. .. ..$ : int 29
  .. .. .. ..$ : int 96
  .. .. .. ..$ : int 98
  .. .. .. ..$ : int 143
  .. .. .. ..$ : int 181
  .. .. .. ..$ : int 197
  .. .. .. ..$ : int 218
  .. .. .. ..$ : int 219
  .. .. .. ..$ : int 221
  .. .. .. ..$ : int 365
  .. .. .. ..$ : int 367
  .. .. .. ..$ : int 370
  .. .. .. ..$ : int 406
  .. .. .. ..$ : int [1:3] 285 319 320
  .. .. .. ..$ : int [1:2] 172 358
  .. .. .. ..$ : int 97
  .. .. .. ..$ : int 428
  .. .. .. ..$ : int 430
  .. .. .. ..$ : int 160
  .. .. .. ..$ : int 273
  .. .. .. ..$ : int 347
  .. .. .. ..$ : int 411
  .. .. .. ..$ : int [1:3] 183 185 329
  .. .. .. ..$ : int [1:2] 22 24
  .. .. .. ..$ : int [1:2] 134 331
  .. .. .. ..$ : Named num 253
  .. .. .. .. ..- attr(*, "names")= chr "253"
  .. .. .. ..$ : Named num 309
  .. .. .. .. ..- attr(*, "names")= chr "309"
  .. .. .. ..$ : int [1:3] 228 306 310
  .. .. .. ..$ : int [1:3] 153 212 213
  .. .. .. ..$ : int [1:2] 201 203
  .. .. .. ..$ : int [1:2] 336 424
  .. .. .. ..$ : int 86
  .. .. .. ..$ : int 88
  .. .. .. ..$ : int 102
  .. .. .. ..$ : int 112
  .. .. .. ..$ : int 231
  .. .. .. .. [list output truncated]
  .. .. ..$ psSamples     : int [1:274] 4 1 10 4 10 10 10 7 3 3 ...
  .. .. ..$ isotopes      :List of 431
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 1
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 1
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 2
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 3
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 2
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 3
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 4
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 4
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 5
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 5
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 6
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 6
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 7
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 7
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 8
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 8
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 9
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 9
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 10
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 10
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 11
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 11
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 12
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 12
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 13
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 13
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 14
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 14
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 15
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 16
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 15
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 16
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 17
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 18
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 18
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 17
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 19
  .. .. .. .. ..$ iso   : num 0
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 4
  .. .. .. .. ..$ y     : num 19
  .. .. .. .. ..$ iso   : Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "iso"
  .. .. .. .. ..$ charge: Named num 1
  .. .. .. .. .. ..- attr(*, "names")= chr "charge"
  .. .. .. .. ..$ val   : Named num 0
  .. .. .. .. .. ..- attr(*, "names")= chr "intrinsic"
  .. .. .. ..$ : NULL
  .. .. .. .. [list output truncated]
  .. .. ..$ ruleset       :'data.frame':	45 obs. of  7 variables:
  .. .. .. ..$ name    : chr [1:45] "[M-H]-" "[M-2H]2-" "[M-2H+Na]-" "[M-H+Cl]2-" ...
  .. .. .. ..$ nmol    : num [1:45] 1 1 1 1 1 1 1 2 2 2 ...
  .. .. .. ..$ charge  : num [1:45] -1 -2 -1 -2 -1 -1 -2 -1 -2 -1 ...
  .. .. .. ..$ massdiff: num [1:45] -1.01 -2.01 20.97 33.96 36.95 ...
  .. .. .. ..$ oidscore: num [1:45] 1 2 4 5 6 7 8 1 2 4 ...
  .. .. .. ..$ quasi   : num [1:45] 1 0 0 0 0 1 0 0 0 0 ...
  .. .. .. ..$ ips     : num [1:45] 1 1 0.25 1 0.25 1 1 0.5 0.5 0.25 ...
  .. .. ..$ derivativeIons:List of 431
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 6
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M+Cl]-"
  .. .. .. .. .. ..$ mass   : Named num 230
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 5
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-2H+K]-"
  .. .. .. .. .. ..$ mass   : Named num 230
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 6
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M+Cl]-"
  .. .. .. .. .. ..$ mass   : Named num 247
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 2
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 5
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-2H+K]-"
  .. .. .. .. .. ..$ mass   : Named num 247
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 6
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M+Cl]-"
  .. .. .. .. .. ..$ mass   : Named num 249
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 5
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-2H+K]-"
  .. .. .. .. .. ..$ mass   : Named num 249
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 2
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 6
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M+Cl]-"
  .. .. .. .. .. ..$ mass   : Named num 254
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 1
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-H]-"
  .. .. .. .. .. ..$ mass   : Named num 290
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 1
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-H]-"
  .. .. .. .. .. ..$ mass   : Named num 292
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 3
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-2H+Na]-"
  .. .. .. .. .. ..$ mass   : Named num 292
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 39
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-H-HCOOH]-"
  .. .. .. .. .. ..$ mass   : Named num 364
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 1
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-H]-"
  .. .. .. .. .. ..$ mass   : Named num 318
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 38
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-H-CO2]-"
  .. .. .. .. .. ..$ mass   : Named num 364
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ : NULL
  .. .. .. ..$ :List of 1
  .. .. .. .. ..$ :List of 5
  .. .. .. .. .. ..$ rule_id: Named num 3
  .. .. .. .. .. .. ..- attr(*, "names")= chr "ruleID"
  .. .. .. .. .. ..$ charge : num -1
  .. .. .. .. .. ..$ nmol   : num 1
  .. .. .. .. .. ..$ name   : chr "[M-2H+Na]-"
  .. .. .. .. .. ..$ mass   : Named num 318
  .. .. .. .. .. .. ..- attr(*, "names")= chr "mass"
  .. .. .. .. [list output truncated]
  .. ..$ peaks              : num [1:5711, 1:11] 430 540 350 554 567 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:5711] "CP0001" "CP0002" "CP0003" "CP0004" ...
  .. .. .. ..$ : chr [1:11] "mz" "mzmin" "mzmax" "rt" ...
  .. ..$ groups             : num [1:431, 1:10] 253 253 254 255 255 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr [1:10] "mzmed" "mzmin" "mzmax" "rtmed" ...
  .. ..$ groupmat           : num [1:431, 1:10] 253 253 254 255 255 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr [1:10] "mzmed" "mzmin" "mzmax" "rtmed" ...
  .. ..$ camera_output      :'data.frame':	431 obs. of  19 variables:
  .. .. ..$ mz        : num [1:431] 253 253 254 255 255 ...
  .. .. ..$ mzmin     : num [1:431] 253 253 254 255 255 ...
  .. .. ..$ mzmax     : num [1:431] 253 253 254 255 255 ...
  .. .. ..$ rt        : num [1:431] 503 498 498 511 577 ...
  .. .. ..$ rtmin     : num [1:431] 502 497 497 511 570 ...
  .. .. ..$ rtmax     : num [1:431] 509 499 499 523 578 ...
  .. .. ..$ CD-6KUCT  : num [1:431] 25638 13603971 1813356 1172629 190357 ...
  .. .. ..$ CD-77FXR  : num [1:431] NA 13982939 1910121 46756 1123058 ...
  .. .. ..$ CD-9OS5Y  : num [1:431] 610779 16205526 2227732 1028948 1067952 ...
  .. .. ..$ CD-9WOBP  : num [1:431] NA 13154606 1876157 NA 2067358 ...
  .. .. ..$ HC-9SNJ4  : num [1:431] 15684 11734941 1632058 592468 NA ...
  .. .. ..$ HC-9X47O  : num [1:431] 253934 9033802 1212304 717266 NA ...
  .. .. ..$ HC-AMR37  : num [1:431] 752963 5449290 541015 14008 40023 ...
  .. .. ..$ HC-AUP8B  : num [1:431] 456517 11455149 1569154 685798 NA ...
  .. .. ..$ QC_PREFA02: num [1:431] 543985 12002263 1599893 5761790 745308 ...
  .. .. ..$ QC_PREFB02: num [1:431] 527491 11919343 1663142 5788800 764796 ...
  .. .. ..$ isotopes  : chr [1:431] "" "[1][M]-" "[1][M+1]-" "[2][M]-" ...
  .. .. ..$ adduct    : chr [1:431] "" "" "" "" ...
  .. .. ..$ pcgroup   : chr [1:431] "239" "11" "138" "197" ...
  .. ..$ peak_result_summary: chr [1:10, 1:6] "CD-6KUCT" "CD-77FXR" "CD-9OS5Y" "CD-9WOBP" ...
  ..@ dataSet         :'data.frame':	432 obs. of  11 variables:
  .. ..$ Sample    : chr [1:432] "Label" "253.0509@502.69" "253.1445@498.26" "254.1478@498.26" ...
  .. ..$ CD-6KUCT  : chr [1:432] "Control" "25638.15531" "13603971.11044" "1813356.13912" ...
  .. ..$ CD-77FXR  : chr [1:432] "Control" NA "13982938.51448" "1910120.95123" ...
  .. ..$ CD-9OS5Y  : chr [1:432] "Control" "610779.40232" "16205526.08603" "2227732.44567" ...
  .. ..$ CD-9WOBP  : chr [1:432] "Control" NA "13154606.30984" "1876157.32445" ...
  .. ..$ HC-9SNJ4  : chr [1:432] "Disease" "15683.72062" "11734941.47374" "1632057.69695" ...
  .. ..$ HC-9X47O  : chr [1:432] "Disease" "253933.79128" "9033801.95777" "1212303.76631" ...
  .. ..$ HC-AMR37  : chr [1:432] "Disease" "752962.62186" "5449290.30062" "541014.59873" ...
  .. ..$ HC-AUP8B  : chr [1:432] "Disease" "456516.62758" "11455148.76817" "1569153.85546" ...
  .. ..$ QC_PREFA02: chr [1:432] "QC" "543984.71889" "12002263.16849" "1599892.93381" ...
  .. ..$ QC_PREFB02: chr [1:432] "QC" "527490.76195" "11919342.76482" "1663142.41353" ...
  ..@ runningplan     : list()
  ..@ msgSet          : list()
  ..@ userpath        : chr(0) 
  ..@ WorkingDir      : chr(0) 
}
\references{
https://www.nature.com/articles/s41586-019-1237-9
}
\examples{
data(mSet)
## maybe str(mSet) ; plot(mSet) ...
}
\keyword{datasets}
