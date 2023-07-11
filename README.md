# Analysis for 
# _Experimental drought reduces the productivity and stability of a recovering calcareous grassland_
#### Jackson, J., Middleton, S. L., Lawson, C. S., Jardine, E., Hawes, N., Maseyk, K., Salguero-Gómez, R., & Hector, A.

#### 2023-06-27
#### Repository created by John Jackson

---

[![DOI](https://zenodo.org/badge/659275253.svg)](https://zenodo.org/badge/latestdoi/659275253)

Further acknowledgement for support of this study: Melanie Stone, David Gowing, Julia Haynes, Lara Clements, Laura McManus, Hannah King, David Encarnation, Abir Patwary, and Lauren Hinchcliffe.

This directory contains scripts and analysis data for our work on the experimental manipulation of precipitation and its impact on biodiversity and productivity. For the manuscript please see {CHANGE} [the bioRxiv entry](https://www.biorxiv.org/content/10.1101/2022.03.08.483493v3). Package version info for this analysis is given below.

Analysis scripts can be found in the `code/` sub-repository, manuscript figures and output in the `output/` sub-repository, and analysis data in the `data/` sub-repository. Raw biodiversity data (which is part of a global network) is available on request from Andrew Hector.

Scripts are labeled A-G in order of the analysis, and are as follows:

1. `A_data_cleaning_aquisition.R` - Compiling raw data sheets, data cleaning and wrangling for the raw biodiversity data.
2. `B_mixed_effects_models.R` - Bayesian hierarchical mixed-effects models to explore how biodiversity measures change with respect to precipitation treatments. Implemented in `brms`
3. `C_community_analyses.R` - Community stability and composition analysis for exploring how the grassland communities changed with precipitation treatments and through time. Implemented in `vegan`.
4. `D_gamma_models.R` - Supplementary Bayesian hierarchical mixed-effects models to explore how productivity changes with respect to precipitation treatments, but using a Gamma distribution to describe productivity.
5. `E_weather_effect_biodiversity.R` - Analysis exploring how biodiversity trends were impacted by local weather variation between 2016-2020, where data was compiled from the Ecological Change Network.
6. `F_manuscript_figures.R` - Script to recreate data figures presented in the manuscript.
7. `G_manuscript_tables.R` - Script to recreate model selection tables used in the modelling framework.

## System Information and Package Versions

<details>
  <summary>Click here to expand</summary>
  
```
R version 4.3.0 (2023-04-21 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19044)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] flextable_0.9.1   psych_2.3.6       MASS_7.3-58.4     gghalves_0.1.4    tidybayes_3.0.4  
 [6] GGally_2.1.2      vegan_2.6-4       lattice_0.21-8    permute_0.9-7     brms_2.19.0      
[11] Rcpp_1.0.10       viridis_0.6.3     viridisLite_0.4.2 patchwork_1.1.2   lubridate_1.9.2  
[16] forcats_1.0.0     stringr_1.5.0     dplyr_1.1.2       purrr_1.0.1       readr_2.1.4      
[21] tidyr_1.3.0       tibble_3.2.1      ggplot2_3.4.2     tidyverse_2.0.0  

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      tensorA_0.36.2          rstudioapi_0.14        
  [4] jsonlite_1.8.5          magrittr_2.0.3          estimability_1.4.1     
  [7] farver_2.1.1            rmarkdown_2.22          ragg_1.2.5             
 [10] vctrs_0.6.2             askpass_1.1             base64enc_0.1-3        
 [13] htmltools_0.5.5         distributional_0.3.2    curl_5.0.0             
 [16] StanHeaders_2.26.27     htmlwidgets_1.6.2       plyr_1.8.8             
 [19] emmeans_1.8.6           zoo_1.8-12              uuid_1.1-0             
 [22] igraph_1.4.3            mime_0.12               lifecycle_1.0.3        
 [25] pkgconfig_2.0.3         colourpicker_1.2.0      Matrix_1.5-4           
 [28] R6_2.5.1                fastmap_1.1.1           shiny_1.7.4            
 [31] digest_0.6.31           colorspace_2.1-0        reshape_0.8.9          
 [34] ps_1.7.5                textshaping_0.3.6       crosstalk_1.2.0        
 [37] fansi_1.0.4             timechange_0.2.0        abind_1.4-5            
 [40] mgcv_1.8-42             compiler_4.3.0          fontquiver_0.2.1       
 [43] bit64_4.0.5             withr_2.5.0             backports_1.4.1        
 [46] inline_0.3.19           shinystan_2.6.0         pkgbuild_1.4.1         
 [49] openssl_2.0.6           gfonts_0.2.0            gtools_3.9.4           
 [52] loo_2.6.0               tools_4.3.0             zip_2.3.0              
 [55] httpuv_1.6.11           threejs_0.3.3           glue_1.6.2             
 [58] callr_3.7.3             nlme_3.1-162            promises_1.2.0.1       
 [61] grid_4.3.0              checkmate_2.2.0         cluster_2.1.4          
 [64] reshape2_1.4.4          generics_0.1.3          gtable_0.3.3           
 [67] tzdb_0.4.0              data.table_1.14.8       hms_1.1.3              
 [70] xml2_1.3.4              utf8_1.2.3              pillar_1.9.0           
 [73] ggdist_3.3.0            markdown_1.7            vroom_1.6.3            
 [76] posterior_1.4.1         later_1.3.1             splines_4.3.0          
 [79] bit_4.0.5               tidyselect_1.2.0        fontLiberation_0.1.0   
 [82] miniUI_0.1.1.1          knitr_1.43              fontBitstreamVera_0.1.1
 [85] arrayhelpers_1.1-0      gridExtra_2.3           V8_4.3.0               
 [88] crul_1.4.0              stats4_4.3.0            xfun_0.39              
 [91] bridgesampling_1.1-2    matrixStats_1.0.0       DT_0.28                
 [94] rstan_2.26.22           stringi_1.7.12          httpcode_0.3.0         
 [97] evaluate_0.21           codetools_0.2-19        officer_0.6.2          
[100] gdtools_0.3.3           cli_3.6.1               RcppParallel_5.1.7     
[103] systemfonts_1.0.4       shinythemes_1.2.0       xtable_1.8-4           
[106] munsell_0.5.0           processx_3.8.1          coda_0.19-4            
[109] svUnit_1.0.6            parallel_4.3.0          rstantools_2.3.1       
[112] ellipsis_0.3.2          prettyunits_1.1.1       dygraphs_1.1.1.6       
[115] bayesplot_1.10.0        Brobdingnag_1.2-9       mvtnorm_1.2-1          
[118] scales_1.2.1            xts_0.13.1              crayon_1.5.2           
[121] rlang_1.1.1             mnormt_2.1.1            shinyjs_2.1.0        
```

</details>
