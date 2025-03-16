# # Drug datasets - apply row-level z-score normalization
# if (exists("ctrp1_drug", envir = .GlobalEnv)) assign("ctrp1_drug", zscore_normalize_drug(base::get("ctrp1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("ctrp2_drug", envir = .GlobalEnv)) assign("ctrp2_drug", zscore_normalize_drug(base::get("ctrp2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("prism_drug", envir = .GlobalEnv)) assign("prism_drug", zscore_normalize_drug(base::get("prism_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("gdsc1_drug", envir = .GlobalEnv)) assign("gdsc1_drug", zscore_normalize_drug(base::get("gdsc1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("gdsc2_drug", envir = .GlobalEnv)) assign("gdsc2_drug", zscore_normalize_drug(base::get("gdsc2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("gCSI_drug", envir = .GlobalEnv)) assign("gCSI_drug", zscore_normalize_drug(base::get("gCSI_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("deng1_drug", envir = .GlobalEnv)) assign("deng1_drug", zscore_normalize_drug(base::get("deng1_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("deng2_drug", envir = .GlobalEnv)) assign("deng2_drug", zscore_normalize_drug(base::get("deng2_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# if (exists("deng3_drug", envir = .GlobalEnv)) assign("deng3_drug", zscore_normalize_drug(base::get("deng3_drug", envir = .GlobalEnv)), envir = .GlobalEnv)
# 
# # Save
# save(
#   prism_drug,
#   ctrp1_drug,
#   ctrp2_drug,
#   gCSI_drug,
#   gdsc1_drug,
#   gdsc2_drug,
#   deng1_drug,
#   deng2_drug,
#   deng3_drug,
#   file = "Input/04/drug_zscore_3016.Rda"
# )
