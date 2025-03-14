# Omics data ----
if(config_list$test_mode == "T"){
  # Test files
  load("Input/03/exp.Rda")
  load("Input/03/cnv.Rda")
  load("Input/03/meth.Rda")
  load("Input/03/protein.Rda")
} else{
  load("Input/01/exp.Rda")
  load("Input/01/cnv.Rda")
  load("Input/01/meth.Rda")
  load("Input/01/protein.Rda")
}

# Small file
load("Input/01/fusion.Rda")
load("Input/01/mut.Rda")

# Drug and annotation data ---- 
# modified new anno and drug 
load("Input/02/drug_0314.Rda")
load("Input/06/anno_0309.Rda") 


