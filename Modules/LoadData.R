# Omics data ----
if(config_list$test_mode == "T"){
  # Test files
  load("Input/04/mRNA.Rda")
  load("Input/04/cnv.Rda")
  load("Input/04/meth.Rda")
  load("Input/04/protein.Rda")
} else{
  load("Input/01/mRNA.Rda")
  load("Input/01/cnv.Rda")
  load("Input/01/meth.Rda")
  load("Input/01/protein.Rda")
}

# Small file
load("Input/01/fusion.Rda")
load("Input/01/mut.Rda")

# Drug and annotation data ---- 
# modified new anno and drug 
load("Input/02/drug.Rda")
load("Input/03/anno.Rda") 


