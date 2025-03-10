# Omics data ----
if(config_list$test_mode == "T"){
  # Test files
  load("Input/03/exp.Rda")
  load("Input/03/cnv.Rda")
  load("Input/03/meth.Rda")
} else{
  load("Input/01/exp.Rda")
  load("Input/01/cnv.Rda")
  load("Input/01/meth.Rda")
}

# Small file
load("Input/01/protein.Rda")
load("Input/01/fusion.Rda")
load("Input/01/mut.Rda")

# Drug and annotation data ----
# load("Input/02/drug2.Rda")
# load("Input/05/pdo_anno.Rda") 

# modified new anno and drug 
load("Input/02/drug_0309.Rda")
# load("Input/05/pdo_anno_new.Rda") 
load("Input/06/anno_0309.Rda") 

# PDO data ----
load("Input/05/pdo_deng_0309.Rda")

# Finished Plot or preplot obj ----

