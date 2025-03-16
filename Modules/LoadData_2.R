if(config_list$test_mode == "T"){
  # Test files
  load("Input/03/exp.Rda")
  load("Input/03/cnv.Rda")
  load("Input/03/meth.Rda")
  load("Input/03/protein.Rda")
  load("Input/02/drug_0314.Rda")
  
} else{
  load("Input/01/exp.Rda")
  load("Input/01/cnv.Rda")
  load("Input/01/meth.Rda")
  load("Input/01/protein.Rda")
  load("Input/02/drug_0314.Rda")
}