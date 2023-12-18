#library(devtools)
library(OOPidemic)

# set up a reference strain with a randomised genome
ref_strain <- ReferenceStrain$new(
  name = "ref_strain",
  g_len = 1000
)

# set up a group 
group <- Group$new(
  id = 1,
  ref_strain = ref_strain
)

# set up a lab to take samples 
lab <- Lab$new()

# run simluation
group$run_simulation(lab)
