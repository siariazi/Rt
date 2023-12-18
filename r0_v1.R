###############################################################################
# library imports 

library(OOPidemic)
# SEIR model estimating R0
###############################################################################

ref_strain <- ReferenceStrain$new(
    name = "ref_strain", # just a name for strain
    g_len = 1000, # genome length
    mut_rate = 1 / 7000 # one mutation be a week
)


runs <- 1000
i <- 1
results <- c()
for (i in seq(1:runs)) {

    # Group is a class
    group <- Group$new(
        id = 1,
        ref_strain = ref_strain,
        init_inf = 3,
        init_sus = 497,
        max_init_dist = 3, # up to 3 mutations in pathogen compared to the ref strain 
        inf_rate = 0.25, # how many people get infected at a particular time 
        # inc_shape = 0, # the incubation shape parm is what the default value is in the package
        rec_shape = 14, rec_rate = 2. # gamma distn 
    )


    index_cases <- group$infectious_hosts(0)
    index_rec_times <- vapply(index_cases, function(h) h$recovery_time, numeric(1L))

    while (group$time < max(index_rec_times)) {
        group$infect()
    }

    cases_infected <- sum(vapply(index_cases, function(i) length(i$infectees), integer(1L)))

    cat("Run:       ", i, "\n")
    cat("   Cases:  ", cases_infected, "\n")
    cat("   R0:     ", cases_infected / group$init_inf, "\n")

    results <- c(results, cases_infected / group$init_inf)
}

print(summary(results))

"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0000  0.3333  0.6667  0.8247  1.0000  3.0000 "
