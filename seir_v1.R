###############################################################################
# library imports 

library(OOPidemic)
# SEIR model generating fasta file and meta data

set.seed(543)

###############################################################################

ref_strain <- ReferenceStrain$new(
    name = "ref_strain",
    g_len = 1000,
    mut_rate = 1 / 7000
)

group <- Group$new(
    id = 1,
    ref_strain = ref_strain,
    init_inf = 3,
    init_sus = 497,
    max_init_dist = 3,
    inf_rate = 0.25,
    # inc_shape = 0,
    rec_shape = 14, rec_rate = 2
)

lab <- Lab$new()

# save index cases
index_cases <- group$infectious_hosts()

lab$sample_hosts(group$hosts_due_for_sampling, group$time)

while (all(
    group$is_outbreak_active,
    group$recovered_size < 100,
    group$time < 90
)) {
    group$infect()
    lab$sample_hosts(group$hosts_due_for_sampling, group$time)
}

cd_plot <- plot_compartment_dynamics(group)
plot(cd_plot)

print(group)
print(lab)

# calculate how many hosts the index_cases infected
num_hosts_inf <- vapply(index_cases, function(i) length(i$infectees), integer(1L)) 
cat("Hosts infected by index cases:   ", sum(num_hosts_inf), "\n")
cat("R0 =                             ", sum(num_hosts_inf) / group$init_inf, "\n")

lab$save("data", "seir_v1")
