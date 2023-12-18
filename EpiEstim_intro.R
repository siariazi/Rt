# EpiEstim demostration

library(EpiEstim)
library(ggplot2)

## load data
data(Flu2009)

## incidence:
head(Flu2009$incidence)

## serial interval (SI) distribution:
Flu2009$si_distr

head(Flu2009$si_data)

library(incidence)
plot(as.incidence(Flu2009$incidence$I, dates = Flu2009$incidence$dates))

res_parametric_si <- estimate_R(Flu2009$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 2.6, 
                                  std_si = 1.5))
)
#> Default config will estimate R on weekly sliding windows.
#>     To change this change the t_start and t_end arguments.

head(res_parametric_si$R)

plot(res_parametric_si, legend = FALSE)

res_non_parametric_si <- estimate_R(Flu2009$incidence, 
                                    method="non_parametric_si",
                                    config = make_config(list(
                                      si_distr = Flu2009$si_distr))
)

plot(res_non_parametric_si, "R")
