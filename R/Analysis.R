library("nhanesdata")

## process the data -- minus the accelerometry
# process_accel()
data("PAXINTEN_C");data("PAXINTEN_D")
process_flags()
process_mort()
process_covar()


library("magrittr");library("dplyr")
## merge covariate, mortality, and accelerometry data
AllAct_C <- left_join(PAXINTEN_C, Mortality_C, by = "SEQN") %>%
                left_join(Covariate_C, by="SEQN")


## subset the data


## re-weight the data
