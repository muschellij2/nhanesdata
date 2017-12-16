# library("nhanesdata")
#
# ## process the data -- minus the accelerometry
# # process_accel()
# data("PAXINTEN_C");data("PAXINTEN_D")
# process_flags()
# process_mort()
# process_covar()
#
#
# library("magrittr");library("dplyr")
# ## apply accelerometry exclusion criteria
#
#
# ## merge covariate, mortality, and accelerometry data
# AllAct_C <- left_join(PAXINTEN_C, Mortality_C, by = "SEQN") %>%
#                 left_join(Covariate_C, by="SEQN")
# AllAct_D <- left_join(PAXINTEN_D, Mortality_D, by = "SEQN") %>%
#                 left_join(Covariate_D, by="SEQN")
# rm(list=c(paste0(c("PAXINTEN_", "Covariate_","Mortality_","Flags_"),rep(LETTERS[3:4],each=4))))
#
# # ## subset the data
# #
# #
# # ## re-weight the data
