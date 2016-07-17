#Let us make a latex Table
library(xtable)

#Lets put in the data
my_df <- data.frame(ID = c("Experiment.1","Experiment.2.Cycle.1","Experiment.2.Cycle.2"),
                    Pu.Recovery =c(76,94,90),
                    Pu.Err = c(0.3,0.9,0.6),
                    U.Recovery=c(25,7,5),
                    U.Err =c(0.03,0.06,0.10))


my_df

print(xtable(my_df), type='latex')
