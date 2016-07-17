#Let us make a latex Table
library(xtable)

#Lets put in the data
my_df <- data.frame(ID = c("U","Rb","Cs","Sr","Ba","Mo","Ru","Pd","Cd","Sn","Ce","Nd","Pm","Sm","Eu"),
                    Atomic.Number =c(92,37,55,38,56,42,44,46,48,50,58,60,61,62,63),
                    Exp.1 = c(6.85,32,146,233.5,344,20.67,49,65,61,7.45,35.24,16.37,10.7,9.94,8.4),
                    Exp.1.err=c(0.46,1.55,7.58,12.74,200,2.03,1.9,14.3,6.6,0.43,1.68,0.65,0.66,0.25,0.49),
                    Exp.2 =c(15.08,1.84,11.92,38.26,0.39,1.19,2.84,3.62,3.5,13.85,3.2,5.94,3.3,2.5,2.6),
                    Exp.2.err=c(0.60,0.26,0.96,2.23,50,0.25,0.111,0.94,0.98,1.29,0.67,2.01,0.5,0.19,0.23),
                    Isotopes.Used=c("238U","85Rb","133Cs","90Sr","Ba","97Mo","101Ru 102Ru 104Ru","108 110Pd","Cd112","119Sn","140,142Ce","Nd143","147Pm","151Sm","154Eu"))



my_df

my_df <- my_df[order(my_df$Atomic.Number),]

my_df

mean(my_df$Exp.2[-15]/my_df$Exp.1[-15])

print(xtable(my_df), type='latex')

