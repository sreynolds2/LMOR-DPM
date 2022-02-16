
inputs<- readRDS("./baseline-parameters/default-parameters.rds")

# FIGURES
## SURVIVAL
par(mar=c(4,4,4,2),
    oma=c(1,1,0,0))
plot(2:40, inputs$phi$vbgf[2:40], pch=19, xlab="Age", ylab="Survival Probability")
points(inputs$phi$ln_vbgf[2:40], pch=21)#=19, col="blue")
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), #col=c("black", "blue"), 
       pch=c(19, 21), bty="n")
#### AVERAGE LENGTH AT AGE
age<-1:inputs$max_age
La_vbgf<- inputs$growth$vbgf$Linf*(1-exp(-inputs$growth$vbgf$k*(age-inputs$growth$vbgf$t0)))
La_ln<- inputs$growth$ln_vbg$Linf*(1-exp(-inputs$growth$ln_vbg$k*(age-inputs$growth$ln_vbg$t0)))
plot(age, La_vbgf, xlab="Age", ylab="Length (mm)", pch=19)
points(age, La_ln, pch=21)#col="blue", pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), #col=c("black", "blue"), 
       pch=c(19, 21), bty="n")
#### AVERAGE FECUNDITY PER SPAWNER
age<-1:41
plot(age, c(rep(-10000,inputs$a_min-1), 
            inputs$eggs$ln_vbgf[inputs$a_min:inputs$max_age]), ylim=c(0,35000),
     xlab="Age", ylab="Mean Number of Eggs per Spawning Female", col="blue", pch=19)
points(age, c(rep(-10000,inputs$a_min-1), 
              inputs$eggs$vbgf[inputs$a_min:inputs$max_age]), pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
       pch=19, bty="n")
#### AVERAGE FECUNDITY PER SPAWNER FIGV2
age<-1:41
plot(age, c(rep(-10000,inputs$a_min-1), 
            inputs$eggs$ln_vbgf[inputs$a_min:inputs$max_age]), ylim=c(0,35000),
     xlab="Age", ylab="Mean Number of Eggs per Spawning Female", pch=21)
points(age, c(rep(-10000,inputs$a_min-1), 
              inputs$eggs$vbgf[inputs$a_min:inputs$max_age]), pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "black"), 
       pch=c(19, 21), bty="n")
# #### AVERAGE FECUNDITY PER FEMALE (NOT ACCOUNTING FOR RR STATUS)
# age<-1:41
# plot(age, c(rep(-10000,inputs$a_min-1), 
#             inputs$allF_eggs$ln_vbgf[inputs$a_min:inputs$max_age]), ylim=c(0,35000),
#      xlab="Age", ylab="Mean Number of Eggs", col="blue", pch=19)
# points(age, c(rep(-10000,inputs$a_min-1), 
#               inputs$allF_eggs$vbgf[inputs$a_min:inputs$max_age]), pch=19)
# legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
#        pch=19, bty="n")
#### MATURATION AGE
age<-1:41
plot(age, inputs$m_i$FSM[1:41], col="green", xlab="Age", ylab="Probability of Maturing", pch=19, type="b")
points(age, inputs$m_i$ln_vbgf[1:41], col="blue", pch=19, type="b")
points(age, inputs$m_i$vbgf[1:41], pch=19, type="b")
legend("topright", c("JAGS VBGF", "TMB LN(VBGF)", "2019 FSM"), 
       col=c("black", "blue", "green"), lty=1, pch=19, bty="n")

#### MATURATION AGE FIGV2
age<-1:41
plot(age, inputs$m_i$vbgf[1:41], xlab="Age", ylab="Probability of Maturing", pch=19, type="b")
points(age, inputs$m_i$ln_vbgf[1:41], pch=21, type="b")
legend("topright", c("JAGS VBGF", "TMB LN(VBGF)"), 
       col="black", lty=1, pch=c(19, 21), bty="n")

#### PROPORTION OF REPRODUCTIVELY READY FEMALES
age<-1:41
plot(age, inputs$psi$FSM[1:41], col="green", xlab="Age", ylab="Proportion of Gravid Females", 
     pch=19, type="b")
points(age, inputs$psi$ln_vbgf[1:41], col="blue", pch=19, type="b")
points(age, inputs$psi$vbgf[1:41], pch=19, type="b")
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)", "2019 FSM"), 
       col=c("black", "blue", "green"), pch=19, lty=1, bty="n")
#### PROPORTION OF REPRODUCTIVELY READY FEMALES FIGV2
age<-1:41
plot(age, inputs$psi$vbgf[1:41], xlab="Age", ylab="Proportion of Gravid Females", 
     pch=19, type="b")
points(age, inputs$psi$ln_vbgf[1:41], pch=21, type="b")
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), 
       col="black", pch=c(19,21), lty=1, bty="n")
#### AVERAGE FECUNDITY PER FEMALE
age<-1:41
plot(age, inputs$psi$FSM[1:41]*inputs$eggs$ln_vbgf[1:41], xlab="Age", ylab="Eggs per Female", 
     col="lightblue", pch=19, type="b")
points(age, inputs$psi$FSM[1:41]*inputs$eggs$vbgf[1:41], col="gray", pch=19, type="b")
points(age, inputs$psi$vbgf[1:41]*inputs$eggs$vbgf[1:41], pch=19, type="b")
points(age, inputs$psi$ln_vbgf[1:41]*inputs$eggs$ln_vbgf[1:41], col="blue", pch=19, type="b")
legend("bottomright", c("JAGS VBGF Maturation & Growth", 
                        "TMB LN(VBGF) Maturation & Growth", 
                        "2019 FSM Maturation + JAGS VBGF Growth", 
                        "2019 FSM Maturation + TMB LN(VBGF) Growth"), 
       col=c("black", "blue", "gray", "lightblue"), pch=19, bty="n")




## FIGURE COMPARISON OF APPROACHES -- NEEDS TO BE RECODED
age<-1:inputs$max_age
### VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ADULTS ONLY
plot(age, fecundity_vbgf$Mean_Eggs_Adults, xlab="Age", ylab="Mean Number of Eggs", col="blue", pch=19)
points(age, fecundity_vbgf2$Mean_Eggs_Adults, pch=19)
points(age, fecundity_vbgf5$Mean_Eggs_Adults, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ALL BY AGE: ADULTS WEIGHTED
plot(age, fecundity_vbgf$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
points(age, fecundity_vbgf2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19)
points(age, fecundity_vbgf5$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### LN_VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ADULTS ONLY
plot(age, fecundity_ln$Mean_Eggs_Adults, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
points(age, fecundity_ln2$Mean_Eggs_Adults, pch=19)
points(age, fecundity_ln5$Mean_Eggs_Adults, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### LN_VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ALL BY AGE: ADULTS WEIGHTED
plot(age, fecundity_ln$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
points(age, fecundity_ln2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19)
points(age, fecundity_ln5$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### VBGF EGGS FROM ALL BY AGE: AVERAGE OF DRAWS VS. ADULTS WEIGHTED
plot(age, fecundity_vbgf2$Mean_Eggs_All, xlab="Age", ylab="Number of Eggs", pch=19)
points(age, fecundity_vbgf2$Mean_Eggs_Adults*p_stage$vbgf$adult, pch=19, col="blue")
legend("bottomright", c("Eggs from All Fish", "Eggs from Adults Weighted by Proportion"),
       col=c("black", "blue"), pch=19, bty="n")

### LN_VBGF EGGS FROM ALL BY AGE: AVERAGE OF DRAWS VS. ADULTS WEIGHTED
plot(age, fecundity_ln2$Mean_Eggs_All, xlab="Age", ylab="Number of Eggs", pch=19)
points(age, fecundity_ln2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="blue")
legend("bottomright", c("Eggs from All Fish", "Eggs from Adults Weighted by Proportion"),
       col=c("black", "blue"), pch=19, bty="n")

rm(dat, mean_fl, sd_fl, fit, intrcpt, slp, disp, n, a, fec_vbgf, fecundity_vbgf,
   fec_ln, fecundity_ln, length_at_age, eggs, age)


#### ADDITIONAL EGGS
E2<- read.csv("./baseline-parameters/fecundity_estimates_by_age_vbgf_2mil.csv")
E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_ln_vbgf_2mil.csv")
inps<- list()
inps$max_age<- 66
inps$a_min<- 8
##### AVERAGE FECUNDITY PER SPAWNER
age<-1:66
plot(age, c(rep(-10000,inps$a_min-1), 
            E$Mean_Eggs_Adults[inps$a_min:inps$max_age]), ylim=c(0,35000),
     xlab="Age", ylab="Mean Number of Eggs per Spawning Female", col="blue", pch=19)
points(age, c(rep(-10000,inps$a_min-1), 
              E2$Mean_Eggs_Adults[inps$a_min:inps$max_age]), pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
       pch=19, bty="n")



