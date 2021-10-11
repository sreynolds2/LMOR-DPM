
source("./R/0_baseline-parameterization.r")

# FIGURES
## SURVIVAL
par(mar=c(4,4,4,2),
    oma=c(1,1,0,0))
plot(2:40, inputs$phi$vbgf[2:40], pch=19, xlab="Age", ylab="Survival Probability")
points(inputs$phi$ln_vbgf[2:40], col="blue", pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
       pch=19, bty="n")
#### AVERAGE LENGTH AT AGE
age<-1:inputs$max_age
La_vbgf<- inputs$growth$vbgf$Linf*(1-exp(-inputs$growth$vbgf$k*(age-inputs$growth$vbgf$t0)))
La_ln<- inputs$growth$ln_vbg$Linf*(1-exp(-inputs$growth$ln_vbg$k*(age-inputs$growth$ln_vbg$t0)))
plot(age, La_vbgf, xlab="Age", ylab="Length (mm)", pch=19)
points(age, La_ln, col="blue", pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
       pch=19, bty="n")
#### AVERAGE FECUNDITY PER SPAWNER
age<-1:41
plot(age, c(rep(-10000,inputs$a_min-1), 
            inputs$eggs$ln_vbgf[inputs$a_min:inputs$max_age]), ylim=c(0,35000),
     xlab="Age", ylab="Mean Number of Eggs per Spawning Female", col="blue", pch=19)
points(age, c(rep(-10000,inputs$a_min-1), 
              inputs$eggs$vbgf[inputs$a_min:inputs$max_age]), pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
       pch=19, bty="n")
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
plot(age, inputs$m_i$FSM, col="green", xlab="Age", ylab="Probability of Maturing", pch=19, type="b")
points(age, inputs$m_i$ln_vbgf, col="blue", pch=19, type="b")
points(age, inputs$m_i$vbgf, pch=19, type="b")
legend("topright", c("JAGS VBGF", "TMB LN(VBGF)", "FSM"), 
       col=c("black", "blue", "green"), lty=1, pch=19, bty="n")
#### PROPORTION OF REPRODUCTIVELY READY FEMALES
age<-1:41
plot(age, inputs$psi$FSM, col="green", xlab="Age", ylab="Proportion of Gravid Females", 
     pch=19, type="b")
points(age, inputs$psi$ln_vbgf, col="blue", pch=19, type="b")
points(age, inputs$psi$vbgf, pch=19, type="b")
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)", "FSM"), 
       col=c("black", "blue", "green"), pch=19, lty=1, bty="n")
#### AVERAGE FECUNDITY PER FEMALE
age<-1:41
plot(age, inputs$psi$FSM*inputs$eggs$ln_vbgf, xlab="Age", ylab="Eggs per Female", 
     col="lightblue", pch=19, type="b")
points(age, inputs$psi$FSM*inputs$eggs$vbgf, col="gray", pch=19, type="b")
points(age, inputs$psi$vbgf*inputs$eggs$vbgf, pch=19, type="b")
points(age, inputs$psi$ln_vbgf*inputs$eggs$ln_vbgf, col="blue", pch=19, type="b")
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)", "FSM Maturation + JAGS VBGF", 
                        "FSM Maturation + TMB LN(VBGF)"), 
       col=c("black", "blue", "gray", "lightblue"), pch=19, bty="n")
