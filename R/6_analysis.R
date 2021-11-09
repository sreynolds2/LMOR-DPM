
source("./R/1_global.r")
source("./R/0_create-inputs.r")
source("./R/2_functions.r")


## BASELINE MATRIX ANALYSIS
ea_vbgf<- matrix_eigen_analysis(inputs=params,
                           growth_type="vbgf")
ea_ln<- matrix_eigen_analysis(inputs=params,
                              growth_type="ln_vbgf")
ea<- list(vbgf=ea_vbgf, ln_vbgf=ea_ln)
saveRDS(ea, "./output/baseline-matrix-results.r")

### SENSITIVITIES AND ELASTICITIES
tbl_vbgf<- sens_elas_table(ea_vbgf)
barplot(tbl_vbgf$Sensitivities$Sensitivity[5:1],
        horiz = TRUE, xlab="Sensitivity", 
        names.arg=tbl_vbgf$Sensitivities$Parameter[5:1], 
        las=1, main="JAGS Von Bertalanffy Growth Fit")

tbl_ln<- sens_elas_table(ea_ln)
barplot(tbl_ln$Sensitivities$Sensitivity[5:1],
        horiz = TRUE, xlab="Sensitivity", 
        names.arg=tbl_ln$Sensitivities$Parameter[5:1], 
        las=1, main="TMB von Bertalanffy Growth Fit Log Scale")

all(tbl_ln$Sensitivities$Parameter[5:1]==tbl_vbgf$Sensitivities$Parameter[5:1])
M<- matrix(c(tbl_ln$Sensitivities$Sensitivity[5:1],
             tbl_vbgf$Sensitivities$Sensitivity[5:1]), 
           nrow = 2, byrow=TRUE)
barplot(M, horiz = TRUE, xlab="Sensitivity", 
        names.arg=tbl_ln$Sensitivities$Parameter[5:1], 
        las=1, beside=TRUE, col=c("lightgray", "black"))
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), 
       fill=c("black", "lightgray"), bty="n")

### AGE-0 SURVIVAL FOR STABILITY
bp_vbgf<- boundary_vals(inputs=params,
                        growth_type="vbgf")
bp_vbgf$boundary
bp_ln<- boundary_vals(inputs=params,
                      growth_type="ln_vbgf")
bp_ln$boundary

c_vbgf<- phi0_phi1_sexratio_curves(bp_vbgf,
                                   probF=0.5)
c_ln<- phi0_phi1_sexratio_curves(bp_ln,
                                 probF=0.5)
c_dat<- rbind(c_vbgf, c_ln)

plot_boundary_curves(c_dat,
                     phi1_upper = 0.96)
## ADD UPPER BOUND FOR AGE-0 SURVIVAL FROM PINE ET AL. 2001 
abline(0.0004, 0, col="lightgray")
## ADD STEFFENSEN ET AL. (2013) STABILITY POINTS
points(0.686, 0.00011, pch=17)
## ADD WILDHABER ET AL. (2017) STABILITY POINTS
points(0.3674,0.00011, pch=19)
legend(0.73, 0.0009, c("Steffensen et al. 2013", "Wildhaber et al. 2017"),
       pch=c(17,19), bty='n')


### UNLINKING SURVIVAL AND MATURATION FROM GROWTH
inps2<- params
inps2$phi$vbgf<- c(0.686, rep(0.922, inps2$max_age))
inps2$phi$ln_vbgf<- c(0.686, rep(0.922, inps2$max_age))
inps2$m_i$vbgf<- inps2$m_i$FSM
inps2$m_i$ln_vbgf<- inps2$m_i$FSM
inps2$psi$vbgf<- inps2$psi$FSM
inps2$psi$ln_vbgf<- inps2$psi$FSM
inps2$gamma<- 1

bp_vbgf2<- boundary_vals(inputs=inps2,
                         growth_type="vbgf")
bp_vbgf2$boundary
bp_ln2<- boundary_vals(inputs=inps2,
                       growth_type="ln_vbgf")
bp_ln2$boundary

abs(bp_vbgf2$boundary$phi0-bp_ln2$boundary$phi0)


E_vbgf<- read.csv("./baseline-parameters/fecundity_estimates_by_age_no_LB_vbgf.csv")
E_ln<- read.csv("./baseline-parameters/fecundity_estimates_by_age_no_LB_ln_vbgf.csv")
inps2$eggs$vbgf<- E_vbgf$Mean_Eggs_Adults
inps2$eggs$ln_vbgf<- E_ln$Mean_Eggs_Adults

bp_vbgf3<- boundary_vals(inputs=inps2,
                         growth_type="vbgf")
bp_vbgf3$boundary
bp_ln3<- boundary_vals(inputs=inps2,
                       growth_type="ln_vbgf")
bp_ln3$boundary

abs(bp_vbgf3$boundary$phi0-bp_ln3$boundary$phi0)
