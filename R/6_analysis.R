
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
## ADD STEFFENSEN ET AL. () STABILITY POINTS
points()
## ADD WILDHABER ET AL. (2016) STABILITY POINTS