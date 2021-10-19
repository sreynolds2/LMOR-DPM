params<- readRDS("./baseline-parameters/default-parameters.rds")

inps<- list()
## MAXIMUM AGE
inps$max_age<- params$max_age
## AGE-0 SURVIVAL
inps$mu$phi0<- 0.000075
inps$sd$phi0<- 0
inps$sd_temp$phi0<- 0
inps$lb$phi0<- 0
inps$ub$phi0<- 0.0004
# ## LARVAL SURVIVAL
# inps$mu$phi_larv<- 0.0510
# inps$sd$phi_larv<- sqrt(1.3055*10^(-4))
# inps$sd_temp$phi_larv<- 0
# inps$lb$phi_larv<- 0
# inps$ub$phi_larv<- 1
## AGE-1+ SURVIVALS
### JAGS VBGF FIT
inps$vbgf$mu$phi<- params$phi$vbgf
inps$vbgf$sd$phi<- rep(0, length(params$phi$vbgf))
inps$vbgf$sd_temp$phi<- rep(0, length(params$phi$vbgf))
inps$vbgf$lb$phi<- rep(0, length(params$phi$vbgf))
inps$vbgf$ub$phi<- rep(1, length(params$phi$vbgf))
### TMB ln(VBGF) FIT
inps$ln_vbgf$mu$phi<- params$phi$ln_vbgf
inps$ln_vbgf$sd$phi<- rep(0, length(params$phi$ln_vbgf))
inps$ln_vbgf$sd_temp$phi<- rep(0, length(params$phi$ln_vbgf))
inps$ln_vbgf$lb$phi<- rep(0, length(params$phi$ln_vbgf))
inps$ln_vbgf$ub$phi<- rep(1, length(params$phi$ln_vbgf))
## SEX RATIO
inps$mu$probF<- 0.5
inps$sd$probF<- 0
inps$sd_temp$probF<- 0
inps$lb$probF<- 0
inps$ub$probF<- 1
## AGE OF FIRST FEMALE REPRODUCTION
inps$mu$mat_age<- params$a_min
inps$sd$mat_age<- 0
inps$sd_temp$mat_age<- 0
inps$lb$mat_age<- 1
inps$ub$mat_age<- params$max_age
  ### IF THIS VARIES, WILL WANT TO RUN NEW SURVIVAL, 
  ### MATURATION, AND PSI VALUES
## PROPORTION SPAWNING
### JAGS VBGF FIT
inps$vbgf$mu$p_spawn<- params$psi$vbgf
inps$vbgf$sd$p_spawn<- rep(0, length(params$psi$vbgf))
inps$vbgf$sd_temp$p_spawn<- rep(0, length(params$psi$vbgf))
inps$vbgf$lb$p_spawn<- rep(0, length(params$psi$vbgf))
inps$vbgf$ub$p_spawn<- rep(1, length(params$psi$vbgf))
### TMB ln(VBGF) FIT
inps$ln_vbgf$mu$p_spawn<- params$psi$ln_vbgf
inps$ln_vbgf$sd$p_spawn<- rep(0, length(params$psi$ln_vbgf))
inps$ln_vbgf$sd_temp$p_spawn<- rep(0, length(params$psi$ln_vbgf))
inps$ln_vbgf$lb$p_spawn<- rep(0, length(params$psi$ln_vbgf))
inps$ln_vbgf$ub$p_spawn<- rep(1, length(params$psi$ln_vbgf))
## GROWTH PARAMETERS --HOW DO WE WANT TO DEAL WITH DIFFERENCES IN EGGS???
### JAGS VBGF FIT
inps$vbgf$Linf<- params$growth$vbgf$Linf
inps$vbgf$k<- params$growth$vbgf$k
inps$vbgf$t0<- params$growth$vbgf$t0
inps$vbgf$fit$sd$l_at_a<- params$growth$vbgf$sigma
inps$vbgf$fit$sd_temp$l_at_a<- 0
#inps$vbgf$params$sd$l_at_a<- #vector of age specific total sd's
#inps$vbgf$params$sd_temp$l_at_a<- #vector of age specific temporal sd's
### TMB ln(VBGF) FIT
inps$ln_vbgf$Linf<- params$growth$ln_vbgf$Linf
inps$ln_vbgf$k<- params$growth$ln_vbgf$k
inps$ln_vbgf$t0<- params$growth$ln_vbgf$t0
inps$ln_vbgf$fit$sd$l_at_a<- params$growth$ln_vbgf$sigma
inps$ln_vbgf$fit$sd_temp$l_at_a<- 0
#inps$ln_vbgf$params$sd$l_at_a<- #vector of age specific total sd's
#inps$ln_vbgf$params$sd_temp$l_at_a<- #vector of age specific temporal sd's
## FECUNDITY PARAMETERS
inps$fec$b0<- params$fecundity$intercept 
inps$fec$b1<- params$fecundity$slope
inps$fec$mean_fl<- params$fecundity$mean_fl
inps$fec$sd_fl<- params$fecundity$sd_fl
inps$fec$fit$sd$ln_eggs_at_l<- params$fecundity$disp
inps$fec$fit$sd_temp$ln_eggs_at_l<- 0
#inps$fec$params$sd$eggs_at_l#vector of length specific total sd's
#inps$fec$params$sd_temp$eggs_at_l<- #vector of length specific temporal sd's
## MEAN EGGS
inps$vbgf$eggs<- params$eggs$vbgf[1:inps$max_age]
inps$ln_vbgf$eggs<- params$eggs$ln_vbgf[1:inps$max_age]
## INITIAL POPULATION
#inps$N0<- #load from file



