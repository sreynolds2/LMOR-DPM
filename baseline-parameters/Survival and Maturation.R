
surv_and_mat_from_stage<- function(max_age=41, # MAXIMUM FEMALE AGE OF REPRODUCTION
                                   a_min=8,  # MINIMUM MATURATION AGE
                                   # MINIMUM STAGE BASED LENGTHS
                                   min_adult_length=800,
                                   min_sa_length=600,
                                   # SURVIVAL PARAMETERS BY STAGE
                                   phi_1=0.151,
                                   phi_juv=0.780,
                                   phi_sa=0.932,
                                   phi_adult=0.976,
                                   # TRANSITION PARAMETERS
                                   j2j=0.893,
                                   sa2sa=0.897)
{
  # INITIALIZE LISTS
  stage<- list()
  growth<- list()
  phi<- list()
  m_i<- list()
  ## MINIMUM STAGE BASED LENGTHS
  stage$min_adult_length<- min_adult_length
  stage$min_sa_length<- min_sa_length
  ## SURVIVAL PARAMETERS BY STAGE
  stage$phi_1<- phi_1
  stage$phi_juv<- phi_juv
  stage$phi_sa<- phi_sa
  stage$phi_adult<- phi_adult
  ## TRANSITION PARAMETERS
  stage$j2j<- j2j
  stage$sa2sa<- sa2sa
  
  # AGE 1+ SURVIVALS
  ## GROWTH PARAMETERS -- NOT SEX SPECIFIC (EFFECTS MATURATION, WHICH DOES DIFFER 
  ## AMONG SEXES)
  ### SR
  vbgf<-readRDS("./growth fits/output/vbgf-known-and-unknown-age-constant-sigma.Rds")
  growth$vbgf$Linf<- as.numeric(vbgf$fit$BUGSoutput$mean$Linf)
  growth$vbgf$k<- as.numeric(vbgf$fit$BUGSoutput$mean$k)
  growth$vbgf$t0<- as.numeric(vbgf$fit$BUGSoutput$mean$t0) 
  growth$vbgf$sigma<- as.numeric(vbgf$fit$BUGSoutput$mean$sigma)
  ### MEC
  model_fit<-readRDS("./growth fits/vbgf-rpma-4-mec/rpma-4-vbgf-re.Rds")
  growth$ln_vbgf$Linf<- exp(model_fit$report$mu[1])
  growth$ln_vbgf$k<- exp(model_fit$report$mu[2])
  growth$ln_vbgf$t0<- model_fit$report$t0 
  growth$ln_vbgf$sigma<- exp(model_fit$report$log_obs_er)
  # WHY IS THIS LOG OBS ERROR???
  rm(vbgf, model_fit)
  ## FIT WITH RANDOM EFFECT OF YEAR
  # J2SA<-plogis(rnorm(n,-3.036152,sqrt(1/4.010753e-01)))
  # SA2A<-plogis(rnorm(n,2.788886e-03 ,sqrt(1/4.796302e-01)))
  # phi_age1<-plogis(rnorm(n,-1.711568,sqrt(1/8.418292e-01 )))
  # phi_junvenile<-plogis(rnorm(n,1.570630,sqrt(1/4.700724e-01)))
  # phi_subadult<-plogis(rnorm(n,1.498391,sqrt(1/1.450544)))
  # phi_adult<-plogis(rnorm(n,1.633617,sqrt(1/4.068083e+01)))
  ## CONVERT TO SURVIVAL BY AGE
  ### PROBABILITY OF BEING IN EACH STAGE BY AGE
  p_stage<- lapply(c("vbgf", "ln_vbgf"), function(x)
  {
    age<- 3:max_age
    if(x=="vbgf")
    {
      La<- growth$vbgf$Linf*(1-exp(-growth$vbgf$k*(age-growth$vbgf$t0)))
      p_juv<- pnorm(stage$min_sa_length, mean=La, sd=growth$vbgf$sigma)
      p_sa<- pnorm(stage$min_adult_length, mean=La, sd=growth$vbgf$sigma)-p_juv
      p_adult<- 1-pnorm(stage$min_adult_length, mean=La, sd=growth$vbgf$sigma)
    }
    if(x=="ln_vbgf")
    {
      La<- growth$ln_vbgf$Linf*(1-exp(-growth$ln_vbgf$k*(age-growth$ln_vbgf$t0)))
      p_juv<- pnorm(log(stage$min_sa_length), mean=log(La), sd=growth$ln_vbgf$sigma)
      p_sa<- pnorm(log(stage$min_adult_length), mean=log(La), sd=growth$ln_vbgf$sigma)-p_juv
      p_adult<- 1-pnorm(log(stage$min_adult_length), mean=log(La), sd=growth$ln_vbgf$sigma)
    }
    # ADD AGES 1 AND 2 AND ADJUST FOR ANNUAL TRANSITIONS BETWEEN STAGES 
    # (COULD ADD A MINIMUM AGE FOR SUBADULT)
    p_juv<- c(0,1,p_juv)
    p_sa<- c(0,0,p_sa[1]+p_adult[1], p_sa[2:length(p_sa)])
    p_adult<-c(0,0,0,p_adult[2:length(p_adult)])
    p_1<- c(1,rep(0, length(p_juv)-1))
    # ADJUST FOR MINIMUM AGE AT MATURATION
    p_sa[1:(a_min-1)]<- p_sa[1:(a_min-1)] + p_adult[1:(a_min-1)]
    p_adult[1:(a_min-1)]<- 0
    prob<- list(age1=p_1, juv=p_juv, subA=p_sa, adult=p_adult)
    return(prob)
  })
  names(p_stage)<- c("vbgf", "ln_vbgf")
  
  ## SURVIVALS FROM SR VBGF FIT
  phi$vbgf<- p_stage$vbgf$age1*stage$phi_1 + 
    p_stage$vbgf$juv*stage$phi_juv + p_stage$vbgf$subA*stage$phi_sa + 
    p_stage$vbgf$adult*stage$phi_adult
  ## SURVIVALS FROM MEC LN_VBGF FIT
  phi$ln_vbgf<- p_stage$ln_vbgf$age1*stage$phi_1 + 
    p_stage$ln_vbgf$juv*stage$phi_juv + p_stage$ln_vbgf$subA*stage$phi_sa + 
    p_stage$ln_vbgf$adult*stage$phi_adult
  
  
  ### MATURATION FROM SR VBGF FIT
  m<- sapply(2:length(p_stage$vbgf$adult), function(i)
  {
    p_stage$vbgf$adult[i]-(p_stage$vbgf$adult[i-1]*stage$phi_adult)/
      (p_stage$vbgf$age1[i-1]*stage$phi_1 + p_stage$vbgf$juv[i-1]*stage$phi_juv + 
         p_stage$vbgf$subA[i-1]*stage$phi_sa + p_stage$vbgf$adult[i-1]*stage$phi_adult)
  })
  m_i$vbgf<- c(0, m)
  ### MATURATION FROM MEC LN_VBGF FIT
  m<- sapply(2:length(p_stage$ln_vbgf$adult), function(i)
  {
    p_stage$ln_vbgf$adult[i]-(p_stage$ln_vbgf$adult[i-1]*stage$phi_adult)/
      (p_stage$ln_vbgf$age1[i-1]*stage$phi_1 + p_stage$ln_vbgf$juv[i-1]*stage$phi_juv + 
         p_stage$ln_vbgf$subA[i-1]*stage$phi_sa + p_stage$ln_vbgf$adult[i-1]*stage$phi_adult)
  })
  #### SET NEGATIVE PROBABILITIES TO 0
  m[m<0]<-0
  m_i$ln_vbgf<- c(0, m)
  #### NOTE NOT ALL FISH ARE MATURE BY AGE 41; CAN ADJUST THIS MANUALLY IF DESIRED
  return(list(max_age=max_age, a_min=a_min, phi=phi, m_i=m_i, growth=growth, 
              stage=stage))
}

if(exists("run_phi_mat"))
{
  if(run_phi_mat)
  {
    inps<- surv_and_mat_from_stage()  
    saveRDS(inps, "./baseline-parameters/phi_mi.rds")
  }
}
