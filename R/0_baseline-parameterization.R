inputs<- list()
# MAXIMUM AGE
inputs$max_age<- 41

# MINIMUM AGE AT MATURATION
inputs$a_min<- 8

# MINIMUM STAGE BASED LENGTHS
inputs$min_adult_length<- 800
inputs$min_sa_length<- 600

# AGE 1+ SURVIVALS
## GROWTH PARAMETERS -- NOT SEX SPECIFIC (EFFECTS MATURATION, WHICH DOES DIFFER 
## AMONG SEXES)
### SR
vbgf<-readRDS("./growth fits/output/vbgf-known-and-unknown-age-constant-sigma.Rds")
inputs$vbgf$Linf<- as.numeric(vbgf$fit$BUGSoutput$mean$Linf)
inputs$vbgf$k<- as.numeric(vbgf$fit$BUGSoutput$mean$k)
inputs$vbgf$t0<- as.numeric(vbgf$fit$BUGSoutput$mean$t0) 
inputs$vbgf$sigma<- as.numeric(vbgf$fit$BUGSoutput$mean$sigma)
### MEC
model_fit<-readRDS("./growth fits/vbgf-rpma-4-mec/rpma-4-vbgf-re.Rds")
inputs$ln_vbgf$Linf<- exp(model_fit$report$mu[1])
inputs$ln_vbgf$k<- exp(model_fit$report$mu[2])
inputs$ln_vbgf$t0<- model_fit$report$t0 
inputs$ln_vbgf$sigma<- exp(model_fit$report$log_obs_er)
  # WHY IS THIS LOG OBS ERROR???
rm(vbgf, model_fit)
## SURVIVAL PARAMETERS BY STAGE
inputs$phi_1<- 0.151
inputs$phi_juv<- 0.780
inputs$phi_sa<- 0.932
inputs$phi_adult<- 0.976
## TRANSITION PARAMETERS
inputs$j2j<- 0.893
inputs$sa2sa<- 0.897
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
  age<- 3:inputs$max_age
  if(x=="vbgf")
  {
    La<- inputs$vbgf$Linf*(1-exp(-inputs$vbgf$k*(age-inputs$vbgf$t0)))
    p_juv<- pnorm(inputs$min_sa_length, mean=La, sd=inputs$vbgf$sigma)
    p_sa<- pnorm(inputs$min_adult_length, mean=La, sd=inputs$vbgf$sigma)-p_juv
    p_adult<- 1-pnorm(inputs$min_adult_length, mean=La, sd=inputs$vbgf$sigma)
  }
  if(x=="ln_vbgf")
  {
    La<- inputs$ln_vbgf$Linf*(1-exp(-inputs$ln_vbgf$k*(age-inputs$ln_vbgf$t0)))
    p_juv<- pnorm(log(inputs$min_sa_length), mean=log(La), sd=inputs$ln_vbgf$sigma)
    p_sa<- pnorm(log(inputs$min_adult_length), mean=log(La), sd=inputs$ln_vbgf$sigma)-p_juv
    p_adult<- 1-pnorm(log(inputs$min_adult_length), mean=log(La), sd=inputs$ln_vbgf$sigma)
  }
  # ADD AGES 1 AND 2 AND ADJUST FOR ANNUAL TRANSITIONS BETWEEN STAGES 
  # (COULD ADD A MINIMUM AGE FOR SUBADULT)
  p_juv<- c(0,1,p_juv)
  p_sa<- c(0,0,p_sa[1]+p_adult[1], p_sa[2:length(p_sa)])
  p_adult<-c(0,0,0,p_adult[2:length(p_adult)])
  p_1<- c(1,rep(0, length(p_juv)-1))
  # ADJUST FOR MINIMUM AGE AT MATURATION
  p_sa[1:(inputs$a_min-1)]<- p_sa[1:(inputs$a_min-1)] + p_adult[1:(inputs$a_min-1)]
  p_adult[1:(inputs$a_min-1)]<- 0
  prob<- list(age1=p_1, juv=p_juv, subA=p_sa, adult=p_adult)
  return(prob)
})
names(p_stage)<- c("vbgf", "ln_vbgf")

## SURVIVALS FROM SR VBGF FIT
inputs$phi$vbgf<- p_stage$vbgf$age1*inputs$phi_1 + 
  p_stage$vbgf$juv*inputs$phi_juv + p_stage$vbgf$subA*inputs$phi_sa + 
  p_stage$vbgf$adult*inputs$phi_adult
## SURVIVALS FROM MEC LN_VBGF FIT
inputs$phi$ln_vbgf<- p_stage$ln_vbgf$age1*inputs$phi_1 + 
  p_stage$ln_vbgf$juv*inputs$phi_juv + p_stage$ln_vbgf$subA*inputs$phi_sa + 
  p_stage$ln_vbgf$adult*inputs$phi_adult


# FERTILITY VALUES
## PROPORTION OF FEMALES THAT ARE REPRODUCTIVELY READY
### MATURATION FROM SR VBGF FIT
m<- sapply(2:length(p_stage$vbgf$adult), function(i)
  {
    p_stage$vbgf$adult[i]-(p_stage$vbgf$adult[i-1]*inputs$phi_adult)/
      (p_stage$vbgf$age1[i-1]*inputs$phi_1 + p_stage$vbgf$juv[i-1]*inputs$phi_juv + 
        p_stage$vbgf$subA[i-1]*inputs$phi_sa + p_stage$vbgf$adult[i-1]*inputs$phi_adult)
  })
inputs$m_i$vbgf<- c(0, m)
### MATURATION FROM MEC LN_VBGF FIT
m<- sapply(2:length(p_stage$ln_vbgf$adult), function(i)
{
  p_stage$ln_vbgf$adult[i]-(p_stage$ln_vbgf$adult[i-1]*inputs$phi_adult)/
    (p_stage$ln_vbgf$age1[i-1]*inputs$phi_1 + p_stage$ln_vbgf$juv[i-1]*inputs$phi_juv + 
       p_stage$ln_vbgf$subA[i-1]*inputs$phi_sa + p_stage$ln_vbgf$adult[i-1]*inputs$phi_adult)
})
#### SET NEGATIVE PROBABILITIES TO 0
m[m<0]<-0
inputs$m_i$ln_vbgf<- c(0, m)
#### NOTE NOT ALL FISH ARE MATURE BY AGE 41; CAN ADJUST THIS MANUALLY IF DESIRED
### MATURATION FROM VALUES PRESENTED AT 2019 FSM
a_m<- 8
a_max<- 15
a_h<- 10
k<- 1.47 #means 3*sigma is approximately 3.7 years
mat_Cdist<-1/(1+exp(-k*(a_m:(a_max-1)-a_h)))
m<-rep(0, a_max-a_m+1)
m[1]<-mat_Cdist[1]
for(i in 2:length(mat_Cdist))
{
  m[i]<- mat_Cdist[i]-mat_Cdist[i-1]
}
m[length(mat_Cdist)+1]<- 1-mat_Cdist[length(mat_Cdist)]
inputs$mat$a_min<- a_m
inputs$mat$a_max<- a_max
inputs$mat$a_h<- a_h
inputs$mat$k<- k
inputs$m_i$FSM<- c(rep(0,a_m-1), m, rep(0, inputs$max_age-a_max))
rm(a_max, a_h, k, mat_Cdist, m, i, a_m)
### REPRODUCTIVE PERIOD
max_period<- 5
#### FROM FULLER ET AL. 2007 AND DELONAY ET AL. 2016:
#### SPECIFIC PERIOD 1-3 PROBS FROM KNOWN DATA WITH THE REMAINING 
#### PROBABILITY SPLIT AMONG PERIODS 4 TO MAX_PERIOD SUCH THAT
#### PROBS[T+1]=PROBS[T]/2 FOR T>4
probs<-c(0, 8/21, 13/21*3/5, 
         13/21*2/5*1/sum(2^(0:(max_period-4)))*2^((max_period-4):0))
#### SAVE REPRODUCTION PERIOD INPUTS
inputs$spawning_period$max_period<- max_period
inputs$spawning_period$p_t<- c(probs, rep(0, inputs$max_age-max_period))
rm(max_period, probs)
### CALCULATE PROPORTION FROM MATURATION PROBABILITIES 
### AND REPRODUCTION PERIODS
#### FROM SR VBGF FIT
inputs$psi$vbgf<- rep(0, inputs$max_age)
inputs$psi$vbgf[inputs$a_min]<- inputs$m_i$vbgf[inputs$a_min]
for(i in (inputs$a_min+1):inputs$max_age)
{
  inputs$psi$vbgf[i]<- inputs$m_i$vbgf[i]+
    sum(inputs$psi$vbgf[inputs$a_min:(i-1)]*inputs$spawning_period$p_t[i-inputs$a_min:(i-1)])
}
rm(i)
#### FROM MEC LN_VBGF FIT
inputs$psi$ln_vbgf<- rep(0, inputs$max_age)
inputs$psi$ln_vbgf[inputs$a_min]<- inputs$m_i$ln_vbgf[inputs$a_min]
for(i in (inputs$a_min+1):inputs$max_age)
{
  inputs$psi$ln_vbgf[i]<- inputs$m_i$ln_vbgf[i]+
    sum(inputs$psi$ln_vbgf[inputs$a_min:(i-1)]*inputs$spawning_period$p_t[i-inputs$a_min:(i-1)])
}
rm(i)
#### FROM UMOR APPROACH AND 2019 FSM VALUES
inputs$psi$FSM<- rep(0, inputs$max_age)
inputs$psi$FSM[inputs$mat$a_min]<- inputs$m_i$FSM[inputs$mat$a_min]
for(i in (inputs$mat$a_min+1):inputs$max_age)
{
  inputs$psi$FSM[i]<- inputs$m_i$FSM[i]+
    sum(inputs$psi$FSM[inputs$mat$a_min:(i-1)]*inputs$spawning_period$p_t[i-inputs$mat$a_min:(i-1)])
}
rm(i)

## SEX RATIO (PROBABILITY OF BEING FEMALE)
inputs$probF<- 0.5

## NUMBER OF EGGS PER SPAWNING FEMALE
# ### AGE-LENGTH RELATIONSHIP
# length_at_age<- function(age=NULL,
#                          reps=NULL,
#                          inputs=NULL,
#                          type="ln_vbgf")
# {
#   if(type=="vbgf")
#   {
#     La<- inputs$vbgf$Linf*(1-exp(-inputs$vbgf$k*(age-inputs$vbgf$t0)))
#     l<-rnorm(reps, La, inputs$vbgf$sigma)
#   }
#   if(type=="ln_vbgf")
#   {
#     La<- inputs$ln_vbgf$Linf*(1-exp(-inputs$ln_vbgf$k*(age-inputs$ln_vbgf$t0)))
#     l<- exp(rnorm(reps, log(La), inputs$ln_vbgf$sigma))
#   }
#   for(i in 1:length(l))
#   {
#     l[i]<- ifelse(l[i]<0, La, l[i])
#     l[i]<- ifelse(l[i]>1800, 1800, l[i])
#     
#   }
#   return(l)
# }
# 
# ### LENGTH-FECUNDITY RELATIONSHIP
# dat<-read.csv("../fecundity/_dat/Fecundity.csv")
# #### FILL MISSING VALUES AS NA
# dat[dat==-99]<-NA
# #### ADD AN INDVIDUAL ID
# dat$id<-1:nrow(dat)
# #### MEAN AND SD TO SCALE VARIABLES
# mean_fl<-mean(na.omit(dat$FL))
# sd_fl<-sd(na.omit(dat$FL))
# #### SCALE DATA 
# dat$fl_std<-scale(dat$FL, center = mean_fl, scale = sd_fl)
# #### FIT MODEL GIVEN DATA (FIT 7 USED IN POP MODEL)
# #### VERIFY ALL ASSUMPTIONS MET
# library(lme4)
# library(merTools)
# fit<-glmer(EGGS~fl_std+(1|id),dat[which(dat$BASIN=="Lower"),],
#            family=poisson(link="log"))
# intrcpt<- unname(fixef(fit)[1])
# slp<- unname(fixef(fit)[2])
# disp<- (sd(unlist(ranef(fit)$id))*(length(unlist(ranef(fit)$id))-1)/length(unlist(ranef(fit)$id))
#         +sd(unlist(ranef(fit)$id)))/2
# # plot(dat$FL, dat$EGGS)
# # points(800:1700, exp(intrcpt+slp*scale(800:1700, center = mean_fl, scale = sd_fl)), 
# #        col="blue", type="l")
# 
# eggs<- function(fork_length=NULL,
#                 a=NULL,
#                 b=NULL,
#                 dispersion_param=NULL,
#                 min_length=NULL)
# {
#   N<-length(fork_length)
#   fl_normalized<- (fork_length - mean_fl)/sd_fl
#   egg_no<- rpois(N,exp(rnorm(N,a + b*fl_normalized,dispersion_param)))
#   if(any(fork_length<inputs$min_adult_length))
#   {
#     egg_no[which(fork_length<inputs$min_adult_length)]<- 0
#   }
#   return(egg_no)
# }
# 
# ## SIMULATE FECUNDITY BY AGE
# n=1000000
# a=1:inputs$max_age
# fec_vbgf<- lapply(a, function(x)
# {
#   lengths<- length_at_age(age=x, reps=n, inputs = inputs, type="vbgf")
#   mn_lgth<- mean(lengths)
#   med_lgth<- median(lengths)
#   fecundity<- eggs(fork_length=lengths, a=intrcpt, b=slp, dispersion_param=disp,
#                    min_length = inputs$min_adult_length)
#   mn_eggs<- mean(fecundity)
#   med_eggs<- median(fecundity) 
#   if(any(lengths<inputs$min_adult_length) & any(lengths>=inputs$min_adult_length))
#   {
#     tmp<- which(lengths>=inputs$min_adult_length)
#     N_adult<- length(tmp)
#     mn_adult_lgth<- mean(lengths[tmp])
#     med_adult_lgth<- median(lengths[tmp])
#     mn_adult_eggs<- mean(fecundity[tmp])
#     med_adult_eggs<- median(fecundity[tmp])
#     min_adult_eggs<- min(fecundity[tmp])
#   }
#   if(all(lengths>=inputs$min_adult_length))
#   {
#     N_adult<- n
#     mn_adult_lgth<- mn_lgth
#     med_adult_lgth<- med_lgth
#     mn_adult_eggs<- mn_eggs
#     med_adult_eggs<- med_eggs
#     min_adult_eggs<- min(fecundity)
#   }
#   if(all(lengths<inputs$min_adult_length))
#   {
#     N_adult<- 0
#     mn_adult_lgth<- NA
#     med_adult_lgth<- NA
#     mn_adult_eggs<- 0
#     med_adult_eggs<- 0
#     min_adult_eggs<- 0
#   }
#   return(list(mean_eggs=mn_eggs, median_eggs=med_eggs,
#               mean_adult_eggs=mn_adult_eggs, median_adult_eggs=med_adult_eggs,
#               max_eggs=max(fecundity), min_eggs=min(fecundity),
#               min_adult_eggs=min_adult_eggs,
#               mean_length=mn_lgth, median_length=med_lgth,
#               mean_adult_length=mn_adult_lgth, median_adult_length=med_adult_lgth,
#               prop_adults=N_adult/n))
# })
# fecundity_vbgf5<- data.frame(Age=a, 
#                             Mean_Eggs_All=sapply(fec_vbgf, "[[", 1),
#                             Median_Eggs_All=sapply(fec_vbgf, "[[", 2),
#                             Mean_Eggs_Adults=sapply(fec_vbgf, "[[", 3),
#                             Median_Eggs_Adults=sapply(fec_vbgf, "[[", 4),
#                             Max_Eggs_Simulated=sapply(fec_vbgf, "[[", 5),
#                             Min_Eggs_Simulated=sapply(fec_vbgf, "[[", 6),
#                             Min_Eggs_Simulated_Adults=sapply(fec_vbgf, "[[", 7),
#                             Mean_Length_All=sapply(fec_vbgf, "[[", 8),
#                             Median_Length_All=sapply(fec_vbgf, "[[", 9),
#                             Mean_Length_Adult=sapply(fec_vbgf, "[[", 10),
#                             Median_Length_Adult=sapply(fec_vbgf, "[[", 11),
#                             Proportion_Adults=sapply(fec_vbgf, "[[", 12))
# write.csv(fecundity_vbgf5, "./baseline-parameters/fecundity_estimates_by_age_vbgf_5mil.csv",
#           row.names = FALSE)
# 
# fec_ln<- lapply(a, function(x)
# {
#   lengths<- length_at_age(age=x, reps=n, inputs=inputs, type="ln_vbgf")
#   mn_lgth<- mean(lengths)
#   med_lgth<- median(lengths)
#   fecundity<- eggs(fork_length=lengths, a=intrcpt, b=slp, dispersion_param=disp,
#                    min_length = inputs$min_adult_length)
#   mn_eggs<- mean(fecundity)
#   med_eggs<- median(fecundity)
#   if(any(lengths<inputs$min_adult_length) & any(lengths>=inputs$min_adult_length))
#   {
#     tmp<- which(lengths>=inputs$min_adult_length)
#     N_adult<- length(tmp)
#     mn_adult_lgth<- mean(lengths[tmp])
#     med_adult_lgth<- median(lengths[tmp])
#     mn_adult_eggs<- mean(fecundity[tmp])
#     med_adult_eggs<- median(fecundity[tmp])
#     min_adult_eggs<- min(fecundity[tmp])
#   }
#   if(all(lengths>=inputs$min_adult_length))
#   {
#     N_adult<- n
#     mn_adult_lgth<- mn_lgth
#     med_adult_lgth<- med_lgth
#     mn_adult_eggs<- mn_eggs
#     med_adult_eggs<- med_eggs
#     min_adult_eggs<- min(fecundity)
#   }
#   if(all(lengths<inputs$min_adult_length))
#   {
#     N_adult<- 0
#     mn_adult_lgth<- NA
#     med_adult_lgth<- NA
#     mn_adult_eggs<- 0
#     med_adult_eggs<- 0
#     min_adult_eggs<- 0
#   }
#   return(list(mean_eggs=mn_eggs, median_eggs=med_eggs,
#               mean_adult_eggs=mn_adult_eggs, median_adult_eggs=med_adult_eggs,
#               max_eggs=max(fecundity), min_eggs=min(fecundity),
#               min_adult_eggs=min_adult_eggs,
#               mean_length=mn_lgth, median_length=med_lgth,
#               mean_adult_length=mn_adult_lgth, median_adult_length=med_adult_lgth,
#               prop_adults=N_adult/n))
# })
# fecundity_ln<- data.frame(Age=a, 
#                           Mean_Eggs_All=sapply(fec_ln, "[[", 1),
#                           Median_Eggs_All=sapply(fec_ln, "[[", 2),
#                           Mean_Eggs_Adults=sapply(fec_ln, "[[", 3),
#                           Median_Eggs_Adults=sapply(fec_ln, "[[", 4),
#                           Max_Eggs_Simulated=sapply(fec_ln, "[[", 5),
#                           Min_Eggs_Simulated=sapply(fec_ln, "[[", 6),
#                           Min_Eggs_Simulated_Adults=sapply(fec_ln, "[[", 7),
#                           Mean_Length_All=sapply(fec_ln, "[[", 8),
#                           Median_Length_All=sapply(fec_ln, "[[", 9),
#                           Mean_Length_Adult=sapply(fec_ln, "[[", 10),
#                           Median_Length_Adult=sapply(fec_ln, "[[", 11),
#                           Proportion_Adults=sapply(fec_ln, "[[", 12))
# write.csv(fecundity_ln, "./baseline-parameters/fecundity_estimates_by_age_ln_vbgf.csv",
#           row.names = FALSE)
# 
# ## FIGURE COMPARISON OF APPROACHES
# age<-1:41
# ### VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ADULTS ONLY
# plot(age, fecundity_vbgf$Mean_Eggs_Adults, xlab="Age", ylab="Mean Number of Eggs", col="blue", pch=19)
# points(age, fecundity_vbgf2$Mean_Eggs_Adults, pch=19)
# points(age, fecundity_vbgf5$Mean_Eggs_Adults, pch=19, col="red")
# legend("bottomright", c("1 million", "2 million", "5 million"), 
#        col=c("black", "blue", "red"), pch=19, bty="n")
# 
# ### VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ALL BY AGE: ADULTS WEIGHTED
# plot(age, fecundity_vbgf$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
# points(age, fecundity_vbgf2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19)
# points(age, fecundity_vbgf5$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="red")
# legend("bottomright", c("1 million", "2 million", "5 million"), 
#        col=c("black", "blue", "red"), pch=19, bty="n")
# 
# ### LN_VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ADULTS ONLY
# plot(age, fecundity_ln$Mean_Eggs_Adults, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
# points(age, fecundity_ln2$Mean_Eggs_Adults, pch=19)
# points(age, fecundity_ln5$Mean_Eggs_Adults, pch=19, col="red")
# legend("bottomright", c("1 million", "2 million", "5 million"), 
#        col=c("black", "blue", "red"), pch=19, bty="n")
# 
# ### LN_VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ALL BY AGE: ADULTS WEIGHTED
# plot(age, fecundity_ln$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
# points(age, fecundity_ln2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19)
# points(age, fecundity_ln5$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="red")
# legend("bottomright", c("1 million", "2 million", "5 million"), 
#        col=c("black", "blue", "red"), pch=19, bty="n")
# 
# ### VBGF EGGS FROM ALL BY AGE: AVERAGE OF DRAWS VS. ADULTS WEIGHTED
# plot(age, fecundity_vbgf2$Mean_Eggs_All, xlab="Age", ylab="Number of Eggs", pch=19)
# points(age, fecundity_vbgf2$Mean_Eggs_Adults*p_stage$vbgf$adult, pch=19, col="blue")
# legend("bottomright", c("Eggs from All Fish", "Eggs from Adults Weighted by Proportion"), 
#        col=c("black", "blue"), pch=19, bty="n")
# 
# ### LN_VBGF EGGS FROM ALL BY AGE: AVERAGE OF DRAWS VS. ADULTS WEIGHTED
# plot(age, fecundity_ln2$Mean_Eggs_All, xlab="Age", ylab="Number of Eggs", pch=19)
# points(age, fecundity_ln2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="blue")
# legend("bottomright", c("Eggs from All Fish", "Eggs from Adults Weighted by Proportion"), 
#        col=c("black", "blue"), pch=19, bty="n")
# 
# rm(dat, mean_fl, sd_fl, fit, intrcpt, slp, disp, n, a, fec_vbgf, fecundity_vbgf,
#    fec_ln, fecundity_ln, length_at_age, eggs)

E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_vbgf_2mil.csv")
inputs$eggs$vbgf<- E$Mean_Eggs_Adults
inputs$eggs$vbgf[1:(inputs$a_min-1)]<- 0
inputs$allF_eggs$vbgf<- E$Mean_Eggs_All
inputs$allF_eggs$vbgf[1:(inputs$a_min-1)]<- 0
E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_ln_vbgf_2mil.csv")
inputs$eggs$ln_vbgf<- E$Mean_Eggs_Adults
inputs$eggs$ln_vbgf[1:(inputs$a_min-1)]<- 0
inputs$allF_eggs$ln_vbgf<- E$Mean_Eggs_All
inputs$allF_eggs$ln_vbgf[1:(inputs$a_min-1)]<- 0

rm(E, p_stage)

### FIGURES
#### SURVIVAL
par(mar=c(4,4,4,2),
    oma=c(1,1,0,0))
plot(2:40, inputs$phi$vbgf[2:40], pch=19, xlab="Age", ylab="Survival Probability")
points(inputs$phi$ln_vbgf[2:40], col="blue", pch=19)
legend("bottomright", c("JAGS VBGF", "TMB LN(VBGF)"), col=c("black", "blue"), 
       pch=19, bty="n")
#### AVERAGE LENGTH AT AGE
age<-1:41
La_vbgf<- inputs$vbgf$Linf*(1-exp(-inputs$vbgf$k*(age-inputs$vbgf$t0)))
La_ln<- inputs$ln_vbgf$Linf*(1-exp(-inputs$ln_vbgf$k*(age-inputs$ln_vbgf$t0)))
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
