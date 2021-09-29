inputs<- list()
# MAXIMUM AGE
inputs$max_age<- 41

# AGE 1+ SURVIVALS
## GROWTH PARAMETERS -- NOT SEX SPECIFIC (EFFECTS MATURATION, WHICH DOES DIFFER 
## AMONG SEXES)
### SR
vbgf<-readRDS("./growth fits/output/vbgf-known-and-unknown-age-constant-sigma.Rds")
inputs$vbgf$Linf<- vbgf$fit$BUGSoutput$mean$Linf
inputs$vbgf$k<- vbgf$fit$BUGSoutput$mean$k
inputs$vbgf$t0<- vbgf$fit$BUGSoutput$mean$t0 
inputs$vbgf$sigma<- vbgf$fit$BUGSoutput$mean$sigma
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
    p_juv<- pnorm(600, mean=La, sd=inputs$vbgf$sigma)
    p_sa<- pnorm(800, mean=La, sd=inputs$vbgf$sigma)-p_juv
    p_adult<- 1-pnorm(800, mean=La, sd=inputs$vbgf$sigma)
  }
  if(x=="ln_vbgf")
  {
    La<- inputs$ln_vbgf$Linf*(1-exp(-inputs$ln_vbgf$k*(age-inputs$ln_vbgf$t0)))
    p_juv<- pnorm(log(600), mean=log(La), sd=inputs$ln_vbgf$sigma)
    p_sa<- pnorm(log(800), mean=log(La), sd=inputs$ln_vbgf$sigma)-p_juv
    p_adult<- 1-pnorm(log(800), mean=log(La), sd=inputs$ln_vbgf$sigma)
  }
  p_juv<- c(0,1,p_juv)
  p_sa<- c(0,0,p_sa[1]+p_adult[1], p_sa[2:length(p_sa)])
  p_adult<-c(0,0,0,p_adult[2:length(p_adult)])
  p_1<- c(1,rep(0, length(p_juv)-1))
  prob<- list(age1=p_1, juv=p_juv, subA=p_sa, adult=p_adult)
  return(prob)
})
names(p_stage)<- c("vbgf", "ln_vbgf")
### MINIMUM AGE AT MATURATION
### PROBABLY WANT TO ADJUST JUST GOING OFF OF THE DISTRIBUTIONS ABOVE
a_min<- 4
inputs$a_min<- a_min
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
inputs$reproduction$max_period<- max_period
inputs$reproduction$p_t<- c(probs, rep(0, inputs$max_age-max_period))
rm(max_period, probs)
### CALCULATE PROPORTION FROM MATURATION PROBABILITIES 
### AND REPRODUCTION PERIODS
#### FROM SR VBGF FIT
inputs$psi$vbgf<- rep(0, inputs$max_age)
inputs$psi$vbgf[a_min]<- inputs$m_i$vbgf[a_min]
for(i in (a_min+1):inputs$max_age)
{
  inputs$psi$vbgf[i]<- inputs$m_i$vbgf[i]+
    sum(inputs$psi$vbgf[a_min:(i-1)]*inputs$reproduction$p_t[i-a_min:(i-1)])
}
rm(i)
#### FROM MEC LN_VBGF FIT
inputs$psi$ln_vbgf<- rep(0, inputs$max_age)
inputs$psi$ln_vbgf[a_min]<- inputs$m_i$ln_vbgf[a_min]
for(i in (a_min+1):inputs$max_age)
{
  inputs$psi$ln_vbgf[i]<- inputs$m_i$ln_vbgf[i]+
    sum(inputs$psi$ln_vbgf[a_min:(i-1)]*inputs$reproduction$p_t[i-a_min:(i-1)])
}
rm(i)
#### FROM UMOR APPROACH AND 2019 FSM VALUES
inputs$psi$FSM<- rep(0, inputs$max_age)
inputs$psi$FSM[inputs$mat$a_min]<- inputs$m_i$FSM[inputs$mat$a_min]
for(i in (inputs$mat$a_min+1):inputs$max_age)
{
  inputs$psi$FSM[i]<- inputs$m_i$FSM[i]+
    sum(inputs$psi$FSM[inputs$mat$a_min:(i-1)]*inputs$reproduction$p_t[i-inputs$mat$a_min:(i-1)])
}
rm(i)

## SEX RATIO (PROBABILITY OF BEING FEMALE)
inputs$probF<- 0.5

## NUMBER OF EGGS PER SPAWNING FEMALE
### AGE-LENGTH RELATIONSHIP
length_at_age<- function(age=NULL,
                         reps=NULL,
                         inputs=NULL,
                         type="ln_vbgf")
{
  if(type=="vbgf")
  {
    La<- inputs$vbgf$Linf*(1-exp(-inputs$vbgf$k*(age-inputs$vbgf$t0)))
    l<-rnorm(reps, La, inputs$vbgf$sigma)
  }
  if(type=="ln_vbgf")
  {
    La<- inputs$ln_vbgf$Linf*(1-exp(-inputs$ln_vbgf$k*(age-inputs$ln_vbgf$t0)))
    l<- exp(rnorm(reps, log(La), inputs$ln_vbgf$sigma))
  }
  for(i in 1:length(l))
  {
    l[i]<- ifelse(l[i]<0, La, l[i])
    l[i]<- ifelse(l[i]>1800, 1800, l[i])
    
  }
  return(l)
}

### LENGTH-FECUNDITY RELATIONSHIP
dat<-read.csv("../fecundity/_dat/Fecundity.csv")
#### FILL MISSING VALUES AS NA
dat[dat==-99]<-NA
#### ADD AN INDVIDUAL ID
dat$id<-1:nrow(dat)
#### MEAN AND SD TO SCALE VARIABLES
mean_fl<-mean(na.omit(dat$FL))
sd_fl<-sd(na.omit(dat$FL))
#### SCALE DATA 
dat$fl_std<-scale(dat$FL, center = mean_fl, scale = sd_fl)
#### FIT MODEL GIVEN DATA (FIT 7 USED IN POP MODEL)
#### VERIFY ALL ASSUMPTIONS MET
library(lme4)
library(merTools)
fit<-glmer(EGGS~fl_std+(1|id),dat[which(dat$BASIN=="Lower"),],
           family=poisson(link="log"))
intrcpt<- unname(fixef(fit)[1])
slp<- unname(fixef(fit)[2])
disp<- (sd(unlist(ranef(fit)$id))*(length(unlist(ranef(fit)$id))-1)/length(unlist(ranef(fit)$id))
        +sd(unlist(ranef(fit)$id)))/2
# plot(dat$FL, dat$EGGS)
# points(800:1700, exp(intrcpt+slp*scale(800:1700, center = mean_fl, scale = sd_fl)), 
#        col="blue", type="l")

eggs<- function(fork_length=NULL,
                a=NULL,
                b=NULL,
                dispersion_param=NULL)
{
  N<-length(fork_length)
  fl_normalized<- (fork_length - mean_fl)/sd_fl
  egg_no<- rpois(N,exp(rnorm(N,a + b*fl_normalized,dispersion_param)))
  return(egg_no)
}

## SIMULATE FECUNDITY BY AGE
n=1000000
a=1:inputs$max_age
fec_vbgf<- lapply(a, function(x)
{
  lengths<- length_at_age(age=x, reps=n, inputs = inputs, type="vbgf")
  mn_lgth<- mean(lengths)
  med_lgth<- median(lengths)
  fecundity<- eggs(fork_length=lengths, a=intrcpt, b=slp, dispersion_param=disp)
  mn_eggs<- mean(fecundity)
  med_eggs<- median(fecundity) 
  return(list(mean_eggs=mn_eggs, median_eggs=med_eggs, 
              max_eggs=max(fecundity), min_eggs=min(fecundity),
              mean_length=mn_lgth, median_length=med_lgth))
})
fecundity_vbgf<- data.frame(Age=a, 
                            Mean_Eggs_Produced=sapply(fec_vbgf, "[[", 1),
                            Median_Eggs_Produced=sapply(fec_vbgf, "[[", 2),
                            Max_Eggs_Simulated=sapply(fec_vbgf, "[[", 3),
                            Min_Eggs_Simulated=sapply(fec_vbgf, "[[", 4),
                            Mean_Length=sapply(fec_vbgf, "[[", 5),
                            Median_Length=sapply(fec_vbgf, "[[", 6))
write.csv(fecundity_vbgf, "./baseline-parameters/fecundity_estimates_by_age_vbgf.csv",
          row.names = FALSE)

fec_ln<- lapply(a, function(x)
{
  lengths<- length_at_age(age=x, reps=n, inputs=inputs, type="ln_vbgf")
  mn_lgth<- mean(lengths)
  med_lgth<- median(lengths)
  fecundity<- eggs(fork_length=lengths, a=intrcpt, b=slp, dispersion_param=disp)
  mn_eggs<- mean(fecundity)
  med_eggs<- median(fecundity) 
  return(list(mean_eggs=mn_eggs, median_eggs=med_eggs, 
              max_eggs=max(fecundity), min_eggs=min(fecundity),
              mean_length=mn_lgth, median_length=med_lgth))
})
fecundity_ln<- data.frame(Age=a, 
                          Mean_Eggs_Produced=sapply(fec_ln, "[[", 1),
                          Median_Eggs_Produced=sapply(fec_ln, "[[", 2),
                          Max_Eggs_Simulated=sapply(fec_ln, "[[", 3),
                          Min_Eggs_Simulated=sapply(fec_ln, "[[", 4),
                          Mean_Length=sapply(fec_ln, "[[", 5),
                          Median_Length=sapply(fec_ln, "[[", 6))
write.csv(fecundity_ln, "./baseline-parameters/fecundity_estimates_by_age_ln_vbgf.csv",
          row.names = FALSE)
rm(dat, mean_fl, sd_fl, fit, intrcpt, slp, disp, n, a, fec_vbgf, fecundity_vbgf,
   fec_ln, fecundity_ln, length_at_age, eggs)


E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_vbgf.csv")
inputs$eggs$vbgf<- E$Mean_Eggs_Produced
inputs$eggs$vbgf[1:(a_min-1)]<- 0
E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_ln_vbgf.csv")
inputs$eggs$ln_vbgf<- E$Mean_Eggs_Produced
inputs$eggs$ln_vbgf[1:(a_min-1)]<- 0

rm(E, a_min, p_stage)
