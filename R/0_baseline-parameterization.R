
# MAXIMUM AGE, MINIMUM AGE AT MATURATION, AGE 1+ SURVIVALS, 
# AND MATURATION AGE PROBABILITIES FROM GROWTH
params<- readRDS("./baseline-parameters/phi_mi_fec.rds")

# FERTILITY VALUES
## PROPORTION OF FEMALES THAT ARE REPRODUCTIVELY READY
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
params$mat$a_min<- a_m
params$mat$a_max<- a_max
params$mat$a_h<- a_h
params$mat$k<- k
params$m_i$FSM<- c(rep(0,a_m-1), m, rep(0, params$max_age-a_max))
rm(a_max, a_h, k, mat_Cdist, m, i, a_m)
### REPRODUCTIVE PERIOD
max_period<- 5
#### FROM FULLER ET AL. 2007 AND DELONAY ET AL. 2016:
#### SPECIFIC PERIOD 1-3 PROBS FROM KNOWN DATA WITH THE REMAINING 
#### PROBABILITY SPLIT AMONG PERIODS 4 TO MAX_PERIOD SUCH THAT
#### PROBS[T+1]=PROBS[T]/2 FOR T>4
probs<-c(0, 8/21, 13/21*3/5, 
         13/21*2/5*1/sum(2^(0:(max_period-4)))*2^((max_period-4):0))
#### SAVE REPRODUCTION PERIOD PARAMETERS
params$spawning_period$max_period<- max_period
params$spawning_period$p_t<- c(probs, rep(0, params$max_age-max_period))
rm(max_period, probs)
### CALCULATE PROPORTION FROM MATURATION PROBABILITIES 
### AND REPRODUCTION PERIODS
#### FROM SR VBGF FIT
params$psi$vbgf<- rep(0, params$max_age)
params$psi$vbgf[params$a_min]<- params$m_i$vbgf[params$a_min]
for(i in (params$a_min+1):params$max_age)
{
  params$psi$vbgf[i]<- params$m_i$vbgf[i]+
    sum(params$psi$vbgf[params$a_min:(i-1)]*params$spawning_period$p_t[i-params$a_min:(i-1)])
}
rm(i)
#### FROM MEC LN_VBGF FIT
params$psi$ln_vbgf<- rep(0, params$max_age)
params$psi$ln_vbgf[params$a_min]<- params$m_i$ln_vbgf[params$a_min]
for(i in (params$a_min+1):params$max_age)
{
  params$psi$ln_vbgf[i]<- params$m_i$ln_vbgf[i]+
    sum(params$psi$ln_vbgf[params$a_min:(i-1)]*params$spawning_period$p_t[i-params$a_min:(i-1)])
}
rm(i)
#### FROM UMOR APPROACH AND 2019 FSM VALUES
params$psi$FSM<- rep(0, params$max_age)
params$psi$FSM[params$mat$a_min]<- params$m_i$FSM[params$mat$a_min]
for(i in (params$mat$a_min+1):params$max_age)
{
  params$psi$FSM[i]<- params$m_i$FSM[i]+
    sum(params$psi$FSM[params$mat$a_min:(i-1)]*params$spawning_period$p_t[i-params$mat$a_min:(i-1)])
}
rm(i)

## SEX RATIO (PROBABILITY OF BEING FEMALE)
params$probF<- 0.5

## AGE-SPECIFIC FECUNDITY
E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_vbgf_2mil.csv")
params$eggs$vbgf<- E$Mean_Eggs_Adults
# params$eggs$vbgf[1:(params$a_min-1)]<- 0
params$allF_eggs$vbgf<- E$Mean_Eggs_All
# params$allF_eggs$vbgf[1:(params$a_min-1)]<- 0
E<- read.csv("./baseline-parameters/fecundity_estimates_by_age_ln_vbgf_2mil.csv")
params$eggs$ln_vbgf<- E$Mean_Eggs_Adults
# params$eggs$ln_vbgf[1:(params$a_min-1)]<- 0
params$allF_eggs$ln_vbgf<- E$Mean_Eggs_All
# params$allF_eggs$ln_vbgf[1:(params$a_min-1)]<- 0

rm(E)

## AGE-0 SURVIVAL
params$phi0<- 0.000075

saveRDS(params, "./output/parameterizations/default-parameters.rds")
