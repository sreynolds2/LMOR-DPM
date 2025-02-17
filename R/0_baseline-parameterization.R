
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
params$m_i$FSM<- c(rep(0,a_m-1), m, rep(0, params$max_age+1-a_max))
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
params$psi$vbgf<- rep(0, params$max_age+1)
params$psi$vbgf[params$a_min]<- params$m_i$vbgf[params$a_min]
for(i in (params$a_min+1):(params$max_age+1))
{
  params$psi$vbgf[i]<- params$m_i$vbgf[i]+
    sum(params$psi$vbgf[params$a_min:(i-1)]*params$spawning_period$p_t[i-params$a_min:(i-1)])
}
rm(i)
#### FROM MEC LN_VBGF FIT
params$psi$ln_vbgf<- rep(0, params$max_age+1)
params$psi$ln_vbgf[params$a_min]<- params$m_i$ln_vbgf[params$a_min]
for(i in (params$a_min+1):(params$max_age+1))
{
  params$psi$ln_vbgf[i]<- params$m_i$ln_vbgf[i]+
    sum(params$psi$ln_vbgf[params$a_min:(i-1)]*params$spawning_period$p_t[i-params$a_min:(i-1)])
}
rm(i)
#### FROM UMOR APPROACH AND 2019 FSM VALUES
params$psi$FSM<- rep(0, params$max_age+1)
params$psi$FSM[params$mat$a_min]<- params$m_i$FSM[params$mat$a_min]
for(i in (params$mat$a_min+1):(params$max_age+1))
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

## SPAWNING PROBABILITY
### 1 - ATRESIA RATE
params$gamma<- 1

saveRDS(params, "./baseline-parameters/default-parameters.rds")
saveRDS(params, "./output/parameterizations/default-parameters.rds")



# EXTEND BASELINE PARAMETERS TO MAX_AGE UP TO 65
## GENERATE ADDITIONAL PARAMETERS FOR TESTING MAX_AGES 
## UP TO 65
source("./R/2_functions.r")

ext<- create_params(max_age=65,
                         growth_fit = list(Linf=1051.872, k=0.1008206, t0=-1.830915, sigma=47.56881),
                         growth_fit_type = "vbgf",
                         input_id="ext_vbgf")
ext_ln<- create_params(max_age = 65,
                       input_id = "ext_ln")

ext$phi$ln_vbgf<- ext_ln$phi$ln_vbgf
ext$m_i$ln_vbgf<- ext_ln$m_i$ln_vbgf
ext$growth$ln_vbgf<- ext_ln$growth$ln_vbgf
ext$psi$ln_vbgf<- ext_ln$psi$ln_vbgf
ext$eggs$ln_vbgf<- ext_ln$eggs$ln_vbgf
rm(ext_ln)

ext$mat<- params$mat 
mat_Cdist<-1/(1+exp(-ext$mat$k*(ext$mat$a_min:(ext$mat$a_max-1)-ext$mat$a_h)))
m<-rep(0, ext$mat$a_max-ext$mat$a_min+1)
m[1]<-mat_Cdist[1]
for(i in 2:length(mat_Cdist))
{
  m[i]<- mat_Cdist[i]-mat_Cdist[i-1]
}
m[length(mat_Cdist)+1]<- 1-mat_Cdist[length(mat_Cdist)]
ext$m_i$FSM<- c(rep(0,ext$mat$a_min-1), m, rep(0, ext$max_age-ext$mat$a_max))
rm(i, m, mat_Cdist) 
ext$psi$FSM<- rep(0, ext$max_age)
ext$psi$FSM[ext$mat$a_min]<- ext$m_i$FSM[ext$mat$a_min]
for(i in (ext$mat$a_min+1):ext$max_age)
{
  ext$psi$FSM[i]<- ext$m_i$FSM[i]+
    sum(ext$psi$FSM[ext$mat$a_min:(i-1)]*ext$spawning_period$p_t[i-ext$mat$a_min:(i-1)])
}
rm(i)

saveRDS(ext, "./baseline-parameters/default-parameters_extended.rds")
saveRDS(ext, "./output/parameterizations/default-parameters_extended.rds")
