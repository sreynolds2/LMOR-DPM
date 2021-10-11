
# NUMBER OF EGGS PER SPAWNING FEMALE
## AGE-LENGTH RELATIONSHIP
length_at_age<- function(age=NULL,
                         reps=NULL,
                         inputs=NULL,
                         type="ln_vbgf")
{
  if(type=="vbgf")
  {
    La<- inputs$growth$vbgf$Linf*(1-exp(-inputs$growth$vbgf$k*(age-inputs$growth$vbgf$t0)))
    l<-rnorm(reps, La, inputs$growth$vbgf$sigma)
  }
  if(type=="ln_vbgf")
  {
    La<- inputs$growth$ln_vbgf$Linf*(1-exp(-inputs$growth$ln_vbgf$k*(age-inputs$growth$ln_vbgf$t0)))
    l<- exp(rnorm(reps, log(La), inputs$growth$ln_vbgf$sigma))
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
rm(dat, fit)

eggs<- function(fork_length=NULL,
                mean_fl=NULL,
                sd_fl=NULL,
                a=NULL,
                b=NULL,
                dispersion_param=NULL,
                min_length=NULL)
{
  N<-length(fork_length)
  fl_normalized<- (fork_length - mean_fl)/sd_fl
  egg_no<- rpois(N,exp(rnorm(N,a + b*fl_normalized,dispersion_param)))
  if(any(fork_length<inputs$stage$min_adult_length))
  {
    egg_no[which(fork_length<inputs$stage$min_adult_length)]<- 0
  }
  return(egg_no)
}

## SIMULATE FECUNDITY BY AGE
fecundity_sim<- function(n=1000000,
                         inputs=NULL,
                         run_type="ln_vbgf")
{
  a=1:inputs$max_age
  fec<- lapply(a, function(x)
  {
    lengths<- length_at_age(age=x, reps=n, inputs = inputs, type=run_type)
    mn_lgth<- mean(lengths)
    med_lgth<- median(lengths)
    fecundity<- eggs(fork_length=lengths, 
                     a=inputs$fecundity$intercept, 
                     b=inputs$fecundity$slope, 
                     dispersion_param=inputs$fecundity$disp,
                     min_length = inputs$stage$min_adult_length)
    mn_eggs<- mean(fecundity)
    med_eggs<- median(fecundity)
    if(any(lengths<inputs$stage$min_adult_length) & any(lengths>=inputs$stage$min_adult_length))
    {
      tmp<- which(lengths>=inputs$stage$min_adult_length)
      N_adult<- length(tmp)
      mn_adult_lgth<- mean(lengths[tmp])
      med_adult_lgth<- median(lengths[tmp])
      mn_adult_eggs<- mean(fecundity[tmp])
      med_adult_eggs<- median(fecundity[tmp])
      min_adult_eggs<- min(fecundity[tmp])
    }
    if(all(lengths>=inputs$stage$min_adult_length))
    {
      N_adult<- n
      mn_adult_lgth<- mn_lgth
      med_adult_lgth<- med_lgth
      mn_adult_eggs<- mn_eggs
      med_adult_eggs<- med_eggs
      min_adult_eggs<- min(fecundity)
    }
    if(all(lengths<inputs$stage$min_adult_length))
    {
      N_adult<- 0
      mn_adult_lgth<- NA
      med_adult_lgth<- NA
      mn_adult_eggs<- 0
      med_adult_eggs<- 0
      min_adult_eggs<- 0
    }
    return(list(mean_eggs=mn_eggs, median_eggs=med_eggs,
                mean_adult_eggs=mn_adult_eggs, median_adult_eggs=med_adult_eggs,
                max_eggs=max(fecundity), min_eggs=min(fecundity),
                min_adult_eggs=min_adult_eggs,
                mean_length=mn_lgth, median_length=med_lgth,
                mean_adult_length=mn_adult_lgth, median_adult_length=med_adult_lgth,
                prop_adults=N_adult/n))
  })
  fecundity<- data.frame(Age=a,
                         Mean_Eggs_All=sapply(fec_vbgf, "[[", 1),
                         Median_Eggs_All=sapply(fec_vbgf, "[[", 2),
                         Mean_Eggs_Adults=sapply(fec_vbgf, "[[", 3),
                         Median_Eggs_Adults=sapply(fec_vbgf, "[[", 4),
                         Max_Eggs_Simulated=sapply(fec_vbgf, "[[", 5),
                         Min_Eggs_Simulated=sapply(fec_vbgf, "[[", 6),
                         Min_Eggs_Simulated_Adults=sapply(fec_vbgf, "[[", 7),
                         Mean_Length_All=sapply(fec_vbgf, "[[", 8),
                         Median_Length_All=sapply(fec_vbgf, "[[", 9),
                         Mean_Length_Adult=sapply(fec_vbgf, "[[", 10),
                         Median_Length_Adult=sapply(fec_vbgf, "[[", 11),
                         Proportion_Adults=sapply(fec_vbgf, "[[", 12))
  return(fecundity)
}

#UPDATE INPUTS
if(exists("fecundity_run"))
{
  if(fecundity_run)
  {
    inputs<- readRDS("./baseline-parameters/phi_mi.rds")
    inputs$fecundity$mean_fl<- mean_fl
    inputs$fecundity$sd_fl<- sd_fl
    inputs$fecundity$intercept<- intrcpt
    inputs$fecundity$slope<- slp
    inputs$fecundity$disp<- disp
    saveRDS(inputs, "./baseline-parameters/phi_mi_fec.rds")
  }
}
rm(mean_fl, sd_fl, intrcpt, slp, disp)


# SIMULATIONS
if(exists("new_egg_sim"))
{
  if(new_egg_sim==TRUE)
  {
    ## SIMULATE VBGF FECUNDITY BY AGE
    fecundity_vbgf<- fecundity_sim(n=1000000,
                                   inputs=inputs,
                                   run_type="vbgf")
    write.csv(fecundity_vbgf, "./baseline-parameters/fecundity_estimates_by_age_vbgf.csv",
              row.names = FALSE)
    
    ## SIMULATE ln(VBGF) FECUNDITY BY AGE
    fecundity_ln<- fecundity_sim(n=1000000,
                                 inputs=inputs,
                                 run_type="ln_vbgf")
    write.csv(fecundity_ln, "./baseline-parameters/fecundity_estimates_by_age_ln_vbgf.csv",
              row.names = FALSE)
  }
}

