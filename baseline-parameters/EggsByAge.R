
inputs<- readRDS("./baseline-parameters/phi_mi.rds")

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

eggs<- function(fork_length=NULL,
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
n=1000
a=1:inputs$max_age
fec_vbgf<- lapply(a, function(x)
{
  lengths<- length_at_age(age=x, reps=n, inputs = inputs, type="vbgf")
  mn_lgth<- mean(lengths)
  med_lgth<- median(lengths)
  fecundity<- eggs(fork_length=lengths, a=intrcpt, b=slp, dispersion_param=disp,
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
fecundity_vbgf<- data.frame(Age=a,
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
write.csv(fecundity_vbgf, "./baseline-parameters/fecundity_estimates_by_age_vbgf.csv",
          row.names = FALSE)

fec_ln<- lapply(a, function(x)
{
  lengths<- length_at_age(age=x, reps=n, inputs=inputs, type="ln_vbgf")
  mn_lgth<- mean(lengths)
  med_lgth<- median(lengths)
  fecundity<- eggs(fork_length=lengths, a=intrcpt, b=slp, dispersion_param=disp,
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
fecundity_ln<- data.frame(Age=a,
                          Mean_Eggs_All=sapply(fec_ln, "[[", 1),
                          Median_Eggs_All=sapply(fec_ln, "[[", 2),
                          Mean_Eggs_Adults=sapply(fec_ln, "[[", 3),
                          Median_Eggs_Adults=sapply(fec_ln, "[[", 4),
                          Max_Eggs_Simulated=sapply(fec_ln, "[[", 5),
                          Min_Eggs_Simulated=sapply(fec_ln, "[[", 6),
                          Min_Eggs_Simulated_Adults=sapply(fec_ln, "[[", 7),
                          Mean_Length_All=sapply(fec_ln, "[[", 8),
                          Median_Length_All=sapply(fec_ln, "[[", 9),
                          Mean_Length_Adult=sapply(fec_ln, "[[", 10),
                          Median_Length_Adult=sapply(fec_ln, "[[", 11),
                          Proportion_Adults=sapply(fec_ln, "[[", 12))
write.csv(fecundity_ln, "./baseline-parameters/fecundity_estimates_by_age_ln_vbgf.csv",
          row.names = FALSE)

## FIGURE COMPARISON OF APPROACHES -- NEEDS TO BE RECODED
age<-1:inputs$max_age
### VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ADULTS ONLY
plot(age, fecundity_vbgf$Mean_Eggs_Adults, xlab="Age", ylab="Mean Number of Eggs", col="blue", pch=19)
points(age, fecundity_vbgf2$Mean_Eggs_Adults, pch=19)
points(age, fecundity_vbgf5$Mean_Eggs_Adults, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ALL BY AGE: ADULTS WEIGHTED
plot(age, fecundity_vbgf$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
points(age, fecundity_vbgf2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19)
points(age, fecundity_vbgf5$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### LN_VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ADULTS ONLY
plot(age, fecundity_ln$Mean_Eggs_Adults, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
points(age, fecundity_ln2$Mean_Eggs_Adults, pch=19)
points(age, fecundity_ln5$Mean_Eggs_Adults, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### LN_VBGF STOCHASTIC DRAW COMPARISON; EGGS FROM ALL BY AGE: ADULTS WEIGHTED
plot(age, fecundity_ln$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, xlab="Age", ylab="Number of Eggs", col="blue", pch=19)
points(age, fecundity_ln2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19)
points(age, fecundity_ln5$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="red")
legend("bottomright", c("1 million", "2 million", "5 million"),
       col=c("black", "blue", "red"), pch=19, bty="n")

### VBGF EGGS FROM ALL BY AGE: AVERAGE OF DRAWS VS. ADULTS WEIGHTED
plot(age, fecundity_vbgf2$Mean_Eggs_All, xlab="Age", ylab="Number of Eggs", pch=19)
points(age, fecundity_vbgf2$Mean_Eggs_Adults*p_stage$vbgf$adult, pch=19, col="blue")
legend("bottomright", c("Eggs from All Fish", "Eggs from Adults Weighted by Proportion"),
       col=c("black", "blue"), pch=19, bty="n")

### LN_VBGF EGGS FROM ALL BY AGE: AVERAGE OF DRAWS VS. ADULTS WEIGHTED
plot(age, fecundity_ln2$Mean_Eggs_All, xlab="Age", ylab="Number of Eggs", pch=19)
points(age, fecundity_ln2$Mean_Eggs_Adults*p_stage$ln_vbgf$adult, pch=19, col="blue")
legend("bottomright", c("Eggs from All Fish", "Eggs from Adults Weighted by Proportion"),
       col=c("black", "blue"), pch=19, bty="n")

rm(dat, mean_fl, sd_fl, fit, intrcpt, slp, disp, n, a, fec_vbgf, fecundity_vbgf,
   fec_ln, fecundity_ln, length_at_age, eggs, age)
