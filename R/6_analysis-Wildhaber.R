
library("truncnorm")

# WILDHABER ET AL 2017 INPUTS
inps<- list()
## MAXIMUM AGE
inps$max_age<- 41
## AGE-0 SURVIVAL
inps$mu$phi0<- 0.00011
inps$sd$phi0<- sqrt(6.0733*10^(-10))
inps$sd_temp$phi0<- 0
inps$lb$phi0<- 0
inps$ub$phi0<- 0.0004
## LARVAL SURVIVAL
inps$mu$phi_larv<- 0.0510
inps$sd$phi_larv<- sqrt(1.3055*10^(-4))
inps$sd_temp$phi_larv<- 0
inps$lb$phi_larv<- 0
inps$ub$phi_larv<- 1
## AGE-1 SURVIVAL
inps$mu$phi1<- 0.3674
inps$sd$phi1<- sqrt(1.7567*10^(-3))
inps$sd_temp$phi1<- 0
inps$lb$phi1<- 0
inps$ub$phi1<- 1
## AGE-2+ SURVIVAL
inps$mu$phi2plus<- 0.9220
inps$sd$phi2plus<- sqrt(2.3573*10^(-3))
inps$sd_temp$phi2plus<- 0
inps$lb$phi2plus<- 0
inps$ub$phi2plus<- 1
## SEX RATIO
inps$mu$probF<- 0.4708
inps$sd$probF<- sqrt(4.6719*10^(-3))
inps$sd_temp$probF<- sqrt(8.5469*10^(-4))
inps$lb$probF<- 0
inps$ub$probF<- 1
## AGE OF FIRST FEMALE REPRODUCTION
inps$mu$mat_age<- 10
inps$sd$mat_age<- sqrt(1.0833)
inps$sd_temp$mat_age<- 0
inps$lb$mat_age<- 1
inps$ub$mat_age<- inps$max_age
## PROPORTION SPAWNING
inps$mu$p_spawn<- 0.2810
inps$sd$p_spawn<- sqrt(5.7906*10^(-3))
inps$sd_temp$p_spawn<- 0
inps$lb$p_spawn<- 0
inps$ub$p_spawn<- 1
## GROWTH PARAMETERS
inps$vbg$Linf<- 1194#exp(6.9750323)
inps$vbg$k<- 0.07419#exp(-2.4073801)
inps$vbg$t0<- -2.672#-2.1603259
#inps$vbg$sd<- #vector of age specific total sd's
#inps$vbg$sd_temp<- #vector of age specific temporal sd's
## FECUNDITY PARAMETERS
inps$fec$b0<- -45224.64
inps$fec$b1<- 83.69
#inps$fec$sd<- #vector of length specific total sd's
#inps$fec$sd_temp<- #vector of length specific temporal sd's
## INITIAL POPULATION
inps$N0<- c(rep(174,29),rep(173,41-29))

project_pop<- function(inputs=NULL,
                       t_spin=20,
                       t_final=38,
                       reps=2000,
                       run_type="Unpartitioned") #"Partitioned" and "Deterministic" also accepted
{
  # CHECK
  chk<- c(all(names(inputs$mu)==names(inputs$sd)),
          all(names(inputs$mu)==names(inputs$sd_temp)),
          all(names(inputs$mu)==names(inputs$lb)),
          all(names(inputs$mu)==names(inputs$ub)))
  if(all(chk)!=TRUE){return(print("Parameter component mismatch."))}
    
  # REPLICATES
  runs<- lapply(1:reps, function(x)
  {
    # INITALIZE POPULATION
    ## AGE-1+
    Nt<- inputs$N0
    ## AGE-0 (APPROXIMATE FROM NUMBER OF AGE-1's)
    age0<- round(inputs$N0[1]/inputs$mu$phi0)
    # PARAMETER MEANS (UPDATED FOR "PARTITIONED" APPROACH)
    mu<- inputs$mu
    ## REPLICATE SPECIFIC PARAMETERS
    if(run_type=="Partitioned")
    {
      # MAY WANT TO CHECK IF LAPPLY OR AS.LIST/UNLIST IS FASTER
      sig_t<- inputs$sd_temp
      sig_p<- lapply(1:length(sig_t), function(p){inputs$sd[[p]]-sig_t[[p]]})
      mu<- lapply(1:length(sig_p), function(p)
      {
        rtruncnorm(1, max(inputs$lb[[p]], inputs$mu[[p]]-1.96*sig_p[[p]]), 
                   min(inputs$ub[[p]], inputs$mu[[p]]+1.96*sig_p[[p]]), 
                   inputs$mu[[p]], sig_p[[p]])
      })
      names(mu)<- names(inputs$mu)
      tau<- as.list(rnorm(length(sig_t), unlist(sig_t),
                          0.05*unlist(sig_t)))
      names(tau)<- names(inputs$mu)
    }
    if(run_type=="Unpartitioned")
    {
      tau<- as.list(rnorm(length(inputs$sd), unlist(inputs$sd), 
                          0.05*unlist(inputs$sd)))
      names(tau)<- names(inputs$mu)
    }
    if(run_type=="Deterministic")
    {
      tau<- as.list(rep(0,length(inputs$mu)))
      names(tau)<- names(inputs$mu)
    }
    # MATRIX OF ANNUAL POPULATION SIZES
    N<- matrix(0, nrow=inputs$max_age+1, ncol=t_spin+t_final+1)
    N[,1]<- c(age0, Nt)
    # PROJECT POPULATION FORWARD
    for(t in 1:(t_spin+t_final))
    {
      ## YEAR SPECIFIC PARAMETERS
      params<-lapply(1:length(mu), function(p)
      {
        ifelse(tau[[p]]==0,
               mu[[p]],
               rtruncnorm(1, max(inputs$lb[[p]], mu[[p]]-1.96*tau[[p]]), 
                          min(inputs$ub[[p]], mu[[p]]+1.96*tau[[p]]), 
                          mu[[p]], tau[[p]]))
      })
      names(params)<- names(mu)
      # SURVIVAL VECTOR
      phi<- c(params$phi1, rep(params$phi2plus, inputs$max_age-2))
      # OPTIONS:
      # A: PULL MEAN FORK LENGTH FROM 95% CI DISTRIBUTION FOR EACH AGE CLASS i
      # B: PULL Nti FORK LENGTHS FOR EACH AGE CLASS i
      # START WITH A SINCE LESS LIKE IBM
      FLi<- inputs$vbg$Linf*(1-exp(-inputs$vbg$k*(1:inputs$max_age)))
      Ei<- inputs$fec$b0+inputs$fec$b1*FLi
      Ei[1:ceiling(params$mat_age-1)]<- 0
      Ft<-  round(Ei*params$probF*rbinom(inputs$max_age, Nt, params$p_spawn)) 
      Nitplus<- rbinom(inputs$max_age,
                       c(age0, Nt[1:(inputs$max_age-1)]), 
                       c(params$phi0, phi))
      age0<- sum(Ft)
      Nt<- Nitplus
      N[,t+1]<- c(age0, Nt)
    }
    return(N)
  })
  return(list(pop_data=runs, run_type=run_type))
}


growth_rate_checks<- function(pop_data=NULL)
{
  Ntot<- colSums(runs[[1]][2:nrow(runs[[1]]),])
  lambda<- sapply(2:length(Ntot), function(x)
    {
      Ntot[x]/Ntot[x-1]
    })
  diff<- sapply(2:length(lambda), function(x)
    {
      lambda[x]-lambda[x-1]
    })
  which(abs(diff)<=0.0001)
  which(abs(diff)<=0.001)
  meanL<- sapply(1:length(lambda), function(x)
    {
      mean(lambda[x:length(lambda)])
    })
  diff2<- sapply(2:length(meanL), function(x)
  {
    meanL[x]-meanL[x-1]
  })
  which(abs(diff2)<=0.0001)
  which(abs(diff2)<=0.001)
  geo_meanL<- sapply(1:(length(Ntot)-1), function(x)
  {
    Ntot[length(Ntot)]/Ntot[x]
  })
  diff3<- sapply(2:length(geo_meanL), function(x)
  {
    geo_meanL[x]-geo_meanL[x-1]
  })
  which(abs(diff3)<=0.0001)
  which(abs(diff3)<=0.001)
  return(list(annual_lambda_diff=diff,
              mean_lambda=list(mean_lambda=meanL, difference=diff2),
              geom_mean_lambda=list(mean_lambda=geo_meanL, difference=diff3)))
}

growth_rate<- function(data=NULL,
                       t_spin=20,
                       t_final=38)
{
  pop_data<- data$pop_data
  # CHECK
  if(t_final > ncol(pop_data[[1]])-1)
  {
    return(print("Population was not projected long enough for given time period."))
  }
  if(t_spin + t_final > ncol(pop_data[[1]])-1)
  {
    return(print("Population was not projected long enough for the given combination of spin-up period and years projected."))
  }
  out<- lapply(1:length(pop_data), function(x)
  {
    Ntot<- colSums(pop_data[[x]][2:nrow(pop_data[[x]]),])
    lambda<- sapply(2:(t_final+1)+t_spin, function(x)
    {
      Ntot[x]/Ntot[x-1]
    })
    return(data.frame(run_type=data$run_type,
                      rep=x,
                      final_size=Ntot[t_final+1], 
                      mean_lambda=mean(lambda),
                      median_lambda=median(lambda),
                      geom_mn_lambda=prod(lambda)^(1/length(lambda)),
                      spin_up_yrs=t_spin,
                      total_yrs=t_final))
  })
  out<- do.call("rbind", out)
  return(out)
}


quasiextinction<- function(data=NULL,
                           threshold=50,
                           years=200,
                           reps=2000,
                           output_yr=200,
                           output_reps=2000)
{
  pop_data<- data$pop_data
  # ERROR CHECK
  if(years>ncol(pop_data[[1]])-1)
  {
    return(print("Population data does not contain enough years for the desired estimate."))
  }
  if(reps>length(pop_data))
  {
    return(print("Population data does not contain enough replicates for the desired estimate."))
  }
  # COMPUTE FRACTION BELOW THRESHOLD
  Ntot<- sapply(1:reps, function(x)
  {
    colSums(pop_data[[x]][2:nrow(pop_data[[x]]),])
  })
  yr<- sapply(1:reps, function(z)
  {
    val<- ifelse(all(Ntot[,z]>=threshold),
                 years*10,
                 min(which(Ntot[,z]<threshold))-1)
    return(val)
  })
  frac<- sapply(1:years, function(y)
  {
    length(which(yr<=y))/reps
  })
  # FRACTION EXTINCT
  if(output_reps==reps)
  {
    output_frac<- frac
  }
  if(output_reps<reps)
  {
    output_frac<- frac*reps/output_reps
  }
  # TIME TO 50%
  time_fifty<- ifelse(output_frac[output_yr]>=0.5,
                      min(which(output_frac>=0.5)),
                      NA)
  # TIME TO 80%
  time_eighty<- ifelse(output_frac[output_yr]>=0.8,
                       min(which(output_frac>=0.8)),
                       NA)
  # EXPECTED TIME
  if(all(yr[1:output_reps]<=output_yr))
  {
    tmp<- data.frame(extinction_yr=yr[1:output_reps], freq=1)
    tmp<- aggregate(freq~extinction_yr, tmp, sum)
    tmp$prob<- tmp$freq/sum(tmp$freq)
    exp_time<- sum(tmp$extinction_yr*tmp$prob)
  }
  if(!all(yr[1:output_reps]<=output_yr))
  {
    exp_time<- NA
  }
  out<- list(ext_dat=list(extinction_yr=yr, frac_extinct=frac),
             ext_table=data.frame(run_type=data$run_type,
                                  frac_ext=output_frac[output_yr],
                                  T_50=time_fifty,
                                  T_80=time_eighty,
                                  E_T=exp_time,
                                  threshold=threshold,
                                  years=output_yr,
                                  reps=output_reps))
  return(out)
}


### TESTS
tm<- Sys.time()
datD<- project_pop(inputs = inps,
                   t_spin=0,
                   t_final=150,
                   reps=2000,
                   run_type = "Deterministic")
Sys.time()-tm
lambdaD<- growth_rate(data = datD,
                      t_spin = 20,
                      t_final = 38)
extD<- quasiextinction(data = datD, reps = 500, output_reps = 500)


tm<- Sys.time()
datP<- project_pop(inputs = inps,
                  t_final=200,
                  reps=500,
                  run_type = "Partitioned")
Sys.time()-tm
lambdaP<- growth_rate(data = datP)
extP<- quasiextinction(data = datP, reps = 500, output_reps = 500)

tm<- Sys.time()
datU<- project_pop(inputs = inps,
                   t_final=200,
                   reps=500,
                   run_type = "Unpartitioned")
Sys.time()-tm
lambdaU<- growth_rate(data = datU)
extU<- quasiextinction(data = datU, reps = 500, output_reps = 500)


lambda<- rbind(lambdaD, lambdaP)
lambda<- rbind(lambda, lambdaU)
lambda<- lambda[order(lambda$rep),]

ext<- rbind(extD$ext_table, extU$ext_table)
ext<- rbind(ext, extP$ext_table)


boundaries<- function(pop_data=NULL)
{
  
}

reproductive_value<- function(pop_data=NULL)
{
  
}

stable_age0<- function(pop_data=NULL)
{
  
}
  #  ## DETERMINISTIC OUTPUTS
  #   Li<- inputs$Linf*(1-exp(-1*inputs$k*(1:inputs$max_age-inputs$t0)))
  #   Ei<- inputs$fec_b0+inputs$fec_b1*Li
  #   Ei[1:(inputs$mat_age-1)]<- 0
  #   Fi<- Ei*inputs$probF*inputs$p_spawn
  #   A<- matrix(0, inputs$max_age, inputs$max_age)
  #   A[1,]<- Fi
  #   for(i in 2:inputs$max_age)
  #   {
  #     
  #     A[i,i-1]<- phi[i-1]
  #   }
  #   N<- matrix(0,nrow=inputs$max_age, ncol=t_final)
  #   N[,1]<- A%*%inputs$N0
  #   for(i in 2:t_final)
  #   {
  #     N[,i]<- A%*%N[,i-1]
  #   }
  #   library(expm) 
  #   N2<- sapply(1:t_final, function(t){A%^%t%*%inputs$N0})
  # }
  # return(deterministic=N, deterministic2=N2, inputs=inputs)
  


# L<- inps$Linf*(1-exp(-1*inps$k*(0:41-inps$t0)))
# f<- inps$fec_b0+inps$fec_b1*775:1058
# L
# f
# (205849.3-inps$fec_b0)/1058
# 205849.3-inps$fec_b1*1058
# 
# (0-inps$fec_b0)/775
# 0-inps$fec_b1*775
# 
# (0-inps$fec_b0)/inps$fec_b1
# (205849.3-inps$fec_b0)/inps$fec_b1
# 
# m<- (205849.3-0)/(1058-775)
# y0<- m*(0-1058)+205849.3
# m*1058+y0
# m*775+y0
# m
# y0

sig_Linf<- 86.9
sig_k<- 0.07419
sig_t0<- 0.2175
vf2<- function(a)
{
  (1-exp(-1*inps$k*(a-inps$t0)))^2*sig_Linf^2 +
    (inps$Linf*exp(-1*inps$k*(a-inps$t0))*(a-inps$t0))^2*sig_k^2 +
    (inps$Linf*exp(-1*inps$k*(a-inps$t0))*inps$k)^2*sig_t0^2 
}
vf2(0:27)
# vf<- function(a)
# {
#   (1-exp(-1*inps$k*(a-inps$t0)))^2*sig_Linf^2 +
#     (inps$Linf*exp(-1*inps$k*(a-inps$t0))*(a-inps$t0))^2*sig_k^2 +
#     (inps$Linf*exp(-1*inps$k*(a-inps$t0))*inps$k)^2*sig_t0^2 + 
#     2*(1-exp(-1*inps$k*(a-inps$t0)))*(inps$Linf*exp(-1*inps$k*(a-inps$t0))*(a-inps$t0))*sig_Linf*sig_k +
#     2*(1-exp(-1*inps$k*(a-inps$t0)))*(inps$Linf*exp(-1*inps$k*(a-inps$t0))*inps$k)*sig_Linf*sig_t0 + 
#     2*(inps$Linf*exp(-1*inps$k*(a-inps$t0))*(a-inps$t0))*(inps$Linf*exp(-1*inps$k*(a-inps$t0))*inps$k)*sig_k*sig_t0
# }