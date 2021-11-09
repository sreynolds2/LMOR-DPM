
create_params<- function(max_age=41, # MAXIMUM FEMALE AGE OF REPRODUCTION
                         a_min=8,  # MINIMUM MATURATION AGE
                         # MINIMUM STAGE BASED LENGTHS
                         min_adult_length=800,
                         min_sa_length=600,
                         # SURVIVAL PARAMETERS BY STAGE
                         phi_1=0.151,
                         phi_juv=0.780,
                         phi_sa=0.932,
                         phi_adult=0.976,
                         # GROWTH MODEL PARAMETERS
                         growth_fit=list(Linf=1065.684, k=0.09023894, t0=-2.166103, sigma=0.1286695), 
                            #list(Linf=1051.872, k=0.1008206, t0=-1.830915, sigma=47.56881),
                         growth_fit_type="ln_vbgf", #"vbgf"
                         # FECUNDITY MODEL PARAMETERS
                         fec_run=FALSE,
                         fec_file=c("fecundity_estimates_by_age_", "_2mil.csv"),
                         mean_fl=1260.167,
                         sd_fl=277.404,
                         b0=10.76789,
                         b1=0.6181171,
                         dispersion_param=0.2986811,
                         reps=2000000,
                         # AGE-0 SURVIVAL
                         phi0=0.000075,
                         # SEX RATIO
                         probF=0.5,
                         # RR SPAWNING PROPORTION
                         gamma=1,
                         # REPRODUCTIVE PERIOD
                         max_period=5,
                         probs=c(0, 8/21, 13/21*3/5),
                         input_id=NULL)
{
  inps<- list()
  # MAXIMUM AGE
  inps$max_age<- max_age
  # MINIMUM MATURATION AGE
  inps$a_min<- a_min
  # AGE-1+ SURVIVALS AND MATURATION FROM STAGE AND GROWTH
  source("./baseline-parameters/Survival and Maturation.r")
  if(growth_fit_type=="vbgf")
  {
    vbgf<- growth_fit
    ln_vbgf=NULL
  }
  if(growth_fit_type=="ln_vbgf")
  {
    vbgf=NULL
    ln_vbgf=growth_fit
  }
  params<- surv_and_mat_from_stage(max_age = max_age+1,
                                   a_min = a_min,
                                   min_adult_length = min_adult_length,
                                   min_sa_length = min_sa_length,
                                   phi_1 = phi_1,
                                   phi_juv = phi_juv,
                                   phi_sa = phi_sa,
                                   phi_adult = phi_adult,
                                   vbgf = vbgf,
                                   ln_vbgf = ln_vbgf)
  if(growth_fit_type=="vbgf")
  {
    params$phi$ln_vbgf<- NULL
    params$m_i$ln_vbgf<- NULL
    params$growth$ln_vbgf<- NULL
    inps$phi$vbgf<- params$phi$vbgf
    inps$m_i$vbgf<- params$m_i$vbgf
    inps$growth$vbgf<- params$growth$vbgf
  }
  if(growth_fit_type=="ln_vbgf")
  {
    params$phi$vbgf<- NULL
    params$m_i$vbgf<- NULL
    params$growth$vbgf<- NULL
    inps$phi$ln_vbgf<- params$phi$ln_vbgf
    inps$m_i$ln_vbgf<- params$m_i$ln_vbgf
    inps$growth$ln_vbgf<- params$growth$ln_vbgf
  }
  saveRDS(params, 
          paste0("./output/parameterizations/phi_mi_",
                 input_id, ".rds"))
  inps$stage<- params$stage
  # FERTILITIES
  ## SEX RATIO
  inps$probF<- probF
  ## REPRODUCTIVE PERIOD
  inps$spawning_period$max_period<- max_period
  if(length(probs)<max_period)
  {
    R<-probs
    R[1]<- 1-probs[1]
    if(length(probs>1))
    {
      for(i in 2:length(probs))
      {
        R[i]<- 1-R[i]/R[i-1]
      }
    }
    probs<-c(probs, 
             prod(R)*1/sum(2^(0:(max_period-4)))*2^((max_period-4):0))
  }
  if(length(probs)>max_period | 
     (length(probs)==max_period & sum(probs)!=1))
  {return(print("Repoductive period distribution error."))}
  inps$spawning_period$p_t<- c(probs, rep(0, params$max_age-max_period+1))
  ## PROPORTION OF FEMALES THAT ARE REPRODUCTIVELY READY TO SPAWN
  psi<- rep(0, max_age)
  psi[a_min]<- inps$m_i[[growth_fit_type]][a_min]
  for(i in (a_min+1):(max_age+1))
  {
    psi[i]<- inps$m_i[[growth_fit_type]][i]+
      sum(psi[a_min:(i-1)]*inps$spawning_period$p_t[i-a_min:(i-1)])
  }
  inps$psi<- list(psi)
  names(inps$psi)<- growth_fit_type
  rm(psi, i)
  ## FECUNDITY
  if(fec_run==TRUE)
  {
    ### PULL SIMULATION FUNCTIONS
    source("./baseline-parameters/EggsByAge.r")
    ### ADD IN FECUNDITY PARAMETERS
    params$fecundity$mean_fl<- mean_fl
    params$fecundity$sd_fl<- sd_fl
    params$fecundity$intercept<- b0
    params$fecundity$slope<- b1
    params$fecundity$disp<- dispersion_param
    saveRDS(params, paste0("./output/parameterizations/phi_mi_fec_",
                           input_id, ".rds"))
    ### SIMULATE NUMBER OF EGGS AND SAVE TO INPUTS
    fec<- fecundity_sim(n=reps,
                        inputs=params,
                        run_type=growth_fit_type)
  }
  if(fec_run==FALSE)
  {
    fec<- read.csv(paste0("./baseline-parameters/", fec_file[1], growth_fit_type, 
                        fec_file[2]))
  }
  inps$fecundity$mean_fl<- mean_fl
  inps$fecundity$sd_fl<- sd_fl
  inps$fecundity$intercept<- b0
  inps$fecundity$slope<- b1
  inps$fecundity$disp<- dispersion_param
  inps$eggs<- list(fec$Mean_Eggs_Adults)
  names(inps$eggs)<- growth_fit_type
  ## AGE-0 SURVIVAL
  inps$phi0<- phi0
  ## PROPORTION OF REPRODUCTIVELY READY FEMALES THAT SPAWN
  inps$gamma<- gamma
  # SAVE PARAMETERS
  saveRDS(inps, paste0("./output/parameterizations/parameters_",
                       input_id, ".rds"))
  return(inps)
}

project_pop<- function(inputs=NULL,
                       growth_type="vbgf",
                       t_spin=20,
                       t_final=38,
                       reps=2000,
                       run_type="Unpartitioned", #"Partitioned" and "Deterministic" also accepted
                       save_path=FALSE)
{
  # SUBSET INPUTS
  inputs$mu<- c(inputs$mu, inputs[[growth_type]]$mu)
  inputs$sd<- c(inputs$sd, inputs[[growth_type]]$sd)
  inputs$sd_temp<- c(inputs$sd_temp, inputs[[growth_type]]$sd_temp)
  inputs$lb<- c(inputs$lb, inputs[[growth_type]]$lb)
  inputs$ub<- c(inputs$ub, inputs[[growth_type]]$ub)
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
        rtruncnorm(length(sig_p[[p]]), 
                   sapply(1:length(sig_p[[p]]), function(x)
                    {
                     ifelse(sig_p[[p]][x]==0,
                            inputs$lb[[p]][x],
                            max(inputs$lb[[p]][x], 
                                inputs$mu[[p]][x]-1.96*sig_p[[p]][x]))
                    }), 
                   sapply(1:length(sig_p[[p]]), function(x)
                    {
                     ifelse(sig_p[[p]][x]==0,
                            inputs$ub[[p]][x],
                            min(inputs$ub[[p]][x], 
                                inputs$mu[[p]][x]+1.96*sig_p[[p]][x]))
                    }), 
                   inputs$mu[[p]], sig_p[[p]])
      })
      names(mu)<- names(inputs$mu)
      if(any(mu$p_spawn=="NaN"))
      {
        mu$p_spawn[which(mu$p_spawn=="NaN")]<- 0
      }
      tau<- lapply(1:length(inputs$sd), function(p)
      {
        rnorm(length(sig_t[[p]]), sig_t[[p]],
                          0.05*sig_t[[p]])
      })
      names(tau)<- names(inputs$mu)
    }
    if(run_type=="Unpartitioned")
    {
      tau<- lapply(1:length(inputs$sd), function(p)
      {
        rnorm(length(inputs$sd[[p]]), inputs$sd[[p]], 
              0.05*inputs$sd[[p]])
      })
      names(tau)<- names(inputs$mu)
    }
    if(run_type=="Deterministic")
    {
      tau<- lapply(1:length(inputs$mu), function(p)
      {
        rep(0,length(inputs$mu[[p]]))
      })
      names(tau)<- names(inputs$mu)
    }
    # MATRIX OF ANNUAL POPULATION SIZES
    N<- matrix(0, nrow=inputs$max_age+1, ncol=t_spin+t_final+1)
    N[,1]<- c(age0, Nt)
    # INITIALIZE PARAMETER PATH
    omega<- NULL
    # PROJECT POPULATION FORWARD
    for(t in 1:(t_spin+t_final))
    {
      ## YEAR SPECIFIC PARAMETERS
      params<-lapply(1:length(mu), function(p)
      {
        out<- sapply(1:length(mu[[p]]), function(x)
        {
          ifelse(tau[[p]][x]==0,
                  mu[[p]][x],
                  rtruncnorm(1, max(inputs$lb[[p]][x], mu[[p]][x]-1.96*tau[[p]][x]), 
                            min(inputs$ub[[p]][x], mu[[p]][x]+1.96*tau[[p]][x]), 
                            mu[[p]][x], tau[[p]][x]))
        })
        return(out)
      })
      names(params)<- names(mu)
      # FERTILITIES
      # OPTIONS:
      # A: PULL MEAN FORK LENGTH FROM 95% CI DISTRIBUTION FOR EACH AGE CLASS i
      # B: PULL Nti FORK LENGTHS FOR EACH AGE CLASS i
      # START WITH A SINCE LESS LIKE IBM
      #
      ## MODIFY THE FOLLOWING TO REFLECT 95% CIs???
      ## YES AND LET's DO LB AND UBs TOO
      #
      # FLi<- rnorm(inputs$max_age, 
      #             inputs[[growth_type]]$Linf*
      #               (1-exp(-inputs[[growth_type]]$k*
      #                        (1:inputs$max_age-inputs[[growth_type]]$t0))),
      #             inputs[[growth_type]]$fit$sd$l_at_a)
      FLi<- inputs[[growth_type]]$Linf*(1-exp(-inputs[[growth_type]]$k*(1:inputs$max_age-inputs[[growth_type]]$t0)))
      FLi_norm<- (FLi-inputs$fec$mean_fl)/inputs$fec$sd_fl
      # Ei<- rpois(length(FLi), 
      #            exp(rnorm(length(FLi),
      #                      inputs$fec$b0+inputs$fec$b1*FLi_norm,
      #                      inputs$fec$fit$sd$ln_eggs_at_l)))
      Ei<- exp(inputs$fec$b0+inputs$fec$b1*FLi_norm)
      Ei[Ei<0]<- 0
      if(save_path==TRUE)
      {
        tmp<- as.list(Ei)
        names(tmp)<- paste0("E", 1:inputs$max_age)
        tmp<- as.data.frame(c(params, tmp))
        omega<- rbind(omega, tmp)
      }
      Ei[1:ceiling(params$mat_age-1)]<- 0
      Ft<-  round(Ei*params$probF*rbinom(inputs$max_age, Nt, params$p_spawn))
      # UPDATE ABUNDANCES
      Nitplus<- rbinom(inputs$max_age,
                       c(age0, Nt[1:(inputs$max_age-1)]), 
                       c(params$phi0, params$phi))
      age0<- sum(Ft)
      Nt<- Nitplus
      N[,t+1]<- c(age0, Nt)
    }
    return(list(N=N, omega=omega))
  })
  omega<- lapply(runs, "[[", 2)
  runs<- lapply(runs, "[[", 1)
  return(list(pop_data=runs, run_type=run_type, 
              input_values=inputs, param_path=omega))
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

## NEED TO REVISIT
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


## CONSTRUCT THE LESLIE MATRIX, PERFORM EIGEN-, SENSITIVITY, & 
## ELASTICITY ANALYSES
matrix_eigen_analysis<- function(inputs,
                                 growth_type="vbgf")
{
  # ERROR HANDLING
  if(!all(sapply(c("max_age", "phi", "psi", "eggs", "probF", "gamma", "phi0"), 
                 exists, where=inputs)))
  {
    return(print("The inputs file must contain all of the following:
                  'max_age', 'phi', 'psi', 'eggs', 'probF', 'gamma', 'phi0_MR'."))
  }
  if(!(exists(growth_type, where=inputs$phi) & 
       exists(growth_type, where=inputs$psi) & 
       exists(growth_type, where=inputs$eggs)))
  {
    return(print("The inputs file must contain the specified growth function
                  type for inputs 'phi', 'psi', and 'eggs'."))
  }
  max_age<- inputs$max_age  
  phi<- inputs$phi[[growth_type]]
  psi<- inputs$psi[[growth_type]]
  eggs<- inputs$eggs[[growth_type]]
  gamma<- inputs$gamma
  sexratio<- inputs$probF
  phi0<- inputs$phi0
  if(!(all(sapply(c(max_age, sexratio, phi0, gamma), 
                  length)==1) & 
       is.numeric(c(max_age, sexratio, phi0, gamma))))
  {
    return(print("Inputs 'max_age', 'probF', 'gamma', and 'phi0' should all be numerical values of length 1."))
  }
  if(!all(c(length(psi), length(eggs))>=max_age+1))
  {
    return(print("Inputs 'psi' and 'eggs' need to be a vector at least 'max_age'+1 in length (one entry for each age plus age max_age+1)."))
  }
  if(length(phi)<max_age)
  {
    return(print("Input 'phi' needs to be a vector of survivals with length max_age."))
  }
  phi_plus<- phi[max_age]
  psi_plus<- psi[max_age+1]
  eggs_plus<- eggs[max_age+1]
  phi<- phi[1:(max_age-1)]
  psi<- psi[1:max_age]
  eggs<- eggs[1:max_age]
  # BUILD LESLIE MATRIX
  A<- matrix(0,max_age,max_age)
  ## SURVIVAL VALUES
  A[cbind(2:max_age,1:(max_age-1))]<- phi
  ## FERTILITY VALUES
  A[1,] <- psi*gamma*eggs*sexratio*phi0
  Aminus<- A[1:(max_age-1), 1:(max_age-1)]
  Aplus<- rbind(A, c(rep(0, max_age-1), phi_plus))
  Aplus<- cbind(Aplus, 
                c(psi_plus*gamma*eggs_plus*sexratio*phi0, rep(0, max_age)))
  
  # EIGENANALYSIS 
  ea<- eigen.analysis(A)
  ea$A<- A
  
  # PARAMETER SENSITIVITIES
  sens<- ea$sensitivities 
  ea$sensitivities<- list()
  ea$sensitivities$entries<- sens
  ea$sensitivities$phi<- rep(0, max_age-1)
  for(i in 1:(max_age-1))
  {
    ea$sensitivities$phi[i]<- sens[i+1,i]
  }
  ea$sensitivities$psi<- sens[1,]*eggs*sexratio*phi0
  # ANALYZE A_MAT_H, K_MAT, AND RR PERIOD PROBS HERE...USE "INPUTS" 
  # AS FUNCTION INPUT
  ea$sensitivities$gamma<- sum(sens[1,]*psi*eggs*sexratio*phi0)
  ea$sensitivities$eggs<- sens[1,]*psi*gamma*sexratio*phi0
  ea$sensitivities$sexratio<- sum(sens[1,]*psi*gamma*eggs*phi0)
  ea$sensitivities$phi0<- sum(sens[1,]*psi*gamma*eggs*sexratio)
  ea_minus<- eigen.analysis(Aminus)
  ea_plus<- eigen.analysis(Aplus)
  ea$sensitivities$max_age<- (ea_plus$lambda1-ea_minus$lambda1)/2
  
  # PARAMETER ELASATICITIES
  elas<- ea$elasticities
  ea$elasticities<- list()
  ea$elasticities$entries<- elas
  rm(elas)
  ea$elasticities$phi<- ea$sensitivities$phi*phi/ea$lambda1
  ea$elasticities$psi<- ea$sensitivities$psi*psi/ea$lambda1
  ea$elasticities$gamma<- ea$sensitivities$gamma*gamma/ea$lambda1
  ea$elasticities$eggs<- ea$sensitivities$eggs*eggs/ea$lambda1
  ea$elasticities$sexratio<- ea$sensitivities$sexratio*sexratio/ea$lambda1
  ea$elasticities$phi0<- ea$sensitivities$phi0*phi0/ea$lambda1
  ea$elasticities$max_age<- ea$sensitivities$max_age*max_age/ea$lambda1
  
  # ADD INPUTS
  ea$inputs<- inputs
  
  # RETURN
  return(ea)
}


## PRODUCE TOP SENSITIVITY AND ELASTICITY TABLES
sens_elas_table<- function(data=NULL,
                           number=10,
                           sensitvities=TRUE,
                           elasticities=TRUE)
{
  sens<-data$sensitivities
  sens$F_i<- sens$entries[1,]
  sens$entries<-NULL
  elas<-data$elasticities
  elas$F_i<- elas$entries[1,] 
  elas$entries<-NULL

  # SENSITIVITIES
  sens_tbl<- data.frame(Parameter=names(unlist(sens)),
                        Sensitivity=unlist(sens))
  sens_tbl$Parameter<-gsub("F_i", "F_", sens_tbl$Parameter)
  sens_tbl$Type<- ifelse(sens_tbl$Parameter %in% paste0("F_", 1:60),
                         "Entry", 
                         ifelse(sens_tbl$Parameter %in% paste0("phi", 1:59),
                                "Both", "Parameter"))
  sens_tbl<- sens_tbl[order(abs(sens_tbl$Sensitivity), decreasing = TRUE),]
  sens_tbl$rank<-1:nrow(sens_tbl)
  rownames(sens_tbl)<- sens_tbl$rank
  sens_tbl$Sensitivity<- round(sens_tbl$Sensitivity,6)
  indx<-max(which(sens_tbl$Sensitivity[number]==sens_tbl$Sensitivity))
  top_sens<- sens_tbl[1:indx,]
  
  # ELASTICITIES
  elas_tbl<- data.frame(Parameter=names(unlist(elas)),
                        Elasticity=unlist(elas))
  elas_tbl$Parameter<-gsub("F_i", "F_", elas_tbl$Parameter)
  elas_tbl$Type<- ifelse(elas_tbl$Parameter %in% paste0("F_", 1:60),
                         "Entry", 
                         ifelse(elas_tbl$Parameter %in% paste0("phi", 1:59),
                                "Both", "Parameter"))
  elas_tbl<- elas_tbl[order(abs(elas_tbl$Elasticity), decreasing = TRUE),]
  elas_tbl$rank<-1:nrow(elas_tbl)
  rownames(elas_tbl)<- elas_tbl$rank
  elas_tbl$Elasticity<- round(elas_tbl$Elasticity,6)
  indx<-max(which(elas_tbl$Elasticity[number]==elas_tbl$Elasticity))
  top_elas<- elas_tbl[1:indx,]
  out<-list(Sensitivities=top_sens, Elasticities=top_elas)
  return(out)
}


## BOUNDARY AGE-0 SURVIVAL & SEXRATIO PRODUCT VALUE
boundary_vals<- function(inputs=NULL,
                         growth_type="vbgf")
{
  # ERROR HANDLING
  if(!all(sapply(c("max_age", "phi", "psi", "gamma", "eggs", "probF"), exists, where=inputs)))
  {
    return(print("The inputs file must contain all of the following:
                 'max_age', 'phi', 'psi', 'gamma', 'eggs', 'probF'."))
    
  }
  if(!(exists(growth_type, where=inputs$phi) & 
       exists(growth_type, where=inputs$psi) & 
       exists(growth_type, where=inputs$eggs)))
  {
    return(print("The inputs file must contain the specified growth function
                  type for inputs 'phi', 'psi', and 'eggs'."))
  }
  max_age<- inputs$max_age  
  phi<- inputs$phi[[growth_type]]
  psi<- inputs$psi[[growth_type]]
  gamma<- inputs$gamma
  eggs<- inputs$eggs[[growth_type]]
  sexratio<- inputs$probF
  if(!(all(sapply(c(max_age, gamma, sexratio),length)==1) &
       is.numeric(c(max_age, gamma, sexratio))))
  {
    return(print("Inputs 'max_age', 'gamma', and 'probF' should be numerical values of length 1."))
  }
  if(!all(c(length(psi), length(eggs))>=max_age))
  {
    return(print("Inputs 'psi' and 'eggs' need to be a vector at least 'max_age' in length (one value for each age)."))
  }
  if(length(phi)<max_age-1)
  {
    return(print("Input 'phi' needs to be a vector of survivals with length at least max_age-1."))
  }
  phi<- phi[1:(max_age-1)]
  psi<- psi[1:max_age]
  eggs<- eggs[1:max_age]
  inputs$phi<- phi
  inputs$psi<- psi
  inputs$eggs<- eggs
  # REARRANGE EULER-LOTKA EQUATION FOR CHARACTERISTIC EQUATION WITH LAMBDA=1
  denom<- sapply(2:length(psi), function(i)
  {
    out<- psi[i]*gamma*eggs[i]*sexratio*prod(phi[1:(i-1)])
    return(out)
  })
  denom<- c(psi[1]*gamma*eggs[1]*sexratio, denom) 
  phi0<- 1/sum(denom)
  ## PRODUCT VALUE FOR AGE0, AGE1, SEXRATIO CURVES
  denom2<- sapply(3:length(psi), function(i)
  {
    out<- psi[i]*eggs[i]*prod(phi[2:(i-1)])
    return(out)
  })
  denom2<- c(psi[2]*eggs[2], denom2)
  prod<- gamma*sexratio*phi0*phi[1]
  B1<- sum(denom2)
  check2<- round(abs(prod-1/B1),14)
  if(psi[1]*gamma*eggs[1]!=0){check2<- "FAILED"}
  
  # CHECK LAMBDA=1 IS THE LARGEST E-VALUE
  ## BUILD LESLIE MATRIX
  A<- matrix(0,max_age,max_age)
  ### SURVIVAL VALUES
  A[cbind(2:max_age,1:(max_age-1))]<- phi
  ### FERTILITY VALUES
  A[1,] <- psi*gamma*eggs*sexratio*phi0
  ## EIGENANALYSIS 
  ea<- eigen.analysis(A)
  check<- round(abs(ea$lambda1-1),14)
  
  # REPORT PHI0 AND CHECK
  inputs$boundary$phi0<- phi0
  inputs$boundary$check<- check
  inputs$boundary$product<- prod
  inputs$boundary$B1<- B1
  inputs$boundary$check2<- check2
  inputs$boundary$growth_type<- growth_type
  return(inputs)
}


## SURVIVAL BOUNDARY CURVES
phi0_phi1_sexratio_curves<- function(boundary_inputs=NULL,
                                     probF=c(0.3,0.5,0.7))
{
  # ERROR HANDLING
  if(!all(sapply(c("product", "check2"), exists, where=boundary_inputs$boundary)))
  {
    return(print("Boundary inputs must contain the list boundary with entries 'product' and 'check2'."))
  }
  if(boundary_inputs$boundary$check2!=0)
  {
    return(print("Boundary check failed."))
  }
  if(!all(sapply(c("max_age", "phi", "psi", "eggs", "probF"), exists, where=boundary_inputs)))
  {
    return(print("The boundary_inputs file must contain all of the following:
                 'max_age', 'phi', 'psi', 'eggs', 'probF'."))
  }
  maxage<- boundary_inputs$max_age
  phi<- boundary_inputs$phi
  psi<- boundary_inputs$psi
  gamma<- boundary_inputs$gamma
  eggs<- boundary_inputs$eggs
  sexratio<- probF
  prod<-boundary_inputs$boundary$product
  if(any(sexratio<=0 | sexratio>1)){sexratio<- sexratio[-which(sexratio<=0 | sexratio>1)]}
  if(length(sexratio)==0)
  {
    return(print("Values of 'probF' greater than 0 and less than or equal to 1
                 must be specified."))
  }
  if(!(all(sapply(c(maxage, gamma, prod), length)==1) &
       is.numeric(c(maxage, gamma, prod))))
  {
    return(print("Boundary inputs 'maxage', 'gamma', and 'boundary$product' should both be numerical values of length 1."))
  }
  if(!all(c(length(psi), length(eggs))==maxage))
  {
    return(print("Boundary inputs 'psi' and 'eggs' need to be a vector of length 'maxage' (one value for each age)."))
  }
  if(length(phi)!=maxage-1)
  {
    return(print("Boundary input 'phi' needs to be a vector of survivals with length maxage-1."))
  }
  curve_dat<- lapply(sexratio, function(r)
  {
    phi1<- c(0.000001,seq(0.01, 1, 0.01))
    phi0<- prod/(gamma*r*phi1)  
    indx<- which(phi0>1)
    if(length(indx)>0)
    {
      phi0<- phi0[-indx]
      phi1<- phi1[-indx]
    }
    #SYMMETRICAL SO CAN EXPAND
    tmp<- phi0
    phi0<- c(phi0, phi1)
    phi1<- c(phi1, tmp)
    out<-NULL
    if(length(phi0>0))
    {
      out<-data.frame(phi1=phi1, phi0=phi0, probF=r, gamma=gamma, 
                      growth_type=boundary_inputs$boundary$growth_type,
                      curve_type="sex ratio")
    }
    return(out)
  })
  curve_dat<- do.call(rbind, curve_dat)
  return(curve_dat)
}


phi0_phi1_atresia_curves<- function(boundary_inputs=NULL,
                                    atresia_prop=c(0.1,1,0.1))
{
  # ERROR HANDLING
  if(!all(sapply(c("product", "check2"), exists, where=boundary_inputs$boundary)))
  {
    return(print("Boundary inputs must contain the list boundary with entries 'product' and 'check2'."))
  }
  if(boundary_inputs$boundary$check2!=0)
  {
    return(print("Boundary check failed."))
  }
  if(!all(sapply(c("max_age", "phi", "psi", "eggs","gamma", "probF"), exists, where=boundary_inputs)))
  {
    return(print("The boundary_inputs file must contain all of the following:
                 'max_age', 'phi', 'psi', 'eggs', 'gamma', 'probF'."))
  }
  maxage<- boundary_inputs$max_age
  phi<- boundary_inputs$phi
  psi<- boundary_inputs$psi
  gamma<- 1-atresia_prop
  eggs<- boundary_inputs$eggs
  probF<- boundary_inputs$probF
  prod<-boundary_inputs$boundary$product
  if(any(gamma<=0 | gamma>1)){gamma<- gamma[-which(gamma<=0 | gamma>1)]}
  if(length(gamma)==0)
  {
    return(print("Values of 'atresia_prop' greater than 0 and less than or equal to 1
                 must be specified."))
  }
  if(!(all(sapply(c(maxage, probF, prod), length)==1) &
       is.numeric(c(maxage, probF, prod))))
  {
    return(print("Boundary inputs 'maxage', 'probF', and 'boundary$product' should both be numerical values of length 1."))
  }
  if(!all(c(length(psi), length(eggs))==maxage))
  {
    return(print("Boundary inputs 'psi' and 'eggs' need to be a vector of length 'maxage' (one value for each age)."))
  }
  if(length(phi)!=maxage-1)
  {
    return(print("Boundary input 'phi' needs to be a vector of survivals with length maxage-1."))
  }
  curve_dat<- lapply(gamma, function(g)
  {
    phi1<- c(0.000001,seq(0.01, 1, 0.01))
    phi0<- prod/(g*probF*phi1)
    indx<- which(phi0>1)
    if(length(indx)>0)
    {
      phi0<- phi0[-indx]
      phi1<- phi1[-indx]
    }
    #SYMMETRICAL SO CAN EXPAND
    tmp<- phi0
    phi0<- c(phi0, phi1)
    phi1<- c(phi1, tmp)
    out<-NULL
    if(length(phi0>0))
    {
      out<-data.frame(phi1=phi1, phi0=phi0, gamma=g, probF=probF, 
                      growth_type=boundary_inputs$boundary$growth_type,
                      curve_type="atresia")
    }
    return(out)
  })
  curve_dat<- do.call(rbind, curve_dat)
  return(curve_dat)
}

## PLOT BOUNDARY CURVES IN SURVIVAL-SEXRATIO SPACE
plot_boundary_curves<- function(curve_dat=NULL,
                                curve_type="atresia",
                                phi0_upper=0.001,
                                phi1_upper=0.8,
                                xlabel=expression(paste("Age-1 Survival Probability  (", phi[1], ")")),
                                ylabel=expression(paste("Age-0 Survival Probability  (  ", phi[0], ")")),
                                xaxis="s")
{
  curve_dat<- subset(curve_dat, curve_type==curve_type)
  if(nrow(curve_dat==0))
  {
    return(print("No curve data of the specified type available."))
  }
  gt<- unique(curve_dat$growth_type)
  cdat<- subset(curve_dat, growth_type==gt[1])
  
  probF<-unique(cdat$probF)
  probF<- probF[order(probF)]
  tmp<-cdat[which(cdat$probF==probF[1]),]
  tmp<-tmp[order(tmp$phi1),]
  plot(tmp$phi1, tmp$phi0, type="l", 
        xlim=c(0,phi1_upper), ylim=c(0,phi0_upper),
        xlab=xlabel, ylab="", xaxt=xaxis,
        tck=0.02, mgp=c(1.5,0.1,0), las=1)
  mtext(ylabel, 2, outer=TRUE, padj=2)
  if(length(probF>1))
  {
    invisible(lapply(2:length(probF), function(i)
    {
      tmp<-cdat[which(cdat$probF==probF[i]),]
      tmp<-tmp[order(tmp$phi1),]
      points(tmp$phi1, tmp$phi0, type="l")
    })) 
  }
  if(length(gt)>1)
  {
    cdat<- subset(curve_dat, growth_type==gt[2])
    probF<-unique(cdat$probF)
    probF<- probF[order(probF)]
    invisible(lapply(1:length(probF), function(i)
    {
      tmp<-cdat[which(cdat$probF==probF[i]),]
      tmp<-tmp[order(tmp$phi1),]
      points(tmp$phi1, tmp$phi0, type="l", lty=2)
    })) 
    g_lab<- ifelse(gt=="vbgf", "JAGS VBGF", "TMB LN(VBGF)")
    legend("topright", g_lab, lty=c(1,2), bty='n')
  }
}



# 7. 
## INITIAL POPULATION
initialize_pop<- function(inputs=NULL,
                          gamma_hist=1,
                          phi0_MR_hist=0.000186,
                          stable_age=FALSE,
                          boom_prob=1/5,
                          bust_prob=1/5,
                          boom_recruits=25000,
                          bust_recruits=0,
                          avg_recruits=1000,
                          initial_adults=5000,
                          exact=TRUE)
{
  # ERROR CHECK
  if(stable_age & (length(initial_adults)!=1 | initial_adults<=0))
  {
    return(print("For a 'stable_age' generation of the initial population, a single
                  positive value for 'initial_adults' must be specified."))
  }
  if(!stable_age)
  {
    if(!all(sapply(c(boom_prob, bust_prob, boom_recruits, bust_recruits,
                     avg_recruits), length)==1))
    {
      return(print("For a random generation of the initial population, a single
                   positive value for each of 'boom_prob', 'bust_prob', 
                   'boom_recruits', 'bust_recruits', 'avg_recruits' must be 
                   specified."))
    }
    if(boom_prob+bust_prob>1)
    {
      return(print("The sum of 'boom_prob' and 'bust_prob' must be between 0 and 1."))
    }
  }
  ## RECORD NEW INPUTS
  inits<- list()
  inits$inputs_id<- inputs$id
  inits$gamma_historical<- gamma_hist
  inits$phi0_historical<- phi0_MR_hist
  inits$stable_age<- stable_age
  inits$exact<- exact
  if(!stable_age)
  {
    inits$boom_prob<- boom_prob
    inits$bust_prob<- bust_prob
    inits$boom_recruits<- boom_recruits
    inits$bust_recruits<- bust_recruits
    inits$avg_recruits<- avg_recruits
  }
  if(stable_age)
  {
    inits$initial_adults<- initial_adults
  }
  ## CREATE INITIAL AGE DISTRIBUTION BY PROJECTIONS OF RANDOM BOOM/BUST YEARS
  if(!stable_age)
  {
    year_type<- rmultinom(inputs$max_age, 1, 
                          c(bust_prob, 1-bust_prob-boom_prob, boom_prob))
    vals<- c(bust_recruits, avg_recruits, boom_recruits)
    N0<- c(vals%*%year_type)
    for(i in 2:60)
    {
      N0[i]<- N0[i]*prod(inputs$phi[1:(i-1)])
    }
    if(!exact)
    {
      N0<- round(N0)
    }
  }
  ## CREATE INITIAL AGE DISTRIBUTION AT ASSUMED HISTORICAL STEADY AGE DISTRIBUTION
  if(stable_age)
  {
    # BUILD HISTORICAL LESLIE MATRIX
    A<- matrix(0,inputs$max_age,inputs$max_age)
    ## SURVIVAL VALUES
    A[cbind(2:inputs$max_age,1:(inputs$max_age-1))]<- inputs$phi
    ## FERTILITY VALUES
    A[1,] <- inputs$psi*gamma_hist*inputs$eggs*inputs$probF*phi0_MR_hist 
    # HISTORICAL EIGENANALYSIS
    ea<- eigen.analysis(A)
    initial_num<- initial_adults/sum(ea$stable.age[inputs$mat$a_h:inputs$max_age])
    if(exact)
    {
      N0<-ea$stable.age*initial_num
    }
    if(!exact)
    {
      N0<- c(rmultinom(1, initial_num-initial_adults, 
                       ea$stable.age[1:(inputs$mat$a_h-1)]),
             rmultinom(1, initial_adults, 
                       ea$stable.age[inputs$mat$a_h:inputs$max_age]))
    }
  }
  inits$N0<- N0
  return(inits)
}  


reproductive_value<- function(pop_data=NULL)
{
  
}

stable_age0<- function(pop_data=NULL)
{
  
}


