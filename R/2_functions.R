
create_inputs<- function(max_age=41,
                         phi0= 0.000075, #AGE-0 (EGG TO AGE-1) SURVIVAL
                         phi_young=c(0.686), #VECTOR OF SURVIVALS FOR AGE-1 TO AGE-N WHERE N IS LENGTH OF VECTOR
                         phi_old=0.92, #SURVIVAL OF AGE-N+1 PLUS FISH
                         growth_model="VBG", #ACCEPTS "VBG", "Gompertz", "Logistic", "Power"
                         growth_params=c(exp(6.9750323), exp(-2.4073801), -2.1603259, 0), #Linf, k, t0, error for all except "Power" where these are a, b, c, error for ax^b+c
                         fec_params=c(),
                         probF=0.5, # PROBABILITY A FISH IS FEMALE (SEX RATIO)
                         
                        )
{
  inps<- list()
  # MAXIMUM AGE
  inps$max_age<- max_age
  # SURVIVALS
  inps$phi<- c(phi_young, rep(phi_old, max_age-1-length(phi_young)))
  inps$phi0<- phi0
  # FERTILITIES
  ## GROWTH FUNCTION (LENGTH-AT-AGE)
  if(growth_model=="VBG")
  {
    L<- growth_params[1]*(1-exp(-1*growth_params[2]*(age-growth_params[3])))
  }
  if(growth_model=="Gompertz")
  {
    #L<- growth_params[1]*(1-exp(-1*growth_params[2]*(age-growth_params[3])))
  }
  if(growth_model=="Logistic")
  {
    #L<- growth_params[1]*(1-exp(-1*growth_params[2]*(age-growth_params[3])))
  }
  if(growth_model=="Power")
  {
    #L<- growth_params[1]*(1-exp(-1*growth_params[2]*(age-growth_params[3])))
  }
  ## FECUNDITY FUNCTION (EGGS-AT-LENGTH) 
  E<- fec_params[1]+fec_params[2]*L 
  ## SEX RATIO
  inps$probF<- probF
}

