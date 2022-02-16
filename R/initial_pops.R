# GENERATE INITIAL POPULATION
init_pop<- function(inputs=NULL,
                    type=NULL) # "Uniform" OR "2020_PSPAP"
{
  inps<- inputs
  ## 2020 AMCR MEDIAN POP EST (APPENDIX TABLE 3-7)
  inps$N_H<- 35857+20582
  inps$N_W<- 6427+2527
  inps$sexratio_H<- 0.5
  inps$sexratio_W<- 0.35
  # SOURCE: THE 63:115 F:M RATIO REPORTED IN STEFFENSEN ET AL. (2013)
  # COULD ALSO USE 0.4708 FROM Wildhaber et al. 2017
  inps$N0_type<- type
  if(type=="Uniform")
  {
    N<- rbinom(1, inps$N_H, inps$sexratio_H) + 
      rbinom(1, inps$N_W, inps$sexratio_W)
    N0<- rep(floor(N/inps$max_age), inps$max_age)
    N0<- N0+c(rep(1, N-sum(N0)), rep(0, inps$max_age-N+sum(N0)))
    inps$N0<- N0
  }
  if(type=="2020_PSPAP")
  {
    f_dat<- fread("../PSPAP-data/fish-processed.csv")
    names(f_dat)<-tolower(names(f_dat))
    f_dat[,rpma:=4]
    # set rpma for segments 5 and 6
    setDT(f_dat)[segment_id%in%c(5,6,23,56),rpma:=3]
    # set rpma for segment 1-4
    setDT(f_dat)[segment_id%in%c(1,2,3,4,21,22,51,52,53,54,55),rpma:=2]# 21 is Milk, 22 is YSR
    # SUBSET OUT RPMA 4 FISH IN MOST RECENT YEARS
    f_dat<- subset(f_dat, rpma==4 & year>=2018)
    # READ IN STOCKING DATA
    stocked<- fread("../PSPAP-data/stocking-processed.csv")
    # REMOVE STOCKINGS WITH MISSING SPAWN DATE
    stocked[,spawn_date:= spawn.date]
    stocked<- subset(stocked, !is.na(spawn_date))
    # PULL RELEVANT COLUMNTS
    indx<- match(c("tag_number","spawn_date"), names(stocked))
    spawn<- stocked[,..indx]
    # REMOVE STOCKINGS WITH MISSING TAG NUMBERS
    tags<- unique(spawn$tag_number)
    tags<- tags[order(tags)]
    setDT(spawn)[tag_number %in% c("",tags[length(tags)], "..........","...",".",
                                   "0000000000","N/A","NOFISHSCAN","unknown",
                                   "XXXXXXXXXX","CAUDAL CL","N PIT"),
                 tag_number:="-99"]
    setDT(spawn)[is.na(tag_number),tag_number:="-99"]    
    setDT(spawn)[tag_number=="NA",tag_number:="-99"] 
    spawn<- subset(spawn, tag_number!="-99")
    # REMOVE DUPLICATE DATA
    spawn<- spawn[!duplicated(spawn),]
    # ADD SPAWNING DATE TO FISH DATA
    f_dat<- merge(f_dat, spawn, all.x=TRUE)
    # SUBSET OUT UNKNOWN AGE FISH
    f_dat<- subset(f_dat, !is.na(spawn_date))
    # KEEP MOST RECENT CAPTURE OF KNOWN AGE FISH
    tags<- f_dat[which(duplicated(f_dat$tag_number)),]$tag_number
    f_ids<- f_dat[which(f_dat$tag_number %in% tags),]$f_id
    indx<- sapply(tags, function(x)
    {
      tmp<- f_dat[which(f_dat$tag_number==x),]
      return(tmp[which.max(tmp$setdate),]$f_id)
      
    })
    f_ids<- setdiff(f_ids, indx)
    f_dat<- f_dat[-which(f_dat$f_id %in% f_ids),]
    rm(f_ids, tags, indx)
    # ADD IN AGE AT CAPTURE
    f_dat$age_capture<- floor(as.numeric(difftime(f_dat$setdate, 
                                                  f_dat$spawn_date,
                                                  units="days")/365.25))
    # COMPUTE AGE AT END OF 2020
    f_dat$age_2020<- as.numeric(difftime(as.POSIXct("2020-12-31"), 
                                         f_dat$spawn_date,
                                         units="days")/365.25)
    # AGGREGATE TO PROPORTIONS
    f_dat$age<- floor(f_dat$age_2020)
    f_dat$yr_diff<- f_dat$age-f_dat$age_capture
    f_dat$freq<- sapply(1:nrow(f_dat), function(i)
    {
      ifelse(f_dat$yr_diff[i]==0, 1,
             prod(inps$phi[f_dat$age_capture[i]:
                             (f_dat$age_capture[i] + f_dat$yr_diff[i]-1)]))
    })
    est<- aggregate(freq~age, f_dat, sum)
    est$prop<- est$freq/sum(est$freq)
    N0_i<- rbinom(inps$N_H, inps$sexratio_H)
    est$N0<- rmultinom(1, N0_i, est$prop)
    est<- merge(est, data.frame(age=1:inps$max_age), all=TRUE)
    est<- est[order(est$age),]
    # NO STOCKING IN 2020; ALL 2018 STOCKING SPAWNED IN 2017
    # ALL 2019 STOCKINGS FROM 2019 SPAWNING WERE DRIFT STUDY FISH
    est$N0[1]<- round(sum(stocked[which(stocked$year.stocked==2019 &
                                          as.numeric(format(stocked$spawn_date, "%Y"))==2019),]$numbers.stocked)*inps$phi0_MR*1)
    est$N0[2]<- round(sum(stocked[which(stocked$year.stocked==2019 &
                                          as.numeric(format(stocked$spawn_date, "%Y"))==2018),]$numbers.stocked)*prod(inps$phi[1]))
    est[is.na(est)]<- 0
    v<- inps$max_age-68+1
    est[which(est$age %in% 68:inps$max_age),]$N0<- rmultinom(1, ceiling(inps$N_W*inps$sexratio_W), rep(1/v, v))
    inps$N0<- est$N0
  }
  return(inps)
}