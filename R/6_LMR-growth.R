#----------------------------------------------------------------------
# 
# IMPORT DATA FOR STOCKED FISH
#
#----------------------------------------------------------------------
stocked<-fread("../PSPAP-data/stocking-processed.csv")
stocked<-stocked[,setdate:=set.date]


#----------------------------------------------------------------------
# 
# read in bend data 
#
#----------------------------------------------------------------------
bends<-fread("../PSPAP-data/bend-data.csv")
bends[BRM_ID==558,]$UPPER_RIVER_MILE<- 753.1
bends[BRM_ID==558,]$LOWER_RIVER_MILE<- 753.0
bends[,id:=.I]

# add segment_id to stocking location
stocked[,segment_id:=NA]
stocked<-split(stocked,by="rpma")
stocked<-lapply(stocked,function(x)
{
  rpma_i<- unique(x$rpma)
  if(rpma_i !=3)
  {
    # make a function to assign bend id
    bend_id<-approxfun(bends[RPMA==rpma_i,]$LOWER_RIVER_MILE,
                       bends[RPMA==rpma_i,]$id,
                       method="constant")
    x[,id:= bend_id(rm)]
  }
  if(rpma_i==3)
  {
    x[,id:=1]
  }   
  return(x)
})
stocked<-rbindlist(stocked)
stocked<-merge(stocked,
               bends[,.SD, .SDcols=c("id","B_SEGMENT")],
               by="id",all.x=TRUE)
setDT(stocked)[rpma==3 ,B_SEGMENT:=5]
setDT(stocked)[rpma==2 & rm>=1761.1 ,B_SEGMENT:=1]
stocked[,segment_id:=B_SEGMENT]
  # NAs IN ID ARE FOR LACK OF RM AND YE AND BH RPMA 2 STOCKING


#----------------------------------------------------------------------
# 
# IMPORT PSPAP AND HAMP DATA [need hamp data added to query]
#
#----------------------------------------------------------------------
dat<-fread("../PSPAP-data/fish-processed.csv")
dat<-subset(dat, SPECIES=="PDSG")
names(dat)<- tolower(names(dat))
dat[,tmp:=1] # FOR COUNTING
# MAKE ALL PIT TAG CHARACTERS UPPER CASE
dat[,tag_number:=toupper(tag_number)]
dat[,origin:=ifelse(origin=="H",1,0)] ## take max to assign origin for known hatchery fish
# add rpma
dat[,rpma:=4]
# set rpma for segments 5 and 6
setDT(dat)[segment_id%in%c(5,6,23),rpma:=3]
# set rpma for segment 1-4
setDT(dat)[segment_id%in%c(1,2,3,4,21,22,51,52,53,54,55),rpma:=2]# 21 is Milk, 22 is YSR
dat[,length_type:=1] # 1 for recovery in system
dat[,setdate:=as.POSIXct(setdate)]

flds<-c("segment_id","tag_number","rpma","setdate",
        "species","length","origin","length_type")
dat<-dat[,.SD,.SDcols=flds]

spawn<- stocked[,.SD, .SDcols=c("tag_number","spawn.date")]
  ## TO OBTAIN ENCODED PART
  # tgs<- unique(spawn$tag_number)
  # tgs<- tgs[order(tgs)]
  # tail(tgs)
  # rm(tgs)
spawn<-subset(spawn,!(tag_number%in%c("â€¦", "..........","")))
spawn<- spawn[-which(duplicated(spawn)),]
  # NOTE THERE ARE SOME STOCKING DATA ISSUES,
  # INCLUDING THESE REPEATS (SOME REPEATS EVEN WITH DIFFERENT PARENTAGE) 
dat<- merge(dat, spawn, 
            all.x=TRUE,
            by="tag_number")

flds<- names(dat)
stock<- stocked[,.SD,.SDcols=flds]
dat<- rbind(dat, stock)


# KNOWN AGE FISH
datA<- subset(dat, !(is.na(spawn.date)))
# calculate age in years
datA[,age:=as.numeric(difftime(setdate,spawn.date,units="days"))/365.25]
## KNOWN AGE FISH WITH LENGTH DATA
L_at_A<- subset(datA, !(is.na(length)) & length!=0)
# restrict to LMR data
L_at_A<- subset(L_at_A, rpma==4)

rm(bends, datA, spawn, stock, stocked, flds)

# UNKNOWN AGE FISH
datUA<- subset(dat, is.na(spawn.date))
datUA<- subset(datUA,!(tag_number%in%c("â€¦", "..........","")))
# drop missing lengths
datUA<- subset(datUA, !(is.na(length)) & length!=0)
## RESTRICT TO LMR
### NOTE THIS REMOVES THE LENGTH DATA OF FISH CAPTURED IN RPMA 3
### THEN RECAPTURED IN RPMA 4 (BUT WILL KEEP MULTIPLE RPMA 4 RECAPTURE LENGTHS)
datUA<-subset(datUA, rpma==4)
## GROWTH DATA OF UNKNOWN AGE FISH
tags<- unique(datUA$tag_number[which(duplicated(datUA$tag_number))])
Gdat<- subset(datUA, tag_number %in% tags)
Gdat<- Gdat[order(tag_number,setdate),]
# add an id for each sequential observation of a fish
Gdat[, length_id := rowid(tag_number)] 
# time at large in years
Gdat[,tal:=difftime(setdate,shift(setdate),units="days"),by="tag_number"]
Gdat[,tal:=as.numeric(tal)/365.25]
# set tal to 0 for 1st captures
setDT(Gdat)[is.na(tal),tal:=0] 
# subset out fish captured again on same day  
Gdat<- Gdat[!(length_id>1& tal==0),] 

# update occasion id
Gdat<- Gdat[order(tag_number,setdate),]
Gdat[, length_id := rowid(tag_number)] 
# change in length
Gdat[,deltaL:=c(0,diff(length))]
Gdat[,deltaL:=ifelse(length_id==1, 0, deltaL)]
length(which(Gdat$deltaL<0)) #205
tags<- unique(Gdat[deltaL<0,]$tag_number) #190
# look at fish with 9 recaps
Gdat[which(Gdat$tag_number=="4626581E0E"),]
# eventually add in a tolerance, but for now, let's call this 
# fish constant over the years
Gdat[which(Gdat$tag_number=="4626581E0E"),]$length<- 
  round(mean(Gdat[which(Gdat$tag_number=="4626581E0E"),]$length))

# RERUN CHANGE IN LENGTH
Gdat[,deltaL:=c(0,diff(length))]
Gdat[,deltaL:=ifelse(length_id==1, 0, deltaL)]
length(which(Gdat$deltaL<0)) #202
tags<- unique(Gdat[deltaL<0,]$tag_number) #189
# remove all remaining fish with ANY negative change in lengths
Gdat<- subset(Gdat, !(tag_number %in% tags))
# NOTE A HANDFUL OF FISH ON THIS LIST HAVE DIFFERING ORIGINS
UALdat<- dcast(Gdat,tag_number~length_id, value.var="length",
               mean, fill=-1)
UALdat[UALdat==-1]<-NA
UALdat$ncap<- sapply(1:nrow(UALdat), function(x){ncol(UALdat)-1-length(which(is.na(UALdat[x,])))})

deltaT<- dcast(Gdat, tag_number~length_id, value.var="tal",
               mean, fill=-1)
deltaT[deltaT==-1]<-NA

rm(dat, datUA, tags, Gdat)


## SET UP MODEL
vbgf_known_and_unknown_age <- function()
{
  #  vbgf model
  for(i in 1:N_known)
  {
    La[i]<- Linf*(1-exp(-k*(Y[i,1]-t0)))
    Y[i,2]~dnorm(La[i], tau[i])
    # precision decreases with age f(delta)
    tau[i] <- 1/((sigma^2)*Y[i,1]^(2*delta))
  }
  
  # fabens model
  for(i in 1:N_unknown)
  {
    for(j in 2:n_cap[i])
    {
      L[i,j]~dnorm(L_hat[i,j],1/(sigma_uknown_age^2))
      L_hat[i,j]<-L[i,j-1]+ ((Linf-L[i,j-1])*(1-exp(-k*dY[i,j])))
    }
  }   
  
  # priors
  t0~dnorm(0,0.0001)
  # growth coefficient k
  lnk~dnorm(0,0.0001)
  log(k)<-lnk
  # L_infinity
  Linf~dunif(0,1800) #MAXIMUM PSPAP PS LENGTH: 1640  
  sigma ~ dgamma(0.001, 0.001)
  # scale precision with age
  delta ~ dnorm(0,1/1000) 
  sigma_uknown_age ~ dgamma(0.001, 0.001)	
}

# FIT MODEL
## GENERATE INPUT DATA FOR JAGS
LMR<- list()
LMR$N_known<- nrow(L_at_A)
LMR$Y<- as.matrix(L_at_A[,.SD, .SDcols=c("age","length")])
LMR$N_unknown<- nrow(UALdat)
LMR$n_cap<- UALdat$ncap
LMR$L<- as.matrix(UALdat[,2:(ncol(UALdat)-1)])
LMR$dY<- as.matrix(deltaT[,2:ncol(deltaT)])
## INITIAL VALUES
inits<- function(t)
{	
  list(
    t0=0,
    lnk=0,
    Linf=max(c(LMR$Y[,2],max(LMR$L,na.rm = TRUE)))+10,
    sigma=0.1,
    delta=0,
    sigma_uknown_age=0.1)
}
## PARAMETERS TO FIT
params<- c("t0","Linf","k","sigma","delta","sigma_uknown_age")	
## RUN MODEL FIT & SAVE
ptm<-proc.time()
fit <- jags.parallel(data=LMR,
                     inits=inits,
                     parameters=params,	
                     model.file=vbgf_known_and_unknown_age,
                     n.chains = 4,	
                     n.iter = 200000,	
                     n.burnin = 100000,
                     export_obj_names=c("LMR"),
                     n.thin=1,
                     working.directory=getwd())
tot<-(proc.time()-ptm)[3]/60
out<- list(fit=fit,dat=LMR,mod=vbgf_known_and_unknown_age,
             params=params,tot=tot)
saveRDS(out, "_output/vbgf-known-and-unknown-age.RDS")

