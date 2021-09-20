#----------------------------------------------------------------------
# 
# IMPORT DATA FOR STOCKED FISH
#
#----------------------------------------------------------------------
stocked<-fread("../PSPAP-data/stocking-processed.csv")
stocked<-stocked[,setdate:=set.date]
# REMOVE LATER TO KEEP LENGTH AT AGE DATA
stocked<-subset(stocked,!(tag_number%in%c("â???¦",
                                          "..........","")))

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
fish<-dat
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
dat<-rbind(dat[,.SD,.SDcols=flds],
           stocked[,.SD,.SDcols=flds])

# ADD IN SPAWN DATE OF KNOWN HATCHERY ORIGIN FISH
# spawn<- stocked[,.SD, .SDcols=c("tag_number","spawn.date","family_lot")]
#   ## # CURRENTLY NOT NEEDED
#   ## spawn<-subset(spawn, !(is.na(tag_number))) 
#   ## spawn<-subset(spawn, tag_number!="")
# # SEE 
# # stocked[which(stocked$tag_number %in% c("6C00090769", "49002A6629")),]
# # FOR SOME STOCKING DATA ISSUES
# spawn<- spawn[-which(duplicated(spawn)),]
# # SEE ADDITIONAL ISSUES HERE; HANDFUL OF FISH WITH 2 PARENTAGES!
# # test<- spawn[c(which(duplicated(spawn$tag_number)),
# #                which(duplicated(spawn$tag_number, fromLast = TRUE))),]
# # test<- test[order(test$tag_number),]
# # test
## FAMILY LOT NOT NEEDED FOR THIS ANALYSIS SO DROP AND REMOVE REPEATS
spawn<- stocked[,.SD, .SDcols=c("tag_number","spawn.date")]
spawn<- spawn[-which(duplicated(spawn)),]
dat<- merge(dat, spawn, 
            all.x=TRUE,
            by="tag_number")
# KNOWN AGE FISH
datA<- subset(dat, !(is.na(spawn.date)))
# calculate age in years
datA[,age:=as.numeric(difftime(setdate,spawn.date,units="days"))/365.25]
## KNOWN AGE FISH WITH LENGTH DATA
L_at_A<- subset(datA, !(is.na(length)) & length!=0)
