###################################################################
# EFDM core
# December 12, 2014
# Authors: Seija Sirkiä, Maria Eronen, Tuula Packalen
# Original Author: Seija Sirkiä, January 25th 2013
#   
#   
#   Copyright 2014 European Union
#   
#   Licensed under the EUPL, Version 1.1 or – as soon they
#   will be approved by the European Commission - subsequent
#   versions of the EUPL (the "Licence");
#   You may not use this work except in compliance with the
#   Licence.
#   You may obtain a copy of the Licence at:
#   
#   http://joinup.ec.europa.eu/software/page/eupl/licence-eupl
#   
#   Unless required by applicable law or agreed to in
#   writing, software distributed under the Licence is
#   distributed on an "AS IS" basis,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
#   express or implied.
#   See the Licence for the specific language governing
#   permissions and limitations under the Licence.
#   
###################################################################
# functions needed:

library(abind) #for function abind

readtransprobs<-function(filename)
  #for reading a ready made transition matrix or several from
  #a text file (file extension .txt), or from an R save file 
  #(file extension .RData)
{
  if(tail(strsplit(filename,"[.]")[[1]],1)=="RData")
  {
    tmp<-load(filename) #load should create a single new object
    #in this local environment, it still needs to be returned
    return(get(tmp))
    #if there were several, if it is the wrong shape...
    #who knows what will happen!
  }
  #else
  rawdata<-readLines(filename)
  if(length(rawdata)==basedim) # if so, 
    #there's just one transition matrix used everywhere
  {
    tmp<-matrix(as.numeric(unlist(strsplit(rawdata,split=" "))),
                ncol=basedim,byrow=TRUE)
    transmat<-rep(tmp,times=prod(otherdim))
    dim(transmat)<-c(basedim,basedim,otherdim)
    return(transmat)
  }
  #otherwise, there are several matrices, for different factor levels
  
  #Which factors are involved? This should be mentioned on the first line
  tmp<-unlist(strsplit(rawdata[1],split=" "))
  tmp<-unlist(strsplit(tmp,split="="))
  fs<-tmp[seq(1,length(tmp),by=2)] #first, third etc are the factor names
  fn<-factdims[fs]
  #rest of the file is the matrices, separated by "header" lines
  frows<-seq(from=1,length=prod(fn),by=basedim+1) #these are the header lines
  tmp<-matrix(as.numeric(unlist(strsplit(rawdata[-frows],split=" "))),
              ncol=basedim,byrow=TRUE)
  #Actually, the header lines are thrown away! It doesn't matter what's on them,
  #but they probably help keeping the file organized as assumed
  dim(tmp)<-c(basedim,fn,basedim)
  tmp<-aperm(tmp,c(1,length(fn)+2,1:length(fn)+1))
  os<-setdiff(otherdimnames,fs)
  transmat<-rep(tmp,times=prod(factdims[os]))
  dim(transmat)<-c(basedim,basedim,fn,factdims[os])
  nr<-1:length(otherdimnames)
  names(nr)<-c(fs,os)
  transmat<-aperm(transmat,c(1,2,nr[otherdimnames]+2))
  return(transmat)
}

mmf<-function(M)
  #simply for doing the necessary "matrix times a vector" multiplications
  #somewhat efficiently
{
  m<-dim(M)[1]
  M[,1:m]%*%M[,m+1]
}

dividebyA<-function(oldstate)
{
  result<-vector("list",nact)
  #if dynamic updating of A, use a function to get a new A:
  #actproblist<<-do.call(updateA,actproblist)
  
  #while(TRUE) { #if iterating A
  for(i in 1:nact) {
    nr<-1:nfact
    names(nr)<-factnames    
    tmp<-oldstate*actproblist[[i]]
    #tmp is the subpopulation receiving the activity i
    result[[i]]<-tmp
  } 
  #if(!iterating) break #if no iteration, we're done here
  #else {
  #tmp<-do.call(iterateA,args=list(result,actproblist))
  #if(is.null(tmp)) break 
  #else if(usethesameAonnextstep) actproblist<<-tmp
  #else actproblist<-tmp #only change the local A
  #} #else, we iterated
  #} while
  names(result)<-actnames
  result
}


newstate<-function(oldstatediv) 
  #for producing the next state from the current one, using the obtained
  #activity and transition probabilities
{
  result<-oldstatediv
  for(i in 1:nact) {
    nr<-1:nfact
    names(nr)<-factnames
    acti<-activities[[i]]
    factnames.here<-c(acti,setdiff(factnames,acti))
    tmp<-aperm(oldstatediv[[actnames[i]]],nr[factnames.here])
    #tmp is the subpopulation receiving the activity i, with appropriate dims
    dim(tmp)<-c(prod(dim(tmp)[1:length(acti)]),1,prod(dim(tmp)[-(1:length(acti))]))
    transmati<-transmats[[i]]
    #transmati is the transition matrix corresponding to activity i
    dim(transmati)<-c(dim(transmati)[1:2],prod(dim(transmati)[-(1:2)]))
    tmp<-apply(abind(transmati,tmp,along=2),3,mmf)
    #tmp now has the new states of this subpopulation
    dim(tmp)<-factdims[factnames.here]
    names(nr)<-factnames.here
    tmp<-aperm(tmp,nr[factnames])
    dimnames(tmp)<-dimnames(tmp)
    result[[actnames[i]]]<-tmp
  }

  result
}


vz <- "C:/Users/regelmann/Documents/GitHub/efdm/toyDataSet/Germany"
setwd(vz)
getwd()

##############################################
# Scenario OLD ###############################
##############################################

#describing the factors:
factors.tmp<-readLines("DEfactors.txt")
nfact<-length(factors.tmp)
factors.tmp<-strsplit(factors.tmp," ")
factnames<-character(nfact)
for(i in 1:nfact) {
  factnames[i]<-factors.tmp[[i]][1]
  factors.tmp[[i]]<-factors.tmp[[i]][-1]
}
names(factors.tmp)<-factnames
factlvls<-factors.tmp
factdims<-sapply(factlvls,length)


#describing the activities
activities.tmp<-readLines("DEactivities.txt")
nact<-length(activities.tmp)
activities<-strsplit(activities.tmp," ")
actnames<-character(nact)
actmethod<-character(nact)
actfiles<-character(nact)
rm(activities.tmp)

for(i in 1:nact) {
  actnames[i]<-activities[[i]][1]
  actmethod[i]<-activities[[i]][2]
  actfiles[i]<-activities[[i]][3]
  activities[[i]]<-activities[[i]][-(1:3)]
}
names(activities)<-actnames

#setting up the transition matrices
transmats<-list()
for(i in 1:length(activities))
{
  #the levels of these factors may change under this activity:
  basedimnames<-activities[[i]] 
  #these are fixed under this activity:
  otherdimnames<-setdiff(factnames,basedimnames)
  basedim<-prod(factdims[basedimnames])
  otherdim<-factdims[otherdimnames]
  
  transmat.forthisround <- readtransprobs(actfiles[i]) 
  
  #the rest is making sure dimensions have correct and sensible names
  lvlnames<-factlvls[[basedimnames[1]]]
  if(length(basedimnames>1))
    for(j in 2:length(basedimnames))
    {
      lvlnames.more<-factlvls[[basedimnames[j]]]
      lvlnames<-paste(rep(lvlnames,times=length(lvlnames.more)),
                      rep(lvlnames.more,each=length(lvlnames)),sep=".")
    }
  dimnames(transmat.forthisround)<-c(list(lvlnames,lvlnames),factlvls[otherdimnames])
  transmats<-c(transmats,list(transmat.forthisround))
  names(transmats)[i]<-actnames[i]
}

rm(basedim,basedimnames,otherdim,otherdimnames,transmat.forthisround,lvlnames,
   lvlnames.more)

#setting up the initial state and the probabilities of activities
statespace<-read.table("DEinitstate.txt",header=TRUE)
actprobtable<-read.table("DEactprobs_nothin.txt",header=TRUE) #neues Einfügen
#actprobtable<-read.table("actprobtable_new.txt",header=TRUE) #neues Einfügen
#this next merge is a bit dumb but the fixing of factor levels and ordering
#needs to be done for both, so I'm doing it at once
statespace<-merge(statespace,actprobtable)
for(fname in factnames)
  statespace[[fname]]<-factor(statespace[[fname]],
                              levels=factlvls[[fname]],ordered=TRUE)
statespace<-statespace[do.call(order,statespace[rev(factnames)]),]

initstate<-statespace[[setdiff(names(statespace),union(factnames,actnames))]]
dim(initstate)<-factdims
dimnames(initstate)<-factlvls


actproblist<-list()
for(i in 1:nact) {
  actproblist[[i]]<-array(statespace[[actnames[i]]],dim=factdims)
  dimnames(actproblist[[i]])<-factlvls
}

#-----------------------------
#finally, the actual simulation

#resultstates<-statespace[factnames]

state<-dividebyA(initstate)
resultstates<-list(step0=state)


#for(i in 1:nact)
#resultstates[[actnames[i]]]<-data.frame(step0=as.vector(tmp[[i]]))

nrofsteps<-as.numeric(20) #how many steps should be taken?
#a silly compatibility issue there with the names, sorry about that

for (i in 1:nrofsteps) {
  nstate<-newstate(state)  
  state<-rowSums(do.call(cbind,lapply(nstate,as.vector)))
  dim(state)<-factdims
  state<-dividebyA(state)
  resultstates[[paste0("step",i)]]<-state
}
#resultstates is now a list of lists, outer list going over the simulation steps, 
#and inner lists going over the activities. The inmost elements are multiway arrays
#in the form of the statespace (that is, factdims)

#Reformulate resultstates to be data frame whose columns are arrays, with nact columns.
#First, vectorize and cbind the multiway arrays:
resultstates<-lapply(resultstates,function(elem) {do.call(cbind,lapply(elem,as.vector))})
#Second, make the data frame:
resultstates<-data.frame(resultstates)

rawoutput<-statespace[factnames]
rawoutput<-data.frame(rawoutput,resultstates)
#rawoutput has 'nact' columns per step, with names suchs as step0.actname
write.table(rawoutput,file="rawoutput.txt")

tmp<-as.matrix(resultstates)
dim(tmp)<-c(dim(tmp)[1],nact,nrofsteps+1)
tmp<-apply(tmp,c(1,3),FUN=sum)
resultstates<-data.frame(tmp)
names(resultstates)<-paste0("step",0:nrofsteps)
resultstates<-data.frame(statespace[factnames],resultstates)
#resultstates is like rawoutput, but with the different columns with name
#such as step0.something summed together and named just step0
write.table(resultstates,file="resultstates.txt")


vol_class <- list(Spruce=c(0,60,153,270,391,503,601,682,748,800),Beech=c(0,36,97,180,272,363,447,522,586,640))
vol_class$Spruce <- vol_class$Spruce+c(diff(vol_class$Spruce)/2,(diff(vol_class$Spruce)/2)[9])
vol_class$Beech <- vol_class$Beech+c(diff(vol_class$Beech)/2,(diff(vol_class$Beech)/2)[9])

resultstates$vol_mean <- 0

for (i in 1:10){
  resultstates[resultstates$species == levels(resultstates$species)[1] & as.numeric(resultstates$vol) == i, ]$vol_mean <- vol_class$Beech[i]
  resultstates[resultstates$species == levels(resultstates$species)[2] & as.numeric(resultstates$vol) == i, ]$vol_mean <- vol_class$Spruce[i]
}
df <- resultstates
df[,6:26] <- df[,6:26]*df[,27]
Bu_Ger <- colSums(df[resultstates$species == levels(df$species)[1],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[1],6:26])

Fi_Ger <- colSums(df[resultstates$species == levels(df$species)[2],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[2],6:26])
plot(seq(0,20), Bu_Ger, ylim=c(0,400),
     main = "Standing Volume Beech",
     xlab = "Simulation Step [n]",
     ylab = expression("Standing volume [m"^3~"ha"^{-1}~"o.b.]"))
plot(seq(0,20), Fi_Ger, ylim=c(0,600),
     main = "Standing Volume Spruce",
     xlab = "Simulation Step [n]",
     ylab = expression("Standing volume [m"^3~"ha"^{-1}~"o.b.]"))

df <- resultstates
df[,6:26] <- df[,6:26]*df[,27]


rawoutput$vol_mean <- 0

for (i in 1:10){
  rawoutput[rawoutput$species == levels(rawoutput$species)[1] & as.numeric(rawoutput$vol) == i, ]$vol_mean <- vol_class$Beech[i]
  rawoutput[rawoutput$species == levels(rawoutput$species)[2] & as.numeric(rawoutput$vol) == i, ]$vol_mean <- vol_class$Spruce[i]
}

df <- rawoutput
df <- df[,-seq(6,46,2)]
df[,6:26] <- df[,6:26]*df[,27]
Bu_ff_Ger <- colSums(df[df$species == levels(df$species)[1],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[1],6:26])
Bu_ff_Ger <- Bu_ff_Ger/5
Bu_ff_Ger

Fi_ff_Ger <- colSums(df[df$species == levels(df$species)[2],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[2],6:26])
Fi_ff_Ger <- Fi_ff_Ger/5
Fi_ff_Ger

plot(seq(0, 20), Bu_ff_Ger, ylim = c(0, round(max(Bu_ff_Ger), 0)), 
     main = "Drain Beech", 
     xlab = "Simulation Step [n]", 
     ylab = expression("Drain [m"^3~"ha"^{-1}~"a"^{-1}~"o.b.]"))
plot(seq(0,20), Fi_ff_Ger, ylim=c(0,round(max(Fi_ff_Ger),0)),
     main = "Drain Spruce", 
     xlab = "Simulation Step [n]", 
     ylab = expression("Drain [m"^3~"ha"^{-1}~"a"^{-1}~"o.b.]"))




#########################################
# Scenario NEW ##########################
#########################################
#describing the factors:
factors.tmp<-readLines("DEfactors.txt")
nfact<-length(factors.tmp)
factors.tmp<-strsplit(factors.tmp," ")
factnames<-character(nfact)
for(i in 1:nfact) {
  factnames[i]<-factors.tmp[[i]][1]
  factors.tmp[[i]]<-factors.tmp[[i]][-1]
}
names(factors.tmp)<-factnames
factlvls<-factors.tmp
factdims<-sapply(factlvls,length)


#describing the activities
activities.tmp<-readLines("DEactivities.txt")
nact<-length(activities.tmp)
activities<-strsplit(activities.tmp," ")
actnames<-character(nact)
actmethod<-character(nact)
actfiles<-character(nact)
rm(activities.tmp)

for(i in 1:nact) {
  actnames[i]<-activities[[i]][1]
  actmethod[i]<-activities[[i]][2]
  actfiles[i]<-activities[[i]][3]
  activities[[i]]<-activities[[i]][-(1:3)]
}
names(activities)<-actnames

#setting up the transition matrices
transmats<-list()
for(i in 1:length(activities))
{
  #the levels of these factors may change under this activity:
  basedimnames<-activities[[i]] 
  #these are fixed under this activity:
  otherdimnames<-setdiff(factnames,basedimnames)
  basedim<-prod(factdims[basedimnames])
  otherdim<-factdims[otherdimnames]
  
  transmat.forthisround <- readtransprobs(actfiles[i]) 
  
  #the rest is making sure dimensions have correct and sensible names
  lvlnames<-factlvls[[basedimnames[1]]]
  if(length(basedimnames>1))
    for(j in 2:length(basedimnames))
    {
      lvlnames.more<-factlvls[[basedimnames[j]]]
      lvlnames<-paste(rep(lvlnames,times=length(lvlnames.more)),
                      rep(lvlnames.more,each=length(lvlnames)),sep=".")
    }
  dimnames(transmat.forthisround)<-c(list(lvlnames,lvlnames),factlvls[otherdimnames])
  transmats<-c(transmats,list(transmat.forthisround))
  names(transmats)[i]<-actnames[i]
}

rm(basedim,basedimnames,otherdim,otherdimnames,transmat.forthisround,lvlnames,
   lvlnames.more)

#setting up the initial state and the probabilities of activities
statespace<-read.table("DEinitstate.txt",header=TRUE)
#actprobtable<-read.table("DEactprobs_nothin.txt",header=TRUE)
actprobtable<-read.table("actprobtable_new.txt",header=TRUE)
#this next merge is a bit dumb but the fixing of factor levels and ordering
#needs to be done for both, so I'm doing it at once
statespace<-merge(statespace,actprobtable)
for(fname in factnames)
  statespace[[fname]]<-factor(statespace[[fname]],
                              levels=factlvls[[fname]],ordered=TRUE)
statespace<-statespace[do.call(order,statespace[rev(factnames)]),]

initstate<-statespace[[setdiff(names(statespace),union(factnames,actnames))]]
dim(initstate)<-factdims
dimnames(initstate)<-factlvls


actproblist<-list()
for(i in 1:nact) {
  actproblist[[i]]<-array(statespace[[actnames[i]]],dim=factdims)
  dimnames(actproblist[[i]])<-factlvls
}

#-----------------------------
#finally, the actual simulation

#resultstates<-statespace[factnames]

state<-dividebyA(initstate)
resultstates<-list(step0=state)


#for(i in 1:nact)
#resultstates[[actnames[i]]]<-data.frame(step0=as.vector(tmp[[i]]))

nrofsteps<-as.numeric(20) #how many steps should be taken?
#a silly compatibility issue there with the names, sorry about that

for (i in 1:nrofsteps) {
  nstate<-newstate(state)  
  state<-rowSums(do.call(cbind,lapply(nstate,as.vector)))
  dim(state)<-factdims
  state<-dividebyA(state)
  resultstates[[paste0("step",i)]]<-state
}
#resultstates is now a list of lists, outer list going over the simulation steps, 
#and inner lists going over the activities. The inmost elements are multiway arrays
#in the form of the statespace (that is, factdims)

#Reformulate resultstates to be data frame whose columns are arrays, with nact columns.
#First, vectorize and cbind the multiway arrays:
resultstates<-lapply(resultstates,function(elem) {do.call(cbind,lapply(elem,as.vector))})
#Second, make the data frame:
resultstates<-data.frame(resultstates)

rawoutput<-statespace[factnames]
rawoutput<-data.frame(rawoutput,resultstates)
#rawoutput has 'nact' columns per step, with names suchs as step0.actname
write.table(rawoutput,file="rawoutput.txt")

tmp<-as.matrix(resultstates)
dim(tmp)<-c(dim(tmp)[1],nact,nrofsteps+1)
tmp<-apply(tmp,c(1,3),FUN=sum)
resultstates<-data.frame(tmp)
names(resultstates)<-paste0("step",0:nrofsteps)
resultstates<-data.frame(statespace[factnames],resultstates)
#resultstates is like rawoutput, but with the different columns with name
#such as step0.something summed together and named just step0
write.table(resultstates,file="resultstates.txt")



######
vol_class <- list(Spruce=c(0,60,153,270,391,503,601,682,748,800),Beech=c(0,36,97,180,272,363,447,522,586,640))
vol_class$Spruce <- vol_class$Spruce+c(diff(vol_class$Spruce)/2,(diff(vol_class$Spruce)/2)[9])
vol_class$Beech <- vol_class$Beech+c(diff(vol_class$Beech)/2,(diff(vol_class$Beech)/2)[9])

resultstates$vol_mean <- 0

for (i in 1:10){
  resultstates[resultstates$species == levels(resultstates$species)[1] & as.numeric(resultstates$vol) == i, ]$vol_mean <- vol_class$Beech[i]
  resultstates[resultstates$species == levels(resultstates$species)[2] & as.numeric(resultstates$vol) == i, ]$vol_mean <- vol_class$Spruce[i]
}
df <- resultstates
df[,6:26] <- df[,6:26]*df[,27]
Bu_Ger <- colSums(df[resultstates$species == levels(df$species)[1],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[1],6:26])

Fi_Ger <- colSums(df[resultstates$species == levels(df$species)[2],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[2],6:26])
plot(seq(0,20), Bu_Ger, ylim=c(0,400),
     main = "Standing Volume Beech",
     xlab = "Simulation Step [n]",
     ylab = expression("Standing volume [m"^3~"ha"^{-1}~"o.b.]"))
plot(seq(0,20), Fi_Ger, ylim=c(0,500),
     main = "Standing Volume Spruce",
     xlab = "Simulation Step [n]",
     ylab = expression("Standing volume [m"^3~"ha"^{-1}~"o.b.]"))

df <- resultstates
df[,6:26] <- df[,6:26]*df[,27]
# Bu_Nuts
Bu_Nuts1 <- colSums(df[df$species == levels(df$species)[1] &
                         df$nuts == levels(df$nuts)[1],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[1] &
                         resultstates$nuts == levels(df$nuts)[1],6:26])
Fi_Nuts1 <- colSums(df[df$species == levels(df$species)[2] &
                         df$nuts == levels(df$nuts)[1],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[2] &
                         df$nuts == levels(df$nuts)[1],6:26])

plot(seq(0,20),Bu_Nuts1,ylim=c(0,400))
plot(seq(0,20),Fi_Nuts1,ylim=c(0,600))


rawoutput$vol_mean <- 0

for (i in 1:10){
  rawoutput[rawoutput$species == levels(rawoutput$species)[1] & as.numeric(rawoutput$vol) == i, ]$vol_mean <- vol_class$Beech[i]
  rawoutput[rawoutput$species == levels(rawoutput$species)[2] & as.numeric(rawoutput$vol) == i, ]$vol_mean <- vol_class$Spruce[i]
}

df <- rawoutput
df <- df[,-seq(6,46,2)]
df[,6:26] <- df[,6:26]*df[,27]
Bu_ff_Ger <- colSums(df[df$species == levels(df$species)[1],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[1],6:26])
Bu_ff_Ger <- Bu_ff_Ger/5
Bu_ff_Ger

Fi_ff_Ger <- colSums(df[df$species == levels(df$species)[2],6:26])/
  colSums(resultstates[resultstates$species == levels(resultstates$species)[2],6:26])
Fi_ff_Ger <- Fi_ff_Ger/5
Fi_ff_Ger

plot(seq(0, 20), Bu_ff_Ger, ylim = c(0, round(max(Bu_ff_Ger), 0)), 
     main = "Drain Beech", 
     xlab = "Simulation Step [n]", 
     ylab = expression("Drain [m"^3~"ha"^{-1}~"a"^{-1}~"o.b.]"))
plot(seq(0,20), Fi_ff_Ger, ylim=c(0,round(max(Fi_ff_Ger),0)),
     main = "Drain Spruce", 
     xlab = "Simulation Step [n]", 
     ylab = expression("Drain [m"^3~"ha"^{-1}~"a"^{-1}~"o.b.]"))

#########
# Overwrite finalfell values based on species and age
#actprobtable_new <- actprobtable
#actprobtable_new <- actprobtable_new %>%
#  mutate(finalfell = case_when(
#    species == "beech" ~ bu[age],
#    species == "spruce" ~ fi[age],
#    TRUE ~ finalfell  # Keep existing values if species doesn't match
#  ))
#actprobtable_new$nomgmt <- 1-actprobtable_new$finalfell
# Define the output file path
#output_file <- "actprobtable_new.txt"

# Write the dataframe to a .txt file with space-separated values
#write.table(actprobtable_new, file = output_file, row.names = FALSE, quote = FALSE, sep = " ")



##################################################
### France #######################################
##################################################
vz <- "C:/Users/regelmann/Documents/GitHub/efdm/toyDataSet/France"
setwd(vz)
getwd()

#describing the factors:
factors_oak.tmp<-readLines("FRfactors_oak.txt")
factors_beech.tmp<-readLines("FRfactors_beech.txt")
factors_conifers.tmp<-readLines("FRfactors_conifers.txt")

nfact_oak<-length(factors_oak.tmp)
nfact_beech<-length(factors_beech.tmp)
nfact_conifers<-length(factors_conifers.tmp)

factors_oak.tmp<-strsplit(factors_oak.tmp," ")
factors_beech.tmp<-strsplit(factors_beech.tmp," ")
factors_conifers.tmp<-strsplit(factors_conifers.tmp," ")

factnames_oak<-character(nfact_oak)
factnames_beech<-character(nfact_beech)
factnames_conifers<-character(nfact_conifers)

for(i in 1:nfact_oak) {
  factnames_oak[i]<-factors_oak.tmp[[i]][1]
  factors_oak.tmp[[i]]<-factors_oak.tmp[[i]][-1]
}
for(i in 1:nfact_beech) {
  factnames_beech[i]<-factors_beech.tmp[[i]][1]
  factors_beech.tmp[[i]]<-factors_beech.tmp[[i]][-1]
}
for(i in 1:nfact_conifers) {
  factnames_conifers[i]<-factors_conifers.tmp[[i]][1]
  factors_conifers.tmp[[i]]<-factors_conifers.tmp[[i]][-1]
}

names(factors_oak.tmp)<-factnames_oak
names(factors_beech.tmp)<-factnames_beech
names(factors_conifers.tmp)<-factnames_conifers

factlvls_oak<-factors_oak.tmp
factlvls_beech<-factors_beech.tmp
factlvls_conifers<-factors_conifers.tmp

factdims_oak<-sapply(factlvls_oak,length)
factdims_beech<-sapply(factlvls_beech,length)
factdims_conifers<-sapply(factlvls_conifers,length)


#describing the activities
activities_oak.tmp<-readLines("FRactivities_oak.txt")
activities_beech.tmp<-readLines("FRactivities_beech.txt")
activities_conifers.tmp<-readLines("FRactivities_conifers.txt")

nact_oak<-length(activities_oak.tmp)
nact_beech<-length(activities_beech.tmp)
nact_conifers<-length(activities_conifers.tmp)

activities_oak<-strsplit(activities_oak.tmp," ")
activities_beech<-strsplit(activities_beech.tmp," ")
activities_conifers<-strsplit(activities_conifers.tmp," ")

actnames_oak<-character(nact_oak)
actnames_beech<-character(nact_beech)
actnames_conifers<-character(nact_conifers)

actmethod_oak<-character(nact_oak)
actmethod_beech<-character(nact_beech)
actmethod_conifers<-character(nact_conifers)

actfiles_oak<-character(nact_oak)
actfiles_beech<-character(nact_beech)
actfiles_conifers<-character(nact_conifers)

rm(activities_oak.tmp,
   activities_beech.tmp,
   activities_conifers.tmp)

for(i in 1:nact_oak) {
  actnames_oak[i]<-activities_oak[[i]][1]
  actmethod_oak[i]<-activities_oak[[i]][2]
  actfiles_oak[i]<-activities_oak[[i]][3]
  activities_oak[[i]]<-activities_oak[[i]][-(1:3)]
}
for(i in 1:nact_beech) {
  actnames_beech[i]<-activities_beech[[i]][1]
  actmethod_beech[i]<-activities_beech[[i]][2]
  actfiles_beech[i]<-activities_beech[[i]][3]
  activities_beech[[i]]<-activities_beech[[i]][-(1:3)]
}
for(i in 1:nact_conifers) {
  actnames_conifers[i]<-activities_conifers[[i]][1]
  actmethod_conifers[i]<-activities_conifers[[i]][2]
  actfiles_conifers[i]<-activities_conifers[[i]][3]
  activities_conifers[[i]]<-activities_conifers[[i]][-(1:3)]
}
names(activities_oak)<-actnames_oak
names(activities_beech)<-actnames_beech
names(activities_conifers)<-actnames_conifers

#setting up the transition matrices
###### hier weiter machen!!!
transmats_oak<-list()
for(i in 1:length(activities_oak))
{
  #the levels of these factors may change under this activity:
  basedimnames<-activities[[i]] 
  #these are fixed under this activity:
  otherdimnames<-setdiff(factnames,basedimnames)
  basedim<-prod(factdims[basedimnames])
  otherdim<-factdims[otherdimnames]
  
  transmat.forthisround <- readtransprobs(actfiles_oak[i]) 
  
  #the rest is making sure dimensions have correct and sensible names
  lvlnames_oak<-factlvls_oak[[basedimnames_oak[1]]]
  if(length(basedimnames_oak>1))
    for(j in 2:length(basedimnames_oak))
    {
      lvlnames.more_oak<-factlvls_oak[[basedimnames_oak[j]]]
      lvlnames_oak<-paste(rep(lvlnames_oak,times=length(lvlnames.more_oak)),
                      rep(lvlnames.more_oak,each=length(lvlnames_oak)),sep=".")
    }
  dimnames(transmat.forthisround_oak)<-c(list(lvlnames_oak,lvlnames_oak),factlvls_oak[otherdimnames_oak])
  transmats_oak<-c(transmats_oak,list(transmat.forthisround_oak))
  names(transmats_oak)[i]<-actnames_oak[i]
}

rm(basedim,basedimnames,otherdim,otherdimnames,transmat.forthisround,lvlnames,
   lvlnames.more)

#setting up the initial state and the probabilities of activities
statespace<-read.table("FRinitstate_oak.txt",header=TRUE)
actprobtable<-read.table("FRactprobs_oak.txt",header=TRUE) #neues Einfügen
#actprobtable<-read.table("actprobtable_new.txt",header=TRUE) #neues Einfügen
#this next merge is a bit dumb but the fixing of factor levels and ordering
#needs to be done for both, so I'm doing it at once
statespace<-merge(statespace,actprobtable)
for(fname in factnames)
  statespace[[fname]]<-factor(statespace[[fname]],
                              levels=factlvls[[fname]],ordered=TRUE)
statespace<-statespace[do.call(order,statespace[rev(factnames)]),]

initstate<-statespace[[setdiff(names(statespace),union(factnames,actnames))]]
dim(initstate)<-factdims
dimnames(initstate)<-factlvls


actproblist<-list()
for(i in 1:nact) {
  actproblist[[i]]<-array(statespace[[actnames[i]]],dim=factdims)
  dimnames(actproblist[[i]])<-factlvls
}

#-----------------------------
#finally, the actual simulation

#resultstates<-statespace[factnames]

state<-dividebyA(initstate)
resultstates<-list(step0=state)


#for(i in 1:nact)
#resultstates[[actnames[i]]]<-data.frame(step0=as.vector(tmp[[i]]))

nrofsteps<-as.numeric(20) #how many steps should be taken?
#a silly compatibility issue there with the names, sorry about that

for (i in 1:nrofsteps) {
  nstate<-newstate(state)  
  state<-rowSums(do.call(cbind,lapply(nstate,as.vector)))
  dim(state)<-factdims
  state<-dividebyA(state)
  resultstates[[paste0("step",i)]]<-state
}
#resultstates is now a list of lists, outer list going over the simulation steps, 
#and inner lists going over the activities. The inmost elements are multiway arrays
#in the form of the statespace (that is, factdims)

#Reformulate resultstates to be data frame whose columns are arrays, with nact columns.
#First, vectorize and cbind the multiway arrays:
resultstates<-lapply(resultstates,function(elem) {do.call(cbind,lapply(elem,as.vector))})
#Second, make the data frame:
resultstates<-data.frame(resultstates)

rawoutput<-statespace[factnames]
rawoutput<-data.frame(rawoutput,resultstates)
#rawoutput has 'nact' columns per step, with names suchs as step0.actname
write.table(rawoutput,file="rawoutput.txt")

tmp<-as.matrix(resultstates)
dim(tmp)<-c(dim(tmp)[1],nact,nrofsteps+1)
tmp<-apply(tmp,c(1,3),FUN=sum)
resultstates<-data.frame(tmp)
names(resultstates)<-paste0("step",0:nrofsteps)
resultstates<-data.frame(statespace[factnames],resultstates)
#resultstates is like rawoutput, but with the different columns with name
#such as step0.something summed together and named just step0
write.table(resultstates,file="resultstates.txt")


vol_class <- list(Oak=(seq(1, 10)*40)-20,
                  Beech=(seq(1, 10)*40)-20,
                  Conifer=(seq(1, 10)*50)-25)
age_class <- list(Oak=(seq(1, 37)*5)-2.5,
                  Beech=(seq(1, 31)*5)-2.5,
                  Conifer=(seq(1, 17)*5)-2.5)

resultstates$vol_mean <- 0

for (i in 1:10){
  resultstates[as.numeric(resultstates$vol) == i, ]$vol_mean <- vol_class$Beech[i]
}
df <- resultstates
df[,4:24] <- df[,4:24]*df[,25]
Oak_Fr <- colSums(df[,4:24])/
  colSums(resultstates[,4:24])

plot(seq(0,20), Oak_Fr, ylim=c(0,300),
     main = "Standing Volume Oak",
     xlab = "Simulation Step [n]",
     ylab = expression("Standing volume [m"^3~"ha"^{-1}~"o.b.]"))

## Drain
rawoutput$vol_mean <- 0

for (i in 1:10){
  rawoutput[as.numeric(rawoutput$vol) == i, ]$vol_mean <- vol_class$Beech[i]
}

df <- rawoutput
df <- df[,-seq(4,64,3)]
df[,4:45] <- df[,4:45]*df[,46]
Oak_ff_Fr <- (colSums(df[seq(4,45,2)])+colSums(df[seq(5,46,2)]))/
  colSums(resultstates[,4:24])
Oak_ff_Fr <- Oak_ff_Fr/5
names(Oak_ff_Fr) <- seq(0,20)
Oak_ff_Fr

plot(seq(0, 20), Oak_ff_Fr, ylim = c(0, round(max(Oak_ff_Fr), 0)), 
     main = "Drain Oak", 
     xlab = "Simulation Step [n]", 
     ylab = expression("Drain [m"^3~"ha"^{-1}~"a"^{-1}~"o.b.]"))
