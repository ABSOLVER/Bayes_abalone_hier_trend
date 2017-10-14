
##R script to fit Bayesian non-linear trend model to Fisheries Victoria Abalone data
##Code by Jim Thomson: jim.thomsson@gmail.com
##Requires that WinBUGS  is installed: available at: http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/
##Also requires Dave Lunn's reversible jump MCMC add on for WinBUGS, available at: http://www.winbugs-development.org.uk
##and all of the R packages listed below.
##Refer to "Abalone Trend Model: User Notes" for instructions

.libPaths("C:/Users/Jim Thomson/Documents/Rlibrary") # Ignore this..

#load packages
library(R2WinBUGS)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(reshape2)
library(scales)

##SET DIRECTORIES
#setwd("~/My Dropbox/Private/Abalone/AbaloneModel") #set working directory (allr req'd files must be in this directory)
wd=getwd()
saveto=paste0(wd,"/Outputs/") # set folder where model outputs should be saved
dir.create(saveto) #create save-to folder if doesn't exist
bugs.directory="C:/WinBUGS14"# Specify where WinBUGS is installed, usually in "C:/Program Files/WinBUGS14"  

##LOAD DATA and SET RESPONSE VARIABLE
alldat=read.csv("RawAbundance.csv") # load Abalone data
##load lookup table for SMUs and reefcodes...the script loops thru all SMUs in the lookup table...see LINE 57
#reefcodes=read.csv("AbaloneReefCodes_AVG.csv") # load reefcode / SMU lookup table..this one for AVG reefs only
reefcodes=read.csv("AbaloneReefCodes.csv") # load reefcode / SMU lookup table..this one for all reefs

Response="TotalCount"  #Name of response variable (Juveniles, Prerecruits,Recruits, TotalCount)
## MODEL OPTIONS
family="poisson" #negbin" #poisson" #error distribution for model..poisson" = poisson (with extra poisson via case-level random error),"negbin" = negative binomial
hierarch.trends=T # T to model trends hierarchically wihin MUs: assumes a common mean trend for each ZONE with MU-specific trends "shrunk" toward mean trend 
                  # F to treat MU-specific trends as completely independent 
armod=T # include (first order) autoregressive errors 
cont=1 # 0 to allow step changes, 1 for continuous trends
base.start=1 #baseline period for each SMU starts at year 1 (used in calculated Pr(>baseline) etc)
base.end=4  #baseline period for each SMU ends at year 4
nits=200 # number of MCMC iterations...start small to check working & speed, then aim for min. of 10000
debug=F # T to manually check MCMC chains etc in WinBUGS after each fit
auto.kmax=T #set max.number of changepoints automatically (to floor(Nyear/5))..if F, set kmax below
kmax=NULL # max. number of change changepoints, ignored if auto.kmax=T
#covariate settings
covar=F # T to include covariates in model: for now, THIS ONLY WORKS WITH AbaloneReefCodes_AVG.csv ..see line 23
prior.prob.cov=0.5 # 0.5 for Bayesian model selection / averaging with standard priors (all models equally probablle a priori)
auto.Kmax.cov=T # should max. no. of covariates in model be automatically set? If FALSE use Kmax.cov to set.
Kmax.cov=NULL #max. no. of covariates in model, by default is set to number of candidates (line), set here to change
addyear=0 #for AVG covariates - 0 to make putative change-point in year of first AVG detection, 1 to make it year after first detection (or any other integer) 

## _save.name is appended to SMU name when saving fitted models and result plots and tables.
#is generated automatically here based on above options, or specifiy your own
save.name=paste(ifelse(hierarch.trends,"hier","ind"),ifelse(armod,"_ar",""),ifelse(cont==0,"_step",""),ifelse(covar,"_covar",""),ifelse(family=="negbin","_negbin",""),sep="") 

#set up a matrix of periods to calc. mean trends for: each row is a period, 1st column is first year, 2nd column last year
#by default  (i.e. if auto.trend.period=T), calc's trends for last 2,3,4, and 5 years (1 year trends are automatically calculated for all years)
auto.trend.period=T #calc's trends for last 2,3,4, and 5 years
trend.period=NULL  # a matrix defining periods to calculate mean trends over (generated automatically if auto.trend.predio=T)


source("AbaloneModel_makeBUGSfile_hier.R")
source("AbaloneModel_makeBUGSfile_indep.R")

alldat$Date=dmy(alldat$Date)

setwd(saveto)
for(lev in levels(reefcodes$ZONE)){   ###ZONE: CHANGED LINE 

reefcode=reefcodes[reefcodes$ZONE==lev,1]  ###ZONE: CHANGED LINE

dat=alldat[alldat$ReefCode%in%reefcode,]

not30=which(dat$SwimLength!=30)
if(length(not30)>0){dat=dat[-not30,]}
#swimlen=ifelse(dat$SwimLength==0,30,dat$SwimLength)
#dat$AdjCount=dat$TotalCount*(30/swimlen) 

melt=melt(dat,id.vars=c(colnames(dat)[1:10])) 

dat=dcast(melt,Date+QuotaYear+Site+ReefCode+Diver~variable,mean)

#add SMU code to data
dat$SMU=reefcodes[match(dat$ReefCode,reefcodes$ReefCode),"Spatial.Management.Units"] #ZONES:NEW LINE
#now change the names, so that SMU's are treated as "ReefCode"s hereafter, and actual ReefCode become "reefblocks" 
names(dat)=gsub("ReefCode","reefblock",names(dat)) #ZONES:NEW LINE
names(dat)=gsub("SMU","ReefCode",names(dat)) #ZONES:NEW LINE

y=dat[,Response]
year=dat$QuotaYear

yearf=as.numeric(factor(year))
site=as.numeric(factor(dat$Site))
reef=as.numeric(factor(dat$ReefCode))
reefyear=as.numeric(ordered(factor(reef):factor(yearf)))
diver=as.numeric(factor(dat$Diver))
reefblock=as.numeric(factor(dat$reefblock)) #ZONES:NEW LINE
Nreef=max(reef)
components=cbind(reef,reefblock,yearf,diver,site) #ZONES:CHANGED LINE...add "reefblock" (i.e. ReefCode) as a random intercept (like site)

Nlevels=apply(components,2,max)
#omitlevs=which(Nlevels==1)

#if(length(omitlevs)>0){components=components[,-omitlevs];Nlevels=apply(components,2,max)}

if(max(reef)>1){components=cbind(components,reefyear);Nlevels=apply(components,2,max)}

#Nmain=ncol(components)
#Nint=nrow(int)
Nbatch=length(Nlevels)

N=length(y)

Nyear=max(yearf)

if(auto.kmax){kmax=max(floor(Nyear/5),1)}

meanyr=min(dat$QuotaYear); sdyr=sd(dat$QuotaYear);yr1=min(dat$QuotaYear)-1

prevcount=NULL
prev.st=prev.end=NULL
noprev=0
if(armod){
  for(i in 1:N){
    prev=which(yearf==yearf[i]-1 & site==site[i])
    prevm=mean(y[prev])
    prevreef=which(yearf==yearf[i]-1 & reef==reef[i])
    prevreefm=mean(y[prevreef])
    ps=min(prev)
    pe=max(prev)
    if(is.nan(prevm)){prevm=prevreefm;ps=min(prevreef);pe=max(prevreef)}
    if(is.nan(prevm)){noprev=noprev+1;ps=N+noprev;pe=ps}
    prevcount=c(prevcount,prev)
    prev.st=c(prev.st,ps);prev.end=c(prev.end,pe)
  }
maxprev=max(prev.end)
}

prevcount=(prevcount-mean(prevcount,na.rm=T))/sd(prevcount,na.rm=T)

yearx=c(1:Nyear)
base.period=c(min(year):max(year))[c(base.start,base.end)]


reefyr.ind=matrix(max(components[,Nbatch])+1,nrow=Nreef+1,ncol=Nyear)
for(i in 1:Nreef){
  for(j in 1:Nyear){
    reefyr.ind[i,j]=components[which(reef==i & yearf==j),Nbatch][1]
    if(is.na(reefyr.ind[i,j])){reefyr.ind[i,j]=max(components[,Nbatch])+1}
  }
}

if(auto.trend.period){trend.period=cbind(Nyear-c(5:2),rep(Nyear,4))}
Ntps=nrow(trend.period)

if(covar){
  if(Nreef==1){components=cbind(components,components[,"yearf"]);colnames(components)[ncol(components)]="reefyear"}
  reefyear=components[,"reefyear"]
  detectyr=reefcodes[match(reefcode,reefcodes$ReefCode),"detectyr"]
  detectyr.no=detectyr-min(year)+1+addyear
  reefyearmat=unique(components[,c("reef","yearf","reefyear")])
  reefyearmat=reefyearmat[order(reefyearmat[,3]),]
  
  Xstep=ifelse(reefyearmat[,"yearf"]>detectyr.no[reefyearmat[,"reef"]],1,0)
  Xtrend=Xstep*reefyearmat[,"yearf"]-detectyr.no[reefyearmat[,"reef"]]*Xstep
  Xcov=NULL
  for(c in 1:Nreef){Xcov=cbind(Xcov,Xstep*ifelse(reefyearmat[,"reef"]==c,1,0))}
  for(c in 1:Nreef){Xcov=cbind(Xcov,Xtrend*ifelse(reefyearmat[,"reef"]==c,1,0))}
  #Xcov=cbind(Xcov,Xstep,Xtrend)
  Qcov=ncol(Xcov)
  if(auto.Kmax.cov){Kmax.cov=Qcov}
  Nyrf=nrow(Xcov)
  Xcov=rbind(Xcov,0)
  colnames(Xcov)=c(paste("Reef", c(1:Nreef),"step",sep=""),paste("Reef", c(1:Nreef),"trend",sep=""))#,"All.step","All.trend")
  if(Nreef==1){components=components[,-ncol(components)]}
}

bugsdata=list("y","yearf","yearx","N","Nreef","reef","components","Nlevels","Nbatch","Nyear","kmax","base.start","base.end","reefyr.ind","trend.period","Ntps")
if(armod){bugsdata=c(bugsdata,"prev.st","prev.end","maxprev")}
beta.ini=matrix(NA,ncol=Nbatch,nrow=max(Nlevels)+ifelse(Nreef>1,1,0))
for(b in 1:Nbatch[1]){beta.ini[1:Nlevels[b],b]=0}
inits=function(){list(sd.rand=rep(0.5,Nbatch),sd=rep(0.1,ifelse(hierarch.trends,3,Nreef+2)),alpha=rnorm(1),beta=beta.ini,mu1=rnorm(N),k=rep(1,Nreef+1))}
if(family=="negbin"){
  inits=function(){list(nb.rho=rep(1,N),logalpha=1,sd.rand=rep(0.5,Nbatch),sd=rep(0.1,ifelse(hierarch.trends,3,Nreef+2)),alpha=rnorm(1),beta=beta.ini,mu1=rnorm(N),k=rep(1,Nreef+1))}
}
params=c("logmean","sd","beta","mu","trend","pch","pinc","fitted","fitted.yr","p.above","p.ok","mean.trend","pmt","ppp")
if(covar){bugsdata=c(bugsdata,"Qcov","Xcov","Nyrf","reefyear","Kmax.cov","prior.prob.cov");params=c(params,"cov.beta","inc.cov")}
filename=paste("Ab_model_",lev,".txt",sep="")
make.model.file.ind(filename,Nreef,AR=armod,covar=covar,family=family)
if(hierarch.trends & Nreef>1){make.model.file.hier(filename,Nreef,AR=armod,covar=covar,family=family)}  

fit=bugs(bugsdata,inits,params,filename,n.iter=nits,debug=debug,bugs.directory=bugs.directory)
save(fit,file=paste(lev,save.name,".R",sep=""))

dat$yadj=y
dat$qyear=dmy(paste("01/01/",dat$QuotaYear,sep=""))
yr=sort(unique(dat$QuotaYear))
yr.date=dmy(paste("01/01/",yr,sep=""))

fitted=data.frame(fit$summary[grep("fitted\\[",rownames(fit$summary)),],yr,yr.date)
fitted.yr=data.frame(fit$summary[grep("fitted.yr\\[",rownames(fit$summary)),],yr,yr.date);fitted.yr=fitted.yr[-grep("log",rownames(fitted)),]
trends=fit$summary[grep("trend",rownames(fit$summary)),];trends[,c(1,3:7)]=trends[,c(1,3:7)]-1
mean.trends=data.frame(trends[grep("mean.trend",rownames(trends)),])
trends=trends[-grep("mean.trend",rownames(trends)),]
trends=data.frame(trends,round(fit$summary[grep("pch",rownames(fit$summary)),1],3),round(fit$summary[grep("pinc",rownames(fit$summary)),1],3),round(fit$summary[grep("p.above",rownames(fit$summary)),1],3),round(fit$summary[grep("p.ok",rownames(fit$summary)),1],3),yr,yr.date)
names(fitted)[3:7]=names(trends)[3:7]=names(mean.trends)[3:7]=paste("pc",gsub("%","",colnames(fit$summary)[3:7]),sep="")
names(trends)[10:13]=c("p.change","p.increase","p.above","p.ok")
mean.trends$pinc=as.vector(t(round(fit$mean$pmt,2)))

prior.ch=kmax*0.5/Nyear
trends$ORs=trends$p.change*(1-prior.ch)/((1-trends$p.change)*prior.ch)

reefnames=c(levels(factor(dat$ReefCode)),"All reefs")

ppps=fit$mean$ppp
names(ppps)=reefnames
write.csv(ppps,file=paste0("PPP_",lev,save.name,".csv"))
  
reefnos=as.numeric(gsub(",.*","",gsub("mean.trend\\[","",rownames(mean.trends))))
yearnos=as.numeric(gsub("\\]","",gsub(".*,","",gsub("mean.trend\\[","",rownames(mean.trends)))))
start.year=trend.period[yearnos,1]+max(year)-Nyear
end.year=trend.period[yearnos,2]+max(year)-Nyear
mean.trends=data.frame(reefnames[reefnos],start.year,end.year,mean.trends)
names(mean.trends)[1]=c("ReefCode")
write.csv(mean.trends,file=paste("mean trends_",lev,save.name,".csv",sep=""))

reefnos=as.numeric(gsub(",.*","",gsub("trend\\[","",rownames(trends))))
yearnos=as.numeric(gsub("\\]","",gsub(".*,","",gsub("trend\\[","",rownames(trends)))))
yearvalue=yearnos +max(year)-Nyear
trends=data.frame(reefnames[reefnos],yearvalue,trends)
names(trends)[1:2]=c("ReefCode","QuotaYear")
write.csv(trends,file=paste("annual trends_",lev,save.name,".csv",sep=""))

if(covar){
  cov.res=cbind(fit$summary[grep("cov.beta",rownames(fit$summary)),],round(fit$mean$inc.cov,2))
  colnames(cov.res)[3:7]=paste("pc",gsub("%","",colnames(cov.res)[3:7]),sep="")
  colnames(cov.res)[ncol(cov.res)]="post.prob"
  rownames(cov.res)=colnames(Xcov)
  cov.prior=Kmax.cov*prior.prob.cov/Qcov
  OR.cov=cov.res$post.prob*(1-cov.prior)/((1-cov.res$post.prob)*cov.prior)
  cov.res=data.frame(cov.res,OR.cov)
  write.csv(cov.res,file=paste("Covariate_",lev,save.name,".csv",sep=""))
}

#plot results
pdf(paste(lev,save.name,".pdf",sep=""))
for(reefnum in 1:(Nreef+1)){

  fitted.r=fitted[grep(paste("\\[",reefnum,",",sep=""),rownames(fitted)),]
  trends.r=trends[grep(paste("\\[",reefnum,",",sep=""),rownames(trends)),]
  mean.trends.r=mean.trends[grep(paste("\\[",reefnum,",",sep=""),rownames(mean.trends)),]
  dat.r=dat[reef==reefnum,]
  if(length(unique(dat.r$ReefCode))>1){stop("check data, >1 reef code in subset")}
  if(reefnum==(Nreef+1)){dat.r=dat}
  p=ggplot()+geom_point(data=dat.r,aes(x=Date,y=yadj))
  CIs=ggplot()+geom_ribbon(data=fitted.r,aes(x=yr.date,ymin=pc2.5,ymax=pc97.5),alpha=0.5,fill="grey80",color=NA) +geom_ribbon(data=fitted.r,aes(x=yr.date,ymin=pc25,ymax=pc75),fill="grey",color="black") 
  count.plot=CIs + geom_point(data=dat.r,aes(x=Date,y=yadj),color="grey30") + geom_line(data=fitted.r,aes(x=yr.date,y=mean)) + ylab("Abalone Count") + xlab("Quota Year")+ ggtitle(paste(lev,": ",reefnames[reefnum],sep="")) + theme_bw()
  p2=ggplot(data=trends.r) + geom_ribbon(aes(x=yr.date,ymin=pc2.5,ymax=pc97.5),alpha=0.5,fill="grey80",color=NA) + geom_ribbon(aes(x=yr.date,ymin=pc25,ymax=pc75),fill="grey",color="black") + geom_line(aes(x=yr.date,y=mean),lwd=2) 
  trend.plot=p2+xlab("Quota Year") + ylab("Trend (proportional change per year)") + geom_abline(intercept=0,slope=0)  + theme_bw()
  
  texty=max(max(dat.r$yadj),max(fitted.r$pc97.5),na.rm=T)*.9
  count.plot=count.plot+ annotate("text",x=yr.date[nrow(trends.r)],y=texty,label=paste("\n","Pr(not declining)=",trends.r[nrow(trends.r),"p.increase"],"\n",
  "Pr(above baseline mean)=",trends.r[nrow(trends.r),"p.above"],"\n",
  "Pr(above baseline & not declining)=",trends.r[nrow(trends.r),"p.ok"]),size=4,hjust=1)#print(count.plot)
  grid.arrange(count.plot,trend.plot,heights=c(1.5,1))
 
}
dev.off()
}
setwd(wd)
