fwhm<-function(x,main) {
   rge<-range(x)
   xa<-floor(rge[1]/10)*10
   xb<-ceiling(rge[2]/10)*10
   nb<-(xb-xa)/10+1     
   d <- density(na.omit(x),n=256,bw="sj")
   xmax <- d$x[d$y==max(d$y)]        
   x1 <- d$x[d$x < xmax][which.min(abs(d$y[d$x < xmax]-max(d$y)/2))]
   x2 <- d$x[d$x > xmax][which.min(abs(d$y[d$x > xmax]-max(d$y)/2))]
   h<-hist(x,breaks=1000,xlim=c(xmax-2*(x2-x1),xmax+2*(x2-x1)),main=main)
   abline(v=c(x1,x2),col="red")   
   return(list(x2-x1,xmax))
}


load("all_data_H.RData")
names(H)
# [1] "HLADR"     "CXCR4"     "CD28"      "CXCR3"     "CD8"       "CXCR5"    
# [7] "CD127"     "CD69"      "TIM3"      "CD3"       "CD45RO"    "PD1"      
#[13] "Tetramer"  "CD4"       "CTLA4"     "CD25"      "ICOS"      "LAG3"     
#[19] "CD62L"     "tetraPos"  "tetraID"   "subjectID" "cohortID" 


features<-c("HLADR","CXCR4","CD28","CXCR3","CD8","CXCR5","CD127","CD69","TIM3","CD3","CD45RO","PD1","CD4","CTLA4","CD25","ICOS","LAG3","CD62L")
tetras<-c("12-20","13-21","CP11","CP18")
#these subjects are ordered by cohort
subjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD009","RAD010","RAD011",
            "TNET006","TNET007","TNET008","TNET009","TNET010","TNET011","TNET012","TNET013","TNET014","TNET015")
cohorts<-c("NBD","ES","JD","MON")

longcohorts<-rep(cohorts,c(5,5,4,10))
names(longcohorts)<-subjects

transf<-function(x,a) asinh(x/a)

ns<-length(subjects)
#but do this for CD4+ cells ONLY
#define cutoffs that define CD4+CD8- population. These have to be set by hand after inspection of scatterplots because the flow cytometer is set differently for every run.
xmin<-rep(10000,ns)
ymax<-rep(5000,ns)
xmin[c(16)]<-30000
xmin[c(24)]<-5000
xmin[c(10,11,13,14,15,18)]<-20000

ymax[c(7,11)]<-6000
ymax[c(13)]<-7000
names(xmin)<-names(ymax)<-subjects

subjectID<-H$subjectID
good<-rep(FALSE,dim(H)[1]) #define mask that will become TRUE for CD4+CD8- cells
pdf("CD4_CD8_scatterplots.pdf")
for (subject in subjects) {
   plot(H[subjectID==subject,c("CD4","CD8")],pch=19,cex=.1,col="#000000EE",xlim=c(-1e4,2e5),ylim=c(-1e4,1e4),main=subject)
   abline(v=xmin[subject],col="grey")
   abline(h=ymax[subject],col="grey")
   good<-good | subjectID==subject & H[,"CD4"]>xmin[subject] & H[,"CD8"]<ymax[subject]
}
dev.off()
sum(good)
#[1] 945892 #this many CD4+CD8- cells total


CD4set<-H[good,c(-5,-14)] #drop CD8, CD4 channels
save(CD4set,file="CD4set.RData")


#set zero peaks at 0 if they are <0
#features with a zero peak:
pdf("control.pdf")
single.features<-c("HLADR","CXCR4","CD28","CXCR3","CXCR5","CD127","CD69","TIM3","PD1","CTLA4","CD25","ICOS","LAG3")
X<-data.matrix(CD4set[,single.features])
for (s in subjects) {
for (f in single.features) {
   selection<-CD4set$subjectID==s
   u<-X[selection,f]
   rge<-quantile(u,probs=c(.01,.99))
   u<-u[u>rge[1] & u<rge[2]] #don't deal with "outliers"
   fw<-fwhm(u,paste(s,f,sep=" "))
   zero<-fw[[2]]
   if (zero<0) X[selection,f]<-X[selection,f]-zero #set mode at 0
   cat(s,f,zero,"\n")
}
}
dev.off()
CD4set[,single.features]<-X #done, zero peaks are at >=0


#now normalize PE using fwhm
Y<-CD4set$Tetramer
nt<-length(cohorts)
pdf("Tetramer_distributions.pdf")
for (j in 1:4) {
   tetra<-tetras[j] 
   for (it in 1:ns) {
      t<-subjects[it]
      cat(tetra,t,"\n")
      selection<-CD4set$subjectID==t & CD4set$tetraID==tetra
      if (sum(selection)==0) next
      u<-Y[selection]
      v<-Y[selection & CD4set$tetraPos] 
      h<-hist(u,breaks=1000,plot=FALSE)
      plot(h$mids,h$counts,log="y",xlim=c(-10000,100000),main=paste(t,tetra,sep=" "))
      abline(v=0,col="grey")
      rug(v,col="red")
   }
}
dev.off()


pdf("fwhm.pdf")
for (j in 1:4) {
   tetra<-tetras[j] 
   for (it in 1:ns) {
      t<-subjects[it]
      selection<-CD4set$subjectID==t & CD4set$tetraID==tetra
      cat(tetra,t,sum(selection),"\n")
      if (sum(selection)==0) next
      u<-Y[selection] 
      rge<-quantile(u,probs=c(.001,.99))
      u<-u[u>rge[1] & u<rge[2]] #don't deal with "outliers"
      fw<-fwhm(u,paste(t,tetra,sep=" "))
      normfactor<-fw[[1]] #this is fwhm
      zeromode<-fw[[2]]
      Y[selection]<-(Y[selection]-zeromode)/normfactor
   }
}
dev.off()
save(Y,file="normPE.RData")






#set the 99%-ile of every PE distribution to that of the normal distribution by scaling. Then, take pos to be PE>5 (5 sd away from 0)!
qnorm99<-qnorm(.99) #2.326348
pdf("tetra_threshold.pdf",width=20,height=8)
par(mfrow=c(1,4),pty="s")
   for (t in subjects) {
      for (tetra in tetras) {
      selection<-CD4set$subjectID==t & CD4set$tetraID==tetra
      if (sum(selection)==0) {
         plot.new()
         next
      }
      x<-Y[selection]
      q99<-quantile(x,.99)
      x<-x/q99*qnorm99 #so that value of 6 would be 5 sd
      Y[selection]<-x
      y<-x[x>0] #only positive
      breaks=seq(0,max(y)+1,by=0.1)
      h<-hist(y,breaks=breaks,plot=FALSE)
      plot(h$mids,h$counts,log="y",main=paste(t,tetra,sep=" "),xlab="transformed intensity",ylab="frequency",xlim=c(0,30))
      pos<-CD4set$tetraPos[selection]
      rug(x[pos],col="red")
      abline(v=c(qnorm99,5),col="grey")
      }
   }
dev.off()
save(Y,file="qnormPE.RData")
CD4set$qTetramer<-Y
CD4set$qtetraPos<-Y>5.
save(CD4set,file="qCD4set.RData")


opos<-CD4set$tetraPos #OG pos
qpos<-CD4set$qtetraPos

pdf("qtetra_threshold.pdf",width=20,height=8)
par(mfrow=c(1,4),pty="s")
   for (t in subjects) {
      for (tetra in tetras) {
      selection<-CD4set$subjectID==t & CD4set$tetraID==tetra
      if (sum(selection)==0) {
         plot.new()
         next
      }
      x<-Y[selection]
      q99<-quantile(x,.99)
      x<-x/q99*qnorm99
      Y[selection]<-x
      y<-x[x>0] #only positive
      breaks=seq(0,max(y)+1,by=0.1)
      h<-hist(y,breaks=breaks,plot=FALSE)
      plot(h$mids,h$counts,log="y",main=paste(t,tetra,sep=" "),xlab="transformed intensity",ylab="frequency",xlim=c(0,50))

      pos<-CD4set$qtetraPos[selection]
      rug(x[pos],col="red")
      abline(v=5,col="grey")
      }
   }
dev.off()



#the following is originally from Step1_separate_into_lineages
load("qCD4set.RData")

features<-names(CD4set)[1:17] #includes PE
nf<-length(features)


for (t in subjects) {
   cat(t, sum(CD4set$subjectID==t & CD4set$tetraID==tetras[1]),sum(CD4set$subjectID==t & CD4set$tetraID==tetras[2]),sum(CD4set$subjectID==t & CD4set$tetraID==tetras[3]),sum(CD4set$subjectID==t & CD4set$tetraID==tetras[4]),"\n")
}
#NBD-BRI001 2475 3671 2891 4204 
#NBD-BRI002 1375 1314 1209 1338 
#NBD-BRI003 2322 2719 2018 2182 
#NBD-BRI004 2766 3488 3199 2970 
#NBD-BRI005 5408 4672 5506 5945 
#T030 31503 30599 31703 29197 
#T031 7861 6979 8032 8419 
#T033 10617 11551 11309 10952 
#T034 12943 14106 11482 0 
#T035 27148 27477 27270 27176 
#RAD007 10409 13253 10633 11905 
#RAD009 6246 6271 6947 6744 
#RAD010 17392 17636 18905 20309 
#RAD011 5172 5212 5167 5675 
#TNET006 13014 12111 15151 14095 
#TNET007 8586 8115 6826 7807 
#TNET008 10413 14401 16308 13437 
#TNET009 26699 21192 13182 21783 
#TNET010 1316 1652 1511 1666 
#TNET011 4100 5204 3467 4689 
#TNET012 7531 5911 8045 8039 
#TNET013 3326 3345 0 7185 
#TNET014 7337 7289 8839 7151 
#TNET015 5988 5063 5217 5561 



#######################################
CD4set<-CD4set[,-12] #drop Tetramer
CD4set<-CD4set[,-9] #drop CD3
#######################################
table(CD4set$cohortID)

#    ES     JD    MON    NBD 
#346324 167876 370020  61672 
#EStablished JustDiagnosed MONitored NormalBloodDonor




X<-data.matrix(CD4set[,1:15]) #no Tetramer
Y<-X
features<-colnames(X)
nf<-length(features) #15

#standardize by subject, use 99% standardization per channel
for (it in 1:ns) {
   t<-subjects[it]
   selection<-CD4set$subjectID== t
   cat(t,sum(selection),"\n")
   for (f in features) {
      u<-X[selection,f]
      u99<-quantile(u,.99)
      u<-u/u99*q99
      Y[selection,f]<-u
   }
}





normCD4set<-CD4set
normCD4set[,1:15]<-Y
save(normCD4set,file="normCD4set.RData")
#these are normalized CD4+CD8- cells

#define quantTetramer - quantile of Tetramer signal
#what is quantile for every subject and tetramer?
normCD4set$quantTetramer<-rep(NA,dim(normCD4set)[1])
for (i in 1:ns) {
   for (tetra in tetras) {
      selection<-normCD4set$subjectID==subjects[i] & normCD4set$tetraID==tetra
      if (length(selection)>0) normCD4set[selection,"quantTetramer"]<-rank(normCD4set[selection,"qTetramer"])/sum(selection)
      cat(subjects[i],tetra,"\n")
   }
}
discard<-is.na(normCD4set$quantTetramer)
normCD4set<-normCD4set[!discard,]
#these are normalized CD4+CD8- cells!!!

dim(normCD4set)
#[1] 918424      22
names(normCD4set)
# [1] "HLADR"         "CXCR4"         "CD28"          "CXCR3"        
# [5] "CXCR5"         "CD127"         "CD69"          "TIM3"         
# [9] "CD45RO"        "PD1"           "CTLA4"         "CD25"         
#[13] "ICOS"          "LAG3"          "CD62L"         "tetraPos"     
#[17] "tetraID"       "subjectID"     "cohortID"      "qTetramer"    
#[21] "qtetraPos"     "quantTetramer"

normCD4set[,"qtetraPos"]<-normCD4set$quantTetramer>0.997
save(normCD4set,file="normCD4set.RData")

#make useful scatterplots
pdf("CD62L_CD45RO.pdf",width=4*4+4,height=ns*4+4)
par(mfrow=c(ns,4),pty="s",mar=c(4, 1, 2, 1),oma = c(4, 4, 4, 2))
for (it in 1:ns) {
   t<-subjects[it]
   for (tetra in tetras) {
      selection<-normCD4set$subjectID==t & normCD4set$tetraID==tetra
      pos<-selection & normCD4set$quantTetramer>0.997
      cat(t,sum(selection),sum(pos),"\n")
      if (sum(selection)>0) { #data exists
         plot(Y[selection,c("CD62L","CD45RO")],pch=".",col="#00000066",main=paste(longcohorts[t],t,tetra,sep=" "),xlim=c(-.2,4),ylim=c(-.2,4))
         abline(h=0,v=0,col="grey")
         points(Y[pos,c("CD62L","CD45RO")],pch=19,cex=1,col="#FF0000")
      } else plot.new()
   }
}
dev.off()

pdf("CD127_CD25.pdf",width=4*4+4,height=ns*4+4)
par(mfrow=c(ns,4),pty="s",mar=c(4, 1, 2, 1),oma = c(4, 4, 4, 2))
for (it in 1:ns) {
   t<-subjects[it]
   for (tetra in tetras) {
      selection<-normCD4set$subjectID==t & normCD4set$tetraID==tetra
      pos<-selection & normCD4set$quantTetramer>0.997
      cat(t,sum(selection),sum(pos),"\n")
      if (sum(selection)>0) { #data exists
         plot(Y[selection,c("CD127","CD25")],pch=".",col="#00000066",main=paste(longcohorts[t],t,tetra,sep=" "),xlim=c(-.2,4),ylim=c(-.2,10))
         abline(h=0,v=0,col="grey")
         points(Y[pos,c("CD127","CD25")],pch=19,cex=1,col="#FF0000")
      } else plot.new()
   }
}
dev.off()

pdf("CXCR3_CXCR4.pdf",width=4*4+4,height=ns*4+4)
par(mfrow=c(ns,4),pty="s",mar=c(4, 1, 2, 1),oma = c(4, 4, 4, 2))
for (it in 1:ns) {
   t<-subjects[it]
   for (tetra in tetras) {
      selection<-normCD4set$subjectID==t & normCD4set$tetraID==tetra
      pos<-selection & normCD4set$quantTetramer>0.997
      cat(t,sum(selection),sum(pos),"\n")
      if (sum(selection)>0) { #data exists
         plot(Y[selection,c("CXCR3","CXCR4")],pch=".",col="#00000066",main=paste(longcohorts[t],t,tetra,sep=" "),xlim=c(-1,4),ylim=c(-.2,4))
         abline(h=0,v=0,col="grey")
         points(Y[pos,c("CXCR3","CXCR4")],pch=19,cex=1,col="#FF0000")
      } else plot.new()
   }
}
dev.off()

pdf("CD25_ICOS.pdf",width=4*4+4,height=ns*4+4)
par(mfrow=c(ns,4),pty="s",mar=c(4, 1, 2, 1),oma = c(4, 4, 4, 2))
for (it in 1:ns) {
   t<-subjects[it]
   for (tetra in tetras) {
      selection<-normCD4set$subjectID==t & normCD4set$tetraID==tetra
      pos<-selection & normCD4set$quantTetramer>0.997
      cat(t,sum(selection),sum(pos),"\n")
      if (sum(selection)>0) { #data exists
         plot(Y[selection,c("CD25","ICOS")],pch=".",col="#00000066",main=paste(longcohorts[t],t,tetra,sep=" "),xlim=c(-1,8),ylim=c(-.2,8))
         abline(h=0,v=0,col="grey")
         points(Y[pos,c("CD25","ICOS")],pch=19,cex=1,col="#FF0000")
      } else plot.new()
   }
}
dev.off()

