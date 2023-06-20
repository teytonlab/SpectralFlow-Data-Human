#save(redCountSummary,fractionRedCountSummary,cohortRedCountSummary,cohortFractionRedCountSummary,file="tetra_Summaries.RData")
load("tetra_Summaries.RData")

##################################begin definitions section######################################################
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}
#my modified heatmap.2 source code, with OLO
source("/media/rs_DATA/s3/20210111_varner_actri_voucher/heatmap.3_original.R")
####################################end definitions section######################################################


features<-c("HLADR","CXCR4","CD28","CXCR3","CXCR5","CD127","CD69","TIM3","CD45RO","PD1","CTLA4","CD25","ICOS","LAG3","CD62L")
tetras<-c("12-20","13-21","CP11","CP18")
#these subjects are ordered by cohort
subjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005","TC017","T012",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD008","RAD009","RAD010","RAD011",
            "TNET006","TNET007","TNET008","TNET009","TNET010","TNET011","TNET012","TNET013","TNET014","TNET015")
cohorts<-c("NBD","ES","JD","MON")

longcohorts<-rep(cohorts,c(7,5,5,10))
names(longcohorts)<-subjects


#define "clean" subjects, etc.   
ccohorts<-c("NBD","ES","JD")
csubjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD008","RAD009","RAD010","RAD011")

ncs<-length(csubjects)
clongcohorts<-rep(ccohorts,c(5,5,5))

colors3<-c("black","red","pink")
names(colors3)<-ccohorts

colors7<-c("grey70","darkolivegreen","brown","purple","magenta","blue","cyan","grey20","orange","red")
clcolors<-makeTransparent(colors7,alpha=0.08) #this covers A, Aa under one color


msubjects<-c("TNET006","TNET007","TNET008","TNET009","TNET010","TNET011","TNET012","TNET013","TNET014","TNET015")
nm<-length(msubjects)


colors3<-c("black","red","pink")
ccolors<-rep(colors3,each=5)
names(ccolors)<-csubjects


load("tetras_V.RData") #loads V
#V is essentially the same as V but with extra UMAPs and tetraclasses
names(V)
 [1] "UMAP 1"        "UMAP 2"        "HLADR"         "CXCR4"        
 [5] "CD28"          "CXCR3"         "CXCR5"         "CD127"        
 [9] "CD69"          "TIM3"          "CD45RO"        "PD1"          
[13] "CTLA4"         "CD25"          "ICOS"          "LAG3"         
[17] "CD62L"         "tetraPos"      "tetraID"       "subjectID"    
[21] "cohortID"      "qTetramer"     "qtetraPos"     "UMAP.1"       
[25] "UMAP.2"        "classes"       "classlabels"   "pointcolors"  
[29] "quantTetramer" "tetraclasses" 


colors4<-c("grey70","deepskyblue4","deeppink3","goldenrod1")
   cl<-unique(V$tetraclasses)
   ncl<-length(cl) #4
   rge<-range(V[,1:2])
labels<-paste0("T",1:ncl)



newfeatures<-c("CD25","ICOS","LAG3","CD127","PD1","HLADR","CD69","CD62L","CD45RO","CD28","TIM3","CXCR4","CXCR3","CXCR5","CTLA4")
#ridgeplots please
trcolors<-makeTransparent(colors4,alpha=0.8)
names(trcolors)<-labels

nf<-length(features)


#I need all samples, their tetras
load("normCD4set.RData")
names(normCD4set)
 [1] "HLADR"         "CXCR4"         "CD28"          "CXCR3"        
 [5] "CXCR5"         "CD127"         "CD69"          "TIM3"         
 [9] "CD45RO"        "PD1"           "CTLA4"         "CD25"         
[13] "ICOS"          "LAG3"          "CD62L"         "tetraPos"     
[17] "tetraID"       "subjectID"     "cohortID"      "qTetramer"    
[21] "qtetraPos"     "quantTetramer"

TL<-normCD4set[normCD4set$quantTetramer>0.997,]
dim(TL)
[1] 3370   22

Nsubjects<-csubjects[1:5]
Esubjects<-csubjects[6:10]
Jsubjects<-csubjects[11:15]


#this will work as a classifier. 
#Jensen-Shannon entropy
library(entropy)

jensen_shannon<-function(x,y) {
   u<-sqrt(entropy((x+y)/2)-entropy(x)/2-entropy(y)/2)
   return(u)
}

d3<-fractionRedCountSummary
#collapse clusters T3+T4
d3<-fractionRedCountSummary[1:3,]
d3[3,]<-colSums(fractionRedCountSummary[3:4,])

D<-matrix(0,nrow=ncs,ncol=ncs)
for (i in 1:(ncs-1)) {
   for (j in (i+1):ncs) {
      D[i,j]<-D[j,i]<-jensen_shannon(d3[,i],d3[,j])
   }
}
extendedD<-matrix(0,nrow=15,ncol=25)
extendedD[1:15,1:15]<-D

library(uwot)
library(plot3D)
set.seed(23571234)
u2<-umap(D,n_components=2,init="random")
pdf("cUMAP.pdf")
   plot(u2,col=ccolors,pch=19,cex=2,xlab="UMAP 1",ylab="UMAP 2")
dev.off()
#pretty nice. A 3nn approach might work well, with 2 mis-classifications, 13% error rate

u3<-umap(D,n_components=3)
save(u3,file="umap3.RData")
scatter3D(u3[,1],u3[,2],u3[,3],pch=19,cex=1,col=ccolors,phi=10,theta=0,colvar=NULL)
#nice view with almost separation
scatter3D(u3[,1],u3[,2],u3[,3],pch=19,cex=1,col=ccolors,phi=-30,theta=50,colvar=NULL)
#especially
pdf("UMAP_training_set_3D.pdf")
   scatter3D(u3[,1],u3[,2],u3[,3],pch=19,cex=1.4,col=ccolors,phi=-55,theta=170,colvar=NULL,xlab="UMAP 1",ylab="UMAP 2",zlab="UMAP 3")
   flr<- min(u3[,2])
   for (i in 1:dim(u3)[1]) {
      lines3D(c(u3[i,1],u3[i,1]),c(flr,u3[i,2]),c(u3[i,3],u3[i,3]),add=TRUE,colvar=NULL,col=ccolors[i])
   }
dev.off()




#leave one out from every cohort, propagate labels, calculate cluster fractions t1..t4, calculate Jensen-Shannon distance to clean samples, csamples
#then do UMAP, with 3nn graph, and find classification by majority vote.

#also, try label propagation into new cells where labels are subjectID? No. use cohortID labels
table(TL$cohortID)

  ES   JD  MON  NBD  OUT 
1063  582 1048  196  481 




#this is collected from old code with label propagation
library(igraph)
library(rnndescent)


pdf("mUMAPs.pdf")
cY<-V[,1:15+2]
MON_fractions<-matrix(0,nrow=4,ncol=length(msubjects))
colnames(MON_fractions)<-msubjects
for (im in 1:10) {
   m<-msubjects[im]
   mY<-normCD4set[normCD4set$subjectID==m & normCD4set$quantTetramer>0.997,1:15] #take tetra+ cells from subject m
   Y<-rbind(cY,mY) #has csubjects on top, and a new msubject at the bottom
   c_m_subjects<-c(V$subjectID,normCD4set$subjectID[normCD4set$subjectID==m & normCD4set$quantTetramer>0.997])

   k<-15
   system.time( gnn<-nnd_knn(Y,k=k,metric="euclidean") )


   el<-gnn$idx[,1:2]
   swap<-el[,1]>el[,2]
   el[swap,]<-el[swap,2:1] #order pairs
   sel<-unique(el)
   g0<-graph_from_edgelist(sel, directed = FALSE)
   for (j in 3:k) {
      el<-gnn$idx[,c(1,j)]
      swap<-el[,1]>el[,2]
      el[swap,]<-el[swap,2:1] #order pairs
      sel<-unique(el)
      g1<-graph_from_edgelist(sel, directed = FALSE)
      g0<-union(g0,g1)
   }
   n<-length(V(g0))
   initial_labels<-rep(-1,n) #negative values mean there is no assigned node label
   names(initial_labels)<-1:n
   originalcomm<-V$tetraclasses 
#  originalcomm[originalcomm==4]<-3 #don't do this if you want to use all 4 fractions. I just think that first 2 + dumpster are better
   initial_labels[1:dim(cY)[1]]<-originalcomm
   all_labels<-cluster_label_prop(g0,initial=initial_labels,fixed=initial_labels>=0)
   final_labels<-all_labels$membership
   table(all_labels$membership)

   #how were the original labels renamed?
   newT1<-unique(final_labels[initial_labels==1]) #is T1
   newT2<-unique(final_labels[initial_labels==2]) #is T2
   newT3<-unique(final_labels[initial_labels==3]) #is T3
   newT4<-unique(final_labels[initial_labels==4]) #is T4
   newnames<-c(newT1,newT2,newT3,newT4)
   
   redCountSummary<-matrix(0,nrow=4,ncol=16) #16 is  clean subjects + 1; 4 is clusters
   rownames(redCountSummary)<-1:4
   colnames(redCountSummary)<-c(csubjects,m)
   for (s in c(csubjects,m)) {
      for (icl in 1:4) {
         redCountSummary[icl,s]<-sum(c_m_subjects==s & final_labels==newnames[icl])
      }
   }
   fractionRedCountSummary<-t(t(redCountSummary)/colSums(redCountSummary))
   MON_fractions[,m]<-fractionRedCountSummary[,16]
   
   d3<-fractionRedCountSummary[1:3,] #combine T3 U T4
   d3[3,]<-colSums(fractionRedCountSummary[3:4,])
   
   D<-matrix(0,nrow=16,ncol=16)
   for (i in 1:(16-1)) {
      for (j in (i+1):16) {
         D[i,j]<-D[j,i]<-jensen_shannon(d3[,i],d3[,j])
      }
   }
   set.seed(2357)
   u2<-umap(D,n_components=2,init="random")
   plot(u2,col=c(ccolors,"cyan"),pch=19,cex=2,xlab="UMAP 1",ylab="UMAP 2",main=m)
   extendedD[,15+im]<-D[1:15,16]
} #m
dev.off()
good_subjects<-c(csubjects,msubjects)
colnames(extendedD)<-good_subjects
gcolors<-c(ccolors,rep("cyan",10))
save(MON_fractions,file="MON_fractions.RData")

pdf("MON_tetra_fractions_by_cluster.pdf",width=10,height=5)
par(mfrow=c(1,2),pty="s",mar=c(8, 1, 2, 1))
barplot(MON_fractions,col=colors4,las=2,border=NA,main=expression(tetramer^"+"~cells),ylab="population fraction")
dev.off()


pdf("tetra_cluster_fractions.pdf",width=5,height=4)
dx<-0:2
dx<-c(0,2,1) #this order was requested
names(dx)<-ccohorts
barplot(t(cohortFractionRedCountSummary[,c(1,3,2)]),beside=TRUE,border=NA,col=colors3[c(1,3,2)],ylab="population fraction",xlab="cluster",ylim=c(0,.9))
legend("topright",legend=ccohorts[c(1,3,2)],pch=15,col=colors3[c(1,3,2)])
for (coh in ccohorts) {
   fr<-fractionRedCountSummary[,clongcohorts==coh]
   for (i in 1:ntc) {
      x<-(i-1)*4+1+dx[coh]+1/2+rnorm(dim(fr)[2],sd=.1)
      y<-fr[i,]
      points(x,y,cex=0.5,col="grey")
   }
}
box()
dev.off()


#plot by distance
library(shape)
symbols<-rep("|",15)
names(symbols)<-csubjects
#symbols["T033"]<-symbols["RAD011"]<-"0"
pdf("JSdistances.pdf")
plot(NA,NA,xlim=c(-.1,.6),ylim=c(0,15+10+1),main="",xlab="Jensen-Shannon distance",ylab="",axes=FALSE)
for (s in good_subjects) {
   d<-extendedD[,s]
   od<-order(d)
   j<-which(good_subjects==s)
   points(d,rep(j,15),col=ccolors,pch=symbols)
   text(-.01,j-.2,labels=s,adj=c(1,0),col=gcolors[j],cex=0.8)
}
Arrows(0,0,.6,0,arr.type="triangle",col="grey",arr.length=0.25)
points(seq(0,.5,.1),rep(0,6),pch="|",cex=.5,col="grey")
text(seq(0,.5,.1),rep(0,6),labels=round(seq(0,.5,.1),1),pos=1,col="grey",cex=0.8)
dev.off()

symbols<-rep(17,15) #triangle, full
names(symbols)<-csubjects
pdf("JS_3nn_distances.pdf")
plot(NA,NA,xlim=c(-.1,.2),ylim=c(0,15+10+1),main="",xlab="Jensen-Shannon distance",ylab="",axes=FALSE)
for (s in msubjects) {
   d<-extendedD[,s]
   od<-order(d)
   j<-which(good_subjects==s) - 15 #that 10 is there to skip the ccubjects 
   lines(c(0,0.205),c(j-.15,j-.15),col="grey",lwd=.2)
   points(d[od[3:1]],rep(j,3),col=ccolors[od[3:1]],pch=symbols) #order 3:1 is there to plot the closest one unobscured
   text(-.01,j-.2,labels=s,adj=c(1,0),col="cyan",cex=0.8)
}
Arrows(0,1-.15,.205,1-.15,arr.type="triangle",col="grey",arr.length=0.25,lwd=0.2)
points(seq(0,.5,.1),rep(1-.15,6),pch="|",cex=.5,col="grey")
text(seq(0,.5,.1),rep(1-.15,6),labels=round(seq(0,.5,.1),1),pos=1,col="grey",cex=0.8)
dev.off()





#classifier power analysis
#withhold samples from the training set and classify them by label propagation to cells from retained training samples, then call the class by majority of 3nn. If undecided, use closest.

cY<-V[,1:15+2]
#at first, withhold one
library(Rfast)

for (outsubject in csubjects) {
   insubjects<-setdiff(csubjects,outsubject)
   inY<-cY[V$subjectID %in% insubjects,] #"training set"
   outY<-cY[V$subjectID==outsubject,] #subject to classify
   atmsubjects<-c(V$subjectID[V$subjectID %in% insubjects],V$subjectID[V$subjectID==outsubject])


   Y<-rbind(inY,outY) #has csubjects on top, and a new msubject at the bottom

   k<-15
   system.time( gnn<-nnd_knn(Y,k=k,metric="euclidean") )


   el<-gnn$idx[,1:2]
   swap<-el[,1]>el[,2]
   el[swap,]<-el[swap,2:1] #order pairs
   sel<-unique(el)
   g0<-graph_from_edgelist(sel, directed = FALSE)
   for (j in 3:k) {
      el<-gnn$idx[,c(1,j)]
      swap<-el[,1]>el[,2]
      el[swap,]<-el[swap,2:1] #order pairs
      sel<-unique(el)
      g1<-graph_from_edgelist(sel, directed = FALSE)
      g0<-union(g0,g1)
   }
   n<-length(V(g0))
   initial_labels<-rep(-1,n) #negative values mean there is no assigned node label
   names(initial_labels)<-1:n
   originalcomm<-V$tetraclasses[V$subjectID %in% insubjects]
   originalcomm[originalcomm==4]<-3 #don't do this if you want to use all 4 fractions. I just think that first 2 + dumpster are better
   initial_labels[1:dim(inY)[1]]<-originalcomm
   all_labels<-cluster_label_prop(g0,initial=initial_labels,fixed=initial_labels>=0)
   final_labels<-all_labels$membership

   redCountSummary<-matrix(0,nrow=3,ncol=15) #16 is  clean subjects + 1; 4 is clusters
   rownames(redCountSummary)<-1:3
   colnames(redCountSummary)<-c(insubjects,outsubject)
   for (s in csubjects) {
      for (icl in 1:3) {
         redCountSummary[icl,s]<-sum(atmsubjects==s & final_labels==icl)
      }
   }
   fractionRedCountSummary<-t(t(redCountSummary)/colSums(redCountSummary))

   d3<-fractionRedCountSummary

   D<-matrix(0,nrow=15,ncol=15)
   for (i in 1:(15-1)) {
      for (j in (i+1):15) {
         D[i,j]<-D[j,i]<-jensen_shannon(d3[,i],d3[,j])
      }
   }
   d<-D[15,-15]
   o<-order(d)
   names(d)<-longcohorts[insubjects]
   #to classify 
   tdf<-as.data.frame(table(names(d)[o[1:3]]))
   mf<-max(tdf$Freq)
   if (mf>1) {
      newclass<-as.character(tdf$Var1)[which.max(tdf$Freq)] #the majority
   } else {
      newclass<-names(d)[o[1]] #the closest
   }
   #cat(outsubject,d[o[1:3]],names(d)[o[1:3]],"\n")
   cat(outsubject,newclass,"\n")
   
} 



#OK, now program all combinations of outsamples such that at most 3 are missing from any cohort - that's between 1 and 9 out
#so just generate all such combinations, N over m, with m between 1 and 9, under constraint that there may be at most 3 missing
for (m in 2:9) {
   #generate combinations of missing ones:
   w<-combn(csubjects,m)
   #remove those that leave fewer than 2 in a cohort
   good<-rep(FALSE,dim(w)[2])
   for (j in 1:dim(w)[2]) {
      outcohorts<-longcohorts[w[,j]] 
      ddf<-as.data.frame(table(outcohorts))
      good[j]<-max(ddf$Freq)<=3
   }
   cat(m,sum(good),mean(good),"\n")
   #now do the work
   gw<-w[,good] #reduce space of combinations to those that leave at least 2 in a cohort
   #what are the original "true" classes?
   trueclasses<-matrix(longcohorts[gw],nrow=m)
   newclasses<-trueclasses #for size
   for (J in 1:dim(gw)[2]) {
      outsubjects<-gw[,J] #these are out, must classify them one by one
      insubjects<-setdiff(csubjects,outsubjects)
      for (l in 1:m) {
         os<-outsubjects[l]

         inY<-cY[V$subjectID %in% insubjects,] #"training set"
         outY<-cY[V$subjectID==os,] #subject to classify
         atmsubjects<-c(V$subjectID[V$subjectID %in% insubjects],V$subjectID[V$subjectID==os])

         Y<-rbind(inY,outY) #has csubjects on top, and a new msubject at the bottom
         k<-15
         gnn<-nnd_knn(Y,k=k,metric="euclidean")
         el<-gnn$idx[,1:2]
         swap<-el[,1]>el[,2]
         el[swap,]<-el[swap,2:1] #order pairs
         sel<-unique(el)
         g0<-graph_from_edgelist(sel, directed = FALSE)
         for (j in 3:k) {
            el<-gnn$idx[,c(1,j)]
            swap<-el[,1]>el[,2]
            el[swap,]<-el[swap,2:1] #order pairs
            sel<-unique(el)
            g1<-graph_from_edgelist(sel, directed = FALSE)
            g0<-union(g0,g1)
         }
         n<-length(V(g0))
         initial_labels<-rep(-1,n) #negative values mean there is no assigned node label
         names(initial_labels)<-1:n
         originalcomm<-V$tetraclasses[V$subjectID %in% insubjects]
         originalcomm[originalcomm==4]<-3 #don't do this if you want to use all 4 fractions. I just think that first 2 + dumpster are better
         initial_labels[1:dim(inY)[1]]<-originalcomm
         all_labels<-cluster_label_prop(g0,initial=initial_labels,fixed=initial_labels>=0)
         final_labels<-all_labels$membership
      
         n15<-15-m+1 #number of subjects in the test
         redCountSummary<-matrix(0,nrow=3,ncol=n15) #16 is  clean subjects + 1; 4 is clusters
         rownames(redCountSummary)<-1:3
         colnames(redCountSummary)<-c(insubjects,os)
         for (s in c(insubjects,os)) {
            for (icl in 1:3) {
               redCountSummary[icl,s]<-sum(atmsubjects==s & final_labels==icl)
            }
         }
         fractionRedCountSummary<-t(t(redCountSummary)/colSums(redCountSummary))
         d3<-fractionRedCountSummary
         D<-rep(0,n15-1)
         for (i in 1:(n15-1)) D[i]<-jensen_shannon(d3[,i],d3[,n15])
         o<-order(D)
         names(D)<-longcohorts[insubjects]
         #to classify 
         tdf<-as.data.frame(table(names(D)[o[1:3]]))
         mf<-max(tdf$Freq)
         if (mf>1) {
            newclass<-as.character(tdf$Var1)[which.max(tdf$Freq)] #the majority
         } else {
            newclass<-names(D)[o[1]] #the closest
         }
         #cat(outsubject,D[o[1:3]],names(D)[o[1:3]],"\n")
         #cat(longcohorts[os],newclass,"\n")
         newclasses[l,J]<-newclass
      }
   }
   cat(m,mean(trueclasses==newclasses),"\n")
}

         
#1 15 1 
#2 105 1 
#3 455 1 
#4 1350 0.989011 
#5 2850 0.9490509 
#6 4300 0.8591409 
#7 4500 0.6993007 
#8 3000 0.4662005 
#9 1000 0.1998002 

#2 105 1 
#2 0.7857143 
#3 455 1 
#3 0.7611722 
#4 1350 0.989011 
#4 0.7355556 
#5 2850 0.9490509 
#5 0.7142456 
#6 4300 0.8591409 
#6 0.6981395 
#7 4500 0.6993007 
#7 0.6906349 
#8 3000 0.4662005 
#8 0.6920417 
#9 1000 0.1998002 
#9 0.7012222 

N<-15-(1:8)
#er<-1-c(.8,0.7857143,0.7611722,0.7355556,0.7142456,0.6981395,0.6906349,0.6920417,0.7012222)
er<-1-c(.8,0.7857143,0.7611722,0.7355556,0.7142456,0.6981395,0.6906349,0.6920417)
pdf("classifier_error_rate.pdf")
plot(N,er,ylim=c(.1,.35),xlim=c(5,25),ylab="error rate",xlab="N-m",pch=19,log="y")
abline(h=.1,col="grey")
dev.off()
