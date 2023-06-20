library(igraph)
library(matrixStats)
library(uwot) #for umap
library(leiden)

##################################begin definitions section##############################################
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

transf<-function(x,a) return(asinh(x/a))
####################################end definitions#######################################################

features<-c("HLADR","CXCR4","CD28","CXCR3","CXCR5","CD127","CD69","TIM3","CD45RO","PD1","CTLA4","CD25","ICOS","LAG3","CD62L")
tetras<-c("12-20","13-21","CP11","CP18")
#these subjects are ordered by cohort
subjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD009","RAD010","RAD011",
            "TNET006","TNET007","TNET008","TNET009","TNET010","TNET011","TNET012","TNET013","TNET014","TNET015")
cohorts<-c("NBD","ES","JD","MON")

longcohorts<-rep(cohorts,c(5,5,4,10))
names(longcohorts)<-subjects

colors4<-c("grey20","red","pink","magenta")
names(colors4)<-cohorts

#define "clean" subjects, etc.
ccohorts<-c("NBD","ES","JD")
csubjects<-c("NBD-BRI001","NBD-BRI002","NBD-BRI003","NBD-BRI004","NBD-BRI005",
            "T030","T031","T033","T034","T035",
            "RAD007","RAD009","RAD010","RAD011")

ncs<-length(csubjects)
clongcohorts<-rep(ccohorts,c(5,5,4))



load("normCD4set.RData")
cnormCD4set<-normCD4set[normCD4set$subjectID %in% csubjects,] #define clean CD4set
normCD4<-cnormCD4set  #use clean CD4set only
n<-dim(normCD4)[1]
#[1] 575872
X<-data.matrix(normCD4[,1:15])

cmass<-rep(0,3) #cohort mass
names(cmass)<-ccohorts
for (coh in ccohorts) cmass[coh]<-mean(normCD4$cohortID==coh)

##############start with global UMAP##########################
k<-15 #k was set to 15 by Roman
system.time( W<-uwot::umap(X=X,n_neighbors=k,nn_method="annoy",ret_extra=c("nn","fgraph")) )
#      user   system  elapsed
#1001.604    6.984  808.293
save(W,file="W.RData")
U<-W$embedding
colnames(U)<-c("UMAP 1","UMAP 2") # g for global UMAP
rownames(U)<-rownames(normCD4)
save(U,file="normCD4_U2.RData")
pdf("UMAP0.pdf")
rge<-range(U)
plot(U,pch=".",col=makeTransparent("black",alpha=0.02),cex=1,xlim=rge,ylim=rge,main="all cohorts")
for (coh in ccohorts) {
   selection<-normCD4$cohort==coh
   plot(U[selection,],pch=".",col=makeTransparent(colors4[coh],alpha=0.04*1/cmass[coh]*0.58),cex=1,xlim=rge,ylim=rge,main=coh) #make plots ~same
}
dev.off()
###############done with global UMAP##########################

colors7<-c("grey70","darkolivegreen","brown","purple","magenta","blue","cyan","grey20","orange","red")
#GLOBAL CLUSTERING
   gnn<-W$nn[[1]][[1]]

   el<-gnn[,c(1,2)] #start edgelist
   for (j in 3:k) el<-rbind(el,gnn[,c(1,j)])
   swap<-el[,1]>el[,2]
   el[swap,]<-el[swap,2:1] #this sorts edge end points so that smaller index goes first
   dup<-duplicated(el) #we will remove these duplicated edges when we build the graph
   elist<-matrix(as.character(el),nrow=dim(el)[1]) #convert node names to character
   g0<-graph_from_edgelist(elist[!dup,], directed = FALSE)
   # save(g0, file=paste0(t,"_g0.RData"))
   save(g0, file=paste0("t_g0.RData"))
   vnames<-as.integer(V(g0)$name)
   o<-order(vnames)

   #in Constant Potts Model, gamma has a direct interpretation: internal density (number of internal edges/choose(nc,2)) of communities is >= gamma
   #see Traag, V. A., Van Dooren, P., & Nesterov, Y. (2011). Narrow scope for resolution-limit-free community detection. Physical Review E, 84(1), 016114.
   #this is also resolution-limit-free! i.e., able to detect smaller clusters in big networks
   system.time( partition <- leiden(g0,partition_type="CPMVertexPartition",resolution_parameter=0.000014) ) #resolution_parameter will depend on the data
   #1991.892    5.124 1996.959
   classes<-partition[o]
   print(table(classes))
#     1      2      3      4      5      6      7      8 
#164333 160182 107076  92893  39414   6360   3075   2539 


   V<-data.frame(U,classes,check.names=FALSE)
   save(V,file="normCD4_V.RData")

pdf("gUMAP.pdf")
   cl<-unique(V$classes)
   ncl<-length(cl)
   rge<-range(V[,1:2])
   plot(V[,1:2],pch=".",cex=1,xlim=rge,ylim=rge,main="all cells",col=makeTransparent("black",alpha=0.02))

   plot(V[,1:2],pch=".",cex=0,xlim=rge,ylim=rge,main="all cells",col=makeTransparent("black",alpha=0.02))
   for (icl in 1:ncl) points(V[V$classes==icl,1:2],pch=".",cex=1,col=makeTransparent(colors7[icl],alpha=0.08))

   for (icl in 1:ncl) plot(V[V$classes==icl,1:2],pch=".",cex=1,col=makeTransparent(colors7[icl],alpha=0.08),main=c("A","B","C","D","E","F","G","H")[icl],xlim=rge,ylim=rge)
dev.off()

cnormCD4setW<-data.frame(cnormCD4set,V,check.names=FALSE)
save(cnormCD4setW,file="cnormCD4set_W.RData")


for (t in ccohorts) {
   cat(t,"\n")
   W<-cnormCD4setW[cnormCD4setW$cohortID==t,]
   save(W,file=paste0(t,"_cnormCD4_W.RData"))
}


trcolors<-makeTransparent(colors7,alpha=.8)
pdf("gUMAP_by_cohort_cluster.pdf")
for (t in ccohorts) {
   cat(t,"\n")
   load(paste0(t,"_cnormCD4_W.RData"))
   cl<-unique(W$classes)
   ncl<-length(cl)
   rge<-range(W[,23:24])
   plot(W[,23:24],pch=".",cex=0,xlim=rge,ylim=rge,main=t,col=makeTransparent("black",alpha=0.06))
   for (icl in 1:ncl) points(W[W$classes==icl,23:24],pch=".",cex=1,col=makeTransparent(colors7[icl],alpha=0.04*1/cmass[t]*0.58))
   print(table(W$classes))
}
dev.off()
#NBD 
#
#    1     2     3     4     5     6     7     8 
#24209 16249  9337  7946  2882   670   139   240 
#ES 
#
#     1      2      3      4      5      6      7      8 
#107704 109378  75592  19829  25406   3760   2893   1762 
#JD 
#
#    1     2     3     4     5     6     7     8 
#32420 34555 22147 65118 11126  1930    43   537 



names(W)
# [1] "HLADR"         "CXCR4"         "CD28"          "CXCR3"        
# [5] "CXCR5"         "CD127"         "CD69"          "TIM3"         
# [9] "CD45RO"        "PD1"           "CTLA4"         "CD25"         
#[13] "ICOS"          "LAG3"          "CD62L"         "tetraPos"     
#[17] "tetraID"       "subjectID"     "cohortID"      "qTetramer"    
#[21] "qtetraPos"     "quantTetramer" "UMAP 1"        "UMAP 2"       
#[25] "classes"      



newfeatures<-c("CD62L","CD25","CD45RO","CD28","CD69","LAG3","TIM3","HLADR","CXCR4","CXCR3","CXCR5","CD127","PD1","CTLA4","ICOS")

#after inspecting clusters and what makes them unique, define this
clcolors<-makeTransparent(colors7,alpha=0.08) #this covers A, Aa under one color

names(clcolors)<-c("A","B","C","D","E","F","G","H","I","J")
transtable.list<-list( 
"NBD"=c("A","B","C","D","E","F","G","H"),
"ES"=c("A","B","C","D","E","F","G","H"),
"JD"=c("A","B","C","D","E","F","G","H")
)

for (t in ccohorts) {
   cat(t,"\n")

   load(paste0(t,"_cnormCD4_W.RData"))
   colnames(W)[23:24]<-c("UMAP 1","UMAP 2")
   labels<-W$classes 
   classlabels<-transtable.list[[t]][W$classes] #this way, clusters outside of transtable will have a NA label
   pointcolors<-clcolors[classlabels]
   W<-data.frame(W,classlabels,pointcolors,check.names=FALSE) #adds global UMAP coordinates
   save(W,file=paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData"))
}


pdf("UMAPs_by_cohort.pdf",width=12,height=4)
   par(mar=c(4, 1, 2, 1), mfrow=c(1,3), oma = c(4, 4, 2, 2),pty="s")
   for (t in ccohorts) {
      cat(t,"\n")
      load(paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData"))
      rge<-range(W[,23:24])
      plot(W[,23:24],pch=".",cex=0,xlim=rge,ylim=rge,main=t,xlab="",ylab="")
      points(W[,23:24],pch=".",cex=1,col=W$pointcolors)
   }
dev.off()

pdf("UMAPs_w_tetras_by_cohort_grey_red.pdf",width=4*4+4,height=4*4)
   par(mar=c(4, 1, 2, 1), mfrow=c(4,3+1), oma = c(4, 4, 4, 2),pty="s")
for (tetra in tetras) {
   plot.new()
   text(0.5,0.5,labels=tetra,cex=4)
   for (t in ccohorts) {
      cat(t,"\n")
      load(paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData"))
   #modify tetraPos a little
   W$tetraPos[is.na(W$tetraPos)]<-""
      rge<-range(W[,23:24])
      plot(W[W$tetraID==tetra,23:24],pch=".",cex=0.1,xlim=rge,ylim=rge,main=t,col="#0000000F",xlab="UMAP 1",ylab="UMAP 2")
      points(W[W$tetraID==tetra & W$quantTetramer>0.997,23:24],pch=21,cex=0.55,col="red3",bg="red",lwd=0.2)
   }
   #title(tetra,line=1,cex=8,outer=TRUE)
}
dev.off()

pdf("UMAPs_w_combined_tetras_by_cohort_grey_red.pdf",width=4*4+4,height=4*4)
   par(mar=c(4, 1, 2, 1), mfrow=c(1,3), oma = c(4, 4, 4, 2),pty="s")
   for (t in ccohorts) {
      cat(t,"\n")
      load(paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData"))
      rge<-range(W[,23:24])
      plot(W[W$tetraID==tetra,23:24],pch=".",cex=0.1,xlim=rge,ylim=rge,main=t,col="#0000000F",xlab="UMAP 1",ylab="UMAP 2")
      points(W[W$quantTetramer>0.997,23:24],pch=21,cex=0.55,col="red3",bg="red",lwd=0.2)
   }
dev.off()

#track colors correspond to point colors, but different opacity
trcolors<-makeTransparent(colors7,alpha=0.8)
names(trcolors)<-c("A","B","C","D","E","F","G","H","I","J")

pdf("cl_ridgeplot.pdf",width=16,height=16)
par(mar=c(4, 1, 2, 1), mfrow=c(4,4), oma = c(4, 4, 2, 2),pty="s")

   #rowFtests
   partition<-c("A","B","C","D","E","F","G","H","I","J")[cnormCD4setW$classes] #here I am assuming that the numbering of clusters is the same in every cohort
   part.df<-as.data.frame(table(partition))
   nbig<-dim(part.df)[1]
   cls<-part.df$partition
   x<-cnormCD4setW
   #the following is good for F test (ANOVA), if desired. Here we just want to plot all genes and their distribution in clusters 
#  Ftest<-colFtests(x=data.matrix(x[,features]),fac=factor(x$classes,levels=bigcl))
#  Fstat<-Ftest$statistic
#  o<-order(-Fstat)
#  Ftest[o,]
      #reshape data into a list
      for (f in newfeatures) {
         dat.list<-list()
         for (cl in cls) dat.list[[cl]]<-x[partition==cl,f]
         #get density values
         vals <- Map(function(x, g, i) { with( density(x,n=1024,bw=.08), data.frame(x,y=y+(i-1), g) ) }, dat.list, names(dat.list), seq_along(dat.list))
         #now plot ridgeline plot
         xrange <- range(unlist(lapply(vals, function(d) range(d$x))))
         yrange <- range(unlist(lapply(vals, function(d) range(d$y))))
         #plot(0,0, type="n", xlim=c(-.4,4.), ylim=yrange, yaxt="n", ylab="", xlab="normalized intensity",main=f)
         plot(0,0, type="n", xlim=c(-1,4), ylim=c(0,10), yaxt="n", ylab="", xlab="normalized intensity",main=f)
         abline(v=0,col="grey")
         for (cl in rev(cls)) { #this guarantees that the clusters are picked in the right order on the y axis of ridgeplot
            #with(vals[[cl]], polygon(x, y, col=trcolors[transtable.list[[t]][cl]],border=NA))
            with(vals[[cl]], polygon(x, y, col=trcolors[cl],border=NA))
         }
         axis(2, at=seq_along(dat.list)-1, names(dat.list))
      }
      #while(!par('page')) plot.new()
dev.off()



#cluster signatures
cl_signatures<-list()
for (t in ccohorts) {
   cat(t,"\n")
   load(paste0(t,"_cnormCD4_W_classlabels_pointcolors.RData")) #W
   comms<-sort(unique(W$classlabels))
   Y<-W[,1:15]

   #store cluster expression signatures
   cl_means<-matrix(0,nrow=length(features),ncol=length(comms),dimnames=list(features,comms)) #add tetramerPos %
   for (l in comms) cl_means[features,l]<- colMeans(Y[W$classlabels==l,features],na.rm=TRUE)
   cl_signatures[[t]]<-cl_means
}
save(cl_signatures,file="CD4_signatures.RData")
