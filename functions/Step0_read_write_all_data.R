####Read in MetaData, establish which fluor measures which protein
MetaData_Hs <- read.csv("/Users/smsharma/Documents/OneDrive - Scripps Research/Fellowships-Grants-Invitations/ScienceTranslMed_paper-2021/Roman_Scripts_and_Manual_for_Teyton-Final/20220427.R/MetaData_Hs.csv", header = TRUE)
  #read.csv("/media/rs_DATA/s3/20210119_teyton_sharma_human_mouse_objective-3-diabetes-multiparametric-spectral-flow-normalization/analysis/human/MetaData_Hs.csv", header = TRUE)
MetaData_Hs
     Protein                    Color
1      HLADR                 BV421.A 
2      CXCR4             eFluor.450.A
3  Viability LIVE.DEAD.Blue.A        
4       CD28                  BV480.A
5      CXCR3                  BV510.A
6        CD8                  BV570.A
7      CXCR5       Super.Bright.645.A
8      CD127                  BV711.A
9       CD69                  BV785.A
10      TIM3                  BB515.A
11       CD3        Alexa.Fluor.532.A
12    CD45RO                  PerCP.A
13       PD1                  BB700.A
14  Tetramer                     PE.A
15       CD4           PE.Texas.Red.A
16     CTLA4                 PE.Cy5.A
17      CD25     PE.Alexa.Fluor.700.A
18      ICOS                 PE.Cy7.A
19      LAG3        Alexa.Fluor.647.A
20     CD62L                APC.Cy7.A



####Fluorophores, put the names into a vector so that you can attempt to loop over the fluorophores
Fluorophores <- as.character(MetaData_Hs$Color)
Proteins <- as.character(MetaData_Hs$Protein)

#define table of proteins with fluor names
fluors<-trimws(Fluorophores[-3])
proteins<-trimws(Proteins[-3])
#names(proteins)<-fluors


####Load csv Files
####there are NO RAD008, TC017 and _T012 files!
datadir<- "/Users/smsharma/Documents/OneDrive - Scripps Research/Laboratory Work/Projects T1D/Log-Aurora_Stained_Exported_Excel_Sheets_Human/Aurora_CSV_Files_Tetramer_Marked/"
  "/media/rs_DATA/s3/20210119_teyton_sharma_human_mouse_objective-3-diabetes-multiparametric-spectral-flow-normalization/analysis/human/csv/"
FileNames<- list.files(datadir, pattern = ".csv")

nf<-length(FileNames)

#from Sidd:
#MON: Monitored -->"TNET" donors
#ES: Established --> Adult donors with T1D
#JD: Just Diagnosed --> Child donors immediately post-diagnosis with T1D
#NBD: Normal Blood Donor --> Donors with risk of T1D (first degree relative) but do not have it.

#AD: Adult
#CH: Child

all_data<-list()
mass<-rep(0,nf)
names(mass)<-FileNames
for (f in FileNames) {
   cat(f,"\n")
   subject<-strsplit(f,"_")[[1]][9]
   t<-strsplit(f,"_")[[1]][7]
   tetra<-strsplit(f,"_")[[1]][10]

   D<-read.csv(file = paste0(datadir, f), header = TRUE)
   tetraPos<-D[,"Tetramer_Positive"]
   tetraPos<-!is.na(tetraPos)
   ncells<-dim(D)[1]
   tetraID<-rep(tetra,ncells)
   subjectID<-rep(subject,ncells)
   cohortID<-rep(t,ncells)
   D2<-D[,fluors]
   colnames(D2)<-proteins
   df<-data.frame(D2,tetraPos,tetraID,subjectID,cohortID,row.names=paste0(subject,"_",t,"_",tetra,".",1:length(tetraPos)),check.names=FALSE)
   all_data[[paste(subject,t,tetra,sep="_")]]<-df

   mass[f]<-dim(df)[1]
}
save(all_data,file=paste0("all_data.RData"))
save(mass,file="mass.RData")

summary(mass)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1788    6113   11212   15833   23314   47819 


H<-do.call("rbind",all_data) #into one data frame
dim(H)
#[1] 1519965      23

save(H,file="all_data_H.RData")
#this ^ includes CD8+ cells
names(H)
 #[1] "HLADR"     "CXCR4"     "CD28"      "CXCR3"     "CD8"       "CXCR5"    
 #[7] "CD127"     "CD69"      "TIM3"      "CD3"       "CD45RO"    "PD1"      
#[13] "Tetramer"  "CD4"       "CTLA4"     "CD25"      "ICOS"      "LAG3"     
#[19] "CD62L"     "tetraPos"  "tetraID"   "subjectID" "cohortID" 

