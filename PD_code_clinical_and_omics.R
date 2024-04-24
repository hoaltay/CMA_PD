#PD  heatmap
library(readxl)
library(VennDiagram)
library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(ggpubr)
############################################################################################################
##################################################heatmap: diff clinical indicators ########################
#heatmap for clinical indicators------------visit 2vs1 and visit 3vs1
rm(list=ls())
rm(list=ls())
rm(list=ls())

sig_th=0.05
disease<-"PD"

path_raw<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/output_from_Arif/Results/Clinical/"
setwd(path_raw)
int_data<-as.data.frame(read_excel('Clinical_2centers_20210201.xlsx', sheet =disease))
int_name<-as.matrix(int_data$Measurements)
log2FC_matrix<-cbind(int_data$`Log2FoldChange (Visit2 vs Visit1 - Active)`,
                    int_data$`Log2FoldChange (Visit3 vs Visit1 - Active)`,
                    int_data$`Log2FoldChange (Visit2 vs Visit1 - Placebo)`,
                    int_data$`Log2FoldChange (Visit3 vs Visit1 - Placebo)`)
colnames(log2FC_matrix)<-c("log2FC_V2vsV1_Active",
                           "log2FC_V3vsV1_Active",
                           "log2FC_V2vsV1_Placebo",
                           "log2FC_V3vsV1_Placebo")
mode(log2FC_matrix)<-"numeric"
index_na_1=which(is.na(log2FC_matrix))
if(length(index_na_1)!=0){
  log2FC_matrix[index_na_1]=0
}

sig_matrix<-cbind(int_data$`P value (Visit2 vs Visit1 - Active)`,
                  int_data$`P value (Visit3 vs Visit1 - Active)`,
                  int_data$`P value (Visit2 vs Visit1 - Placebo)`,
                  int_data$`P value (Visit3 vs Visit1 - Placebo)`)

colnames(sig_matrix)<-c("P_V2vsV1_Active",
                        "P_V3vsV1_Active",
                        "P_V2vsV1_Placebo",
                        "P_V3vsV1_Placebo")
mode(sig_matrix)<-"numeric"
index_na_2=which(is.na(sig_matrix))
if(length(index_na_2)!=0){
  sig_matrix[index_na_2]=1
}

index_1<-which(sig_matrix[,"P_V2vsV1_Active"]<sig_th)
index_2<-which(sig_matrix[,"P_V3vsV1_Active"]<sig_th)
index_3<-which(sig_matrix[,"P_V2vsV1_Placebo"]<sig_th)
index_4<-which(sig_matrix[,"P_V3vsV1_Placebo"]<sig_th)
index<-unique(c(index_1,index_2,index_3,index_4))

int_name<-as.matrix(int_name[index,])
#int_name[which(int_name=="Bel ?evresi")]="Waist circumference"# for AD
int_name[which(int_name=="Kalça Çevresi")]="Hip circumference"# for PD
log2FC_matrix<-log2FC_matrix[index,]
sig_matrix<-sig_matrix[index,]

label_matrix<-matrix("",dim(sig_matrix)[1],dim(sig_matrix)[2])
index_sig<-which(sig_matrix<sig_th)
label_matrix[index_sig]<-"*"

range(log2FC_matrix)
my.breaks <- c(seq(-0.5, -0.0001, by=0.01), seq(0, 0.5, by=0.01))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors = c("white","red"))(length(my.breaks)/2))


path_out<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/Figures/clinical_data/heatmap_all/"
setwd(path_out)

tiff(filename = paste0("heatmap_clincial_data_",disease,"2024.tiff"),units="in",width=4.5,height=7,res=400)
pheatmap(log2FC_matrix,
         cluster_cols = F,
         cluster_rows = F,
         #clustering_method = "ward.D2",
         color =  my.colors, 
         breaks = my.breaks,
         legend_breaks = seq(-0.5,0.5,by=0.5),
         display_numbers = label_matrix,
         fontsize_number = 20,
         fontsize = 15,
         treeheight_row =90,
         cellwidth = 20,
         cellheight = 20,
         gaps_col = 2,
         labels_col = c("CMA: Visit 2vs1","CMA: Visit 3vs1","Placebo: Visit 2vs1","Placebo: Visit 3vs1"),
         labels_row = int_name)
dev.off()

##################################################################################################
##################################################heatmap: diff metabolites########################
#PD
#extract significantly diff met----------amino acid, lipid ans others are separately
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
fdr_th=0.05

setwd("C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Kinetic_project\\raw_data\\")
trans_info<-data.frame(read_excel("KGCO-03-20ML DATA TABLES.XLSX",sheet="Chemical Annotation"))

path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\output_from_Arif\\Results\\Metabolomics\\"
setwd(path_raw)
int_data<-as.matrix(read.csv(paste0(disease,"_ALL_mets_Visit.txt"),header=F,row.names = 1,sep="\t"))
colnames(int_data)<-int_data[1,]
int_data<-int_data[-1,]
mode(int_data)="numeric"

loc<-match(rownames(int_data),trans_info$CHEMICAL_NAME)
which(is.na(loc))
trans_info<-cbind(trans_info$SUPER_PATHWAY[loc],
                  trans_info$SUB_PATHWAY[loc],
                  trans_info$CHEMICAL_NAME[loc])
colnames(trans_info)<-c("SUPER_PATHWAY","SUB_PATHWAY","CHEMICAL_NAME")

int_data<-int_data[,c("LogFC (Active - Visit 3 vs Active - Visit 1)",
                      "LogFC (Placebo - Visit 3 vs Placebo - Visit 1)",
                      "FDR (Active - Visit 3 vs Active - Visit 1)",
                      "FDR (Placebo - Visit 3 vs Placebo - Visit 1)")]
index_sig<-which(int_data[,"FDR (Active - Visit 3 vs Active - Visit 1)"]<fdr_th|int_data[,"FDR (Placebo - Visit 3 vs Placebo - Visit 1)"]<fdr_th)

#index_common<-which(int_data[,"FDR (Active - Visit 3 vs Active - Visit 1)"]<fdr_th&int_data[,"FDR (Placebo - Visit 3 vs Placebo - Visit 1)"]<fdr_th)
#index_ext_2<-which(int_data[,"FDR (Active - Visit 3 vs Active - Visit 1)"]>=fdr_th&int_data[,"FDR (Placebo - Visit 3 vs Placebo - Visit 1)"]<fdr_th)
#index_ext_3<-which(int_data[,"FDR (Active - Visit 3 vs Active - Visit 1)"]<fdr_th&int_data[,"FDR (Placebo - Visit 3 vs Placebo - Visit 1)"]>=fdr_th)

int_data<-int_data[index_sig,]
trans_info<-trans_info[index_sig,]

index_amino<-which(trans_info[,"SUPER_PATHWAY"]=="Amino Acid")
index_lipid<-which(trans_info[,"SUPER_PATHWAY"]=="Lipid")
index_others<-setdiff(c(1:dim(trans_info)[1]),c(index_amino,index_lipid))

data_amino<-int_data[index_amino,]
data_lipid<-int_data[index_lipid,]
data_others<-int_data[index_others,]

info_amino<-trans_info[index_amino,]
info_lipid<-trans_info[index_lipid,]
info_others<-trans_info[index_others,]

path_out<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\sig_met_all\\"
path_out<-paste0(path_out,disease,"\\")
setwd(path_out)
save(file=paste0(disease,"_sig_met_amino_lipid_others.Rdata"),data_amino,data_lipid,data_others,info_amino,info_lipid,info_others)
########################################################
####################
#heatmap for PD---amino acid
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
fdr_th=0.05
int_super_pathway<-"amino_acid"

path_raw<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/Figures/metabolics/sig_met_all/"
path_raw<-paste0(path_raw,disease,"/")
setwd(path_raw)
load(paste0(disease,"_sig_met_amino_lipid_others.Rdata"))
int_data<-data_amino#----------change

log2FC_matrix<-int_data[,c("LogFC (Active - Visit 3 vs Active - Visit 1)",
                           "LogFC (Placebo - Visit 3 vs Placebo - Visit 1)")]
sig_matrix<-int_data[,c("FDR (Active - Visit 3 vs Active - Visit 1)",
                        "FDR (Placebo - Visit 3 vs Placebo - Visit 1)")]

label_matrix<-matrix("",dim(sig_matrix)[1],dim(sig_matrix)[2])
index_sig<-which(sig_matrix<fdr_th)
label_matrix[index_sig]<-"*"

range(log2FC_matrix)
my.breaks <- c(seq(-2, -0.0001, by=0.01), seq(0, 2, by=0.01))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors = c("white","red"))(length(my.breaks)/2))


tiff(filename = paste0(int_super_pathway,"_heatmap_",disease,"_visit3vs1_2024.tiff"),width=7.3,height=11,units="in",res=400)
pheatmap(log2FC_matrix,
         cluster_cols = F,
         cluster_rows = T,
         clustering_method = "ward.D2",
         color =  my.colors, 
         breaks = my.breaks,
         legend_breaks = seq(-2,2,by=1),
         display_numbers = label_matrix,
         fontsize_number = 22,
         fontsize = 14,
         treeheight_row =90,
         cellwidth = 25,
         cellheight = 18,
         #gaps_col = 2,
         labels_col = c("CMA","Placebo"),
         labels_row = rownames(log2FC_matrix))
dev.off()
####################
#heatmap for PD---lipid
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
fdr_th=0.05
int_super_pathway<-"lipid"

#path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\sig_met_all\\"
path_raw<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/Figures/metabolics/sig_met_all/"
path_raw<-paste0(path_raw,disease,"/")
setwd(path_raw)
load(paste0(disease,"_sig_met_amino_lipid_others.Rdata"))
int_data<-data_lipid#------change

log2FC_matrix<-int_data[,c("LogFC (Active - Visit 3 vs Active - Visit 1)",
                           "LogFC (Placebo - Visit 3 vs Placebo - Visit 1)")]
sig_matrix<-int_data[,c("FDR (Active - Visit 3 vs Active - Visit 1)",
                        "FDR (Placebo - Visit 3 vs Placebo - Visit 1)")]

label_matrix<-matrix("",dim(sig_matrix)[1],dim(sig_matrix)[2])
index_sig<-which(sig_matrix<fdr_th)
label_matrix[index_sig]<-"*"

range(log2FC_matrix)
my.breaks <- c(seq(-2, -0.0001, by=0.01), seq(0, 2, by=0.01))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors = c("white","red"))(length(my.breaks)/2))


tiff(filename = paste0(int_super_pathway,"_heatmap_",disease,"_visit3vs1_2024.tiff"),width=7.3,height=6,units="in",res=400)
pheatmap(log2FC_matrix,
         cluster_cols = F,
         cluster_rows = T,
         clustering_method = "ward.D2",
         color =  my.colors, 
         breaks = my.breaks,
         legend_breaks = seq(-2,2,by=1),
         display_numbers = label_matrix,
         fontsize_number = 22,
         fontsize = 14,
         treeheight_row =90,
         cellwidth = 25,
         cellheight = 18,
         #gaps_col = 2,
         labels_col = c("CMA","Placebo"),
         labels_row = rownames(log2FC_matrix))
dev.off()
##########################################
#heatmap for PD---others
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
fdr_th=0.05
int_super_pathway<-"others"

#path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\sig_met_all\\"
path_raw<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/Figures/metabolics/sig_met_all/"
path_raw<-paste0(path_raw,disease,"/")
setwd(path_raw)
load(paste0(disease,"_sig_met_amino_lipid_others.Rdata"))
int_data<-data_others#------change

log2FC_matrix<-int_data[,c("LogFC (Active - Visit 3 vs Active - Visit 1)",
                           "LogFC (Placebo - Visit 3 vs Placebo - Visit 1)")]
sig_matrix<-int_data[,c("FDR (Active - Visit 3 vs Active - Visit 1)",
                        "FDR (Placebo - Visit 3 vs Placebo - Visit 1)")]

label_matrix<-matrix("",dim(sig_matrix)[1],dim(sig_matrix)[2])
index_sig<-which(sig_matrix<fdr_th)
label_matrix[index_sig]<-"*"

range(log2FC_matrix)
my.breaks <- c(seq(-2, -0.0001, by=0.01), seq(0, 2, by=0.01))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors = c("white","red"))(length(my.breaks)/2))


tiff(filename = paste0(int_super_pathway,"_heatmap_",disease,"_visit3vs1_2024.tiff"),width=8.5,height=8,units="in",res=400)
pheatmap(log2FC_matrix,
         cluster_cols = F,
         cluster_rows = T,
         clustering_method = "ward.D2",
         color =  my.colors, 
         breaks = my.breaks,
         legend_breaks = seq(-2,2,by=1),
         display_numbers = label_matrix,
         fontsize_number = 22,
         fontsize = 14,
         treeheight_row =90,
         cellwidth = 25,
         cellheight = 18,
         #gaps_col = 2,
         labels_col = c("CMA","Placebo"),
         labels_row = rownames(log2FC_matrix))
dev.off()
################################################################################
##############################heatpmap: proteomics #############################
################################################################################
#PD------protein
rm(list=ls())
rm(list=ls())
rm(list=ls())

p_th<-0.01
disease<-"PD"

path_raw<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/output_from_Arif/Results/Proteomics/"
setwd(path_raw)
int_data<-as.matrix(read.csv(paste0(disease,"_Visit.txt"),header=F,row.names = 1,sep="\t"))
colnames(int_data)<-int_data[1,]
int_data<-int_data[-1,]
mode(int_data)="numeric"

index_sig<-which(int_data[,"P-Value (Active - Visit 3 vs Active - Visit 1)"]<p_th|int_data[,"P-Value (Placebo - Visit 3 vs Placebo - Visit 1)"]<p_th)
int_data<-int_data[index_sig,]

log2FC_matrix<-int_data[,c("LogFC (Active - Visit 3 vs Active - Visit 1)",
                           "LogFC (Placebo - Visit 3 vs Placebo - Visit 1)")]

sig_matrix<-int_data[,c("P-Value (Active - Visit 3 vs Active - Visit 1)",
                        "P-Value (Placebo - Visit 3 vs Placebo - Visit 1)")]


label_matrix<-matrix("",dim(sig_matrix)[1],dim(sig_matrix)[2])
index_sig<-which(sig_matrix<p_th)
label_matrix[index_sig]<-"*"

range(log2FC_matrix)
my.breaks <- c(seq(-1, -0.0001, by=0.01), seq(0, 1, by=0.01))
my.colors <- c(colorRampPalette(colors = c("blue", "white"))(length(my.breaks)/2),
               colorRampPalette(colors = c("white","red"))(length(my.breaks)/2))

path_out<-"/Users/lixiangyu/Library/CloudStorage/OneDrive-个人/Computer/Other_project/Clinical_AD_PD/Figures/proteomics/"
path_out<-paste0(path_out,disease,"_all/")
setwd(path_out)
tiff(filename = "heatmap_protein_visit3vs1_2024.tiff",width=10,height=5,units = "in",res=400)
log2FC_matrix<-t(log2FC_matrix)
label_matrix<-t(label_matrix)
pheatmap(log2FC_matrix,
         cluster_cols = T,
         cluster_rows = F,
         clustering_method = "ward.D2",
         color =  my.colors, 
         breaks = my.breaks,
         legend_breaks = seq(-1,1,by=1),
         display_numbers = label_matrix,
         fontsize_number = 18,
         fontsize = 15,
         treeheight_row =50,
         cellwidth = 20,
         cellheight = 20,
         #gaps_col = 2,
         labels_col = colnames(log2FC_matrix),
         labels_row = c("CMA","Placebo"))
dev.off()
###############################
