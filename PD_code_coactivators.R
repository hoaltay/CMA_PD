#AD PD based line box plot---based on all samples
library(readxl)
library(VennDiagram)
library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(ggplot2)
library(ggpubr)
require(reshape2)
library(ggbeeswarm)
############################################################################################
########################box plot:clinical indicators##############################################
#MoCA-------PD-
rm(list=ls())
rm(list=ls())
rm(list=ls())

library(ggbeeswarm)

disease="PD"
int_par<-"MoCA"
yylabel="MoCA"

path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\output_from_Arif\\Data\\Clinical\\"
setwd(path_raw)
int_data_1<-as.data.frame(read_excel('20210421 ScandiBio AD-PD Both_Final.xlsx', sheet =paste0(disease,"_Alanya")))
int_data_2<-as.data.frame(read_excel('20210421 ScandiBio AD-PD Both_Final.xlsx', sheet =paste0(disease,"_Medipol")))
int_data_1<-int_data_1[,c("Patient ID","Visit NO","Group",int_par)]
int_data_2<-int_data_2[,c("Patient ID","Visit NO","Group",int_par)]

int_data<-rbind(int_data_1,int_data_2)
colnames(int_data)[4]="value"
mode(int_data$value)="numeric"
int_data$type<-paste0(int_data$Group,"_",int_data$`Visit NO`)


path_out<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\clinical_data\\box_plot_all\\"
path_out<-paste0(path_out,disease,"\\")
setwd(path_out)


int_data$type <- factor(int_data$type,levels = c("Active_Visit 1","Active_Visit 2","Active_Visit 3","Placebo_Visit 1","Placebo_Visit 2","Placebo_Visit 3"))
mycolor<-c("#636363","#6a3d9a","#e31a1c","#bdbdbd","#cab2d6","#fb9a99")

sig_file<- data.frame(read_excel('MoCA_stat_visit123.xlsx', sheet = 1))#-----------change
y_pos<-as.matrix(c(30,34,30,34))#----------change
sig_file$y_pos=y_pos

ggplot(int_data,aes(x=type,y=value))+
  geom_boxplot(fill=mycolor,outlier.size = -1,lwd=0.5)+
  geom_quasirandom(cex=0.8,color="black")+
  scale_colour_manual(values=mycolor)+
  scale_y_continuous(limits = c(5,35))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="",y=yylabel,size=14)+
  font("xy.text", size = 12, color = "black",angle=90)+
  font("ylab", size = 20)+
  stat_pvalue_manual(sig_file, label = "p_label",xmin="group1",xmax="group2",y.position = "y_pos",tip.length = 0.02,bracket.size = 0.7,vjust=-0.3,label.size = 4)


ggsave(paste0(int_par,"_box_beeswarm_visit123_CV2to3",".pdf"),width=4,height=5,dpi=300)#----change
####################################################################################################
###################################################box plot: four cofactors#########################
#for PD
#serine ------PD
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
int_met<-"serine"

path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\output_from_Arif\\Data\\Metabolomics\\"
setwd(path_raw)
int_data<-as.matrix(read.csv(paste0(disease,"_ALL_mets_data.txt"),header=F,sep="\t"))
colnames(int_data)<-int_data[1,]
int_data<-int_data[-1,]

int_data<-int_data[,c("sampleID",int_met)]
colnames(int_data)[2]="value"
int_data<-as.data.frame(int_data)
mode(int_data$value)="numeric"

split_result<-strsplit(int_data$sampleID,"_")
visit<-as.matrix(unlist(lapply(split_result,function(x) x[[2]])))
group<-as.matrix(unlist(lapply(split_result,function(x) x[[3]])))
int_data$type<-as.matrix(paste0(visit,"_",group))

path_out<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\four_cofactors\\"
path_out<-paste0(path_out,disease,"\\")
setwd(path_out)

int_data$type <- factor(int_data$type,levels = c("Visit 1_Active","Visit 3_Active","Visit 1_Placebo","Visit 3_Placebo"))
mycolor<-c("#636363","#e31a1c","#bdbdbd","#fb9a99")

sig_file<- data.frame(read_excel(paste0(int_met,'_stat_fdr.xlsx'), sheet = 1))#-----------change
y_pos<-as.matrix(c(2.1,1.1,2.5,2.9))#----------change
sig_file$y_pos=y_pos

ggplot(int_data,aes(x=type,y=value))+
  geom_boxplot(fill=mycolor,outlier.size =-1 ,lwd=0.5)+
  geom_quasirandom(cex=0.8,color="black")+
  scale_colour_manual(values=mycolor)+
  scale_y_continuous(limits = c(0,3))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="",y=int_met,size=14)+
  font("xy.text", size = 12, color = "black",angle=90)+
  font("ylab", size = 20)+
  stat_pvalue_manual(sig_file, label = "fdr",xmin="group1",xmax="group2",y.position = "y_pos",tip.length = 0.018,bracket.size = 0.7,vjust=-0.3,label.size = 4)

ggsave(paste0(int_met,"_box_beeswarm_fdr",".pdf"),width=4,height=5,dpi=300)#----change
####################################################################
#carnitine ------PD
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
int_met<-"carnitine"

path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\output_from_Arif\\Data\\Metabolomics\\"
setwd(path_raw)
int_data<-as.matrix(read.csv(paste0(disease,"_ALL_mets_data.txt"),header=F,sep="\t"))
colnames(int_data)<-int_data[1,]
int_data<-int_data[-1,]

int_data<-int_data[,c("sampleID",int_met)]
colnames(int_data)[2]="value"
int_data<-as.data.frame(int_data)
mode(int_data$value)="numeric"

split_result<-strsplit(int_data$sampleID,"_")
visit<-as.matrix(unlist(lapply(split_result,function(x) x[[2]])))
group<-as.matrix(unlist(lapply(split_result,function(x) x[[3]])))
int_data$type<-as.matrix(paste0(visit,"_",group))

path_out<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\four_cofactors\\"
path_out<-paste0(path_out,disease,"\\")
setwd(path_out)

int_data$type <- factor(int_data$type,levels = c("Visit 1_Active","Visit 3_Active","Visit 1_Placebo","Visit 3_Placebo"))
mycolor<-c("#636363","#e31a1c","#bdbdbd","#fb9a99")

sig_file<- data.frame(read_excel(paste0(int_met,'_stat_fdr.xlsx'), sheet = 1))#-----------change
y_pos<-as.matrix(c(1.51,1.3,1.63,1.76))#----------change
sig_file$y_pos=y_pos

ggplot(int_data,aes(x=type,y=value))+
  geom_boxplot(fill=mycolor,outlier.size =-1 ,lwd=0.5)+
  geom_quasirandom(cex=0.8,color="black")+
  scale_colour_manual(values=mycolor)+
  scale_y_continuous(limits = c(0.6,1.8))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="",y=int_met,size=14)+
  font("xy.text", size = 12, color = "black",angle=90)+
  font("ylab", size = 20)+
  stat_pvalue_manual(sig_file, label = "fdr",xmin="group1",xmax="group2",y.position = "y_pos",tip.length = 0.018,bracket.size = 0.7,vjust=-0.3,label.size = 4)

ggsave(paste0(int_met,"_box_beeswarm_fdr",".pdf"),width=4,height=5,dpi=300)#----change
####################################################################
#cysteine ------PD
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
int_met<-"cysteine"

path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\output_from_Arif\\Data\\Metabolomics\\"
setwd(path_raw)
int_data<-as.matrix(read.csv(paste0(disease,"_ALL_mets_data.txt"),header=F,sep="\t"))
colnames(int_data)<-int_data[1,]
int_data<-int_data[-1,]

int_data<-int_data[,c("sampleID",int_met)]
colnames(int_data)[2]="value"
int_data<-as.data.frame(int_data)
mode(int_data$value)="numeric"

split_result<-strsplit(int_data$sampleID,"_")
visit<-as.matrix(unlist(lapply(split_result,function(x) x[[2]])))
group<-as.matrix(unlist(lapply(split_result,function(x) x[[3]])))
int_data$type<-as.matrix(paste0(visit,"_",group))

path_out<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\four_cofactors\\"
path_out<-paste0(path_out,disease,"\\")
setwd(path_out)

int_data$type <- factor(int_data$type,levels = c("Visit 1_Active","Visit 3_Active","Visit 1_Placebo","Visit 3_Placebo"))
mycolor<-c("#636363","#e31a1c","#bdbdbd","#fb9a99")

sig_file<- data.frame(read_excel(paste0(int_met,'_stat_fdr.xlsx'), sheet = 1))#-----------change
y_pos<-as.matrix(c(1.55,1.55,1.75,1.95))#----------change
sig_file$y_pos=y_pos

ggplot(int_data,aes(x=type,y=value))+
  geom_boxplot(fill=mycolor,outlier.size =-1 ,lwd=0.5)+
  geom_quasirandom(cex=0.8,color="black")+
  scale_colour_manual(values=mycolor)+
  #scale_y_continuous(limits = c(0.6,1.8))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="",y=int_met,size=14)+
  font("xy.text", size = 12, color = "black",angle=90)+
  font("ylab", size = 20)+
  stat_pvalue_manual(sig_file, label = "fdr",xmin="group1",xmax="group2",y.position = "y_pos",tip.length = 0.018,bracket.size = 0.7,vjust=-0.3,label.size = 4)

ggsave(paste0(int_met,"_box_beeswarm_fdr",".pdf"),width=4,height=5,dpi=300)#----change
####################################################################
#nicotinamide ------PD
rm(list=ls())
rm(list=ls())
rm(list=ls())

disease="PD"
int_met<-"nicotinamide"

path_raw<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\output_from_Arif\\Data\\Metabolomics\\"
setwd(path_raw)
int_data<-as.matrix(read.csv(paste0(disease,"_ALL_mets_data.txt"),header=F,sep="\t"))
colnames(int_data)<-int_data[1,]
int_data<-int_data[-1,]

int_data<-int_data[,c("sampleID",int_met)]
colnames(int_data)[2]="value"
int_data<-as.data.frame(int_data)
mode(int_data$value)="numeric"

split_result<-strsplit(int_data$sampleID,"_")
visit<-as.matrix(unlist(lapply(split_result,function(x) x[[2]])))
group<-as.matrix(unlist(lapply(split_result,function(x) x[[3]])))
int_data$type<-as.matrix(paste0(visit,"_",group))

path_out<-"C:\\Users\\xiangyu.li\\OneDrive\\Computer\\Other_project\\Clinical_AD_PD\\Figures\\metabolics\\four_cofactors\\"
path_out<-paste0(path_out,disease,"\\")
setwd(path_out)

int_data$type <- factor(int_data$type,levels = c("Visit 1_Active","Visit 3_Active","Visit 1_Placebo","Visit 3_Placebo"))
mycolor<-c("#636363","#e31a1c","#bdbdbd","#fb9a99")

sig_file<- data.frame(read_excel(paste0(int_met,'_stat_fdr.xlsx'), sheet = 1))#-----------change
y_pos<-as.matrix(c(1.3,1,1.55,1.8))#----------change
sig_file$y_pos=y_pos

ggplot(int_data,aes(x=type,y=value))+
  geom_boxplot(fill=mycolor,outlier.size =-1 ,lwd=0.5)+
  geom_quasirandom(cex=0.8,color="black")+
  scale_colour_manual(values=mycolor)+
  #scale_y_continuous(limits = c(0,1.5))+
  theme_classic()+
  theme(legend.position = "none")+
  labs(x="",y=int_met,size=14)+
  font("xy.text", size = 12, color = "black",angle=90)+
  font("ylab", size = 20)+
  stat_pvalue_manual(sig_file, label = "fdr",xmin="group1",xmax="group2",y.position = "y_pos",tip.length = 0.018,bracket.size = 0.7,vjust=-0.3,label.size = 4)

ggsave(paste0(int_met,"_box_beeswarm_fdr",".pdf"),width=4,height=5,dpi=300)#----change

