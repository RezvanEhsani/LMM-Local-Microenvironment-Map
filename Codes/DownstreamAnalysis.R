
#' This is an example to clustering patients from keren et al. paper:
#'                https://www.cell.com/cell/pdf/S0092-8674(18)31100-0.pdf
#' Here we just use Macrophage cell type as cells of interest.


source("C:/Users/reh021/Downloads/Github/ExtractFeatures.R")
source("C:/Users/reh021/Downloads/Github/ClusteringCells.R")
source("C:/Users/reh021/Downloads/Github/ClusteringSamples.R")

##
# libraries
library(data.table)
library(Rfast)
library(dplyr)
library(caramellar)
library(progressr)
library(parallel)
library(tidyverse)
library(pbmcapply)
library(Rtsne)
#
library(ConsensusClusterPlus)
library(pals)
library(pheatmap)
library(ggalluvial)
library(ggpubr)
library(ggrepel) 
library(RColorBrewer)
library(NLP)
#
library(ggplot2)
library(ggfortify)
library(survival)
library(survminer)
library("readxl")
###############
# load keren data
KerenData <- read.csv("C:/Users/reh021/Downloads/Github/KerenData.csv",
                      header=T,
                      sep=",",
                      stringsAsFactors = FALSE)

# Setting parameters:
setwd("C:/Users/reh021/Downloads/Github/Macrophages")
Path = "C:/Users/reh021/Downloads/Github/Macrophages/"
#
Data = KerenData
SampleID_col = 2
CellID_col = 1
CellType_col = 63
CellTypesOfInterest = "Macrophages"
MarkersOfInterest_col = NULL
Zscore = TRUE
PositivityCutoff = 0.5
MetaCellType_col = 53
X_col = 58
Y_col = 59
r = 100
minCount = 20
ExistingClass_col = 64
#
NCellTypes <- 11


##########################################################
#                                                        #
#               Feature extraction                       #
#                                                        #
##########################################################

# Load feature data if you save it before, if not call Extract_Features function
KerenFeatures <- read.csv("C:/Users/reh021/Downloads/Github/Macrophages/MacrophagesFeatures.csv",
                        header=T,
                        sep=",",
                        stringsAsFactors = FALSE)
# Extract Features
start_time <- Sys.time()
KerenFeatures <- Extract_Features(

    Data  = Data,
    SampleID_col = SampleID_col,
    CellID_col = CellID_col,
    CellType_col = CellType_col,
    CellTypesOfInterest = CellTypesOfInterest,
    MarkersOfInterest_col = MarkersOfInterest_col,
    Zscore = Zscore,
    PositivityCutoff = PositivityCutoff,
    MetaCellType_col = MetaCellType_col,
    X_col = X_col,
    Y_col = Y_col,
    r = r,
    minCount = minCount,
    ExistingClass_col = ExistingClass_col
    )
end_time <- Sys.time()
print(end_time - start_time)

# Save Feature dataframe
write.table(
    
    KerenFeatures,
    "C:/Users/reh021/Downloads/Github/Macrophages/MacrophagesFeatures.csv",
    append = FALSE,
    sep = ",",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE,
    quote=FALSE
    )


#Filter "cold" samples, see the keren paper
SamplesKeren <- c(35,28,16,37,40,4,41,36,3,5,34,32,6,9,10,13,39,29,17,23,1,
                  33,12,27,8,2,38,20,7,14,11,21,31,18)
KerenFeatures <- KerenFeatures[(KerenFeatures$sample_id %in% SamplesKeren),]

##########################################################
#                                                        #
#               clustering in Cell level                 #
#                                                        #
##########################################################

#optimal k for clustering the Macrophage cells
Nclusters_CellLevel(KerenFeatures,
                    NCellTypes = NCellTypes,
                    Nseed = 1234,
                    PlotName="Keren_Features_Macrophages")

# Clustering Macrophage cells
ResCellLevel <- ClusteringCellLevel(Dat_Feature = KerenFeatures,
                                    NCellTypes = 11,
                                    k = 8
                                    )

Dat_Cells = ResCellLevel$CellClusterData

# save plots
pdf(paste(Path,"CellHeatmap.pdf",sep=""))
print(ResCellLevel$CellHeatmap)
dev.off()
#
pdf(paste(Path,"CellTsne.pdf",sep=""))
print(ResCellLevel$CellTsne)
dev.off()
#
pdf(paste(Path,"CellBoxplot.pdf",sep=""))
print(ResCellLevel$CellBoxplot)
dev.off()
#
pdf(paste(Path,"CellBarplot.pdf",sep=""))
print(ResCellLevel$CellBarplot)
dev.off()


##########################################################
#                                                        #
#               clustering in sample level               #
#                                                        #
##########################################################

Nclusters_SampleLevel(Dat_Cells,
                      maxK = 5,
                      PlotName = "Keren_Macrophages_Sample")


ResSampleLevel <- ClusteringSampleLevel(Dat_Cells,
                                        NCellTypes = NCellTypes,
                                        k_sample = 4,
                                        Nseed = 1234,
                                        ExistClasses = TRUE
                                        )
#
SampleClusterData <- ResSampleLevel$SampleClusterData
# save plots
pdf(paste(Path,"StackBar.pdf",sep=""))
print(ResSampleLevel$StackBar)
dev.off()
#
pdf(paste(Path,"SampleHeatmap.pdf",sep=""))
print(ResSampleLevel$SampleHeatmap)
dev.off()
#
pdf(paste(Path,"SampleAlluvium.pdf",sep=""))
print(ResSampleLevel$SampleAlluvium)
dev.off()
#
pdf(paste(Path,"SampleTsne.pdf",sep=""))
print(ResSampleLevel$SampleTsne)
dev.off()



##########################################################
#                                                        #
#               Extra downstream analysis                #
#                                                        #
##########################################################

library(patchwork)
library(aplot)

# 1) immunregulatory proteins analysis from keren paper figure5. 

keren_immunregulatory_proteins <- function(Samples, Scores, Titel, Dat_Cells){

        df <- data.frame(Samples, Scores)
        colnames(df) <- c("sample_id", "Keren_Score")
        df <- df[(df$sample_id %in% unique(Dat_Cells$sample_id)),]
        df$Class[df$sample_id %in% c(3, 4, 5, 6, 9, 10, 16, 28, 32, 35, 36, 37, 40,34, 41)] <- 'Compartmentalized'
        df$Class[df$sample_id %in% c(13,39,29,17,23,1,33,12,27,8,2,38,20,7,14,11,21,31,18)] <- 'Mixed'

        df$sample_id <- factor(df$sample_id,levels=df$sample_id)
    
        sample_classes <- SampleClusterData
        sample_classes <- sample_classes[(sample_classes$sample_id %in% unique(df$sample_id)), ]
        sample_classes <- merge(df, sample_classes, by="sample_id")
        sample_classes <- sample_classes[order(sample_classes$SampleCluster, decreasing = TRUE), ]
        
        sample_classes$sample_id <- factor(sample_classes$sample_id,levels=sample_classes$sample_id)

        Keren_plot <- ggplot(data=sample_classes, aes(x=sample_id, y=Keren_Score, fill=SampleCluster)) +
                      geom_bar(stat="identity")+
                      scale_fill_manual(values=Sample_Colors)+
                      labs(title="", y=Titel, x="")+
                      theme(
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                            axis.title = element_text(size = 10),
                            panel.background = element_rect(fill = "white")
                            )+
                      coord_flip()

        Dat_Sample <- dcast(Dat_Cells, sample_id  ~ Clusters)
        Dat_Sample <- Dat_Sample[(Dat_Sample$sample_id %in% unique(sample_classes$sample_id)),]
        L <- length(unique(Dat_Cells$Clusters))+1

        Dat_Sample[,2:L] <- Dat_Sample[,2:L]/rowSums(Dat_Sample[,2:L])

        Dat_Sample$sample_id <- factor(Dat_Sample$sample_id,unique(Dat_Sample$sample_id))

        Dat_Sample$sample_id <- factor(Dat_Sample$sample_id,levels=sample_classes$sample_id)
        melted <- melt(Dat_Sample[1:L], id = c("sample_id"))

        # Stackbar
        colnames(melted)[2] <- "CellCluster"
        pa <- ggplot(data=melted, aes(x=sample_id, y=value,fill=CellCluster)) +
                 geom_bar(stat="identity", width=0.83)+
                 labs(title="",y="", x="")+
                 scale_fill_manual(values=Cell_Colors)+
                theme(axis.title.x = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      panel.background = element_rect(fill = "white")
                      )+
          coord_flip()

     
        g <- ggplot(sample_classes, aes(sample_id,y=1,fill=Class)) +
                geom_tile(colour="white")+
                scale_fill_manual(values=c("#990099","#009900"))+
                labs(title="", y="", x="Samples")+
                coord_equal()+
                theme(
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks = element_blank(),
                      panel.background = element_rect(fill = "white")
                      )+
                  coord_flip()

    cowplot::plot_grid(g, Keren_plot, pa,  nrow=3)
    plot_all_3 <- pa %>% insert_left(Keren_plot) %>% insert_left(g, width=0.08) 

    return(plot_all_3)
}



Nfeatues <- 6 + 2*NCellTypes + 12
Dat_Cells1 <- Dat_Cells[-c(2:Nfeatues)]

###########
# For PDL1+
Samples <- c(8,33,27,31,18,23,38,21,14,29,36,20,3,17,2,11,13,40,12,28,39,16,41,4,10,32,5,35,9,37,6)
Scores <- c(7,6,6,5,4,3.7,3.4,3,2.5,2.3,2.1,2.1,1.4,1,1,0.8,0.4,-0.1,-0.1,-1.4,-1.5,-1.6,-1.9,-2,-2.1,
             -2.2,-3.7,-4,-5.8,-6,-6)

Titel="Ratio of PDL1+ cells \n [log2(Tumor/Immune)]"
PPDL1 <- keren_immunregulatory_proteins(Samples, Scores, Titel, Dat_Cells1)
PPDL1
pdf(paste(Path,"Plot_PDL1.pdf",sep=""))
print(PPDL1)
dev.off()

#######
#PD1+
Titel="Ratio of PD1+ cells \n [log2(CD8+/CD4+)]"
Samples <- c(14,20,2,29,13,11,18,10,32,37,39,16,9,40,12,5,27,38,6,4,17,3,35,28)
Scores <- c(3.5,2.4,2.3,1.8,1.6,1.5,1.3,0.6,0.5,0.2,0.2,-0.2,-0.3,-0.5,-0.7,-0.7,-0.8,
            -1,-1.1,-1.8,-1.8,-2.7,-3.2,-3.5)

PPD1 <- keren_immunregulatory_proteins(Samples, Scores, Titel, Dat_Cells1)
PPD1

pdf(paste(Path,"Plot_PD1.pdf",sep=""))
print(PPD1)
dev.off()
#######
#IDO+
Titel="Ratio of IDO+ cells \n [log2(Tumor/Immune)]"
Samples <- c(20,14,39,12,13,32,11,40,17,10,2,37,6,9,29,5,41,3,16,4,28,35)
Scores <- c(2,1.7,0.5,0.1,0.1,-0.3,-1,-1.7,-2.4,-2.8,-3,-3.1,-3.3,-3.3,-3.4,-3.7,-4,-4,-5.2,-6.8,-7,-7)

#
PIDO <- keren_immunregulatory_proteins(Samples, Scores, Titel, Dat_Cells1)
pdf(paste(Path,"Plot_IDO.pdf",sep=""))
print(PIDO)
dev.off()

# 2) mixing Score figure S5. B.
MIXINGPLOT <- function(Samples, Scores, Titel, Dat_Cells){

        df <- data.frame(Samples, Scores)
        colnames(df) <- c("sample_id", "Keren_Score")
        df <- df[(df$sample_id %in% unique(Dat_Cells$sample_id)),]
        df$Class[df$sample_id %in% c(3, 4, 5, 6, 9, 10, 16, 28, 32, 35, 36, 37, 40,34, 41)] <- 'Compartmentalized'
        df$Class[df$sample_id %in% c(13,39,29,17,23,1,33,12,27,8,2,38,20,7,14,11,21,31,18)] <- 'Mixed'

        df$sample_id <- factor(df$sample_id,levels=df$sample_id)
        sample_classes <- SampleClusterData
        sample_classes <- sample_classes[(sample_classes$sample_id %in% unique(df$sample_id)), ]
        sample_classes <- merge(df, sample_classes, by="sample_id")
        sample_classes <- sample_classes[order(sample_classes$Keren_Score, decreasing = TRUE), ]
        
        sample_classes$sample_id <- factor(sample_classes$sample_id,levels=sample_classes$sample_id)

        Keren_plot <- ggplot(data=sample_classes, aes(x=sample_id, y=Keren_Score, fill=SampleCluster)) +
                      geom_bar(stat="identity")+
                      scale_fill_manual(values=Sample_Colors)+
                      labs(title="", y=Titel, x="")+
                      theme(axis.title.y = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.title = element_text(size = 10),
                            panel.background = element_rect(fill = "white")
                            )+
                      coord_flip()

        Dat_Sample <- dcast(Dat_Cells, sample_id  ~ Clusters)
        Dat_Sample <- Dat_Sample[(Dat_Sample$sample_id %in% unique(df$sample_id)),]
        L <- length(unique(Dat_Cells$Clusters))+1
        Dat_Sample[,2:L] <- Dat_Sample[,2:L]/rowSums(Dat_Sample[,2:L])
        Dat_Sample$sample_id <- factor(Dat_Sample$sample_id,unique(Dat_Sample$sample_id))

        Dat_Sample$sample_id <- factor(Dat_Sample$sample_id,levels=sample_classes$sample_id)
        melted <- melt(Dat_Sample, id = c("sample_id"))


        # Stackbar
        colnames(melted)[2] <- "CellCluster"
        pa <- ggplot(data=melted, aes(x=sample_id, y=value,fill=CellCluster)) +
                 geom_bar(stat="identity", width=0.83)+
                 labs(title="",y="", x="")+
                 scale_fill_manual(values=Cell_Colors)+
                theme(axis.title.x = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      panel.background = element_rect(fill = "white")
                      )+
          coord_flip()

        
        g <- ggplot(sample_classes, aes(sample_id,y=1,fill=MixingClass)) +
                geom_tile(colour="white")+
                scale_fill_manual(values=c("#990099","#009900"))+
                labs(title="", y="", x="Samples")+
                coord_equal()+
                theme(
                      axis.title.x = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks = element_blank(),
                      panel.background = element_rect(fill = "white")
                      )+
                  coord_flip()

    cowplot::plot_grid(g, Keren_plot, pa,  nrow=3)
    plot_all_3 <- pa %>% insert_left(Keren_plot) %>% insert_left(g, width=0.08) 

    return(plot_all_3)
}

Titel="Mixing Score"
Scores <- c(.01,.03,.05,.06,.08,.08,.1,.13,.15,.16,.16,.2,.2,.21,.21,.23,.28,.29,.32,.37,.39
                  ,.47,.50,.5,.505,.595,.605,.605,.62,.66,.67,.72,.76,.87)
Samples <- c(35,28,16,37,40,4,41,36,3,5,34,32,6,9,10,13,39,29,17,23,1,
               33,12,27,8,2,38,20,7,14,11,21,31,18)
PMixing <- MIXINGPLOT(Samples, Scores, Titel, Dat_Cells1)

pdf(paste(Path,"Plot_Mixing.pdf",sep=""))
print(PMixing)
dev.off()



# 3) survival analysis figure7 B.
Surv <- read_excel("C:/Users/reh021/Downloads/Github/mibitof_cohort_metadata.xlsx")
Surv <- data.frame(Surv)
Surv <- Surv[2:40,c(1,28,29)]

colnames(Surv) <- c("sample_id", "Survival_days", "Censored")
Surv <- Surv[(Surv$sample_id %in% c(35,28,16,37,40,4,41,36,3,5,34,32,6,9,10,
                                   13,39,29,17,23,1,33,12,27,8,2,38,20,7,14,11,21,31,18)),]

L1 <- length(unique(Dat_Cells$Clusters))
DatSample <- SampleClusterData[,c(1,(L1+2),(L1+3))]

Surv_Dat <- merge(Surv, DatSample, by="sample_id")
Surv_Dat$Survival_days <- as.numeric(Surv_Dat$Survival_days)/365
Surv_Dat$Censored <- as.numeric(Surv_Dat$Censored)
Surv_Dat$Censored1[Surv_Dat$Censored==0] <- 1
Surv_Dat$Censored1[Surv_Dat$Censored==1] <- 0

colnames(Surv_Dat)[4] <- "SC"
Surv_Dat$col[Surv_Dat$SC=="SC1"] <- Sample_Colors[1]
Surv_Dat$col[Surv_Dat$SC=="SC2"] <- Sample_Colors[2]
Surv_Dat$col[Surv_Dat$SC=="SC3"] <- Sample_Colors[3]
Surv_Dat$col[Surv_Dat$SC=="SC4"] <- Sample_Colors[4]
#Surv_Dat$col[Surv_Dat$SC=="SC5"] <- Sample_Colors[5]

model_fit <- survfit(Surv(Survival_days, Censored1) ~ SC, data = Surv_Dat)
surv1 <- ggsurvplot(model_fit,  size = 1,
           break.time.by = 2, 
           palette = Sample_Colors[1:4],
           risk.table = TRUE, risk.table.y.text.col = TRUE)

pdf(paste(Path,"SurvivalPlot.pdf",sep=""))
print(surv1)
dev.off()

