library(tidyverse)
library(dplyr)
library(edgeR)
library(gplots)
library(topGO)
library(WGCNA)
library(dendextend)
library(rlist)
library(plyr)
library(gridExtra)
library(ggExtra)
library(data.table)
library(viridis)
library(scales)
library(bioDist)
library(cowplot)
library(RColorBrewer)
library(rlang)
library(ggpubr)
library(ggforce)

#########################################
#Read in data as a matrix (output of featureCounts)
dat <- read_csv("Poplar_counts.csv") %>% 
  mutate(Geneid = str_sub(Geneid, 1, -6)) %>% 
  dplyr::arrange(Geneid) 
countdata <- dat %>% 
  dplyr::select(-c(1:6)) %>% 
  data.matrix() 
colnames(countdata) <- c("LM30A","LM80C","LM30B","LM30C","LM50A","LM50B","LM50C","LM80A","LM80B",
                         "EE80B","EE30A","EE80C","EE30B","EE30C","EE50A","EE50B","EE50C","EE80A")
rownames(countdata) <- dat$Geneid

#Group the samples and filter based on average log CPM
group <- factor(c("LM30","LM80","LM30","LM30","LM50","LM50","LM50","LM80","LM80",
                  "EE80","EE30","EE80","EE30","EE30","EE50","EE50","EE50","EE80"),
                levels = c("LM30", "LM50", "LM80", "EE30", "EE50", "EE80"))
x <- DGEList(countdata, group = group)
y <- x[aveLogCPM(x, normalized.lib.sizes = F) > -2 & aveLogCPM(x, normalized.lib.sizes = F) < 9, , 
       keep.lib.sizes=FALSE]

#Normalize and estimate dispersion 
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- calcNormFactors(y, method = "TMM")
y <- estimateGLMRobustDisp(y, design, record = T)

#######################################
#DGE analysis

#DE test (p-value < 0.05, lfc = 0) for each pair
#the results are combined into a single dataframe
#1 = upregulated, -1 = downregualted
DGE <- data.frame(decideTestsDGE(exactTest(y, pair = c(1,2))),
                  decideTestsDGE(exactTest(y, pair = c(1,3))),
                  decideTestsDGE(exactTest(y, pair = c(2,3))),
                  decideTestsDGE(exactTest(y, pair = c(4,5))),
                  decideTestsDGE(exactTest(y, pair = c(4,6))),
                  decideTestsDGE(exactTest(y, pair = c(5,6))),
                  decideTestsDGE(exactTest(y, pair = c(1,4))),
                  decideTestsDGE(exactTest(y, pair = c(2,5))),
                  decideTestsDGE(exactTest(y, pair = c(3,6))))

#turning this into a function so I don't have to write it out 12 times
#takes the rownames (genes) that have a value of 1 or -1 in column x
sel <- function(x,y){
  rownames(filter(DGE, DGE[,x] == y))
}

#list of genes in each group; split into up and down
DGEgeneList <- list(LM30.50up = sel(1,1), LM30.50down = sel(1,-1),
                    LM30.80up = sel(2,1), LM30.80down = sel(2,-1),
                    LM50.80up = sel(3,1), LM50.80down = sel(3,-1),
                    EE30.50up = sel(4,1), EE30.50down = sel(4,-1),
                    EE30.80up = sel(5,1), EE30.80down = sel(5,-1),
                    EE50.80up = sel(6,1), EE50.80down = sel(6,-1),
                    LM30.EE30up = sel(7,1), LM30.EE30down = sel(7,-1),
                    LM50.EE50up = sel(8,1), LM50.EE50down = sel(8,-1),
                    LM80.EE80up = sel(9,1), LM80.EE80down = sel(9,-1))

#write as csv (don't have to run this every time)
list.save(DGEgeneList, "DGEgeneList.rdata")

#venn that shows both up and down genes
#use DGE[,1:3] for LM, [4:6] for EE, [7:9] for time of day
vennDiagram(DGE[,1:3],
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))
dev.off()

#######################################
#load DGE data
DGEgeneList <- list.load("DGEgeneList.rdata")

#Convert to logCPM
z <- cpm(y, prior.count=1, log=T)

#or load csv so you don't have to run all the previous stuff
z <- as.matrix(read.csv("z.csv", row.names = "X"))
  
#Filtering (if only LM or EE is wanted)
#need to refilter aveLogCPM due to time-of-day differences
#run either LM or EE to set all downstream analyses to one

#subset LM samples 
z <- as.data.frame(z[,1:9]) %>% 
  mutate(rowM = rowMeans(z)) %>% 
  filter(rowM > -2 & rowM < 9) %>% 
  dplyr::select(-rowM) %>% 
  as.matrix()
time <- "LM"

#subset EE samples
z <- as.data.frame(z[,10:18]) %>% 
  mutate(rowM = rowMeans(z)) %>% 
  filter(rowM > -2 & rowM < 9) %>% 
  dplyr::select(-rowM) %>% 
  as.matrix()
time <- "EE"

#######################################

#stacked bar plot
#code modified from https://static-content.springer.com/esm/art%3A10.1038%2Fs41477-021-00874-5/MediaObjects/41477_2021_874_MOESM1_ESM.pdf

#load list of differentially expressed genes
DGEgeneList <- list.load("DGEgeneList.rdata")

#alternative DGE list for any pval
#will be used to rank top 1000 genes of each category
DGE <- data.frame(rep(NA, 1000))
DGEgeneList <- list()
for(i in 1:6){
  temp <- data.frame(exactTest(y, pair = data.frame(c(1,2),c(1,3),c(2,3),c(4,5),c(4,6),c(5,6))[,i]))
  tempUp <- temp %>% 
    filter(logFC > 0) %>% 
    dplyr::arrange(PValue)
  tempDown <- temp %>% 
    filter(logFC < 0) %>% 
    dplyr::arrange(PValue)
  DGE <- cbind(DGE, rownames(tempUp[1:1000,]), rownames(tempDown[1:1000,]))
}

#list of genes in each group; split into up and down
DGEgeneList <- list(LM30.50up = DGE[,2], LM30.50down = DGE[,3],
                    LM30.80up = DGE[,4], LM30.80down = DGE[,5],
                    LM50.80up = DGE[,6], LM50.80down = DGE[,7],
                    EE30.50up = DGE[,8], EE30.50down = DGE[,9],
                    EE30.80up = DGE[,10], EE30.80down = DGE[,11],
                    EE50.80up = DGE[,12], EE50.80down = DGE[,13])

####################################
#function to make list of up or down regulated genes
SingleTimeTable <- function(direction){
  LM30.50 <- unlist(DGEgeneList[paste0("LM30.50", direction)])
  LM30.80 <- unlist(DGEgeneList[paste0("LM30.80", direction)])
  LM50.80 <- unlist(DGEgeneList[paste0("LM50.80", direction)])
  EE30.50 <- unlist(DGEgeneList[paste0("EE30.50", direction)])
  EE30.80 <- unlist(DGEgeneList[paste0("EE30.80", direction)])
  EE50.80 <- unlist(DGEgeneList[paste0("EE50.80", direction)])
  
  totalgenes <- c(LM30.50,LM30.80,LM50.80,EE30.50,EE30.80,EE50.80)
  
  timetable <- data.frame(ATG=totalgenes, LM30.50=0,LM30.80=0,LM50.80=0,EE30.50=0,EE30.80=0,EE50.80=0) 
  timetable$LM30.50[timetable$ATG %in% LM30.50] <- 1
  timetable$LM30.80[timetable$ATG %in% LM30.80] <- 1
  timetable$LM50.80[timetable$ATG %in% LM50.80] <- 1
  timetable$EE30.50[timetable$ATG %in% EE30.50] <- 1 
  timetable$EE30.80[timetable$ATG %in% EE30.80] <- 1 
  timetable$EE50.80[timetable$ATG %in% EE50.80] <- 1 
  
  timetable <- timetable[!duplicated(timetable),]
  rownames(timetable) <- timetable$ATG
  timetable <- timetable[-1]
  return(timetable) }

#make the up/down tables
up <- SingleTimeTable("up")
down <- SingleTimeTable("down")

#create dataframe with all data needed for plot
Direction <- c(rep("up",36), rep("down",36))
Comparison <- rep(c(rep("LM30.50",6), rep("LM30.80",6), rep("LM50.80",6), rep("EE30.50",6), rep("EE30.80",6), rep("EE50.80",6)),2)
Degree <- rep(seq(1:6),12)
Cardinality <- rep(0,72)
bar <- data.frame(Direction, Comparison, Degree, Cardinality)

#Filling in that table
comps <- c("LM30.50","LM30.80","LM50.80","EE30.50","EE30.80","EE50.80")[2]
FillInBar <- function(timetable, direction, pass){ 
  bar <- pass

  for(h in 1:6){
    temp <- timetable[timetable[,h]==1,]
    for(i in 1:7){
      bar$Cardinality[bar$Comparison== c("LM30.50","LM30.80","LM50.80","EE30.50","EE30.80","EE50.80")[h] & bar$Degree==i & bar$Direction== direction] <- length(rowSums(temp)[rowSums(temp)==i])
    }
  }

  return(bar)
}

bar <- FillInBar(up, "up", bar)
bar <- FillInBar(down, "down", bar)

#Make the downregulated DEG bars extend downwards
bar$Cardinality[bar$Direction=="down"] <- -bar$Cardinality[bar$Direction=="down"]

#Some parameter wrangling for better graphing
bar$Degree <- as.factor(bar$Degree)
bar$Comparison <- factor(bar$Comparison, levels = c("LM30.80","LM30.50","LM50.80","EE30.80","EE30.50","EE50.80"))
definedvibrant <- c("0"="white", "LM30.50"="#33BBEE", "LM30.80"="#0077BB", "LM50.80"="#009988", 
                    "EE30.50"="#EE7733", "EE30.80"="#DDCC77", "EE50.80"="#CC3311")

#actual plot
barplot <- ggplot(bar, aes(x=Comparison, y=Cardinality, group=Degree)) +
  geom_bar(aes(fill=Comparison), stat="identity", position=position_stack()) +
  scale_fill_manual(values=definedvibrant, guide="none") +
  geom_bar(aes(alpha=Degree), fill="black", stat="identity", position=position_stack()) +
  geom_bar(aes(alpha=as.factor(7-as.numeric(Degree))), fill="white", stat="identity", position=position_stack())  +
  scale_alpha_manual(values=c(0,0.05,0.1,0.17,0.27,0.35))+
  scale_x_discrete(labels = c("LM30v80","LM30v50","LM50v80","EE30v80","EE30v50","EE50v80")) +
  guides(alpha=guide_legend(ncol=6)) +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(hjust=0.5, size=8), 
        panel.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        panel.border = element_blank(), 
        legend.box=NULL, legend.key.size=unit(0.4, "cm"), 
        legend.position ="top", 
        legend.text=element_text(size=10), 
        legend.title = element_text(size=10), 
        axis.title=element_text(size=10)) +
  labs(y= "Count differentially regulated", x= "", alpha="# Patterns") +
  scale_y_continuous(label=comma, limits=c(-1000, 1000)) 
barplot
ggsave("StackedBarPlotNew.pdf", height = 8, width = 4, units = "in")

#############################################
#load original DGE data
DGEgeneList <- list.load("DGEgeneList.rdata")

#scale the data and average the columns
scaleAvg <- function(x){
  c <- t(scale(t(x)))
  colnames(c) <- str_sub(colnames(c), 1, -2)
  c <- t(apply(c, 1, function(x) tapply(x, colnames(c), mean)))
}

#subset the counts to LM (1:6) or EE (7:12)
#or 3:6 / 9:12 for no 30v50
hm <- z %>% 
  subset(rownames(z) %in% unlist(DGEgeneList[9:12])) %>% 
  scaleAvg()
#hierarchally cluster and divide into clusters with cutreeDynamic
sd <- dist(hm)
hmc <- hclust(sd)
clusterSize <- 20
deepSplit <- 1
cut <- cutreeDynamic(hmc, minClusterSize = clusterSize,  
                     method = "hybrid", distM = as.matrix(sd), deepSplit = deepSplit,
                     pamStage = F, pamRespectsDendro = T)

#Plot dendrogram with colored clusters to see if there are too many/few clusters
plotDendroAndColors(
  hmc, labels2colors(cut), paste0(""),
  dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, 
  main = paste0(time," DEGs (n = ",nrow(hm),")"), marAll = c(1,5,4,2),
  autoColorHeight = F, colorHeight = 0.1)

#plot cluster to see which ones to remove
Clusters <- data.frame("Gene" = rownames(hm), "Cluster" = cut)[hmc$order,] %>% 
  filter(Cluster != 0)
  order <- unique(Clusters$Cluster)
  logCPM <- data.frame(hm) %>% 
    setDT(keep.rownames = "Gene") %>% 
    merge(Clusters, by = "Gene") %>% 
    pivot_longer(!c("Gene","Cluster"), names_to = "Sample", values_to = "CPM") %>% 
    arrange(Cluster)
  logCPM2 <- data.frame()
  for(i in 1:length(order)){
    temp <- logCPM %>% 
      filter(Cluster == order[i]) %>% 
      mutate(Cluster = i)
    logCPM2 <- rbind(logCPM2, temp)
  }
  logCPM2 <- logCPM
  logCPM2$Cluster <- as.character(logCPM2$Cluster)
  clustSum <- c()
  for(i in seq(1:length(unique(logCPM2$Cluster)))){
    clustSum[i] = sum(logCPM2$Cluster == i)/3
  }
  for(i in seq(1:nrow(logCPM2))){
    logCPM2[i,2] <- paste0("Cluster ",logCPM2[i,2]," (n=",clustSum[as.numeric(logCPM2[i,2])],")")
  }
  logCPM2$Sample <- factor(logCPM2$Sample, levels = c("LM30","LM50","LM80","EE30","EE50","EE80"))
  logCPM2$Cluster <- factor(logCPM2$Cluster, levels = unique(logCPM2$Cluster)[order])  
  ClustersMeans <- logCPM2 %>% 
    group_by(Cluster,Sample) %>% 
    dplyr::summarise(mean = mean(CPM))
  ggplot() +
    geom_line(data = logCPM2, aes(Sample, CPM, group = Gene)) +
    geom_line(data = ClustersMeans, aes(Sample, mean, group = Cluster), color = "red") +
    labs(title = paste0(time, " DEGs (p<0.05,lfc=0), ",length(unique(Clusters$Gene))," genes"), y = "z-score") +
    facet_wrap(~Cluster)
  ggsave("EEClusters9.16.21.png", width = 8, height = 8, units = "in")
  

###################################
#LM
#run this section for LM figure 
  
#remove unwanted clusters and assign others to patterns
Clusters <- data.frame("Gene" = rownames(hm), "Cluster" = cut)[hmc$order,]

Clusters$Group <- ifelse(Clusters$Cluster %in% c(4), "Down",
                         ifelse(Clusters$Cluster %in% c(10,5,6,1), "Down Transient",
                                ifelse(Clusters$Cluster %in% c(3,8,9), "Up",
                                       ifelse(Clusters$Cluster %in% c(7,2), "Up Transient", "Remove"))))
Clusters <- dplyr::arrange(Clusters, Gene)  

#scale and average clustered genes and combine with pattern assignment
logCPM <- t(scale(t(hm)))
logCPM <- t(apply(logCPM, 1, function(x) tapply(x, colnames(logCPM), mean)))
logCPM <- logCPM %>% 
  data.frame("Gene" = rownames(logCPM)) %>% 
  merge(Clusters, by = "Gene") %>% 
  pivot_longer(!c("Gene","Cluster","Group"), names_to = "Sample", values_to = "CPM") %>% 
  arrange(Cluster) 

#find the mean points of each cluster, assign labels
CPMmeans <- logCPM %>% 
  group_by(Cluster,Sample, Group) %>% 
  dplyr::summarise(mean = mean(CPM)) %>% 
  mutate(Group = as.character(Group)) %>% 
  filter(Group != "Remove")
CPMmeans$Cluster <- factor(CPMmeans$Cluster, levels = c(unique(CPMmeans$Cluster)))
CPMmeans$Sample <- factor(CPMmeans$Sample, levels = c("LM30", "LM50", "LM80"))
labs <- c("Up-regulated Only in Low\nWater Availability (n=226)",
          "Up-regulated with Decreasing\nWater Availability (n=227)",
          "Down-regulated Only in Low\nWater Availability (n=389)",
          "Down-regulated with Decreasing\nWater Availability (n=118)")
CPMmeans$Group <- revalue(CPMmeans$Group, c("Up Transient" = labs[1],
                                            "Up" = labs[2],
                                            "Down Transient" = labs[3],
                                            "Down" = labs[4]))
CPMmeans$Group <- factor(CPMmeans$Group, levels = labs[1:4])
text <- data.frame(labels = c("I","II","III","IV"), Group = factor(labs))

###################################
#EE
#run this section for EE figure

Clusters <- data.frame("Gene" = rownames(hm), "Cluster" = cut)[hmc$order,]
Clusters$Cluster <- revalue(as.character(Clusters$Cluster), c("13" = "0",
                                                              "12" = "0",
                                                              "10" = "0",
                                                              "5" = "0"))
Clusters$Group <- ifelse(Clusters$Cluster %in% c(6,4,7,2), "Down",
                         ifelse(Clusters$Cluster %in% c(9), "Down Transient",
                                ifelse(Clusters$Cluster %in% c(3,1), "Up",
                                       ifelse(Clusters$Cluster %in% c(8,11), "Up Transient", "Remove"))))
Clusters <- dplyr::arrange(Clusters, Gene)  

logCPM <- t(scale(t(hm)))
logCPM <- t(apply(logCPM, 1, function(x) tapply(x, colnames(logCPM), mean)))
logCPM <- logCPM %>% 
  data.frame("Gene" = rownames(logCPM)) %>% 
  merge(Clusters, by = "Gene") %>% 
  pivot_longer(!c("Gene","Cluster","Group"), names_to = "Sample", values_to = "CPM") %>% 
  arrange(Cluster) 

CPMmeans <- logCPM %>% 
  group_by(Cluster,Sample, Group) %>% 
  dplyr::summarise(mean = mean(CPM)) %>% 
  mutate(Group = as.character(Group)) %>% 
  filter(Group != "Remove")
CPMmeans$Cluster <- factor(CPMmeans$Cluster, levels = c(unique(CPMmeans$Cluster)))
CPMmeans$Sample <- factor(CPMmeans$Sample, levels = c("EE30", "EE50", "EE80"))
labs <- c("Up-regulated Only in Low\nWater Availability (n=104)",
          "Up-regulated with Decreasing\nWater Availability (n=490)",
          "Down-regulated Only in Low\nWater Availability (n=53)",
          "Down-regulated with Decreasing\nWater Availability (n=631)")
CPMmeans$Group <- revalue(CPMmeans$Group, c("Up Transient" = labs[1],
                                            "Up" = labs[2],
                                            "Down Transient" = labs[3],
                                            "Down" = labs[4]))
CPMmeans$Group <- factor(CPMmeans$Group, levels = labs[1:4])
text <- data.frame(labels = c("V","VI","VII","VIII"), Group = factor(labs))

###################################
#continue here for either time

#make dendrogram of clustered genes and assign colors
dend <- as.dendrogram(rev(hmc)) %>% 
  dendextend::set("branches_lwd", 0.15)  
groups <- Clusters$Group[order.dendrogram(dend)]
groups <- replace(groups, groups=="Remove", 0)
groups <- replace(groups, groups=="Down", 1)
groups <- replace(groups, groups=="Down Transient", 2)
groups <- replace(groups, groups=="Up", 3)
groups <- replace(groups, groups=="Up Transient", 4)
groups <- as.integer(groups)

cols <- brewer.pal(length(unique(groups)), "Set1")[-c(5:6)]

dend <- branches_attr_by_clusters(dend, groups, cols[c(2,1,3,4)]) 
order <- labels(dend)

#set up clustered genes for heatmap
hmlong <- as.data.frame(hm) %>% 
  mutate(GeneID = rownames(hm)) %>% 
  pivot_longer(colnames(hm), names_to = "sample", values_to = "L2FC") %>% 
  mutate(sample = factor(sample, levels = c("LM30","LM50","LM80","EE30","EE50","EE80")))

#plot expression patterns
expression <- ggplot(CPMmeans, aes(Sample, mean, group = Cluster,color = Group)) +
  geom_line(size = 1) +
  scale_color_manual(values = cols) +
  scale_y_continuous(position = "right") +
  facet_wrap(~Group, ncol = 1) +
  geom_label(data = text, mapping = aes(Inf,-Inf, label = labels), color = cols, fontface = "bold",
            hjust = 1, vjust = 0, inherit.aes = FALSE) +
  labs(y = "",x = "") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, size = 8),
        axis.ticks = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = "#D2D2D2"),
        strip.text = element_text(size = 10),
        panel.border = element_rect(color = "#D2D2D2"),
        panel.grid.minor = element_blank()) +
  coord_cartesian(clip = "off")
expression

#############################
#GO term enrichment

#read in GO term data
GO2geneID <- readMappings("GOList.txt")
geneID2GO <- inverseList(GO2geneID)
geneNames <- names(geneID2GO)

#for loop to get results from all clusters
#will return GO terms for MF and BP ontologies
Res <- data.frame()
Clusters.Removed <- filter(Clusters, Group != "Remove")
GroupIndex <- unique(Clusters.Removed$Group)
for(i in 1:length(GroupIndex)){
  #Read in list of genes from clusters
  geneList <- factor(as.integer(geneNames %in% filter(Clusters.Removed, Group == GroupIndex[i])$Gene))
  names(geneList) <- geneNames
  
  for(h in 1:2){
    ont <- c("MF","BP")
    #GO function; change "ontology" to MF, BP, or CC
    GOdata <- new("topGOdata", ontology = ont[h], allGenes = geneList, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
    
    #Run fisher test with default weighted algorithm and filter by pval
    resultFis <- runTest(GOdata, statistic = "fisher")
    sig <- data.frame(resultFis@score) #%>% 
    #filter(resultFis.score < 0.5)
    
    #Generate results table
    temp <- GenTable(GOdata, weightFisher = resultFis, orderBy = resultFis, topNodes = nrow(sig), numChar = 1000)
    temp$Group <- GroupIndex[i]
    Res <- rbind(Res, temp)
  }
}

allRes <- data.frame("Term" = Res$Term,
                     "Group" = Res$Group,
                     "pval" = as.numeric(Res$weightFisher))
allRes <- allRes[!duplicated(allRes[,1:2]),]

#write.csv(allRes, "allResLM9.16.21.csv", row.names = F)
#allRes <- read.csv("allResLM9.16.21.csv")
#write.csv(allRes, "allRes6.3.21EE.csv", row.names = F)

#load in GoTerms with unwanted terms removed
keep <- read.csv("LMGosKeep2.csv") %>% 
  filter(x != "")
#colnames(remove) <- c("Term", "NewTerm")

#set up GO term data for heatmap
distRes <- pivot_wider(allRes,  names_from = Group, values_from = pval)
#distRes <- filter(distRes, !Term %in% remove$x)
distResBASE <- distRes
distResBASE$min <- do.call(pmin, distResBASE)
distResBASE <- filter(distResBASE, min <= 0.05)
distResBASE <- distResBASE[,1:(ncol(distResBASE)-1)]
distResBASE[cbind(Term = FALSE, distResBASE[2:ncol(distResBASE)] >0.05)] <- 0.05
distResBase <- distResBASE %>% 
  dplyr::select(Term, "Up Transient", "Up", "Down Transient","Down") #%>% 
  #filter(Term %in% keep$x) #%>% 
  #merge(remove, by="Term") %>% 
  #mutate(Term = NewTerm) %>% 
  #dplyr::select(-NewTerm)
distResBase <- dplyr::arrange(distResBase, distResBase[,2:ncol(distResBase)]) 
distResBase$Term <- factor(distResBase$Term, levels = rev(distResBase$Term))
BaseUp <- pivot_longer(distResBase, -Term, names_to = "Group", values_to = "pval")
BaseUp$Group <- factor(BaseUp$Group, levels = c("Up Transient", "Up", "Down Transient","Down"))

#write.csv(BaseUp, "LMgos2.csv", row.names = F)
#BaseUp <- read.csv("LMgos2.csv")
#write.csv(BaseUp, "EEgos2.csv", row.names = F)
#BaseUp <- read.csv("LMgos2.csv")
#write.csv(unique(BaseUp$Term), "LMGosRemove9.17.21.csv", row.names = F)
#write.csv(unique(BaseUp$Term), "EEGosRemove3.csv", row.names = F)

#heatmap with ggplot
gos <- ggplot(BaseUp, aes(Group,Term)) +
  geom_tile(aes(fill = pval)) +
  scale_fill_viridis(name = "p-value",limit = c(0, 0.05), direction = -1,
                     guide=guide_colorbar(barheight=unit(0.4, "cm"),
                                          barwidth=unit(4,"cm"),
                                          title.position = "top")) + 
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = text$labels) +
  labs(y = "", x = "",title = "") +
  theme_minimal(base_size = 12) +
  theme(plot.margin = unit(c(0,0,0,1), "cm"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 12, face = "bold",color = cols),
        plot.title = element_text(hjust = 0.5),
        rect=element_rect(fill="transparent",color=NA),
        panel.grid.major=element_blank(),
        legend.position = "top",
        legend.title.align = 0.5)
gos
ggsave("LMgos9.24.21ALL.png", height = 12, width = 8, units = "in")


####################################
#put the figure together

#tree plot
tree <- ggplot(rev(dend), horiz=T, labels=F)+
  theme(plot.margin = margin(0, 0, 0,0, "pt") )+
  scale_x_continuous(expand=c(0.01,0.01))+
  theme(rect=element_rect(fill="transparent", color=NA),
        panel.border=element_blank())

#heatmap
heatmap <- ggplot(hmlong, aes(sample, factor(GeneID, levels=rev(order)), fill=L2FC))+
  geom_tile(height = 1.05)+
  scale_fill_viridis_c(guide=guide_colorbar(barheight=unit(0.4, "cm"),
                                            barwidth=unit(4,"cm"),
                                            title.position = "top"), 
                       name = "z-score")+
  scale_x_discrete(expand=c(0,0))+ scale_y_discrete(expand=c(0,0))+
  theme_minimal(base_size=12)+
  theme(axis.ticks=element_blank(), axis.title = element_blank(),
        axis.text.y=element_blank(), rect=element_rect(fill="transparent", color=NA),
        panel.border=element_blank(),legend.position="top",
        legend.title.align = 0.5)


#LM figure
combo <- ggdraw(xlim = c(0,0.91))+
  draw_plot(tree, x=0, width=0.126, y=0.0301, height=0.8635)+
  draw_plot(heatmap, x=0.1, width=0.2, y=0.03, height=0.985)+
  draw_plot(expression, x=0.32, width=0.22, y=0.008, height=0.888)+
  draw_plot(gos, x=0.5, width=0.396, y=0.012, height=1.03)
combo
ggsave("LMHeat9.24.21.png", plot=combo, height = 8, width = 13.4, units = "in")

#EE figure
combo <- ggdraw(xlim = c(0,0.91))+
  draw_plot(tree, x=0, width=0.126, y=0.0301, height=0.8645)+
  draw_plot(heatmap, x=0.1, width=0.2, y=0.03, height=0.985)+
  draw_plot(expression, x=0.32, width=0.22, y=0.008, height=0.888)+
  draw_plot(gos, x=0.5, width=0.4, y=0.012, height=1.03)
ggsave("EEHeat6.9.21.png", plot=combo, height = 8, width = 13.4, units = "in")

#####################################
#get genes from GO terms
library(GO.db)

#select pattern
pattern <- BaseUp %>% 
  filter(Group == "Up") %>% 
  filter(pval < 0.05)

#get GO IDs
goterms <- Term(GOTERM)
goids <- names(goterms[goterms %in% pattern$Term])

#genes from GO terms
GoGenes <- unlist(GO2geneID[names(GO2geneID) %in% goids])

GoGenes2 <- data.frame(GoGenes) %>% 
  dplyr::rename(Gene = GoGenes)

func <- read_delim("geneFunctions.txt", delim = "\t") %>% 
  distinct(locusName, .keep_all = TRUE) %>% 
  filter(locusName %in% GoGenes2$Gene) %>% 
  dplyr::select(locusName, best_arabi_defline) 
colnames(func) <- c("Gene", "Function")

GoGenes3 <- GoGenes2 %>% 
  merge(func, by = "Gene")
#######################################

#plot line graphs of clock genes
#Convert to logCPM
cl <- read.csv("clock_genes.csv")

z <- cpm(y, prior.count=1, log=T)
logCPM <- subset(z, rownames(z) %in% cl$Gene)
logCPM <- t(scale(t(logCPM)))
colnames(logCPM) <- str_sub(colnames(logCPM), 1, -2)
logCPM <- t(apply(logCPM, 1, function(x) tapply(x, colnames(logCPM), mean)))
logCPM <- data.frame(Gene = rownames(logCPM), logCPM)

logCPM <- logCPM %>% 
  merge(read.csv("clock_genes.csv"), by = "Gene")

#morning clock genes
logLM <- filter(logCPM, EE30 < 0) %>%
  as.data.frame() %>% 
  pivot_longer(-c(Gene, Function), names_to = "Sample", values_to = "score") 
logLM$Sample <- factor(logLM$Sample, levels = c("LM30","LM50","LM80",
                                                  "EE30","EE50","EE80"))
morning <- ggplot(logLM, aes(Sample, score, group = Gene)) +
  geom_line() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank()) +
  labs(y = "z-score", x = "", title = "Morning Clock Genes (n = 8)") +
  scale_colour_discrete("Gene")


#evening clock genes
logEE <- filter(logCPM, EE30 > 0) %>%
  as.data.frame() %>% 
  pivot_longer(-c(Gene, Function), names_to = "Sample", values_to = "score") 
logEE$Sample <- factor(logEE$Sample, levels = c("LM30","LM50","LM80",
                                                  "EE30","EE50","EE80"))
evening <- ggplot(logEE, aes(Sample, score, group = Gene)) +
  geom_line() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y = "z-score", x = "", title = "Evening Clock Genes (n = 6)") +
  scale_colour_discrete("Gene")

h <- grid.arrange(morning, evening) 
ggsave("clockGenes2.pdf", h, height = 2, width = 6.7, units = "in")

#combine the two
logLM$Time <- "Morning-phased (n = 8)"
logEE$Time <- "Evening-phased (n = 6)"
clock <- rbind(logLM, logEE)

ggplot(clock, aes(Sample, score, group = Gene, color = Time)) +
  geom_line() +
  labs(y = "z-score", x = "") +
  scale_color_manual(values = c("#E74C3C", "#3498DB"),"") +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.4,1.05), legend.direction = "horizontal",
        legend.key.height = unit(0.1, "cm"),
        legend.margin=margin(t = 0, unit='cm'), plot.margin = unit(x = c(0.5, 0, -0.2, 0), units = "cm"),
        legend.text = element_text(size=7)) +
  guides(color = guide_legend(reverse = TRUE))
ggsave("clockGenesNew.pdf", height = 2, width = 3.35, units = "in")

###################################
#find slope for all clusters

clustersDat <- unique(ClustersMeans$Cluster)
slope <- data.frame(Cluster = clustersDat,
                    Slope1 = NA,
                    Slope2 = NA,
                    SlopeSum = NA)
for(i in 1:length(clustersDat)){
  temp <- as.data.frame(filter(ClustersMeans, Cluster == clustersDat[i]))
  slope[i,2] <- lm(c(1,temp[1,3]) ~ c(2,temp[2,3]))$coeff[[2]] 
  slope[i,3] <- lm(c(2,temp[2,3]) ~ c(3,temp[3,3]))$coeff[[2]]
  slope[i,4] <- slope[i,2] + slope[i,3]
}

#find difference between last 2 points
  
clustersDat <- unique(ClustersMeans$Cluster)
slope <- data.frame(Cluster = clustersDat,
                    zDifference1.2 = NA,
                    zDifference2.3 = NA)
for(i in 1:length(clustersDat)){
  temp <- as.data.frame(filter(ClustersMeans, Cluster == clustersDat[i]))
  slope[i,2] <- temp[1,3] - temp[2,3]
  slope[i,3] <- temp[3,3] - temp[2,3]
}

#find correlations of each group
library(DescTools)

#get the average pearson correlation value for each cluster
#find cor between each gene, converts to Fisher z score, averages, then converts back to cor
cors <- data.frame("Clusters" = unique(Clusters$Cluster),
                   "Cor" = NA)
for(i in 1:length(unique(Clusters$Cluster))){
  clus <- filter(Clusters, Clusters$Cluster == unique(Clusters$Cluster)[i])
  clus2 <- hm[rownames(hm) %in% clus$Gene,]
  corl <- stats::cor(t(clus2))
  corlTrans <- FisherZ(corl)
  corlTrans[sapply(corlTrans, is.infinite)] <- NA
  meanCorl <- mean(corlTrans, na.rm = TRUE)
  meanCorl <- FisherZInv(meanCorl)
  
  cors[i,2] <- meanCorl
}
cors <- dplyr::arrange(cors, Cor)

##############################
#code for Figure 5

#alternative DGE list for any pval
#For top 1000
#will be used to rank top 1000 genes of each category
DGE <- data.frame(rep(NA, 1000))
DGEgeneList <- list()
for(i in 1:6){
  temp <- data.frame(exactTest(y, pair = data.frame(c(1,2),c(1,3),c(2,3),c(4,5),c(4,6),c(5,6))[,i]))
  tempUp <- temp %>% 
    filter(logFC > 0) %>% 
    dplyr::arrange(PValue)
  tempDown <- temp %>% 
    filter(logFC < 0) %>% 
    dplyr::arrange(PValue)
  DGE <- cbind(DGE, rownames(tempUp[1:1000,]), rownames(tempDown[1:1000,]))
}

#list of genes in each group; split into up and down
DGEgeneList <- list(
                    LM30.80up = DGE[,4], LM30.80down = DGE[,5],
                    LM50.80up = DGE[,6], LM50.80down = DGE[,7],
                    EE30.80up = DGE[,10], EE30.80down = DGE[,11],
                    EE50.80up = DGE[,12], EE50.80down = DGE[,13])

DGEgeneList <- list(
  LM30.80up = DGE[DGE[,4] %in% unlist(Quantitativetop1000[1]),4], LM30.80down = DGE[DGE[,5] %in% unlist(Quantitativetop1000[1]),5],
  LM50.80up = DGE[DGE[,6] %in% unlist(Quantitativetop1000[2]),6], LM50.80down = DGE[DGE[,7] %in% unlist(Quantitativetop1000[2]),7],
  EE30.80up = DGE[DGE[,10] %in% unlist(Quantitativetop1000[3]),10], EE30.80down = DGE[DGE[,11] %in% unlist(Quantitativetop1000[3]),11],
  EE50.80up = DGE[DGE[,12] %in% unlist(Quantitativetop1000[4]),12], EE50.80down = DGE[DGE[,13] %in% unlist(Quantitativetop1000[4]),13]
)


###find total number of unique genes
allDGE <- unlist(DGEgeneList)
allDGE <- unique(allDGE)


#for all DEGS
SingleTimeTable <- function(direction){
  LM30.80 <- unlist(DGEgeneList[paste0("LM30.80", direction)])
  LM50.80 <- unlist(DGEgeneList[paste0("LM50.80", direction)])
  EE30.80 <- unlist(DGEgeneList[paste0("EE30.80", direction)])
  EE50.80 <- unlist(DGEgeneList[paste0("EE50.80", direction)])
  
  totalgenes <- c(LM30.80,LM50.80,EE30.80,EE50.80)
  
  timetable <- data.frame(ATG=totalgenes, LM30.80=0,LM50.80=0,EE30.80=0,EE50.80=0) 
  timetable$LM30.80[timetable$ATG %in% LM30.80] <- 1
  timetable$LM50.80[timetable$ATG %in% LM50.80] <- 1
  timetable$EE30.80[timetable$ATG %in% EE30.80] <- 1 
  timetable$EE50.80[timetable$ATG %in% EE50.80] <- 1 
  
  timetable <- timetable[!duplicated(timetable),]
  rownames(timetable) <- timetable$ATG
  timetable <- timetable[-1]
  return(timetable) }

#make the up/down tables
up <- SingleTimeTable("up")
down <- SingleTimeTable("down")

#find common genes in 2 & 3
comb <- rbind(up, down)
comb$sum <- rowSums(comb)
table(comb$sum)

#see if the 2s are more common in same time of day
comb2 <- filter(comb, sum == 2)
combDup <- aggregate(cbind(comb2[0],dup=1), comb2, length)

#create dataframe with all data needed for plot
Direction <- c(rep("up",16), rep("down",16))
Comparison <- rep(c(rep("LM30.80",4), rep("LM50.80",4), rep("EE30.80",4), rep("EE50.80",4)),2)
Degree <- rep(seq(1:4),8)
Cardinality <- rep(0,32)
bar <- data.frame(Direction, Comparison, Degree, Cardinality)

#Filling in that table
comps <- c("LM30.80","LM50.80","EE30.80","EE50.80")[2]
FillInBar <- function(timetable, direction, pass){ 
  bar <- pass
  
  for(h in 1:4){
    temp <- timetable[timetable[,h]==1,]
    for(i in 1:7){
      bar$Cardinality[bar$Comparison== c("LM30.80","LM50.80","EE30.80","EE50.80")[h] & bar$Degree==i & bar$Direction== direction] <- length(rowSums(temp)[rowSums(temp)==i])
    }
  }
  
  return(bar)
}

bar <- FillInBar(up, "up", bar)
bar <- FillInBar(down, "down", bar)

#Make the downregulated DEG bars extend downwards
bar$Cardinality[bar$Direction=="down"] <- -bar$Cardinality[bar$Direction=="down"]

#Some parameter wrangling for better graphing
bar$Degree <- as.factor(bar$Degree)
bar$Comparison <- factor(bar$Comparison, levels = c("LM30.80","LM50.80","EE30.80","EE50.80"))
definedvibrant <- c("0"="white",  "LM30.80"="#0077BB", "LM50.80"="#009988", 
                     "EE30.80"="#DDCC77", "EE50.80"="#CC3311")

#actual plot
barplot <- ggplot(bar, aes(x=Comparison, y=Cardinality, group=Degree)) +
  geom_bar(aes(fill=Comparison), stat="identity", position=position_stack()) +
  scale_fill_manual(values=definedvibrant, guide="none") +
  geom_bar(aes(alpha=Degree), fill="black", stat="identity", position=position_stack()) +
  geom_bar(aes(alpha=as.factor(5-as.numeric(Degree))), fill="white", stat="identity", position=position_stack())  +
  scale_alpha_manual(values=c(0.05,0.1,0.2,0.4))+
  scale_x_discrete(labels = c("LM30v80","LM50v80","EE30v80","EE50v80")) +
  guides(alpha=guide_legend(ncol=4)) +
  theme_minimal(base_size=16) +
  theme(#axis.text.x=element_text(hjust=0.5, size=8), 
        panel.background = element_blank(), 
        plot.background = element_blank(), 
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        panel.border = element_blank(), 
        legend.box=NULL, legend.key.size=unit(0.4, "cm"), 
        legend.position ="top",
        #legend.text=element_text(size=10), 
        legend.title = element_text(size=12), 
        #axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y= "Count differentially regulated", x= "", alpha="Common Patterns") +
  scale_y_continuous(label=comma, limits=c(-500, 700)) 
barplot
ggsave("StackedBarPlotTop1000.pdf", height = 8, width = 4, units = "in")


#top 1000 genes in each comparison
top1000 <- data.frame("LM.30v80" = unlist(Quantitativetop1000[1]),
                      "LM.50v80" = unlist(Quantitativetop1000[2]),
                      "EE.30v80" = unlist(Quantitativetop1000[3]),
                      "EE.50v80" = unlist(Quantitativetop1000[4]))
write.csv(top1000, "top1000DEGs.csv", row.names = F)

##############################3
#code for PCA
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdataMatrix), Drought_Level, Time_of_Day))
design1 =~ Drought_Level + Time_of_Day
dds <- DESeqDataSetFromMatrix(countData=countdataMatrix, colData=coldata, design= design1)
dds

##Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)

# Modified plotPCA from DESeq2 package. Shows the Names of the Samples (the first col of SampleTable), and uses ggrepel pkg to plot them conveniently.
# @SA 10.02.2017 
library(genefilter)
library(ggplot2)
library(ggrepel)

plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(vsd)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}

newpca <-plotPCA.san(vsd, intgroup=c("Drought_Level", "Time_of_Day"))

k <- newpca + theme_bw()
k
col11 <- "#E74C3C"
col12 <- "#2ECC71"
col13 <- "#3498DB"
colvariation <- c( "30"=col11, "50"=col12, "80"=col13)

pcaData <-plotPCA.san(vsd, intgroup=c("Drought_Level", "Time_of_Day"), returnData=TRUE)
pcaData
df_out <- as.data.frame(pcaData)
df_out$Time_of_Day <- as.character(df_out$Time_of_Day)
df_out$Time_of_Day[df_out$Time_of_Day=="LM"] <- "LM (10:00)"
df_out$Time_of_Day[df_out$Time_of_Day=="EE"] <- "EE (16:00)"
df_out$Time_of_Day <- factor(df_out$Time_of_Day, levels=c("LM (10:00)","EE (16:00)")) 
df_out$Drought_Level <- factor(df_out$Drought_Level, levels = c("80", "50", "30"))
percentVar <- round(100 * attr(pcaData, "percentVar"))

PanelPC1and2 <- ggplot(df_out, aes(x=PC1, y=-PC2, shape=Time_of_Day, color= Drought_Level )) + geom_point(size=0.7, stroke=1)+
  scale_shape_manual(values=c(5,2)) +
  scale_color_manual(values = colvariation)+
  theme_minimal(base_size=12) +labs(x=paste0("PC1:", percentVar[1],"% variance"), y=paste0("PC2:", percentVar[2],"% variance"), shape = "Time of day",color= expression(paste(theta,"g")))+
  theme(rect=element_rect(fill="white"),legend.key.size=unit(0.35,"cm"), legend.spacing=unit(0.055,"cm"))+
  guides(colour=guide_legend(reverse = TRUE, order = 2))

PanelPC1and2
ggsave(PanelPC1and2, filename="2021sept11PanelPC1and2.pdf", height=7, width=9, units="cm", bg="white")


PanelPC1and3 <- ggplot(df_out, aes(x=PC1, y=-PC3, color= Drought_Level , shape=Time_of_Day)) + geom_point(size=0.7, stroke=1)+
  scale_shape_manual(values=c(5,2,3,1,0,4)) +
  scale_color_manual(values = colvariation)+
  theme_minimal(base_size=12) +labs(x=paste0("PC1:", percentVar[1],"% variance"), y=paste0("PC3:", percentVar[3],"% variance"), shape = "Time of day",color= expression(paste(theta,"g")))+
  theme(rect=element_rect(fill="white"),legend.key.size=unit(0.35,"cm"), legend.spacing=unit(0.055,"cm"))+
  guides(colour=guide_legend(reverse = TRUE))

PanelPC1and3
ggsave(PanelPC1and3, filename="2021sept11PanelPC1and3.pdf", height=8, width=16, units="cm", bg="white")


PanelPC2and3 <- ggplot(df_out, aes(x=PC2, y=-PC3, color= Drought_Level , shape=Time_of_Day)) + geom_point(size=0.7, stroke=1)+
  scale_shape_manual(values=c(5,2,3,1,0,4)) +
  scale_color_manual(values = colvariation)+
  theme_minimal(base_size=12) +labs(x=paste0("PC2:", percentVar[2],"% variance"), y=paste0("PC3:", percentVar[3],"% variance"), shape = "Time of day",color= expression(paste(theta,"g")))+
  theme(rect=element_rect(fill="white"),legend.key.size=unit(0.35,"cm"), legend.spacing=unit(0.055,"cm"))+
  guides(colour=guide_legend(reverse = TRUE))

PanelPC2and3
ggsave(PanelPC2and3, filename="2021sept11PanelPC2and3.pdf", height=8, width=16, units="cm", bg="white")

########################
#Quantitative shift plots
#LM30.80
is.de <-decideTestsDGE(exactTest(y, pair = c(1,3)))
summary(is.de)
tab <- topTags(exactTest(y, pair = c(1,3)), n=Inf)$table
dim(tab)
LM30.80 <- tab %>% dplyr::slice(1:1000) %>% 
  rownames_to_column("Geneid")
tab <- tab %>%  
  rownames_to_column("Geneid")
#write.csv(tab, "LM30v80_500.csv")

#LM50.80
is.de1 <-decideTestsDGE(exactTest(y, pair = c(2,3)))
summary(is.de1)
tab1 <- topTags(exactTest(y, pair = c(2,3)), n=Inf)$table
dim(tab1)
LM50.80 <- tab1 %>% dplyr::slice(1:1000) %>% 
  rownames_to_column("Geneid")
tab1 <- tab1 %>%  
  rownames_to_column("Geneid")
#write.csv(tab1, "LM50v80_500.csv")

#EE30.80
is.de2 <-decideTestsDGE(exactTest(y, pair = c(4,6)))
summary(is.de2)
tab2 <- topTags(exactTest(y, pair = c(4,6)), n=Inf)$table
dim(tab2)
EE30.80 <- tab2 %>% dplyr::slice(1:1000) %>% 
  rownames_to_column("Geneid")
tab2 <- tab2 %>%  
  rownames_to_column("Geneid")
#write.csv(tab2, "EE30v80.csv")

#EE50.80
is.de3 <-decideTestsDGE(exactTest(y, pair = c(5,6)))
summary(is.de3)
tab3 <- topTags(exactTest(y, pair = c(5,6)), n=Inf)$table
dim(tab3)
EE50.80 <- tab3 %>% dplyr::slice(1:1000) %>% 
  rownames_to_column("Geneid")
tab3 <- tab3 %>%  
  rownames_to_column("Geneid")
#write.csv(tab3, "EE50v80.csv")

Quantitativetop1000 <- list("LM30.80" = LM30.80$Geneid,
                            "LM50.80" = LM50.80$Geneid,
                            "EE30.80" = EE30.80$Geneid,
                            "EE50.80" = EE50.80$Geneid)

combined <- unique(unlist(Quantitativetop1000))

####Plot 1 LM30.80 vs LM50.80
plotdata1 <- merge(tab[1:2],tab1[1:2], by ="Geneid")
colnames(plotdata1) <- c("Geneid", "LM30v80logFC", "LM50v80logFC")

plotlfc1 <- plotdata1 %>%
  subset(plotdata1$Geneid %in% unlist(Quantitativetop1000[1:2]))
plotlfc1$LM30v80logFC <- as.numeric(plotlfc1$LM30v80logFC)
plotlfc1$LM50v80logFC <- as.numeric(plotlfc1$LM50v80logFC)

m0 <- lm(LM30v80logFC ~ LM50v80logFC, data=plotlfc1)
summary(m0)
mm0 <- cor.test(plotlfc1$LM50v80logFC, plotlfc1$LM30v80logFC,
                 data=plotlfc1,
                 method = "spearman",
                 continuity = FALSE,
                 conf.level = 0.95)
mm0

plotk21 <- ggplot(plotlfc1, aes(x=LM50v80logFC, y=LM30v80logFC))+
  geom_point(size=0.2)+geom_smooth(method="lm", formula=y~x, color="gray")+
  xlim(-5, 7)+
  ylim(-5, 7)+
  labs(x=expression('LM50v80 log'[2]*'FC'),y=expression('LM30v80 log'[2]*'FC'))+
  #annotate(geom="text", label=bquote('rho'*' = '*.(round((mm0)$estimate, 3))), x=0.8*min(plotlfc1$LM30v80logFC, na.rm=T), y=0.8*max(plotlfc1$LM30v80logFC, na.rm=T),hjust=0, size=1.5)+
  annotate(geom="text", label=c(paste('rho =', round(mm0$estimate, 2)), paste('Slope =', round(unname(m0$coefficients[2]), 2)), paste('n =', nrow(plotlfc1))) , x=Inf, y=c(-2.5, -3.5, -4.5), size= 4, fontface = "bold",hjust = 1, vjust = 0)+
  theme_minimal(base_size=16) +
  geom_hline(yintercept = -5, color = "#009988", size = 1) +
  geom_vline(xintercept = -5, color = "#0077BB", size = 1)
plotk21


####Plot 2 EE30.80 vs EE50.80
plotdata2 <- merge(tab2[1:2],tab3[1:2], by ="Geneid")
colnames(plotdata2) <- c("Geneid", "EE30v80logFC", "EE50v80logFC")

plotlfc2 <- plotdata2 %>%
  subset(plotdata2$Geneid %in% unlist(Quantitativetop1000[3:4]))
plotlfc2$EE30v80logFC <- as.numeric(plotlfc2$EE30v80logFC)
plotlfc2$EE50v80logFC <- as.numeric(plotlfc2$EE50v80logFC)


#lm test (y~x)
m1 <- lm(EE30v80logFC ~ EE50v80logFC, data=plotlfc2)
summary(m1)

#cor.test(x,y)
mm1 <- cor.test(plotlfc2$EE50v80logFC, plotlfc2$EE30v80logFC,
                 data=plotlfc2,
                 method = "spearman",
                 continuity = FALSE,
                 conf.level = 0.95)
mm1

plotk22 <- ggplot(plotlfc2, aes(x=EE50v80logFC, y=EE30v80logFC))+
  geom_point(size=0.2)+geom_smooth(method="lm", formula=y~x, color="gray")+
  xlim(-5, 7)+
  ylim(-5, 7)+
  labs(x=expression('EE50v80 log'[2]*'FC'),y=expression('EE30v80 log'[2]*'FC'))+
  #annotate(geom="text", label=bquote('rho'*' = '*.(round((mm1)$estimate, 3))), x=0.8*min(plotlfc2$EE30v80logFC, na.rm=T), y=0.8*max(plotlfc2$EE30v80logFC, na.rm=T),hjust=0, size=1.5)+
  annotate(geom="text", label=c(paste('rho =', round(mm1$estimate, 2)), paste('Slope =', round(unname(m1$coefficients[2]), 2)), paste('n =', nrow(plotlfc2))) , x=Inf, y=c(-2.5, -3.5, -4.5), size= 4, fontface = "bold",hjust = 1, vjust = 0)+
  theme_minimal(base_size=16) +
  geom_hline(yintercept = -5, color = "#CC3311", size = 1) +
  geom_vline(xintercept = -5, color = "#DDCC77", size = 1)
plotk22

#############################################################
####Plot 3 LM30.80 vs EE30.80

plotdata3 <- merge(tab[1:2],tab2[1:2], by ="Geneid")
colnames(plotdata3) <- c("Geneid", "LM30v80logFC", "EE30v80logFC")

plotlfc3 <- plotdata3 %>%
  subset(plotdata3$Geneid %in% unlist(Quantitativetop1000[c(1,3)]))
plotlfc3$LM30v80logFC <- as.numeric(plotlfc3$LM30v80logFC)
plotlfc3$EE30v80logFC <- as.numeric(plotlfc3$EE30v80logFC)

#lm test (y~x)
m2 <- lm(LM30v80logFC ~ EE30v80logFC, data=plotlfc3)
summary(m2)

#cor.test(x,y)
mm2 <- cor.test(plotlfc3$EE30v80logFC, plotlfc3$LM30v80logFC,
                data=plotlfc3,
                method = "spearman",
                continuity = FALSE,
                conf.level = 0.95)
mm2

plotk23 <- ggplot(plotlfc3, aes(x=EE30v80logFC, y=LM30v80logFC))+
  geom_point(size=0.2)+geom_smooth(method="lm", formula=y~x, color="gray")+
  xlim(-5, 7)+
  ylim(-5, 7)+
  labs(x=expression('EE30v80 log'[2]*'FC'),y=expression('LM30v80 log'[2]*'FC'))+
  #annotate(geom="text", label=bquote('rho'*' = '*.(round((mm2)$estimate, 3))), x=0.8*min(plotlfc3$EE30v80logFC, na.rm=T), y=0.8*max(plotlfc3$EE30v80logFC, na.rm=T),hjust=0, size=1.5)+
  annotate(geom="text", label=c(paste('rho =', round(mm2$estimate, 2)), paste('Slope =', round(unname(m2$coefficients[2]), 2)), paste('n =', nrow(plotlfc3))) , x=Inf, y=c(-2.5, -3.5, -4.5), size= 4, fontface = "bold",hjust = 1, vjust = 0)+
  theme_minimal(base_size=16) +
  geom_hline(yintercept = -5, color = "#DDCC77", size = 1) +
  geom_vline(xintercept = -5, color = "#0077BB", size = 1)
plotk23

####Plot 4 EE30.80 vs EE50.80

plotdata4 <- merge(tab1[1:2],tab3[1:2], by ="Geneid")
colnames(plotdata4) <- c("Geneid", "LM50v80logFC", "EE50v80logFC")

plotlfc4 <- plotdata4 %>%
  subset(plotdata4$Geneid %in% unlist(Quantitativetop1000[c(2,4)]))
plotlfc4$LM50v80logFC <- as.numeric(plotlfc4$LM50v80logFC)
plotlfc4$EE50v80logFC <- as.numeric(plotlfc4$EE50v80logFC)

#lm test (y~x)
m4 <- lm(LM50v80logFC ~ EE50v80logFC, data=plotlfc4)
summary(m4)

#cor.test(x,y)
mm4 <- cor.test(plotlfc4$EE50v80logFC, plotlfc4$LM50v80logFC,
                data=plotlfc4,
                method = "spearman",
                continuity = FALSE,
                conf.level = 0.95)
mm4

plotk24 <- ggplot(plotlfc4, aes(x=EE50v80logFC, y=LM50v80logFC))+
  geom_point(size=0.4, alpha = 0.2)+geom_smooth(method="lm", formula=y~x, color="gray")+
  xlim(-5, 7)+
  ylim(-5, 7)+
  labs(x=expression('EE50v80 log'[2]*'FC'),y=expression('LM50v80 log'[2]*'FC'))+
  #annotate(geom="text", label=bquote('rho'*' = '*.(round((mm4)$estimate, 3))), x=0.8*min(plotlfc4$EE30v80logFC, na.rm=T), y=0.8*max(plotlfc4$EE30v80logFC, na.rm=T),hjust=0, size=1.5)+
  annotate(geom="text", label=c(paste('rho =', round(mm4$estimate, 2)), paste('Slope =', round(unname(m4$coefficients[2]), 2)), paste('n =', nrow(plotlfc4))) , x=Inf, y=c(-2.5, -3.5, -4.5), size= 4, fontface = "bold",hjust = 1, vjust = 0)+
  theme_minimal(base_size=16) +
  geom_hline(yintercept = -5, color = "#CC3311", size = 1) +
  geom_vline(xintercept = -5, color = "#009988", size = 1)
plotk24

#library(cowplot)
top1000genes <- plot_grid(
  plotk21, plotk22, plotk23, plotk24, ncol = 2)
top1000genes


ggsave("sep14top1000genesnew2.pdf", top1000genes,height=8, width=8, units="in")
