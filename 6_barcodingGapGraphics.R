## REQUIRED LIBRARIES

library(readxl)
library(tidyverse)
library(usefun)
library(stats)
library(ggpubr)
library(hrbrthemes)
library(NLP)
library(tm)
library(proustr)
library(viridisLite)
library(viridis)
library(forcats)
library(purrr)
library(gridExtra)

#set working directory
setwd("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/")

#load data 
Distance1<- read.delim("Sequences_Folmer_pairs.dist", header=FALSE)
colnames(Distance1) <- c("Comparison","Taxonomy_sp1","sp1","Taxonomy_sp2","sp2","Similarity")
Distance1$Barcode <- rep(c("Folmer"),nrow(Distance1))

Distance2<- read.delim("Sequences_Leray_pairs.dist", header=FALSE)
colnames(Distance2) <- c("Comparison","Taxonomy_sp1","sp1","Taxonomy_sp2","sp2","Similarity")
Distance2$Barcode <- rep(c("Leray"),nrow(Distance2))

Distance3<- read.delim("Sequences_MiFish_pairs.dist", header=FALSE)
colnames(Distance3) <- c("Comparison","Taxonomy_sp1","sp1","Taxonomy_sp2","sp2","Similarity")
Distance3$Barcode <- rep(c("MiFish"),nrow(Distance3))

Distance4<- read.delim("Sequences_Teleo_pairs.dist", header=FALSE)
colnames(Distance4) <- c("Comparison","Taxonomy_sp1","sp1","Taxonomy_sp2","sp2","Similarity")
Distance4$Barcode <- rep(c("Teleo"),nrow(Distance4))

#Calculate descriptive statistics 
Folmer.data <- Distance1 %>% dplyr::select(Comparison, Similarity)
Leray.data <-  Distance2 %>% dplyr::select(Comparison, Similarity)
Mifish.data <- Distance3 %>% dplyr::select(Comparison, Similarity)
Teleo.data <-  Distance4 %>% dplyr::select(Comparison, Similarity)

data <- c()
for (barcode in c("Folmer", "Leray", "Mifish", "Teleo")) {
  eval(parse(text=paste0("i<-", barcode,".data[",barcode,".data$Comparison == 'SP','Similarity']")))
  eval(parse(text = paste0("y<-min(",barcode,".data$'Similarity')")))
  data <- rbind(data, c(unlist(boxplot.stats(i, coef = 1.5)[1]), y,(length(unlist(boxplot.stats(i, coef = 1.5)[4]))/length(i)*100)))
}
data <- as.data.frame(data)
colnames(data) <- c("min", "Q1", "median", "Q3", "max", "min_value", "%_outliers") #outliers 7% of the data for normal distribution
data$barcode <- c("Folmer", "Leray", "Mifish", "Teleo")
data #boxplot values for intraspecific distances 

#Set  ranges for each level 
ranges<-c()
for (barcode in c("Folmer", "Leray", "Mifish", "Teleo")){
  V2<-round(data$min_value + (data$min - data$min_value)* 3/4)
  V3<-round(data$min_value + (data$min - data$min_value)* 2/4)
  V4<-round(data$min_value + (data$min - data$min_value)* 1/4)
  Barcode<- c("Folmer", "Leray","MiFish","Teleo")
  ranges <- data.frame(cbind(round(data$min),V2, V3,V4,round(data$min_value), Barcode))
}
sapply(ranges, class)
ranges<- transform(ranges, V1=as.numeric(V1)) %>% transform(ranges, V2=as.numeric(V2)) %>% transform(ranges, V3=as.numeric(V3)) %>% transform(ranges, V4=as.numeric(V4))  %>% transform(ranges, V5=as.numeric(V5)) 
ranges<-ranges[,1:6]
colnnames(ranges)
ranges #Lower limits for each level 

# PLOT THE DATA

#BOXPLOT
data.boxplot <- do.call("rbind", list(Distance1, Distance2, Distance3, Distance4))
data.boxplot$Barcode <- factor(data.boxplot$Barcode, levels=c("Folmer","MiFish","Leray","Teleo"))

data.boxplot %>%
  mutate(Comparison = fct_relevel(Comparison, "SP", "GE", "FA", "OR", "CL", "PH")) %>%
  ggplot( aes(x=Comparison, y=Similarity)) +
  geom_hline(data = ranges, aes(yintercept=V1, linetype="dashed"))+ #level 1
  geom_hline(data = ranges, aes(yintercept=V2, linetype="solid"))+ #level 2:5
  geom_hline(data = ranges, aes(yintercept=V3, linetype="solid"))+
  geom_hline(data = ranges, aes(yintercept=V4, linetype="solid"))+
  geom_hline(data = ranges, aes(yintercept=V5, linetype="solid"))+
  geom_boxplot(
    color="blue",
    fill="blue",
    alpha=0.2,
    notch=F,
    outlier.colour="red",
    outlier.fill=NULL,
    outlier.shape = 1,
    outlier.size=0.5)+
  facet_wrap(~ Barcode, ncol=2)+
  theme_minimal()

#HEATMAPS
data2<- c()
for (comp in c("SP", "GE", "FA", "OR", "CL", "PH")) {
  for (barcode in c("Folmer", "Leray", "Mifish", "Teleo")) {
    eval(parse(text = paste0("dists <- ", barcode, ".data[", barcode, ".data$Comparison == comp, 'Similarity']")))
    e1 <- length(dists[dists <=100 & dists >= data[data$barcode == barcode, 1]])# Level 1 
    e2 <- length(dists[dists < data[data$barcode == barcode, 1] & dists >= (data[data$barcode == barcode, 6] + ((data[data$barcode == barcode, 1]-data[data$barcode == barcode, 6])*3/4))]) # Level 2
    e3 <- length(dists[dists < (data[data$barcode == barcode, 6] + ((data[data$barcode == barcode, 1]-data[data$barcode == barcode, 6])*3/4)) & dists >= (data[data$barcode == barcode, 6] + ((data[data$barcode == barcode, 1]-data[data$barcode == barcode, 6])*2/4))]) # Level 3
    e4 <- length(dists[dists < (data[data$barcode == barcode, 6] + ((data[data$barcode == barcode, 1]-data[data$barcode == barcode, 6])*2/4)) & dists >=  (data[data$barcode == barcode, 6] + ((data[data$barcode == barcode, 1]-data[data$barcode == barcode, 6])*1/4))]) # Level 4
    e5 <- length(dists[dists < (data[data$barcode == barcode, 6] + ((data[data$barcode == barcode, 1]-data[data$barcode == barcode, 6])*1/4)) & dists >= data[data$barcode == barcode, 6]]) # Level 5
    data2 <- rbind(data2, c(e1, e2, e3, e4, e5, comp, barcode))
  }
}
data2 <- as.data.frame(data2)
colnames(data2) <- c("Level 1", "Level 2", "Level 3", "Level 4", "Level 5", "Comparison", "Barcode")
data2 #number of comparisons for each level for each barcode
data2.elongated <- as.data.frame(pivot_longer(data2, cols=1:5, names_to = "Range", values_to = "Values"))
data2.elongated$Values <- as.numeric(as.character(data2.elongated$Values))
data2.elongated$Comparison <- fct_relevel(data2.elongated$Comparison,"SP","GE","FA","OR","CL","PH")
data2.elongated$Range <- fct_relevel(data2.elongated$Range, "Level 5", "Level 4","Level 3", "Level 2", "Level 1")
data2.elongated[data2.elongated == 0] <- NA
order <- c("Folmer", "Mifish", "Leray", "Teleo")
data2.elongated<-data2.elongated %>% arrange(sapply(Barcode, function(y) which(y == order)))

lapply(unique(data2.elongated$Barcode), function(cc) {
  ggplot(filter(data2.elongated, Barcode==cc),
         (aes(x=Comparison,y=Range, fill= Values, frame=Barcode))) +
    geom_tile(color="white", size=0.1)+
    labs(x=NULL, y=NULL, fill="Comparisons",
         title=sprintf("%s", cc))+
    geom_text(aes(label = Values))+
    scale_fill_gradient(low="gray95", high="cornflowerblue", na.value = "white")+
    xlab("")+ ylab("")+
    theme_ipsum()
}) -> cclist
cclist[["ncol"]] <- 2
do.call(grid.arrange, cclist) #recommended size 1400x770


