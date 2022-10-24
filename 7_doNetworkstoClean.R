
## REQUIRED LIBRARIES
library(ape)
library(dplyr)
library(tidyverse)
library(igraph)
library(RColorBrewer)
library(visNetwork)
library(qgraph)
library(tidyverse)
library(usefun)
library(Biostrings)

## Function to create networks
do_network <- function(data, name){
  edges <- data %>% select(sp1,sp2, Comparison, Similarity)
  names(edges) <- c("from", "to", "Comparison", "Similarity")
  edges$label <- as.character(round(edges$Similarity,2))
  nodes1 <- data %>% select(species1, sp1)
  names(nodes1)<- c("Species","Accnumber")
  nodes2 <- data %>% select(species2, sp2)
  names(nodes2)<- c("Species", "Accnumber")
  nodes <- rbind(nodes1, nodes2)
  rm(nodes1);rm(nodes2)
  nodes<- unique(nodes)
  nodes <- nodes[, c(2, 1)]
  links <- data %>% select(sp1,sp2) 
  network <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  deg <- degree(network, mode="all")
  nodes$deg <- deg
  names(nodes) <- c("id", "Species", "Connection")
  nodes$label = str_replace(nodes$Species, "_", " ")
  nodes$title = paste0("<p><b>", nodes$id,"</b><br>",nodes$Species,"</p>")
  nb.cols<-length(unique(nodes$Species))
  col <- colorRampPalette(brewer.pal(nlevels(as.factor(nodes$Species)),"Set3"))(nb.cols)
  my_color <- col[as.numeric(as.factor(V(network)$Species))]
  nodes$color = my_color
  plot <- visNetwork(nodes, edges, height = "500px", main = name, width = "100%") %>% 
    visOptions(highlightNearest = list(enabled = T, degree = 2, hover = T), selectedBy = "Connection", nodesIdSelection = TRUE) %>%
    visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>% visLayout(randomSeed = 1) %>%
    visEdges(smooth = FALSE) %>%   visPhysics(solver = "barnesHut", stabilization = FALSE) %>% visIgraphLayout()
  visSave(plot, file = paste0(name))
  clusters<- components(network) %>% groups()
  return(list("nodes"=nodes, "clusters"=clusters, "edges"=edges))
}

# Set working directory

setwd("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/")

#Load ranges 
ranges <- read.delim("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/Barcode_ranges.txt", header=F)

#Load reference database to be curated. Load Teleo, Mifish or DB to wish (both tax and dist files)
### 1: Teleo DB 
DB_raw <- read.delim("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/Sequences_Teleo.tax", header=FALSE)
Teleo <- read.delim("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/Sequences_Teleo_pairs.dist", header=FALSE)
data <- Teleo
name <- "Teleo"

### 2. Mifish DB
DB_raw <- read.delim("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/Sequences_Mifish.tax", header=FALSE)
MiFish <- read.delim("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/Sequences_Mifish_pairs.dist", header=FALSE)
data <- MiFish
name <- "MiFish"

# Prepare df 
colnames(data) <- c("Comparison","Taxonomy_sp1","sp1","Taxonomy_sp2","sp2","Similarity")
for (i in name) {
  rngs <- ranges %>% filter(ranges$V6 == i)
  data <- separate(data, Taxonomy_sp1, c(NA,"phylum1","class1","order1","family1","genus1","species1"), sep = ";" )
  data <- separate(data, Taxonomy_sp2, c(NA,"phylum2","class2","order2","family2","genus2","species2"),sep = ";" )
  data$Level<-ifelse(data$Similarity >= rngs[1,1], "Level_1",
                            ifelse(data$Similarity < rngs[1,1] & data$Similarity >= rngs[1,2], "Level_2",
                                   ifelse(data$Similarity < rngs[1,2] & data$Similarity >= rngs[1,3], "Level_3",
                                          ifelse(data$Similarity < rngs[1,3] & data$Similarity >= rngs[1,4], "Level_4", "Level_5"))))
  intra = data[data$Comparison =="SP" ,]
  intra = intra %>% select("Comparison","class1","order1","family1","genus1","species1","sp1","species2","sp2","Similarity","Level")
  inter = data[data$Comparison !="SP" ,]
  DB = DB_raw %>% separate(V2, c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
  summary <- DB %>% group_by(Species) %>% summarise(total_sequences=n())
  tax = unique(DB[,3:8])
  summary <- merge(tax, summary, by = "Species", all = T)
}


## NETWORK ANALYSIS AND POTENTIAL ERROR DETECTION 

#IMPORTANT! Set parameters to wish according to heatmap boxes to be included . 
intralevels <-3:5 #intraspecific levels to be analysed

interlevels <- c("Level_1", "Level_2") #interspecific levels to be analysed
intertax <- c("PH","CL", "OR")

#create empty df to be filled 
inter.net <- c()
intra.net <- c()
erroneous.sequences <-c()
problematic.sequences <-c()

{for (m in intralevels) {
  eval(parse(text=paste0("l='Level_",m,"'")))
  intra.net <- intra[intra$Level == l, ]
  if (dim(intra.net)[1]>0) {
    print(paste0("I N T R A S P E C I F I C   A N A L Y S I S: ",l))
    filename <- paste0(name,"_SP_",l,"_network.html")
    output <- do_network(data = intra.net, name = filename)
    for (i in 1:dim(output$clusters)) {
      eval(parse(text=paste0("n<-c(output$clusters$'",i,"')")))
      cluster <- output$nodes %>% select("id","Species", "Connection")
      cluster <- cluster %>% filter(id %in% n)
      a <- merge(cluster, summary, by= "Species", all=F)
      if (length(n)>2) {
        b <- a[a$Connection == nrow(a)-1,]
        if (dim(b)[1]>0) {
          print(paste("The cluster",i,"is generated by central sequence"))
          seq_central <- separate(inter, Level, c(NA, "Level"), sep ="_") %>% filter(sp1==b$id | sp2== b$id) %>% filter(Level < m)
          if (dim(seq_central)[1] > 0) {
            names = unique(c(seq_central$species1, seq_central$species2))
            names = names[!names %in% unique(b$Species)]
            print(paste("Central sequence",b$id,"is more similar to other species than to sequences of its own species", c(seq_central$id)))
            erroneous.sequences <- rbind(erroneous.sequences, c(b[1,], paste("Sequence more similar to:", do.call(paste, c(as.list(names), sep = ", ")))))
          }
        } else {
          print(paste("There is no central sequence in cluster", i))
          print("Network visualization is desirable")
          }
      } else {
        print(paste("The cluster",i,"is composed by a single pair (n=2), and no central sequence exists"))
        }
      }
    }
  }
  for (e in interlevels) {
    print(paste0("I N T E R S P E C I F I C   A N A L Y S I S: ",e))
    for (t in intertax) {
      print(paste("T A X O N O M I C  L E V E L:", t))
      inter.net <- inter[inter$Level == e, ]
      inter.net <- inter.net[inter.net$Comparison == t, ]
      if (dim(inter.net)[1]>0){
        filename <- paste0(name,"_",t,"_",e,"_network.html")
        output <- do_network(data = inter.net, name = filename)
        for (i in 1:dim(output$clusters)){
          eval(parse(text=paste0("n<-c(output$clusters$'",i,"')")))
          cluster <- output$nodes %>% select("id","Species", "Connection") %>% filter(id %in% n)
          a <- merge(cluster, summary, by= "Species", all=F)
          if (length(n)>2) {
            b <- a[a$Connection == nrow(a)-1,]
            if (dim(b)[1]>0) {
              print(paste("The cluster",i,"is generated by central sequence"))
              if (length(n)==a$total_sequences[1]) {
                print(paste0("There are NO more sequences for species in cluster ", i))
                problematic.sequences <- rbind(problematic.sequences, c(b, "The central sequence"))
              } else {
                print(paste0("There are more sequences for species in cluster ", i))
                erroneous.sequences <- rbind(erroneous.sequences, c(b, "The central sequence"))
              }
            } else {
              print(paste("There is NO central sequence in cluster",i,", and complex structure occurs"))
              sp <- unique(a$Species)
              if (length(sp) == 2) {
                print("Two species compose the complex strucutre")
                for (f in sp) {
                  x <- a[a$Species == f, ]
                  ids = a %>% filter(Species == f) %>% select(id)
                  N = a %>% filter(Species == f) %>% select(total_sequences) %>% unique() %>% as.numeric()
                  if (as.numeric(count(ids)) == N) {
                    print("There are no more sequences of this species out of the cluster")
                    g <- summary %>% filter(Species == f) %>% select(Genus)
                    y <- DB %>% filter(Genus == g$Genus) %>% select(Species) %>% unique()
                    if (dim(y)[1] > 1) {
                      print(paste0("There are/is ", (dim(y)[1]-1) ," more species (sequenced) in the genus"))
                      gg <- inter %>% filter(genus1 == g$Genus) %>% filter(Comparison == "GE")
                      if (dim(gg)[1] == 0) {
                        print("Sequences are too diferent to others from the genus to be compared")
                        for (q in 1:dim(x)[1]) {
                          s <- x[q, ]
                          erroneous.sequences <- rbind(erroneous.sequences, c(s, "Intra-genus distances are bigger than expected A"))
                        } 
                      } else {
                        ggg <- gg %>% filter(species1 == f | species2 == f)
                        if (dim(ggg)[1] == 0) {
                          print("Sequences are too diferent to others from the genus to be compared")
                          for (q in 1:dim(x)[1]) {
                            s <- x[q, ]
                            erroneous.sequences <- rbind(erroneous.sequences, c(s, "Intra-genus distance are bigger than expected B"))
                          } 
                        }
                      }
                    } else {
                      print("There are no more species in the genus or they are not sequenced")
                      problematic.sequences <- rbind(problematic.sequences, c(x[1,], "There are no more sequences in the genus to conclude (id non relevant)"))
                    }
                  } else {
                    print(paste0("There are/is ", (N - as.numeric(count(ids))) ," more sequence of this species"))
                    u <- intra %>% filter(species1 == f) %>% subset(sp1 %in% ids$id | sp2 %in% ids$id) %>% select(Level)
                    if (any(u$Level %in% c("Level_1", "Level_2"))) {
                      print("Intra-0sp distances are within the normal range")
                    } else {
                      print("Intra-sp distances are bigger than expected")
                      for (q in 1:dim(x)[1]) {
                        s <- x[q, ]
                        erroneous.sequences <- rbind(erroneous.sequences, c(s, "Intra-sp distances are bigger than expected C"))
                      }
                    }
                  }
                }
              }
            }
          } else {
            print(paste("The cluster",i,"is a pair (n = 2)"))
            for (f in a$Species){
              s <- a %>% filter(Species == f)
              if (summary[summary$Species == f, 7] == 1) {
                print(paste("There are no more sequences of", f))
                g <- summary %>% filter(Species == f) %>% select(Genus)
                y <- DB %>% filter(Genus == g$Genus) %>% select(Species) %>% unique()
                if (dim(y)[1] > 1) {
                  print(paste0("There are/is ", (dim(y)[1]-1) ," more species (sequenced) in the genus"))
                  gg <- inter %>% filter(genus1 == g$Genus) %>% filter(Comparison == "GE")
                  if (dim(gg)[1] == 0) {
                    print("Sequences are too diferent to others from the genus to be compared (Level 6)")
                    erroneous.sequences <- rbind(erroneous.sequences, c(s, "Intra-genus distances are bigger than expected"))
                  }
                } else {
                  print("There are no more species in the genus, or they are not sequenced")
                  problematic.sequences <- rbind(problematic.sequences, c(s, "There are no more sequences in the genus to reach a conclusion"))
                }
              } else {
                print(paste("There is/are", summary[summary$Species == f, 7] -1 ,"more sequences of", f))
                u <- intra %>% filter(species1 == f) %>% subset(sp1 %in% ids$id | sp2 %in% ids$id) %>% select(Level)
                if (any(u$Level %in% c("Level_1", "Level_2"))) {
                  print("Intra-genus distances are within the normal range")
                } else {
                  print("Intra-genus distances are bigger than expected")
                  erroneous.sequences <- rbind(erroneous.sequences, c(s, "Intra-genus distances are bigger than expected"))
                }
              }
            }
          }
        }
      }
    }
  }
erroneous.sequences <- unique(erroneous.sequences)
write.table(problematic.sequences, file = paste0("Problematic_sequences_",name,"_database.txt"), sep="\t", row.names = F)
write.table(erroneous.sequences, file = paste0("Removed_sequences_",name,"_database.txt"), sep="\t", row.names = F)
delete = data.frame(erroneous.sequences) %>% select(id) %>% unlist() %>% unique()
print(paste0("IN THE ",name," DATASET, ",length(delete)," ERRONEOUS SEQUENCES ARE DETECTED AND REMOVED"))
DB_clean = DB %>% subset(!V1 %in% delete)
summary_clean <- DB_clean %>% group_by(Species) %>% summarise(total_sequences=n())
tax = unique(DB_clean[,3:8])
summary_clean <- merge(tax, summary_clean, by = "Species", all = T)
write.table(summary_clean, file = paste0("Summary_file_",name,"_database.txt"),sep="\t", row.names = F)
DB = DB_raw %>% subset(V1 %in% DB_clean$V1)
write.table(DB, file = paste0("Sequences_", name, "_curated.tax"), sep = "\t", row.names = F, col.names = F) #tax file
Sequences = readDNAStringSet(paste0("~/OneDrive - AZTI/DB-Cristina/scripts/output_files/Sequences_",name,".fasta"))
Sequences = Sequences[c(which(names(Sequences) %in% DB_clean$V1))]
writeXStringSet(Sequences, filepath = paste0("Sequences_",name,"_curated.fasta"), append=FALSE, compress=FALSE, compression_level=NA, format="fasta") #fasta file
}

#This script will automatically delete sequences tagged as erroneous from the final version of the reference database (i.e., DB_clean). For manual inspection of tagged sequences and removal, a list of erroenously tagged sequences is outputted
#Output files include: Potentially erroneous sequence list, problematic sequence list, summary file of cleaned reference database, tax file of cleaned reference database, fasta file of cleaned reference database