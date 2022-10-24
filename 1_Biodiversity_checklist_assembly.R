# REQUIRED LIBRARIES
library(readxl)
library(taxize)
library(Orcs)
library(tidyverse)
library(usefun)
library(RSQLite)
library(rfishbase) #if problems with rfishbase, try unlink(rfishbase:::db_dir())
library(DBI)

#########################################################################
###########################################################################
# CREATION OF THE FISH BIODIVERSITY CHECKLIST 

setwd("~/Escritorio/Cristina/Datos/2021/Pipeline/Step 1/")

# Import data from fishbase
Baltic_Sea <-species_by_ecosystem(ecosystem = "Baltic Sea", server = "fishbase")
Barents_Sea <- species_by_ecosystem(ecosystem = "Barents Sea", server = "fishbase")
Black_Sea <-species_by_ecosystem(ecosystem = "Black Sea", server = "fishbase")
Canary_Current <-species_by_ecosystem(ecosystem = "Canary Current", server = "fishbase")
Celtic_Biscay_Shelf <-species_by_ecosystem(ecosystem = "Celtic-Biscay Shelf", server = "fishbase")
Faroe_Plateau <-species_by_ecosystem(ecosystem = "Faroe Plateau", server = "fishbase")
Greenland_Sea <-species_by_ecosystem(ecosystem = "Greenland Sea", server = "fishbase")
Iberian_Coastal <-species_by_ecosystem(ecosystem = "Iberian Coastal", server = "fishbase")
Iceland_Shelf_and_Sea <-species_by_ecosystem(ecosystem = "Iceland Shelf and Sea", server = "fishbase")
North_Sea <-species_by_ecosystem(ecosystem = "North Sea", server = "fishbase")
Norwegian_Sea <-species_by_ecosystem(ecosystem = "Norwegian Sea", server = "fishbase")
Mediterranean_Sea<-species_by_ecosystem(ecosystem = "Mediterranean Sea", server = "fishbase")

# Merge data
df <- list(data.frame(Baltic_Sea), data.frame(Barents_Sea), data.frame(Black_Sea), data.frame(Canary_Current),
            data.frame(Celtic_Biscay_Shelf),data.frame(Faroe_Plateau), data.frame(Greenland_Sea), data.frame(Iberian_Coastal),
            data.frame(Iceland_Shelf_and_Sea), data.frame(North_Sea), data.frame(Norwegian_Sea), data.frame(Mediterranean_Sea))
species_fishbase <- merge(df, by = "Species", all = TRUE)
species_fishbase <- data.frame(unique(species_fishbase$Species))
names(species_fishbase)<-"Species"

# Remove subspecies
species_fishbase <- species_fishbase %>% 
  separate(col=Species, into = c("genus","species","subspecies"), sep = " ") %>%
  unite(col = "Species", genus, species, sep = " ") %>% dplyr::select(Species)
species_fishbase <- unique(species_fishbase)

# Retrieve taxonomy from WORMS 
taxonomy <- classification(species_fishbase$Species, db="worms", marine_only=FALSE)

#Turn classification object into dataframe
nombre <- c()
output <- array(dim = c(dim(species_fishbase)[1],7))
colnames(output) <- c("Kingdom", "Phylum", "Class_Gigaclass", "Order", "Family", "Genus", "Species")
for (i in seq(1,length(taxonomy))) {
  nombre <- as.character(names(taxonomy[i]))
  if (!is.na(taxonomy[nombre])) { 
    if (sum(as.integer(eval(parse(text = paste0("taxonomy$`",nombre,"`$rank == 'Class'"))))) == 1) {
      output[i,] <- t(eval(parse(text = paste0("taxonomy$`",nombre, "`$name[taxonomy$`", nombre, "`$rank == 'Kingdom' |",
                                               "taxonomy$`", nombre, "`$rank == 'Phylum' |",
                                               "taxonomy$`", nombre, "`$rank == 'Class' |", 
                                               "taxonomy$`", nombre, "`$rank == 'Order' |",
                                               "taxonomy$`", nombre, "`$rank == 'Family' |",
                                               "taxonomy$`", nombre, "`$rank == 'Genus' |", 
                                               "taxonomy$`", nombre, "`$rank == 'Species']" ))))
    } else {
      output[i,] <- t(eval(parse(text = paste0("taxonomy$`",nombre, "`$name[taxonomy$`", nombre, "`$rank == 'Kingdom' |",
                                               "taxonomy$`", nombre, "`$rank == 'Phylum' |",
                                               "taxonomy$`", nombre, "`$rank == 'Gigaclass' |",
                                               "taxonomy$`", nombre, "`$rank == 'Order' |",
                                               "taxonomy$`", nombre, "`$rank == 'Family' |",
                                               "taxonomy$`", nombre, "`$rank == 'Genus' |",
                                               "taxonomy$`", nombre, "`$rank == 'Species']" ))))
    }
    
  }
  rm(nombre)
}
output <- data.frame(output)
order_clean <- data.frame(gsub("\\ .*","",output$Order))
output<-cbind(output, order_clean)
output<-output[, c(1,2,3,8,4,5,6,7)]
output$Order<-NULL
names(output)[4]<- "Order"
output<-output[complete.cases(output),]# LIST OF MARINE EUROPEAN FISH SPECIES

# Export data
write.table(output, "Species_list.txt", sep = " ",
            row.names = FALSE, col.names = TRUE)

# Create species list to retrieve sequence information from NCBI 
sp_names <- data.frame(output$Species)
write.table(sp_names, "sp_names.txt",
            row.names = FALSE, col.names = FALSE)
