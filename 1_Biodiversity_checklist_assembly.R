###########################################################################

# REQUIRED LIBRARIES
library(readxl)
library(taxize)
library(Orcs)
library(tidyverse)
library(usefun)

#########################################################################
###########################################################################
# CREATION OF THE FISH BIODIVERSITY CHECKLIST 
setwd()
setwd("~/Escritorio/Cristina/Datos/2021/")

#Load data
Baltic_Sea <- read_excel("FishBase_LME.xlsx", sheet = "Baltic Sea")
Barents_Sea <-read_excel("FishBase_LME.xlsx", sheet = "Barents Sea")
Black_Sea <-read_excel("FishBase_LME.xlsx", sheet = "Black Sea")
Canary_Current <-read_excel("FishBase_LME.xlsx", sheet = "Canary current")
Celtic_Biscay_Shelf <-read_excel("FishBase_LME.xlsx", sheet = "Celtic Biscay Shelf")
Faroe_Plateau <-read_excel("FishBase_LME.xlsx", sheet = "Faroe Plateau")
Greenland_Sea <-read_excel("FishBase_LME.xlsx", sheet = "Greenland Sea")
Iberian_Coastal <-read_excel("FishBase_LME.xlsx", sheet = "Iberian Coastal")
Iceland_Shelf_and_Sea <-read_excel("FishBase_LME.xlsx", sheet = "Iceland Shelf and Sea")
North_Sea <-read_excel("FishBase_LME.xlsx", sheet = "North Sea")
Norwegian_Sea <-read_excel("FishBase_LME.xlsx", sheet = "Norwegian Sea")
Mediterranean_Sea <-read_excel("FishBase_LME.xlsx", sheet = "Mediterranean Sea")

# Merge data
df <- list(data.frame(Baltic_Sea), data.frame(Barents_Sea), data.frame(Black_Sea), data.frame(Canary_Current),
            data.frame(Celtic_Biscay_Shelf),data.frame(Faroe_Plateau), data.frame(Greenland_Sea), data.frame(Iberian_Coastal),
            data.frame(Iceland_Shelf_and_Sea), data.frame(North_Sea), data.frame(Norwegian_Sea), data.frame(Mediterranean_Sea))
species_fishbase <- merge(df, by = "Species", all = TRUE)
species_fishbase <- species_fishbase$Species

# Retrieve taxonomy from WORMS 
taxonomy <- classification(species_fishbase, db="worms")
class(taxonomy)
#Warning: Not found 12 species*

#Turn classification object into dataframe
nombre <- c()
output <- array(dim = c(length(species_fishbase),7))
colnames(output) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
for (i in seq(1,length(taxonomy))) {
  nombre <- as.character(names(taxonomy[i]))
  if (!is.na(taxonomy[nombre])) {
    output[i,] <- t(eval(parse(text = paste0("taxonomy$`",nombre, "`$name[taxonomy$`", nombre, "`$rank == 'Kingdom' |",
                                             "taxonomy$`", nombre, "`$rank == 'Phylum' |",
                                             "taxonomy$`", nombre, "`$rank == 'Class' |", 
                                             "taxonomy$`", nombre, "`$rank == 'Order' |",
                                             "taxonomy$`", nombre, "`$rank == 'Family' |",
                                             "taxonomy$`", nombre, "`$rank == 'Genus' |", 
                                             "taxonomy$`", nombre, "`$rank == 'Species']" ))))
  }
  rm(nombre)
}
output <- data.frame(output)
head(output)

# Check those "Not found records"
nrow(distinct(as.data.frame(species_fishbase))) - nrow(distinct(output)) #input - output

found.sp <- distinct(output) 
found.sp <- found.sp$Species # species that were found in WORMS and whose taxonomy is compiled
notfound.sp <- as.data.frame(outersect(species_fishbase, found.sp))
names(notfound.sp) <- "Species"
head(notfound.sp) #there are some subspecies
notfound.sp <- separate(notfound.sp, Species, into = c("genus","species","subspecies"), sep =" ")
notfound.sp <- as.data.frame(paste(notfound.sp$genus,"",notfound.sp$species))
names(notfound.sp)<- "Species"
notfound.sp <- distinct(notfound.sp) #subspecies removed, now check again in worms
notfound.sp <- notfound.sp$Species
notfound.sp <- as.character(notfound.sp)
notfound.sp <- str_squish(notfound.sp)
taxonomy2 <- classification(notfound.sp, db="worms")
# Again, we get 12 not found species, which are not marine. e.g. Barbatula barbatula is freschwater.

# Check those found 8 species 
found2 <- map(taxonomy2,~data.frame(.))
found2 <- map_dfr(found2,~mutate_all(.,as.character)) 
found2 <- found2[found2$rank =="Species",]
found2 <- found2$name
found2<-found2[!is.na(found2)]
length(intersect(output$Species, found2))
# This species were found in the first classification() but in our species_fishbase file we had both subspecies
#and species e.g. Stomias boa, Stomias boa boa, Stomias boa ferox being multiplpe records for the same species.

output <- unique(output) # LIST OF MARINE EUROPEAN FISH SPECIES
write.csv(output, file = "Species_list.csv")
write.table(output, "Species_list.txt", sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
