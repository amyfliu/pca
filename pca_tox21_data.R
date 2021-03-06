#load the needed libraries
library(readxl)
library(kohonen)
library(tidyverse)
library(factoextra)

#read in the data
tox21 <- read.delim("tox21_data.txt") #8,971 distinct compounds/drugs

#Change rows with 'x' to NA and filter out drugs/chemicals with incomplete data using the complete.case()
new_df <- tox21
new_df[new_df == 'x'] <- NA  
final_data <- new_df[complete.cases(new_df),] 

#all columns are character, need to convert all columns to numerical for analysis
#data frames are not "character" or "numerical", but COLUMNS are...
final_data_2 <- as.data.frame(sapply(final_data[-c(1,2)], as.numeric)) 
#NOTE: column 1 and 2 is excluded

#rejoin the CAS ID and bioassay data into a data frame 
final_data_3 =
  data.frame(final_data$Structure.ID, final_data_2) %>% 
  janitor::clean_names() %>% 
  rename(structure_id = final_data_structure_id) %>% 
  mutate(num_zero = rowSums(. == 0)) %>% 
  filter(num_zero < 162) %>% 
  select(-num_zero)

########################## END OF DATA CLEANING ###############################

#PCA ANALSYSIS 
pca <- prcomp(final_data_3[,-1], scale = TRUE)

#plot pc1 and pc2
plot(pca$x[,1], pca$x[,2]) #NOT VERY HELPFUL... 

#make scree plot 
pca.var <- pca$sdev^2 #variance
pca.var.percent <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.percent, xlab = "Principal component", ylab = "% variation", main= "Scree plot")

pca.var.percent[1:10]


#Use ggplot and make a fancier plot 
pca.data <- data.frame(Sample=final_data_3$structure_id,
                       X=pca$x[,1], #gives us the pca scores (i.e, the linear combinations)
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 (", pca.var.percent[1], "%)", sep="")) +
  ylab(paste("PC2 (", pca.var.percent[2], "%)", sep="")) + 
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 bioassay points that contribute most to PC1.
loading_scores <- pca$rotation[,1]
assay_scores <- abs(loading_scores) ## get the magnitudes
assay_score_ranked <- sort(assay_scores, decreasing=TRUE)
top_20_assay <- names(assay_score_ranked[1:20]) #top 20

top_20_assay ## show the names of the top bioassay that dominate the 1st principle component 
pca$rotation[top_20_assay, 1] #what 'variable/bioassay' dominates the first principle component?? and what are the weights

#write the results to a file
write.csv(top_20_assay, file = "top10_assay_PC3.csv")



###### More exploring ####### 

#find which bioassay points can be 'dropped'
tox <- read_excel("~/Documents/Axle/clustering/SOM/axle_clustering/data/tox21_data.xlsx")
tox[tox == 'x'] <- NA  #convert all missing x's to "NA"

sum(is.na(tox$Structure_ID))
sum(is.na(tox$`tox21-ache-p1_ratio`) | tox$`tox21-ache-p1_ratio`== 0)
sum(is.na(tox$`tox21-ache-p3_ratio`) | tox$`tox21-ache-p3_ratio`== 0)

#my_vec <- vector('num_missing'); #create empty vector 
#my_vec <- c("num_missing") #character vector 
my_vec <- vector()
for(i in 2:ncol(tox)) {
  new_value <- sum(is.na(tox[i]) | tox[i] == 0) #gives # of null or 0's
  my_vec <- c(my_vec, new_value) #use c() to combine the vectors
}
my_vec

#new_tox <- rbind(tox, my_vec) #adds new row with the # of missing/0's
#tail(new_tox)

my_vec_2 <- my_vec/8971
my_vec_2

test3 <- sort(my_vec_2)
test3

#check how many 0 or missing data there are
sum(test3 > 0.75) #gives us the % of missing or null values 
sum(test3 > 0.80)
sum(test3 > 0.85)
sum(test3 > 0.90)
