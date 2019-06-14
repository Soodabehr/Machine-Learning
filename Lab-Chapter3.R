

#1.Load data in peptide_peak30.csv, and peptide_peak40.csv, into 2 separate data frames in R.
peptide_peak30<- data.frame(read.csv("C:/Users/soodr/Documents/CS5310-Data Mining/Labs/Chapter 3/peptide_peak30.csv", stringsAsFactors = FALSE))
peptide_peak40<- data.frame(read.csv("C:/Users/soodr/Documents/CS5310-Data Mining/Labs/Chapter 3/peptide_peak40.csv", stringsAsFactors = FALSE))


#2.Process each dataset to remove the columns we don't need, 
 peptide_peak30$Fraction<- NULL
peptide_peak30$Source.File<-NULL
peptide_peak30$X.Spec.Sample.7<-NULL
peptide_peak30$PTM<- NULL
peptide_peak30$AScore<-NULL


peptide_peak40$Fraction<- NULL
peptide_peak40$Source.File<-NULL
peptide_peak40$X.Spec.Sample.9<-NULL
peptide_peak40$PTM<- NULL
peptide_peak40$AScore<-NULL



library(tidyverse)
library(tidyr)

#Split the 1/K0 Range column into 3 columns
peptide_peak30_R <- peptide_peak30 %>% separate(X1.K0.Range,c('StartingRange','EndingRange'), sep = "-")
peptide_peak30_R$Range <- as.numeric(peptide_peak30_R$EndingRange)-as.numeric(peptide_peak30_R$StartingRange)

peptide_peak40_R <- peptide_peak40 %>% separate(X1.K0.Range,c('StartingRange','EndingRange'), sep = "-")
peptide_peak40_R$Range <- as.numeric(peptide_peak40_R$EndingRange)-as.numeric(peptide_peak40_R$StartingRange)


#Normalize the columns with numeric data
# create a function 
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


numericals_cols_peak30 <- subset(peptide_peak30_R,select = X.10lgP:Intensity.Sample.7)
numeric_df_peak30 <- mutate_all(numericals_cols_peak30, function(x) as.numeric(as.character(x)))
dfNormal_peak30 <- as.data.frame(lapply(numeric_df_peak30, normalize))
final_peptide_peak30 <- cbind(peptide_peak30_R$Peptide,dfNormal_peak30,subset(peptide_peak30_R,select = Precursor.Id:Range))
write.csv(final_peptide_peak30,"peptide_peak30_clean.csv")

numericals_cols_peak40 <- subset(peptide_peak40_R,select = X.10lgP:Intensity.Sample.9)
numeric_df_peak40 <- mutate_all(numericals_cols_peak40, function(x) as.numeric(as.character(x)))
dfNormal_peak40 <- as.data.frame(lapply(numeric_df_peak40, normalize))
final_peptide_peak40 <- cbind(peptide_peak40_R$Peptide,dfNormal_peak40,subset(peptide_peak40_R,select = Precursor.Id:Range))
write.csv(final_peptide_peak40,"peptide_peak40_clean.csv")

#3.Compare peptides in 2 datasets and create the target feature to indicate whether the peptide appears in both datasets or only one of them.
#library(tidyverse)
final_peptide_peak40$common<-final_peptide_peak40$peptide %in% final_peptide_peak30$peptide
write.csv(final_peptide_peak40,"peptide_peak40_clean.csv")

final_peptide_peak30$common<-final_peptide_peak30$peptide %in% final_peptide_peak40$peptide
write.csv(final_peptide_peak30,"peptide_peak30_clean.csv")


#Note--This Common variable will be used as the target feature. 

#4.Combine two datasets into one.
names(final_peptide_peak30)[names(final_peptide_peak40) == 'Intensity.Sample.7'] <- 'Intensity.Sample.9'
df_peptide_peak30=as.data.frame(final_peptide_peak30)
df_peptide_peak40=as.data.frame(final_peptide_peak40)
#install.packages("plyr")
library(plyr)
combind_dataset<- rbind.fill(final_peptide_peak30,final_peptide_peak40)
#combind_dataset<-(rbind(df_peptide_peak30,df_peptide_peak40))
numbers_dataset<-as.data.frame(combind_dataset[2:10])

write.csv(combind_dataset,"peptide_combined.csv")
#eighty_percent_data=round((80*length(combind_dataset$Peptide))/100)


#Randomly sample 80% of the data as the training dataset and the remaining 20% as the testing datase
n=dim(combind_dataset)[1]
v<-sample(1:n,0.8*n)


#Creating training and test datasets
df_peptide_train<-numbers_dataset[v,] 
df_peptide_test<-numbers_dataset[-v,]

#Store class lables in factor vectors
df_peptide_train_lables<- combind_dataset[v,15] 
df_peptide_test_lables<-combind_dataset[-v,15]

#5.Train the kNN model using the training dataset. Note that only numeric features should be used in the training process.
#install.packages("class")
library(class)
k=round(sqrt(dim(df_peptide_train)[1]))
peptide_test_pred<- knn(train=df_peptide_train[-9], test=df_peptide_test[-9], cl=df_peptide_train_lables,k)

#6&7.Use the kNN model to predict the target feature of the testing dataset.
#Use CrossTable() function to compare the predicted state values to the saved true values and analyze the results.
#install.packages("gmodels")
library(gmodels)
CrossTable(x=df_peptide_test_lables, y=peptide_test_pred,prop.chisq = FALSE)




