# Data Wrangling

learning <- read.csv("Data_Cortex_Nuclear.csv", TRUE, ",", stringsAsFactors=FALSE)
library(tidyr)
library(dplyr)
library(stringr)

# Add class_number column
learning1 <- learning %>% mutate(class_number = class)

# Assign Classes to #1-8
learning1$class_number[learning1$class_number == "c-CS-m"] <- "1"
learning1$class_number[learning1$class_number == "c-SC-m"] <- "2"
learning1$class_number[learning1$class_number == "c-CS-s"] <- "3"
learning1$class_number[learning1$class_number == "c-SC-s"] <- "4"
learning1$class_number[learning1$class_number == "t-CS-m"] <- "5"
learning1$class_number[learning1$class_number == "t-SC-m"] <- "6"
learning1$class_number[learning1$class_number == "t-CS-s"] <- "7"
learning1$class_number[learning1$class_number == "t-SC-s"] <- "8"

# Add learning column
learning1 <- learning1 %>% mutate (learning = class_number)

# Assign 1 = learn 0 = do not learn
learning1$learning[learning1$learning == "1"] <- "1"
learning1$learning[learning1$learning == "3"] <- "1"
learning1$learning[learning1$learning == "5"] <- "1"
learning1$learning[learning1$learning == "2"] <- "0"
learning1$learning[learning1$learning == "4"] <- "0"
learning1$learning[learning1$learning == "6"] <- "0"
learning1$learning[learning1$learning == "7"] <- "0"
learning1$learning[learning1$learning == "8"] <- "0"

# replace NA values with mean of the corresponding class/class_number
learning2 <- learning1 %>% group_by(class) %>%
  mutate_each(funs(replace(., which(is.na(.)), mean(., na.rm=TRUE))), ends_with('_N'))

write.csv(learning2, "learning_capstone_DW.csv")

# Statistical Analyses

learning_stats <- read.csv("learning_capstone_DW.csv", TRUE, ",", stringsAsFactors=FALSE)
library(ggplot2)
library(dplyr)

# Determine the mean protein levels for each class
learning_prot_mean <- aggregate(learning_stats[, 3:79], list(learning_stats$class_number), mean)
learning_prot_mean

# Use transpose() to switch rows and columns
learning_prot_mean_trans <- as.data.frame(t(as.matrix(learning_prot_mean)))
learning_prot_mean_trans

# Remove first row labeled 'Group.1' which is the class number not a protein
learning_prot_mean_trans2 <- learning_prot_mean_trans[-c(1),]
learning_prot_mean_trans2

# Rename column headings to class description
learning_prot_mean_trans3 <- learning_prot_mean_trans2 %>%
  rename(
    class1 = V1,
    class2 = V2,
    class3 = V3,
    class4 = V4,
    class5 = V5,
    class6 = V6,
    class7 = V7,
    class8 = V8
  )

# plot relationships of mean protein levels between each class in scatter plot matrices
# points outside of linear relationship may be result of experimental class
pairs(learning_prot_mean_trans3, upper.panel = panel.smooth, lower.panel = NULL)

# generate correlation matrix of mean protein levels between each class

install.packages("Hmisc")
library(Hmisc)

flattenCorMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

learning_prot_mean_matrix <- rcorr(as.matrix(learning_prot_mean[, 2:78]))
learning_prot_matrix2 <- flattenCorMatrix(learning_prot_mean_matrix$r, learning_prot_mean_matrix$P)

# extract combinations of proteins that have correlation >= 0.95
learning_prot_matrix2[learning_prot_matrix2$cor>=0.95, ]

# note: 14 protein combinations meet the >=0.95 criteria

# remove 9 proteins with high correlation
proteins_removed <- c('DYRK1A_N', 'ITSN1_N', 'pERK_N', 'NUMB_N', 'CDK5_N', 'pS6_N', 'pPKCAB_N', 'pAKT_N', 'NR2B_N')
learning_stats2 <- learning_stats[ , !(names(learning_stats) %in% proteins_removed)]
learning_prot_mean_removed <- learning_prot_mean[ , !(names(learning_prot_mean) %in% proteins_removed)]

# Were 9 columns removed?
str(learning_stats)
str(learning_stats2)
str(learning_prot_mean)
str(learning_prot_mean_removed)

# Yes, they were removed

# plot correlations of mean protein levels
install.packages("corrplot")
library(corrplot)
correlations <- cor(learning_protein_mean_removed[, 2:69])
corrplot(correlations, method="circle", tl.cex = 0.6)

# remove all columns except class_number and each protein tested
columns_removed <- c('X', 'MouseID', 'Genotype', 'Treatment', 'Behavior', 'class', 'learning')
learning_stats3 <- learning_stats2[ , !(names(learning_stats2) %in% columns_removed)]

# move class_number variable from last to first
column_reorder <- grep("class_number", names(learning_stats3))
learning_stats4 <- learning_stats3[, c(column_reorder, (1:ncol(learning_stats3))[-column_reorder])]
head(learning_stats4)

# perform repeated measures ANOVA to see if the protein levels between classes 1-8 are different

learning_stats4$class_number <- as.factor(learning_stats4$class_number)
i= 3

for (i in 2:dim(learning_stats4)[2]) {
  fm1 <- aov(learning_stats4[,i]~learning_stats4$class_number)
  summary(fm1)
  print(paste(i, colnames(learning_stats4)[i]))
  print(summary(fm1))
  boxplot(learning_stats4[,i] ~ learning_stats4$class_number,
          col=rainbow(length(unique(learning_stats4$class_number)), alpha=0.3),
          medcol = rainbow(length(unique(learning_stats4$class_number)), v=0.7),
          whiskcol = rainbow(length(unique(learning_stats4$class_number))),
          staplecol = rainbow(length(unique(learning_stats4$class_number)), v=0.7),
          boxcol = rainbow(length(unique(learning_stats4$class_number)), v=0.7),
          outcol = rainbow(length(unique(learning_stats4$class_number)), v=0.7),
          ylab = paste("Protein", colnames(learning_stats4)[i]),
          xlab = "class")
}

# Change colors for representative plots that show potential trends with mean protein levels based on classes
# change colors for APP_N classes by genotype WT=red mutant=blue
boxplot(learning_stats4[,27] ~ learning_stats4$class_number,
        col=alpha(c("red", "red", "red", "red", "blue", "blue", "blue", "blue"),0.3),
        medcol = c("red", "red", "red", "red", "blue", "blue", "blue", "blue"),
        whiskcol = c("red", "red", "red", "red", "blue", "blue", "blue", "blue"),
        staplecol = c("red", "red", "red", "red", "blue", "blue", "blue", "blue"),
        boxcol = c("red", "red", "red", "red", "blue", "blue", "blue", "blue"),
        outcol = c("red", "red", "red", "red", "blue", "blue", "blue", "blue"),
        ylab = paste("Protein", colnames(learning_stats4)[27]),
        xlab = "class")

legend("topright", c("Wild-type", "Mutant"), lty=1, col=c("red", "blue"), bty='n', cex=0.8, title="Genotype", text.font = 4)

# change colors for SOD1_N classes by behavior test given CS=red SC=blue
boxplot(learning_stats4[,29] ~ learning_stats4$class_number,
        col=alpha(c("red", "blue", "red", "blue", "red", "blue", "red", "blue"),0.3),
        medcol = c("red", "blue", "red", "blue", "red", "blue", "red", "blue"),
        whiskcol = c("red", "blue", "red", "blue", "red", "blue", "red", "blue"),
        staplecol = c("red", "blue", "red", "blue", "red", "blue", "red", "blue"),
        boxcol = c("red", "blue", "red", "blue", "red", "blue", "red", "blue"),
        outcol = c("red", "blue", "red", "blue", "red", "blue", "red", "blue"),
        ylab = paste("Protein", colnames(learning_stats4)[29]),
        xlab = "class")

legend("topright", c("CS", "SC"), lty=1, col=c("red", "blue"), bty='n', cex=0.8, title="Behavior test", text.font = 4)

# change colors for pPKCG_N classes by drug given drug=red no drug=blue
boxplot(learning_stats4[,41] ~ learning_stats4$class_number,
        col=alpha(c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),0.3),
        medcol = c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),
        whiskcol = c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),
        staplecol = c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),
        boxcol = c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),
        outcol = c("red", "red", "blue", "blue", "red", "red", "blue", "blue"),
        ylab = paste("Protein", colnames(learning_stats4)[41]),
        xlab = "class")

legend("topright", c("drug", "no drug"), lty=1, col=c("red", "blue"), bty='n', cex=0.8, title="Drug given", text.font = 4)

# change colors for GluR4_N classes so all have same color
boxplot(learning_stats4[,53] ~ learning_stats4$class_number,
        col=alpha(c("purple", "purple", "purple", "purple", "purple", "purple", "purple", "purple"),0.3),
        medcol = c("purple", "purple", "purple", "purple", "purple", "purple", "purple", "purple"),
        whiskcol = c("purple", "purple", "purple", "purple", "purple", "purple", "purple", "purple"),
        staplecol = c("purple", "purple", "purple", "purple", "purple", "purple", "purple", "purple"),
        boxcol = c("purple", "purple", "purple", "purple", "purple", "purple", "purple", "purple"),
        outcol = c("purple", "purple", "purple", "purple", "purple", "purple", "purple", "purple"),
        ylab = paste("Protein", colnames(learning_stats4)[53]),
        xlab = "class")


# Identifying the best 2 parameter model
temp = expand.grid(a = colnames(learning_stats4)[2:69],
                   b = colnames(learning_stats4)[2:69])

# Remove models that use the same protein. ex. SHH_N + SHH_N
temp = temp[which(temp$a != temp$b),]

# Create an AIC column to store the scoring
# AIC measures goodness-of-fit of the model
# the lower the AIC value, the better the model is
temp$AIC = 999999999
temp$coef = 999999999

# create an empty data frame to put the final data set into
out_final = vector(length = 0)

## Looping through each class - fitting a logistic regression
for(k in 1:length(unique(learning_stats4$class_number))){
  print(paste("we are on class", k))
  
  ## Take a temp dataset to fill out
  out = temp
  
  ## dummy variable coding
  # First making a vector of all 0s
  class_number1 = rep(0, dim(learning_stats4)[1])
  
  # one hot encoding - change the 0s to 1s if the class is equal to k, where k is 1 through 8
  class_number1[which(learning_stats4$class_number==k)]=1
  out$class = k
  
  ## Testing the 2 possible logistic regression models for each dummy variable coded classes
  for(i in 1:dim(out)[1]){
    print(paste(i, as.character(out[i,1]), as.character(out[i,2])) )
    model = glm(class_number1 ~ learning_stats4[,as.character(out[i,1])] +
                  learning_stats4[,as.character(out[i,2])],
                family = binomial(link='logit') )
    
    ## fill in empty vector with AIC and coefficient values
    out$AIC[i] = model$aic
    rm(model)
  }
  out_final = rbind(out_final, out)
}

# Find lowest AIC for each class
library(data.table)

DT <- data.table(out_final)
DT[ , .SD[which.min(AIC)], by = class]

# Run models with the lowest AIC for each class

###class 1
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==1)]=1
model1 = glm(class_number1 ~ learning_stats4[,"SOD1_N"] + learning_stats4[,"MTOR_N"],  family = binomial(link='logit'))
predicted_percent = round(predict(model1, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy
# accuracy = 0.887

## class 2
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==2)]=1

model2 = glm(class_number1 ~ learning_stats4[,"APP_N"] + learning_stats4[,"pCAMKII_N"],  family = binomial(link='logit'))

predicted_percent = round(predict(model2, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy
# accuracy = 0.944

## class 3
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==3)]=1

model3 = glm(class_number1 ~ learning_stats4[,"Ubiquitin_N"] + learning_stats4[,"pGSK3B_N"],  family = binomial(link='logit'))

predicted_percent = round(predict(model3, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy

# accuracy = 0.877

## class 4
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==4)]=1

model4 = glm(class_number1 ~ learning_stats4[,"pNUMB_N"] + learning_stats4[,"MTOR_N"],  family = binomial(link='logit'))

predicted_percent = round(predict(model4, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy

# accuracy = 0.957

## class 5
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==5)]=1

model5 = glm(class_number1 ~ learning_stats4[,"SOD1_N"] + learning_stats4[,"APP_N"],  family = binomial(link='logit'))

predicted_percent = round(predict(model5, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy

# accuracy = 0.877

## class 6
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==6)]=1

model6 = glm(class_number1 ~ learning_stats4[,"CaNA_N"] + learning_stats4[,"ERBB4_N"],  family = binomial(link='logit'))

predicted_percent = round(predict(model6, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy

# accuracy = 0.888

## class 7
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==7)]=1

model7 = glm(class_number1 ~ learning_stats4[,"GluR3_N"] + learning_stats4[,"P38_N"],  family = binomial(link='logit'), maxit = 100)

predicted_percent = round(predict(model7, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy

# accuracy = 0.924

## class 8
class_number1 = rep(0, dim(learning_stats4)[1])
class_number1[which(learning_stats4$class_number==8)]=1

model8 = glm(class_number1 ~ learning_stats4[,"pPKCG_N"] + learning_stats4[,"BRAF_N"],  family = binomial(link='logit'), maxit = 100)

predicted_percent = round(predict(model8, type="response"),2)
out1 = predicted_percent
out1[out1>=0.5] = 1
out1[out1<0.5] = 0

y_output = data.frame(predicted_percent = predicted_percent,
                      predicted = round(out1,1),
                      observed = class_number1)

dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])[1] ## true positive
dim(y_output[which(y_output$predicted == 1 & y_output$observed == 0),])[1] ## false positive
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 1),])[1] ## false negative
dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1] ## true negative
dim(y_output[which(y_output$predicted != 0 & y_output$predicted != 1),])[1]

plot(x = predicted_percent,
     y = class_number1,
     ylab = "observed",
     xlab = "predicted")

accuracy <- ((dim(y_output[which(y_output$predicted == 1 & y_output$observed == 1),])) + (dim(y_output[which(y_output$predicted == 0 & y_output$observed == 0),])[1]))/1080
accuracy

# accuracy = 0.908


## Identifying proteins whose levels differ in mutant mice given the learning behavior with drug or no drug

# subset data to include only classes 5 and 7 

learning_stats5 = learning_stats4[learning_stats4$class_number %in% c("5", "7"), ]


learning_stats5$class_number <- as.factor(learning_stats5$class_number)
i= 3


# Determine if there are significant differences in protein levels between classes 5 and 7
t.test_results <- t(sapply(learning_stats5[-1], function(x)
  unlist(t.test(x~learning_stats5$class_number)[c("estimate", "p.value", "statistic", "conf.int")])))

t.test_results2 <- as.data.frame(t.test_results)

str(t.test_results2)

# extract proteins that display p value <=0.05
t.test_results2_sig <- t.test_results2[t.test_results2$p.value<=0.05, ]

str(t.test_results2_sig)

# 36 proteins display significant difference between class 5 (drug) and class 7 (no drug)

# subset data to include only classes 1 and 3

learning_stats6 = learning_stats4[learning_stats4$class_number %in% c("1", "3"), ]


learning_stats6$class_number <- as.factor(learning_stats6$class_number)
i= 3

# Determine if there are significant differences in protein levels between classes 1 and 3
t.test_resultsb <- t(sapply(learning_stats6[-1], function(x)
  unlist(t.test(x~learning_stats6$class_number)[c("protein", "estimate", "p.value", "statistic", "conf.int")])))

t.test_resultsb2 <- as.data.frame(t.test_resultsb)

str(t.test_resultsb2)

# extract proteins that display p value <=0.05
t.test_resultsb2_sig <- t.test_resultsb2[t.test_resultsb2$p.value<=0.05, ]

str(t.test_resultsb2_sig)

# 30 proteins display significant difference between class 1 and class 3

# From these 2 lists of proteins, 19 are unique to being different in class 5 and 7

# subset data to include only classes 1, 3, 5, and 7

learning_stats7 = learning_stats4[learning_stats4$class_number %in% c("1", "3", "5", "7"), ]

learning_stats7$class_number <- as.factor(learning_stats7$class_number)
i= 3

# create box plot for each protein only for classes 1, 3, 5, and 7
for (i in 2:dim(learning_stats7)[2]) {
  fm4 <- aov(learning_stats7[,i]~learning_stats7$class_number)
  summary(fm4)
  print(paste(i, colnames(learning_stats7)[i]))
  boxplot(learning_stats7[,i] ~ learning_stats7$class_number,
          col=alpha(c("red", "red", "orange", "red", "blue", "blue", "purple", "blue"),0.4),
          medcol = c("red", "red", "orange", "red", "blue", "blue", "purple", "blue"),
          whiskcol = c("red", "red", "orange", "red", "blue", "blue", "purple", "blue"),
          staplecol = c("red", "red", "orange", "red", "blue", "blue", "purple", "blue"),
          boxcol = c("red", "red", "orange", "red", "blue", "blue", "purple", "blue"),
          outcol = c("red", "red", "orange", "red", "blue", "blue", "purple", "blue"),
          ylab = paste("Protein", colnames(learning_stats7)[i]),
          xlab = "class")
}
legend("topright", c("WT-drug", "WT-no drug", "Ts65Dn-drug", "Ts65Dn-no drug"), lty=1, col=c("red", "orange", "blue", "purple"), bty='n', cex=0.8, title="Genotype and Drug given", text.font = 4)
