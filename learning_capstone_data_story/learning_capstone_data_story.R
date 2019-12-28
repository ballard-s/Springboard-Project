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

# remove all columns except class_number and each protein tested
columns_removed <- c('X', 'MouseID', 'Genotype', 'Treatment', 'Behavior', 'class', 'learning')
learning_stats3 <- learning_stats2[ , !(names(learning_stats2) %in% columns_removed)]

# move class_number column from last column to first column
column_reorder <- grep("class_number", names(learning_stats3))
learning_stats4 <- learning_stats3[, c(column_reorder, (1:ncol(learning_stats3))[-column_reorder])]
head(learning_stats4) # class_number has been moved to first column

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
