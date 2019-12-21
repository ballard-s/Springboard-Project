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

# Create new '.csv' file
write.csv(learning2, "learning_capstone_DW.csv")


