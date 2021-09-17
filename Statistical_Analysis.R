library(reshape2)
library(data.table)
# Reading data

csv_data <- read.csv(file = 'FDs_FirstSecond.csv') # read from working directory

# Data extraction and transpose
Brain_experts<- c(5, 6, 7, 9, 10, 12, 18, 19, 24, 26, 28, 30, 32, 35, 39)
Other_experts<- c(11, 15, 20, 21, 25, 27, 29, 33, 34, 36, 37)  
Registrars<- c(13, 14, 17, 23, 31)
Naives<-40:65
All_radiologists<-c(Brain_experts,Other_experts,Registrars)

Group_A<-Brain_experts # need to change group to get analysis among different groups
Group_B<-Naives

Pathological_Stimuli<- 2:21
Normal_Stimuli<- 22:41

extracted_data <- abs(csv_data[Pathological_Stimuli,c(Group_A, Group_B)]) # need to change stimuli type to get analysis for pathological and normal stimuli separately 

# Long format conversion

long_data <- melt(setDT(extracted_data), measure.vars = 1:ncol(extracted_data), variable.name = "Subject",value.name = "FD")
#print(long_data)
#summary(long)

# Add a new column
Final_long_data <- cbind(long_data, Group = rep(c("Brain-experts","Registrars"), times = c(length(Group_A)*20,length(Group_B)*20)))

#############
### Linear model  ON  the MEANS BY SUBJECT 
#############
### get the new structure of the data meaning mean of FD for each subject
#########################################
library(dplyr)
library(effectsize)
new.data <- Final_long_data %>%
  group_by(Subject) %>%
  summarize(mean_FD = mean(FD, na.rm = TRUE),Group=unique(Group))
##################

#### To get the mean by group and standard error by group
sd.data <- new.data %>%
  group_by(Group)  %>%
  summarize(FD.mean = mean(mean_FD, na.rm = TRUE),Sd.error=sd(mean_FD)/sqrt(n()))
sd.data
#################

effectsize::cohens_d(mean_FD ~ Group, data = new.data)   ###  use this one which is the pooled version (classical one)

### degree freedom 
lm.model <- lm(mean_FD~Group,data=new.data)
res.analysis <- summary(lm.model)
DF <- res.analysis$df[2]   
Pvalue <- res.analysis$coefficients[2,4]
tvalue <- abs(res.analysis$coefficients[2,3])
DF
Pvalue
tvalue
#summary(lm.model)  ## for double check 


### pvalue correction
pvalues <- c(0.4855,0.0826,0.7827,0.8333,0.2971) # p values are corrected after getting values for all tests
p.adjust(pvalues,method="fdr")


