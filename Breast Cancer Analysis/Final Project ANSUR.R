# --- LIBRARY --- #
library(dplyr)
library(tidyverse)
library(ggplot2)
library(survival) 
library(ggsurvfit) 
library(KMsurv)


# --- READ DATA --- #
data <- read.delim("breast cancer data.tsv", header = TRUE, sep = "\t")
data
str(data)


# --- VARIABLES SELECTION --- #
data <- data %>% 
  select(Overall.Survival.Status,
         Overall.Survival..Months.,
         Diagnosis.Age,
         Race.Category,
         Subtype,
         Radiation.Therapy,
         Person.Neoplasm.Cancer.Status)
str(data)
names(data)


# --- CHECK MISSING VALUE --- #
sum(is.na(data))
data <- data %>% drop_na()
sum(is.na(data))


# --- CHECK DUPLICATE DATA --- #
duplicates <- duplicated(data)
sum(duplicates)
data <- data %>% distinct()
str(data)
summary(data)


# --- CHECK DATA TYPE --- #
unique(data$Race.Category)
data$Race.Category <- factor(data$Race.Category, 
                         levels = c("White", 
                                    "Black or African American", 
                                    "Asian", 
                                    "American Indian or Alaska Native"))
data$Race.Category <- as.numeric(data$Race.Category)

unique(data$Subtype)
data$Subtype <- factor(data$Subtype, 
                       levels = c("BRCA_LumA", 
                                  "BRCA_Her2", 
                                  "BRCA_LumB", 
                                  "BRCA_Normal", 
                                  "BRCA_Basal"))
data$Subtype <- as.numeric(data$Subtype)

unique(data$Overall.Survival.Status)
data$Overall.Survival.Status <- ifelse(data$Overall.Survival.Status == "0:LIVING", 0, 1)

unique(data$Radiation.Therapy)
data$Radiation.Therapy<- ifelse(data$Radiation.Therapy == "Yes", 1, 0)

unique(data$Person.Neoplasm.Cancer.Status)
data$Person.Neoplasm.Cancer.Status<- ifelse(data$Person.Neoplasm.Cancer.Status == "With Tumor", 1, 0)

str(data)


# --- BOXPLOT DIAGNOSIS AGE --- #
boxplot(data$Diagnosis.Age,  
        main = "Boxplot Diagnosis Age",  
        ylab = "Age",  
        col = "lightblue",  
        border = "darkblue",  
        horizontal = FALSE)


# --- BOXPLOT OVERALL SURVIVAL (MONTHS) ---#
boxplot(data$Overall.Survival..Months.,  
        main = "Boxplot Overall Survival (Months)",  
        ylab = "Duration (Months)",  
        col = "lightblue",  
        border = "darkblue",  
        horizontal = FALSE)


# --- PREPROCESSING DATA --- #
data$Survival_Status <- ifelse(data$Overall.Survival.Status == "1:DECEASED", 
                               1, 0)

data$Age_Category <- cut(data$Diagnosis.Age,  
                         breaks = c(-Inf, 40, 60, Inf),  
                         labels = c("<40", "40-60", ">60"))

data$Age_Category <- as.factor(data$Age_Category)

data <- data[data$Race.Category != "American Indian or Alaska Native", ] 
data$Race.Category <- as.factor(data$Race.Category) 
data$Subtype <- as.factor(data$Subtype) 
data$Radiation.Therapy <- as.factor(data$Radiation.Therapy) 
data$Person.Neoplasm.Cancer.Status <- as.factor(data$Person.Neoplasm.Cancer.Status)

str(data)


# --- CATEGORY LEVELS --- #
levels(data$Age_Category)
levels(data$Race.Category) 
data$Race.Category <- droplevels(data$Race.Category) 
levels(data$Race.Category)
levels(data$Subtype) 
levels(data$Radiation.Therapy)
levels(data$Person.Neoplasm.Cancer.Status) 


# --- DATA EXPLORATORY --- #
summary(data$Overall.Survival..Months.)
summary(data$Diagnosis.Age)


# --- BARPLOT FOR AGE CATEGORY --- #
age_freq <- table(data$Age_Category)
age_freq_sorted <- sort(age_freq, decreasing = FALSE) 
age_df <- as.data.frame(age_freq_sorted) 
ggplot(age_df, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", fill = "lightcoral") + 
  labs(title = "Distribusi Age Category", x = "Age Category", y = "Frekuensi"
  )


# --- BARPLOT FOR RACE CATEGORY --- #
race_freq <- table(data$Race.Category) 
race_freq_sorted <- sort(race_freq, decreasing = FALSE) 
race_df <- as.data.frame(race_freq_sorted) 
ggplot(race_df, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", fill = "lightgreen") + 
  labs(title = "Distribusi Race Category", x = "Race Category", y = "Frekuens
 i")


# --- BARPLOT FOR SUBTYPE --- #
sub_freq <- table(data$Subtype) 
sub_freq_sorted <- sort(sub_freq, decreasing = FALSE) 
sub_df <- as.data.frame(sub_freq_sorted) 
ggplot(sub_df, aes(x = Var1, y = Freq)) + 
  geom_bar(stat = "identity", fill = "lightblue") + 
  labs(title = "Distribusi Subtype", x = "Subtype", y = "Frekuensi")


# --- PIE CHART FOR RADIATION THERAPY --- #
radiation_df <- as.data.frame(radiation_freq) 
radiation_df$percentage <- 100 * radiation_df$Freq / sum(radiation_df$Freq) 
radiation_df$label <- paste0(round(radiation_df$percentage, 1), "%") 
ggplot(radiation_df, aes(x = "", y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", width = 1) +  
  coord_polar(theta = "y") +  
  labs(title = "Proporsi Radiation Therapy") + 
  scale_fill_manual(values = c("tomato", "gold")) + 
  theme_void() +  
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color 
            = "white", size = 5)  


# PIE CHART FOR CANCER STATUS --- #
neo_freq <- table(data$Person.Neoplasm.Cancer.Status) 
neo_df <- as.data.frame(neo_freq) 
neo_df$percentage <- 100 * neo_df$Freq / sum(neo_df$Freq) 
neo_df$label <- paste0(round(neo_df$percentage, 1), "%") 
ggplot(neo_df, aes(x = "", y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity", width = 1) +  
  coord_polar(theta = "y") +  
  labs(title = "Proporsi Cancer Status") + 
  scale_fill_manual(values = c("pink", "skyblue")) +  
  theme_void() 
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color 
            = "white", size = 5)  


# --- KAPLAN-MEIER ANALYSIS --- #
Y <- Surv(time = data$Overall.Survival..Months., event = data$Survival_Status)

kmfit1 <- survfit(Y~1) 
kmfit1

plot(kmfit1, 
     xlab = "Time (day)",  
     ylab = "Survival Probabilities",  
     main = "Kurva Kaplan-Meier")


# --- LOG RANK TEST & KM CURVES --- #
# AGE VARIABLE
lr1 <- survdiff(Y ~ data$Age_Category) ; lr1

plot(survfit(Y ~ data$Age_Category), col = c("red", "purple", "brown"), lty = 
       c("solid", "dashed", "dotted"), 
     xlab = "Time", 
     ylab = "Survival Probabilities", 
     main = "Kurva Kaplan-Meier\nusia") 
legend("bottomleft", c("<40 Tahun", "40-60 Tahun", ">60 Tahun"),  
       lty = c("solid", "dashed", "dotted"), 
       col = c("red", "purple", "brown"))

# RACE CATEGORY
lr2 <- survdiff(Y ~ data$Race.Category) ; lr2
plot(survfit(Y ~ data$Race.Category), col = c( "green", "blue", "red"),  
     lty = c("solid", "dashed", "dotted", "dotdash"), 
     xlab = "Time", 
     ylab = "Survival Probabilities", 
     main = "Kurva Kaplan-Meier\nRas") 
legend("bottomleft", levels(data$Race.Category),  
       lty = c("solid", "dashed", "dotted", "dotdash"), 
       col = c( "green", "blue", "red")) 

# SUBTYPE
lr3 <- survdiff(Y ~ data$Subtype) ; lr3
plot(survfit(Y ~ data$Subtype), col = c("darkblue", "gold", "purple", "green"
                                        ,"red"),  
     lty = c("solid", "dashed", "dotted", "dotdash", "twodash"), 
     xlab = "Time", 
     ylab = "Survival Probabilities", 
     main = "Kurva Kaplan-Meier\nSubtipe Kanker") 
legend("bottomleft", levels(data$Subtype),  
       lty = c("solid", "dashed", "dotted", "dotdash", "twodash"), 
       col = c("darkblue", "gold", "purple", "green", "red"))

# RADIATION THERAPY
lr4 <- survdiff(Y ~ data$Radiation.Therapy); lr4
plot(survfit(Y ~ data$Radiation.Therapy), col = c("blue", "magenta"),  
     lty = c("solid", "dashed"), 
     xlab = "Time", 
     ylab = "Survival Probabilities", 
     main = "Kurva Kaplan-Meier\nTerapi Radioterapi") 
legend("bottomleft", levels(data$Radiation.Therapy),  
       lty = c("solid", "dashed"), 
       col = c("blue", "magenta"))

#NEOPLASM CANCER STATUS
lr5 <- survdiff(Y ~ data$Person.Neoplasm.Cancer.Status) ; lr5
plot(survfit(Y ~ data$Person.Neoplasm.Cancer.Status), col = c("orange", "dark
 green"),  
     lty = c("solid", "dashed"), 
     xlab = "Time", 
     ylab = "Survival Probabilities", 
     main = "Kurva Kaplan-Meier\nStatus Kanker") 
legend("bottomleft", levels(data$Person.Neoplasm.Cancer.Status),  
       lty = c("solid", "dashed"), 
       col = c("orange", "darkgreen"))


# --- COX-PH MODEL --- #
Y <- Surv(data$Overall.Survival..Months., data$Survival_Status) 
cox_model <- coxph(Y ~ data$Age_Category + data$Race.Category + data$Subtype 
                   + data$Radiation.Therapy + data$Person.Neoplasm.Cancer.Status) 
cox_model 

cox_model2 <- coxph(Y ~ data$Radiation.Therapy + data$Person.Neoplasm.Cancer.Status) 
cox_model2
summary(cox_model2)
cox_model2$loglik


# --- PH ASSUMPTION TEST --- #
minusloglog <- function(x) -log(-log(x)) 

kmfit_radiation <- survfit(Y ~ data$Radiation.Therapy) 
plot(kmfit_radiation, fun = minusloglog, col = c("#b83143", "#1483A0"),  
     lwd = 1.8, xlab = "Waktu dalam Skala Log",  
     ylab = "Log-Log Survival", main = "Kurva Log-Log\nRadioterapi") 
legend("topright", c("Yes", "No"),  
       lty = 1:1, lwd = 2, col = c("#b83143", "#1483A0"), cex = 0.8)

kmfit_cancer_status <- survfit(Y ~ data$Person.Neoplasm.Cancer.Status) 
plot(kmfit_cancer_status, fun = minusloglog, col = c("#b83143", "#1483A0"),  
     lwd = 1.8, xlab = "Waktu dalam Skala Log",  
     ylab = "Log-Log Survival", main = "Kurva Log-Log\nStatus Tumor") 
legend("topright", c("With Tumor", "Tumor Free"),  
       lty = 1:1, lwd = 2, col = c("#b83143", "#1483A0"), cex = 0.8)


# --- PH ASSUMPTION TEST (OBS vs EXP) --- #
plot(kmfit_radiation, 
     main = "Observed vs Expected Survival Probability\nfor Radioterapi", 
     xlab = "Time (day)", ylab = "Survival Probability", 
     lty = 1, col = c("brown", "darkgreen")) 
log_kmfit_rad <- coxph(formula = Y ~ Radiation.Therapy , data = data) 
Radiation.Therapy2 <- data.frame(Radiation.Therapy = 0 : 1) 
lines(survfit(log_kmfit_rad, Radiation.Therapy2), col = c("brown", "darkgreen"), lty = 2)
legend("bottomleft", legend = c("0 (Obs)", "0 (Eks)", "1 (Obs)", "1 (Eks)"),  
       col = c("brown", "brown", "darkgreen", "darkgreen"), lty = c(1, 2, 1, 
                                                                    2), title = "Radioterapi", cex = 0.8) 
plot(kmfit_cancer_status, 
     main = "Observed vs Expected Survival Probability\nStatus Tumor", 
     xlab = "Time (day)", ylab = "Survival Probability", 
     lty = 1, col = c("brown", "darkgreen")) 
log_kmfit_can <- coxph(formula = Y ~ Person.Neoplasm.Cancer.Status , data = data) 
Person.Neoplasm.Cancer.Status2 <- data.frame(Person.Neoplasm.Cancer.Status = 
                                               0 : 1) 
lines(survfit(log_kmfit_can, Person.Neoplasm.Cancer.Status2), col = c("brown", "darkgreen"), lty = 2)
legend("bottomleft", legend = c("0 (Obs)", "0 (Eks)", "1 (Obs)", "1 (Eks)"),  
       col = c("brown", "brown", "darkgreen", "darkgreen"), lty = c(1, 2, 1, 
                                                                    2), title = "Status Tumor", cex = 0.8) 


# --- PH ASSUMPTION TEST (GOODNES OF FIT) --- #
check_ph <- cox.zph(cox_model2, transform = rank) 
check_ph


# --- BEST MODEL --- #
# 01: without any interaction
y <- Surv(time = data$Overall.Survival..Months., event = data$Survival_Status) 
model1 <- coxph(y~data$Radiation.Therapy+data$Person.Neoplasm.Cancer.Status) 
summary(model1)

# 02: with interaction
model2<-coxph(y~data$Radiation.Therapy+data$Person.Neoplasm.Cancer.Status+data$Radiation.Therapy*data$Person.Neoplasm.Cancer.Status) 
summary(model2)


# --- LRT TEST --- #
model2$loglik[2]
model1$loglik[2]

LR_stat2 <- -2*model1$loglik[2]-(-2*model2$loglik[2]) 
LR_stat2 

chisq_table2 <- qchisq(0.95,1) #chisq tabel 
pvalue <- 1-pchisq(LR_stat2,1) #pvalue 
pvalue

if(LR_stat2 > chisq_table2){ 
  print("Reject H0") 
  print("Conclusion: With 95% confidence level, model 1 is better.") 
}else{ 
  print("Do not reject H0") 
  print("Conclusion: With 95% confidence level, model 2 is adequate.") 
}

LR_stat2 