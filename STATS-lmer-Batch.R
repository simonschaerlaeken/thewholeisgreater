#Lancer les packages:
library(Matrix)
library(lme4)
library(MuMIn)
library(carData)
library(effects)
library(ggplot2)
library(pracma)
library(car)
library(phia)
library(lmerTest)
#library(plyr) 
library(dplyr)

rm(list=ls())

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                   ----------------                 VARIABLES              ------------------
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

filenameLog <- 'Log_iLead_STATS_LMER.txt'
filenameOutput <- 'Summary_iLead_STATS_LMER.csv'
directoryInput <-
  "/YourPath_input/"
filenameInput <-
  "ALL_Data_iLead_Zscore_SBAC_9Var_6Clu.csv"
directoryOutput <-
  "/YourPath_output/"

##------------- PARAMETERS -------------
analysis_type <- "cat" 
histogram_flag <- FALSE
y_prefix <- ""
y_suffix <- ""
x_prefix <- ""
x_suffix <- ""

##------------- LM demographic -------------

# All Y
y_list <- c(
 'zscore_MATH_FLUENCY.acc_sum.overall',
 'zscore_READING_COMPREHENSION.acc_score.overall',
 'zscore_combined.Math',
 'zscore_combined.ELA'
)

# All X
x_list <- c(
  'Clusters + Language.Fluency + Gender + Parent.Ed.Lvl + Ethnicity + (1|pid) + (1|School)') 

effect_variable <- 'Clusters'

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
#                   ----------------                 LINEAR MODELS            ------------------
#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


# Start writing to an output file
options(warn=1)

# OPEN FILE
dat<-read.csv(paste(directoryInput,filenameInput,sep=""),header=TRUE, sep=',', na.strings = "NaN")
names(dat)
summary(dat)

# --- Modify data
dat$Clusters <- factor(dat$Clusters)
dat$Gender <- factor(dat$Gender)
dat$Grade <- factor(dat$Grade)
dat$Cohort <- factor(dat$Cohort)
dat$Time <- factor(dat$Time)

# CHECK RESULTS FOLDER EXIST
if (!file.exists(directoryOutput)){
  dir.create(directoryOutput)
}

# Log
fileLog <- file(paste(directoryOutput,filenameLog,sep = ""), open="wt")
sink(paste(directoryOutput,filenameLog,sep = ""))

# Batch stats
for (i in 1:length(y_list)){
  warning(paste('--- VI: ', y_list[i]))
  print("----------------------------------------")
  print("----------------------------------------")
  print(paste('    VI: ', y_list[i]))
  print("----------------------------------------")
  print("----------------------------------------")
  
  for (j in 1:length(x_list)){
    warning(paste('VD: ', x_list[j]))
    print("                                        ")
    print("##################################################")
    print(paste('  VD: ', x_list[j]))
    print("##################################################")
    print("                                        ")
    
    y_name = y_list[i]
    x_name = x_list[j]
    
    # Select Rows based on preprocessing criteria
    # 1 - No English Learners
    dat_el <- dat %>% filter(is.na(Language.Fluency)|Language.Fluency!="English Learner")
    # 2 - No SpEd
    dat_sp <- dat_el %>% filter(is.na(SpEd)|SpEd!="Yes")
    # 3 - Only T2 and T4
    dat_T2T4 <- dat_sp %>% filter(is.na(Time)|Time == "T2" | Time == "T4")
    # 4 - Above Chance level for Math and Reading
    if (strcmp(y_name, "zscore_MATH_FLUENCY.acc_sum.overall")){
      dat_math <- dat_T2T4 %>% mutate(Math_Accuracy = MATH_FLUENCY.acc_sum.overall/MATH_FLUENCY.acc_length.overall)
      dat_selected <- dat_math %>% filter(is.na(Math_Accuracy)|Math_Accuracy >= 0.5)
    }
    else if (strcmp(y_name, "zscore_READING_COMPREHENSION.acc_score.overall")){
      dat_reading <- dat_T2T4 %>% mutate(Reading_Accuracy = READING_COMPREHENSION.acc_sum.overall/READING_COMPREHENSION.acc_length.overall)
      dat_selected <- dat_reading %>% filter(is.na(Reading_Accuracy)|Reading_Accuracy >= 0.5)
    }
    else {
      dat_selected <- dat_T2T4
    }
    
    # Create text
    model_name_text <- paste(y_name,"~",x_name,sep = "")
    
    # --- LMER
    lm_text <- paste("lmer(",model_name_text,", data=dat_selected)",sep = "")
    model.lm<-eval(parse(text = lm_text))
    warning("Model Finished")
    
    # Print Result Model
    print(summary(model.lm))
    print(anova(model.lm))
    print(r.squaredGLMM(model.lm))
    
    # Print Contrast
    cont_text <- paste("testInteractions(model.lm, pairwise = '",effect_variable,"',adjustment = 'fdr')",sep = "")
    print(eval(parse(text = cont_text)))
    warning("Contrasts Finished")

    # PLOT EFFECT AND 95%
    #Plot Model
    plot_name <- paste(directoryOutput,y_list[i],"__",effect_variable, "_STATS_EFF_Line",".jpg", sep="")
    eff_text <- paste("effect(\'",effect_variable,"\', model.lm,typical=mean, confidence.level=0.95)",sep = "")
    eff <- print(eval(parse(text = eff_text)))
    eff_frame_text <- paste("data.frame(lower = eff$lower, upper = eff$upper, fit = eff$fit, ", effect_variable," = eff$x$",effect_variable,")",sep = "")
    eff_frame <- eval(parse(text = eff_frame_text))
    jpeg(plot_name)
    print(plot(eff))
    dev.off()
    
    # GGPLOT bar
    plot_name_pdf <- paste(directoryOutput,y_name,"__",effect_variable, "_STATS_EFF_bar",".pdf", sep="")
    pdf(plot_name_pdf)
    #png(plot_name_pdf,width = 960, height = 960, units = "px",res=150)
    plot_text <- paste("ggplot(eff_frame, aes(x =", effect_variable,", y = fit, color = ",effect_variable,", fill = ",effect_variable,")) +",
                       "geom_bar(stat='identity') +",
                       "geom_errorbar(aes(ymin= lower, ymax= upper), width=.2) +",
                       "scale_color_brewer(palette='Dark2') +",
                       "scale_fill_brewer(palette='Pastel2') +",
                       #"coord_cartesian(ylim = c(min(dat_selected$",y_name,"),max(dat_selected$",y_name,"))) +",
                       "theme_light()",sep = "")
    print(eval(parse(text = plot_text)))
    dev.off()
    warning("Plot EFF Finished")
  }
}
closeAllConnections()