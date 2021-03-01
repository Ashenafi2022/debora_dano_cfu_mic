# debora_dano_cfu_mic

# Date: Feb. 28, 2021
# Created by Ashenafi Beyi

## R packages
install.pckages("ggplot2")
library("ggpubr")
install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
install.packages("FSA") # Dunn test, paired wise nonparametric test

library(ggplot2)
library(ggpubr)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)   # alternatively, this also loads %>%
library(FSA)



# uploading dataset to R, check 
dano_cfu_mic <- read.csv("C:/Data analysis/Debora/MIC_CFU4analysis.csv", header=T)
dim(dano_cfu_mic)
str(dano_cfu_mic)
dano_cfu_mic$TimePoint <- as.factor(dano_cfu_mic$TimePoint) # convert timpoint from

##______________ CFU ____________________________________________
# Summary by group and time point
cfu_summary = group_by(dano_cfu_mic,Group_Time)%>%summarize(mean=mean(log10_CipR), 
                                            se=sd(log10_CipR)/sqrt(length(log10_CipR)))
                                            
# Kruskal Wallis test: cfu by group and sampling time point
kruskal.test(dano_cfu_mic$log10_CipR ~ dano_cfu_mic$Group_Time)
  
  #pair wise comparion using Dun test
dunnTest(dano_cfu_mic$log10_CipR ~ dano_cfu_mic$Group_Time, method="bh")

# Line graphs
mic_summaries <- group_by(dano_cfu_mic, Group, TimePoint) %>%
            dplyr::summarise(
              count = n(),
              mean_cfuR = mean(log10_CipR, na.rm = TRUE),sd_cfuR = sd(log10_CipR, na.rm = TRUE),
              median_cfuR = median(log10_CipR, na.rm = TRUE), IQR_cfuR = IQR(log10_CipR, na.rm = TRUE),
              
              mean_cfuT = mean(log10_total, na.rm = TRUE),sd_cfU = sd(log10_total, na.rm = TRUE),
              median_cfuT = median(log10_total, na.rm = TRUE), IQR_cfuT = IQR(log10_total, na.rm = TRUE))

ggline(mic_summaries, x="TimePoint", y= "median_cfuR", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "Ciprofloxacin resistant CFU (median, log10)", xlab = "Sampling day")
ggline(mic_summaries, x="TimePoint", y= "median_cfuT", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "Total CFU (median, log10)", xlab = "Sampling day")
       
       
       ## ..............Newdata sets by TREATMENT.................
dano_cfu_mic_TRT <- subset(dano_cfu_mic, Treatment=="Treatment",
                        select=TimePoint:Tetracycline_BP)
dano_cfu_mic_Control <- subset(dano_cfu_mic, Treatment=="Control",
                            select=TimePoint:Tetracycline_B)
# ..........Boxplot.....................
par(mfrow=c(2,5))
# Margins - Avoiding overlapping labels axes 
# https://stackoverflow.com/questions/6778425/avoid-overlapping-axis-labels-in-r
par(mar = c(3, 3, 3, 3), mgp = c(3, 1, 0))

  # Treatment
boxplot(dano_cfu_mic_TRT$log10_CipR ~ dano_cfu_mic_TRT$TimePoint,
        xlab = "Sampling day", ylab = "Ciprofloxacin resistant CFU (log10)", col= "10")
boxplot(dano_cfu_mic_Control$log10_total ~ dano_cfu_mic_Control$TimePoint,
        xlab = "Sampling day", ylab = "Total CFU (log10)", col= "10")
  
  #Control      
boxplot(dano_cfu_mic_Control$log10_CipR ~ dano_cfu_mic_Control$TimePoint,
        xlab = "Sampling day", ylab = "Ciprofloxacin resistant CFU (log10)", col= "10")
boxplot(dano_cfu_mic_Control$log10_total ~ dano_cfu_mic_Control$TimePoint,
        xlab = "Sampling day", ylab = "Total CFU (log10)", col= "10")



## ____________________ MIC - ciprofloxacin______________________________________
# Summary by group and time point
MIC_Cipro_summary = group_by(dano_cfu_mic,Group_Time)%>%summarize(mean=mean(Ciprofloxacin_MIC), 
                                            se=sd(Ciprofloxacin_MIC)/sqrt(length(Ciprofloxacin_MIC)))
                                            
# Kruskal Wallis test: cfu by group and sampling time point
kruskal.test(dano_cfu_mic$Ciprofloxacin_MIC ~ dano_cfu_mic$Group_Time)
   #pair wise comparion using Dun test
dunnTest(dano_cfu_mic$Ciprofloxacin_MIC ~ dano_cfu_mic$Group_Time, method="bh")


#......................Line graphs................................
# 01 Line graphs groups by sampling time points
mic_summaries <- group_by(dano_cfu_mic, Group, TimePoint) %>%
            dplyr::summarise(
              count = n(),
              mean_cipro = mean(Ciprofloxacin_MIC, na.rm = TRUE),sd_cipro = sd(Ciprofloxacin_MIC, na.rm = TRUE),
              median_cipro = median(Ciprofloxacin_MIC, na.rm = TRUE), IQR_cipro = IQR(Ciprofloxacin_MIC, na.rm = TRUE),
              
              mean_clinda = mean(Clindamycin_MIC, na.rm = TRUE),sd_clinda = sd(Clindamycin_MIC, na.rm = TRUE),
              median_clinda = median(Clindamycin_MIC, na.rm = TRUE), IQR_clinda = IQR(Clindamycin_MIC, na.rm = TRUE),
              
              mean_erythro = mean(Erythromicin_MIC, na.rm = TRUE),sd_erythro = sd(Erythromicin_MIC, na.rm = TRUE),
              median_erythro = median(Erythromicin_MIC, na.rm = TRUE), IQR_erythro = IQR(Erythromicin_MIC, na.rm = TRUE),
              
              mean_flor = mean(Florfenicol_MIC, na.rm = TRUE),sd_flor = sd(Florfenicol_MIC, na.rm = TRUE),
              median_flor = median(Florfenicol_MIC, na.rm = TRUE), IQR_flor = IQR(Florfenicol_MIC, na.rm = TRUE),
              
              mean_genta = mean(Gentamicin_MIC, na.rm = TRUE),sd_genta = sd(Gentamicin_MIC, na.rm = TRUE),
              median_genta = median(Gentamicin_MIC, na.rm = TRUE), IQR_genta = IQR(Gentamicin_MIC, na.rm = TRUE))
              
              mean_nalidic = mean(Nalidixic.acid_MIC, na.rm = TRUE),sd_nalidic = sd(Nalidixic.acid_MIC, na.rm = TRUE),
              median_nalidic = median(Nalidixic.acid_MIC, na.rm = TRUE), IQR_nalidic = IQR(Nalidixic.acid_MIC, na.rm = TRUE),
              
              mean_telithro = mean(Telithromycin_MIC, na.rm = TRUE),sd_telithro = sd(Telithromycin_MIC, na.rm = TRUE),
              median_telithro = median(Telithromycin_MIC, na.rm = TRUE), IQR_telithro = IQR(Telithromycin_MIC, na.rm = TRUE),
              
              mean_ttc = mean(Tetracycline_MIC, na.rm = TRUE),sd_ttc = sd(Tetracycline_MIC, na.rm = TRUE),
              median_ttc = median(Tetracycline_MIC, na.rm = TRUE), IQR_ttc = IQR(Tetracycline_MIC, na.rm = TRUE))
     

library(ggpubr)
ggline(mic_summaries, x="TimePoint", y= "median_cipro", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Ciprofloxacin (median)", xlab = "Sampling day")
ggline(mic_summaries, x="TimePoint", y= "median_clinda", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Clindamycin (median)", xlab = "Sampling day")
ggline(mic_summaries, x="TimePoint", y= "median_erythro", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Erythromicin (median)", xlab = "Sampling day")
ggline(mic_summaries, x="TimePoint", y= "median_flor", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Florfenicol (median)", xlab = "Sampling day")
ggline(mic_summaries, x="TimePoint", y= "median_genta", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Gentamicin (median)", xlab = "Sampling day")
ggline(mic_summaries, x="TimePoint", y= "median_nalidic", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Nalidixic Acid (median)", xlab = "Sampling day")       
ggline(mic_summaries, x="TimePoint", y= "median_telithro", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Telithromycin (median)", xlab = "Sampling day")
 ggline(mic_summaries, x="TimePoint", y= "median_ttc", linetype = "Group", 
       shape = "Group", color = "Group", siz=1.0,
       ylab = "MIC of Tetracycline (median)", xlab = "Sampling day")


# ..........Boxplot.....................
par(mfrow=c(2,1))
# Margins - Avoiding overlapping labels axes 
# https://stackoverflow.com/questions/6778425/avoid-overlapping-axis-labels-in-r
par(mar = c(3, 3, 3, 3), mgp = c(3, 1, 0))

  # Treatment
boxplot(dano_cfu_mic_TRT$Ciprofloxacin_MIC ~ dano_cfu_mic_TRT$TimePoint,
        xlab = "Sampling day", ylab = "Ciprofloxacin MIC", col= "10")
  
  #Control      
boxplot(dano_cfu_mic_Control$Ciprofloxacin_MIC ~ dano_cfu_mic_Control$TimePoint,
        xlab = "Sampling day", ylab = "Ciprofloxacin MIC (log10)", col= "10")




#___________________ Correlation of MIC's among different antibiotics____________________
    # Multiple comparisons
# https://rcompanion.org/rcompanion/f_01.html

    dim(Genus4RFTrt) # 138 X 240, 228 genera (Total 298 genera in 209 samples, zero inflatted onces are removed)
    
    which( colnames(Genus4RFTrt)=="Proteobacteria_Campylobacterales_Campylobacteraceae_g_Campylobacter" ) # 205
    match("Proteobacteria_g_Campylobacter", names(Genus4RFTrt)) # 205
    colnames(Genus4RFTrt[205]) # "Campylobacter"
       
    #Pearson .....False discovery rate adjusted p-value
    empty_micTRT <- data.frame(Corr = rep(0, 11),  p_value = rep(0, 11)) # 11 is the #columns compared
    k = c(8:18) #the MIC columns re found in columns 10 to 18
    corr.output.micTRT <- data.frame(rep(empty_micTRT, ))
    for(i in k){
      cor.assoc <- cor.test(dano_cfu_mic_TRT[8:18], dano_cfu_mic_TRT[8:18, i], type="pearson")
      corr.output.micTRT[i,] <- cbind(corr=cor.assoc$estimate,  p_value=cor.assoc$p.value)
    }
    corr.pear.micTRT
    corr.pear.micTRT$FDR = p.adjust(corr.pear.micTRT$p_value, method = "fdr") # False discovery rate adjusted p-value
    corr.pear.micTRT <- as.data.frame(corr.pear.micTRT)
    
    # Spearman.........False discovery rate adjusted p-value
    corr.spear.micTRT <- data.frame(rep(empty_micTRT, ))
    for(i in k){
      cor.assoc <- cor.test(dano_cfu_mic_TRT[8:18], dano_cfu_mic_TRT[8:18, i], type="spearman")
      corr.output.micTRT[i,] <- cbind(corr=cor.assoc$estimate,  p_value=cor.assoc$p.value)
     }
    corr.spear.micTRT
    corr.spear.micTRT$FDR = p.adjust(corr.spear.micTRT$p_value, method = "fdr") # False discovery rate adjusted p-value
    corr.spear.micTRT <- as.data.frame(corr.spear.micTRT)
    
    #Cobine both pearson and spearman correlation files and export 
    CorrAdj.micTRT <- cbind.data.frame(pear.r = corr.pear.micTRT$Corr, pear.p_raw = corr.pear.micTRT$p_value, pear.p_fdr = corr.pear.micTRT$FDR, 
                                         spear.r = corr.spear.micTRT$Corr, spear.p_raw = corr.spear.micTRT$p_value, spear.p_fdr = corr.spear.micTRT$FDR)
    
    #EXPORT FILE
    install.packages("rio")
    library("rio")
    export(CorrAdj.genusTrt, file="C:/Data analysis/DanofloxacinData/Debora/Corr.micTRT.xlsx")
    
    cor.test(corr.pear.micTRT$Corr, corr.spear.micTRT$Corr) # 1, Pearson and Pearson correlation values has 1.0 correlation cofficient
    


# Repeated measure correlation 
install.packages("rmcorr")
library(rmcorr)

mic_rmcorr <- rmcorr(TimePoint, Ciprofloxacin_MIC, Clindamycin_MIC,Erythromicin_MIC, 
Florfenicol_MIC, Gentamicin_MIC, Nalidixic.acid_MIC, Telithromycin_MIC,Tetracycline_MIC,
dano_cfu_mic, CIs = c("analytic,", "bootsrap"), nreps = 100, 
bstrap.out = F)
