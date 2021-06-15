#Purpose: To convert MS signal intensities to the corresponding QC loading amounts according to PH-PA selection and then perform
#         quantitative comparison in untargeted metabolomics

#1. Make notation for metabolic features in method blank
#2. Determine to use PH or PA for each metabolic feature 
#3. Using PH and PA, build calibration curves using QC serial injections respectively
#4. Perform signal calibration to convert MS signal intensities to QC loading amounts
#5. Generate the new feature table.

# Created: 2020-07-06
# Edited : 2020-07-31  
# Created by: Huaxu Yu
# Edited by : Huaxu Yu
# Version: 0.0.1
# copy right 2019 Huan lab @ University of British Columbia

##################################################
#Parameter setting
 #File input
Calibration_datapath = "F:/Results/PH&PA/Application/Application_test_20200707/Final_test"                                                              # the data path for the metabolite intensity table
Peak_area_FileName = "Leukemia_Peak_area.csv"                                          # the file name of peak area intensity table                                 
Peak_height_FileName = "Leukemia_Peak_height.csv"                                      # the file name of peak height intensity table
QC_PA_FileName = "QC_Peak_area.csv"                                                    # the file name of QC intensity table (peak area)
QC_PH_FileName = "QC_Peak_height.csv"                                                  # the file name of QC intensity table (peak height)
 #QC information
QC_conc = c(0.1,0.3,0.5,0.7,1)                                                         # volumn points of QC serial injections
rep_QC = 3                                                                             # number of the injection replicates at each QC volumn point
 #Calibration settings
QC_zero_num = 1                                                                        # number of QC calibration points that have zero values
R2_threshold = 0.8                                                                     # the feature with R2 value lager than this threshold will be qualified for signal calibration
k_threshold = 0                                                                        # the feature with k value lager than this threshold will be qualified for signal calibration

##################################################

#Load files
setwd(Calibration_datapath)
PA_sample = read.csv(Peak_area_FileName, stringsAsFactors = FALSE)
PH_sample = read.csv(Peak_height_FileName, stringsAsFactors = FALSE)
QC_PA = read.csv(QC_PA_FileName, stringsAsFactors = FALSE)
QC_PH = read.csv(QC_PH_FileName, stringsAsFactors = FALSE)

New_IntensityTable = data.frame(matrix(nrow = nrow(PA_sample), ncol = (ncol(PA_sample)+1)))
colnames(New_IntensityTable)[1:ncol(PA_sample)] = colnames(PA_sample)
colnames(New_IntensityTable)[(ncol(PA_sample)+1):ncol(New_IntensityTable)] = c("Notation")
New_IntensityTable[,1:3] = PA_sample[,1:3]

######################################################

#Calculate average and RSD for each QC serial injection point by using peak height and peak area respectively
QC_PA_summary = data.frame(matrix(ncol = 2*length(QC_conc)+3, nrow = nrow(PA_sample)))
colnames(QC_PA_summary)[(2*length(QC_conc)+1) : (2*length(QC_conc)+3)] = c("k","b","R2")
QC_PH_summary = QC_PA_summary
for (i in 1:length(QC_conc)) {
  PA_matrix = QC_PA[,(4+(i-1)*rep_QC):(4+i*rep_QC-1)]
  PH_matrix = QC_PH[,(4+(i-1)*rep_QC):(4+i*rep_QC-1)]
  QC_PA_summary[,i] = rowMeans(PA_matrix)
  QC_PH_summary[,i] = rowMeans(PH_matrix)
  QC_PA_summary[,length(QC_conc)+i] = apply(PA_matrix, 1, sd) / rowMeans(PA_matrix)
  QC_PH_summary[,length(QC_conc)+i] = apply(PH_matrix, 1, sd) / rowMeans(PH_matrix)
}

# Calculate k, b, R2 in QC serial injections by using peak height and peak area respectively.
for (i in 1:nrow(QC_PA_summary)) {
  
  PA_regression = lm(as.numeric(QC_PA_summary[i,1:length(QC_conc)]) ~ QC_conc)
  QC_PA_summary$k[i] = as.numeric(PA_regression$coefficients[2])
  QC_PA_summary$b[i] = as.numeric(PA_regression$coefficients[1])
  QC_PA_summary$R2[i] = as.numeric(summary(PA_regression)[8])
  
  PH_regression = lm(as.numeric(QC_PH_summary[i,1:length(QC_conc)]) ~ QC_conc)
  QC_PH_summary$k[i] = as.numeric(PH_regression$coefficients[2])
  QC_PH_summary$b[i] = as.numeric(PH_regression$coefficients[1])
  QC_PH_summary$R2[i] = as.numeric(summary(PH_regression)[8])
  
  if(sum(as.numeric(QC_PH_summary[i,1:length(QC_conc)] == 0)) > QC_zero_num){
    New_IntensityTable$Notation[i] = "Insufficient QC data points"
  } else if(QC_PH_summary$k[i] < k_threshold | QC_PH_summary$R2[i] < R2_threshold) {
    New_IntensityTable$Notation[i] = "Poor regression"
  } else {
    New_IntensityTable$Notation[i] = "Good regression"
  }
  
}

#For each metabolic feature, find the closest intensity in QC serial injections
#Make the decesion of using peak height or peak area for each feature in each sample
Decision_table = data.frame(matrix(ncol = ncol(PA_sample), nrow = nrow(PA_sample)))
colnames(Decision_table) = colnames(PA_sample)
Decision_table[,1:3] = PA_sample[,1:3]

for (i in 1:nrow(PH_sample)) {
  for (j in 1:(ncol(PH_sample)-3)) {
    int_diff = abs(PH_sample[i,3+j]-QC_PH_summary[i,1:length(QC_conc)])
    cor_QC_position = match(min(int_diff),int_diff)
    sum_position = cor_QC_position + length(QC_conc)
    if(QC_PH_summary[i,sum_position] < QC_PA_summary[i,sum_position] | is.nan(QC_PH_summary[i,sum_position]) | is.nan(QC_PA_summary[i,sum_position])){
      Decision_table[i,3+j] = "PH"
      New_IntensityTable[i,3+j] = (PH_sample[i,3+j]-QC_PH_summary$b[i])/QC_PH_summary$k[i]
    } else{
      Decision_table[i,3+j] = "PA"
      New_IntensityTable[i,3+j] = (PA_sample[i,3+j]-QC_PA_summary$b[i])/QC_PA_summary$k[i]
    }
  }
}

for (i in 1:nrow(PH_sample)) {
  if(New_IntensityTable$Notation[i] != "Good regression"){
    New_IntensityTable[i,4:ncol(PA_sample)] = PH_sample[i,4:ncol(PA_sample)]
  }
}

for (i in 1:nrow(New_IntensityTable)) {
  for (j in 4:17) {
    if(New_IntensityTable[i,j] < 0){
      New_IntensityTable[i,j] = 0
    }
  }
}
Decision_table = cbind(Decision_table,New_IntensityTable$Notation)
write.csv(New_IntensityTable, file = "Calibrated intensity table.csv",row.names = FALSE)
write.csv(Decision_table, file = "decision table.csv",row.names = FALSE)
