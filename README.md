# HDX-Metabolomics
Underlying R script for database reduction and experimental data analysis

# Background/Purpose
The workflow described here outlines data processing for hydrogen/deuterium exchange (HDX) mass spectrometry based compound identification. Briefly, HDX aids in compound identification because it highlights the exchangeable protons in chemical structures. By incorporating deuterium into LC mobile phase composition, exchangeable protons for all compounds in the sample should be replaced with deuterium, resulting in a 1 Da mass shift. By investigating both MS1 and MS/MS data, it is possible to determine the number and positions of exchangeable protons, which can greatly aid in filtering isomer candidates so that only high confidence potential structural matches remain. 

## What is contained in this Repository:
- Experimental data analysis for hydrogen/deuterium exchange-based compound ID (outlined below)
  - Sample Data: Data_NoLabel.xlsx and Data_withLabel.xlsx
  - Script: HDX_Matching.R
  - Sample output file: Data_NoLabel_HDX_output.txt

The script outlined below solves one of the underlying issues in HDX data processing, linking unlabeled and labeled features across data files. To do this, two example datasets were generated using a ThermoScientific LC-Q Exactive instrument operated in data dependent-MSMS mode with a BEH Amide column. One data set (Data_NoLabel.xlsx) was acquired using mobile phases consisting of 0.125% formic acid and 10 mM ammonium formate dissolved in water (A) or acetonitrile/water (95:5) (B). The other data set (Data_withLabel.xlsx) was acquired using 0.125% d2-Formic acid and 10mM d5-Ammonium Formate dissolved in D2O (A) or acetonitrile/D2O (95:5) (B). No matter what acquisition parameters are used, this script will only work if the same chromatography is used across different runs, with the only difference being deuterium oxide replacing water and deuterated buffers/mobile phase modifiers replacing their unlabeled counterparts. In other words, the retention times should be comparable across the two data sets, or at least the retention time differences should be relatively predictable.

# Step 0: Data Preparation
Save the two sample datasheets into the same directory, for this example, we will call it "C:\HDX". These data sheets were prepared by first processing the raw data generated via LC-Q Exactive in MS-DIAL. From MS-DIAL, the data was exported in .txt format, and all columns other than those needed for this scriopt were removed. It should be noted that if left in the datasheets, those columns would be untouched and keep their original data, they were removed simply to avoid confusion. Prior to processing data, set the working directory to the saved folder:
```
setwd("C:/HDX")
```

# Step 1: Data input
Two datasheets containing the HDX and non-HDX runs are needed to link unlabeled features to their labeled counterparts.
These sheets should be inputted in .xlsx format with the required fields being: ID, retention time (RT), ms1 mass (MZ), and MS/MS pattern (MSMS; as m/z1:intensity1 mz2:intensity2 ...).
```
#Data Input
Data_NoLabel <- readxl::read_excel("Data_NoLabel.xlsx", sheet = "DataDictionary")
Data_withLabel <- readxl::read_excel("Data_withLabel.xlsx", sheet = "DataDictionary")

Data_NoLabel_MSMS <- Data_NoLabel[which(!is.na(Data_NoLabel$MSMS)),]
Data_withLabel_MSMS <- Data_withLabel[which(!is.na(Data_withLabel$MSMS)),]
```

# Step 2: For loop to identify deuterium exchanged features connected to unlabeled features
We used three criteria for the matching:
1. Retention time difference between labeled and unlabeled features is less than 0.5 minutes
2. Parent mass difference is equal to a multiple of the difference between hydrogen and deuterium masses (1.006277 Da), within a 1 mDa (0.001 Da) window.
3. The five most abundant unlabeled MSMS features must be present in the labeled MSMS with mass differences coming from different deuterium incorporation (they must also have a mass difference equal to a multiple of 1.006277 within a 1 mDa window).
```
for (i in 1:nrow(Data_NoLabel_MSMS)) {
  RT1 <- as.numeric(Data_NoLabel_MSMS$RT[i])
  mz1 <- Data_NoLabel_MSMS$MZ[i]
  Label_Features <- Data_withLabel_MSMS[which(Data_withLabel_MSMS$RT > (RT1 - 0.2) & Data_withLabel_MSMS$RT < (RT1 + 0.2) & 
    Data_withLabel_MSMS$MZ >= mz1 & (Data_withLabel_MSMS$MZ - mz1) < 20 & (((Data_withLabel_MSMS$MZ - mz1) %% 1.006277) < 0.001 |
    ((Data_withLabel_MSMS$MZ - mz1) %% 1.006277) > 1.005277)),]
  if(nrow(Label_Features) < 1){
    Data_NoLabel_MSMS$LabeledID[i] <- NA
    Data_NoLabel_MSMS$ExchangeNumber[i] <- NA
    next
  }
  Label_Features.list <- lapply(1:nrow(Label_Features), function(l){
    MSMS.df <- data.frame(do.call(rbind,lapply(strsplit(Label_Features$MSMS[l], " ")[[1]], function(x) {
      strsplit(x, ":")[[1]]
    })),stringsAsFactors = F)
    MSMS.df$X1 <- as.numeric(MSMS.df$X1)
    MSMS.df$X2 <- as.numeric(MSMS.df$X2)
    MSMS.df$ID <- Label_Features$ID[l]
    MSMS.df$HDX <- round(Label_Features$MZ[l] - mz1)
    MSMS.df
  })
  Label_Features.list.df <- data.frame(do.call(rbind, Label_Features.list),stringsAsFactors = F)
  
  MSMS.df.1 <- data.frame(do.call(rbind,lapply(strsplit(Data_NoLabel_MSMS$MSMS[i], " ")[[1]], function(x) {
    strsplit(x, ":")[[1]]
  })),stringsAsFactors = F)
  MSMS.df.1$X1 <- as.numeric(MSMS.df.1$X1)
  MSMS.df.1$X2 <- as.numeric(MSMS.df.1$X2)
  MSMS.df.1 <- MSMS.df.1[order(MSMS.df.1$X2, decreasing = T),]
  MSMS.df.1$NoLabel_ID <- Data_NoLabel_MSMS$ID[i]
  MSMS.df.1$NoLabel_MS1 <- Data_NoLabel_MSMS$MZ[i]
  MSMS.df.1$NoLabel_RT <- Data_NoLabel_MSMS$RT[i]
  
  MS2_Matching <- lapply(1:3, function(j) {
    mz2 <- MSMS.df.1$X1[j]
    MSMS_match <- Label_Features.list.df[which((((Label_Features.list.df$X1 - mz2) %% 1.006277) < 0.001 | ((Label_Features.list.df$X1 - 
                                                                                                              mz2) %% 1.006277) > 1.005277) & round(Label_Features.list.df$X1 - mz2) <= Label_Features.list.df$HDX),]
    if(nrow(MSMS_match) == 0){
      MSMS_match <- data.frame(NULL)
    } else {
      names(MSMS_match) <- c("MZ","Intensity","ID","HDX")
      cbind(MSMS.df.1[j,], MSMS_match)
    }
  })
  MS2_Matching.df <- data.frame(do.call(rbind, MS2_Matching),stringsAsFactors = F)
  if(nrow(MS2_Matching.df) == 0){
    Data_NoLabel_MSMS$LabeledID[i] <- NA
    Data_NoLabel_MSMS$ExchangeNumber[i] <- NA
    next
  }
  MS2_Matching.df$MS2_HDX <- round(MS2_Matching.df$MZ - MS2_Matching.df$X1)
  names(MS2_Matching.df) <- c("NoLabel_MS2","NoLabel_Frag_intensity","NoLabel_ID","NoLabel_MS1","NoLabel_RT","Labeled_MS2", "Labeled_Frag_intensity", "LabeledID","MS1_HDX","MS2_HDX")
  Data_NoLabel_MSMS$LabeledID[i] <- paste(unique(MS2_Matching.df$LabeledID), collapse = ";")
  Data_NoLabel_MSMS$ExchangeNumber[i] <- paste(round((Label_Features$MZ) - mz1), collapse = ";")
}
```
# Step 3: Writing Output file
The output file should look exactly the same as the input file with some additional fields appended to the end. These fields will be "Labeled_ID" and "ExchangeNumber".
Labeled_ID: This field will contain the IDs from the labeled data file that match the three criteria outlined in step 2. If there are multiple matches, they will be separated by semicolons (';').
ExchangeNumber: This field is calculated by rounding the mass difference between labeled and unlabeled features. This field will give the number of deuterium atoms exchanged for protons in the chemical structure for each of the matched IDs. Again, if there are multiple matches, the exchange numbers will be separated by semicolons.
For multiple matches the order is maintained, for example: if IDs 110, 145, and 170 are all found to match one unlabeled feature with 4, 5, and 6 H/D exchanges respectively, the "Labeled_ID" field will look like "110;145;170" and the "ExchangeNumber" field will be "4;5;6".
```
write.table(Data_NoLabel, file = "Data_NoLabel_HDX.txt", sep = "\t")
```
