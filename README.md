# HDX-Metabolomics
Underlying R script for database reduction and experimental data analysis

# Step 1: Data input
Two datasheets containing the HDX and non-HDX runs are needed to link unlabeled features to their labeled counterparts.
These sheets should be inputted in .xlsx format with fields for ID, retention time (RT), ms1 mass (MZ), and MS/MS pattern (MSMS; as m/z1:intensity1 mz2:intensity2 ...)
```
#Data Input
Data_NoLabel <- readxl::read_excel("Data_NoLabel.xlsx", sheet = "DataDictionary")
Data_withLabel <- readxl::read_excel("Data_withLabel.xlsx", sheet = "DataDictionary")
```
(optional) Only include features for which MS/MS was collected
```
Data_NoLabel_MSMS <- Data_NoLabel[which(!is.na(Data_NoLabel$MSMS)),]
Data_withLabel_MSMS <- Data_withLabel[which(!is.na(Data_withLabel$MSMS)),]
```

# Step 2: For loop to identify deuterium exchanged features connected to unlabeled features
We used three criteria for the matching:
1. Retention time difference between labeled and unlabeled features is less than 0.5 minutes
2. Parent mass difference is equal to a multiple of the difference between hydrogen and deuterium masses (1.006277 Da), within a 1 mDa (0.001 Da) window.
3. The five most abundant unlabeled MSMS features must be present in the labeled MSMS with mass differences coming from different deuterium incorporation (they must also have a mass difference equal to a multiple of 1.006277 within a 1 mDa window).
```
for (i in 1:nrow(Data_NoLabel)) {
  RT1 <- as.numeric(Data_NoLabel$RT[i])
  mz1 <- Data_NoLabel$MZ[i]
  Label_Features <- Data_withLabel[which(Data_withLabel$RT > (RT1 - 0.2) & Data_withLabel$RT < (RT1 + 0.2) & Data_withLabel$MZ >= mz1 &
    (Data_withLabel$MZ - mz1) < 20 & (((Data_withLabel$MZ - mz1) %% 1.006277) < 0.001 |((Data_withLabel$MZ - mz1) %% 1.006277) >
    1.005277)),]
  if(nrow(Label_Features) < 1){
    Data_NoLabel$LabeledID <- NA
    Data_NoLabel$ExchangeNumber <- NA
    next
  }
  Label_Features.list <- lapply(1:nrow(Label_Features), function(l){
    MSMS.df <- data.frame(do.call(rbind,lapply(strsplit(Label_Features$MSMS, " ")[[1]], function(x) {
      strsplit(x, ":")[[1]]
    })),stringsAsFactors = F)
    MSMS.df$ID <- Label_Features$ID[l]
    MSMS.df$HDX <- round(Label_Features$MZ[l] - mz1)
    MSMS.df
  })
  Label_Features.list.df <- data.frame(do.call(rbind, Label_Features.list),stringsAsFactors = F)
  
  MSMS.df.1 <- data.frame(do.call(rbind,lapply(strsplit(Data_NoLabel$MSMS[i], " ")[[1]], function(x) {
    strsplit(x, ":")[[1]]
  })),stringsAsFactors = F)
  MSMS.df.1$X1 <- as.numeric(MSMS.df.1$X1)
  MSMS.df.1$X2 <- as.numeric(MSMS.df.1$X2)
  MSMS.df.1 <- MSMS.df.1[order(MSMS.df.1$X2, decreasing = T),]
  MSMS.df.1$NoLabel_ID <- TVC_NL_MSMS$ID[i]
  MSMS.df.1$NoLabel_MS1 <- TVC_NL_MSMS$MZ[i]
  MSMS.df.1$NoLabel_RT <- TVC_NL_MSMS$RT[i]
  
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
    Data_NoLabel$LabeledID[i] <- NA
    TVC_NL_MSMS$ExchangeNumber[i] <- NA
    next
  }
  MS2_Matching.df$MS2_HDX <- round(MS2_Matching.df$MZ - MS2_Matching.df$X1)
  names(MS2_Matching.df) <- c("NoLabel_MS2","NoLabel_Frag_intensity","NoLabel_ID","NoLabel_MS1","NoLabel_RT","Labeled_MS2", "Labeled_Frag_intensity", "Labeled_ID","MS1_HDX","MS2_HDX")
  TVC_NL_MSMS$Labeled_ID[i] <- paste(unique(MS2_Matching.df$Labeled_ID), collapse = ";")
  TVC_NL_MSMS$ExchangeNumber[i] <- paste(round((Label_Features$MZ) - mz1), collapse = ";")
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
