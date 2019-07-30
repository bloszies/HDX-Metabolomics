setwd("F:/Clay/D2O_Study/Github")
#Data Input
Data_NoLabel <- readxl::read_excel("Data_NoLabel.xlsx", sheet = "DataDictionary")
Data_withLabel <- readxl::read_excel("Data_withLabel.xlsx", sheet = "DataDictionary")

Data_NoLabel_MSMS <- Data_NoLabel[which(!is.na(Data_NoLabel$MSMS)),]
Data_withLabel_MSMS <- Data_withLabel[which(!is.na(Data_withLabel$MSMS)),]

for (i in 1:nrow(Data_NoLabel_MSMS)) {
  RT1 <- as.numeric(Data_NoLabel_MSMS$RT[i])
  mz1 <- Data_NoLabel_MSMS$MZ[i]
  Label_Features <- Data_withLabel_MSMS[which(Data_withLabel_MSMS$RT > (RT1 - 0.2) & Data_withLabel_MSMS$RT < (RT1 + 0.2) & Data_withLabel_MSMS$MZ >= mz1 &
    (Data_withLabel_MSMS$MZ - mz1) < 20 & (((Data_withLabel_MSMS$MZ - mz1) %% 1.006277) < 0.001 |((Data_withLabel_MSMS$MZ - mz1) %% 1.006277) > 1.005277)),]
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

write.table(Data_NoLabel_MSMS, row.names = F, file = "Data_NoLabel_HDX.txt", sep = "\t")
