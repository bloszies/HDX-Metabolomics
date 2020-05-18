setwd("F:/Clay/D2O_Study/Github")
#Data Input
HILIC <- readxl::read_excel("Urine2_posCombined.xlsx", sheet = "NL_DataDict")
HDX <- readxl::read_excel("Urine2_posCombined.xlsx", sheet = "wL_DataDict")
HDX$`Average Rt(min)` <- as.numeric(HDX$`Average Rt(min)`)

HILIC_MSMS <- HILIC[which(HILIC$`MS/MS assigned` == T),]
HDX_MSMS <- HDX[which(HDX$`MS/MS assigned` == T),]

for (i in 1:nrow(HILIC_MSMS)) {
  RT1 <- as.numeric(HILIC_MSMS$`Average Rt(min)`[i])
  mz1 <- HILIC_MSMS$`Average Mz`[i]
  Lab_Features <- HDX_MSMS[which(HDX_MSMS$`Average Rt(min)` > (RT1 - 0.2) & HDX_MSMS$`Average Rt(min)` < (RT1 + 0.2) & HDX_MSMS$`Average Mz` >= mz1 & 
                                      (HDX_MSMS$`Average Mz` - mz1) < 20 & (((HDX_MSMS$`Average Mz` - mz1) %% 1.006277) < 0.001 | ((HDX_MSMS$`Average Mz` - mz1) %% 1.006277) > 1.005277)),]
  if(nrow(Lab_Features) < 1){
    HILIC_MSMS$`wL Alignment IDs`[i] <- NA
    HILIC_MSMS$`wL Exchange(s)`[i] <- NA
    next
  }
  Lab_Features.list <- lapply(1:nrow(Lab_Features), function(l){
    MSMS.df <- data.frame(do.call(rbind,lapply(strsplit(Lab_Features$`MS/MS spectrum`[l], " ")[[1]], function(x) {
      strsplit(x, ":")[[1]]
    })),stringsAsFactors = F)
    MSMS.df$X1 <- as.numeric(MSMS.df$X1)
    MSMS.df$X2 <- as.numeric(MSMS.df$X2)
    MSMS.df$ID <- Lab_Features$`Alignment ID`[l]
    MSMS.df$HDX <- round(Lab_Features$`Average Mz`[l] - mz1)
    MSMS.df
  })
  
  Lab_Features.list.df <- data.frame(do.call(rbind, Lab_Features.list),stringsAsFactors = F)
  
  
  MSMS.df.1 <- data.frame(do.call(rbind,lapply(strsplit(HILIC_MSMS$`MS/MS spectrum`[i], " ")[[1]], function(x) {
    strsplit(x, ":")[[1]]
  })),stringsAsFactors = F)
  MSMS.df.1$X1 <- as.numeric(MSMS.df.1$X1)
  MSMS.df.1$X2 <- as.numeric(MSMS.df.1$X2)
  MSMS.df.1 <- MSMS.df.1[order(MSMS.df.1$X2, decreasing = T),]
  MSMS.df.1$NL_AlignmentID <- HILIC_MSMS$`Alignment ID`[i]
  MSMS.df.1$NL_MS1 <- HILIC_MSMS$`Average Mz`[i]
  MSMS.df.1$NL_RT <- HILIC_MSMS$`Average Rt(min)`[i]
  
  MS2_Matching <- lapply(1:3, function(j) {
    mz2 <- MSMS.df.1$X1[j]
    MSMS_match <- Lab_Features.list.df[which((((Lab_Features.list.df$X1 - mz2) %% 1.006277) < 0.001 | ((Lab_Features.list.df$X1 - mz2) %% 1.006277) > 1.005277) & 
                                               round(Lab_Features.list.df$X1 - mz2) <= Lab_Features.list.df$HDX),]
    if(nrow(MSMS_match) == 0){
      MSMS_match <- data.frame(NULL)
    } else {
      names(MSMS_match) <- c("m/z","Intensity","ID","HDX")
      cbind(MSMS.df.1[j,], MSMS_match)
    }
  })
  MS2_Matching.df <- data.frame(do.call(rbind, MS2_Matching),stringsAsFactors = F)
  if(nrow(MS2_Matching.df) == 0){
    HILIC_MSMS$`wL Alignment IDs`[i] <- NA
    HILIC_MSMS$`wL Exchange(s)`[i] <- NA
    next
  }
  MS2_Matching.df$MS2_HDX <- round(MS2_Matching.df$m.z - MS2_Matching.df$X1)
  names(MS2_Matching.df) <- c("NL_MS2","NL_Frag_intensity","NL_AlignmentID","NL_MS1","NL_RT","wL_MS2", "wL_Frag_intensity", "wL_ID","MS1_HDX","MS2_HDX")
  HILIC_MSMS$`wL Alignment IDs`[i] <- paste(unique(MS2_Matching.df$wL_ID), collapse = ";")
  HILIC_MSMS$`wL Exchange(s)`[i] <- paste(unique(MS2_Matching.df$MS1_HDX), collapse = ";")
}

write.table(HILIC_MSMS, file = "MSMS_LabelMatching.txt", sep = "\t")
