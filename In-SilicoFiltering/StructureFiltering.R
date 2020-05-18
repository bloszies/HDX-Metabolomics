setwd("F:/Clay/D2O_Study/GitHub")
data <- readxl::read_excel("F:/Clay/D2O_Study/Mass Frontier/FragmentsAnalysis_0919_ver3.xlsx")

data_test <- FillDown(data, 'AlignmentID') # Start: 11:34
data_test <- FillDown(data_test, 'InChIKey')
data_test <- FillDown(data_test, 'True Adduct')

MF_data <- data_test
write.table(MF_data, file = "MFFragments.txt")
save(MF_data, file = "MFFragments.Rdata")
load("MFFragments.Rdata")

MF_DSCounts <- read.delim("F:/Clay/D2O_Study/Mass Frontier/Fragments_2ndAttempt/FragmentData.txt", stringsAsFactors = F)
colnames(MF_DSCounts) <- c("Frag_InChIKey", "DSCount", "ExactWeight")
library("dplyr")
MF_DSCountsNoDups <- distinct(MF_DSCounts, Frag_InChIKey, .keep_all = TRUE)

library("data.table")
MF_FragsDScounts <- merge(MF_data, MF_DSCountsNoDups, by = "Frag_InChIKey")
MF_FragsDScounts <- as.data.frame(c(MF_FragsDScounts[,1:7],MF_FragsDScounts[,10:12], MF_FragsDScounts[,27:29]))
write.table(MF_FragsDScounts, file = "MF_FragsMetaData.txt")
save(MF_FragsDScounts, file = "MF_FragsDScounts.RData")
load("MF_FragsDScounts.RData")

####### Correcting the Parent DSCounts to reflect Protonation #######
MF_FragsDScounts$DSCount <- as.numeric(MF_FragsDScounts$DSCount)
MF_FragsDScounts$DSCount[which(grepl("+H", MF_FragsDScounts$MF_Formula, fixed = T) == T)] <- 
  MF_FragsDScounts$DSCount[which(grepl("+H", MF_FragsDScounts$MF_Formula, fixed = T) == T)] + 1

####### Preparing the experimental data w/ Identifiers ##########
ExpHDX <- readxl::read_excel("F:/Clay/D2O_Study/Mass Frontier/FragmentsAnalysis_0919_ver4.xlsx", sheet = "Sheet3")
ExpHDX <- ExpHDX[2:18]

##########Started for loop @ 9:20
MatchList1 <- list()
MatchList2 <- list()
MatchList3 <- list()
MatchList4 <- list()
MatchList5 <- list()
for (i in 1:nrow(ExpHDX)) {
  MF_Matches <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i]),]
  CorrectHDX1 <- numeric(5)
  CorrectHDX2 <- numeric(5)
  CorrectHDX3 <- numeric(5)
  CorrectHDX4 <- numeric(5)
  CorrectHDX5 <- numeric(5)
  for (j in 1:nrow(MF_Matches)) {
    CorrectHDX1[j] <- list(MF_Matches[which((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak1[i]) < 0.1) & (MF_Matches$DSCount[j] == ExpHDX$HDX1[i])),])
    CorrectHDX2[j] <- list(MF_Matches[which((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak2[i]) < 0.1) & (MF_Matches$DSCount[j] == ExpHDX$HDX2[i])),])
    CorrectHDX3[j] <- list(MF_Matches[which((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak3[i]) < 0.1) & (MF_Matches$DSCount[j] == ExpHDX$HDX3[i])),])
    CorrectHDX4[j] <- list(MF_Matches[which((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak4[i]) < 0.1) & (MF_Matches$DSCount[j] == ExpHDX$HDX4[i])),])
    CorrectHDX5[j] <- list(MF_Matches[which((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak5[i]) < 0.1) & (MF_Matches$DSCount[j] == ExpHDX$HDX5[i])),])
  }
  MatchList1[[i]] <- CorrectHDX1
  MatchList2[[i]] <- CorrectHDX2
  MatchList3[[i]] <- CorrectHDX3
  MatchList4[[i]] <- CorrectHDX4
  MatchList5[[i]] <- CorrectHDX5
}

####### For Loop that works, just can't figure out how to paste the data into MF_FragsDScounts (big spreadsheet) #####
ExpHDX[is.na(ExpHDX)] <- 0
for (i in 1:nrow(ExpHDX)) {
  MF_Matches <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i]),]
  MF_Matches$Match <- NA
  for (j in 1:nrow(MF_Matches)) {
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak1[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX1[i])) {
      MF_Matches$Match[j] <- "Peak1"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak2[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX2[i])) {
      MF_Matches$Match[j] <- "Peak2"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak3[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX3[i])) {
      MF_Matches$Match[j] <- "Peak3"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak4[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX4[i])) {
      MF_Matches$Match[j] <- "Peak4"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak5[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX5[i])) {
      MF_Matches$Match[j] <- "Peak5"
    }
  }
}

####### For Loop that pastes matches into ExpHDX #####
ExpHDX[is.na(ExpHDX)] <- 0
ExpHDX$Match1 <- NA
ExpHDX$Match2 <- NA
ExpHDX$Match3 <- NA
ExpHDX$Match4 <- NA
ExpHDX$Match5 <- NA
ExpHDX$Match6 <- NA
for (i in 1:nrow(ExpHDX)) {
  MF_Matches <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i]),]
  if (nrow(MF_Matches) < 1) {
    next
  }
  for (j in 1:nrow(MF_Matches)) {
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak1[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX1[i])) {
      ExpHDX$Match1[i] <- "TRUE"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak2[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX2[i])) {
      ExpHDX$Match2[i] <- "TRUE"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak3[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX3[i])) {
      ExpHDX$Match3[i] <- "TRUE"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak4[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX4[i])) {
      ExpHDX$Match4[i] <- "TRUE"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak5[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX5[i])) {
      ExpHDX$Match5[i] <- "TRUE"
    }
    if ((abs(MF_Matches$MonoisotopicMass[j] - ExpHDX$Peak6[i]) < 0.05) & (MF_Matches$DSCount[j] == ExpHDX$HDX6[i])) {
      ExpHDX$Match6[i] <- "TRUE"
    }
  }
}
write.csv(ExpHDX, file = "Processed HDX.csv")





colnames(data)
for (j in 7:(length(ExpHDX)-1)) {
  if (MF_Matches$MonoisotopicMass - ExpHDX$Peak1[i] < 0.05) {
  }
}
Match1 <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i] & (abs(MF_FragsDScounts$MonoisotopicMass - ExpHDX$Peak1[i]) < 0.05) &  
                                   (MF_FragsDScounts$DSCount == ExpHDX$HDX1[1])),]
Match2 <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i] & (abs(MF_FragsDScounts$MonoisotopicMass - ExpHDX$Peak2[i]) < 0.05) &  
                                   (MF_FragsDScounts$DSCount == ExpHDX$HDX2[1])),]
Match3 <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i] & (abs(MF_FragsDScounts$MonoisotopicMass - ExpHDX$Peak3[i]) < 0.05) & 
                                   (MF_FragsDScounts$DSCount == ExpHDX$HDX3[1])),]
Match4 <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i] & (abs(MF_FragsDScounts$MonoisotopicMass - ExpHDX$Peak4[i]) < 0.05) & 
                                   (MF_FragsDScounts$DSCount == ExpHDX$HDX4[1])),]
Match5 <- MF_FragsDScounts[which(MF_FragsDScounts$InChIKey == ExpHDX$InChIKey[i] & (abs(MF_FragsDScounts$MonoisotopicMass - ExpHDX$Peak5[i]) < 0.05) & 
                                   (MF_FragsDScounts$DSCount == ExpHDX$HDX5[1])),]
