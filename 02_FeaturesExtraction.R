# Initial stuff
library(Biostrings)
library(Peptides)
library(ggplot2)
library(optparse)

# Optparse
option_list <- list(
  make_option("--file1", action = "store", dest = "file1",
              help = "Input stringset ['name.fa'] for the positive condition."),
  make_option("--file0", action = "store", dest = "file0",
              help = "Input stringset ['name.fa'] for the negative condition."),
  make_option(c("-d", "--dframe"), action = "store", dest = "dframe",
              help = "Give a name to the data.frame with features. [Example: DF]"),
  make_option(c("-s", "--stdf"), action = "store_true", dest = "stdf",
              default = FALSE, help = "Make a stadistic data.frame.
              [default: %default]"),
  make_option(c("-b", "--boxp"), action = "store_true", dest = "boxplot",
              default = FALSE, help = "Make boxplot to compare features.
              [default: %default]"),
  make_option(c("-p", "--dplot"), action = "store_true", dest = "denplot",
              default = FALSE, help = "Make density to plots to compare features.
              [default: %default")
)
opt <- parse_args(OptionParser(option_list=option_list))


# Helper function to be called from within getProtFeatures
getAAComposition <- function(protSeqs) {
  AADataFrame <- data.frame(Peptides::aaComp(protSeqs))
  AADataFrame <- AADataFrame[, c(FALSE, TRUE)]
  AADataFrame <- data.frame(t(AADataFrame))
  rownames(AADataFrame) <- NULL
  colnames(AADataFrame) <- c("aaTiny_pc", "aaSmall_pc", "aaAliphatic_pc", "aaAromatic_pc",
                             "aaNonPolar_pc", "aaPolar_pc", "aaCharged_pc", "aaBasic_pc",
                             "aaAcidic_pc")
  return(AADataFrame)
}

# Function to obtain annotation
getProtFeatures <- function(protSeqDF) {
  protSeqs <- protSeqDF[, "sequence"] #Extrae la secuencia
  molecular_weight <- Peptides::mw(protSeqs, monoisotopic = FALSE) 
  pepLength <- Peptides::lengthpep(protSeqs) 
  isoelectric_point <- Peptides::pI(protSeqs, pKscale = "EMBOSS") 
  instability <- Peptides::instaIndex(protSeqs) 
  aliphaticIndex <- Peptides::aIndex(protSeqs)
  bomanIndex <- Peptides::boman(protSeqs) 
  AADataFrame <- getAAComposition(protSeqs) 
  protInfo <- as.data.frame(t(as.matrix(data.frame(strsplit(protSeqDF[, "seq_accession"], split='[|]', fixed=FALSE)))))
  if(dim(protInfo)[2] != 1){ #If the fasta has more than name
    protClass <- protInfo[,1]
    protNames <- protInfo[,2] 
    protFeatsDF <- data.frame(protClass,protNames, molecular_weight, pepLength, 
                              isoelectric_point, instability, aliphaticIndex, bomanIndex, 
                              AADataFrame)
  }else{
    protNames <- protInfo[,1] 
    protFeatsDF <- data.frame(protNames, molecular_weight, pepLength, 
                              isoelectric_point, instability, aliphaticIndex, bomanIndex, 
                              AADataFrame)
  }
  
  
  return(protFeatsDF)
}

# Function to make a statistic data.frame
presentStatistic <- function(DF) {
  col_names_to_test <- names(DF)[-c(1:2)]
  
  tTest_list <- list()
  for (col_name in col_names_to_test) {
    tTest_list[[col_name]] <- t.test(DF[[col_name]] ~ protClass, data = DF)
  }
  mwTest_list <- list()
  for (col_name in col_names_to_test) {
    mwTest_list[[col_name]] <- wilcox.test(DF[[col_name]] ~ protClass, data = DF)
  }
  test_list <- list("tTest_list" = tTest_list, "mwTest_list" = mwTest_list)
  
  # Define feature matrix - for each feature (row) we will put the results of
  # the statistical tests
  featureM <- matrix(nrow=dim(DF)[2]-2, ncol = 8)
  colnames(featureM) <- c("Allergen mean", "Non allergen mean", 
                          "Log2 means", "P value mean", "Allergen median",
                          "Non allergen median", "Log2 medians",
                          "P value median")
  i <- 1
  for (protType in unique(DF$protClass)) {
    featuresDF <- DF[DF$protClass == protType, -c(1:2)]
    meanV <- colMeans(featuresDF)
    featureM[, i] <- meanV
    i <- i + 1
    j <- 1
    medianV <- vector()
    for (feat in names(featuresDF)) {
      medianV[j] <- median(featuresDF[, feat])
      j <- j + 1
    }
    featureM[, i] <- medianV
    i <- i + 3
  }
  
  
  pValueMean <- vector(length=dim(featuresDF)[2])
  names(pValueMean) <- names(featuresDF)
  for (feat in names(featuresDF)) {
    pValueMean[feat] <- test_list[["tTest_list"]][[feat]][[3]]
  }
  
  pValueMedian <- vector(length=dim(featuresDF)[2])
  names(pValueMedian) <- names(featuresDF)
  for (feat in names(featuresDF)) {
    pValueMedian[[feat]] <- c(test_list[["mwTest_list"]][[feat]][[3]])
  }
  
  signMeans <- vector(length=dim(featuresDF)[2])
  for (i in 1:length(pValueMean)){
    if(pValueMean[i]<=0.05){
      signMeans[i] <- TRUE
    }else{
      signMeans[i] <- FALSE
    }
  }
  
  signMedians <- vector(length=dim(featuresDF)[2])
  for (i in 1:length(pValueMedian)){
    if(pValueMean[i]<=0.05){
      signMedians[i] <- TRUE
    }else{
      signMedians[i] <- FALSE
    }
  }
  
  featureM[, "P value mean"] <- pValueMean
  featureM[, "P value median"] <- pValueMedian
  log2Means <- log2(featureM[, "Allergen mean"]/featureM[, "Non allergen mean"])
  featureM[, "Log2 means"] <- log2Means
  log2Medians <- log2(featureM[, "Allergen median"]/featureM[, "Non allergen median"])
  featureM[, "Log2 medians"] <- log2Medians
  
  featureM <- cbind(featureM,signMeans,signMedians)
  featureM  <-featureM [, c(1,2,3,4,9,5,6,7,8,10)]
  
  row.names(featureM) <- names(featuresDF)
  
  featureM <- as.data.frame(featureM)
  names(featureM)[names(featureM) == "signMeans"] <- "Significative means"
  names(featureM)[names(featureM) == "signMedians"] <- "Significative medians"
  featureM <- format(featureM,digits=3,scientific = 4)

  return(featureM)
}
# Function to make density curves for each feature
makeDensityPlot <- function(DF, feature, xName) {
  DF$protClass <- as.factor(DF$protClass)
  return(ggplot(DF, aes(x = feature, color = protClass)) + geom_density()+ xlab(xName))
}
getAllDensityPlots <- function(DF) {
  col_names_to_plot <- names(DF[-c(1:2)])
  plots <- list()
  
  for (i in 1:length(col_names_to_plot)) {
    if (i == 1 | i == 2) {
      plots[[i]] <- eval(substitute(makeDensityPlot(DF,feature = log2(DF[[col_names_to_plot[i]]]), xName = col_names_to_plot[i]),list(i = i)))
    } else {
      plots[[i]] <- eval(substitute(makeDensityPlot(DF, feature = DF[[col_names_to_plot[i]]], xName = col_names_to_plot[i]),list(i = i)))
    }
  }
  return(plots)
}
