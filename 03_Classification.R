
#library(tidyverse)
library(rpart)
library(RSNNS)
library(naivebayes)
library(LogicReg)
library(caret)
library(MLmetrics)

#Function to normalize the data with the mean and standard deviation of train data
scaleData <-function(data,train){ 
  m <- apply(train[,-c(1:2)],2,mean)
  sd <- apply(train[,-c(1:2)],2,sd)
  s <- data.frame(matrix(ncol = ncol(data), nrow = 0))
  colnames(s) <- names(data)
  k <- grep("protNames", colnames(data)) #If you have class and name it starts at 3 and if you only have the name it starts at 2
  for (i in 1:length(m)) { 
    for (j in 1:dim(data)[1]) {
      s[j,i+k] <- (data[j,i+k]-m[i])/sd[i]
    }
  }
  if(k==1){
    s[,1] <- data[,1]
  }else{
    s[,1:2] <- data[,1:2]
  }
  return(s) 
}

#Cross validation k-Folds

train_cv <- function(kfolds, classifier, training_set){

  control<- trainControl(method="cv", number=kfolds, savePredictions = "final")
  metric <- "Accuracy"
  set.seed(29)
  
  if(classifier == "Decision Tree"){
    mdl <- train(protClass~., data=training_set[-c(2)], method="rpart", metric=metric, trControl=control)
  }
  
  if(classifier == "KNN"){
    mdl <- train(protClass~., data=training_set[-c(2)], method="knn", metric=metric, trControl=control)
  }
  
  if(classifier == "MLP"){
    mdl <- train(protClass~., data=training_set[-c(2)], method="mlp", metric=metric, trControl=control)
  }
  
  if(classifier == "NB Gaussian"){
    mdl <- train(protClass~., data=training_set[-c(2)], method="naive_bayes", metric=metric, trControl=control)
  }
  
  if(classifier == "SVM"){
    mdl <- train(protClass~., data=training_set[-c(2)], method="svmLinear", metric=metric, trControl=control)
  }
  
  return(mdl)
}

#To plot the validation results in each iteration and its average
train_results <- function(model,kfolds){
  p <- model$pred
  acc_train <- data.frame(matrix(ncol = 2, nrow = kfolds))
  
  if (kfolds == 5){
    for (i in 1:kfolds) {
      fold <- paste("Fold",i, sep = "")
      y_pred <- p$pred[p$Resample == fold]
      y_true <- p$obs[p$Resample == fold]
      acc_train[i,2] <- round(Accuracy(y_pred, y_true)*100, 2)
    }
  }else{
  
    for (i in 1:kfolds) {
      if(i < 10){
        fold <- paste("Fold0",i, sep = "")
      }else{
        fold <- paste("Fold",i, sep = "")
      }
     
      y_pred <- p$pred[p$Resample == fold]
      y_true <- p$obs[p$Resample == fold]
      acc_train[i,2] <- round(Accuracy(y_pred, y_true)*100, 2)
    }
  }
  acc_train[,1] <- c(1:kfolds)
  colnames(acc_train) <- c('folds','acc')
  return(ggplot(acc_train, aes(x = folds, y = acc)) + geom_point()+ geom_line()+ xlab('Folds') + ylab('Accuracy (%)') +
                                               geom_hline(yintercept = mean(acc_train$acc), color="blue"))
  
}

#Classifies and generates a confusion matrix if it is a performance test or generates a data.frame with the results if it is a test of samples
classify <- function(model, test){
  predicted <- predict(model, test)
  
  if("protClass" %in% colnames(test)) #To test and create a confusion matrix
  {
    truth <- factor(test$protClass)
    cm <- confusionMatrix(truth, predicted, positive = '1') #1 means allergen
    return(cm)

  }else  #To test new data (without the actually class)
  { 
    scores <- predict(model, test,type = "prob")
    predClass <- data.frame(matrix(ncol = 1, nrow = length(predicted)))
    colnames(predClass) <- c("Class")
    
    probClass <- data.frame(matrix(ncol = 1, nrow = length(predicted)))
    colnames(probClass) <- c("Probability")
    probClass[,1] <-  round(pmax(scores[,1],scores[,2])*100,digits = 2)
    
    for (i in 1:length(predicted)) {
      if(predicted[i]==1){
        predClass[i,1] <- "Allergen"
      }else{
        predClass[i,1] <- "Non-allergen"
      }
    }
    results <- data.frame(test$protNames,predClass,probClass)
    return(results)
  }
  
}

draw_confusion_matrix <- function(cm,classifier,kfolds) {
  
  total <- sum(cm$table)
  res <- as.numeric(cm$table)
  TP <- res[4] 
  TN <- res[1] 
  FP <- res[2] 
  FN <- res[3]
  
  if (TP+FN == 0 || TP+FP == 0 || TN+FP == 0 || TN+FN == 0){
    mcc <- TP*TN-FP*FN
  }else{
    mcc <- (TP*TN-FP*FN)/(sqrt((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN)))
  }
  
  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }
  
  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(285, 440), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  tit <- toupper(paste(classifier, "with", as.character(kfolds),"folds", sep = " "))
  title(tit, cex.main=2)
  
  # create the matrix 
  #classes = colnames(cm$table)
  classes <- c("Non allergen", "Allergen")
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  #text(195, 435, classes[1], cex=1.2)
  text(195, 295, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  #text(295, 435, classes[2], cex=1.2)
  text(295, 295, classes[2], cex=1.2)
  text(125, 370, 'Predicted Class', cex=1.3, srt=90, font=2)
  #text(245, 450, 'True class', cex=1.3, font=2)
  text(245, 285, 'True class', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)
  
  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(20, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(20, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(50, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(50, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
  text(80, 35, "MCC", cex=1.5, font=2)
  text(80, 20, round(mcc, 3), cex=1.4)
}