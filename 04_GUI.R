
rm(list=ls()) #clear the environment

library(caTools)
library(dplyr) 
library(rvest)
library(magrittr)
library(Biostrings)
library(seqRFLP)
library(xml2)
library(xlsx)
library(ggplot2)


source("00_Functions.R")  
source("01_AllergenSearch.R")
source("02_FeaturesExtraction.R")
source("03_Classification.R")


library(gWidgets2)
library(gWidgets2RGtk2)

options(guiToolkit="RGtk2")

str_extract <- stringr::str_extract  

#CREATE THE GUI
if(interactive()) {
  
  w <- gwindow("Allergen Classification", visible=FALSE) #Title of the window
  g <- ggroup(cont=w, horizontal=FALSE) #Provides a vertical box container for packing in child components
  
  # Declaration of the different components
  
  l1 <- glabel("Allergen family to search in WHO database", container=g)
  e1<- gedit(text = "Enter the biochemical names separated by ; ", width = 25, cont=g)
  b1 <- gbutton("Search", cont=g, handler=function(...) {
    allergen_names <<-search_allergens_WHO(svalue(e1)) #Extracts the allergen names and its links
    l3$widget$setText(length(unlist(allergen_names[2]))) #Shows the number of allergens found
    })
  
  l2 <- glabel("Number of allergens found", container=g)
  l3 <- glabel("", container=g)
  
  
  #Pop up window to save the list of ID's allergen found as csv file
  b2 <- gbutton("Save Allergen IDs List", cont=g, handler=function(...) {
   subw <- gwindow("Save Allergen IDs List", visible=FALSE) 
   subg <- ggroup(cont=subw, horizontal=FALSE) 
    l4 <- glabel("Path to save", container=subg)
    e2<- gedit(text = "", width = 25, cont=subg)
    l5 <- glabel("File name", container=subg)
    e3<- gedit(text = "", width = 25, cont=subg)
    b21 <- gbutton("Save", cont=subg, handler=function(...) {
      if(svalue(e2)!="" & svalue(e3)!="") #If the user enter both values
      path = paste(svalue(e2), "\\", svalue(e3),".csv", sep="")
      write.csv(allergen_names[2],path,row.names =FALSE)
      visible(subw) <- FALSE
    })
    b22 <- gbutton("Cancel", cont=subg, handler=function(...) {
      visible(subw) <- FALSE
    })
   visible(subw) <- TRUE
  })
  
  l6 <- glabel("Protein family to search in InterPro database", container=g) #Only in case of do not have previously downloaded the FASTA file
  e4<- gedit(text = "Enter the accessions separated by ; ", width = 25, cont=g)
  b3 <- gbutton("Load the given accessions into a browser", cont=g, handler=function(...) {
    search_database_interpro(svalue(e4))
  })
  
  #To open the FASTA files downloaded before
  b4 <- gbutton("Open InterPro FASTA files", cont=g, handler=function(...) {
    fnames <- choose.files()
    db <<-extract_seq(fnames)
    l8$widget$setText(dim(db)[1]) ##Shows the number of proteins found
    allergen_names_DF <<- data.frame(allergen_names[2])
    seq_alg <<- subset(db, db$seq_accession %in% allergen_names_DF[,1]) #Extracts the sequences of the allergen list found in WHO
    num_alg <<- dim(seq_alg)[1]
    seq_alg$seq_accession <<- paste(as.character(rep(1,num_alg)),seq_alg$seq_accession, sep = '|')
    num <<- num_alg*10 #To have a 1:10 proportion
    seq_nalg_total <<- db[ !(db$seq_accession %in% seq_alg$seq_accession), ] #Creates a non-allergen list from the allergen list found in WHO
    num_nalg <<- dim(seq_nalg_total)[1]
    seq_nalg_total$seq_accession <<- paste(as.character(rep(0,num_nalg)),seq_nalg_total$seq_accession, sep = '|')
    seq_nalg <<- seq_nalg_total[sample(nrow(seq_nalg_total), num), ] #Extracts the desired number of proteins randomly.
    l10$widget$setText(dim(seq_alg)[1])
    l12$widget$setText(dim(seq_nalg_total)[1])
  }) 
  
  l7 <- glabel("Number of proteins found in FASTA files", container=g)
  l8 <- glabel("", container=g)
  l9 <- glabel("Number of ALLERGEN proteins found in FASTA files", container=g)
  l10 <- glabel("", container=g)
  l11 <- glabel("Number of NON ALLERGEN proteins found in FASTA files", container=g)
  l12 <- glabel("", container=g)
  
  #Pop up window to save the list of ID's allergen found as csv file
  b5 <- gbutton("Save sets found as FASTA files", cont=g, handler=function(...) {
    subw <- gwindow("Save sets found as FASTA files", visible=FALSE) 
    subg <- ggroup(cont=subw, horizontal=FALSE) 
    l13 <- glabel("Path to save", container=subg)
    e5<- gedit(text = "", width = 25, cont=subg)
    l14 <- glabel("ALLERGEN file name", container=subg)
    e6<- gedit(text = "", width = 25, cont=subg)
    l15 <- glabel("NON-ALLERGEN file name", container=subg)
    e7<- gedit(text = "", width = 25, cont=subg)
    l16 <- glabel("Subset NON-ALLERGEN file name", container=subg)
    e8<- gedit(text = "", width = 25, cont=subg)
    b51 <- gbutton("Save", cont=subg, handler=function(...) { #Saves the set which you want
      if(svalue(e5)!="" & svalue(e6)!="") {#If the user enter both values TO SAVE ALG SET save it
        alg_fa <- dataframe2fas(seq_alg %>% select(seq_accession,sequence))
        path <- paste(svalue(e5), "\\", svalue(e6),".fa", sep="")
        write.fasta(alg_fa, path)
        visible(subw) <- FALSE
    } 
      if(svalue(e5)!="" & svalue(e7)!="") {#If the user enter both values TO SAVE NON ALG TOTAL SET save it
        nalg_total_fa <- dataframe2fas(seq_nalg_total %>% select(seq_accession, sequence))
        path = paste(svalue(e5), "\\", svalue(e7),".fa", sep="")
        write.fasta(nalg_total_fa, path)
        visible(subw) <- FALSE
    }  
      if(svalue(e5)!="" & svalue(e8)!="") {#If the user enter both values TO SAVE NON ALG SUBSET save it
        nalg_fa <- dataframe2fas(seq_nalg %>% select(seq_accession, sequence))
        path <- paste(svalue(e5), "\\", svalue(e8),".fa", sep="")
        write.fasta(nalg_fa, path)
        visible(subw) <- FALSE
      }    
    })
    b52 <- gbutton("Cancel", cont=subg, handler=function(...) {
      visible(subw) <- FALSE
    })
    visible(subw) <- TRUE
  })
  
  #Windows to extract and show the features and statistical analysis
  b6 <- gbutton("Extract features", cont=g, handler=function(...) {
    #We must differentiate between if you are doing all the process (search allergens in WHO, etc.) or if you ...
    #have the fasta files of allergens and non allergens, that is, you already have the subsets or not.
    if(exists("seq_alg")){
      alg_features <<- getProtFeatures(protSeqDF = seq_alg)
    }else {
      fnames <- choose.files(caption = "Select allergen fasta file", multi = FALSE)
      seq_alg <<- makeProtSeqDF2(fnames)
      alg_features <<- getProtFeatures(protSeqDF = seq_alg)
    }
    if(exists("seq_nalg")){
      nalg_features <<- getProtFeatures(protSeqDF = seq_nalg)
    }else {
      fnames <- choose.files(caption = "Select subset non allergen fasta file", multi = FALSE)
      seq_nalg <<- makeProtSeqDF2(fnames)
      nalg_features <<- getProtFeatures(protSeqDF = seq_nalg)
    }
    
    DF <<- rbind(alg_features, nalg_features)
    subw <- gwindow("Extracted features", visible=FALSE)
    subg <- gvbox(cont=subw)
    tbl <- gtable(DF, cont=subg, expand=TRUE, fill=TRUE)
    addHandlerClicked(tbl, handler=function(h,...) sprintf("You selected %s", svalue(h$obj)))
    
    b61 <- gbutton("Export table", cont=subg, handler=function(...) {
      w2 <- gwindow("Export table", visible=FALSE) 
      g2 <- ggroup(cont=w2, horizontal=FALSE)
      
      fl <- gformlayout(cont=g2)
      formats<- c("CSV","TXT","EXCEL")
      gformat <-gradio(formats, selected=1, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
      
      l4 <- glabel("Path to save", container=g2)
      e2<- gedit(text = "", width = 25, cont=g2)
      l5 <- glabel("File name", container=g2)
      e3<- gedit(text = "", width = 25, cont=g2)
      
      b611 <- gbutton("Save", cont=g2, handler=function(...) {
        if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
          
          if(svalue(gformat) == "TXT"){
            path <- paste(svalue(e2), "\\", svalue(e3),".txt", sep="")
            write.table(DF,path,sep="\t",row.names=FALSE)
          }
          if (svalue(gformat) == "EXCEL"){
            path <- paste(svalue(e2), "\\", svalue(e3),".xlsx", sep="")
            write.xlsx(DF,path,col.names = TRUE, row.names = FALSE)
          }
          if(svalue(gformat) == "CSV"){
            path <- paste(svalue(e2), "\\", svalue(e3),".csv", sep="")
            write.csv(DF,path,row.names =FALSE)
          }
          visible(w2) <- FALSE
        }
      })
      
      b612 <- gbutton("Cancel", cont=g2, handler=function(...) {
        visible(w2) <- FALSE
      })
      visible(w2) <- TRUE
    })
    
    b62 <- gbutton("Statistical analysis", cont=subg, handler=function(...) {
      featureAnalysis <<- presentStatistic(DF)
      rowNames <- as.data.frame(rownames(featureAnalysis))
      colnames(rowNames) <- c("Features")
      tblFeatureAnalysis <- cbind(rowNames, featureAnalysis)
      
      subw <- gwindow("Statistical analysis", visible=FALSE)
      subg <- gvbox(cont=subw)
      
      tbl <- gtable(tblFeatureAnalysis, cont=subg, expand=TRUE, fill=TRUE)
      
      b621 <- gbutton("Export table", cont=subg, handler=function(...) {
        w2 <- gwindow("Export table", visible=FALSE) 
        g2 <- ggroup(cont=w2, horizontal=FALSE) 
        
        fl <- gformlayout(cont=g2)
        formats<- c("CSV","TXT","EXCEL")
        gformat <-gradio(formats, selected=1, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
        
        
        l4 <- glabel("Path to save", container=g2)
        e2<- gedit(text = "", width = 25, cont=g2)
        l5 <- glabel("File name", container=g2)
        e3<- gedit(text = "", width = 25, cont=g2)
        
        b6211 <- gbutton("Save", cont=g2, handler=function(...) {
          if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
            if(svalue(gformat) == "TXT"){
              path = paste(svalue(e2), "\\", svalue(e3),".txt", sep="")
              write.table(featureAnalysis,path,sep="\t",row.names=TRUE)
            }
            if (svalue(gformat) == "EXCEL"){
              path = paste(svalue(e2), "\\", svalue(e3),".xlsx", sep="")
              write.xlsx(featureAnalysis,path,col.names = TRUE, row.names = TRUE)
            }
            if(svalue(gformat) == "CSV"){
              path = paste(svalue(e2), "\\", svalue(e3),".csv", sep="")
              write.csv(featureAnalysis,path,row.names =TRUE)
            }
            visible(w2) <- FALSE
          }
        })
        
        b6212 <- gbutton("Cancel", cont=g2, handler=function(...) {
          visible(w2) <- FALSE
        })
        visible(w2) <- TRUE
      })
      
      b622 <- gbutton("Show density plots", cont=subg, handler=function(...) {
      
        w3 <- gwindow("Density plots",visible = FALSE)
        g3 <- ggroup(cont=w3, horizontal=FALSE) 
        
        gg <- ggraphics(cont=g3,label = "Density plots",visible=FALSE)
        plots <<- getAllDensityPlots(DF)
        
        visible(w3) <- TRUE
        visible(gg) <- TRUE
        
        multiplot(plotlist = plots, cols = 3)
        
        b6221 <- gbutton("Export", cont=g3, handler=function(...) {
          w4 <- gwindow("Export plots", visible=FALSE) 
          g4 <- ggroup(cont=w4, horizontal=FALSE) 
          fl <- gformlayout(cont=g4)
          formats<- c("PDF","PNG","EPS")
          gformat <-gradio(formats, selected=2, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
          visible(w4) <- TRUE
          
          
          l4 <- glabel("Path to save", container=g4)
          e2<- gedit(text = "", width = 25, cont=g4)
          l5 <- glabel("File name", container=g4)
          e3<- gedit(text = "", width = 25, cont=g4)
          
          b62211 <- gbutton("Save", cont=g4, handler=function(...) {
            if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
              if(svalue(gformat) == "PDF"){
                path <- paste(svalue(e2), "\\", svalue(e3),".pdf", sep="")
                pdf(path)
                for (i in 1:length(plots)) {
                  print(plots[[i]])
                }
                dev.off() 
              }
              if (svalue(gformat) == "PNG"){
                path <- paste(svalue(e2), "\\", svalue(e3),".png", sep="")
                png(path, units="px", width=1600, height=1100, res=160)
                multiplot(plotlist = plots, cols = 3)
                dev.off() 
              }
              if(svalue(gformat) == "EPS"){
                path <- paste(svalue(e2), "\\", svalue(e3),".eps", sep="")
                setEPS()
                postscript(path,width = 8,height =12)
                multiplot(plotlist = plots, cols = 2)
                dev.off() 
              }
              visible(w4) <- FALSE
            }
          
        })
          
          b62212 <- gbutton("Cancel", cont=g4, handler=function(...) {
            visible(w4) <- FALSE
          })
          visible(w4) <- TRUE
          
        })
        b6222 <- gbutton("Close", cont=g3, handler=function(...) {
          visible(w3) <- FALSE
        })
      
      })
      b623 <- gbutton("Close", cont=subg, handler=function(...) {
        visible(subw) <- FALSE
      })
      
      visible(subw) <- TRUE
    
    })   
    
    b63 <- gbutton("Close", cont=subg, handler=function(...) {
      visible(subw) <- FALSE
    })
    
    visible(subw) <- TRUE
  })
    
  #Windows to classify the dataset 
  b7 <- gbutton("Classification", cont=g, handler=function(...) {
      subw <- gwindow("Classification", visible=FALSE) 
      subg <- ggroup(cont=subw, horizontal=FALSE) 
      
      fl <- gformlayout(cont=subg)
      features<- c("All features","Features significantly different")
      gfeatures <-gradio(features, selected=1, horizontal=TRUE, cont=fl,index = TRUE, label="Select a input data   ")
      
      fl2 <- gformlayout(cont=subg)
      t<- c("Performance test","Classify new data")
      gtest <-gradio(t, selected=1, horizontal=TRUE, cont=fl2,index = TRUE, label="Test samples   ")
      #You can choose if you want to extract the test sample(s) from the database just searched or if you have some
      #sample(s) to classify only search its protein family and classify.
      
      fl3 <- gformlayout(cont=subg)
      numFolds<- c("5","10","15")
      gnumFolds <-gradio(numFolds, selected=2, horizontal=TRUE, cont=fl2,index = TRUE, label="Number of folds in cross validation   ")
    
      fl4<- gformlayout(cont=subg)
      classifiers<- c("Decision Tree","KNN","MLP","NB Gaussian","SVM")
      gclassifiers <-gradio(classifiers, selected=5, horizontal=TRUE, cont=fl2,index = TRUE, label="Classifier   ")
      
      
      if(!exists("DF")){
          w2 <- gwindow("Message", visible=FALSE) 
          g2 <- ggroup(cont=w2, horizontal=FALSE) 
          
          l4 <- glabel("Please, load allergen and non allergen data in \"Extract features\"", container=g2)
          b7111 <- gbutton("Ok", cont=g2, handler=function(...) {
            visible(w2) <- FALSE
          })
          visible(w2)<-TRUE
      }
      
      if(!exists("featureAnalysis")){
        featureAnalysis <<- presentStatistic(DF)
      }
      
      b71 <- gbutton("Classify", cont=subg, handler=function(...) { #Select features as input according to the previous selection,
        #normalizes the data and create K=10 subsets randomly to use in cross validation
        
        if(svalue(gtest) == "Classify new data"){
          fnames <- choose.files(caption = "Select test fasta file", multi = FALSE) #Selects a fasta file
          seq_test <<- makeProtSeqDF2(fnames) #Generates a dataframe from this fasta file
          test_set <<- getProtFeatures(protSeqDF = seq_test) #Extracts the features of this sequences
          training_set <<- subset(DF, !(DF$protNames %in% test_set$protNames)) #In this case, the train data is all samples 
          #of the database except the test set coincidences (no matter how many samples are).
        }else{
          #In this case, test set is a proportion (10%) of the database.
          set.seed(29)
          split <- sample.split(DF$protClass, SplitRatio = 0.90)#To extract randomly 10% of each class
          training_set <<- subset(DF, split == TRUE)
          test_set <<- subset(DF, split == FALSE)
        }
        
        if(svalue(gfeatures) == "Features significantly different"){
          non_signif_features <- subset(featureAnalysis, featureAnalysis$`Significative means` == 0 & featureAnalysis$`Significative medians` == 0)
          test_set <<- test_set[, !names(test_set) %in% names(test_set[rownames(non_signif_features)])]
          training_set <<- training_set[, !names(training_set) %in% names(training_set[rownames(non_signif_features)])]
        }
        
       normalized_test_set <<- scaleData(data = test_set,train = training_set)
       normalized_training_set <<- scaleData(data = training_set,train = training_set)
       kfolds <<- as.numeric(svalue(gnumFolds))
       
       w2 <- gwindow("Training", visible=FALSE) 
       g2 <- ggroup(cont=w2, horizontal=FALSE) 
       
       gg <- ggraphics(cont=g2,label = "Results",visible=FALSE)
       
       mdl <<- train_cv(kfolds = kfolds,classifier = svalue(gclassifiers), training_set = normalized_training_set)
       
       visible(w2) <- TRUE
       visible(gg) <- TRUE
       
       trainResults <<- train_results(mdl,kfolds)
       print(trainResults)
       
       b7111 <- gbutton("Export figure", cont=g2, handler=function(...) {
         w3 <- gwindow("Export figure", visible=FALSE) 
         g3 <- ggroup(cont=w3, horizontal=FALSE) 
         fl <- gformlayout(cont=g3)
         formats<- c("PDF","PNG","EPS")
         gformat <-gradio(formats, selected=2, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
         visible(w3) <- TRUE
         
         
         l4 <- glabel("Path to save", container=g3)
         e2<- gedit(text = "", width = 25, cont=g3)
         l5 <- glabel("File name", container=g3)
         e3<- gedit(text = "", width = 25, cont=g3)
         
         b71111 <- gbutton("Save", cont=g3, handler=function(...) {
           if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
             if(svalue(gformat) == "PDF"){
               path <- paste(svalue(e2), "\\", svalue(e3),".pdf", sep="")
               pdf(path)
               print(trainResults)
               dev.off() 
             }
             if (svalue(gformat) == "PNG"){
               path <- paste(svalue(e2), "\\", svalue(e3),".png", sep="")
               png(path, units="px", width=1600, height=1100, res=160)
               print(trainResults)
               dev.off() 
             }
             if(svalue(gformat) == "EPS"){
               path <- paste(svalue(e2), "\\", svalue(e3),".eps", sep="")
               setEPS()
               postscript(path,width = 8,height =12)
               print(trainResults)
               dev.off() 
             }
             visible(w3) <- FALSE
           }
           
         })
         
         b71112 <- gbutton("Cancel", cont=g3, handler=function(...) {
           visible(w3) <- FALSE
         })
         
       })
       
       b7112 <- gbutton("Close", cont=g2, handler=function(...) {
          visible(w2) <- FALSE
       })
       
       if(svalue(gtest) == "Classify new data"){
         results <<- classify(model = mdl, test = normalized_test_set)
         
         subw <- gwindow("Results", visible=FALSE)
         subg <- gvbox(cont=subw)
         tbl <- gtable(results, cont=subg, expand=TRUE, fill=TRUE)
         addHandlerClicked(tbl, handler=function(h,...) sprintf("You selected %s", svalue(h$obj)))
         
         b71121 <- gbutton("Export table", cont=subg, handler=function(...) {
           w2 <- gwindow("Export table", visible=FALSE) 
           g2 <- ggroup(cont=w2, horizontal=FALSE)
           
           fl <- gformlayout(cont=g2)
           formats<- c("CSV","TXT","EXCEL")
           gformat <-gradio(formats, selected=1, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
           
           l4 <- glabel("Path to save", container=g2)
           e2<- gedit(text = "", width = 25, cont=g2)
           l5 <- glabel("File name", container=g2)
           e3<- gedit(text = "", width = 25, cont=g2)
           
           b711211 <- gbutton("Save", cont=g2, handler=function(...) {
             if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
               
               if(svalue(gformat) == "TXT"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".txt", sep="")
                 write.table(results,path,sep="\t",row.names=FALSE)
               }
               if (svalue(gformat) == "EXCEL"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".xlsx", sep="")
                 write.xlsx(results,path,col.names = TRUE, row.names = FALSE)
               }
               if(svalue(gformat) == "CSV"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".csv", sep="")
                 write.csv(results,path,row.names =FALSE)
               }
               visible(w2) <- FALSE
             }
           })
           
           
           b711212 <- gbutton("Cancel", cont=g2, handler=function(...) {
             visible(w2) <- FALSE
           })
           visible(w2) <- TRUE
         })
         
         b71122 <- gbutton("Close", cont=subg, handler=function(...) {
           visible(subw) <- FALSE
         })
         visible(subw) <- TRUE
         
       }else{
         confm<<- classify(model = mdl, test = normalized_test_set)
         w3 <- gwindow("Performance test", visible=FALSE) 
         g3 <- ggroup(cont=w3, horizontal=FALSE) 
         
         gg <- ggraphics(cont=g3,label = "Results",visible=FALSE)
         
         visible(w3) <- TRUE
         visible(gg) <- TRUE
         
         draw_confusion_matrix(confm,svalue(gclassifiers),kfolds)
         
         b71123 <- gbutton("Export figure", cont=g3, handler=function(...) {
           w4 <- gwindow("Export figure", visible=FALSE) 
           g4 <- ggroup(cont=w4, horizontal=FALSE) 
           fl <- gformlayout(cont=g4)
           formats<- c("PDF","PNG","EPS")
           gformat <-gradio(formats, selected=2, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
           visible(w4) <- TRUE
           
           
           l4 <- glabel("Path to save", container=g4)
           e2<- gedit(text = "", width = 25, cont=g4)
           l5 <- glabel("File name", container=g4)
           e3<- gedit(text = "", width = 25, cont=g4)
           
           b71111 <- gbutton("Save", cont=g4, handler=function(...) {
             if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
               if(svalue(gformat) == "PDF"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".pdf", sep="")
                 pdf(path)
                 draw_confusion_matrix(confm,svalue(gclassifiers),kfolds)
                 dev.off() 
               }
               if (svalue(gformat) == "PNG"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".png", sep="")
                 png(path, units="px", width=1600, height=1100, res=160)
                 draw_confusion_matrix(confm,svalue(gclassifiers),kfolds)
                 dev.off() 
               }
               if(svalue(gformat) == "EPS"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".eps", sep="")
                 setEPS()
                 postscript(path,width = 8,height =12)
                 draw_confusion_matrix(confm,svalue(gclassifiers),kfolds)
                 dev.off() 
               }
               
             }
             visible(w4) <- FALSE
             
           })
           
           b612 <- gbutton("Cancel", cont=g4, handler=function(...) {
             visible(w4) <- FALSE
           })
           
           visible(w4) <- TRUE
         })
         
         b71124 <- gbutton("Export results", cont=g3, handler=function(...) {
           w4 <- gwindow("Export results", visible=FALSE) 
           g4 <- ggroup(cont=w4, horizontal=FALSE) 
           
           l4 <- glabel("Path to save", container=g4)
           e2<- gedit(text = "", width = 25, cont=g4)
           l5 <- glabel("File name", container=g4)
           e3<- gedit(text = "", width = 25, cont=g4)
           
           fl <- gformlayout(cont=g4)
           formats<- c("CSV","EXCEL")
           gformat <-gradio(formats, selected=1, horizontal=TRUE, cont=fl,index = TRUE, label="Select a format         ")
           
           res <- as.numeric(confm$table)
           TP <- res[4] 
           TN <- res[1] 
           FP <- res[2] 
           FN <- res[3]
           
           if (TP+FN == 0 || TP+FP == 0 || TN+FP == 0 || TN+FN == 0){
             MCC <- TP*TN-FP*FN
           }else{
             MCC <- (TP*TN-FP*FN)/(sqrt((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN)))
           }
           
           t_cm <- data.frame(rbind(TN,FN,FP,TP,MCC))
           names(t_cm)[1] <- "Metrics"
           t_results <- data.frame(confm$byClass) 
           names(t_results)[1] <- "Metrics"
           test_results <<- round(rbind(t_cm,t_results),3)
           
           
           b71121 <- gbutton("Save", cont=g4, handler=function(...) {
             if(svalue(e2)!="" & svalue(e3)!="") {#If the user enter both values
               if (svalue(gformat) == "EXCEL"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".xlsx", sep="")
                 write.xlsx(test_results,path,col.names = TRUE, row.names = TRUE)
               }
               if(svalue(gformat) == "CSV"){
                 path <- paste(svalue(e2), "\\", svalue(e3),".csv", sep="")
                 write.csv(test_results,path,row.names =TRUE)
               }
               visible(w4) <- FALSE
             }
           })
           
           b612 <- gbutton("Cancel", cont=g4, handler=function(...) {
             visible(w4) <- FALSE
           })
           
           visible(w4) <- TRUE
         })
         
         b71125 <- gbutton("Close", cont=g3, handler=function(...) {
           visible(w3) <- FALSE
         })
         
       }
       
      })
      
      b72 <- gbutton("Close", cont=subg, handler=function(...) {
        visible(subw) <- FALSE
      })
      
      visible(subw) <- TRUE
      
      })
  
  b8 <- gbutton("Close", cont=g, handler=function(...) {
    visible(w) <- FALSE
  })
  
  visible(w) <- TRUE
  
  
}
