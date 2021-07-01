# Function to make data.frame from FASTA, separate the information of the sequence
makeProtSeqDF <- function(stringset) {
  first_DF <- readAAStringSet(stringset) #Creates an AAStringSet (container for storing AA sequences separated from names)
  seq_info_DF <- data.frame(strsplit(names(first_DF), split='[|]', fixed=FALSE))
  seq_info_transpose <- as.data.frame(t(as.matrix(seq_info_DF)))
  seq_accession <- seq_info_transpose$V1
  seq_state <- seq_info_transpose$V2
  seq_name<- seq_info_transpose$V3
  seq_taxID <- seq_info_transpose$V4
  sequence <- paste(first_DF)
  protSeqDF <- data.frame(seq_accession,seq_state,seq_name,seq_taxID, sequence, stringsAsFactors = FALSE)
  return(protSeqDF)
}
makeProtSeqDF2 <- function(stringset) {
  first_DF <- readAAStringSet(stringset) #Creates an AAStringSet (container for storing AA sequences separated from names)
  seq_accession <- names(first_DF)
  sequence <- paste(first_DF)
  protSeqDF <- data.frame(seq_accession, sequence, stringsAsFactors = FALSE)
  return(protSeqDF)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Function to normalize the data with the mean and standard deviation of train data
scaleData <-function(data,train){ 
  m <- apply(train[,-c(1:2)],2,mean)
  sd <- apply(train[,-c(1:2)],2,sd)
  s <- data.frame(matrix(ncol = ncol(data), nrow = 0))
  colnames(s) <- names(data)
  k <- dim(data)[2] -14
  for (i in k:dim(data)[2]) { #To go through all the features
    for (j in 1:dim(data)[1]) { #To go through all the data which we want to normalize.
      if(k==2){
        s[j,i] <- (data[j,i]-m[i-1])/sd[i-1]
      }else{
        s[j,i-1] <- (data[j,i]-m[i-2])/sd[i-2] 
      }
    }
  }
  if(k==2){
    s[,1] <- data[,1]
  }else{
    s[,1:2] <- data[,1:2]
  }
  return(s) 
}