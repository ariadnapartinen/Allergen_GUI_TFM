#To search the allergen from WHO by biochemical names and uniprot ID.

# Initial stuff
library(rvest)
library(magrittr)
library(Biostrings)
library(seqRFLP)
library(xml2)

str_extract <- stringr::str_extract 

######################################################################################

search_allergens_WHO <- function(allergen_names) {

  alrgNames <- c() #empty array to use later
  links <- c()
  
  url <- "http://www.allergen.org/search.php?allergenname=&allergensource=&TaxSource=&TaxOrder=&foodallerg=all&bioname=" #list of the biochemical names to search in the allergen WHO website, separated by ";"
  
  split_allergen_names <- c(unlist(strsplit(allergen_names,";"))) #separates the different biochemical names from the list
  
  
  #### TO SEARCH THE LIST OF ALLERGEN PROTEIN REGISTERED IN THE FAMILY GIVEN AS INPUT ####
  
  for (i in 1:length(split_allergen_names)) { #loop to search in the website each biochemical name
    
    complete_url <- paste(url,split_allergen_names[i],sep="") #Concatenates the biochemical name to the url
    
    html <- read_html(complete_url) #To do the search in the website
    
    # Get the names for the different allergens from the table
    population <- html %>% html_nodes(xpath='//*[@id="example"]') %>% html_table %>% extract2(1) #Creates a table with the information from the website
    
    alrgNames <- unique(c(alrgNames, population[population[,1] == "", 2])) #Concatenates and extracts the allergen name (only the population's rows with the species value empty (the first row is the species repited))
    
    links <- union(links,paste0("http://www.allergen.org/", html %>% html_nodes(xpath = "//td/a") %>%
               html_attr("href"))) #concatenates and extracts the URLs of each entry (allergen)
    
  }
  
  # Follow these links to get the Uniprot IDs
  alrg_uniprot_list <- sapply(unlist(links),function(u) u %>% read_html %>% html_nodes(xpath = '//*[@id="isotable"]') %>% 
                                html_table() %>%  extract2(1) %>% extract("UniProt")) #create a list of allergen
  
  # Convert to vector and rename to allergen names, to do it extract the allergen name from the url lengths(alrg_uniprot_list)
  # and replace the url of uniprot for the name of WHO
  alrg_uniprot <- setNames(unlist(alrg_uniprot_list, use.names=F),rep(unlist(alrgNames), lengths(alrg_uniprot_list)))
  
  # Clean up NAN and empty name
  alrg_uniprot <- alrg_uniprot[! is.na(alrg_uniprot)] 
  alrg_uniprot <- alrg_uniprot[! alrg_uniprot == ""]
  
  alg_names <- list(alrgNames,alrg_uniprot)
  return(alg_names)
}

######################################################################

search_database_interpro <- function(interPro_accessions) {
  
  split_interPro_accessions <- c(unlist(strsplit(interPro_accessions,";")))
  
  #### TO SEARCH THE LIST OF non-ALLERGEN PROTEIN REGISTERED IN THE FAMILY GIVEN AS INPUT ####
  
  for (i in 1:length(split_interPro_accessions)) {
    inicial_url <- "http://www.ebi.ac.uk/interpro/entry/"
    end_url <-"/protein/UniProt/#table"
    complete_url <- paste(inicial_url, split_interPro_accessions[i],end_url,sep="")
    browseURL(complete_url) #Loads a given URL into an HTML browser
  }
  
}

######################################################################
 #To extract and bind the sequences from FASTA files chosen
extract_seq <- function(fnames){
 
  db <- data.frame(seq_accession = "seq_accession",seq_state = "seq_state",seq_name = "seq_name",seq_taxID = "seq_taxID",sequence = "sequence") #Create the db data frame

  for (i in 1:length(fnames)){
    seqs <- makeProtSeqDF(fnames[i]) #Separates the sequences of the proteins' information
    db <- rbind(db,seqs) #Concatenates all the sequences chosen in FASTA files
  }
  db <- db[-1,] #To delete the first row
  db <- unique(db) #To verify that any protein is repeated
  
  return(db)
}


