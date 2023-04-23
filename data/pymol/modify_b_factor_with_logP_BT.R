
# Script to modify the B factor column of a pdb file with the logP value of the BayesTraits result

pdb_filename <- './3b8e_chainA.pdb'
new_pdb_filename111 <- './3b8e_chainA_New_B_Factor_logP_111.pdb'
new_pdb_filename122 <- './3b8e_chainA_New_B_Factor_logP_122.pdb'

#read files
data111_BT <- read.table("../BayesTraits/data111_BT_restrictedModels_allParalogs.txt",header =T)
data122_BT <- read.table("../BayesTraits/data122_BT_restrictedModels_allParalogs.txt",header =T)


modify_b_factor <- function(pdb_filename,BTfile){
  #read pdb file
  pdb_file <- readLines(pdb_filename, warn=FALSE)
  
  # modify b-factor column in pdb_file, iterating through pdb_file line by line
  for(i in 1:length(pdb_file)) {
    if(grepl('ATOM', pdb_file[i]) == TRUE) {
      
      # get residue number for current line
      residue <- substr(pdb_file[i], 23, 26)
      residue <- as.integer(residue)
      
      # check if residue is in the site column in the sitewise_diversity file
      if(residue %in% BTfile$site) {
        
        # get the sitewise amino acid diversity for the current residue
        logP <- dplyr::filter(BTfile,site==residue) %>% select(.,logP)
        logP <- logP[1,]
        
        # replace current data in b-factor column with site.count
        substr(pdb_file[i], 61, 66) <- paste0('  ', format(round(logP, 2), nsmall=2))

        # replace remaining b-factor columns with NA
      } else {
        substr(pdb_file[i], 61, 66) <- '    NA'
      }
      
    } else if(grepl('HETATM', pdb_file[i]) == TRUE) {
      substr(pdb_file[i], 61, 66) <- '    NA'
    }
  }
  return(pdb_file)
}

pdb_111 <- modify_b_factor(pdb_filename,data111_BT)
pdb_122 <- modify_b_factor(pdb_filename,data122_BT)

# write modified pdb data to new file
writeLines(pdb_111, new_pdb_filename111)
writeLines(pdb_122, new_pdb_filename122)


