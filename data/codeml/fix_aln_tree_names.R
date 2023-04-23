
# PReparing data to run codeml - PAML
# This script fixes species names in the tree and alignment so that they match in order to run codeml

library(phytools)
library(bioseq) # https://fkeck.github.io/bioseq/reference/index.html
library(dplyr)
library(seqinr)

# Import tree and alignment
# The 'tree' is the ML tree from RAxML without branch lengths
tr <- read.tree("./Tconstrained_amniota_RAxML_bestTree.result_NoBL.tre")
tr$edge.length <- NULL
tr$root.edge <- NULL
aln <- read.fasta("../alignments/A1-A3_AA_alignment_convergence_modified.fasta")

tip_names <- tr$tip.label
seq_names <- aln$seq.name

# Fix species names
new_tip_names <- stringr::str_extract(tip_names,"(lungfish|fish|caecilian|frog|lizard|snake|bird|crocodilian|mammal|turtle)(.*?)_(A1|A2|A3)") %>% 
  gsub("Neophocaena_asiaeorientalis_asiaeorientalis","Neophocaena_asiaeorientalis",.)
new_seq_names <- stringr::str_extract(seq_names,"(lungfish|fish|caecilian|frog|lizard|snake|bird|crocodilian|mammal|turtle)(.*?)_(A1|A2|A3)") %>% 
  gsub("Neophocaena_asiaeorientalis_asiaeorientalis","Neophocaena_asiaeorientalis",.)

#check match
new_tip_names[!(new_tip_names %in% new_seq_names)]
new_seq_names[!(new_seq_names %in% new_tip_names)]

#change names in tree and alignment
tr2 <- tr
aln2 <- aln

tr2$tip.label <- new_tip_names
aln2$seq.name <- new_seq_names

# Export new fixed tree for codeml
phylotools::dat2phylip(aln2,outfile = "A1-A3_AA_alignment_convergence_modified_codeml.phy")
write.tree(tr2,"Tconstrained_amniota_RAxML_bestTree.result_NoBL_codeml.tre")

# check duplicate names and manually fix them
new_tip_names[duplicated(new_tip_names)]
