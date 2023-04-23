
# Create backbone tree to run RAxML
# Constrain topology using main clades of tetrapods

library(phytools)
library(bioseq) # https://fkeck.github.io/bioseq/reference/index.html
library(dplyr)
library(seqinr)

#Use alignment seq.names to build a backbone phylogeny
seqs <- as_tibble.AAbin(read_fasta("../alignments/A1-A3_AA_alignment_convergence_modified_codeml.fasta",type="AA"))

sp_names <- seqs$label

# Extract clades for each ATP1A paralog

# A1
frogA1 <- sp_names[grepl('frog',sp_names) & grepl('_A1',sp_names)]
caecilianA1 <- sp_names[grepl('caecilian',sp_names) & grepl('_A1',sp_names)]
lizardA1 <- sp_names[grepl('lizard',sp_names) & grepl('_A1',sp_names)]
snakeA1 <- sp_names[grepl('snake',sp_names) & grepl('_A1',sp_names)]
turtleA1 <- sp_names[grepl('turtle',sp_names) & grepl('_A1',sp_names)]
mammalA1 <- sp_names[grepl('mammal',sp_names) & grepl('_A1',sp_names)]
birdA1 <- sp_names[grepl('bird',sp_names) & grepl('_A1',sp_names)]
crocodilianA1 <- sp_names[grepl('crocodilian',sp_names) & grepl('_A1',sp_names)]

# A1
frogA2 <- sp_names[grepl('frog',sp_names) & grepl('_A2',sp_names)]
caecilianA2 <- sp_names[grepl('caecilian',sp_names) & grepl('_A2',sp_names)]
lizardA2 <- sp_names[grepl('lizard',sp_names) & grepl('_A2',sp_names)]
snakeA2 <- sp_names[grepl('snake',sp_names) & grepl('_A2',sp_names)]
turtleA2 <- sp_names[grepl('turtle',sp_names) & grepl('_A2',sp_names)]
mammalA2 <- sp_names[grepl('mammal',sp_names) & grepl('_A2',sp_names)]
birdA2 <- sp_names[grepl('bird',sp_names) & grepl('_A2',sp_names)]
crocodilianA2 <- sp_names[grepl('crocodilian',sp_names) & grepl('_A2',sp_names)]

# A3
frogA3 <- sp_names[grepl('frog',sp_names) & grepl('_A3',sp_names)]
caecilianA3 <- sp_names[grepl('caecilian',sp_names) & grepl('_A3',sp_names)]
lizardA3 <- sp_names[grepl('lizard',sp_names) & grepl('_A3',sp_names)]
snakeA3 <- sp_names[grepl('snake',sp_names) & grepl('_A3',sp_names)]
turtleA3 <- sp_names[grepl('turtle',sp_names) & grepl('_A3',sp_names)]
mammalA3 <- sp_names[grepl('mammal',sp_names) & grepl('_A3',sp_names)]
birdA3 <- sp_names[grepl('bird',sp_names) & grepl('_A3',sp_names)]
crocodilianA3 <- sp_names[grepl('crocodilian',sp_names) & grepl('_A3',sp_names)]

groups <- c(frogA1,caecilianA1,lizardA1,snakeA1,turtleA1,mammalA1,birdA1,frogA2,caecilianA2,lizardA2,snakeA2,turtleA2,mammalA2,birdA2,frogA3,caecilianA3,lizardA3,snakeA3,turtleA3,mammalA3,birdA3,
            crocodilianA1,crocodilianA2,crocodilianA3)

# make sure sum == length(sp_names)
sum <- 0
for(i in 1:length(groups)){
  sum <- sum + length(i)
}

sp_names[!(sp_names %in% groups)] # detect which species were not captured (should only be the fish)

# make backbone trees --> constrain main groups to be monophyletic --> 'backbone_tree.tre' file.
# *Used this to create a newick tree file manually.
cat(frogA3,"\n",sep=",")
cat(caecilianA3,"\n",sep=",")
cat(lizardA3,"\n",sep=",")
cat(snakeA3,"\n",sep=",")
cat(turtleA3,"\n",sep=",")
cat(mammalA3,"\n",sep=",")
cat(birdA3,"\n",sep=",")
cat(crocodilianA3,"\n",sep=",")


## Check that names in alignment and backbone tree match
t <- read.tree("./backbone_tree.tre")
tip_names <- t$tip.label
seq_names <- seqs$label

tip_names[!(tip_names %in% seq_names)]
seq_names[!(seq_names %in% tip_names)]
