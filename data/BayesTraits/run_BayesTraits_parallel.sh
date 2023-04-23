#!/bin/bash

# Running BayesTraits for N files:
# 
# Global Variables
TREE="/Users/santiagoherrera/Desktop/PROYECTOS/ATPase_Project/ASR/a1-3_v2/Epistasis_analysis/BayesTraits_data/codeml_brlens_tree_nexus.tre"

NULL="/Users/santiagoherrera/Desktop/PROYECTOS/ATPase_Project/ASR/a1-3_v2/Epistasis_analysis/BayesTraits_data/null_discrete_independent_v2.txt"

ALT="/Users/santiagoherrera/Desktop/PROYECTOS/ATPase_Project/ASR/a1-3_v2/Epistasis_analysis/BayesTraits_data/alt_discrete_dependent_v2.txt"


#Function to run BayesTraits on a single file
# $1 is $TREE
# $2 is input file
# $3 is $NULL
# $4 is $ALT

runBT () {
	# Likelihood of null model per site
	LIK_NULL=$(../BayesTraitsV3 "$1" "$2" < "$3" | tail -n 2 | head -n 1 | cut -f 2)

	# Likelihood of alt model per site
	LIK_ALT=$(../BayesTraitsV3 "$1" "$2" < "$4" | tail -n 2 | head -n 1 | cut -f 2)

	#Write results to file
	echo "$2":$LIK_NULL:$LIK_ALT | column -t -s ":" >> out_BT.txt
}
export -f runBT

#runBT $TREE V134_V1.txt $NULL $ALT # This works

#Run the function in parallel using 10 jobs
parallel -j10 runBT $TREE {} $NULL $ALT ::: V1*.txt