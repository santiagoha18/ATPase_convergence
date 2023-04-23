from Bio import PDB
import numpy
import os, sys

# RUN: python pairwise_distance.py $PWD 3b8e.pdb

# This script computes the 3D distance between the alpha carbons of all the pairs of sites with a given
# focal site.

try:
    sourcedir=sys.argv[1]
    print(sourcedir + '...OK!')
except:
    print('Please pass directory_name')


try:
    PDBfile=sys.argv[2]
    print(PDBfile + '...OK!')
except:
    print('Please pass ATP1A1 PDB file')    


#Import ATP1A1 structure
parser = PDB.PDBParser()
pdb_path = os.path.join(sourcedir,PDBfile)
pdb_id = "3b8e"
struct = parser.get_structure(pdb_id, pdb_path)
atp1a1 = struct[0]["A"]

#Select alpha carbon of residue 111
res111 = atp1a1[111]["CA"]

#Select alpha carbon of residue 122
res122 = atp1a1[122]["CA"]


def pairwise_dist(protein, site, focal_residue):
	site1 = []
	site2 = []
	s2_resname = []
	dist = []
	for residue in protein:
		if focal_residue != residue:
			site1.append(site)
			site2.append(residue.id[1])
			s2_resname.append(residue.resname)
        	
        	# compute distance between CA atoms
			try:
				distance = focal_residue - residue["CA"]
				dist.append(distance)
			except KeyError:
            ## no CA atom, e.g. for H_NAG
				dist.append("NA")
	df=numpy.vstack((site1,site2,s2_resname,dist))
	return(df)


pairs111 = pairwise_dist(atp1a1,"111",res111)
pairs122 = pairwise_dist(atp1a1,"122",res122)


#Write output for site 111
filepath1=os.path.join(sourcedir,'distance_pairs111.out')
out_file=open(filepath1,'w')

out_file.write("site1" + "\t" + "site2" + "\t" + "site2_aa" + "\t" + "distance" + "\n")
for i in range(len(pairs111[0])):
    out_file.write(pairs111[0][i] + '\t' + str(pairs111[1][i]) + '\t' + pairs111[2][i] + "\t" + str(pairs111[3][i]) + '\n')
out_file.close()

#Write output for site 122
filepath2=os.path.join(sourcedir,'distance_pairs122.out')
out_file=open(filepath2,'w')

out_file.write("site1" + "\t" + "site2" + "\t" + "site2_aa" + "\t" + "distance" + "\n")
for i in range(len(pairs122[0])):
    out_file.write(pairs122[0][i] + '\t' + str(pairs122[1][i]) + '\t' + pairs122[2][i] + "\t" + str(pairs122[3][i]) + '\n')
out_file.close()



