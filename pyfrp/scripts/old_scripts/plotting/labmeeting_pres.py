#Script to sum up molecules 

#Molecule folder
mol_folder=""

#Molecules
mols=["2015_05_24_BMP2_D2_10min","2015_05_25_BMP2_D2_10min","2015_05_26_BMP2_D2_10min","2015_05_30_BMP2_D2_10min","2015_05_31_BMP2_D2_10min","2015_06_01_BMP2_D2_10min"]

#Embryos
embs=[[2,3],[4],[3,4],[2,3],[2,3,4,5],[4,6]]

for i,mol in enumerate(mols):
	
	#Load molecule file
	mol=molecule("dummy")
	mol.load_molecule(mol_folder+"/"+mol+".pk")
	
	for j in embs[i]:
		
		emb=mol.embryos[j]
		
		emb.fits=[]
		
		emb.add_fit("ext_F_cNM_below",0,"default")
		emb.add_fit("ext_F_cNM_below",1,"default")
		emb.add_fit("ext_F_cNM_below",2,"default")
		
		for fit in emb.fits:
			


