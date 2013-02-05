"""
These functions and constants are for creating distance and dihedral constraints based on secondary structure for an EVFold reimplementaion.

The original paper doesn't spell everything out, even in the supplementary material. Specifically, the dihedral angles had to be reverse engineered from the output files provided by the EVFold authors.
Actually, there are places where the paper says there is something in supplementary material, and it isn't even there.

@author: blunt@cs.ucsd.edu



www.evfold.org

Marks, Debora S, Thomas a Hopf, and Chris Sander. "Protein Structure Prediction from Sequence Variation." Nature Biotechnology 30, no. 11 (November 8, 2012): 1072-1080. http://www.nature.com/doifinder/10.1038/nbt.2419.

Hopf, Thomas A., Lucy J. Colwell, Robert Sheridan, Burkhard Rost, Chris Sander, and Debora S. Marks. "Three-Dimensional Structures of Membrane Proteins from Genomic Sequencing." Cell (May 2012): 
1-15. http://linkinghub.elsevier.com/retrieve/pii/S0092867412005090.


"""

SSDistanceStr = 'assign (resid {0} and name {2:s})  (resid {1} and name {3:s})  {4:.1f} {5:.1f} {6:.1f} weight {7:.5f}\n'

helixConstraints = [
	[('O','O', 3.07, 0.2, 0.2),	('N','N',2.82, 0.2, 0.2),				('CA','CA',3.82,0.2,0.2),	('CA','O', 4.45, 0.2, 0.2),	('CB','CB',3.6,0.4,0.4)],
	[('O','O',4.65,0.4,0.4),	('N','N',4.4,0.4,0.4),		('O','N',3.3,0.3,0.3),	('CA','CA',5.5,0.3,0.3),	('CA','O',5.6,0.4,0.4),		('CB','CB',5.15,0.4,0.4)],
	[('O','O',5.05,0.65,0.65),	('N','N',5.0,0.6,0.6),		('O','N',3.95,0.5,0.5),	('CA','CA',5.3,0.65,0.65),	('CA','O',5.9,0.8,0.8),		('CB','CB',5.2,0.8,0.8)],
	[('O','O',6.3,0.65,0.65),	('N','N',6.25,0.8,0.8),		('O','N',4.3,0.7,0.7),	('CA','CA',6.35,0.7,0.7),	('CA','O',7.5,0.7,0.7),		('CB','CB',6.35,0.8,0.8)],
	[('O','O',8.3,0.55,0.55),	('N','N',8.2,0.5,0.5),		('O','N',6.1,0.6,0.6),	('CA','CA',8.7,0.6,0.6),	('CA','O',9.55,0.6,0.6),	('CB','CB',8.55,0.65,0.65)],
	[('O','O',9.7,0.6,0.6),		('N','N',9.6,0.55,0.55),	('O','N',7.95,0.6,0.6),	('CA','CA',10.05,0.6,0.6),	('CA','O',10.65,0.65,0.65),	('CB','CB',9.9,0.75,0.75)],
	[('O','O',10.75,0.75,0.75),	('N','N',10.75,0.65,0.65),	('O','N',9.05,0.7,0.7),	('CA','CA',10.8,0.75,0.75),	('CA','O',11.7,0.75,0.75),	('CB','CB',10.8,1.0,1.0)],
	[('O','O',12.3,0.8,0.8),	('N','N',12.3,0.7,0.7),		('O','N',10.3,0.75,0.75),('CA','CA',12.45,0.8,0.8),	('CA','O',13.5,0.75,0.75),	('CB','CB',12.45,1.0,1.0)]
]

strandConstraints = [
	[('O','O',3.4,0.3,0.3),		('N','N',3.4,0.3,0.3),									('CA','O',4.6,0.2,0.2),		('CB','CB',4.4,0.5,0.5)],
	[('O','O',6.45,0.6,0.6),	('N','N',6.45,0.6,0.6),		('O','N',4.2,0.5,0.5),	('CA','CA',6.6,0.5,0.5),	('CA','O',6.6,0.5,0.5),		('CB','CB',7.6,0.7,0.7)],
	[('O','O',9.5,1.2,1.2),		('N','N',9.5,1.2,1.2),		('O','N',7.3,0.8,0.8),	('CA','CA',9.7,1.2,1.2),	('CA','O',10.6,1.3,1.3),	('CB','CB',9.9,1.3,1.3)],
	[('O','O',12.5,1.6,1.6),	('N','N',12.5,1.6,1.6),		('O','N',10.3,1.3,1.3),	('CA','CA',12.6,1.6,1.6),	('CA','O',13.6,1.7,1.7),	('CB','CB',12.7,1.7,1.7)],
	[('O','O',15.4,2.2,2.2),	('N','N',15.4,2.2,2.2),		('O','N',13.3,1.9,1.9),	('CA','CA',15.5,2.3,2.3),	('CA','O',16.4,2.3,2.3),	('CB','CB',15.6,2.3,2.3)],
	[('O','O',18.0,3.0,3.0),	('N','N',18.0,3.0,3.0),		('O','N',16.0,2.6,2.6),	('CA','CA',18.1,3.0,3.0),	('CA','O',19.0,3.2,3.2),	('CB','CB',18.2,3.1,3.1)],
	[('O','O',20.5,4.0,4.0),	('N','N',20.5,4.0,4.0),		('O','N',18.6,3.6,3.6),	('CA','CA',20.5,4.1,4.1),	('CA','O',21.3,4.2,4.2),	('CB','CB',20.6,4.1,4.1)]
]

loopConstraints = []

secondaryDict = {'H':helixConstraints,'E':strandConstraints,'C':loopConstraints}

def make_SS_dists(start,length,type='H',weight=5.0):
	end = start+length

	SS_constraints = secondaryDict[type]

	retstr = ""

	for i in range(start,end):
		for j,offset in zip(range(i+1,end),range(len(SS_constraints))):
			for single in SS_constraints[offset]:
				retstr += SSDistanceStr.format(*((i,j) + single + (weight,)))
	
	return retstr


#TODO : Figure out about dihedral angles in CNS, is the last argument a float or an int? shouldnt it be the weight or something? (It seems to be an int in the paper...)

#dihedralStrs = ["assign (resid {0} and name c) (resid {1} and name n) (resid {1} and name ca) (resid {1} and name c)  {2:.1f} {3:.1f} {4:.1f} {5:.1f}\n",
#		"assign (resid {0} and name n) (resid {0} and name ca) (resid {0} and name c) (resid {1} and name n)  {2:.1f} {3:.1f} {4:.1f} {5:.1f}\n"]

dihedralStrs = ["assign (resid {0} and name c) (resid {1} and name n) (resid {1} and name ca) (resid {1} and name c)  {2:.1f} {3:.1f} {4:.1f} {5:.0f}\n",
		"assign (resid {0} and name n) (resid {0} and name ca) (resid {0} and name c) (resid {1} and name n)  {2:.1f} {3:.1f} {4:.1f} {5:.0f}\n"]

alphaDihedrals = [(5.0,-57.0,7.0,2),(5.0, -47.0, 7.0, 2)]
betaDihedrals  = [(5.0,-120.0,7.0,2),(5.0,120.0,7.0,2)]
loopDihedrals = []

secondaryAngleDict = {'H':alphaDihedrals, 'E':betaDihedrals, 'C':loopDihedrals } 

def make_SS_angles(start,length,type="H",weight=2.0):
	retstr = ""

	for i in range(start,start+length):
		for oneDihedralStr,theAngles in zip(dihedralStrs,secondaryAngleDict[type]):
			retstr += oneDihedralStr.format(*((i,i+1) + theAngles))
	
	return retstr

if __name__ == "__main__":
	pass
	#TODO : Put some unit tests here.
