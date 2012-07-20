import re


def read_jnet(file):
	lines = [l for l in file][1:-1]
	letters = ''.join([l[2] for l in lines])
	return letters

def read_hmmtop(file):
	for line in HMMTOPFILE:
		if line.startswith('Transmembrane helices'):
			TMstr = map(lambda x: tuple(map(int,x.split('-')))+("H",), line.split(':')[1].strip().split(' '))
	return [(i[0],i[1],'H') for i in TMstr]

def make_ranges(SS_string,index=0):
	"""
	The Secondary Structure string must be: C or - for coil/loop, E for Beta-Sheet, H for helix.
	Don't ask me why.

	returns a list of tuples (Start_int, End_int, "C" or "E" or "H")
	"""

	SS_string = SS_string.replace('-','C')

	matches = re.findall('H+|E+|C+',SS_string)

	regions = list()

	i=0
	for amatch in matches:
		regions.append((i+index,i+len(amatch)-1+index,amatch[0]))
		i += len(amatch) + index

	return regions

if __name__ == '__main__':
	foo = '-----EEEEE-----EEEE----EEEEHHHHHHHHHHH------EEEEHHHHHHHHHHHHHHHHHHHH----EEEEEEHHHHHHH------HHHHHHHHHHHHHHHHHHHHHHHHHHHEEEE---EEEEE----E----HHHH---------HHHHHHHHHHHHHHH-------HHHHHHHHHHH-------HHHHHHHH----E-------HHHHHHHHHHHHHHHH--------HHHH-----HHHHHHHHHHHHHHHHH----'

	ranges = make_ranges(foo)
	print foo
	print ranges