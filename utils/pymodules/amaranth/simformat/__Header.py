'''
Created on Sep 21, 2010

@author: lunt
'''

import re

#TODO refactor this and all client-code so that it supports any varying number of domains, triplets/quadruplets, etc
class Header:
	"""This object represents the information contained in the header of a simformat file that has been read in, or a grouped-alignment that has been generated in memory.
	it provides an easy to use API for finding the similarity comparison cut-offs that were used to generate an alignment, and the names and sizes of domains appearing in a grouped-alignment.
	"""
	#
	def __init__(self,header=""):
		self.domainNames = list()
		self.domainLengths = list()
		self.cutoffs = list()
		self.origString = header

	def __str__(self):
		retstring = ""
		for i,j in zip(self.domainNames,self.domainLengths):
			retstring += "# %s : %i\n" % (i,j)
		if self.cutoffs is not None and len(self.cutoffs) > 0:
			retstring += "# " + ",".join([format(i) for i in self.cutoffs]) + "\n"
		retstring += "######BEGIN FASTA######\n"
		return retstring
	
	def __add__(self,other):#
		assert isinstance(other,Header) , "only two headers can be added together, type of other is : %r , class of other is %r" % (type(other),other.__class__)
		retval = Header("")
		retval.domainNames = list(self.domainNames)
		retval.domainNames.extend(other.domainNames)
		retval.domainLengths = list(self.domainLengths)
		retval.domainLengths.extend(other.domainLengths)
		return retval
	
	@staticmethod
	def join(severalHeads):
		#assert all([isinstance(i,Header) for i in severalHeads]) , "only headers can be added together"
		retval = Header("")
		for i in severalHeads:
			retval += i
		return retval

	def addDomain(self,nameLenTuple):
		"""Add a domain named nameLenTuple[0] of length nameLenTuple[1] to this header.
		nameLenTuple must be a tuple of the form: (domain_name, domain_length)
		"""
		self.domainNames.append(nameLenTuple[0])
		self.domainLengths.append(nameLenTuple[1])


def readHeader(infile):
	"""This function reads a header from our file-format and instantiates it.
	This function leaves the file handle with its iterator/cursor set to the beginning of the line just after the header.
	"""
	header = ""
	while True:#this loop is needed because there is no peek function.
		beforepos = infile.tell()
		line = infile.readline()
		if line.startswith('#'):
			header += line
		else:
			infile.seek(beforepos)
			break
	if header == "":
		return None
	else:
		return parseHeader(header)
	

def parseHeader(header):
	retHeader = Header(header)
	
	if header is not "" and header is not None:
			lines = header.splitlines()
			assert lines[-1].find("BEGIN FASTA") > -1 , "Is this Header malformed?"
			end = -1#unless we find a line with similarity cutoffs in the next bit, we want to read all but the last line of the header.
			#this bit seems like it could be very error-prone
			cutoffLineRe = "# *(?:(\d+\.\d+) ?, *)+(\d+\.\d+) ?,? ?"
			if re.match(cutoffLineRe,lines[-2]):
				end = -2#and ignore the last TWO lines of the header for domain names
				cutoffre = "(\d+\.\d+)"
				cutoffmatch = re.findall(cutoffre,lines[-2])#check the second to last line of the header to see if it has similarity cutoffs
				if(len(cutoffmatch) > 0):#if there are cutoffs
					retHeader.cutoffs = [float(i) for i in cutoffmatch]#set the cutoffs
			domainRe = "# *(\S*) *: *(\d+)"
			for oneline in lines[:end]:#for every header line that should have a domain-name and length
				aMatch = re.match(domainRe,oneline)
				if aMatch:
					name,length = aMatch.groups()
					length = int(length)
					retHeader.domainNames.append(name)
					retHeader.domainLengths.append(length)
				else:
					assert False, "Header Lines Malformed!"
	
	
	
	return retHeader