"""A library for dealing with the similarity file that we have created, including classes for representing their headers, and a subclass of Bio.Align.Generic.Alignment that maps sequence IDs to SeqRecord s."""
import os.path

import Bio.SeqIO
import Bio.AlignIO
import Bio.Align.Generic
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import *
from Bio.SeqRecord import *
import re
import copy

from amaranth.utils import flatten
from itertools import izip
import scipy as S

from __Header import *

class AmaranthAlignment(Bio.Align.MultipleSeqAlignment):
	"""This object represents an alignment who's members are grouped sequences, such as protein pairings.
	In addition to attributes inherited from Bio.Align.Generic.Alignment, it also contains a dictionary mapping sequence names to the SeqRecord objects that contain them, and an amaranth.simformat.Header .
	"""
	def __init__(self,alphabet=Alphabet.Gapped(IUPAC.protein),header=None,dictKey=None):
		super(AmaranthAlignment,self).__init__(None,alphabet)
		self.namedict = dict()
		"""A dictionary mapping the sequence name, or some portion of it (such as GI or accession) to the SeqRecord object."""
		if header == None:
			self.header = Header()
		elif isinstance(header,str):
			self.header = Header(header)
		else:
			self.header = copy.deepcopy(header)
		"""A header representing information about the names of domains that are grouped in this alignment, their widths, and similarity cutoffs."""
		self.dictKey = dictKey
		"""Store the name of the key used to index SeqRecords in the namedict."""

	def addSeqReq(self,seqReq):
		self.append(seqReq)
		
	def append(self,seqReq):
		super(AmaranthAlignment,self).append(seqReq)
		if self.dictKey:
			seqKey = seqReq.annotations.get(self.dictKey,None)
			if seqKey:
				self.namedict[seqKey] = seqReq
		
	def get_alignment_length(self):
		if len(self._records) > 0:
			return len(self._records[0])
		else:
			return 0

	def get_weights(self,cutoff=0.7):
		"""
		Gets the inverse of the normalized similarities on an alignment in our format.
		That is, if there are N other sequences above the similarity threshold for the current entry, it's weight will be 1/(N+1), an entry with no similar entries is weighted 1, 1 similar entry gets a weight of 1/2, etc.
		"""
		
		index = self.header.cutoffs.index(cutoff)
		return [(1.0/(X.annotations["weights"][index]+1.0)) for X in self._records]

def read(infile,format="fasta",converter=lambda x:x):
	"""
	Read an alignment in our format, with header and potentially similarity data, this is not guaranteed to work on anything but fasta, in fact, I'd be surprised if it didn't blow up on everything else.
	"""
	aheader = readHeader(infile)
	
	alignment = AmaranthAlignment(header = aheader)
	for record in Bio.SeqIO.parse(infile,format):
		annotateSeqRecord(record)
		alignment.append(record)
	
	if aheader == None: #if the file has no header, we create one, giving the filename and width as one domain.
		fileBasenameNoExtention = os.path.basename(infile.name).rsplit(".",1)[0]
		alignment.header.addDomain( (fileBasenameNoExtention, alignment.get_alignment_length()) )
	else:
		assert alignment.get_alignment_length() == sum(alignment.header.domainLengths), "Alignment not of the width expected from the header."
	
	return alignment

def __writerTitleGen(seqReq):
	title = seqReq.name + " "
	if seqReq.annotations.has_key("Eval"):
		title += " E={%s}" % ",".join(seqReq.annotations["Eval"])
	if seqReq.annotations.has_key("bits"):
		title += " B={%s}" % ",".join(seqReq.annotations["bits"])
	if seqReq.annotations.has_key("organism"):
		title += " [%s]" % seqReq.annotations["organism"]
	elif seqReq.annotations.has_key("source"):
		title += " [%s]" % seqReq.annotations["source"]

	if seqReq.annotations.has_key("weights"):
		title += " <%s>" % ",".join([str(k) for k in seqReq.annotations["weights"]])
	return title

def write(outfile,align,header=None,format="fasta"):
	"""
	Write a file in our format, if header is provided, it will override align.header, if it is not provided and align.header exists, it will be written as header to the file.
	Fasta data is written after that.
	"""
	if header is not None:
		outfile.write(str(header))
	elif hasattr(align,"header"):
		outfile.write(str(align.header))
	
	if format is "fasta":
		fWriter = Bio.SeqIO.FastaIO.FastaWriter(outfile,wrap=0,record2title=__writerTitleGen)
		fWriter.write_file(align)
		outfile.close()
	else:
		Bio.AlignIO.write([align],outfile,format)
	
REJECT = 0
EVAL_GREATEST = 1
GROUP = 2
KEEP_FIRST = 3
def alignmentToNamedict(analign,key="gi",collision=REJECT):
	"""
	DEPRECATED, Just Use An AmaranthAlignment
	Run over an Alignment and create entries in its namedict attribute.
	The key variable defines the annotation from each SeqRecord's annotations dictionary that will be used as that SeqRecord's name.
	Be WARNED that this function implicitly calls annotateSeqRecord on every SeqRecord in analign, which adds entries to their annotations dictionary.
	
	There are currently four actions that can be taken on sequences that would end up having colliding names.
	By providing the appropriate flag, they can be Rejected (REJECT (Default)), kept according to the one with the greatest EVal (EVAL_GREATEST), grouped so as not to be lost (GROUP (**Unimplemented**)), or just the first that was encountered can be kept (KEEP_FIRST (**Unimplemented**)).
	
	"""
	import warnings
	warnings.warn("simformat.alignmentToNamedict is DEPRECATED!")
	namedict = analign.namedict = dict()#assign a namedict to analign if it does not already have one. DOES this WORK?
	
	#This maintains backwards compatibility with string-typed keys by converting them to a key lambda
	if isinstance(key,str):
		originalKey = key#Remember that there are no closures in python
		key = lambda x:x.annotations[originalKey]
	
	for i in analign:
		annotateSeqRecord(i)
		dictKey = key(i)
		if not namedict.has_key(dictKey):
			namedict[dictKey] = i
		else:
			if collision == REJECT:
				namedict.pop(dictKey)
			elif collision == EVAL_GREATEST:
				if float(namedict[dictKey].annotations["Eval"]) < float(i.annotations["Eval"]):
					namedict[dictKey] = i

def annotateAlignment(analign):
	"""
	Run over all members of an Alignment and annotate them.
	"""
	for i in analign:
		annotateSeqRecord(i)

giMatcher=re.compile("gi\|(\d*)")
specMatcher=re.compile("\[(.*)\]")
floatRE = "([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
floatMatcher = re.compile(floatRE)
eValMatcher=re.compile("E=\S+")
bitMatcher=re.compile("B=\S+")
subSeqMatcher=re.compile("/(\d+)(?:~|-)(\d+)")
weightsMatcher=re.compile("<(\d+(,\d+)*)>")
def annotateSeqRecord(seqreq):
	"""
	Scan a SeqRecord taken from a fasta alignment generated by the MCHammer tool, and add annotations "gi", "organism", "source", "Eval", and "begin" and "end" if available.
	If this is a grouped alignment, gi will be a tuple of GIs, in all cases, GIs are held as strings.
	"""
	giMatch = giMatcher.findall(seqreq.name)
	if len(giMatch) > 0:
		seqreq.annotations["gi"] = tuple(giMatch)

	subSeqMatch = subSeqMatcher.search(seqreq.name)
	if subSeqMatch:
		seqreq.annotations["start"] = subSeqMatch.groups()[0]
		seqreq.annotations["end"] = subSeqMatch.groups()[1]
		seqreq.annotations["region"] = tuple(subSeqMatch.groups())
	specMatch = specMatcher.search(seqreq.description)
	if specMatch:
		seqreq.annotations["organism"] = specMatch.groups()[0]
		seqreq.annotations["source"] = specMatch.groups()[0]
	eValMatch = eValMatcher.search(seqreq.description)
	if eValMatch:
		eVals = floatMatcher.findall(eValMatch.group())
		seqreq.annotations["Eval"] = tuple(eVals)
	bitsMatch = bitMatcher.search(seqreq.description)
	if bitsMatch:
		bits = floatMatcher.findall(bitsMatch.group())
		seqreq.annotations["bits"] = tuple(bits)
	weightsMatch = weightsMatcher.search(seqreq.description)
	if weightsMatch:
		weightList = floatMatcher.findall(weightsMatch.group())
		seqreq.annotations["weights"] = [int(i) for i in weightList]


def groupReqs(seqRecs):
	"""
	Join SeqRecord objects according to our pairing format and return the new SeqReq object that gets created.
	SeqRecs will be joined in the order that they were provided in the argument seqRecs, which should be a list.
	"""
	#seqReqs = flatten(seqReqs)
	#TODO fix this so that seqReqs can either be a tuple with only a list of seqRecs, or can be a tuple of seqRecs
	idString = ":".join([i.name.split(" ")[0] for i in seqRecs])
	mySeq = seqRecs[0]
	for i in seqRecs[1:]:
		mySeq+=i
	#handle merging of annotations
	mySeq.id = idString
	mySeq.name = idString
	mySeq.annotations["gi"] = tuple(flatten([k.annotations.get("gi","") for k in seqRecs]))
	mySeq.annotations["organism"] = seqRecs[0].annotations.get("organism","")
	mySeq.annotations["source"] = seqRecs[0].annotations.get("source","")
	mySeq.annotations["Eval"] = tuple(flatten([i.annotations.get("Eval","") for i in seqRecs]))
	mySeq.annotations["bits"] = tuple(flatten([i.annotations.get("bits","") for i in seqRecs]))

	return mySeq


def listWindow(list,n=2):
	for i in range(len(list)-n+1):
		yield tuple(list[i:i+n])

def splitReqs(seqRec,lengths):
	names = seqRec.name.split(":")
	splitpoints=S.cumsum([0]+lengths)
	retval = list()
	for aname,arange in izip(names,listWindow(splitpoints)):
		arec = SeqRecord(seqRec.seq[arange[0]:arange[1]],name=aname,id=aname)
		retval.append(arec)
		
	orgname = seqRec.annotations.get("organism",
									seqRec.annotations.get("source",None))
	if orgname:
		for i in retval:
			i.annotations["organism"] = orgname
			i.annotations["source"] = orgname
			
	for splitKey in ["gi","Eval","bits"]:
		evals = seqRec.annotations.get(splitKey,None)
		if evals:
			for i,j in izip(retval,evals):
				i.annotations[splitKey] = (j,)

	
	return retval
	
	
	

def splitAligns(align,head=None):
	"""
	Read a paired alignment in our format, and split it into separate alignment objects, one for each of the original domains used to form the pair/triplet/etc file.
	"""
	if head is None and hasattr(align,"header"):
		head = align.header
	splitpoints = head.domainLengths
	splitList = map(lambda x: splitReqs(x,splitpoints),align._records)
	outAligns = [AmaranthAlignment() for i in head.domainLengths]
	for h,n,l in izip(outAligns,head.domainNames,head.domainLengths):
		h.header.domainNames.append(n)
		h.header.domainLengths.append(l)
	
	for i in splitList:
		for j,k in zip(i,outAligns):
			k._records.append(j)
	return outAligns


def getinvnormsim(theAlign,index=3):
	"""
	DEPRECATED, use Align.get_weights(CUTOFF_SIMILARITY)
	"""
	return [(1.0/(X.annotations["weights"][index]+1.0)) for X in theAlign]


