offset =33 # 33 for most newer reads, 64 for old Illumina reads

class Sequence:
	def __init__(self,name="",seq=""):
		self.name = name
		self.seq = seq
		self.qualstr = "" #string of quality scores in characters
		self.qualarr = [] #list of quality scores in ASCII numbers

	#input a string of quality characters, set both qualstr and qualarr
	def set_qualstr(self,qual):
		self.qualstr = qual
		if len(self.qualarr) == 0:
			for j in self.qualstr:
				quality_score = ord(j) - offset
				self.qualarr.append(quality_score)
				assert quality_score <= 41 and quality_score > 0, \
					"Change the offset in seq.py\nqual"

	#output in four-line fastq format
	def get_fastq(self):
		retstr = "@"+self.name+"\n"
		retstr += self.seq+"\n"
		retstr += "+\n"
		retstr += self.qualstr+"\n"
		return retstr

	def get_fasta(self):
		return ">"+self.name+"\n"+self.seq+"\n"

	def rev_comp(self):
		tseq = ""
		for i in self.seq[::-1]:
			if i.lower() == "a":
				tseq += "t"
			elif i.lower() == "t":
				tseq += "a"
			elif i.lower() == "c":
				tseq += "g"
			elif i.lower() == "g":
				tseq += "c"
			else:
				tseq += i
		self.seq = tseq

def fastq_generator(infile):
	line = infile.readline()
	while len(line) > 0: #no empty string at end-of-file in readline
		if line[0] == "@":
			name = line[1:].strip() #name of sequence, minus "@"
			seq = infile.readline().strip() #the actual sequence
			line = infile.readline().strip() #"+"
			qual = infile.readline().strip() #the quality string
			tseq = Sequence(name=name,seq=seq)
			tseq.set_qualstr(qual)
			yield tseq
		line = infile.readline()



def read_fasta_file(filename):
	"""given the path to a fasta file, retern a list of seq objects"""
	fl = open(filename,"r")
	seqlist = [] #list of sequence objects
	templab = ""
	tempseq = ""
	first = True
	for i in fl:
		if i[0] == ">":
			if first == True:
				first = False
			else:#if not first store lastseq read
				seqlist.append(Sequence(templab,tempseq))
			templab = i.strip()[1:].split(" ")[0]
			tempseq = ""
		else:
			tempseq = tempseq + i.strip()
	fl.close()
	seqlist.append(Sequence(templab,tempseq))
	return seqlist
