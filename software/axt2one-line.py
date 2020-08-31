#!/usr/bin/python
import sys
import re
mydict = {}
args = sys.argv
reader = args[1]
writer = open(args[2],"w")
with open(reader) as fh:
	for line in fh:
		if re.search(r"\d+",line):
			header = line.strip()
			mydict[header] = []
		else:
			mydict[header].append(line.strip())
for key,value in mydict.items():
	seqs = "".join(value)
	length = len(seqs)
	mid = int(length/2)
	seq1 = seqs[0:mid]
	seq2 = seqs[mid:]
	writer.write(key + "\n" + "\n".join([seq1,seq2]))

