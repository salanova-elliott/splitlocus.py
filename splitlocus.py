#!/usr/local/Cluster-Apps/python/3.6.8/bin/python3
#Script is modified from 00_SortbyLocus_FastA_ALL.py

from fuzzywuzzy import fuzz
import sys

#Forward primer dictionary
#Actual ZBJ primer used was AGATATTGGAACWTTATATTTTATTTTTGG; using T for the W wobble below b/c it was more abundant in test sample
primerdict = {
	"zbj" : "AGATATTGGAACTTTATATTTTATTTTTGG",
	"anml": "GGTCAACAAATCATAAAGATATTGG",
	"plants3" : "CTAAATTGGGATTATCCGCT",
	"plants5" : "TTTCACTCAAGATTGGGTTTCT",
	"plants7" : "CTCCTGAATATGAAACCAAAGA",
	"plantsc" : "CGAAATCGGTAGACGCTACG",
	"plantsITS1" : "AGAAGTCGTAACAAGGTTTCCGTAGG"
}

#Output dictionary with each primer (and unknown category) associated with a list
outdict = {primer : [] for primer in primerdict}
outdict["unknown"] = []

#Import sequences
sequencelist = []
with open(sys.argv[1]) as f:
	for line in f:

		if line.startswith(">"):
			header = line.rstrip()
		elif line.startswith(("A","C","G","T")):
			sequencelist.append((header, line.rstrip()))

#All console output is here
def lengthDisplay():
	print(sys.argv[1], file=sys.stderr)
	print("Number of sequences is: {}".format(len(sequencelist)), file=sys.stderr)
	for primer in outdict:
		print("{} {}".format(primer, len(outdict[primer])), file=sys.stderr)

#File writer
def writeFile(primer):
	with open("{}_{}.fasta".format(sys.argv[2], primer), "w") as outf:
		for sequence in outdict[primer]:
			outf.write("{}\n{}\n".format(sequence[0],sequence[1]))

#Finds the best primer match for each sequence
for sequence in sequencelist:
	currentmax = ("", 0)

	#Scores each primer and keeps track of the current max
	for primer in primerdict:
		score = fuzz.partial_ratio(sequence[1][0:39], primerdict[primer])
		if score >= 90 and score > currentmax[1]:
			currentmax = (primer, score)

	#Adds sequence to appropriate list
	if currentmax[1] >= 90:
		outdict[currentmax[0]].append(sequence)
	else:
		outdict["unknown"].append(sequence)

#Writes out files and console data
#lengthDisplay()
for primer in outdict:
	writeFile(primer)


