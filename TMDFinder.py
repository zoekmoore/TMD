'''
Zoe Moore 9/15/2019

A program that reads RNA sequences from a text file and translates
them to their corresponding amino acid sequences.
'''
from readfasta import readfasta
from RNATranslate import translate
from TMDFilter import assignHydrophobicity
from TMDFilter import filterHydrophobicity

# Read in the RNA sequences from a file specified by user input
filename = input("Please enter the input file name: ")
rnainfo = readfasta(filename)

# Prepare to re-write the RNA sequences to an output file specified by user input
outfilename = input("Please enter the rna output file name: ")
handle = open(outfilename, mode="w")

# Separate out three RNA sequences and write them to separate lines of a .txt file
seqone = rnainfo[0][2]
handle.write(seqone + "\n\n")

seqtwo = rnainfo[1][2]
handle.write(seqtwo + "\n\n")

seqthree = rnainfo[2][2]
handle.write(seqthree + "\n\n")

handle.close()

# Translate RNA sequences to their single-letter amino acid sequences
aaseqone=translate(seqone)
aaseqtwo=translate(seqtwo)
aaseqthree=translate(seqthree)

# Save each amino acid sequence to an output file

# Prepare to write the amino acid sequences to an output file specified by user input
outfilename = input("Please enter the amino acid output file name: ")
handle = open(outfilename, mode="w")

# Save each amino acid sequences by writing them to separate lines of a .txt file
handle.write(aaseqone + "\n\n")

handle.write(aaseqtwo + "\n\n")

handle.write(aaseqthree + "\n\n")

handle.close()

# For each amino acid sequence, determine the regions under our window size of 11
# that have a hydrophobicity sum under our threshold of 5
hydroone = assignHydrophobicity(aaseqone)
hydrotwo = assignHydrophobicity(aaseqtwo)
hydrothree = assignHydrophobicity(aaseqthree)

print(hydroone)

# Prompt user for where to save transmembrane domain results
tmdfilename = input("Please enter the a hydrophobicity score output file name: ")

tmdone = filterHydrophobicity(hydroone)

print(tmdone)