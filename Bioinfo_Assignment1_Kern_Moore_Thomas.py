'''
Zoe Moore
Updated: 9/25/2019

Reads RNA sequences from an input .txt file and determines whether or
not they contain transmembrane domains
'''
from readfasta import readfasta

'''
main - streamlines the process of determining transmembrane domains by
calling necessary functions for each step of the process; handles formatting
output and files
'''
def main():
    # Read in the RNA sequences from a file specified by user input
    filename = input("Please enter the input file name: ")
    rnainfo = readfasta(filename)

    # Prepare to re-write the RNA sequences to an output file specified by user input
    outfilename = input("Please enter the output file name: ")
    handle = open(outfilename, mode="w")

    # Iterate through each RNA sequence in the input file
    for i in range(len(rnainfo)):
        # Specify gene that is being evaluated
        handle.write("Gene " + str(i+1) + ":\n" + rnainfo[i][2] + "\n\n")
        
        # Translate the RNA Sequence to its corresponding single-letter amino acid sequence
        # Write information to the output file
        translatedseq = translate(rnainfo[i][2])
        handle.write("Protein Sequence " + str(i+1) + ": " + translatedseq + "\n\n")
        
        # Scan the single-letter amino acid sequence for transmembrane domains
        # Write results to the output file
        findTMD(translatedseq, handle)

    # Close file
    handle.close()

'''
translate - a subroutine that translates an RNA sequence into its
corresponding single-letter amino acid sequence.

Parameter: an RNA sequence that includes only one line of nucleotides
formatted through readfasta.py. Must include only the characters: 'A', 'C',
'G', and 'U'.
'''
def translate(sequence):
    # Create an empty amino acid sequence list to copy translations into
    aasequence = ''

    # Establish a code of trios of nucleotides and their corresponding single-letter amino acids
    code = {'UCA' : 'S',    # Serine
            'UCC' : 'S',    # Serine
            'UCG' : 'S',    # Serine
            'UCU' : 'S',    # Serine
            'UUC' : 'F',    # Phenylalanine
            'UUU' : 'F',    # Phenylalanine
            'UUA' : 'L',    # Leucine
            'UUG' : 'L',    # Leucine
            'UAC' : 'Y',    # Tyrosine
            'UAU' : 'Y',    # Tyrosine
            'UAA' : '_',    # Stop
            'UAG' : '_',    # Stop
            'UGC' : 'C',    # Cysteine
            'UGU' : 'C',    # Cysteine
            'UGA' : '_',    # Stop
            'UGG' : 'W',    # Tryptophan
            'CUA' : 'L',    # Leucine
            'CUC' : 'L',    # Leucine
            'CUG' : 'L',    # Leucine
            'CUU' : 'L',    # Leucine
            'CCA' : 'P',    # Proline
            'CCC' : 'P',    # Proline
            'CCG' : 'P',    # Proline
            'CCU' : 'P',    # Proline
            'CAC' : 'H',    # Histidine
            'CAU' : 'H',    # Histidine
            'CAA' : 'Q',    # Glutamine
            'CAG' : 'Q',    # Glutamine
            'CGA' : 'R',    # Arginine
            'CGC' : 'R',    # Arginine
            'CGG' : 'R',    # Arginine
            'CGU' : 'R',    # Arginine
            'AUA' : 'I',    # Isoleucine
            'AUC' : 'I',    # Isoleucine
            'AUU' : 'I',    # Isoleucine
            'AUG' : 'M',    # Methionine
            'ACA' : 'T',    # Threonine
            'ACC' : 'T',    # Threonine
            'ACG' : 'T',    # Threonine
            'ACU' : 'T',    # Threonine
            'AAC' : 'N',    # Asparagine
            'AAU' : 'N',    # Asparagine
            'AAA' : 'K',    # Lysine
            'AAG' : 'K',    # Lysine
            'AGC' : 'S',    # Serine
            'AGU' : 'S',    # Serine
            'AGA' : 'R',    # Arginine
            'AGG' : 'R',    # Arginine
            'GUA' : 'V',    # Valine
            'GUC' : 'V',    # Valine
            'GUG' : 'V',    # Valine
            'GUU' : 'V',    # Valine
            'GCA' : 'A',    # Alanine
            'GCC' : 'A',    # Alanine
            'GCG' : 'A',    # Alanine
            'GCU' : 'A',    # Alanine
            'GAC' : 'D',    # Aspartic Acid
            'GAU' : 'D',    # Aspartic Acid
            'GAA' : 'E',    # Glutamic Acid
            'GAG' : 'E',    # Glutamic Acid
            'GGA' : 'G',    # Glycine
            'GGC' : 'G',    # Glycine
            'GGG' : 'G',    # Glycine
            'GGU' : 'G'}    # Glycine

    # Go through each trio of nucleotides in the RNA sequence
    for i in range(0, len(sequence), 3):
        # Add the corresponding single-letter amino acid to the amino acid sequence
        aasequence += code[sequence[i:i+3]]

        # If the stop codon is reached, the translation is finished
        if code[sequence[i:i+3]] == '_':
            break

    # Return the resulting amino acid sequence
    return aasequence

'''
assignHydrophobicity - uses hydrophobicities based on the Hessa
biological hydrophobicity scale assigned to each amino acid. Sums
the values on sliding windows of size 11 amino acids, and it records
the indices of hydrophobicity sums below the threshold of value 5.

Parameter: a single-letter amino acid sequence translated from an RNA
sequence.
'''
def assignHydrophobicity(sequence):
    # Create hydrophobicity assignments per amino acids
    # Based on Hessa Scale
    code = {'C' : -0.13,
            'N' : 2.05,
            'Q' : 2.36,
            'S' : 0.84,
            'T' : 0.52,
            'D' : 3.49,
            'E' : 2.68,
            'K' : 2.71,
            'R' : 2.58,
            'A' : 0.11,
            'G' : 0.74,
            'I' : -0.60,
            'L' : -0.55,
            'M' : -0.10,
            'P' : 2.23,
            'V' : -0.31,
            'F' : -0.32,
            'H' : 2.06,
            'W' : 0.30,
            'Y' : 0.68}

    # Create a list to store the hydrophobicity scores
    hydrosequence = list()

    # Initialize a sum variable for the window size
    windowSum = 0

    # Iterate through the amino acid sequence
    # Sum the hydrophobicity scores in window sizes of 11 amino acids
    # If the sum is lower than our threshold of 5, store the start index i
    # of that sum [taken from index i to i + 11] into the hydrosequence list
    for i in range(0, len(sequence) - 11, 1):       # iterate through the amino acid sequence
        for j in range(i, i + 11, 1):               # iterate through 11 amino acids, per the window size
            # sum the hydrophobicity values of the window
            windowSum += code[sequence[j]]         
        if(windowSum < 5):                          # if the sum is less than our threshold of 5 . . .
            # store the index of that low hydrophobicity region
            hydrosequence.append(i)                 
        windowSum = 0
        
    # Return the list of indices containing hydrophobicity value sums less than our threshold
    return hydrosequence
    
'''
filterHydrophobicity - filters through the windows of hydrophobicity
sums below the threshold of 5, looking for sequences between length 18 to 30 amino
acids (or indices of 7 to 19 in terms of the length 11 sliding windows), per the
requirement for transmembrane domain length, that do not contain any breaks.
Stores the indices of the continuous sequences, if they exist.

Parameter: the sequence of hydrophobicity value sums on the window size of length 11
amino acids.
'''
def filterHydrophobicity(hydrosequence):
    # Initialize the possible start of a transmembrane domain sequence
    tempstart = -1

    # Initialize the length of the transmembrane domain sequence
    length = 0

    # Initialize a list of the transmembrane domain locations
    tmdlocations = ''

    # Iterate through the hydrophobicity sums to look for consecutive
    # sequences of length 18 to 30 amino acids
    for i in range(len(hydrosequence) - 1):
        if ((hydrosequence[i+1] - hydrosequence[i]) == 1 and i != len(hydrosequence) - 2):      # If there is not a gap in hydrophobicity sums
            if tempstart == -1:                                                                 # If the sequence start has been reset
                # Store the possible starting index
                tempstart = hydrosequence[i]                                                    
            # Increment the length of the sequence
            length += 1
        else:                                                                                   # If a gap has been found in the sequence
            if length > 6 and length < 20:                                                      # If the sequence is between 18 to 30 amino acids
                # Add the interval of the sequence start and sequence end to the transmembrane domain list
                tmdlocations += str(tempstart) + ':' + str(hydrosequence[i] + 11) + '\n'
            
            # Reset the sequence start and sequence length trackers
            tempstart = -1
            length = 0
            
    return tmdlocations

'''
findTMD - calls assignHydrophobicity and filterHydrophobicity in succession,
then it determines from the results if transmembrane domains were found and indicates
those results accordingly on the output file.

Parameter: the single-letter amino acid sequence that results from translating the
input RNA sequence.

Parameter: access to the output file to write the discovered results.
'''
def findTMD(sequence, handle):
    # Assign hydrophobicity scores to the amino acid sequence
    hydrosequence = assignHydrophobicity(sequence)

    # Find consecutive sequences of 18-30 amino acids from the ideal
    # hydrophobicity scores
    tmdlocations = filterHydrophobicity(hydrosequence)

    if tmdlocations == '':                                                   # If no transmembrane domain locations are found
        # Indicate no transmembrane domains were found
        handle.write("Transmembrane Protein: NO\n\n\n\n")
    else:                                                                    # Otherwise, if transmembrane domain locations are found
        # Indicate transmembrane domains were found
        handle.write("Transmembrane Protein: YES\n\n")
        
        # Specify the indices of the regions where the domains were likely found
        handle.write("Transmembrane domains found at amino acid indices: \n" + tmdlocations + "\n\n\n\n")

if __name__ == "__main__":
    main()
