
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
            windowSum += code[sequence[j]]          # sum the hydrophobicity values of the window
        if(windowSum < 5):                          # if the sum is less than our threshold of 5 . . .
            hydrosequence.append(i)                 # store the index of that low hydrophobicity region
        windowSum = 0
        
    # Return the list of indices containing hydrophobicity value sums less than our threshold
    return hydrosequence
    
