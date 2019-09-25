'''
DEFINE FUNCTION HERE
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
            windowSum += code[sequence[j]]          # sum the hydrophobicity values of the window
        if(windowSum < 5):                          # if the sum is less than our threshold of 5 . . .
            hydrosequence.append(i)                 # store the index of that low hydrophobicity region
        windowSum = 0
        
    # Return the list of indices containing hydrophobicity value sums less than our threshold
    return hydrosequence
    
'''
DEFINE FUNCTION HERE
'''
def filterHydrophobicity(hydrosequence):
    tempstart = 0
    length = 0

    for i in range(len(hydrosequence) - 1):
        if ((hydrosequence[i+1] - hydrosequence[i]) == 1 and i != len(hydrosequence) - 2):    # If there is not a gap in hydrophobicity scores
            if tempstart == 0:
                tempstart = hydrosequence[i]                                        # Store the possible starting index
            length += 1
            print(tempstart)
        else:
            if length > 6 and length < 20:
                print(tempstart, hydrosequence[i] + 11)
            tempstart = 0
            length = 0
            
# THIS NEEDS TO BE FIXED: you don't have it correctly selecting the chunks. the business of getting it to select the gap-less intervals
# appears to work fine with the mechanism of restarting, but the business of getting things within the ideal range of 18-30 amino acids
# doesn't appear to work how we want it to so need to fix that
            
        