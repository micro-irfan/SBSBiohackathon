from collections import namedtuple

#query: PATTERN
#reference: SEQUENCE

#Generating the bitmask B[Tj] for the elements found in reference
def _generateAlphabet(reference, query): #defining function, reference is the seq and query is the pattern
	alphabet = list(set(reference)) #seperate the ref into iterable UNIQUE elements
	bitap_dict = {} #opening up a bitmask dictionary B[Tj] for each unique letter
	for letter in alphabet: #for every single unique element in the reference
		letterPositionInQuery = 0 #start off with 0
		for symbol in query:
			letterPositionInQuery = letterPositionInQuery << 1 #left shift to move specific letter across binary string (from left to right) to tally with the letter position in the query
			letterPositionInQuery |= int(letter != symbol) ## adds 1/0 to the end of binary string based on boolean value #check if the specified character matches corresponding character in the query in the specific position
		bitap_dict[letter] = letterPositionInQuery #repeat above block till last character in the query and update dictionary with the obtained bit mask of len (equals to query length)
	return bitap_dict #dictionary of <all possible letters (from seq/ref) in query/pattern> (key) and <letterPositionInQuery> (value)

    
placeholder = namedtuple('placeholder','query seq start end mismatch')
#placeholder is the name, and query seq start end mismatch is the associated value; better than dictionary since it is immutable and less writing is involved

## Shift-Or Function
def bitapSearch(reference, query, mismatch = 1): #NOTE: mismatch refers to a mismatch of a single character between query and reference [number of mismatches defined refers to number of erroneous characters tolerated for a pattern to be found]
	referenceLen = len(reference) #sequence
	queryLen = len(query) #pattern

	alphabet = _generateAlphabet(reference, query)

	matrix = [] 
	emptyColumn = (2 << (queryLen - 1)) - 1 #binary code of 2 (10) left shift (queryLen-1) spaces (the minus 1 to convert all binaries to 1 and reduce length of binary string by 1)
    #e.g. if queryLen is 4, 
    #2 has a binary string of 10 
    #2<<(4-1) --> 0b10 000
    # 0b10 000 -1 --> 0b1111 (binary string length reduce by 1 --> same length as queryLen)
    
	underground = [emptyColumn for i in range(referenceLen + 1)] #underground is a list of 11111 (D) inital bitmasks for every letter in the reference (sequence)
	matrix.append(underground)
	gRNAs = []
	skip = []

	for k in range(1, mismatch + 2): 
		matrix.append([emptyColumn])
        #in the matrix:
            #index [0]: underground
            #index [1]: no mismatch (thats why the inexact matching only occurs for k>1, when k=1, the if clause does not run)
            #index [2]: 1 mismatch 
            #index [k]: k-1 mismatches
            #NOTE: k has a range of 1, mismatch+2 because k starts at 1, since index 0 of matrix is used for underground, and mismatch specified plus 2 because the range function does not run the max number stated
            #e.g. mismatch = 1, so k runs in range (1,3) --> k = 1, matrix[1]: stores no mismatch, k = 2, matrix[2]: stores 1 mismatch
            
		for columnNum in range(1, referenceLen + 1):
			prevColumn = (matrix[k][columnNum - 1]) >> 1 #D'#right shift to move on to next letter 
			letterPattern = alphabet[reference[columnNum - 1]] #bitmask B[Tj] #Note that the letter index in reference (seq) will always be 1 less than columnNum (since columnNum starts with 1 to denote first letter in seq, which has index 0)
			curColumn = prevColumn | letterPattern

			if k > 1: #INEXACT STRING SEARCH
				## Mismatch 
				curColumn = curColumn & (matrix[k - 1][columnNum - 1] >> 1)
				#(matrix[k - 1][columnNum - 1] >> 1) looks at the previous result for the same character in the same position, AND conditional on the current column and no mismatch D (matrix[1][columnNum-1])
                
                ##IMPORTANT NOTE: 
                    #OR: favours 1 (where 1,1  0,1  1,0 --> gives 1 and only 0,0 --> gives 0)
                    #OR function is more sensitive to mismatches --> hence used in exact string matching
                    #AND: favours 0 (where 0,0  0,1  1,0 --> gives 0 and only 1,1 --> gives 1)
                    #AND function tolerates mismatches --> hence used in inexact string matching
                    
			matrix[k].append(curColumn)

			if (curColumn & 0x1) == 0: #0x1 represents hexadecimal 1 (binary string of 1) #in a right shift matching --> 0 in first position implies a match in first character (note that & 0x1 --> implies that curColumn must have a 0 in first posiiton to offset the 1 in 0x1)
				startPos = max(0, columnNum - queryLen) # taking in account Replace operation #ensures that the minimum length of a match is the query length
				if startPos in skip: continue
				endPos = min(columnNum, referenceLen) # taking in account Replace operation #finds the shortest possible instance of a match
				place = reference[startPos:endPos]
				temp = placeholder(query, place, startPos, endPos, k - 1)
				gRNAs.append(temp)
				skip.append(startPos) #prevents repeated reporting of same strain
                #NOTE: inexact string search will also capture previous patterns where there are no mismatches and instances of fewer mismatches --> repeated reporting
			
	return gRNAs

#Implementation of defined functions:
reference  = 'GGGCNCTGCTGAGAATGNACTGAATATAAACTTGTGGTAGTTGGANGCTGGTGGCGTAGGCTTGTGGTTGTGGGANGCTGGTGGCGAAGAGTGCCTTGACGATACAGNCTANATTNCAGAATNCATTTTGTGGNACGAATATGATCCANACAATAGNAGGATTC'
string_search = 'TTGTGGTAGTTGGANGCTGGTGGCG' #pattern we are looking for
errors = 2 #allowing mismatch of 2; subsitution mutations to take place
gRNAs = bitapSearch(reference, string_search, errors)
for i, g in enumerate(gRNAs):
	print (g)

# Will Print: 
# placeholder(query='TTGTGGTAGTTGGANGCTGGTGGCG', seq='TTGTGGTAGTTGGANGCTGGTGGCG', start=31, end=56, mismatch=0)
# placeholder(query='TTGTGGTAGTTGGANGCTGGTGGCG', seq='TTGTGGTTGTGGGANGCTGGTGGCG', start=61, end=86, mismatch=2)
