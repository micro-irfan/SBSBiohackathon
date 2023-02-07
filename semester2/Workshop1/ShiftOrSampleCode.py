from collections import namedtuple

def _generateAlphabet(reference, query):
	alphabet = list(set(reference))
	bitap_dict = {}
	for letter in alphabet:
		letterPositionInQuery = 0
		for symbol in query:
			letterPositionInQuery = letterPositionInQuery << 1
			letterPositionInQuery |= int(letter != symbol)
		bitap_dict[letter] = letterPositionInQuery
	return bitap_dict

placeholder = namedtuple('placeholder','query seq start end mismatch')

## Shift-Or Function
def bitapSearch(reference, query, mismatch = 1):
	referenceLen = len(reference)
	queryLen = len(query)

	alphabet = _generateAlphabet(reference, query)

	matrix = [] 
	emptyColumn = (2 << (queryLen - 1)) - 1
	underground = [emptyColumn for i in range(referenceLen + 1)]
	matrix.append(underground)
	gRNAs = []
	skip = []

	for k in range(1, mismatch + 2):
		matrix.append([emptyColumn])

		for columnNum in range(1, referenceLen + 1):
			prevColumn = (matrix[k][columnNum - 1]) >> 1
			letterPattern = alphabet[reference[columnNum - 1]]
			curColumn = prevColumn | letterPattern

			if k > 1:
				#insertColumn = curColumn & (matrix[k - 1][columnNum - 1])
				#deleteColumn = curColumn & (matrix[k - 1][columnNum] >> 1)
				#replaceColumn = curColumn & (matrix[k - 1][columnNum - 1] >> 1)
				#curColumn = insertColumn & deleteColumn & replaceColumn
				
				## Mismatch 
				curColumn = curColumn & (matrix[k - 1][columnNum - 1] >> 1)
				
			matrix[k].append(curColumn)

			if (curColumn & 0x1) == 0:
				startPos = max(0, columnNum - queryLen) # taking in account Replace operation
				if startPos in skip: continue
				endPos = min(columnNum, referenceLen) # taking in account Replace operation
				place = reference[startPos:endPos]
				temp = placeholder(query, place, startPos, endPos, k - 1)
				gRNAs.append(temp)
				skip.append(startPos)
			
	return gRNAs

reference  = 'GGGCNCTGCTGAGAATGNACTGAATATAAACTTGTGGTAGTTGGANGCTGGTGGCGTAGGCTTGTGGTTGTGGGANGCTGGTGGCGAAGAGTGCCTTGACGATACAGNCTANATTNCAGAATNCATTTTGTGGNACGAATATGATCCANACAATAGNAGGATTC'
string_search = 'TTGTGGTAGTTGGANGCTGGTGGCG'
errors = 2
gRNAs = bitapSearch(reference, string_search, errors)
for i, g in enumerate(gRNAs):
	print (g)

# Will Print: 
# placeholder(query='TTGTGGTAGTTGGANGCTGGTGGCG', seq='TTGTGGTAGTTGGANGCTGGTGGCG', start=31, end=56, mismatch=0)
# placeholder(query='TTGTGGTAGTTGGANGCTGGTGGCG', seq='TTGTGGTTGTGGGANGCTGGTGGCG', start=61, end=86, mismatch=2)
