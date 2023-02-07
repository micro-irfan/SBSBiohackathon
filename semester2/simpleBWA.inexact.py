#!/usr/bin/env python
# coding: utf-8

from collections import namedtuple

## ------------------------------------------
## Burrows-Wheeler Aligner - Inexact Matching
## Author : Muhammad Irfan
## ------------------------------------------

Suffix = namedtuple('Suffix', 'text pos')

class simpleBWA:
    """
    A Burrows-Wheeler Aligner class. 
    Contains the 4 core datastructures used in the algorithm
    SA, BWT, C and Occ which are created in the constructor which is passed the reference string as an argument
    """

    # Initialize the Datastructure 
    def __init__(self, reference):
        reference = reference.upper()

        self.alphabet = list(set(reference)) 
        for_reference = reference + '$'
        rev_reference = reference[::-1] + '$'

        ## Initialize 2 auxillary datastructures
        C = {k:0 for k in self.alphabet}
        for_OCC = {k:[] for k in self.alphabet}
        rev_OCC = {k:[] for k in self.alphabet}

        ## create all the rotation/suffix combinations of the reference and reverse reference, and their starting index positions
        rotation_list = list()
        rotation_list_reverse = list()
        for i in range(len(for_reference)):
            new_rotation = "%s%s" % (for_reference[i:], for_reference[0:i])
            struct = Suffix(new_rotation, i)
            rotation_list.append(struct)

            new_rotation_reverse = "%s%s" % (rev_reference[i:], rev_reference[0:i])
            struct_rev = Suffix(new_rotation_reverse,i)
            rotation_list_reverse.append(struct_rev)
        
            ## create the C datastructure
            if for_reference[i] !='$':
                for char in self.alphabet:
                    if for_reference[i] < char:
                        C[char] = C[char] + 1  

        ## Sort the rotations/suffixes lexicographically using the suffix 
        rotation_list.sort(key = lambda x: x.text)
        rotation_list_reverse.sort(key = lambda x: x.text)
        
        ## Initialize DataStrucutre
        suffix_array, bwt = [list() for i in range(2)]

        ## BWT and SA construction
        for i in rotation_list:
            ## Position of the reordered suffixes forms the Suffix Array elements
            suffix_array.append(i.pos)

            ## Last character in each rotation (in the new order) forms the BWT string elements
            bwt.append(i.text[-1])
        
            ## Occ datastructure Construction
            for char in self.alphabet:
                if len(for_OCC[char]) == 0:
                    prev = 0
                else:
                    prev = for_OCC[char][-1]
                
                if i.text[-1:] == char:
                    for_OCC[char].append(prev + 1)
                else:
                    for_OCC[char].append(prev)

        for i in rotation_list_reverse:
            ## Occ datastructure Construction for reverse reference
            for char in self.alphabet:
                if len(rev_OCC[char]) == 0:
                    prev = 0
                else:
                    prev = rev_OCC[char][-1]
                if i.text[-1:] == char:
                    rev_OCC[char].append(prev+1)
                else:
                    rev_OCC[char].append(prev)  

        ## Store all the useful datastructures as class variables for easy future access
        self.SA  = suffix_array
        self.BWT = bwt
        self.C   = C
        self.Occ = for_OCC
        self.n   = len(for_reference)

        ## Occ datastructure for the reverse reference, using to construct the D array (the lower bound on the number of differences allowed), to speed up alignments 
        self.Occ_reverse = rev_OCC 

        ## empty list for later use
        self.D = list()

    def printDataStructure(self):
        print ('-' * 10)
        for i in range(self.n):
            print (self.SA[i], self.BWT[i])
        print ('-' * 10)

    ## Get the position(s) of the query in the reference
    def find_match(self, query, mismatch):
        query = query.upper()
        if mismatch == 0:
            print(self.exact_match(query))
        else:
            print(self.inexact_match(query, mismatch))

    ## Exact matching - no indels or mismatches allowed
    def exact_match(self, query):
        query = query.upper()
        i = 0
        j = self.n - 1
        
        for x in range(len(query)):
            newChar = query[-x-1]
            newI = self.C[newChar] + self.OCC(newChar, i-1) + 1
            newJ = self.C[newChar] + self.OCC(newChar, j)
            i = newI
            j = newJ

        matches = self.SA[i:j+1]
        print(matches)

    ## Inexct matching - based on algorithm explained in bwa paper
    ## inexact matching, z is the max threshold for allowed edits
    def inexact_match(self, query, z):
        self.calculate_d(query)
        SA_indices = self.inexact_recursion(query, len(query)-1, z, 0, self.n-1)
        ## Return the values in the SA
        return [self.SA[x] for x in SA_indices]

    ## recursion function that effectively "walks" through the suffix tree using the SA, BWT, Occ and C datastructures
    def inexact_recursion(self, query, i, z, k, l):
        tempset = set()
        use_lower_bound_tree_pruning = True
        indels_allowed = False
        insertion_penalty = 1
        deletion_penalty = 1
        mismatch_penalty = 1

        ## 2 stop conditions, one when too many differences have been encountered, another when the entire query has been matched, terminating in success
        ## Reached the limit of differences at this stage, terminate this path traversal
        if (z < self.get_D(i) and use_lower_bound_tree_pruning) or (z < 0 and not use_lower_bound_tree_pruning):
            ## Return empty set   
            return set()
        
        ## Empty query string, entire query has been matched, return SA indexes k:l
        if i < 0:
            for m in range(k,l+1):
                tempset.add(m)
            return tempset
            
        result = set()
        if indels_allowed: result = result.union(self.inexact_recursion(query,i-1,z-insertion_penalty,k,l))#without finding a match or altering k or l, move on down the query string. Insertion
        ## For each character in the alphabet
        for char in self.alphabet:
            ## Find the SA interval for the char
            newK = self.C[char] + self.OCC(char,k-1) + 1 
            newL = self.C[char] + self.OCC(char,l)
            
            ## If the substring was found
            if newK <= newL:
                if indels_allowed: result = result.union(self.inexact_recursion(query,i,z-deletion_penalty,newK,newL))# Deletion
                
                ## If the char was correctly aligned, then continue without decrementing z (differences)
                if char == query[i]: result = result.union(self.inexact_recursion(query,i-1,z,newK,newL))
                
                ## Continue but decrement z, to indicate that this was a difference/unalignment
                else: result = result.union(self.inexact_recursion(query,i-1,z-mismatch_penalty,newK,newL))
        
        return result

    #calculates the D array for a query, used to prune the tree walk and increase speed for inexact searching
    def calculate_d(self,query):
        k = 0
        l = self.n-1
        z = 0

        ## Empty the D array
        self.D = list()
        for i in range(len(query)):
            k = self.C[query[i]] + self.OCC(query[i], k-1, reverse=True) + 1
            l = self.C[query[i]] + self.OCC(query[i], l, reverse=True)
            
            ## If this character has NOT been found
            if k > l:
                k = 0
                l = self.n - 1
                z = z + 1
            self.D.append(z)

    ## returns normal Occ value, otherwise returns the reverse Occ if explicitly passed as an argument
    ## NOTE Occ('a',-1) = 0 for all 'a'
    def OCC(self, char, index, reverse=False):
        if index < 0:
            return 0
        else:
            if reverse:
                return self.Occ_reverse[char][index]
            else:
                return self.Occ[char][index]
    
    ## gets values from the D array
    ## NOTE D(-1) = 0
    def get_D(self, index):
        if index < 0:
            return 0
        else:
            return self.D[index]

if __name__ == '__main__':

    reference = 'MISSISSIPPI'
    pattern   = 'IMS'

    reference  = 'GGGCCTGCTGAGAATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAGAATCATTTTGTGGNACGAATATGATCCAACAATAGAGGATTC'
    pattern = 'TAAACTTGTGGTAGTTGGGGCTGG'

    bwa = simpleBWA(reference)
    
    bwa.printDataStructure()
    bwa.find_match(pattern, 1)