#!/usr/bin/env python
# coding: utf-8

from collections import namedtuple

## ----------------------------------------
## Burrows-Wheeler Aligner - Exact Matching
## Author : Muhammad Irfan
## ----------------------------------------

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
        
        # print ("C DataStructure", sorted(list(for_reference)))

        ## Initialize 2 auxillary datastructures
        C = {k:0 for k in self.alphabet}
        for_OCC = {k:[] for k in self.alphabet}

        ## create all the rotation/suffix combinations of the reference and reverse reference, and their starting index positions
        rotation_list = list()
        for i in range(len(for_reference)):
            new_rotation = "%s%s" % (for_reference[i:],for_reference[0:i])
            struct = Suffix(new_rotation, i)
            rotation_list.append(struct)
        
            ## create the C datastructure
            if for_reference[i] !='$':
                for char in self.alphabet:
                    if for_reference[i] < char:
                        C[char] = C[char] + 1  

        ## Sort the rotations/suffixes lexicographically using the suffix 
        rotation_list.sort(key = lambda x: x.text)
        # for r in rotation_list: 
        #     print (r.text)
        
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

        ## Store all the useful datastructures as class variables for easy future access
        self.SA  = suffix_array
        self.BWT = bwt
        self.C   = C
        self.Occ = for_OCC
        self.n   = len(for_reference)

    def printDataStructure(self):
        print ('-' * 10)
        for i in range(self.n):
            print (self.SA[i], self.BWT[i])
        print ('-' * 10)

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

    def OCC(self, char, index):
        return self.Occ[char][index] if index > 0 else 0

if __name__ == '__main__':

    reference = 'GGATCGATGT'
    pattern   = 'GAT'

    bwa = simpleBWA(reference)
    bwa.exact_match(pattern)
    bwa.printDataStructure()