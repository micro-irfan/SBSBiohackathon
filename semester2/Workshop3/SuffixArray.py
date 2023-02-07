#!/usr/bin/env python
# coding: utf-8

## -----------
## SuffixArray
## Author : Muhammad Irfan
## -----------

wordOfTheDay = 'CAMEL'

## Store (index, sliced word)
listOfSuffix = []

## Loop & Slice Words
for i in range(len(wordOfTheDay)):
	slicedWord = wordOfTheDay[i:]
	listOfSuffix.append((i, slicedWord))

## Sort Lexicographically
listOfSuffix.sort(reverse = False, key = lambda x: x[1])  

## Retrieve the indices
suffixArray = [index for index,suffix in listOfSuffix]

print ('-' * 10)
for i in listOfSuffix:
	print (i)
print ('-' * 10)
print (suffixArray)

'''
Expected Output for Camel:
----------
(1, 'AMEL')
(0, 'CAMEL')
(3, 'EL')
(4, 'L')
(2, 'MEL')
----------
[1, 0, 3, 4, 2]

'''