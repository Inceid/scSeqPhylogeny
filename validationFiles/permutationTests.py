###############################################################################
#####                                                                     #####
#####     STATISTICAL PERMUTATION TESTS FOR ALL TUMOR SAMPLE CLUSTERS     #####
#####                                                                     #####
###############################################################################

'''
updated 7-26-18

Validation Strategy:

(1) Call convert (function imported from dot_to_nw) on the tree.
This converts the tree from dot to newick format.

(2) Call Tree.get_leaves (function imported from ete3) on the newick-formatted 
tree. This collects all leaves (strings) from the tree into a list.

(3) Call getAllClusters on the leaves list. This sorts leaves into six clusters 
(lists), each cluster based on the combination of patient and region sampled 
from (e.g. 'GBM07SP1', where 'GBM07' represents the patient number, and 'SP1' 
represents the region sampled.

(4) Concatenate all clusters together into a larger list (concat). 

(5) Call getPwClusterDist on all clusters. This function calculates the sum of 
all pairwise distances between profiles within each cluster, and returns the 
cumulative sum of these cluster sums.

(6) Input a number of trials and prompt the user to proceed or quit.

(7) Shuffle concat. This randomizes all labels within concat.

(8) Call getPwClusterDist on concat, subdividing concat into six equally sized 
sublists. This recreates six completely random leaf clusters and collects new 
pairwise cluster sums for each random cluster.

(9) Append the sum of these cluster sums into a new list called 'sums'. 
Do this for all trials. 

(10) Manually calculate a p-value for actualSum given the distribution in 
'sums'. The null hypothesis is that actualSum is greater than the mean cluster
sum in 'sums'. Whereas the alternate hypothesis is that actualSum is less than
the mean cluster sum in 'sums'.

We take 0.05 to be the threshold for P-value significance.
'''

# First import dependencies

import os
import imp
from ete3 import Tree
import random as rd

# get Newick conversion function as well
from dot_to_nw import *


###############################################################################
##### HELPER FUNCTIONS FOR STATISTICAL COMPARISONS VIA PAIRWISE DISTANCES #####
###############################################################################


def getAllClusters(leaves):

    ''' Takes individual tree leaves from the original tree, and separates  '''
    '''  them into six clusters, each cluster labelled based on the prefix  '''
    '''       present in the sample. returns each cluster at the end.       '''

    # initialize empty clusters
    clu1 = clu2 = clu3 = clu4 = clu5 = clu6 = []

    # now sort each leaf into its appropriate cluster
    for i in range(len(leaves)):
        if 'G07SP1' in str(leaves[i]):
            clu1.append(leaves[i])
        elif 'G07SP2' in str(leaves[i]):
            clu2.append(leaves[i])
        elif 'G07SP3' in str(leaves[i]):
            clu3.append(leaves[i])
        elif 'G33SP1' in str(leaves[i]):
            clu4.append(leaves[i])
        elif 'G33SP2' in str(leaves[i]):
            clu5.append(leaves[i])
        else: # only thing left is 'G33SP3'
            clu6.append(leaves[i])
    
    return (clu1, clu2, clu3, clu4, clu5, clu6)


def getPwClusterDist(l1, l2, l3, l4, l5, l6):

    ''' Takes in all 6 clusters of CNP samples, each sample coming from its '''
    '''  own site and differentiated by label; returns all of the pairwise  '''
    '''       distances among the CNP's within each list of samples.        '''

    # initialize all distances
    cc1 = cc2 = cc3 = cc4 = cc5 = cc6 = 0

    for i in range(len(l1)):
        for j in range(i+1,len(l1)):
            cc1 += l1[i].get_distance(l1[j])

    for i in range(len(l2)):
        for j in range(i+1,len(l2)):
            cc2 += l2[i].get_distance(l2[j])

    for i in range(len(l3)):
        for j in range(i+1,len(l3)):
            cc3 += l3[i].get_distance(l3[j])

    for i in range(len(l4)):
        for j in range(i+1,len(l4)):
            cc4 += l4[i].get_distance(l4[j])

    for i in range(len(l5)):
        for j in range(i+1,len(l5)):
            cc5 += l5[i].get_distance(l5[j])

    for i in range(len(l6)):
        for j in range(i+1,len(l6)):
            cc6 += l6[i].get_distance(l6[j])

    return (cc1, cc2, cc3, cc4, cc5, cc6)


def getPValue(val, distribution):

    ''' manually calculates and returns P-value for val given distribution '''

    count = 0
    for i in range(len(distribution)):
        if val > distribution:
            count += 1
    numSamples = len(distribution)

    return count / numSamples


###############################################################################
################################ MAIN FUNCTION ################################
###############################################################################


def validation(testFile):
    
    deanTree = convert(testFile)
    try: deanTree = convert(testFile)
    except: raise Exception('Something wrong with the file\'s path! Exiting.')

    print 'Validation Report from %s' % testFile

    # extract all clusters and their distances from the tree
    leaves = Tree.get_leaves(deanTree)
    clu1, clu2, clu3, clu4, clu5, clu6 = getAllClusters(leaves)
    
    # concatenate clusters into a single list
    concat = clu1 + clu2 + clu3 + clu4 + clu5 + clu6
    actualDists = getPwClusterDist(clu1, clu2, clu3, clu4, clu5, clu6)
    
    # we want the sum of our actual distances for comparison
    actualSum = sum(actualDists)
    sums = []

    # ask for the number of trials for the permutation test
    trials = input('How many trials for permutation test? ')

    print 'sum of actual distances: %d' % actualSum
    runTrial = raw_input('Proceed? ')
    if runTrial == 'yes' or runTrial == 'Yes':
        pass
    else:
        print 'exited. '
        return

    increment = trials // 10 # we will print 10 trial iterations out in total

    a = len(clu1)
    b = a + len(clu2)
    c = b + len(clu3)
    d = c + len(clu4)
    e = d + len(clu5)
    f = e + len(clu6)

    for k in range(trials):
        # randomly shuffle the whole list of clusters
        rd.shuffle(concat)
        # then partition into new, randomized clusters
        (c1, c2, c3, c4, c5, c6) = getPwClusterDist(concat[0:a], concat[a:b], 
                                                    concat[b:c], concat[c:d], 
                                                    concat[d:e], concat[e:f])
        # find sums for each new concatenated array
        sums.append(c1 + c2 + c3 + c4 + c5 + c6)

        if k % increment == 0: # print each significant iteration
            print 'finished trial: ', (k + 1)

    # calculate statistics
    pValue = getPValue(actualSum, sums)

    # now display them
    print 'FINAL values below.' 
    print '=========='
    print 'pairwise distance cluster sum of experimental tree = ', actualSum
    print 'P-value = ', pValue
    print 'trials = ', trials
    print '=========='

#EOF