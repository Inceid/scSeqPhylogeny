# Suraj Joshi (2017)
# 10/3/2017

###############################################################################
####                   SUBSAMPLE TESTING FOR ALL METHODS                   ####
###############################################################################

# Importing Dependencies :

import random

###############################################################################
###                                ALL FUNCTIONS                            ###
###############################################################################


def createRandomSample(profiles, cellNames, samples, subSize):
    ''' Take in proiles, names, number of total samples, and size of '''
    ''' desired subsample. '''

    print "Creating a random subsample of %s profiles... " % subSize
    randomElems, subsample, newNames = [], [], []

    print "samples = ", samples
    allElems = range(0, samples)
    randomElems = random.sample(allElems, subSize)
    print "RANDOM ELEMS: ", randomElems

    for index in randomElems: # then add profiles / names into subsample list
        # any profile should be matched by index to its name
        profileToAdd = profiles[index]
        nameToAdd = cellNames[index]
        subsample.append(profileToAdd)
        newNames.append(nameToAdd)

    print "done."
    return (subsample, newNames)

#EOF