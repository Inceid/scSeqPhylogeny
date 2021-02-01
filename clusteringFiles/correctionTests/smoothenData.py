# Suraj Joshi (2017)
# 9/24/2017

###############################################################################
####         DATA-SMOOTHING : SIMPLE MOVING AVERAGE IMPLEMENTATION         ####
###############################################################################

# Importing Dependencies :

import random

###############################################################################
###                                ALL FUNCTIONS                            ###
###############################################################################

def movingAvg(profiles, width): 

    '''   Takes in a data set of profiles, and even width (width % 2 == 0) '''
    '''   and returns data set where each profile is a weighted average of '''
    '''     it and the profiles neighboring it, over the entire data set.  '''

    print "Calculating moving-average of data set... "
    newSet = [0] * len(profiles)
    stretch = width / 2
    maxIndex = len(profiles) - 1
    # make sure the segment accesses are all safe:
    minBound = 0
    maxBound = len(profiles)

    for index in range(maxBound):
        # get the range / window of profiles to calculate avg
        currSlice = profiles[max(minBound, index - stretch) : 
                             min(index + stretch, maxBound) + 1]
        # add +1 to right range of slice (slice ranges are [min, max)!)
        # and then calculate avg, place into new set
        newSet[index] = getAverageProfile(currSlice)

    return newSet


def getAverageProfile(profiles):
    ''' Takes in a list of copy number profiles, returns a profile where  '''
    ''' every index i contains the integer average of all copy numbers at '''
    '''                  i from the list of profiles                      '''
    averageProfile = []
    profileLength, numProfiles = len(profiles[0]), len(profiles)

    for i in xrange(profileLength):
        average = 0
        for profile in profiles:
            average += profile[i]
        average //= numProfiles
        averageProfile.append(average)

    return averageProfile


def testGetAverageProfile():

    profiles1 = [ [1,2,1], [2,2,1], [4,7,10] ]
    assert(getAverageProfile(profiles1) == [2, 3, 4])

    profiles2 = [ [1, 4, 21, 21, 9], [2, 3, 0, 12, 9], [1, 5, 0, 6, 9], 
                  [3, 1, 6, 8, 9], [2, 4, 0, 5, 11] ]
    assert(getAverageProfile(profiles2) == [1, 3, 5, 10, 9])
    
    print "testGetAverageProfile passed!"


def testMovingAvg():

    oldSet = [ [1, 0, 1], [5, 5, 5], [3, 1, 2], 
               [2, 1, 3], [5, 5, 5], [0, 0, 0] ]
    window = 3
    newSet = movingAvg(oldSet, window)
    assert newSet == [ [3, 2, 3], [3, 2, 2], [3, 2, 3], 
                       [3, 2, 3], [2, 2, 2], [2, 2, 2] ]

    print "testMovingAvg passed!"

#EOF