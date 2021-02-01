###############################################################################
####                                                                       ####
####              MODIFYING DATA TO FIT FISHTREES PARAMETERS               ####
####                                                                       ####
###############################################################################

# Importing Dependencies :

from __future__ import division # we need float division!
from cntp_final import parseDeanTiming, getCopyNumberDirec

import os
import time
import copy
import sys
import traceback
import random
import string
from collections import Counter


def getData():

    ''' retrieves and returns all copy number data. '''

    print 'Getting data... \n'
    file = getCopyNumberDirec()
    profiles, cells, regions, names = parseDeanTiming(file)
    return profiles, cells, regions, names


def pickRandomProbes(profiles, cells, regions, names, marked, n):

    ''' takes data from parseDeanTiming & list of marked profile indices, '''
    ''' and 'n' number of indices to select from data. A marked index has '''
    ''' already been picked on a previous trial. '''
    ''' Returns modified profile set limited to 'n' selected indices,   '''
    ''' a list of the chosen indices, and updated list of marked indices. '''

    possibleChoices = []
    for index in range(len(regions)):
        if index not in marked:
            possibleChoices.append(index)

    # need to deal with possibleChoices being too small
    if possibleChoices == []: return 'all indices have already been picked!'
    elif len(possibleChoices) <= n: randIndices = possibleChoices
    else: # size is large enough; sample randomly from possibleChoices.
        randIndices = random.sample(possibleChoices, n)

    sampleProfiles = []
    for profile in profiles:
        # build the new profile
        newProfile = []
        for randIndex in randIndices:
            newProfile.append(profile[randIndex])
        sampleProfiles.append(newProfile)

    marked.append(randIndices) # update marked destructively
    return sampleProfiles, randIndices


def initializeHeader(numRandomProbes):

    ''' function that returns initial header line for FISHtrees text files. '''
    ''' Single cell data is assumed to have no chromosome probes, therefore '''
    ''' we take chromosome probe count to be '1'. This probe represents the '''
    ''' 'PLOIDY' of each cell. Then we add the number of random probes that '''
    '''                 were sampled from pickRandomProbes.                 '''

    return ['1', 'PLOIDY', str(numRandomProbes)]


def createHeader(numRandomProbes):

    ''' function that creates header text for each FISHtrees data file. '''

    header = initializeHeader(numRandomProbes) # create the header string
    for probe in range(0, numRandomProbes): # append each probe name
        header.append(string.ascii_uppercase[probe])
    return header


def convertProfilesToTuples(sampleProfiles):

    ''' converts all profiles in sampleProfiles to tuples of themselves. '''
    ''' e.g. [2, 2, 2] becomes (2, 2, 2). '''

    for i in range(len(sampleProfiles)):
        profile = sampleProfiles[i] # grab profile
        for j in range(len(profile)):
            profile[j] = str(profile[j]) # replace copy numbers with strings
        sampleProfiles[i] = tuple(profile) # replace with corresponding tuple


def writeTextFile(sampleProfiles, numCells, numRandomProbes, name):

    ''' function that outputs a .txt file named 'name'. '''

    print 'writing a text file... \n'
    fileName = "%s.txt" % name # create new file
    newFile = open(fileName,'w')

    # create and write header to file
    header = createHeader(numRandomProbes)
    newFile.write('\t'.join(header) + '\n')

    # convert each profile in sampleProfiles into a tuple. This makes the
    # profiles immutable, so that they can become dictionary keys as below.
    convertProfilesToTuples(sampleProfiles)

    # retrieve dictionary, each element being (profile, count) where 'count'
    # is the number of occurrences of 'profile' in sampleProfiles.
    uniqueProfilesDict = Counter(sampleProfiles)
    numUniqueProfiles = len(uniqueProfilesDict)

    # now we can write the second header line
    newFile.write(str(numUniqueProfiles) + '\t' + str(numCells))

    for profile in uniqueProfilesDict:
        count = uniqueProfilesDict[profile]
        newFile.write('\n')
        newFile.write('2\t') # assume ploidy is always 2 for single-cell data
        # now write profile + tab + count
        newFile.write('\t'.join(profile) + '\t' + str(count))
    newFile.close()

    return fileName


def runPreProcessing():

    ''' runs all preprocessing trials, setting up FISHtrees data '''

    marked = []
    profiles, cells, regions, names = getData()
    # create new directories to fill with sample profiles
    # each directory will eventually represent a patient, but just have
    # one directory for now for all cell samples from both patients
    separated_data = 'processed_data'
    os.mkdir(separated_data)
    numRandomProbes = 10 # number of probes to select each trial.
                         # hardcoded to 10 for now
    # total number of preprocessing trials:
    randTrials = 10 # change later to #int(round(regions / numRandomProbes))

    oldPath = os.getcwd() # first, get old path
    print 'Current path: %s' % oldPath
    # now switch to directory where you want to write text files
    os.chdir(separated_data)

    for trial in range(1, randTrials+1):
        print 'trial: ', trial
        sampleProfiles, indices = pickRandomProbes(profiles, cells,
                                                   names, regions, marked,
                                                   numRandomProbes)
        # now convert to FISHtrees format
        name = 'subset_%s' % str(trial)
        print 'number of profiles: ', len(sampleProfiles)
        writeTextFile(sampleProfiles, cells, numRandomProbes, name)
        # note : need to deal with case where len(sampleProfile) < 10 !
    os.chdir(oldPath)
    print 'done'


def main():
    runPreProcessing()


if __name__ == '__main__':
    main()