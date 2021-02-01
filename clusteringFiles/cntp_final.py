# Suraj Joshi (2017)
# 6/28/2017

###############################################################################
#####                                                                     #####
#####         IMPLEMENTING SHAMIR ET AL CNT ALGORITHM ON FISH DATA        #####
#####                                                                     #####
###############################################################################

# # # # #
# L O G #
# # # # #

'''
6/18/17 - Completed! Will run to see what we get!

6/19/17 - runtime is a serious issue. within 5 minutes the distance matrix
was not able to be calculated. Must fix!

6/26/17 - Changed distance function such that it works on all test cases; gives
out results similar to the polytimeDistance function

6/28/17 - Cleaned up code, replaced hardcoded paths and file names with
user-and-OS-flexible names. This should make the overall program more robust.
Sent to Prof. Schwartz.

 === took a break ===

8/29/17 - Returned to lab, cleaning up code

'''

#   Variable names correspond to Shamir, Zehavi, and Zeira's paper from TAU   #
# =========================================================================== #
#                  My variable names -> their variable names                  #
# =========================================================================== #
#                                Src -> S                                     #
#                                Tgt -> T                                     #
#                                  n -> n                                     #
#                                  N -> N                                     #
#                           deletion -> d                                     #
#                             mTable -> M                                     #
#                                amp -> a[i, d] for all indices i             #
#                        prevIndices -> prev(i) for all indices i             #
#                             qArray -> Qi for all indices i                  #
#                               difs -> u                                     #
#                           baseVals -> base[i] for all indices i             #
#                               a, b -> a, b                                  #
# =========================================================================== #

# The paper itself is in scSeqPhylogeny/_other/papers

# Importing Dependencies :

import os
import time
import copy
import sys
import traceback

# To use cmd's dot-to-png conversion command from python
try: from subprocess import check_call
except: print 'warning: failed to import check_call'

# We would like to use functions from Tyler's script
import final_nj as nj

# We also want to import our correction tests
from correctionTests.subsampleTests import *
from correctionTests.smoothenData import *

# [SHAMIR'S SAMPLE PROBLEM]
SRC_EX = [3, 1, 2, 3, 2, 1, 4]
TGT_EX = [2, 0, 0, 0, 0, 0, 2]

# [MY SAMPLE PROBLEM]
SRC_EX2 = [3, 1, 2, 3, 2, 1, 4]
TGT_EX2 = [2, 0, 1, 0, 1, 1, 3]


###############################################################################
###                         ALL DATA-PARSING FUNCTIONS                      ###
###############################################################################


def getCopyNumberDirec():

    '''   Finds, prints, and returns the path containing our raw CN data   '''
    ''' moves to scSeqPhylogeny direc, then navigates to \allcp from there '''

    # topDirec is scSeqPhylogeny, two levels up from current file
    topDirec = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    cpDirec = topDirec + '\\deanFiles\\allcp\\'
    print 'location of all copy number files: ', cpDirec
    return cpDirec


def parseDeanTiming(direc):

    ''' slightly edited parseDean ; this prints out total time taken '''

    profiles, cells, names = [], 0, []
    startTime = time.time()

    for filename in os.listdir(direc):

        start = filename.index('SC')
        end = filename.index('.')
        cellname = 'G'+filename[3:5]+'SP'+filename[8:9]+filename[start:end]
        names.append(cellname)
        
        f = open(direc + filename)
        fl = f.readlines()
        cpArr, copynumbers, regions = [], fl[1:], 0

        for i in xrange(len(copynumbers)):
            cpArr.append(int(copynumbers[i].replace('\n','')))
            regions += 1

        profiles.append(cpArr)
        f.close()
        cells += 1

    endTime = time.time()
    print 'parseDeanTiming took %d seconds.' % (endTime-startTime)
    return (profiles, cells, regions, names)


def writeCNTPDotFile(t):

    ''' function that outputs a dotfile. Asks for output filename. '''
    '''        filename should be a string such as "test01"        '''
    '''              dot file can be read with GVEdit              '''

    dataname = raw_input("output filename? ")
    fileName = "CNTP_%s.dot" % dataname # create new file
    f = open(fileName,'w')

    # read and concatenate entries of instructions
    text = "digraph tree_of_%s {\n" % dataname + "".join(t) + "}"
    f.write(text)
    f.close()

    return fileName


def runCNTP(profiles, cells, regions, names):
    
    ''' takes file path, performs CNTP tree-joining and dotFile writing. '''

    print 'running CNTP... \n',   
    distMatrix, taxaL, taxaL2 = linearDistanceMatrix(profiles, cells, regions,
                                                     names)
    t = nj.neighborJoin(distMatrix, cells, taxaL, taxaL2, names)
    fileName = writeCNTPDotFile(t) # create and write into a new dotFile
    fileName = fileName.replace('.dot', '')
    
    try: #now try to convert to png within program
        print 'command line is converting your .dot file... ' 
        check_call(['dot','-Tpng', fileName + '.dot','-o', fileName + '.png'])
    except: print 'command-line call failed. Convert the .dot file manually.'
    
    print 'Your tree (.dot and .png) files are in cntp_final\'s directory.'


def polytimeDistanceMatrix(profiles, cells, regions, names):

    ''' very similar to nj.distanceMatrixL1 in final_nj_Suraj, except that '''
    '''       the distance is calculated using the CNTP approach,       '''
    '''            as opposed to calculating pairwise sums.             '''

    print 'calculating polynomial time distance matrix... '
    distMatrix = dict()
    taxaL = []

    for i in xrange(cells):
        taxaL.append(str(names[i]))
        for j in xrange(i+1,cells):
            # not counting distance between L[i] and itself
            distance = polytimeDistance(profiles[i], profiles[j])
            distMatrix[str(names[i]),str(names[j])] = distance
    print 'got the matrix!'
    return (distMatrix, taxaL, copy.copy(taxaL))


def linearDistanceMatrix(profiles, cells, regions, names):

    ''' linear time version of polytimeDistanceMatrix '''

    print 'calculating linear time distance matrix... '
    beginning = time.time()
    distMatrix = dict()
    taxaL = []
    for i in xrange(cells):
        # keep track of our progress
        print 'finished %d out of %d indices' % (i+1, cells)
        taxaL.append(str(names[i]))
        for j in xrange(i+1,cells):
            # not counting distance between L[i] and itself
            # here is the only difference: use linear instead of polytime
            distance = linearDistance(profiles[i], profiles[j])
            distMatrix[str(names[i]),str(names[j])] = distance
    print 'got the matrix!'
    end = time.time()
    print 'this took %d seconds. ' % (end - beginning)
    return (distMatrix, taxaL, copy.copy(taxaL))


def runALL(profiles, cells, regions, names):

    ''' Runs all distance methods for neighbor joining, returns tree files '''

    print "Running all tests... "

    print "Doing L1... "
    LD = 'L1'
    (dM, taxaL, taxaL2) = nj.distanceMatrixL1(profiles, cells, regions, names)
    t = nj.neighborJoin(dM, cells, taxaL, taxaL2, names)
    fileName = nj.writeDotFile(t, LD) # create and write into a new dotFile
    fileName = fileName.replace('.dot', '')
    try: # now try to convert to png within program
        print 'command line is converting your .dot file... '
        check_call(['dot','-Tpng', fileName + '.dot','-o', fileName + '.png'])
    except: 'command-line call failed. Convert the .dot file manually.'

    print "Doing L2... "
    LD = 'L2'
    (dM, taxaL, taxaL2) = nj.distanceMatrixL2(profiles, cells, regions, names) 
    t = nj.neighborJoin(dM, cells, taxaL, taxaL2, names)
    fileName = nj.writeDotFile(t, LD) # create and write into a new dotFile
    fileName = fileName.replace('.dot', '')
    try: # now try to convert to png within program        
        print 'command line is converting your .dot file... '
        check_call(['dot','-Tpng', fileName + '.dot','-o', fileName + '.png'])
    except: 'command-line call failed. Convert the .dot file manually.'

    print "Doing CNTP... "
    runCNTP(profiles, cells, regions, names)


def runTrialPrompts(profiles, cells, regions, names):

    (profiles, cells, 
     regions, names) = chooseDataSet(profiles, cells, 
                                     regions, names)
    chooseSmoothing(profiles)
    print cells
    chooseTest(profiles, cells, regions, names)


def chooseDataSet(profiles, cells, regions, names):

    # ask what part of data to use
    print("Use all data or a subset? ")
    dataToUse = raw_input("Type in \'all\' or \'sub\', then press Enter. ")

    if dataToUse == 'sub':
        subSize = input("How many profiles in subsample? ")
        print "evaluating subsample... "
        # regions remains the same as before, but everything else changes
        profiles, names = createRandomSample(profiles, names, cells, subSize)
        cells = subSize # change cells only after creating random sample
    elif dataToUse == 'all': 
        # keep data set the same
        print "evaluating full data set... "
    else: raise ValueError("Error: should either choose \'sub\' or \'all\'!")
    return profiles, cells, regions, names


def chooseSmoothing(profiles):

    choice = raw_input("Smooth data? ")

    if choice == "yes":
        # choose the width
        width = input("Choose window size for smoothing. ")
        profiles = movingAvg(profiles, width)


def chooseTest(profiles, cells, regions, names):

    print "Which test? "
    choice = raw_input("Type ALL, L1, L2, or CNTP, then press Enter. ")
    if choice == 'CNTP': runCNTP(profiles, cells, regions, names)
    elif choice == 'ALL': runALL(profiles, cells, regions, names)
    elif choice in ['L1', 'L2']: nj.runL1L2(profiles, cells, regions, 
                                            names, choice)
    else: raise ValueError("Error: test needs to be L1, L2, or CNTP!")


###############################################################################
###                     DISTANCE-EVALUATING FUNCTIONS!                      ###
###############################################################################


def polytimeDistance(Src, Tgt):

    '''         takes in source, target profiles        '''
    '''  finds and returns CNTP-distance between them   '''
    '''  Uses first O(n x N^2) method in Shamir et al.  '''

    n, N, prevIndices, qArray, difs = calcInitVars(Src, Tgt)
    tableInit = mTableInit(Src, Tgt, difs, n, N, prevIndices, qArray)
    tableFull = mTableFull(Src, Tgt, difs, n, N, prevIndices, 
                           qArray, tableInit)
    cntpDistance = getFinalValue(tableFull)
    return cntpDistance


''' Note: the function linearDistance(), which is the counterpart to this '''
'''   function above, is in the complexity reduction section as of now.   '''


###############################################################################
###                          ALL HELPER FUNCTIONS!                          ###
###############################################################################


def calcInitVars(Src, Tgt):

    '''   takes in source, target profiles, calculates and returns  '''
    ''' iteration bounds, raw differences, prevIndices, and qArray. '''
    '''             these variables are explained below.            '''

    # iteration bounds
    n = len(Src)
    assert(n == len(Tgt))
    N = max( max(Src), max(Tgt)) # max copy number

    # pre-computed arrays
    prevIndices = getPrevIndices(Tgt)
    qArray = getQArray(prevIndices, Src)
    difs = [ (Src[i] - Tgt[i]) for i in xrange(len(Src)) ]

    return (n, N, prevIndices, qArray, difs)


def mTableFull(Src, Tgt, difs, n, N, prevIndices, qArray, tableInit):

    '''             temporary function only for demonstration            '''
    '''       whereas the algorithm stores only two rows of mTable,      '''
    '''               this function shows the whole mTable.              '''

    table = [copy.copy(tableInit)] #full table is 2d list
    for i in xrange(1, n): 
        if Tgt[i] == 0: 
            col = [float('inf')] * N #simply append infs all the way
        else:
            col = []
            for deletion in xrange(N+1):
                # second iterator equivalent to deletion
                segment = []
                for deletion2 in xrange(N+1):
                    # find each amplification value
                    ampDel = getAmp(i, deletion, Src, Tgt, difs)
                    ampDel2 = getAmp(prevIndices[i], deletion2, Src, Tgt, difs)
                    # if any amp value is infinite, solution is infinite length
                    if ampDel == float('inf') or ampDel2 == float('inf'):
                        segmentVal = float('inf')
                    else:
                        #see Recursion formula from Shamir et al
                        segmentVal = ( table[prevIndices[i]][deletion2] + 
                                       max(deletion - deletion2, 0) + 
                                       max(ampDel - ampDel2, 0) + 
                                       max(qArray[i] - 
                                           max(deletion, deletion2), 0) )
                    segment.append(segmentVal)
                # the value of M[i,d] at this specific point
                M_i_d_val = min(segment)
                col.append(M_i_d_val)
        table.append(col)
    return table


def mTableInit(Src, Tgt, difs, n, N, prevIndices, qArray):

    '''     takes profiles, differences, lengths, prev(i) and Qi     '''
    ''' returns first column of table M containing optimal solutions '''

    mTable = []
    # initial table computation
    for deletion in xrange(N+1): #include N
        amp = getAmp(0, deletion, Src, Tgt, difs)
        val = deletion + amp
        mTable.append(val)
    return mTable


def getFinalValue(mTable):

    ''' takes in solution table M, returns final minimum length solution '''

    segment = mTable[-1]
    return min(segment)


def getPrevIndices(Tgt):

    ''' takes in Tgt profile, returns prev(i) values '''
    ''' prev(i) for any i > 0 is the largest index less '''
    ''' than i such that Tgt[that index] > 0. '''

    prevIndices = []
    # first elem will be none, since prev(i) is offset by 1 elem from profiles
    prevIndices.append(None) 
    for i in xrange(1, len(Tgt)):
        if Tgt[i-1] == 0: # adopt the previous prevIndex
            # edge-case: i = 1
            if i == 1:
                valToAdd = 0 # [prev(1) set to 0 for now]
            else:
                valToAdd = prevIndices[i-1]
        else: # make i-1 the prevIndex
            valToAdd = i-1
        prevIndices.append(valToAdd)
    return prevIndices


def getQArray(prevIndices, Src):
    
    '''      takes in prevIndices, Src profile      '''
    '''  returns qArray, where qArray[i] signifies  '''
    '''  the number of deletions to perform over a  '''
    '''  contiguous segment of zeros if one arises. '''

    qArray = []
    qArray.append(None) # see similar comment in getPrevIndices()
    for i in xrange(1, len(Src)):
        if prevIndices[i] == i-1:
            # no contiguous zero segment here; perform 0 deletions
            valToAdd = 0
        else:
            #use the max copy number in the src between prev(i) and i
            segment = Src[prevIndices[i]+1 : i]
            valToAdd = max(segment)
        qArray.append(valToAdd)
    return qArray


def getAmp(i, deletion, Src, Tgt, difs):

    ''' takes in index i, deletion, profiles / their difs '''
    ''' returns # of amplifications needed to offset current deletion '''
    
    if Tgt[i] > 0 and max(difs[i], 0) <= deletion < Src[i]: #assumption 12
        ampVal = -difs[i] + deletion
    else: 
        ampVal = float('inf')
    return ampVal


###############################################################################
#############         SAME ALGORITHM, REDUCED COMPLEXITY          #############
###############################################################################


'''    Many array variables ex. qArray, prevIndices can be expressed as     '''
'''                 recursive values, with base values at 0                 '''


def memoizer(fn):

    ''' takes and modifies input function (fn) to memoize repeating answers '''

    answers = {}
    
    def helper(*args): # check for previously calculated answers
        argsList = list(args)
        for (i, val) in enumerate(argsList):
            if isinstance(val, list): # make it a tuple to then hash it
                argsList[i] = tuple(val)

        argKey = tuple(argsList) # use argsList for hashing (args has mutables)
        if argKey not in answers:
            answers[argKey] = fn(*args) # add the answer in
        return answers[argKey]
    
    return helper


def prev(i, Tgt):

    ''' Takes in index, target profile and returns maximum index less than  '''
    '''       i at which the target profile at that index is nonzero.       '''

    if i == 0 or i == 1:
        val = 0
    elif Tgt[i - 1] == 0:
        val = prev(i - 1, Tgt)
    else:
        val = i - 1
    return val


def Q(i, j, Src, Tgt):

    if i == j == 0: # edge-case for prev(i) of 0
        val = 0
    elif j == i - 1: 
        val = 0
    else:
        segment = Src[j + 1: i]
        val = max(segment)
    return val


def f(d_i_min, d_i_max, a_i, b_i):

    ''' helper 'piecewise' function for the O(nN) piecewise algorithm below '''
    ''' takes in minimum and maximum allowed deletions and placeholder vars '''
    ''' returns current distance value at index i according to these values '''

    if d_i_min <= a_i:
        return 0
    elif a_i < d_i_min <= b_i:
        return d_i_min - a_i
    elif b_i < d_i_min <= d_i_max:
        return 2*d_i_min - a_i - b_i 
    else: return 0


def linearDistance(Src, Tgt):

    ''' Takes in src, tgt profiles. Returns CNT distance from src -> tgt or ''' 
    '''   the number of edit operations necessary to transform src -> tgt.  '''
    '''  To make sense of each of the variables, see Shamir et al's paper.  '''

    assert(len(Src) == len(Tgt)) # otherwise we cannot compare!

    # initialize a dynamic programming vector as zero-dictionaries
    placeholder_a, placeholder_b, bases = {0:0}, {0:0}, {0:0}
    N = max( max(Src), max(Tgt) ) # max copy number from src / tgt
    # to justify initial zeroes in the dicts (doesn't affect total dist):
    prelimEdits(Src, Tgt, N) # destructively modify src and tgt

    for i in range(len(Src)):
        if Tgt[i] == 0: continue

        j = prev(i, Tgt) # maximum j < i at which tgt[j] != 0
        Q_i = Q(i, j, Src, Tgt) # deletions to perform from prev(i) to i
        diff_j, diff_i = Src[j] - Tgt[j], Src[i] - Tgt[i]
        R_i = diff_j - diff_i
        d_i_min, d_i_max = max(diff_i, 0), max(Src[i] - 1, 0)

        if R_i >= 0:
            if Q_i <= placeholder_a[j]:
                placeholder_a[i] = placeholder_a[j] - R_i
                placeholder_b[i] = placeholder_b[j]
            elif placeholder_a[j] < Q_i <= placeholder_b[j]:
                placeholder_a[i] = placeholder_a[j] - R_i
                placeholder_b[i] = placeholder_b[j]
            elif placeholder_b[j] < Q_i:
                placeholder_a[i] = placeholder_b[j] - R_i
                placeholder_b[i] = Q_i

        elif R_i < 0:
            if Q_i <= placeholder_a[j]:
                placeholder_a[i] = placeholder_a[j]
                placeholder_b[i] = placeholder_b[j] - R_i
            elif placeholder_a[j] < Q_i <= placeholder_b[j]:
                placeholder_a[i] = Q_i
                placeholder_b[i] = placeholder_b[j] - R_i
            elif placeholder_b[j] < Q_i:
                placeholder_a[i] = min(Q_i, placeholder_b[j] - R_i)
                placeholder_b[i] = max(Q_i, placeholder_b[j] - R_i)

        # now calculate variables that are returned
        bases[i] = ( bases[j] + max(Q_i - placeholder_a[j], 0) + 
                     f(d_i_min, d_i_max, placeholder_a[i], placeholder_b[i]) )

        # limit the ranges of delimiters according to d_i_min, d_i_max
        placeholder_a[i] = max(d_i_min, min(placeholder_a[i], d_i_max))
        placeholder_b[i] = max(placeholder_a[i], 
                               min(placeholder_b[i], d_i_max))  
    # return cumulative distance; reverse initial modifications
    finalIndex = max(bases.keys())
    undoPrelimEdits(Src, Tgt, N)
    return bases[finalIndex]


def prelimEdits(Src, Tgt, N):

    ''' Add in placeholder copy numbers so that recursive induction starts '''
    ''' and ends at unimportant values. Does not affect the final distance '''

    Src.insert(0, N+1)
    Tgt.insert(0, N+1)
    Src.insert(len(Src), N+1)
    Tgt.insert(len(Tgt), N+1)


def undoPrelimEdits(Src, Tgt, N):

    ''' An undo function for prelimEdits. We don't want CN profiles to be '''
    ''' twice as big as before distance analysis, so we need this fxn. '''

    Src.pop(0)
    Tgt.pop(0)
    Src.pop()
    Tgt.pop()


###############################################################################
############################### TEST FUNCTIONS! ###############################
###############################################################################


def testDistances():

    ''' verifies that our polynomial time function will generate the same '''
    '''  distance result as its linear time counterpart, on all inputs.   '''

    src = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    tgt = [1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1]
    assert(polytimeDistance(src, tgt) == 2)
    assert(linearDistance(src, tgt) == 2)

    src, tgt = [1, 1, 1], [2, 0, 2]
    assert(polytimeDistance(src, tgt) == 2)
    assert(linearDistance(src, tgt) == 2)

    src = [3, 1, 2, 3, 2, 1, 4]
    tgt = [2, 0, 0, 0, 0, 0, 2]
    assert(polytimeDistance(src, tgt) == 3)
    assert(linearDistance(src, tgt) == 3)

    print 'all distance tests passed!'


###############################################################################
#################################### MAIN! ####################################
###############################################################################


def main():
    print 'running main... \n',

    file = getCopyNumberDirec()
    print 'getting values... '
    profiles, cells, regions, names = parseDeanTiming(file)
    print 'done.'

    response = raw_input('''Type one of the following and press Enter:
                \'trial\' to get a tree from copy number data.
                See the README for more info on using these globals.
                \'exit\' to exit the prompt entirely.\n''')

    try:
        if response in ['trial', 'sub']:
            runTrialPrompts(profiles, cells, regions, names)
        elif response == 'exit':
            return # exit prompt
        else:
            print 'invalid response.'
            return
    
    except:
        getErrorMessage()
        return


def getErrorMessage():

    print '\n'
    print '#' * 70
    print 'ERROR OCCURRED !!!'
    print '#' * 70
    print ' '
    traceback.print_exc(file=sys.stdout) 
    print ' '
    print '#' * 70


###############################################################################
###                            EVERYTHING THAT RUNS                         ###
###############################################################################


if __name__ == '__main__':
    main()

#EOF