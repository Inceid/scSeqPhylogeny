# Tae Yoon (Tyler) Park (2016-2017)
# Suraj Joshi (2017)
# 10/14/2017

###############################################################################
####  IMPLEMENTATION OF APPLYING NEIGHBOR-JOINING ALGORITHMS TO FISH DATA  ####
###############################################################################

'''
LOG: Tyler

 1-7 : 02/19/16 fix x & y, joined identical rows
 1-8 : 04/11/16 add dummy node, add legend
 1-9 : added parse function for csv files
 1-9_2 : modified node names to real names from imported data file

LOG: Suraj

 1-9_2_Suraj : repaired tree generation
               re-labeled internal nodes for legibility
 1-9_3_Suraj : documented code, cleaned up comments
 1-9_4_Suraj : deleted all defunct functions.
 final_nj_Suraj : edited some comments, cleaned up code, "robustified".
 final_nj : modified labelling, took out hardcoding of highest merge step

'''

# first import our dependencies

import copy
import os
import profile
import random
import time

# to use cmd's dot-to-png conversion command from python

try: from subprocess import check_call
except: print 'warning: failed to import check_call'


###############################################################################
####    HELPER FUNCTIONS FOR PARSING AND CALCULATING PAIRWISE DISTANCES    ####
###############################################################################


def getCopyNumberDirec():

    ''' Finds, prints, and returns the path containing our raw CN data '''

    topDirec = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    cpDirec = topDirec + '\\deanFiles\\allcp\\'
    print 'location of all copy number files: ', cpDirec

    return cpDirec


def parseDean(direc):

    ''' Reads all cell data from deanFiles directory 'direc', '''
    ''' returns a list of copy number profiles, the number of cells, '''
    ''' the number of regions sampled, and a list of cell names as strings. '''
    ''' NOTE: cells = number of profiles ; regions = length of one profile. '''

    profiles, cells, names = [], 0, []
    
    for filename in os.listdir(direc):

        start = filename.index('SC')
        end = filename.index('.')
        cellname = 'G'+filename[3:5]+'SP'+filename[8:9]+filename[start:end]

        names.append(cellname)
        f = open(direc + filename)
        fl = f.readlines()
        profile = []
        copynumbers = fl[1:]
        regions = 0
        for i in range(len(copynumbers)):
            profile.append(int(copynumbers[i].replace('\n','')))
            regions += 1
        profiles.append(profile)
        f.close()
        cells += 1

    return (profiles, cells, regions, names)


def distanceMatrixL1(profiles, cells, regions, names):

    ''' Takes list of profiles (l), number of profiles (cells), len of each '''
    '''     profile (regions), list of names of samples (names).    ''' 
    '''  Calculates and returns a distance matrix (dM), a dict that stores  ''' 
    '''   distances between each src / tgt profile, and two copies of the   '''
    '''                       of taxa to be clustered.                      '''

    matrix, taxaList = dict(), []

    for i in range(cells):
        taxaList.append(str(names[i]))
        for j in range(i + 1, cells):
            dist = 0
            for k in range(regions):
                dist += abs(profiles[i][k] - profiles[j][k])
            matrix[str(names[i]), str(names[j])] = dist

    return (matrix, taxaList, copy.copy(taxaList))


def distanceMatrixL2(profiles, cells, regions, names):

    '''   takes same inputs as distanceMatrixL1, returns a distance matrix  '''
    '''       similar to dM1, except that L2 distances are calculated.      '''

    matrix, taxaList = dict(), []
    for i in range(cells):
        taxaList.append(str(names[i]))
        for j in range(i + 1, cells):
            dist = 0
            for k in range(regions):
                dist += abs(profiles[i][k] - profiles[j][k])**2
            # take square root of the sum of squares, add into matrix
            matrix[str(names[i]), str(names[j])] = dist ** 0.5

    return (matrix, taxaList, copy.copy(taxaList))


def calcRowSums(dM):

    ''' takes input distance matrix, finds and returns row sum '''
    '''               input row should be in string            '''
    
    rowSumDict = dict()
    for (x,y) in dM:

        try: rowSumDict[x] += dM[x,y]
        except: rowSumDict[x] = dM[x,y]
        
        try: rowSumDict[y] += dM[x,y]
        except: rowSumDict[y] = dM[x,y]

    return rowSumDict


def calcQm(dM, rsD, taxa):

    ''' calculates Q matrix needed for neighborJoin from dist matrix(dM) '''
    '''         dict(distmat)*int(row size) -> dict(Qmat) '''

    qM = dM.copy()
    for (i, j) in dM:
        qM[i,j] = (taxa-2) * dM[i,j] - (rsD[i]) - (rsD[j])

    return qM


def initTree(dM, taxaL):

    '''    initializes the completely unresolved tree    '''
    ''' the center node in unresolved tree is labeled -1 '''

    seen = ["-1"]
    t = ["\t-1 [shape = box];\n"] # initial dummy node
    for i in taxaL: # i represents the taxum string here

        if i not in seen: # add nodes
            color = selectColor(i)
            t.append("\t%s [shape=box, style=filled, color=%s];\n" % (i, color))
            t.append("\t-1 -> %s;\n" % i)
            seen.append(i)

    return t


def selectColor(taxaStr):
    if 'G07SP1' in taxaStr: color = 'crimson'
    elif 'G07SP2' in taxaStr: color = 'orange'
    elif 'G07SP3' in taxaStr: color = 'gold'
    elif 'G33SP1' in taxaStr: color = 'forestgreen'
    elif 'G33SP2' in taxaStr: color = 'lightskyblue'
    else:
        assert('G33SP3' in taxaStr)
        color = 'darkviolet'
    return color


def updateDm(dM, u, mi, mj, rsD, taxa, taxaL):

    ''' destructively updates the input dM after merging two nodes '''
    '''  calculates u to i and u to j as well as u to other nodes  '''
    '''   where u is i and j merged, and (i,j) are the min pair    '''

    rowIsum = rsD[mi]
    rowJsum = rsD[mj]
    d = dM.copy()

    u_mi = (0.5) * (dM[mi,mj]) + (1.0 / (2 * (taxa-2))) * (rowIsum-rowJsum)
    u_mj = dM[mi,mj] - u_mi

    for x in taxaL:
        if (x != mi and x != mj):
            mimj = dM[mi,mj]
            if (mi,x) in dM:
                mix = dM[mi,x]
            else:
                mix = dM[x,mi]
            if (mj,x) in dM:
                mjx = dM[mj,x]
            else:
                mjx = dM[x,mj]
            
            dM[u,x] = (0.5) * ( mix + mjx - mimj )

    #pop out members of merged node
    for (i,j) in d:
        if (i == mi or i == mj or j == mi or j == mj):
            del dM[i,j]

    taxaL.append(u)
    taxaL.remove(mi)
    taxaL.remove(mj)

    return (dM, u_mi, u_mj,taxaL)


def updateGraph(t, u, mi, mj, u_mi, u_mj):

    '''           tree updating function          '''
    ''' manually adds in the next data-containing '''
    '''   node given the appropriate string data  '''
    ''' u,mi,mj: string, u_mi,u_mj: int(distance) '''
    '''  returns t, list of instructions for tree '''

    t.append("\t%s [shape = box];\n" % u)  #add the merged node OLD LINE

    t.remove("\t-1 -> %s;\n" % mi)
    t.remove("\t-1 -> %s;\n" % mj)
    t.append("\t%s -> %s [label = %f];\n" % (u, mi, u_mi))
    t.append("\t%s -> %s [label = %f];\n" % (u, mj, u_mj))
    t.append("\t-1 -> %s;\n" % u)

    return t
    

def isParent(label):

    ''' Takes in a label, returns whether label denotes a parent node. '''
    return ('_' in label)


def labelRoot(t, label1, label2, dist):

    ''' Takes in tree and data from final distance matrix tuple, labels '''
    '''  the tree with the root node and corresponding final distance.  '''

    t.append("\tROOT [shape = box];\n")
    if isParent(label1):
        # then label distance on label2
        t.append("\tROOT -> %s ;\n" % label1)
        t.append("\tROOT -> %s [label = %f];\n" % (label2, dist))
    else: 
        assert(isParent(label2)) # this is the only other possibility
        t.append("\tROOT -> %s [label = %f];\n" % (label1, dist))
        t.append("\tROOT -> %s;\n" % label2)


def neighborJoin(matrix, taxa, taxaL, taxaL2, names):

    ''' Takes in distance matrix (matrix), number of taxa (taxa), '''
    ''' two lists of all samples to cluster (taxaL, taxaL2), list of sample '''
    ''' names, and a bool that says whether the tree's internal nodes are '''
    ''' labelled. returns instructions (t) for drawing the tree.'''

    # code for debugging: 
    print "our matrix is using... \n"
    for taxum in taxaL: print taxum

    t = initTree(matrix, taxaL)
    dM = matrix.copy() # non-destructive

    parentNum = 1 # initialize label for each parent node; apply when needed

    while (taxa > 2):

        rsD = calcRowSums(dM)
        qM = calcQm(dM, rsD, taxa)
        (mi, mj) = min(qM, key=qM.get)

        useParentLabel = False

        if isParent(mi): 
            useParentLabel = True
            new_mi = 'Parent%d' % parentNum
            parentNum += 1
        else: 
            new_mi = mi
        
        if  isParent(mj): 
            if not useParentLabel:
                parentNum += 1
            useParentLabel = True
            new_mj = 'Parent%d' % parentNum
        else:
            new_mj = mj

        if useParentLabel:
            u = ( 'Parent%d' + '_' + new_mi + '_' + new_mj ) % parentNum
        else:
            u = mi + '_' + mj

        (dM, u_mi, u_mj, taxaL) = updateDm(dM, u, mi, mj, rsD, taxa, taxaL)
        taxa -= 1
        t = updateGraph(t, u, mi, mj, u_mi, u_mj)

    if (taxa == 2):
        ( (final_mi, final_mj), dist ) = dM.popitem()
        t.remove("\t-1 -> %s;\n" % final_mi)
        t.remove("\t-1 -> %s;\n" % final_mj)
        t.remove("\t-1 [shape = box];\n")
        # add default root node in
        labelRoot(t, final_mi, final_mj, dist)

    return t


def writeDotFile(t, LD):

    ''' Takes in instructions for writing tree (t), L-distance being used, '''
    '''   and outputs a .dot file of the tree. Asks for output filename.   '''
    '''            The output .dot file can be read with GVEdit            '''

    dataname = raw_input("output filename? ")
    fileName = "%s_%s.dot" % (LD, dataname) # create new file

    f = open(fileName,'w')
    # read and concatenate entries of instructions
    text = "digraph tree_of_%s {\n" % dataname + "".join(t) + "}"
    f.write(text)
    f.close()

    return fileName


def runL1L2(l, cells, probes, names, LD):

    startTime = time.time()

    if LD == 'L1': # calculate L1 Distance tree
        (dM, taxaL, taxaL2) = distanceMatrixL1(l,cells,probes,names)

    elif LD == 'L2': # calculate L2 Distance tree 
        (dM, taxaL, taxaL2) = distanceMatrixL2(l,cells,probes,names) 
        
    # get tree-drawing instructions as string t
    t = neighborJoin(dM, cells, taxaL, taxaL2, names)
    fileName = writeDotFile(t, LD)
    fileName = fileName.replace('.dot', '')
    
    try: # now try to convert to png within program
        print 'command line is converting your .dot file... '
        check_call(['dot','-Tpng', fileName + '.dot','-o', fileName + '.png'])
    except: 'command-line call failed. Convert the .dot file manually.'

    endTime = time.time()
    print 'done. time elapsed: %d seconds' % (endTime - startTime)
    print 'Your tree (.dot and/or .png) should be in this file\'s directory.'


###############################################################################
####################################  MAIN ####################################
###############################################################################


def main():

    print 'running main...'
    cpDirec = getCopyNumberDirec()
    (l, cells, probes, names) = parseDean(cpDirec)
    runL1L2(l, cells, probes, names)


###############################################################################
############################ EVERYTHING THAT RUNS! ############################
###############################################################################


if __name__ == '__main__':
    main()
    exit()

#EOF