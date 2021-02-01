# Suraj Joshi (2017)
# 10/14/2017

###############################################################################
####                        MAIN FOR ALL OTHER MAINS                       ####
###############################################################################


# import dependencies
import os

from validationFiles.tree_cluster_forall_final import validation

try: from subprocess import check_call
except: print 'warning: failed to import check_call'

def main():
    # runs every file needed for analysis
    print "Would you like to do analysis or validation? Type one, press Enter."
    choice = raw_input("") # no prompt needed
    topDirec = os.path.abspath(os.path.dirname(__file__))

    if choice == "analysis":
        clusterPath = (topDirec + '/clusteringFiles/cntp_final.py')
        try: check_call(['python2','-i', clusterPath])
        except: raise Exception("Program crashed or closed.")

    elif choice == "validation":
        setupValidation() # run this setup function before tree_cluster_forall

    else: print "Invalid choice! Re-call main."


def setupValidation():

    topDirec = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    # access pre-generated phylogenies
    deanFile_L1 = topDirec + '/outputs/finalTrees/L1_all_samples.dot'
    deanFile_L2 = topDirec + '/outputs/finalTrees/L2_all_samples.dot'
    deanFile_CNTP = topDirec + '/outputs/finalTrees/CNTP_all_samples.dot'

    testFile = raw_input('''Which file to test? 
            Type in one of the following and press Enter,
            \'L1\', \'L2\', or \'CNTP\',
            to validate the full-sample phylogenies
            that we have already generated, or
            \'NEW\' to validate a newly generated phylogeny. You will be
            asked for the phylogeny\'s .dot file\'s relative path. \n''')

    if testFile == 'L1': testFile = deanFile_L1
    elif testFile == 'L2': testFile = deanFile_L2
    elif testFile == 'CNTP': testFile = deanFile_CNTP
    elif testFile == 'NEW':

        print 'example path: scSeqPhylogeny/myPath/myTree.dot'
        fileStr = raw_input('''            Type in your tree\'s relative path 
            according to the example path above, then press Enter. 
            Use single-slash (/) separaters! \n''')
        testFile = fileStr.replace('scSeqPhylogeny/', '')
        validation(testFile)

    else: raise ValueError('invalid choice: pick L1, L2, CNTP, or NEW!')


if __name__ == '__main__':
    main()
    exit()

#EOF