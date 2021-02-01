===============================================================================
                    SINGLE-CELL SEQUENCING PHYLOGENY PROJECT                    
===============================================================================
Last Updated 07/26/18
===============================================================================
Tae Yoon (Tyler) Park (Implemented Neigbor-Joining Algorithms, 2016 - 2017)
===============================================================================
Suraj Joshi (Implemented New Methods for Phylogeny Construction, 2017 - 2018)
===============================================================================
Russell Schwartz (Principal Investigator)
===============================================================================

################################################
###                                          ###  
###      This file contains:                 ###
###                                          ###
###   1. OVERVIEW                            ###
###   2. NECESSARY LANGUAGES / LIBRARIES     ###
###   3. FILES (GENERAL DESCRIPTIONS)        ###
###   4. HOW TO RUN (SPECIFIC INSTRUCTIONS)  ###
###   5. FUTURE DIRECTIONS                   ###
###   6. KNOWN ISSUES                        ###
###   7. RESULTS                             ###
###                                          ###
################################################

===============================================================================
 1. OVERVIEW : CONCEPTION, EXPERIMENTAL DATA, METHODS, GOALS
===============================================================================

This project uses copy number data obtained from single-cell sequencing of
multiple tumor sites from two patients to generate phylogenies of each of
the tumor samples' sequences.

The sequences we are comparing are lists of copy numbers, or 
"copy number profiles".

We use various measures of pairwise evolutionary distance for every possible
pair of profiles to cluster our single cell data into complete phylogenies. 

Goals:

1) Implement a robust program for clustering any input list of copy number
   profiles.
2) Use a realistic distance measure that mimics mutation / genetic variation
   among tumor cells during cancer progression.
3) Implement an algorithm that accomplishes these goals with appropriate
   efficiency (timewise / memorywise).

We use three measures of evolutionary distance: 

1) L1 (Minkowski) Distance
2) L2 (Euclidean) Distance
3) CNT (Copy Number Transformation) Distance

We want to compare the accuracies of these methods in phylogeny generation, 
with hopes of determining a method statistically superior to the rest.

The third method is a direct implementation of an algorithm proposed by
Ron Shamir, Meira Zehavi, and Ron Zeira from Tel Aviv University in 2016.
So far, this method shows a much better accuracy in validation tests compared
to both L1 and L2.

Complexity:

The L1 and L2 scripts follow O(n^3) neighbor joining, where n is the number 
of profiles. Each distance calculation is O(m), where m is the length of each 
profile.

The CNTP script also follows O(n^3), but my distance calculation is on the order
of O(mN), where N is the maximum copy number among all profiles.


===============================================================================
 2. NECESSARY LANGUAGES / LIBRARIES
===============================================================================

Python 2 is required for running this code. We use Python 2.7.12 specifically.

To download Python 2, visit this website and download the appropriate version.

https://www.python.org/downloads/

===============================================================================

Some additional non-standard Python 2 libraries are needed to run this code:

1) ETE Toolkit (ete3), a data visualization library that allows phylogenetic 
trees to be visualized and analyzed for topology and accuracy.

   There are two ways to install ete3:

   A. Visit this website:
      http://etetoolkit.org/
      And install the appropriate version (ete3 is the specific one we use).

   B. Open the command line prompt and type in:
      python2 -m pip install ete3
      This method assumes that you have already added Python 2 to your default
      PATH.

2) Numpy, a scientific computing library that allows manipulation of matrices
   and statistical variables.

   There are two ways to install numpy:

   A. Visit this website and follow the directions for adding the Numpy binary
      to Python 2:
      https://docs.scipy.org/doc/numpy-1.10.1/user/install.html

   B. Open the command line prompt and type in:
      python2 -m pip install numpy

3) MatPlotLib, a library for plotting and analyzing graphs, also allowing
   manipulation of statistical variables.

   Visit this page for a comprehensive guide on how to install matplotlib for
   all operating systems:

   https://matplotlib.org/users/installing.html

4) Scipy, another scientific computing software that also allows manipulation
   of statistical variables.

   There are two ways to install scipy:

   A. Visit this website and follow the directions for installing scipy for
      your OS:

      https://www.scipy.org/install.html

   B. Open the command line prompt and type in:

   python2 -m pip install scipy

   NOTE: Some people have issues installing scipy if both numpy and mkl are
   not installed (a package called numpy+mkl). If this error occurs, please
   uninstall numpy using "python2 -m pip uninstall numpy" and visit this page
   to install both packages simultaneously before attempting to install Scipy
   again:

   https://pypi.python.org/pypi/numpy-mkl

===============================================================================

Some data visualization programs are also required to run this code.

1) Graphviz

    Our code outputs .dot files containing phylogenies. These must be converted 
    to .png files for viewing. Graphviz is able to perform the conversion.

    How to install: 
        1) http://www.graphviz.org/Download..php
           Follow the directions on this page.

        If you are on Windows:
        2) Add Graphviz to your default PATH, so that you can 
           run the .dot conversion from the command line prompt.

        If not on Windows:
        3) You may have to do the conversion manually from Graphviz. 
	   Open Graphviz, navigate to the outputted .dot file (usually located 
           in scSeqPhylogeny/clusteringFiles), and convert it to the desired 
           file format (.png is recommended).


===============================================================================
 3. FILES : GENERAL DESCRIPTIONS
===============================================================================

The top-level directory is scSeqPhylogeny. Files within this directory include:

1) _other

This contains files unnecessary to the functioning of the program:

2) clusteringFiles

clusteringFiles includes:

    A. cntp_final.py - The main file that runs on all copy number profiles and
       constructs a complete phylogeny of each of these profiles.

    B. final_nj_Suraj.py - The main neighbor-joining script.

       This file runs the neighbor-joining algorithm on any given distance 
       matrix containing pairwise distances for all relevant profile pairs.
       The distance matrix is calculated based on the input distance
       measure (L1 or L2) when the file itself is run, and a .dot file of the
       tree is generated within the same directory. The program will attempt
       to convert the dotFile into a .png file of the phylogenetic tree, but
       if this fails, the user is obliged to use the following command to
       complete the conversion manually:

       dot -Tpng <<filename>>.dot -o <<outputname>>.png

    C. trees - A directory containing all .dot and .png files generated from
       the algorithms in clusteringFiles.

    D. shamir_CNTP.pdf - The paper from which the CNTP method was derived and
       implemented.

    E. older versions - contains older versions cntp_final and final_nj.

2) deanFiles

deanFiles includes allcp, a directory of all of the copy number profiles, .dot 
files of older clustering trials from the entire data set and of sections 
of the data set separated by patient / region sampled.

4) outputs

outputs contains all outputs regarding the phylogenies created:

    A. analyses, a directory containing the outputs from the validation tests
       conducted over all of our trees.
    B. finalTrees, a directory containing the .dot and .png versions of all of
       our created trees.
    C. scSeqPhylogeny-REPORT, a directory containing the most recent
       all project proposals and results.

5) SCSEQ

This is a directory containing sample single cell data from SRA.

6) tree_cluster_forall_final.py

tree_cluster_forall_final.py performs validation tests for all samples in a
given input tree.

7) dot_to_nw.py

This file converts .dot files to Newick formatted graphs.

8) undir_to_dir.py

This file converts unrooted trees to rooted trees.


===============================================================================
 4. HOW TO RUN
===============================================================================

    A. To run L1 / L2:

        i. Open Python 2 in an IDE.
        ii. Navigate to final_nj_Suraj.py, open, and run it.
        iii. Follow the program's prompts. The tree takes 2-5 minutes to
             generate.
        iv. When naming your .dot file, please do not include numerical / 
	    arithmetic symbols or numbers in the name.

    B. To run CNT:

        i. Open Python 2 in an IDE.
        ii. Navigate to cntp_final.py, open, and run it.
        iii. Follow the program's prompts. The tree takes about 20-30 minutes
	     to generate. 
             Everytime a row of pairwise distances is calculated to completion, 
             an index will be printed. This means the program is still running.

    C. To validate trees:

        i. Open Python 2 in an IDE.
        ii. Navigate to tree_cluster_forall_final.py, open, and run it.
        iii. If you are validating one of our phylogenies, follow the program's
             directions.
        iv. If you are validating a phylogeny you generated, type in the 
	    relative path (from the top level project directory) to your 
	    phylogeny, and press Enter.
            Use '\\' as your separater.
            Ex. scSeqPhylogeny\\rest_of_my_path\\my_tree.dot
        v. Input the number of desired trials.

    D. To access copy number variables:

        i. Open Python 2 in an IDE.
        ii. Navigate to cntp_final.py, open, and run it.
        iii. Press enter immediately upon the program's initial prompt. Wait
        about 4-6 seconds, and all of the copy number data will be accesible
        from the shell / IDE.
        iv. The variable names you can access:
            
            profiles: all copy number profiles collected into a list
            
            cells: the number of profiles / samples.
        
            regions: the number of copy numbers per profile. All profiles are
            of the same length.

            names: the names of each profile's sample, which follows a specific
            string encoding. If you would like to see the specific nature of
            this encoding, Navigate to deanFiles/allcp from the top-level
            directory.


===============================================================================
 5. FUTURE DIRECTIONS
===============================================================================

As of now, the motivation is to investigate the inefficiencies in the CNT
algorithm, find a fourth method that we can implement to a similar or improved 
level of accuracy, improve the efficiency of the methods we are currently 
using, or create our own method by combining some or all of the methods we 
have used thus far.

The general idea is to allow for more possibilities of mutations over a given
CNP than simply segmental deletions and amplifications.


===============================================================================
 6. KNOWN ISSUES
===============================================================================

1) While Shamir, Zehavi, and Zeira's method is noted to be computable in linear
time on the length of a profile for any pair of profiles, our implementation
of it takes considerably more time to generate a complete phylogeny than the
standard L1 / L2 clustering methods. 


===============================================================================
 7. RESULTS
===============================================================================

Please navigate to scSeqPhylogeny/outputs.

===============================================================================