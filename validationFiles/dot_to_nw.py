###############################################################################

from ete3 import Tree

def convert(s):

    ''' Takes in .dot filename as a string, reads through file and converts '''
    '''               into Newick, which is readable by ete3.               '''

    nwlist = []
    f = open(s)
    fl = f.readlines()

    for i in fl:
        i=i.replace(']','')
        i=i.replace(';','')

        if ('->' in i.split()): #for GT with no branch lengths
            nwlist.append((i.split()[0],i.split()[2])) 

        elif ('--' in i.split()):
            nwlist.append((i.split()[0],i.split()[2]))

    tree = Tree.from_parent_child_table(nwlist)
    return tree

#EOF