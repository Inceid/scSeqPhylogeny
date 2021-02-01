###############################################################################

def convert(s):

    ''' Takes in .dot filename as string, converts undirected tree into   '''
    '''  directed tree and returns the instructions for drawing the tree. '''

    dglist = []
    nlist = []
    f = open(s)
    fl = f.readlines()
    for i in fl:
        i=i.replace('"','')
        i=i.replace('|','_')
        if '--' in i.split():
            isP = True
            for j in range(0,len(i.split()[0]),2):
                if i.split()[0][j] not in i.split()[2]:
                    isP = False
            dglist.append('\t')
            dglist.append(i.split()[0].replace(';',''))
            dglist.append(' [shape = box];\n')
            dglist.append('\t')
            dglist.append(i.split()[2].replace(';',''))
            dglist.append(' [shape = box];\n')
            if not isP:
                dglist.append('\t')
                dglist.append(i.split()[0].replace(';',''))
                dglist.append(' -> ')
                dglist.append(i.split()[2])
                dglist.append('\n')
            else:
                dglist.append('\t')
                dglist.append(i.split()[2].replace(';',''))
                dglist.append(' -> ')
                dglist.append(i.split()[0])
                dglist.append('\n')

        elif 'strict' in i:
            i = i.replace('strict graph','digraph g1')
            dglist.append(i)
        else:
            dglist.append(i)
                        
    return dglist

#EOF