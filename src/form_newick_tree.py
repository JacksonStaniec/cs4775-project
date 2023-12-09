#!/usr/bin/env python3

''' Outputs the Newick formatted tree after performing the neighbor-joining
    algorithm on an arbitrary number of species.

Arguments:
    -f: distances file (symmetric matrix with 0 on the diagonal)
        (default is dist10.txt)
Outputs:
    Newick formatted tree after neighbor-joining

Example usage:
    python 1a.py -f dist10.txt
'''

import argparse
import numpy as np
import pandas as pd

''' Reads the input file of distances between the sequences

Arguments:
    distances_file: file name of distances between sequences
Returns:
    D: matrix of distances (map of maps)
    mapping: index to name mapping (dictionary)
'''
def read_data(distances_file):
    with open(distances_file, "r") as f:
        lines = [l.strip().split() for l in f.readlines()]
        mapping = {i: s for i, s in enumerate(lines[0])}
        lines = [l[1:] for l in lines[1:]]
        D = {i: {} for i in range(len(lines))}
        for i, l in enumerate(lines):
            for j, sval in enumerate(l):
                D[i][j] = float(sval)
    return D, mapping


''' Performs the neighbor joining algorithm on a given set of sequences.

Arguments:
    D: map of maps, defining distances between the sequences
       (initially n x n, symmetric, 0 on the diagonal)
       (index -> index -> distance)
Returns:
    To be determined
'''

def neighbor_join(D,mapping):
    #the tree joining algo kinda ended up being rooted already
    d=pd.DataFrame(D).to_numpy()
    n = np.shape(d)[0]
    idxs = pd.DataFrame(np.array(list(mapping.items())))[1].to_numpy()
    while n>=3:
        Q = np.zeros((n,n))
        rowsums = np.sum(d,1)
        #compute Q
        #Q(i,j)=(n-2)d(i,j)-sum(d(i,k))-sum(d(j,k))
        for j in range(n):
            for i in range(j):
                Q[i][j]=(n-2)*d[i][j]-rowsums[i]-rowsums[j]
        mincoord=np.unravel_index(np.argmin(Q),np.shape(Q))
        #we have to merge mincoord[0] and mincoord[1] into one node
        #branch length
        d1=0.5*d[mincoord[0]][mincoord[1]]+(1/(2*(n-2)))*(rowsums[mincoord[0]]-rowsums[mincoord[1]])
        d2=d[mincoord[0]][mincoord[1]]-d1
        mc1=idxs[mincoord[0]]
        mc2=idxs[mincoord[1]]
        idxs=np.delete(idxs,mincoord)
        idxs = np.append('({}:{:.6f},{}:{:.6f})'.format(mc1,d1,mc2,d2),idxs)
        #update D
        n-=1
        dNew=np.zeros((n,n))
        dNew[1:,1:]=np.delete(np.delete(d,mincoord,axis=1),mincoord,axis=0)
        dNew[1:,0]=np.delete(0.5*(d[mincoord[0]]+d[mincoord[1]]-d[mincoord[0],mincoord[1]]),mincoord)
        dNew[0,1:]=np.delete(0.5*(d[mincoord[0]]+d[mincoord[1]]-d[mincoord[0],mincoord[1]]),mincoord)
        d=dNew
    returned ='({}:{},{}:{});'.format(idxs[0],0.5*d[1,0],idxs[1],0.5*d[1,0])
    return returned


''' Helper function for defining a tree data structure.
    First finds the edge to add a root node to and then generates binary tree.
    Root node should be at the midpoint of the last added edge.

Arguments:
    To be determined
Returns:
    To be determined
'''
def assemble_tree(edges):
    ''' Complete this function. '''
    pass


''' Returns a string of the Newick tree format for the tree rooted at `root`.

Arguments:
    To be determined
Returns:
    output: rooted tree in Newick tree format (string)
'''
def generate_newick(root, tree_map, D, mapping = None):
    ''' Complete this function. '''
    # TODO: print with 6 decimal digits such as
    # string = '%s:%.6f' % (name, length)
    pass
 

def main():
    parser = argparse.ArgumentParser(
        description='Neighbor-joining algorithm on a set of n sequences')
    parser.add_argument('-f', action="store", dest="f", type=str, default='dist10.txt')
    args = parser.parse_args()
    distances_file = args.f

    D, mapping = read_data(distances_file)
    nwk_str = neighbor_join(D,mapping) # TODO: decide what to return here
    #assemble_tree() # TODO: decide arguments to pass and what to return here if needed
    #nwk_str = generate_newick() # TODO: decide arguments to pass here if needed
    
    # Print and save the Newick string.
    print(nwk_str)
    with open("tree.nwk", "w") as nwk_file:
        print(nwk_str, file=nwk_file)


if __name__ == "__main__":
    main()
