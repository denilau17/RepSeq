import sys
import os
import csv
import sqlite3

import clusteringcore
import numpy as np
import scipy as sp
import scipy.cluster

import datetime
import multiprocessing as mp
import itertools

def pdist(X,metric):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = metric(X[i], X[j])
            k += 1
    return dm

def cluster_seqs(seqs,cutoff,linkage='single'):
    if len(seqs) == 0:
        return (np.array([]),{})

    #collapses identical seqs
    #unique_seqs = list(set(seqs))
    #seq_idxs = dict( [(j,i) for (i,j) in enumerate(unique_seqs)] )

    #checks if there is only 1 unique seq
    if len(seqs) == 1:
        T = np.array([1]*len(seqs))
        return T

    #compute distance matrix
    #Y = pdist( unique_seqs, clusteringcore.levenshtein )
    Y = pdist(seqs, clusteringcore.levenshtein)
    
    #compute linkage
    Z = sp.cluster.hierarchy.linkage(Y,method=linkage)

    # determine the clusters at level cutoff
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')
    #return (T,seq_idxs)
    return T

#get list of subgroups for each pool of clonal assignments
def get_subgroups(c, subject):
    query = "SELECT subgroup FROM " + subject + ";"
    results = c.execute(query).fetchall()
    subgroup = [x[0].encode('ascii', 'ignore') for x in results]
    return subgroup

#get sequence, celltype and CDR3 len info for clustering and post-clustering analysis
def get_subgroup_seqs(c, subgroup):
    query = "SELECT Sequence, cell_type, CDR3_len FROM " + subject + " WHERE subgroup = '" + subgroup + "';"
    results = c.execute(query).fetchall()
    seqs = [x[0].encode('ascii', 'ignore') for x in results]
    seqs = list(set(seqs))
    celltype = [x[1].encode('ascii', 'ignore') for x in results]
    cdr3_len = results[0][2]
    return [seqs, celltype, cdr3_len]

#clusters a particular subgroup into clones with max edit distance of the CDR3 len 
def assign_clones(c, subgroup):
    seqs, celltype, cdr3_len = get_subgroup_seqs(c, subgroup)
    seqs = list(set(seqs))
    results = cluster_seqs(seqs, cdr3_len)
    t = [int(x) for x in results]
        
    for i in range(len(t)):
        rv.append([subgroup_list[1], t[i], celltype[i], seqs[i], cdr3_len])

def clones(data):
    results = cluster_seqs(data[0], data[2])
    t = [int(x) for x in results]
    return t

def format_data(subgroup_list, data, results):
    if len(data) != len(results):
        return []
    rv = []
    for i in range(len(data)):
        subgroup = subgroup_list[i]
        seqs = data[i][0]
        celltype = data[i][1]
        cdr3_len = data[i][2]
        clone_assignments = results[i]
        for j in range(len(seqs)):
            if len(seqs) != len(clone_assignments):
                print("not correct order!!")
                return []
            row = [subgroup, cdr3_len, seqs[j], celltype[j], clone_assignments[j]]
            rv.append(row)
    
    return rv
        
    
def main(db, subject, outfile):

    connection = sqlite3.connect(db)
    c = connection.cursor()
    
    subgroup_list = get_subgroups(c, subject)

    subgroup_list = list(set(subgroup_list))

    subgroup_list = subgroup_list[0:200]
    data = []
    
    for subgroup in subgroup_list:
        x = get_subgroup_seqs(c, subgroup)
        data.append(x)

    pool = mp.Pool(processes=4)

    results = pool.map(clones, data)

    rv = format_data(subgroup_list, data, results)

    out = open(outfile, 'wb')
    csv_out = csv.writer(out)
    csv_out.writerows(rv)

    connection.close()

if __name__ == "__main__":
    db = '/Users/denise/Documents/RepSeq2/IMGT_parsed.sqlite'
    subject = 'IMGT_007'
    outfile = '/Users/denise/Documents/RepSeq2/clones_007.csv'

    main(db, subject, outfile)
