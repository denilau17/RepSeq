import sqlite3
import csv
import time

#retrive list of all clone family names
def get_clone_list(c, table):
    query = 'select distinct cloneID from ' + table + ';'
    results = c.execute(query).fetchall()
    results = [y[0] for y in results]
    return results

#get all sequences and sequence IDs for a clonal family
def get_seqs(c, cloneID, table):
    query = 'select Sequence, seqID from ' + table + ' where cloneID =?'
    clone = (cloneID,)
    results = c.execute(query, clone).fetchall()
    return results

#make a FASTA file
def write_csv(cloneID, seqs, out_folder):
    outfile = out_folder + str(cloneID) + ".fasta"
    csvfile = open(outfile, "wb")
    for seq in seqs:
        fasta_id = seq[1].split(' ')[-1]
        csvfile.write('>' + fasta_id + '\n')
        csvfile.write(seq[0] + '\n')

def main(db, table, out_folder):
    conn = sqlite3.connect(db)
    conn.text_factory = str
    c = conn.cursor()

    clone_list = get_clone_list(c, table)
    clone_list= clone_list[1:10]
    seq_list = [get_seqs(c, clone, table) for clone in clone_list]

    for i in range(len(clone_list)):
        write_csv(clone_list[i], seq_list[i], out_folder)
    
    return seq_list

if __name__ == "__main__":
    db = 'clone_assigned.sqlite'
    table = 'IMGT_007'
    out_folder = 'FASTA_007/'
    out = main(db, table, out_folder)
