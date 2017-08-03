#!/usr/bin/python

"""Processes raw fasta files to obtain the reverse complement of raw sequences
"""

import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO


raw_file_names = ['../data/01107PB.fasta', '../data/01207PB.fasta', '../data/SFAPB.fasta']

complement = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

test_sequence = Seq('AAGG', generic_dna)
assert str(test_sequence.reverse_complement()) == 'CCTT'

def main(argv):

    for file_name in raw_file_names:

        print 'Processing ' + file_name.replace('../data/','')

        processed_file_name = file_name.replace('.fasta','_processed.fasta')

        first_line = True

        with open(file_name, 'rU') as raw_file, open(processed_file_name,'w') as processed_file:

            fasta_object = list(SeqIO.parse(raw_file, "fasta"))

            for record in fasta_object:
                record.seq = record.seq.reverse_complement()
                record.description = ''

            SeqIO.write(fasta_object, processed_file, 'fasta')


if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)






