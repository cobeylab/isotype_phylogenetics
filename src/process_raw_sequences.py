#!/usr/bin/python

"""Processes raw fasta files to obtain the reverse complement of raw sequences
"""

import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


raw_file_names = ['../data/01107PB.fasta', '../data/01207PB.fasta', '../data/SFAPB.fasta']

complement = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

test_sequence = Seq('AAGG', generic_dna)
assert str(test_sequence.reverse_complement()) == 'CCTT'

def main(argv):

    for file_name in raw_file_names:

        print 'Processing ' + file_name

        processed_file_name = file_name.replace('.fasta','_processed.fasta')

        first_line = True

        with open(file_name, 'r') as raw_file, open(processed_file_name,'w') as processed_file:
            for line in raw_file:
                #If line contains a sequence identifier:
                if line.find('>') > -1:
                    seq_id = re.search(r'>[A-Z|0-9]*', line)
                    seq_id = seq_id.group()

                    if first_line == True:
                        processed_file.write(seq_id + '\n')
                        first_line = False

                    else:
                        processed_file.write('\n' + seq_id + '\n')

                #If not, then it contains a sequence:
                else:
                    sequence = line.replace('\n','')
                    sequence = Seq(sequence, generic_dna)
                    sequence = sequence.reverse_complement()

                    #sequence.translate(table=1)

                    processed_file.write(str(sequence))


if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)






