# copyright (c) 2014-2016 Ksenia Khelik
#
# This Source Code Form is subject to the 
# terms of the Mozilla Public License, v. 2.0. 
# If a copy of the MPL was not distributed 
# with this file, You can obtain one at 
# https://mozilla.org/MPL/2.0/.
#
#-------------------------------------------------------------------------------

from __future__ import print_function

from Bio.Seq import Seq
from Bio import SeqIO

import os


def READ_FASTA_ENTRY(file_name):
    fasta_sequences=[]
    sequence_name=[]
    full_names_dict={}

    sequences_dict={}

    
    if os.stat(file_name)[6]!=0: #not empty
       
        fh = open(file_name, "r")
        for record in SeqIO.parse(fh, "fasta"):
            short_name=str(record.id).split(' ')[0]
            sequence_name.append(short_name)
            full_names_dict[short_name]=str(record.id)

            fasta_sequences.append(str(record.seq))
            sequences_dict[short_name]=str(record.seq)

        
    return sequences_dict,fasta_sequences, sequence_name, full_names_dict


def COMPL_STRING(line):

    my_seq = Seq(line)
    compl_line=str(my_seq.reverse_complement())

    return compl_line




