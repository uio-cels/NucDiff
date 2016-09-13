# copyright (c) 2014-2016 Ksenia Khelik
#
# This Source Code Form is subject to the 
# terms of the Mozilla Public License, v. 2.0. 
# If a copy of the MPL was not distributed 
# with this file, You can obtain one at 
# https://mozilla.org/MPL/2.0/.
#
#-------------------------------------------------------------------------------



import os


def READ_FASTA_ENTRY(file_name):
    fasta_sequences=[]
    sequence_name=[]
    full_names_dict={}

    sequences_dict={}

    
    if os.stat(file_name)[6]!=0: #not empty
        


        f=open(file_name,'r')

        line=f.readline()
        
        seq=''
        
        while line:
            if line.startswith('>'):
                if not seq=='':
                    fasta_sequences.append(seq)
                    seq=''
                sequence_name.append(line[1:-1].split(' ')[0])
                full_names_dict[line[1:-1].split(' ')[0]]=line[1:-1]
                
            else:
                
                seq+=line[:-1]
                
            
            last_simb=line[-1]                
            line=f.readline()

        

        if last_simb!='\n' and last_simb!=' ':
            seq+=last_simb
        fasta_sequences.append(seq)


        for i in range(len(sequence_name)):
            sequences_dict[sequence_name[i]]=fasta_sequences[i]
        
    return sequences_dict,fasta_sequences, sequence_name, full_names_dict


def COMPL_STRING(line):
    subst_dict={'A':'T','T':'A','G':'C','C':'G','N':'N','-':'-','a':'T','t':'A','g':'C','c':'G','n':'N'}

    
    compl_line=''
    for i in range(len(line)-1,-1,-1):
        if subst_dict.has_key(line[i]):
            compl_line+=subst_dict[line[i]]
        else:
            compl_line+=line[i]
    

    return compl_line




