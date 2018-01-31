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



def FIND_STATISTICS(err_ref_cont_coord_errors_list):
    

    insertions_dict={'insertion':0, 'insertion-multiple_copy':0,'insertion-tandem_multiple_copy':0,
                     'wrong_beginning':0,'wrong_end':0,'wrong_gap':0,'wrong_gap-overestimated':0}

    deletions_dict={'deletion':0,'deletion-collapsed_repeat':0,'deletion-collapsed_tandem_repeat':0,'deletion-gap_underestimated':0}

    substitutions_dict={'substitution':0,'gap':0}

    misjoin_dict={'circular_genome_start':0,'misjoin':0,'misjoin-insertion':0,'misjoin-mixed_fragments':0,'misjoin-wrong_scaffolding':0,'misjoin-overlap':0}
    translocation_dict={'translocation':0,'translocation-insertion':0,'translocation-mixed_fragments':0,'translocation-wrong_scaffolding':0,'translocation-overlap':0}
    wrong_seq_num=0
    reshufflings_block_num=0
    reshufflings_num=0
    inversions_num=0
    uncov_ref_frag_num=0
    uncov_ref_frag_len=0
   
    for entry in err_ref_cont_coord_errors_list:
        

        if entry[8]=='r':
            if entry[6]=='contig_st' or entry[6]=='contig_end':
                a='contig_end, contig_st. Do nothing.'
            elif entry[6]=='uncovered_region':
                uncov_ref_frag_num+=1
                uncov_ref_frag_len+=entry[7]
           
           
        elif entry[8]=='a':

            if entry[6]=='wrong_sequence':
                wrong_seq_num+=1
            
        
        
        elif entry[8]=='struct' or entry[8]=='snps':

            if entry[6] in insertions_dict:
                insertions_dict[entry[6]]+=1
            elif entry[6] in deletions_dict:
                deletions_dict[entry[6]]+=1
            elif entry[6] in substitutions_dict:
                substitutions_dict[entry[6]]+=1
            elif entry[6].startswith('reshuffling'):
                
                temp=entry[6].split('reshuffling-part_')
                
                if temp[1][:2]=='1_':
                    reshufflings_num+=1
                    
                reshufflings_block_num+=1
            elif entry[6]=='inversion':
                inversions_num+=1
            
           
        elif entry[8]=='struct2':
            
            if entry[3] in misjoin_dict:
                
                misjoin_dict[entry[3]]+=1

                if entry[3]in ['misjoin-insertion','misjoin-mixed_fragments'] :
                    insertions_dict['insertion']+=1

                elif entry[3]in ['misjoin-wrong_scaffolding'] :
                    insertions_dict['wrong_gap']+=1
                        
            elif entry[3] in translocation_dict:
                translocation_dict[entry[3]]+=1

                if entry[3]in ['translocation-insertion','translocation-mixed_fragments'] :
                    insertions_dict['insertion']+=1

                elif entry[3]in ['translocation-wrong_scaffolding'] :
                    insertions_dict['wrong_gap']+=1
                
            
        

    insertion_num=sum(insertions_dict.values())
    deletion_num=sum(deletions_dict.values())
    substitution_num=sum(substitutions_dict.values())
    misjoin_num=sum(misjoin_dict.values())
    translocation_num=sum(translocation_dict.values())

    total_num=insertion_num+deletion_num+substitution_num+misjoin_num+translocation_num+wrong_seq_num+reshufflings_num+inversions_num

    statistics_output_lines=[]
    statistics_output_lines.append('Total number\t'+str(total_num)+'\n')
    statistics_output_lines.append('Insertions\t'+str(insertion_num)+'\n')
    statistics_output_lines.append('Deletions\t'+str(deletion_num)+'\n')
    statistics_output_lines.append('Substitutions\t'+str(substitution_num)+'\n')
    statistics_output_lines.append('Translocations\t'+str(translocation_num)+'\n')
    statistics_output_lines.append('Relocations\t'+str(misjoin_num - misjoin_dict['circular_genome_start'])+'\n')
    statistics_output_lines.append('Reshufflings\t'+str(reshufflings_num)+'\n')
    statistics_output_lines.append('Reshuffled blocks\t'+str(reshufflings_block_num)+'\n')
    statistics_output_lines.append('Inversions\t'+str(inversions_num)+'\n')
    statistics_output_lines.append('Unaligned sequences\t'+str(wrong_seq_num)+'\n')
    statistics_output_lines.append('\n')
    statistics_output_lines.append('Uncovered ref regions num\t'+str(uncov_ref_frag_num)+'\n')
    statistics_output_lines.append('Uncovered ref regions len\t'+str(uncov_ref_frag_len)+'\n')
    statistics_output_lines.append('\n')
    statistics_output_lines.append('DETAILED INFORMATION:\t\n')
    statistics_output_lines.append('substitution\t'+str(substitutions_dict['substitution'])+'\n')
    statistics_output_lines.append('gap\t'+str(substitutions_dict['gap'])+'\n')
    statistics_output_lines.append('\n')
    statistics_output_lines.append('insertion\t'+str(insertions_dict['insertion'])+'\n')
    statistics_output_lines.append('duplication\t'+str(insertions_dict['insertion-multiple_copy'])+'\n')
    statistics_output_lines.append('tandem_duplication\t'+str(insertions_dict['insertion-tandem_multiple_copy'])+'\n')
    statistics_output_lines.append('unaligned_beginning\t'+str(insertions_dict['wrong_beginning'])+'\n')
    statistics_output_lines.append('unaligned_end\t'+str(insertions_dict['wrong_end'])+'\n')
    statistics_output_lines.append('inserted_gap\t'+str(insertions_dict['wrong_gap']+insertions_dict['wrong_gap-overestimated'])+'\n')
    
    statistics_output_lines.append('\n')
    statistics_output_lines.append('deletion\t'+str(deletions_dict['deletion']+deletions_dict['deletion-gap_underestimated'])+'\n')
    statistics_output_lines.append('collapsed_repeat\t'+str(deletions_dict['deletion-collapsed_repeat'])+'\n')
    statistics_output_lines.append('tandem_collapsed_repeat\t'+str(deletions_dict['deletion-collapsed_tandem_repeat'])+'\n')
    
    statistics_output_lines.append('\n')
    statistics_output_lines.append('translocation\t'+str(translocation_dict['translocation'])+'\n')
    statistics_output_lines.append('translocation-insertion\t'+str(translocation_dict['translocation-insertion'])+'\n')
    statistics_output_lines.append('translocation-insertion_ATGCN\t'+str(translocation_dict['translocation-mixed_fragments'])+'\n')
    statistics_output_lines.append('translocation-inserted_gap\t'+str(translocation_dict['translocation-wrong_scaffolding'])+'\n')
    statistics_output_lines.append('translocation-overlap\t'+str(translocation_dict['translocation-overlap'])+'\n')
    statistics_output_lines.append('\n')
    statistics_output_lines.append('circular_genome_start\t'+str(misjoin_dict['circular_genome_start'])+'\n')
    statistics_output_lines.append('relocation\t'+str(misjoin_dict['misjoin'])+'\n')
    statistics_output_lines.append('relocation-insertion\t'+str(misjoin_dict['misjoin-insertion'])+'\n')
    statistics_output_lines.append('relocation-insertion_ATGCN\t'+str(misjoin_dict['misjoin-mixed_fragments'])+'\n')
    statistics_output_lines.append('relocation-inserted_gap\t'+str(misjoin_dict['misjoin-wrong_scaffolding'])+'\n')
    statistics_output_lines.append('relocation-overlap\t'+str(misjoin_dict['misjoin-overlap'])+'\n')
    
    

    return statistics_output_lines
        
