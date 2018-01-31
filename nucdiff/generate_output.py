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

from . import general

err_new_names_dict={'substitution':'substitution', 'gap':'gap',

                      'insertion':'insertion','insertion-multiple_copy':'duplication',
                      'insertion-tandem_multiple_copy':'tandem_duplication','wrong_beginning':'unaligned_beginning',
                      'wrong_end':'unaligned_end', 'wrong_gap':'inserted_gap', 'wrong_gap-overestimated':'inserted_gap',

                      'deletion':'deletion', 'deletion-collapsed_repeat':'collapsed_repeat',
                      'deletion-collapsed_tandem_repeat':'collapsed_tandem_repeat','deletion-gap_underestimated':'deletion',

                      'translocation':'translocation','translocation-insertion':'translocation-insertion',
                      'translocation-mixed_fragments':'translocation-insertion_ATGCN',
                      'translocation-wrong_scaffolding':'translocation-inserted_gap',
                      'translocation-overlap':'translocation-overlap',

                      'circular_genome_start':'circular_genome_start', 'misjoin':'relocation',
                      'misjoin-insertion':'relocation-insertion',
                      'misjoin-mixed_fragments':'relocation-insertion_ATGCN',
                      'misjoin-wrong_scaffolding':'relocation-inserted_gap',
                      'misjoin-overlap':'relocation-overlap',

                      'reshuffling':'reshuffling', 'inversion':'inversion','forward':'forward',
                      'wrong_sequence':'unaligned_sequence',

                      'translocation_st':'translocation_st','translocation-insertion_st':'translocation-insertion_st',
                      'translocation-mixed_fragments_st':'translocation-insertion_ATGCN_st',
                      'translocation-wrong_scaffolding_st':'translocation-inserted_gap_st',
                      'translocation-overlap_st':'translocation-overlap_st',

                      'translocation_end':'translocation_end','translocation-insertion_end':'translocation-insertion_end',
                      'translocation-mixed_fragments_end':'translocation-insertion_ATGCN_end',
                      'translocation-wrong_scaffolding_end':'translocation-inserted_gap_end',
                      'translocation-overlap_end':'translocation-overlap_end',


                      'misjoin_st':'relocation_st', 'misjoin-insertion_st':'relocation-insertion_st',
                      'misjoin-mixed_fragments_st':'relocation-insertion_ATGCN_st',
                      'misjoin-wrong_scaffolding_st':'relocation-inserted_gap_st',
                      'misjoin-overlap_st':'relocation-overlap_st',

                      'misjoin_end':'relocation_end','misjoin-insertion_end':'relocation-insertion_end',
                      'misjoin-mixed_fragments_end':'relocation-insertion_ATGCN_end',
                      'misjoin-wrong_scaffolding_end':'relocation-inserted_gap_end',
                      'misjoin-overlap_end':'relocation-overlap_end',     

                      'contig_st':'query_st', 'contig_end':'query_end', 'uncovered_region':'uncovered_region',
                      'clipped_repeated_region':'uncovered_repeated_region'
                      }




def OUTPUT_LOCAL_PART(f,structure_dict_local, cont_name, level_line):

    
        if structure_dict_local!=[]: #local_differences
            for entry in structure_dict_local:
                    f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

   

def OUTPUT_INVERS_PART(f,structure_dict_block, cont_name,  level_line):
    
        if structure_dict_block['invers_output']==[]: #no inversion
            OUTPUT_LOCAL_PART(f,structure_dict_block['local_output'], cont_name, level_line)

        else: #inversion
            if structure_dict_block['transp_output']!=[]:
            
                    entry=structure_dict_block['invers_output'][0]
                    f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

                    
                    OUTPUT_LOCAL_PART(f,structure_dict_block['local_output'], cont_name,  '\t'+level_line)

            else:
                if structure_dict_block['between_output']!=[] :
                    if structure_dict_block['invers_output'][0][0]>structure_dict_block['between_output'][0][0]:
                        entry=structure_dict_block['between_output'][0]
                        f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

                        entry=structure_dict_block['invers_output'][0]
                        f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

                    
                        OUTPUT_LOCAL_PART(f,structure_dict_block['local_output'], cont_name, '\t'+level_line)

                        if len(structure_dict_block['between_output'])>1:
                            for entry in structure_dict_block['between_output'][1:]:
                                f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
                    else:

                        entry=structure_dict_block['invers_output'][0]
                        f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

                    
                        OUTPUT_LOCAL_PART(f,structure_dict_block['local_output'], cont_name, '\t'+level_line)

                        for entry in structure_dict_block['between_output']:
                                f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
                    
                        
                
                else:

                    entry=structure_dict_block['invers_output'][0]
                    f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

                    
                    OUTPUT_LOCAL_PART(f,structure_dict_block['local_output'], cont_name, '\t'+level_line)

                       
    
def OUTPUT_TRANSP_PART(f,structure_dict_misj, cont_name, level_line):

    
        if len(list(structure_dict_misj.keys()))==1: #no transposition
            block_name=list(structure_dict_misj.keys())[0]
            OUTPUT_INVERS_PART(f,structure_dict_misj[block_name], cont_name, level_line)

        else: #transposition
            for block_name in sorted(structure_dict_misj.keys()):
                if structure_dict_misj[block_name]['transp_output']!=[]:

                    if structure_dict_misj[block_name]['between_output']!=[]:
                        if structure_dict_misj[block_name]['transp_output'][0][0]>structure_dict_misj[block_name]['between_output'][0][0]:
                            entry=structure_dict_misj[block_name]['between_output'][0]
                            f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')

                            entry=structure_dict_misj[block_name]['transp_output'][0]
                            f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+entry[2]+'\t'+str(entry[3])+'\n')

                            OUTPUT_INVERS_PART(f,structure_dict_misj[block_name], cont_name,  '\t'+level_line)

                            if len(structure_dict_misj[block_name]['between_output'])>1:
                                for entry in structure_dict_misj[block_name]['between_output'][1:]:
                                    f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
                        else:
                            entry=structure_dict_misj[block_name]['transp_output'][0]
                            f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+entry[2]+'\t'+str(entry[3])+'\n') 

                            OUTPUT_INVERS_PART(f,structure_dict_misj[block_name], cont_name, '\t'+level_line)

                            for entry in structure_dict_misj[block_name]['between_output']:
                                    f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
                        


                            
                    else:        
                        
                        entry=structure_dict_misj[block_name]['transp_output'][0]
                        f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+entry[2]+'\t'+str(entry[3])+'\n')

                        OUTPUT_INVERS_PART(f,structure_dict_misj[block_name], cont_name, '\t'+level_line)

                        
                else:
                    OUTPUT_INVERS_PART(f,structure_dict_misj[block_name], cont_name,  level_line)

                       
        
            
    
    

def OUTPUT_MISJOIN_PART(f,structure_dict_transl, cont_name, level_line):

    
        if len(list(structure_dict_transl.keys()))==1: #no misjoin
            misj_group=list(structure_dict_transl.keys())[0]
            OUTPUT_TRANSP_PART(f,structure_dict_transl[misj_group]['blocks'], cont_name,  level_line)

        else: #misjoin
            for misjoin_group in sorted(structure_dict_transl.keys()):
                    entry=structure_dict_transl[misjoin_group]['output_line']
                    f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict['misjoin']+'_'+str(entry[2])+'\t'+str(entry[1]-entry[0]+1)+'\n')

                    OUTPUT_TRANSP_PART(f,structure_dict_transl[misjoin_group]['blocks'], cont_name, '\t'+level_line)

                    if structure_dict_transl[misjoin_group]['reason']!=[]:
                        for entry in structure_dict_transl[misjoin_group]['reason']:
                            f.write(level_line+'|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
            
    

def OUTPUT_TRANSLOC_PART(f,structure_dict_cont, cont_name):

    
        if len(list(structure_dict_cont.keys()))==1: #no translocation
            OUTPUT_MISJOIN_PART(f,structure_dict_cont[0]['blocks'], cont_name, '')

        else: #translocation
            for transl_group in sorted(structure_dict_cont.keys()):
                    entry=structure_dict_cont[transl_group]['output_line']
                    f.write('|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict['translocation']+'_'+str(entry[2])+'\t'+str(entry[1]-entry[0]+1)+'\n')

                    OUTPUT_MISJOIN_PART(f,structure_dict_cont[transl_group]['blocks'], cont_name, '\t')

                    if structure_dict_cont[transl_group]['reason']!=[]:
                        for entry in structure_dict_cont[transl_group]['reason']:
                            f.write('|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
            
                    

   
        

def OUTPUT_READABLE(structure_dict,end_err_dict,unmapped_list, cont_dict,contig_names, contig_full_names_dict, ref_dict,ref_full_names_dict,working_dir, prefix):

        

        file_name=working_dir+prefix+'_errors_readable.out'

        f=open(file_name, 'w')

        for cont_name in contig_names:
            cont_name=cont_name.split(' ')[0]
            if cont_name in end_err_dict and (end_err_dict[cont_name]['wrong_beginning']!=[] or end_err_dict[cont_name]['wrong_end']!=[]):
                f.write(contig_full_names_dict[cont_name]+'\t'+str(len(cont_dict[cont_name]))+'\n' )

                if  end_err_dict[cont_name]['wrong_beginning']!=[]:
                        entry=end_err_dict[cont_name]['wrong_beginning']
                        f.write('|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
                        
                        
                        

                if len(list(structure_dict[cont_name].keys()))==1 and \
                       len(list(structure_dict[cont_name][0]['blocks'].keys()))==1 and \
                       len(list(structure_dict[cont_name][0]['blocks'][0]['blocks'].keys()))==1 and \
                       structure_dict[cont_name][0]['blocks'][0]['blocks'][0]['block'][10]==[]:
                            a='no errors inside'
                else:
                    OUTPUT_TRANSLOC_PART(f,structure_dict[cont_name], cont_name)


                if  end_err_dict[cont_name]['wrong_end']!=[]:
                        entry=end_err_dict[cont_name]['wrong_end']
                        f.write('|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')


                f.write('\n\n')

            else:
                if cont_name in unmapped_list:
                    f.write(contig_full_names_dict[cont_name]+'\t'+str(len(cont_dict[cont_name]))+'\n' )
                    f.write('|1\t'+str(len(cont_dict[cont_name]))+'\t'+err_new_names_dict['wrong_sequence']+'\n')
                    f.write('\n\n')

                else:
                    if len(list(structure_dict[cont_name].keys()))==1 and \
                       len(list(structure_dict[cont_name][0]['blocks'].keys()))==1 and \
                       len(list(structure_dict[cont_name][0]['blocks'][0]['blocks'].keys()))==1 and \
                       structure_dict[cont_name][0]['blocks'][0]['blocks'][0]['block'][10]==[]:
                            a='correct_query'
                    else:
                        f.write(contig_full_names_dict[cont_name]+'\t'+str(len(cont_dict[cont_name]))+'\n' )
                        OUTPUT_TRANSLOC_PART(f,structure_dict[cont_name], cont_name)
                        f.write('\n\n')
                        
        f.close()


    
def OUTPUT_VCF_FILE_QUERY(asmb_name_full,contig_full_names_dict, contig_names,snp_err_dict,output_file,cont_dict,ref_dict):

    fq=open(output_file,'w')
    
    fq.write('##fileformat=VCFv4.2\n')
    fq.write('##source=NucDiffv2.0\n')
    fq.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n')
   
    for cont_name in contig_names:
        if asmb_name_full=='yes':
            c_name=contig_full_names_dict[cont_name]
        else:
            c_name=cont_name
                

        if snp_err_dict[cont_name]!=[]:
                for entry in snp_err_dict[cont_name]:
                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap':
                            fq.write(c_name+'\t'+str(entry[4]-1)+'\t.\t'+cont_dict[cont_name][entry[4]-2:entry[5]].upper()+
                                         '\t'+cont_dict[cont_name][entry[4]-2:entry[4]-1].upper()+'\t.\tPASS\n')
        
                            
                        elif entry[6]=='deletion':
                            
                            if entry[9]==1:
                                 ref_bases=ref_dict[entry[0]][entry[1]-1:entry[2]].upper()
                            else:
                                 ref_bases=general.COMPL_STRING(ref_dict[entry[0]][entry[1]-1:entry[2]]).upper()

                            fq.write(c_name+'\t'+str(entry[4])+'\t.\t'+cont_dict[cont_name][entry[4]-1:entry[4]].upper()+
                                         '\t'+cont_dict[cont_name][entry[4]-1:entry[4]].upper()+ref_bases+'\t.\tPASS\n')
                         
                            
                           
                        else: #subst and gap
                            ref_coord=str(entry[1])+'-'+str(entry[2])
                            if entry[9]==1:
                                 ref_bases=ref_dict[entry[0]][entry[1]-1:entry[2]].upper()
                            else:
                                 ref_bases=general.COMPL_STRING(ref_dict[entry[0]][entry[1]-1:entry[2]]).upper()

                                
                            fq.write(c_name+'\t'+str(entry[4])+'\t.\t'+cont_dict[cont_name][entry[4]-1:entry[5]].upper()+
                                         '\t'+ref_bases+'\t.\tPASS\n')
                         
                           
                       
                        
    
    fq.close()

        
def OUTPUT_VCF_FILE_REF(ref_name_full,ref_full_names_dict, ref_names,ref_snp_err_dict,output_file,ref_dict,cont_dict):
    fr=open(output_file,'w')

    fr.write('##fileformat=VCFv4.2\n')
    fr.write('##source=NucDiffv2.0\n')
    fr.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n')

    for ref_name in ref_names:
        if ref_name_full=='yes':
            r_name=ref_full_names_dict[ref_name]
        else:
            r_name=ref_name
                

        if ref_snp_err_dict[ref_name]!=[]:
                
                for entry in ref_snp_err_dict[ref_name]:
                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap':
                            if entry[9]==1:
                                 query_bases=cont_dict[entry[3]][entry[4]-1:entry[5]].upper()
                            else:
                                 query_bases=general.COMPL_STRING(cont_dict[entry[3]][entry[4]-1:entry[5]]).upper()

                            fr.write(r_name+'\t'+str(entry[1])+'\t.\t'+ref_dict[r_name][entry[1]-1:entry[1]].upper()+
                                         '\t'+ref_dict[r_name][entry[1]-1:entry[1]].upper()+query_bases+'\t.\tPASS\n')
                         

                         

                        elif entry[6]=='deletion':
                            fr.write(r_name+'\t'+str(entry[1]-1)+'\t.\t'+ref_dict[r_name][entry[1]-2:entry[2]].upper()+
                                         '\t'+ref_dict[r_name][entry[1]-2:entry[1]-1].upper()+'\t.\tPASS\n')
        

                            
                            
                        
                        else: #subst and gap
                            if entry[9]==1:
                                query_bases=cont_dict[entry[3]][entry[4]-1:entry[5]].upper()
                            else:
                                query_bases=general.COMPL_STRING(cont_dict[entry[3]][entry[4]-1:entry[5]]).upper()
                         
                            fr.write(r_name+'\t'+str(entry[1])+'\t.\t'+ref_dict[r_name][entry[1]-1:entry[2]].upper()+
                                         '\t'+query_bases+'\t.\tPASS\n')
                         
                           
    
    
    fr.close()


def OUTPUT_REF_ASSEM_TABLE(err_ref_cont_coord_errors_list, ref_dict,ref_names,ref_full_names_dict,cont_dict,contig_names,contig_full_names_dict,working_dir, prefix,asmb_name_full,ref_name_full,mapped_blocks_dict,vcf_flag):

   #query_snps.gff
    temp_entry=[]
    
    snp_err_dict={}
    struct_err_dict={}
    repeat_region_dict={}
    for cont_name in list(cont_dict.keys()):
        snp_err_dict[cont_name]=[]
        struct_err_dict[cont_name]=[]
        repeat_region_dict[cont_name]=[]
        

   
        
    for entry in err_ref_cont_coord_errors_list:
        
        if entry[8]=='a':
            c_name=entry[3]
            struct_err_dict[c_name].append(entry)
        elif entry[8]=='r':
            a='contig_end contig_st'
        elif entry[8]=='snps':
            c_name=entry[3]
            snp_err_dict[c_name].append(entry)
        elif entry[8]=='struct':
            c_name=entry[3]
            struct_err_dict[c_name].append(entry)
        else:
            c_name=entry[0]
            struct_err_dict[c_name].append(['.','.','.']+entry)

     
            
        
    for c_name in list(snp_err_dict.keys()):
        snp_err_dict[c_name]= sorted(snp_err_dict[c_name],key=lambda inter:inter[4], reverse=False)
            
    for c_name in list(struct_err_dict.keys()):
        struct_err_dict[c_name]= sorted(struct_err_dict[c_name],key=lambda inter:inter[4], reverse=False)

        
        for entry in struct_err_dict[c_name]:
            if entry[:3]==['.','.','.']:
                for u in range(3):
                     entry.pop(0)
                
       
    if vcf_flag=='yes':
            OUTPUT_VCF_FILE_QUERY(asmb_name_full,contig_full_names_dict, contig_names,snp_err_dict,working_dir+prefix+'_query_snps.vcf',cont_dict,ref_dict)

    fq=open(working_dir+prefix+'_query_snps.gff','w')
    
    fq.write('##gff-version 3\n')
    

    ID_cur=1
    for cont_name in contig_names:
        if asmb_name_full=='yes':
            c_name=contig_full_names_dict[cont_name]
        else:
            c_name=cont_name
                

        if snp_err_dict[cont_name]!=[]:
                fq.write('##sequence-region\t'+c_name+'\t1\t'+str(len(cont_dict[cont_name]))+'\n')

                for i in range(len(snp_err_dict[cont_name])):
                    snp_err_dict[cont_name][i].append('SNP_'+str(ID_cur))
                    ID_cur+=1

                for entry in snp_err_dict[cont_name]:
                        if ref_name_full=='yes':
                            r_name=ref_full_names_dict[entry[0]]
                        else:
                            r_name=entry[0]

                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap':
                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000667'+'\t'+str(entry[4])+'\t'+str(entry[5])+
                                         '\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+
                                         ';query_bases='+cont_dict[cont_name][entry[4]-1:entry[5]].lower()+\
                                        ';ref_bases=-'+';color=#EE0000'+'\n')

                        
                        elif entry[6]=='deletion':
                            ref_coord=str(entry[1])+'-'+str(entry[2])
                            if entry[9]==1:
                                 ref_bases=ref_dict[entry[0]][entry[1]-1:entry[2]]
                            else:
                                 ref_bases=general.COMPL_STRING(ref_dict[entry[0]][entry[1]-1:entry[2]])
                         
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+ref_coord+';query_bases=-;ref_bases='+ref_bases.lower()+\
                                     ';color=#0000EE'+'\n')
                        
                        else: #subst and gap
                            ref_coord=str(entry[1])+'-'+str(entry[2])
                            if entry[9]==1:
                                 ref_bases=ref_dict[entry[0]][entry[1]-1:entry[2]]
                            else:
                                 ref_bases=general.COMPL_STRING(ref_dict[entry[0]][entry[1]-1:entry[2]])
                         
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[4])+'\t'+str(entry[5])+
                                     '\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';subst_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+
                                     ';ref_coord='+ref_coord+';query_bases='+cont_dict[cont_name][entry[4]-1:entry[5]].lower()+\
                                     ';ref_bases='+ref_bases.lower()+';color=#42C042'+'\n')
                         
                        if entry[4]>entry[5] or entry[7]<0:
                            print(entry)
                            print('ERROR: a wrong query coordinate')

                    
                        
    
    fq.close()

    #query struct
    
    fq=open(working_dir+prefix+'_query_struct.gff','w')
    
    fq.write('##gff-version 3\n')
    

    ID_cur=1
    for cont_name in contig_names:
        if asmb_name_full=='yes':
            c_name=contig_full_names_dict[cont_name]
        else:
            c_name=cont_name
                

        if struct_err_dict[cont_name]!=[]:
                fq.write('##sequence-region\t'+c_name+'\t1\t'+str(len(cont_dict[cont_name]))+'\n')

                for i in range(len(struct_err_dict[cont_name])):
                    struct_err_dict[cont_name][i].append('SV_'+str(ID_cur))
                    ID_cur+=1

                for entry in struct_err_dict[cont_name]:
                    
                    
                    if entry[8]=='struct':
                        if ref_name_full=='yes':
                            r_name=ref_full_names_dict[entry[0]]
                        else:
                            r_name=entry[0]
                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap' or entry[6]=='wrong_end' or entry[6]=='wrong_beginning' or entry[6]=='wrong_gap-overestimated':

                            if len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'
                            elif len(entry)==12:
                                ID_name=entry[11]
                                query_dir=str(entry[9])
                            elif len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                            
                            
                            
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000667'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+query_dir+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=#EE0000'+'\n')
                                
                        elif entry[6]=='insertion-multiple_copy':
                            
                            
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+\
                                        ';query_repeated_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=#EE0000'+'\n')

                                    
                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=#EE0000'+'\n')
                            elif len(entry)==12:
                                ID_name=entry[11]
                            
                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=#EE0000'+'\n')

                            else:
                                print(entry)
                                input('ku')
    

                        elif entry[6]=='insertion-tandem_multiple_copy':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+\
                                        ';query_repeated_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=#EE0000'+'\n')

                                elif len(entry[11])==2:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+\
                                        ';query_repeated_region_1='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+\
                                        ';query_repeated_region_2='+str(entry[11][1][4])+'-'+str(entry[11][1][5])+';color=#EE0000'+'\n')

                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=#EE0000'+'\n')
    
                            elif len(entry)==12:
                                ID_name=entry[11]
                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=#EE0000'+'\n')
                            
                            

                        elif entry[6]=='deletion' or entry[6]=='deletion-gap_underestimated':

                            if len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'
                            elif len(entry)==12:
                                ID_name=entry[11]
                                query_dir=str(entry[9])
                            elif len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                            
                            
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+query_dir+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#0000EE'+'\n')

                        elif entry[6]=='deletion-collapsed_repeat':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+\
                                        ';query_repeated_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=#0000EE'+'\n')

                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#0000EE'+'\n')
    
                            elif len(entry)==12:
                                ID_name=entry[11]

                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#0000EE'+'\n')

                            
                        elif entry[6]=='deletion-collapsed_tandem_repeat':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+\
                                        ';query_repeated_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=#0000EE'+'\n')

                                elif len(entry[11])==2:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+\
                                        ';query_repeated_region_1='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+\
                                        ';query_repeated_region_2='+str(entry[11][1][4])+'-'+str(entry[11][1][5])+';color=#0000EE'+'\n')

                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#0000EE'+'\n')
    
                            elif len(entry)==12:
                                ID_name=entry[11]
                                
                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#0000EE'+'\n')
                            

                        elif entry[6]=='substitution' or entry[6]=='gap':
                            
                            if len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'
                            elif len(entry)==12:
                                ID_name=entry[11]
                                query_dir=str(entry[9])
                            elif len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                        
                                                        
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';subst_len='+str(entry[7])+';query_dir='+query_dir+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#42C042'+'\n')

                        elif entry[6].startswith('reshuf'):
                            ID_name=entry[11]

                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+entry[6]+\
                                        ';blk_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#404040'+'\n')

                        elif entry[6]=='inversion':
                            ID_name=entry[11]

                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';blk_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=#01DFD7'+'\n')
                            
                        else:
                            print(entry)
                            #raw_input('kfj')
                            a='reshuffling inversion'
                            
                    elif entry[8]=='a':
                        ID_name=entry[9]

                        fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';color=#990000'+'\n')

                        

                    else:# entry[8]=='struct2'
                        
                        if entry[3]=='circular_genome_start' or entry[3]=='misjoin':
                            ID_name=entry[9]

                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';ref_sequence='+entry[6][0][4]+';blk_1_query='+str(entry[6][0][0])+'-'+str(entry[6][0][1])+\
                                     ';blk_1_ref='+str(entry[6][0][2])+'-'+str(entry[6][0][3])+\
                                     ';blk_1_query_len='+str(entry[6][0][1]-entry[6][0][0]+1)+\
                                     ';blk_1_ref_len='+str(entry[6][0][3]-entry[6][0][2]+1)+\
                                     ';blk_1_st_query='+str(entry[6][0][5][0])+';blk_1_st_ref='+str(entry[6][0][5][1])+\
                                     ';blk_1_end_query='+str(entry[6][0][6][0])+';blk_1_end_ref='+str(entry[6][0][6][1])+\
                                     ';blk_2_query='+str(entry[7][0][0])+'-'+str(entry[7][0][1])+\
                                     ';blk_2_ref='+str(entry[7][0][2])+'-'+str(entry[7][0][3])+\
                                     ';blk_2_query_len='+str(entry[7][0][1]-entry[7][0][0]+1)+\
                                     ';blk_2_ref_len='+str(entry[7][0][3]-entry[7][0][2]+1)+\
                                     ';blk_2_st_query='+str(entry[7][0][5][0])+';blk_2_st_ref='+str(entry[7][0][5][1])+\
                                     ';blk_2_end_query='+str(entry[7][0][6][0])+';blk_2_end_ref='+str(entry[7][0][6][1])+\
                                     ';color=#990099'+'\n')

                        elif entry[3]=='translocation':
                            ID_name=entry[9]

                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001873'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';ref_sequence_1='+entry[6][0][4]+\
                                     ';blk_1_query='+str(entry[6][0][0])+'-'+str(entry[6][0][1])+\
                                     ';blk_1_ref='+str(entry[6][0][2])+'-'+str(entry[6][0][3])+\
                                     ';blk_1_query_len='+str(entry[6][0][1]-entry[6][0][0]+1)+\
                                     ';blk_1_ref_len='+str(entry[6][0][3]-entry[6][0][2]+1)+\
                                     ';blk_1_st_query='+str(entry[6][0][5][0])+';blk_1_st_ref='+str(entry[6][0][5][1])+\
                                     ';blk_1_end_query='+str(entry[6][0][6][0])+';blk_1_end_ref='+str(entry[6][0][6][1])+\
                                     ';ref_sequence_2='+entry[7][0][4]+\
                                     ';blk_2_query='+str(entry[7][0][0])+'-'+str(entry[7][0][1])+\
                                     ';blk_2_ref='+str(entry[7][0][2])+'-'+str(entry[7][0][3])+\
                                     ';blk_2_query_len='+str(entry[7][0][1]-entry[7][0][0]+1)+\
                                     ';blk_2_ref_len='+str(entry[7][0][3]-entry[7][0][2]+1)+\
                                     ';blk_2_st_query='+str(entry[7][0][5][0])+';blk_2_st_ref='+str(entry[7][0][5][1])+\
                                     ';blk_2_end_query='+str(entry[7][0][6][0])+';blk_2_end_ref='+str(entry[7][0][6][1])+\
                                     ';color=#A0A0A0'+'\n')

                        elif entry[3]=='misjoin-insertion' or entry[3]=='misjoin-mixed_fragments' or entry[3]=='misjoin-wrong_scaffolding':
                            ID_name=entry[9]
                            
                            
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';ins_len='+str(entry[4])+\
                                     ';ref_sequence='+entry[6][0][4]+';blk_1_query='+str(entry[6][0][0])+'-'+str(entry[6][0][1])+\
                                     ';blk_1_ref='+str(entry[6][0][2])+'-'+str(entry[6][0][3])+\
                                     ';blk_1_query_len='+str(entry[6][0][1]-entry[6][0][0]+1)+\
                                     ';blk_1_ref_len='+str(entry[6][0][3]-entry[6][0][2]+1)+\
                                     ';blk_1_st_query='+str(entry[6][0][5][0])+';blk_1_st_ref='+str(entry[6][0][5][1])+\
                                     ';blk_1_end_query='+str(entry[6][0][6][0])+';blk_1_end_ref='+str(entry[6][0][6][1])+\
                                     ';blk_2_query='+str(entry[7][0][0])+'-'+str(entry[7][0][1])+\
                                     ';blk_2_ref='+str(entry[7][0][2])+'-'+str(entry[7][0][3])+\
                                     ';blk_2_query_len='+str(entry[7][0][1]-entry[7][0][0]+1)+\
                                     ';blk_2_ref_len='+str(entry[7][0][3]-entry[7][0][2]+1)+\
                                     ';blk_2_st_query='+str(entry[7][0][5][0])+';blk_2_st_ref='+str(entry[7][0][5][1])+\
                                     ';blk_2_end_query='+str(entry[7][0][6][0])+';blk_2_end_ref='+str(entry[7][0][6][1])+\
                                     ';color=#990099'+'\n')


                        elif entry[3]=='translocation-insertion' or entry[3]=='translocation-mixed_fragments' or entry[3]=='translocation-wrong_scaffolding':
                            ID_name=entry[9]

                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001873'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';ins_len='+str(entry[4])+\
                                     ';ref_sequence_1='+entry[6][0][4]+\
                                     ';blk_1_query='+str(entry[6][0][0])+'-'+str(entry[6][0][1])+\
                                     ';blk_1_ref='+str(entry[6][0][2])+'-'+str(entry[6][0][3])+\
                                     ';blk_1_query_len='+str(entry[6][0][1]-entry[6][0][0]+1)+\
                                     ';blk_1_ref_len='+str(entry[6][0][3]-entry[6][0][2]+1)+\
                                     ';blk_1_st_query='+str(entry[6][0][5][0])+';blk_1_st_ref='+str(entry[6][0][5][1])+\
                                     ';blk_1_end_query='+str(entry[6][0][6][0])+';blk_1_end_ref='+str(entry[6][0][6][1])+\
                                     ';ref_sequence_2='+entry[7][0][4]+\
                                     ';blk_2_query='+str(entry[7][0][0])+'-'+str(entry[7][0][1])+\
                                     ';blk_2_ref='+str(entry[7][0][2])+'-'+str(entry[7][0][3])+\
                                     ';blk_2_query_len='+str(entry[7][0][1]-entry[7][0][0]+1)+\
                                     ';blk_2_ref_len='+str(entry[7][0][3]-entry[7][0][2]+1)+\
                                     ';blk_2_st_query='+str(entry[7][0][5][0])+';blk_2_st_ref='+str(entry[7][0][5][1])+\
                                     ';blk_2_end_query='+str(entry[7][0][6][0])+';blk_2_end_ref='+str(entry[7][0][6][1])+\
                                     ';color=#A0A0A0'+'\n')


                           
                        elif entry[3]=='misjoin-overlap':
                             ID_name=entry[10]
                            
                             fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';overlap_len='+str(entry[4])+\
                                     ';ref_sequence='+entry[6][0][4]+';blk_1_query='+str(entry[6][0][0])+'-'+str(entry[6][0][1])+\
                                     ';blk_1_ref='+str(entry[6][0][2])+'-'+str(entry[6][0][3])+\
                                     ';blk_1_query_len='+str(entry[6][0][1]-entry[6][0][0]+1)+\
                                     ';blk_1_ref_len='+str(entry[6][0][3]-entry[6][0][2]+1)+\
                                     ';blk_1_st_query='+str(entry[6][0][5][0])+';blk_1_st_ref='+str(entry[6][0][5][1])+\
                                     ';blk_1_end_query='+str(entry[6][0][6][0])+';blk_1_end_ref='+str(entry[6][0][6][1])+\
                                     ';blk_2_query='+str(entry[7][0][0])+'-'+str(entry[7][0][1])+\
                                     ';blk_2_ref='+str(entry[7][0][2])+'-'+str(entry[7][0][3])+\
                                     ';blk_2_query_len='+str(entry[7][0][1]-entry[7][0][0]+1)+\
                                     ';blk_2_ref_len='+str(entry[7][0][3]-entry[7][0][2]+1)+\
                                     ';blk_2_st_query='+str(entry[7][0][5][0])+';blk_2_st_ref='+str(entry[7][0][5][1])+\
                                     ';blk_2_end_query='+str(entry[7][0][6][0])+';blk_2_end_ref='+str(entry[7][0][6][1])+\
                                     ';color=#990099'+'\n')


                        elif entry[3]=='translocation-overlap':
                            ID_name=entry[10]
                            
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0001873'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';overlap_len='+str(entry[4])+\
                                     ';ref_sequence_1='+entry[6][0][4]+\
                                     ';blk_1_query='+str(entry[6][0][0])+'-'+str(entry[6][0][1])+\
                                     ';blk_1_ref='+str(entry[6][0][2])+'-'+str(entry[6][0][3])+\
                                     ';blk_1_query_len='+str(entry[6][0][1]-entry[6][0][0]+1)+\
                                     ';blk_1_ref_len='+str(entry[6][0][3]-entry[6][0][2]+1)+\
                                     ';blk_1_st_query='+str(entry[6][0][5][0])+';blk_1_st_ref='+str(entry[6][0][5][1])+\
                                     ';blk_1_end_query='+str(entry[6][0][6][0])+';blk_1_end_ref='+str(entry[6][0][6][1])+\
                                     ';ref_sequence_2='+entry[7][0][4]+\
                                     ';blk_2_query='+str(entry[7][0][0])+'-'+str(entry[7][0][1])+\
                                     ';blk_2_ref='+str(entry[7][0][2])+'-'+str(entry[7][0][3])+\
                                     ';blk_2_query_len='+str(entry[7][0][1]-entry[7][0][0]+1)+\
                                     ';blk_2_ref_len='+str(entry[7][0][3]-entry[7][0][2]+1)+\
                                     ';blk_2_st_query='+str(entry[7][0][5][0])+';blk_2_st_ref='+str(entry[7][0][5][1])+\
                                     ';blk_2_end_query='+str(entry[7][0][6][0])+';blk_2_end_ref='+str(entry[7][0][6][1])+\
                                     ';color=#A0A0A0'+'\n')


                           
                        
                       
                            

    fq.close()
    
    #reference
    ref_snp_err_dict={}
    ref_struct_err_dict={}
    
    for ref_name in list(ref_dict.keys()):
        ref_snp_err_dict[ref_name]=[]
        ref_struct_err_dict[ref_name]=[]
        

    for cont_name in list(snp_err_dict.keys()):
        for entry in snp_err_dict[cont_name]:
               ref_snp_err_dict[entry[0]].append(entry)

    for cont_name in list(struct_err_dict.keys()):
        for entry in struct_err_dict[cont_name]:
           
            if entry[8]!='a':
                if entry[8]!='struct2':
                    ref_struct_err_dict[entry[0]].append(entry)

               
                else:
                    #if entry[3]=='misjoin' or entry[3]=='circular_genome_start' or entry[3]=='translocation':
                    if entry[3]=='misjoin-overlap' or entry[3]=='translocation-overlap':
                        ref_name=entry[6][0][4]

                        ref_struct_err_dict[ref_name].append([ref_name,entry[6][0][6][1],entry[6][0][6][1],entry[3],entry[6][0][6][2],entry[0],entry[6][0][6][0],
                                                              [entry[1],entry[2],entry[4],'breakpoint'],
                                                              [entry[6][0][0],entry[6][0][1],entry[6][0][2],entry[6][0][3],'block_coord'],
                                                              [entry[6][0][5][0],entry[6][0][5][1],entry[6][0][5][2],'block_start'],entry[10]+'.1'])
                        ref_name=entry[7][0][4]

                        ref_struct_err_dict[ref_name].append([ref_name,entry[7][0][5][1],entry[7][0][5][1],entry[3],entry[7][0][5][2],entry[0],entry[7][0][5][0],
                                                              [entry[1],entry[2],entry[4],'breakpoint'],
                                                              [entry[7][0][0],entry[7][0][1],entry[7][0][2],entry[7][0][3],'block_coord'],
                                                              [entry[7][0][6][0],entry[7][0][6][1],entry[7][0][6][2],'block_end'],entry[10]+'.2'])
                        

                    else:
                        ref_name=entry[6][0][4]

                        ref_struct_err_dict[ref_name].append([ref_name,entry[6][0][6][1],entry[6][0][6][1],entry[3],entry[6][0][6][2],entry[0],entry[6][0][6][0],
                                                              [entry[1],entry[2],entry[4],'breakpoint'],
                                                              [entry[6][0][0],entry[6][0][1],entry[6][0][2],entry[6][0][3],'block_coord'],
                                                              [entry[6][0][5][0],entry[6][0][5][1],entry[6][0][5][2],'block_start'],entry[9]+'.1'])
                        ref_name=entry[7][0][4]

                        ref_struct_err_dict[ref_name].append([ref_name,entry[7][0][5][1],entry[7][0][5][1],entry[3],entry[7][0][5][2],entry[0],entry[7][0][5][0],
                                                              [entry[1],entry[2],entry[4],'breakpoint'],
                                                              [entry[7][0][0],entry[7][0][1],entry[7][0][2],entry[7][0][3],'block_coord'],
                                                              [entry[7][0][6][0],entry[7][0][6][1],entry[7][0][6][2],'block_end'],entry[9]+'.2'])
                        
                    
                        
    for entry in err_ref_cont_coord_errors_list:
        
        if entry[8]=='r':
                
                if entry[6]!='contig_st' and entry[6]!='contig_end' and entry[6]!='uncovered_region':
                    ref_name=entry[0]
                    ref_struct_err_dict[ref_name].append(entry)
                    ref_struct_err_dict[ref_name][-1].append('SV_'+str(ID_cur))
                    ID_cur+=1
                    


    for r_name in list(ref_snp_err_dict.keys()):
        ref_snp_err_dict[r_name]= sorted(ref_snp_err_dict[r_name],key=lambda inter:inter[1], reverse=False)
            
    for r_name in list(ref_struct_err_dict.keys()):
        ref_struct_err_dict[r_name]= sorted(ref_struct_err_dict[r_name],key=lambda inter:inter[1], reverse=False)

    if vcf_flag=='yes':
            OUTPUT_VCF_FILE_REF(ref_name_full,ref_full_names_dict, ref_names,ref_snp_err_dict,working_dir+prefix+'_ref_snps.vcf',ref_dict,cont_dict)  

    #ref snps
    fr=open(working_dir+prefix+'_ref_snps.gff','w')
    fr.write('##gff-version 3\n')

    for ref_name in ref_names:
        if ref_name_full=='yes':
            r_name=ref_full_names_dict[ref_name]
        else:
            r_name=ref_name
                

        if ref_snp_err_dict[ref_name]!=[]:
                fr.write('##sequence-region\t'+r_name+'\t1\t'+str(len(ref_dict[ref_name]))+'\n')

               
                for entry in ref_snp_err_dict[ref_name]:
                        if asmb_name_full=='yes':
                            c_name=contig_full_names_dict[entry[3]]
                        else:
                            c_name=entry[3]

                            
                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap':
                            if entry[9]==1:
                                 query_coord=str(entry[4])+'-'+str(entry[5])
                                 query_bases=cont_dict[entry[3]][entry[4]-1:entry[5]]
                            else:
                                 query_coord=str(entry[5])+'-'+str(entry[4])
                                 query_bases=general.COMPL_STRING(cont_dict[entry[3]][entry[4]-1:entry[5]])
                         

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000667'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+query_coord+';query_bases='+query_bases+\
                                        ';ref_bases=-'+';color=#EE0000'+'\n')
                        elif entry[6]=='deletion':
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+';query_bases=-;ref_bases='+\
                                        ref_dict[entry[0]][entry[1]-1:entry[2]]+';color=#0000EE'+'\n')
                        
                        else: #subst and gap
                            if entry[9]==1:
                                 query_coord=str(entry[4])+'-'+str(entry[5])
                                 query_bases=cont_dict[entry[3]][entry[4]-1:entry[5]]
                            else:
                                 query_coord=str(entry[5])+'-'+str(entry[4])
                                 query_bases=general.COMPL_STRING(cont_dict[entry[3]][entry[4]-1:entry[5]])
                         
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';subst_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+query_coord+\
                                     ';query_bases='+query_bases+';ref_bases='+ref_dict[entry[0]][entry[1]-1:entry[2]]+';color=#42C042'+'\n')
                         
                        if entry[4]>entry[5] or entry[7]<0:
                            print(entry)
                            print('ERROR: a wrong query coordinate')
    
    
    fr.close()


    #ref_struct.gff        
    fr=open(working_dir+prefix+'_ref_struct.gff','w')
    
    fr.write('##gff-version 3\n')
    

    
    for ref_name  in ref_names:
        if ref_name_full=='yes':
            r_name=ref_full_names_dict[ref_name]
        else:
            r_name=ref_name
                

        if ref_struct_err_dict[ref_name]!=[]:
                fr.write('##sequence-region\t'+r_name+'\t1\t'+str(len(ref_dict[ref_name]))+'\n')

                for entry in ref_struct_err_dict[ref_name]:
                    

                    if entry[8]=='struct':
                        
                        if asmb_name_full=='yes':
                            c_name=contig_full_names_dict[entry[3]]
                        else:
                            c_name=entry[3]
                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap' or entry[6]=='wrong_end' or entry[6]=='wrong_beginning' or entry[6]=='wrong_gap-overestimated':

                            if len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'
                            elif len(entry)==12:
                                ID_name=entry[11]
                                query_dir=str(entry[9])
                            elif len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                            
                            
                            
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000667'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+query_dir+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#EE0000'+'\n')
                            
                            
                            
                                    
                        elif entry[6]=='insertion-multiple_copy':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    #repeat_region_dict[cont_name].append([entry])
                                    
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+\
                                        ';ref_repeated_region='+str(entry[11][0][1])+'-'+str(entry[11][0][2])+';color=#EE0000'+'\n')

                                    
                                elif entry[11]==[]:
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#EE0000'+'\n')

                                    
                            elif len(entry)==12:
                                ID_name=entry[11]
                            
                                fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir=1;query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#EE0000'+'\n')

                                

                                

                        elif entry[6]=='insertion-tandem_multiple_copy':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                   
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+\
                                        ';ref_repeated_region='+str(entry[11][0][1])+'-'+str(entry[11][0][2])+';color=#EE0000'+'\n')

                                elif len(entry[11])==2:
                                                                       
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+\
                                        ';ref_repeated_region_1='+str(entry[11][0][1])+'-'+str(entry[11][0][2])+\
                                        ';ref_repeated_region_2='+str(entry[11][1][1])+'-'+str(entry[11][1][2])+';color=#EE0000'+'\n')


                                    
                                elif entry[11]==[]:
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#EE0000'+'\n')
    
                            elif len(entry)==12:
                                ID_name=entry[11]
                                fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir=1;query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#EE0000'+'\n')
                            
                            
                                
                        elif entry[6]=='deletion' or entry[6]=='deletion-gap_underestimated':

                            if len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'
                            elif len(entry)==12:
                                ID_name=entry[11]
                                query_dir=str(entry[9])
                            elif len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                            
                            
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+query_dir+';query_sequence='+c_name+';query_coord='+str(entry[4])+';color=#0000EE'+'\n')

                           
                           
                        elif entry[6]=='deletion-collapsed_repeat':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                   
                                    
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+\
                                        ';ref_repeated_region='+str(entry[11][0][1])+'-'+str(entry[11][0][2])+';color=#0000EE'+'\n')

                                    
                                elif entry[11]==[]:
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+';color=#0000EE'+'\n')

                                    
                                    
                            elif len(entry)==12:
                                ID_name=entry[11]

                                fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir=1;query_sequence='+c_name+';query_coord='+str(entry[4])+';color=#0000EE'+'\n')
                                
                        

                                
                        elif entry[6]=='deletion-collapsed_tandem_repeat':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    
                                    
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+\
                                        ';ref_repeated_region='+str(entry[11][0][1])+'-'+str(entry[11][0][2])+';color=#0000EE'+'\n')

                                    

                                elif len(entry[11])==2:
                                    
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+\
                                        ';ref_repeated_region_1='+str(entry[11][0][1])+'-'+str(entry[11][0][2])+\
                                        ';ref_repeated_region_2='+str(entry[11][1][1])+'-'+str(entry[11][1][2])+';color=#0000EE'+'\n')

                                    
                                elif entry[11]==[]:
                                    fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+';color=#0000EE'+'\n')
    
                            elif len(entry)==12:
                                ID_name=entry[11]
                                
                                fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir=1;query_sequence='+c_name+';query_coord='+str(entry[4])+';color=#0000EE'+'\n')
                            
                                
                        elif entry[6]=='substitution' or entry[6]=='gap':
                            
                            if len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'
                            elif len(entry)==12:
                                ID_name=entry[11]
                                query_dir=str(entry[9])
                            elif len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                        
                                                        
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';subst_len='+str(entry[7])+';query_dir='+query_dir+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#42C042'+'\n')

                            

                            
                        elif entry[6].startswith('reshuf'):
                            ID_name=entry[11]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+entry[6]+\
                                        ';blk_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#404040'+'\n')
                            
                            
                        elif entry[6]=='inversion':
                            ID_name=entry[11]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';blk_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+'-'+str(entry[5])+';color=#01DFD7'+'\n')

                            
                            '''        
                    elif entry[8]=='r':

                        
                        
                        if entry[6]=='uncovered_region' :
                            ID_name=entry[9]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';region_len='+str(entry[2]-entry[1]+1)+';color=#990000'+'\n')

                            '''   
                            

                            
                    else:# entry[8]=='struct2'
                        if entry[3]=='circular_genome_start' or entry[3]=='misjoin':
                            ID_name=entry[10]


                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';query_sequence='+str(entry[5])+';query_coord='+str(entry[6])+\
                                     ';breakpoint_query='+str(entry[7][0])+'-'+str(entry[7][1])+\
                                     ';blk_query='+str(entry[8][0])+'-'+str(entry[8][1])+\
                                     ';blk_ref='+str(entry[8][2])+'-'+str(entry[8][3])+\
                                     ';blk_query_len='+str(entry[8][1]-entry[8][0]+1)+\
                                     ';blk_ref_len='+str(entry[8][3]-entry[8][2]+1)+\
                                     ';color=#990099'+'\n'
                                     )
                            
                        elif entry[3]=='translocation':
                            ID_name=entry[10]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001873'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';query_sequence='+str(entry[5])+';query_coord='+str(entry[6])+\
                                     ';breakpoint_query='+str(entry[7][0])+'-'+str(entry[7][1])+\
                                     ';blk_query='+str(entry[8][0])+'-'+str(entry[8][1])+\
                                     ';blk_ref='+str(entry[8][2])+'-'+str(entry[8][3])+\
                                     ';blk_query_len='+str(entry[8][1]-entry[8][0]+1)+\
                                     ';blk_ref_len='+str(entry[8][3]-entry[8][2]+1)+\
                                     ';color=#A0A0A0'+'\n'
                                     )
                            

                            
                        elif entry[3]=='misjoin-insertion' or entry[3]=='misjoin-mixed_fragments' or entry[3]=='misjoin-wrong_scaffolding':
                            
                            
                            ID_name=entry[10]

                            
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';ins_len='+str(entry[7][2])+\
                                     ';query_sequence='+str(entry[5])+';query_coord='+str(entry[6])+\
                                     ';breakpoint_query='+str(entry[7][0])+'-'+str(entry[7][1])+\
                                     ';blk_query='+str(entry[8][0])+'-'+str(entry[8][1])+\
                                     ';blk_ref='+str(entry[8][2])+'-'+str(entry[8][3])+\
                                     ';blk_query_len='+str(entry[8][1]-entry[8][0]+1)+\
                                     ';blk_ref_len='+str(entry[8][3]-entry[8][2]+1)+\
                                     ';color=#990099'+'\n'
                                     )

                            
                        elif entry[3]=='translocation-insertion' or entry[3]=='translocation-mixed_fragments' or entry[3]=='translocation-wrong_scaffolding':
                            ID_name=entry[10]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001873'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';ins_len='+str(entry[7][2])+\
                                     ';query_sequence='+str(entry[5])+';query_coord='+str(entry[6])+\
                                     ';breakpoint_query='+str(entry[7][0])+'-'+str(entry[7][1])+\
                                     ';blk_query='+str(entry[8][0])+'-'+str(entry[8][1])+\
                                     ';blk_ref='+str(entry[8][2])+'-'+str(entry[8][3])+\
                                     ';blk_query_len='+str(entry[8][1]-entry[8][0]+1)+\
                                     ';blk_ref_len='+str(entry[8][3]-entry[8][2]+1)+\
                                     ';color=#A0A0A0'+'\n'
                                     )


                                                        
                        elif entry[3]=='misjoin-overlap':
                            
                            
                            ID_name=entry[10]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001874'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';overlap_len='+str(entry[7][2])+\
                                     ';query_sequence='+str(entry[5])+';query_coord='+str(entry[6])+\
                                     ';breakpoint_query='+str(entry[7][0])+'-'+str(entry[7][1])+\
                                     ';blk_query='+str(entry[8][0])+'-'+str(entry[8][1])+\
                                     ';blk_ref='+str(entry[8][2])+'-'+str(entry[8][3])+\
                                     ';blk_query_len='+str(entry[8][1]-entry[8][0]+1)+\
                                     ';blk_ref_len='+str(entry[8][3]-entry[8][2]+1)+\
                                     ';color=#990099'+'\n'
                                     )

                           
                            
                        elif entry[3]=='translocation-overlap':
                            ID_name=entry[10]

                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0001873'+'\t'+str(entry[1])+'\t'+str(entry[2])+\
                                     '\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[3]]+\
                                     ';overlap_len='+str(entry[7][2])+\
                                     ';query_sequence='+str(entry[5])+';query_coord='+str(entry[6])+\
                                     ';breakpoint_query='+str(entry[7][0])+'-'+str(entry[7][1])+\
                                     ';blk_query='+str(entry[8][0])+'-'+str(entry[8][1])+\
                                     ';blk_ref='+str(entry[8][2])+'-'+str(entry[8][3])+\
                                     ';blk_query_len='+str(entry[8][1]-entry[8][0]+1)+\
                                     ';blk_ref_len='+str(entry[8][3]-entry[8][2]+1)+\
                                     ';color=#A0A0A0'+'\n'
                                     )


                           
                       
                                

    fr.close()
    

    #additional files

    
    ref_duplication_dict={}
    for ref_name in list(mapped_blocks_dict.keys()):
        ref_duplication_list=FIND_REF_DUPLICATION(mapped_blocks_dict[ref_name])
        ref_duplication_dict[ref_name]=ref_duplication_list

       
    for entry in err_ref_cont_coord_errors_list:
        if entry[6]=='uncovered_region':
            ref_name=entry[0]

            if ref_name not in ref_duplication_dict:
                ref_duplication_dict[ref_name]=[]

            ref_duplication_dict[ref_name].append([entry[1],entry[2],'Uncovered_region'])    

             

    cont_repeat_region_dict={}
    for cont_name in list(repeat_region_dict.keys()):
        cont_repeat_region_dict[cont_name]=[]
        for entry in repeat_region_dict[cont_name]:
            entry=entry[0]
            
            ref_name=entry[0]
            cont_repeat_region_dict[cont_name].append([entry[11][0][4],entry[11][0][5],'Repeated_region',entry[11][0][6],[entry[4],entry[5],entry[6],entry[7]]])

            if ref_name not in ref_duplication_dict:
                ref_duplication_dict[ref_name]=[]

            ref_duplication_dict[ref_name].append([entry[11][0][1],entry[11][0][2],'Repeated_region',entry[11][0][6],[entry[1],entry[2],entry[6],entry[7]]])

            if len(entry[11])==2:
                cont_repeat_region_dict[cont_name].append([entry[11][1][4],entry[11][1][5],'Repeated_region',entry[11][1][6],[entry[4],entry[5],entry[6],entry[7]]])
                ref_duplication_dict[ref_name].append([entry[11][1][1],entry[11][1][2],'Repeated_region',entry[11][1][6],[entry[1],entry[2],entry[6],entry[7]]])

                

            
    for c_name in list(cont_repeat_region_dict.keys()):
        cont_repeat_region_dict[c_name]= sorted(cont_repeat_region_dict[c_name],key=lambda inter:inter[0], reverse=False)

    for r_name in list(ref_duplication_dict.keys()):
        ref_duplication_dict[r_name]= sorted(ref_duplication_dict[r_name],key=lambda inter:inter[0], reverse=False)
            
    
    for entry in err_ref_cont_coord_errors_list:
        if entry[8]=='struct2':
            if entry[3]=='misjoin-overlap': 
                ref_name=entry[9][0][0]
                if ref_name not in ref_duplication_dict:
                    ref_duplication_dict[ref_name]=[]

                ref_duplication_dict[ref_name].append([entry[9][0][1],entry[9][0][2],'Relocation_overlap_region'])

                ref_name=entry[9][1][0]
                if ref_name not in ref_duplication_dict:
                    ref_duplication_dict[ref_name]=[]

                ref_duplication_dict[ref_name].append([entry[9][1][1],entry[9][1][2],'Relocation_overlap_region'])
    
            
            elif entry[3]=='translocation-overlap': 
                ref_name=entry[9][0][0]
                if ref_name not in ref_duplication_dict:
                    ref_duplication_dict[ref_name]=[]

                ref_duplication_dict[ref_name].append([entry[9][0][1],entry[9][0][2],'Translocation_overlap_region'])

                ref_name=entry[9][1][0]
                if ref_name not in ref_duplication_dict:
                    ref_duplication_dict[ref_name]=[]

                ref_duplication_dict[ref_name].append([entry[9][1][1],entry[9][1][2],'Translocation_overlap_region'])
    
            

    for ref_name in list(ref_duplication_dict.keys()):
        ref_duplication_dict[ref_name]=sorted(ref_duplication_dict[ref_name],key=lambda inter:inter[0], reverse=False)    

    f=open(working_dir+prefix+'_ref_additional.gff','w')
    
    f.write('##gff-version 3\n')

    ID=1
    for r_name in ref_names:
        if ref_name_full=='yes':
            ref_name=ref_full_names_dict[r_name]
        else:
            ref_name=r_name
        

        if r_name in ref_duplication_dict:
            if ref_duplication_dict[r_name]!=[]:
                    f.write('##sequence-region\t'+ref_name+'\t1\t'+str(len(ref_dict[r_name]))+'\n')

                    for entry in ref_duplication_dict[r_name]:
                        
                        ID_cur='Region_'+str(ID)
                        ID+=1

                        entry.append(ID_cur)

                        if entry[2]=='Repeated_region':
                            f.write(ref_name+'\tNucDiff_v2.0\t'+'SO:0000657'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t.\t.\tID='+entry[5]+';Name='+entry[2]+\
                                ';repeat_len='+str(entry[1]-entry[0]+1)+';difference_type='+err_new_names_dict[entry[4][2]]+';difference_coord_ref='+str(entry[4][0])+'-'+str(entry[4][1])+\
                                    ';difference_len='+str(entry[4][3])+';color=#DB0101'+'\n')

                        elif entry[2]=='Relocation_overlap_region':
                            f.write(ref_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t.\t.\tID='+entry[3]+';Name='+entry[2]+\
                                ';overlap_len='+str(entry[1]-entry[0]+1)+';color=#00A123'+'\n')

                        elif  entry[2]=='Translocation_overlap_region':
                            f.write(ref_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t.\t.\tID='+entry[3]+';Name='+entry[2]+\
                                ';overlap_len='+str(entry[1]-entry[0]+1)+';color=#8E00A1'+'\n')


                        elif entry[2]=='Uncovered_region':  
                            f.write(ref_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t.\t.\tID='+entry[3]+';Name='+entry[2]+\
                                        ';region_len='+str(entry[1]-entry[0]+1)+';color=#000000'+'\n')
   
                        else:
                            
                            
                            f.write(ref_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t.\t.\tID='+entry[3]+';Name='+'Ref_duplication'+\
                                ';duplic_len='+str(entry[1]-entry[0]+1)+';color=#4005BF'+'\n')
                        
    f.close()


    f=open(working_dir+prefix+'_query_additional.gff','w')
    
    f.write('##gff-version 3\n')

    ID=1
    for c_name in contig_names:
        if asmb_name_full=='yes':
            asmb_name=contig_full_names_dict[c_name]
        else:
            asmb_name=c_name

        if c_name in cont_repeat_region_dict:

            if cont_repeat_region_dict[c_name]!=[]:

                f.write('##sequence-region\t'+asmb_name+'\t1\t'+str(len(cont_dict[c_name]))+'\n')

                for entry in cont_repeat_region_dict[c_name]:
                    entry.append('Region_'+str(ID))
                    ID+=1

                    f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000657'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[5]+';Name='+entry[2]+\
                                ';repeat_len='+str(entry[1]-entry[0]+1)+';difference_type='+err_new_names_dict[entry[4][2]]+';difference_coord_query='+str(entry[4][0])+'-'+str(entry[4][1])+\
                                    ';difference_len='+str(entry[4][3])+'\n')


                        
   
                    
                        
                        
    f.close()
    
    
    

def OUTPUT_STAT(statistics_output_lines,working_dir,prefix,cont_num,ref_num):
    f=open(working_dir+prefix+'_stat.out','w')
    for entry in statistics_output_lines:
        f.write(entry)
    f.write('\n\n')
    f.write('ADDITIONAL INFORMATION:\t\n')
    f.write('query sequences\t'+str(cont_num)+'\n')
    f.write('reference sequences\t'+str(ref_num)+'\n')
    f.close()


def FIND_REF_DUPLICATION(input_list):
    ref_duplication_list=[]


    

    temp_dict={}
    for entry in input_list:
        if entry[0] not in temp_dict:
            temp_dict[entry[0]]={'st':0,'end':0}
        temp_dict[entry[0]]['st']+=1

        if entry[1]+1 not in temp_dict:
            temp_dict[entry[1]+1]={'st':0,'end':0}
        temp_dict[entry[1]+1]['end']+=1

    cur_num=0
    flag=0
    for pos in sorted(temp_dict.keys()):
        cur_num=cur_num+temp_dict[pos]['st']-temp_dict[pos]['end']

        if flag==0:
            if cur_num>1:
                int_st=pos
                flag=1

        else:
            
            if cur_num<2:
                ref_duplication_list.append([int_st,pos-1,'Ref_duplication'])
                flag=0

            
  
        
    return ref_duplication_list       
def OUTPUT_MAPPED_BLOCKS_TO_REF(struct_dict,ref_dict,ref_names,ref_full_names_dict,cont_dict,contig_full_names_dict,working_dir, prefix,asmb_name_full,ref_name_full):

    

    mapped_blocks_dict={}
    for r_name in list(ref_dict.keys()):
        mapped_blocks_dict[r_name]=[]

   
    for cont in list(struct_dict.keys()):
        for trl in list(struct_dict[cont].keys()):
            for msj in list(struct_dict[cont][trl]['blocks'].keys()):
                for bl in list(struct_dict[cont][trl]['blocks'][msj]['blocks'].keys()):
                    blk_coord=struct_dict[cont][trl]['blocks'][msj]['blocks'][bl]['block']

                    ref_name=blk_coord[8]
                    mapped_blocks_dict[ref_name].append([blk_coord[2],blk_coord[3],cont,blk_coord[0],blk_coord[1], blk_coord[4]])
                    
    for ref in list(mapped_blocks_dict.keys()):
        mapped_blocks_dict[ref]=sorted(mapped_blocks_dict[ref],key=lambda inter:inter[0], reverse=False)


    
    f=open(working_dir+prefix+'_ref_blocks.gff','w')
    
    f.write('##gff-version 3\n')

    ID=1
    for r_name in ref_names:
        if ref_name_full=='yes':
            ref_name=ref_full_names_dict[r_name]
        else:
            ref_name=r_name
        

        if r_name in mapped_blocks_dict:
            if mapped_blocks_dict[r_name]!=[]:
                    f.write('##sequence-region\t'+ref_name+'\t1\t'+str(len(ref_dict[r_name]))+'\n')

                    for entry in mapped_blocks_dict[r_name]:
                        if asmb_name_full=='yes':
                            c_name=contig_full_names_dict[entry[2]]
                        else:
                            c_name=entry[2]

                        
                        ID_cur='Blk_'+str(ID)
                        ID+=1

                        entry.append(ID_cur)

                        f.write(ref_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t.\t.\tID='+entry[6]+';Name='+c_name+'_block'+\
                                ';blk_length='+str(entry[1]-entry[0]+1)+';query_sequence_length='+str(len(cont_dict[entry[2]]))+\
                                ';blk_coord_query='+str(entry[3])+'-'+str(entry[4])+';query_dir='+str(entry[5])+'\n')
                        
    f.close()


    return mapped_blocks_dict


    
def OUTPUT_BLOCKS_TO_QUERY(cont_blocks_dict,ref_dict,contig_names,ref_full_names_dict,cont_dict,contig_full_names_dict,working_dir, prefix,asmb_name_full,ref_name_full):

    f=open(working_dir+prefix+'_query_blocks.gff','w')
    
    f.write('##gff-version 3\n')

    
    

    ID=1
    for c_name in contig_names:
        if asmb_name_full=='yes':
            asmb_name=contig_full_names_dict[c_name]
        else:
            asmb_name=c_name

        if c_name in cont_blocks_dict:

            if cont_blocks_dict[c_name]!=[]:

                f.write('##sequence-region\t'+asmb_name+'\t1\t'+str(len(cont_dict[c_name]))+'\n')

                for entry in cont_blocks_dict[c_name]:
                    entry.append('Blk_'+str(ID))
                    ID+=1

                    if ref_name_full=='yes':
                            if not entry[2] in ['Circular_genome_block','Relocation_block','Translocation_block']:
                                
                                rf_name=ref_full_names_dict[entry[4]]
                            else:
                                rf_name=entry[4]
                    else:
                            rf_name=entry[4]

                    if entry[2]=='Inversion': 
                        f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[8]+';Name='+entry[2]+\
                                ';blk_length='+str(entry[3])+';query_dir='+str(entry[7])+';ref_sequence='+rf_name+';ref_coord='+str(entry[5])+'-'+str(entry[6])+';color=#DF0101'+'\n')

                    elif entry[2]=='Block':
                        f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[8]+';Name='+entry[2]+\
                                ';blk_length='+str(entry[3])+';query_dir='+str(entry[7])+';ref_sequence='+rf_name+';ref_coord='+str(entry[5])+'-'+str(entry[6])+';color=#000000'+'\n')

                    elif entry[2].startswith('Resh'):
                         f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[8]+';Name='+entry[2]+\
                                ';blk_length='+str(entry[3])+';query_dir='+str(entry[7])+';ref_sequence='+rf_name+';ref_coord='+str(entry[5])+'-'+str(entry[6])+';color=#04B404'+'\n')

                    elif entry[2]=='Relocation_block':
                        f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[4]+';Name='+entry[2]+\
                                ';blk_length='+str(entry[3])+';color=#01DFD7'+'\n')

                    elif entry[2]=='Translocation_block':
                        f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[4]+';Name='+entry[2]+\
                                ';blk_length='+str(entry[3])+';color=#0404B4'+'\n')

                    elif entry[2]=='Circular_genome_block':
                        f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[4]+';Name='+entry[2]+\
                                ';blk_length='+str(entry[3])+';color=#990099'+'\n')
   
                    else:  
                        if len(entry)==5:
                            f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[4]+';Name='+entry[2]+\
                                    ';blk_length='+str(entry[3])+'\n')

                            

                        elif len(entry)==9:
                            f.write(asmb_name+'\tNucDiff_v2.0\t'+'SO:0000001'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tID='+entry[8]+';Name='+entry[2]+\
                                    ';blk_length='+str(entry[3])+';query_dir='+str(entry[7])+';ref_sequence='+rf_name+';ref_coord='+str(entry[5])+'-'+str(entry[6])+'\n')

                        
                    
                        
                        
    f.close()
    

def GENERATE_OUTPUT(struct_dict,end_err_dict,unmapped_list, file_ref, file_contigs,working_dir, prefix,err_ref_cont_coord_errors_list,statistics_output_lines,asmb_name_full,ref_name_full,cont_blocks_dict,vcf_flag):

    contigs_dict, contig_seqs, contig_names, contig_full_names_dict=general.READ_FASTA_ENTRY(file_contigs)
    ref_dict, ref_seqs, ref_names,ref_full_names_dict=general.READ_FASTA_ENTRY(file_ref)

    
    #OUTPUT_READABLE(struct_dict,end_err_dict,unmapped_list, contigs_dict,contig_names, contig_full_names_dict, ref_dict,ref_full_names_dict,working_dir, 'results/'+prefix)

    mapped_blocks_dict=OUTPUT_MAPPED_BLOCKS_TO_REF(struct_dict,ref_dict,ref_names,ref_full_names_dict,contigs_dict,contig_full_names_dict,working_dir, 'results/'+prefix,asmb_name_full,ref_name_full)
    OUTPUT_REF_ASSEM_TABLE(err_ref_cont_coord_errors_list, ref_dict,ref_names,ref_full_names_dict,contigs_dict,contig_names,contig_full_names_dict,working_dir, 'results/'+prefix,asmb_name_full,ref_name_full,mapped_blocks_dict,vcf_flag)
    OUTPUT_STAT(statistics_output_lines,working_dir, 'results/'+prefix, len(list(contigs_dict.keys())),len(list(ref_dict.keys())))
    OUTPUT_BLOCKS_TO_QUERY(cont_blocks_dict, ref_dict,contig_names,ref_full_names_dict,contigs_dict,contig_full_names_dict,working_dir, 'results/'+prefix,asmb_name_full,ref_name_full)
