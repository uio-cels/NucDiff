# copyright (c) 2014-2016 Ksenia Khelik
#
# This Source Code Form is subject to the 
# terms of the Mozilla Public License, v. 2.0. 
# If a copy of the MPL was not distributed 
# with this file, You can obtain one at 
# https://mozilla.org/MPL/2.0/.
#
#-------------------------------------------------------------------------------



import general

err_new_names_dict={'substitution':'substitution', 'gap':'gap',

                      'insertion':'insertion','insertion-multiple_copy':'duplication',
                      'insertion-tandem_multiple_copy':'tandem_duplication','wrong_beginning':'unaligned_beginning',
                      'wrong_end':'unaligned_end', 'wrong_gap':'inserted_gap', 'wrong_gap-overestimated':'inserted_gap',

                      'deletion':'deletion', 'deletion-collapsed_repeat':'collapsed_repeat',
                      'deletion-collapsed_tandem_repeat':'tandem_collapsed_repeat','deletion-gap_underestimated':'deletion',

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

                      'contig_st':'query_st', 'contig_end':'query_end', 'uncovered_region':'uncovered_region'
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

    
        if len(structure_dict_misj.keys())==1: #no transposition
            block_name=structure_dict_misj.keys()[0]
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

    
        if len(structure_dict_transl.keys())==1: #no misjoin
            misj_group=structure_dict_transl.keys()[0]
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

    
        if len(structure_dict_cont.keys())==1: #no translocation
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
            if end_err_dict.has_key(cont_name) and (end_err_dict[cont_name]['wrong_beginning']!=[] or end_err_dict[cont_name]['wrong_end']!=[]):
                f.write(contig_full_names_dict[cont_name]+'\t'+str(len(cont_dict[cont_name]))+'\n' )

                if  end_err_dict[cont_name]['wrong_beginning']!=[]:
                        entry=end_err_dict[cont_name]['wrong_beginning']
                        f.write('|'+str(entry[0])+'\t'+str(entry[1])+'\t'+err_new_names_dict[entry[2]]+'\t'+str(entry[3])+'\n')
                        
                        
                        

                if len(structure_dict[cont_name].keys())==1 and \
                       len(structure_dict[cont_name][0]['blocks'].keys())==1 and \
                       len(structure_dict[cont_name][0]['blocks'][0]['blocks'].keys())==1 and \
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
                    if len(structure_dict[cont_name].keys())==1 and \
                       len(structure_dict[cont_name][0]['blocks'].keys())==1 and \
                       len(structure_dict[cont_name][0]['blocks'][0]['blocks'].keys())==1 and \
                       structure_dict[cont_name][0]['blocks'][0]['blocks'][0]['block'][10]==[]:
                            a='correct_query'
                    else:
                        f.write(contig_full_names_dict[cont_name]+'\t'+str(len(cont_dict[cont_name]))+'\n' )
                        OUTPUT_TRANSLOC_PART(f,structure_dict[cont_name], cont_name)
                        f.write('\n\n')
                        
        f.close()


    



def OUTPUT_REF_ASSEM_TABLE(err_ref_cont_coord_errors_list, ref_dict,ref_names,ref_full_names_dict,cont_dict,contig_names,contig_full_names_dict,working_dir, prefix,asmb_name_full,ref_name_full):

   #query_snps.gff

    
    snp_err_dict={}
    struct_err_dict={}
    repeat_region_dict={}
    for cont_name in cont_dict.keys():
        snp_err_dict[cont_name]=[]
        struct_err_dict[cont_name]=[]
        repeat_region_dict[cont_name]=[]

   
        
    for entry in err_ref_cont_coord_errors_list:
        if entry[8]=='b':
           
            if len(entry)==11 and entry[10]=='snps':
                c_name=entry[3]
                snp_err_dict[c_name].append(entry)
            else:
                c_name=entry[3]
                struct_err_dict[c_name].append(entry)
       # elif entry[8]=='a':
       #     print entry
       #    raw_input('jdf')
       
    
    for c_name in snp_err_dict.keys():
        snp_err_dict[c_name]= sorted(snp_err_dict[c_name],key=lambda inter:inter[4], reverse=False)
            
    for c_name in struct_err_dict.keys():
        struct_err_dict[c_name]= sorted(struct_err_dict[c_name],key=lambda inter:inter[4], reverse=False)
        


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
                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000667'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';query_bases='+cont_dict[cont_name][entry[4]-1:entry[5]]+\
                                        ';ref_bases=-'+';color=EE0000'+'\n')
                        elif entry[6]=='deletion':
                            if entry[9]==1:
                                 ref_coord=str(entry[1])+'-'+str(entry[2])
                                 ref_bases=ref_dict[entry[0]][entry[1]-1:entry[2]]
                            else:
                                 ref_coord=str(entry[2])+'-'+str(entry[1])
                                 ref_bases=general.COMPL_STRING(ref_dict[entry[0]][entry[1]-1:entry[2]])
                         
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+ref_coord+';query_bases=-;ref_bases='+ref_bases+\
                                     ';color=0000EE'+'\n')
                        
                        else: #subst and gap
                            if entry[9]==1:
                                 ref_coord=str(entry[1])+'-'+str(entry[2])
                                 ref_bases=ref_dict[entry[0]][entry[1]-1:entry[2]]
                            else:
                                 ref_coord=str(entry[2])+'-'+str(entry[1])
                                 ref_bases=general.COMPL_STRING(ref_dict[entry[0]][entry[1]-1:entry[2]])
                         
                            fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';subst_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+ref_coord+';query_bases='+cont_dict[cont_name][entry[4]-1:entry[5]]+\
                                     ';ref_bases='+ref_bases+';color=42C042'+'\n')
                         
                        if entry[4]>entry[5] or entry[7]<0:
                            print entry
                            print 'ERROR: a wrong query coordinate'
    

    
    fq.close()


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
                        if ref_name_full=='yes':
                            r_name=ref_full_names_dict[entry[0]]
                        else:
                            r_name=entry[0]
                                
                        if entry[6]=='insertion' or entry[6]=='wrong_gap' or entry[6]=='wrong_end' or entry[6]=='wrong_beginning':
                            if len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                                
                            elif len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'

                            
                            
                            #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000667'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                            #            ';ins_len='+str(entry[7])+';query_dir='+query_dir+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=EE0000'+'\n')
                                
                        elif entry[6]=='insertion-multiple_copy':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                    #    ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+\
                                    #    ';query_repeat_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=EE0000'+'\n')

                                    
                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=EE0000'+'\n')
    
                            elif len(entry)==10:
                                ID_name=entry[9]
                            
                                #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000035'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                #        ';ins_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=EE0000'+'\n')
    

                        elif entry[6]=='insertion-tandem_multiple_copy':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                    #    ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+\
                                    #    ';query_repeat_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=EE0000'+'\n')

                                elif len(entry[11])==2:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                    #    ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+\
                                    #    ';query_repeat_region_1='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+\
                                    #    ';query_repeat_region_2='+str(entry[11][1][4])+'-'+str(entry[11][1][5])+';color=EE0000'+'\n')

                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=EE0000'+'\n')
    
                            elif len(entry)==10:
                                ID_name=entry[9]

                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000173'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';ins_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+';color=EE0000'+'\n')
                            

                        if entry[6]=='deletion':
                            if len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                                
                            elif len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'

                            
                            
                            #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                            #            ';del_len='+str(entry[7])+';query_dir='+query_dir+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=0000EE'+'\n')

                        elif entry[6]=='deletion-collapsed_repeat':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                    #    ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+\
                                    #    ';query_repeat_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=0000EE'+'\n')

                                    
                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=0000EE'+'\n')
    
                            elif len(entry)==10:
                                ID_name=entry[9]

                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=0000EE'+'\n')
    
                        elif entry[6]=='deletion-collapsed_tandem_repeat':
                            if len(entry)==13:
                                ID_name=entry[12]

                                if len(entry[11])==1:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                    #    ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+\
                                    #    ';query_repeat_region='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+';color=0000EE'+'\n')

                                elif len(entry[11])==2:
                                    repeat_region_dict[cont_name].append([entry])
                                    
                                    #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                    #    ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+\
                                    #    ';query_repeat_region_1='+str(entry[11][0][4])+'-'+str(entry[11][0][5])+\
                                    #    ';query_repeat_region_2='+str(entry[11][1][4])+'-'+str(entry[11][1][5])+';color=0000EE'+'\n')

                                    
                                elif entry[11]==[]:
                                    fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=0000EE'+'\n')
    
                            elif len(entry)==10:
                                ID_name=entry[9]

                               

                                fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir=1;ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=0000EE'+'\n')

                        elif entry[6]=='substitution' or entry[6]=='gap':
                            if len(entry)==13:
                                ID_name=entry[12]
                                query_dir=str(entry[9])
                                
                            elif len(entry)==10:
                                ID_name=entry[9]
                                query_dir='1'

                                                        
                            #fq.write(c_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t.\t.\tID='+ID_name+';Name='+err_new_names_dict[entry[6]]+\
                            #            ';subst_len='+str(entry[7])+';query_dir='+query_dir+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+';color=42C042'+'\n')
                                            
                        elif entry[6]=='circular_genome_start':
                            print entry
                            raw_input('hd')

                            #parse afterwards!!!!!
                        #elif entry[6]=='misjoin':
                            

    fq.close()
    

    ref_snp_err_dict={}
    for ref_name in ref_dict.keys():
        ref_snp_err_dict[ref_name]=[]
        

    for cont_name in snp_err_dict.keys():
        for entry in snp_err_dict[cont_name]:
           ref_snp_err_dict[entry[0]].append(entry)

    for r_name in ref_snp_err_dict.keys():
        ref_snp_err_dict[r_name]= sorted(ref_snp_err_dict[r_name],key=lambda inter:inter[1], reverse=False)
            
       


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
                            c_name=ref_full_names_dict[entry[3]]
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
                                        ';ref_bases=-'+';color=EE0000'+'\n')
                        elif entry[6]=='deletion':
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:0000159'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';del_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+str(entry[4])+';query_bases=-;ref_bases='+\
                                        ref_dict[entry[0]][entry[1]-1:entry[2]]+';color=0000EE'+'\n')
                        
                        else: #subst and gap
                            if entry[9]==1:
                                 query_coord=str(entry[4])+'-'+str(entry[5])
                                 query_bases=cont_dict[entry[3]][entry[4]-1:entry[5]]
                            else:
                                 query_coord=str(entry[5])+'-'+str(entry[4])
                                 query_bases=general.COMPL_STRING(cont_dict[entry[3]][entry[4]-1:entry[5]])
                         
                            fr.write(r_name+'\tNucDiff_v2.0\t'+'SO:1000002'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t.\t.\tID='+entry[11]+';Name='+err_new_names_dict[entry[6]]+\
                                        ';subst_len='+str(entry[7])+';query_dir='+str(entry[9])+';query_sequence='+c_name+';query_coord='+query_coord+\
                                     ';query_bases='+query_bases+';ref_bases='+ref_dict[entry[0]][entry[1]-1:entry[2]]+';color=42C042'+'\n')
                         
                        if entry[4]>entry[5] or entry[7]<0:
                            print entry
                            print 'ERROR: a wrong query coordinate'
    
    
    fr.close()
                         
    

    '''
    #copy modified
    f=open(working_dir+prefix+'_query_snps.gff','w')
    
    f.write('##gff-version 3\n')

    for cont_name in contig_names:
        if asmb_name_full=='yes':
            c_name=contig_full_names_dict[cont_name]
        else:
            c_name=cont_name
                

        if cont_err_dict.has_key(cont_name):
                if cont_err_dict[cont_name]!=[]:
                    f.write('##sequence-region\t'+contig_full_names_dict[cont_name]+'\t1\t'+str(len(cont_dict[cont_name]))+'\n')
                    for entry in cont_err_dict[cont_name]:
                        if entry[0]=='-':
                            f.write(contig_full_names_dict[cont_name]+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+'\n')
                        else:
                            if ref_name_full=='yes':
                                r_name=ref_full_names_dict[entry[0]]
                            else:
                                r_name=entry[0]
                                
                            if entry[6].startswith('reshuffling'):
                                f.write(c_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                            else:
                                f.write(c_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';ref_sequence='+r_name+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                
                        if entry[4]>entry[5] or entry[7]<0:
                            print entry
                            print 'ERROR: a wrong query coordinate'
    

    f.close()
    '''

    
    '''
    cont_err_dict={}
    for cont_name in cont_dict.keys():
        cont_err_dict[cont_name]=[]
        
    for entry in err_ref_cont_coord_errors_list:
        if entry[8]=='b' or entry[8]=='a':
            c_name=entry[3]
            cont_err_dict[c_name].append(entry)

    for c_name in cont_err_dict.keys():
        cont_err_dict[c_name]= sorted(cont_err_dict[c_name],key=lambda inter:inter[4], reverse=False)
            
    

    
    f=open(working_dir+prefix+'_query_coord.gff','w')
    
    f.write('##gff-version 3\n')

    if asmb_name_full=='yes':
         for cont_name in contig_names:
                if cont_err_dict.has_key(cont_name):
            
                    if cont_err_dict[cont_name]!=[]:
                        f.write('##sequence-region\t'+contig_full_names_dict[cont_name]+'\t1\t'+str(len(cont_dict[cont_name]))+'\n')
                        for entry in cont_err_dict[cont_name]:
                            if entry[0]=='-':
                                f.write(contig_full_names_dict[cont_name]+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+'\n')
                            else:
                                if ref_name_full=='yes':
                                    if entry[6].startswith('reshuffling'):
                                        f.write(contig_full_names_dict[cont_name]+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';ref_sequence='+ref_full_names_dict[entry[0]]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                    else:
                                        f.write(contig_full_names_dict[cont_name]+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';ref_sequence='+ref_full_names_dict[entry[0]]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                else:
                                    if entry[6].startswith('reshuffling'):    
                                        f.write(contig_full_names_dict[cont_name]+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';ref_sequence='+entry[0]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                    else:
                                        f.write(contig_full_names_dict[cont_name]+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';ref_sequence='+entry[0]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')

                            if entry[4]>entry[5] or entry[7]<0:
                                print entry
                                print 'ERROR: a wrong query coordinate'
    else:
          for cont_name in contig_names:
                if cont_err_dict.has_key(cont_name):
            
                    if cont_err_dict[cont_name]!=[]:
                        f.write('##sequence-region\t'+cont_name+'\t1\t'+str(len(cont_dict[cont_name]))+'\n')
                        for entry in cont_err_dict[cont_name]:
                            if entry[0]=='-':
                                f.write(cont_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+'\n')
                            else:
                                if ref_name_full=='yes':
                                    if entry[6].startswith('reshuffling'):
                                        f.write(cont_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';ref_sequence='+ref_full_names_dict[entry[0]]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                    else:
                                        f.write(cont_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';ref_sequence='+ref_full_names_dict[entry[0]]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                else:
                                    if entry[6].startswith('reshuffling'):
                                        f.write(cont_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';ref_sequence='+entry[0]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')
                                    else:    
                                        f.write(cont_name+'\t.\t'+'Difference'+'\t'+str(entry[4])+'\t'+str(entry[5])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';ref_sequence='+entry[0]+';ref_coord='+str(entry[1])+'-'+str(entry[2])+'\n')

                            if entry[4]>entry[5] or entry[7]<0:
                                print entry
                                print 'ERROR: a wrong query coordinate'


    f.close()
    


   
    

    ref_err_dict={}
    for r_name in ref_dict.keys():
        ref_err_dict[r_name]=[]
        
    for entry in err_ref_cont_coord_errors_list:
        if entry[8]=='b' or entry[8]=='r':
            r_name=entry[0]
            ref_err_dict[r_name].append(entry)

    for r_name in ref_err_dict.keys():
        ref_err_dict[r_name]= sorted(ref_err_dict[r_name],key=lambda inter:inter[1], reverse=False)
            
    
    
    f=open(working_dir+prefix+'_ref_coord.gff','w')
    
    f.write('##gff-version 3\n')

    if ref_name_full=='yes':

        for r_name in ref_names:
            if ref_err_dict.has_key(r_name):

        
                if ref_err_dict[r_name]!=[]:
                    f.write('##sequence-region\t'+ref_full_names_dict[r_name]+'\t1\t'+str(len(ref_dict[r_name]))+'\n')
                    for entry in ref_err_dict[r_name]:
                        if entry[3]=='-':
                            f.write(ref_full_names_dict[r_name]+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+'\n')
                        else:
                            if asmb_name_full=='yes':
                                if entry[6].startswith('reshuffling'):
                                        f.write(ref_full_names_dict[r_name]+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';query_sequence='+contig_full_names_dict[entry[3]]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')
                                else:        
                                        f.write(ref_full_names_dict[r_name]+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';query_sequence='+contig_full_names_dict[entry[3]]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')
                            else:
                                if entry[6].startswith('reshuffling'):
                                        f.write(ref_full_names_dict[r_name]+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';query_sequence='+entry[3]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')
                                else:
                                        f.write(ref_full_names_dict[r_name]+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';query_sequence='+entry[3]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')

                        if entry[1]>entry[2] or entry[7]<0:
                            print entry
                            print 'ERROR: a wrong reference coordinate'
    else:
        for r_name in ref_names:
            if ref_err_dict.has_key(r_name):

        
                if ref_err_dict[r_name]!=[]:
                    f.write('##sequence-region\t'+r_name+'\t1\t'+str(len(ref_dict[r_name]))+'\n')
                    for entry in ref_err_dict[r_name]:
                        if entry[3]=='-':
                            f.write(r_name+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+'\n')
                        else:
                            if asmb_name_full=='yes':
                                if entry[6].startswith('reshuffling'):
                                        f.write(r_name+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';query_sequence='+contig_full_names_dict[entry[3]]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')
                                else:
                                        f.write(r_name+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';query_sequence='+contig_full_names_dict[entry[3]]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')
                            else:
                                if entry[6].startswith('reshuffling'):
                                     f.write(r_name+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+entry[6]+';length='+str(entry[7])+';query_sequence='+entry[3]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')   
                                else:
                                        f.write(r_name+'\t.\t'+'Difference'+'\t'+str(entry[1])+'\t'+str(entry[2])+'\t.\t+\t.\tName='+err_new_names_dict[entry[6]]+';length='+str(entry[7])+';query_sequence='+entry[3]+';query_coord='+str(entry[4])+'-'+str(entry[5])+'\n')

                        if entry[1]>entry[2] or entry[7]<0:
                            print entry
                            print 'ERROR: a wrong reference coordinate'
                        

    f.close()
    '''
    
def OUTPUT_STAT(statistics_output_lines,working_dir,prefix,cont_num,ref_num):
    f=open(working_dir+prefix+'_stat.out','w')
    for entry in statistics_output_lines:
        f.write(entry)
    f.write('\n\n')
    f.write('ADDITIONAL INFORMATION:\t\n')
    f.write('query sequences\t'+str(cont_num)+'\n')
    f.write('reference sequences\t'+str(ref_num)+'\n')
    f.close()


def OUTPUT_MAPPED_BLOCKS_TO_REF(struct_dict,ref_dict,ref_names,ref_full_names_dict,cont_dict,contig_full_names_dict,working_dir, prefix,asmb_name_full,ref_name_full):

    

    mapped_blocks_dict={}
    for r_name in ref_dict.keys():
        mapped_blocks_dict[r_name]=[]

    for cont in struct_dict.keys():
        for trl in struct_dict[cont].keys():
            for msj in struct_dict[cont][trl]['blocks'].keys():
                for bl in struct_dict[cont][trl]['blocks'][msj]['blocks'].keys():
                    blk_coord=struct_dict[cont][trl]['blocks'][msj]['blocks'][bl]['block']

                    ref_name=blk_coord[8]
                    mapped_blocks_dict[ref_name].append([blk_coord[2],blk_coord[3],cont,blk_coord[0],blk_coord[1]])

    for ref in mapped_blocks_dict.keys():
        mapped_blocks_dict[ref]=sorted(mapped_blocks_dict[ref],key=lambda inter:inter[0], reverse=False)

    
    f=open(working_dir+prefix+'_mapped_blocks.gff','w')
    
    f.write('##gff-version 3\n')

    if ref_name_full=='yes':

        for r_name in ref_names:
            
            
            if mapped_blocks_dict.has_key(r_name):

        
                if mapped_blocks_dict[r_name]!=[]:
                    f.write('##sequence-region\t'+ref_full_names_dict[r_name]+'\t1\t'+str(len(ref_dict[r_name]))+'\n')

                    for entry in mapped_blocks_dict[r_name]:
                        if asmb_name_full=='yes':
                            f.write(ref_full_names_dict[r_name]+'\t.\t'+'MappedBlock'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t+\t.\tName='+contig_full_names_dict[entry[2]]+';length='+str(entry[1]-entry[0]+1)+';query_length='+str(len(cont_dict[entry[2]]))+';query_coord='+str(entry[3])+'-'+str(entry[4])+'\n')
                        else:
                            f.write(ref_full_names_dict[r_name]+'\t.\t'+'MappedBlock'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t+\t.\tName='+entry[2]+';length='+str(entry[1]-entry[0]+1)+';query_length='+str(len(cont_dict[entry[2]]))+';query_coord='+str(entry[3])+'-'+str(entry[4])+'\n')
    else:
        for r_name in ref_names:
            
            
            if mapped_blocks_dict.has_key(r_name):

        
                if mapped_blocks_dict[r_name]!=[]:
                    f.write('##sequence-region\t'+r_name+'\t1\t'+str(len(ref_dict[r_name]))+'\n')

                    for entry in mapped_blocks_dict[r_name]:
                        if asmb_name_full=='yes':
                            f.write(r_name+'\t.\t'+'MappedBlock'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t+\t.\tName='+contig_full_names_dict[entry[2]]+';length='+str(entry[1]-entry[0]+1)+';query_length='+str(len(cont_dict[entry[2]]))+';query_coord='+str(entry[3])+'-'+str(entry[4])+'\n')
                        else:
                            f.write(r_name+'\t.\t'+'MappedBlock'+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t+\t.\tName='+entry[2]+';length='+str(entry[1]-entry[0]+1)+';query_length='+str(len(cont_dict[entry[2]]))+';query_coord='+str(entry[3])+'-'+str(entry[4])+'\n')
    
    f.close()
    
def OUTPUT_BLOCKS_TO_QUERY(cont_blocks_dict,ref_dict,contig_names,ref_full_names_dict,cont_dict,contig_full_names_dict,working_dir, prefix,asmb_name_full,ref_name_full):

    f=open(working_dir+prefix+'_query_blocks.gff','w')
    
    f.write('##gff-version 3\n')

    dir_dict={-1:'-',1:'+'}

    gff_dict={'block':'Block', 'translocation_block':'TRL','relocation_block':'RLC', 'inversion':'INV'}

    for c_name in contig_names:
        if asmb_name_full=='yes':
            asmb_name=contig_full_names_dict[c_name]
        else:
            asmb_name=c_name

        if cont_blocks_dict.has_key(c_name):

            if cont_blocks_dict[c_name]!=[]:

                f.write('##sequence-region\t'+asmb_name+'\t1\t'+str(len(cont_dict[c_name]))+'\n')

                for entry in cont_blocks_dict[c_name]:
                    if len(entry)==8:
                    
                        if ref_name_full=='yes':
                            rf_name=ref_full_names_dict[entry[4]]
                        else:
                            rf_name=entry[4]

                        if entry[2].startswith('reshuf'):
                            field_2='RSH'
                        else:
                            field_2=gff_dict[entry[2]]

                        f.write(asmb_name+'\t.\t'+field_2+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+dir_dict[entry[7]]+'\t.\tName='+entry[2]+';length='+str(entry[3])+';ref_name='+rf_name+';ref_length='+str(len(ref_dict[entry[4]]))+';ref_coord='+str(entry[5])+'-'+str(entry[6])+'\n')

                    else:
                        field_2=gff_dict[entry[2]]
                        f.write(asmb_name+'\t.\t'+field_2+'\t'+str(entry[0])+'\t'+str(entry[1])+'\t.\t'+'.'+'\t.\tName='+entry[2]+';length='+str(entry[3])+'\n')

        
                        
    f.close()
    

def GENERATE_OUTPUT(struct_dict,end_err_dict,unmapped_list, file_ref, file_contigs,working_dir, prefix,err_ref_cont_coord_errors_list,statistics_output_lines,asmb_name_full,ref_name_full,cont_blocks_dict):

    contigs_dict, contig_seqs, contig_names, contig_full_names_dict=general.READ_FASTA_ENTRY(file_contigs)
    ref_dict, ref_seqs, ref_names,ref_full_names_dict=general.READ_FASTA_ENTRY(file_ref)

    
    #OUTPUT_READABLE(struct_dict,end_err_dict,unmapped_list, contigs_dict,contig_names, contig_full_names_dict, ref_dict,ref_full_names_dict,working_dir, 'results/'+prefix)
    OUTPUT_REF_ASSEM_TABLE(err_ref_cont_coord_errors_list, ref_dict,ref_names,ref_full_names_dict,contigs_dict,contig_names,contig_full_names_dict,working_dir, 'results/'+prefix,asmb_name_full,ref_name_full)
    OUTPUT_STAT(statistics_output_lines,working_dir, 'results/'+prefix, len(contigs_dict.keys()),len(ref_dict.keys()))
      
    OUTPUT_MAPPED_BLOCKS_TO_REF(struct_dict,ref_dict,ref_names,ref_full_names_dict,contigs_dict,contig_full_names_dict,working_dir, 'results/'+prefix,asmb_name_full,ref_name_full)

    OUTPUT_BLOCKS_TO_QUERY(cont_blocks_dict, ref_dict,contig_names,ref_full_names_dict,contigs_dict,contig_full_names_dict,working_dir, 'results/'+prefix,asmb_name_full,ref_name_full)
