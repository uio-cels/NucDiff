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

#-----------------------------------------------------------------------------------------------------------------------check functions----------------------------------------
def RESTORE_FRAGMENT_LOCAL_ERRORS(old_contig_subseq, c_st,local_errors_list):

    correctness_dict={}

    
    #1. mark each position with difference type and length

    for i in range(len(old_contig_subseq)):
        correctness_dict[i]=['correct',0]


    for err in local_errors_list:
        
      err_type=err[2]
      err_st=err[0]
      err_end=err[1]
      

      if err_st>=c_st:
      

      
         
          if err_type.startswith('insertion') or err_type.startswith('wrong_gap') :
                
                for i in range(err_st, err_end+1):

                    
                    if correctness_dict[i-c_st][0]=='correct':
                        correctness_dict[i-c_st]=['insertion',1]
                    
                    else:
                        print correctness_dict[i-c_st]
                        print  i
                        print local_errors_list
                        #raw_input('check_local_error: insertion or wrong_gap after strange error. Press any key...')
                        print 'check_local_error: insertion or wrong_gap after strange error'
                
          elif err_type=='gap':
                for i in range(err_st, err_end+1):
                    if correctness_dict[i-c_st][0]=='correct':
                        correctness_dict[i-c_st][0]='correct'
                    else:
                        print correctness_dict[i-c_st]
                        print  i
                        print local_errors_list
                        #raw_input('check_local_error: gap after strange error. Press any key...')
                        print 'check_local_error: gap after strange error'

          
            
              
                        
          elif err_type.startswith('deletion'):
                err_len=err[3]
                if correctness_dict[err_st-c_st][0]=='correct':
                    correctness_dict[err_st-c_st]=['deletion',err_len]
                elif correctness_dict[err_st-c_st][0]=='deletion':
                    correctness_dict[err_st-c_st][1]+=err_len
                elif correctness_dict[err_st-c_st][0]=='substitution':
                    correctness_dict[err_st-c_st][0]='subdel'
                    correctness_dict[err_st-c_st].append(err_len)
                elif correctness_dict[err_st-c_st][0]=='insertion':
                    correctness_dict[err_st-c_st][0]='insdel'
                    correctness_dict[err_st-c_st].append(err_len)
                
                elif correctness_dict[err_st-c_st][0]=='subdel':
                    correctness_dict[err_st-c_st][2]+=err_len

                elif correctness_dict[err_st-c_st][0]=='insdel':
                    correctness_dict[err_st-c_st][2]+=err_len
                else:
                    print correctness_dict[err_st-c_st]
                    print err_st
                    print local_errors_list
                    #raw_input('check_local_error: deletion after strange error. Press any key...')
                    print 'check_local_error: deletion after strange error'
            
                    

          elif err_type=='substitution' :
                for i in range(err_st, err_end+1):
                    
                    if correctness_dict[i-c_st][0]=='correct':
                        correctness_dict[i-c_st]=['substitution',1]
                    elif correctness_dict[i-c_st][0]=='deletion':
                        correctness_dict[i-c_st][0]='subdel'
                        correctness_dict[i-c_st].append(correctness_dict[i-c_st][1])
                        correctness_dict[i-c_st][1]=1
                    else:
                        print correctness_dict[i-c_st]
                        print  i
                        print local_errors_list
                        #raw_input('check_local_error: substitution after strange error. Press any key...')
                        print 'check_local_error: substitution after strange error'
            
                
    #create a correct query sequence
        
       
    new_contig=''
    
    for err in local_errors_list:
        
      err_type=err[2]
      err_st=err[0]
      err_end=err[1]
      
      

      if err_st<c_st:
          if err_type.startswith('deletion'):
              err_len=err[3]

              for i in range(err_len):
                  new_contig+='N'
      
      else:
          break
    
    for i in range(len(old_contig_subseq)):
            
                    
            if correctness_dict[i][0]=='correct':
                new_contig+=old_contig_subseq[i]
                
                
            elif correctness_dict[i][0]=='insertion':
                a='do nothing'

            elif correctness_dict[i][0]=='deletion':
                line=''
                for j in range(correctness_dict[i][1]):
                    line+='N'

                new_contig+=old_contig_subseq[i]
                new_contig+=line
            elif correctness_dict[i][0]=='substitution':
                new_contig+='N'

            elif correctness_dict[i][0]=='subdel':
                line=''
                for j in range(correctness_dict[i][2]+1):
                    line+='N'
                new_contig+=line
            elif correctness_dict[i][0]=='insdel':
                line=''
                for j in range(correctness_dict[i][2]):
                    line+='N'
                new_contig+=line
                
                
            else:
                #raw_input('check_local_error: unknown case during contig restoring. Press any key...')
                print 'check_local_error: unknown case during contig restoring' 

    
    
    
    return new_contig


def TEMP_CHECK_CORRECTNESS_ALL_ERRORS(structure_dict, contigs_dict, ref_dict):


    for cont_name in structure_dict.keys():
        cont_seq=contigs_dict[cont_name]

        print cont_name ,'check correctness'

        for transl_group in structure_dict[cont_name].keys():
            for misj_group in structure_dict[cont_name][transl_group]['blocks'].keys():
                for block_name in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys():
                    
                    c_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][0]
                    c_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][1]
                    r_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][2]
                    r_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][3]
                    c_dir=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][4]
                    r_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][8]
                    local_errors_list=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][10]

                    
                    
                    
                    for entry in local_errors_list:
                        if entry[0]>entry[1]:
                            print entry
                            #raw_input('ERROR: start > end')
                            print 'ERROR: start > end'
                        else:
                            if entry[3]<=0:
                                print entry
                                #raw_input('ERROR: error length < or = 0')
                                print 'ERROR: error length < or = 0'
                    
                    
                    old_contig_subseq=cont_seq[c_st  -1:c_end+1  -1]
                    ref_subseq=ref_dict[r_name][r_st  -1:r_end+1  -1]
                    
                    
                    new_contig_subseq=RESTORE_FRAGMENT_LOCAL_ERRORS(old_contig_subseq, c_st,local_errors_list)

                    #save restored query for quality check of transpositions
                    if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys())>1:
                        if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['invers_output']==[]:
                            structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['new_c_subseq']=new_contig_subseq
                        else:
                            if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['invers_output'][0][2]=='inversion':
                                structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['new_c_subseq']=general.COMPL_STRING(new_contig_subseq)
                            else:
                                 structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['new_c_subseq']=new_contig_subseq
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['ref_coord']=[r_name, r_st, r_end]    
                            
                    
                    #check if lengths are equal 
                    if len(new_contig_subseq)!=len(ref_subseq):
                        print '#######'
                        print cont_name, r_name
                        print c_st,c_end, c_dir
                        print r_st,r_end
                        print len(new_contig_subseq), len(ref_subseq) , len(old_contig_subseq)
                        print
                        
                        for entry in local_errors_list:
                            print entry

                        
                        print '#######'

                        
                        #raw_input('check correctness1: wrong length. Press any key...')
                        print 'check correctness1: wrong length'
                        
                    #check if bases are equal
                    else:
                    
                        if c_dir==1:
                            for i in range(len(ref_subseq)):
                                if (new_contig_subseq[i] in 'ATGCatgc') and (ref_subseq[i] in 'ATGCatgc'):
                                    if new_contig_subseq[i].upper()!=ref_subseq[i].upper():
                                        print '#######'
                                        print cont_name, r_name
                                        print c_st,c_end, c_dir
                                        print r_st,r_end
                                        print i,c_st+i-1, i+r_st-1,new_contig_subseq[i],ref_subseq[i] 
                                        print
                                        
                                        for entry in local_errors_list:
                                            print entry

                                        print
                                        
                                        print
                                        
                                        print

                                        

                                        print '#######'
                                        #raw_input('check correctness2: wrong base. Press any key...')
                                        print 'check correctness2: wrong base'
                                    break
                        else:
                            new_contig_compl=general.COMPL_STRING(new_contig_subseq)
                            
                            
                            for i in range(len(ref_subseq)):
                                if (new_contig_compl[i] in 'ATGCatgc') and (ref_subseq[i] in 'ATGCatgc'):
                                    if new_contig_compl[i].upper()!=ref_subseq[i].upper():
                                        print '#######'
                                        print cont_name, r_name
                                        print c_st,c_end, c_dir
                                        print r_st,r_end
                                        print i, c_end-i+1,i+r_st-1,new_contig_compl[i],ref_subseq[i] 
                                        print
                                        
                                        for entry in local_errors_list:
                                            print entry

                                        print
                                        
                                        print
                                        
                                        print

                                        print '#######'
                                        #raw_input('check correctness3: wrong base. Press any key...')
                                        print 'check correctness3: wrong base'
                                    break

                #create fragments consisting of transposed or inverted fragments            
                if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys())>1:
                    new_fragment=''
                    for i in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])):
                        for  block_name in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys():
                            if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][9]==i:
                                new_fragment+=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['new_c_subseq']

                                if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output']!=[]:
                                    for entry in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output']:
                                        if entry[2]=='deletion':
                                            for j in range(entry[3]):
                                                new_fragment+='N'
                                if i==0:
                                    r_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['ref_coord'][1]
                                    r_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['ref_coord'][0]
                                elif i==len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])-1:
                                    r_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['ref_coord'][2]

                                break
                                    

                    
                    new_ref_fragment=ref_dict[r_name][r_st  -1:r_end+1  -1]

                    #check lengths of restored fragments
                    if len(new_fragment)!=len(new_ref_fragment):
                            print '#######'
                            print cont_name, r_name
                            
                            print r_st,r_end
                            print len(new_fragment), len(new_ref_fragment)
                            print
                            
                            print '#######'

                            
                            #raw_input('check correctness: wrong transp length. Press any key...')
                            print 'check correctness: wrong transp length'
                    else: #check if bases are equal
                        for i in range(len(new_ref_fragment)):
                                if (new_fragment[i] in 'ATGCatgc') and (new_ref_fragment[i] in 'ATGCatgc'):
                                    if new_fragment[i].upper()!=new_ref_fragment[i].upper():
                                        print '#######'
                                        print cont_name, r_name
                                        
                                        print r_st,r_end
                                        print i,i+r_st-1,new_fragment[i],new_ref_fragment[i] 
                                        print
                                        
                                        for entry in local_errors_list:
                                            print entry

                                        print '#######'
                                        #raw_input('check correctness: wrong transp base. Press any key...')
                                        print 'check correctness: wrong transp base'
                                    break
                            
            
    

    for cont_name in structure_dict.keys():
        for transl_group in structure_dict[cont_name].keys():
            for misj_group in structure_dict[cont_name][transl_group]['blocks'].keys():
                for block_name in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys():
                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['new_c_subseq']=''

    
                     
                         

