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

import sys
import os

from Bio import SeqIO




def READ_FASTA_ENTRY(file_name):
    fasta_sequences=[]
    sequence_name=[]
    full_names_dict={}

    sequences_dict={}

    
    if os.stat(file_name)[6]!=0: #not empty
       
        fh = open(file_name, "r")
        for record in SeqIO.parse(fh, "fasta"):
            short_name=str(record.id).split(' ')[0]
            sequences_dict[short_name]=str(record.seq)

        
    return sequences_dict



#------------------------------------------------------------------------------------------





def FIND_CORR_INTERVAL_AND_ERRORS(entry, cont_name):

    interv_list=[]

    c_st=entry[0]
    end_c=entry[1]
    start_r=entry[2]
    end_r=entry[3]
    ref_name=entry[8]
    c_dir=entry[4]
    local_errors_list=sorted(entry[10],key=lambda inter:inter[0], reverse=False)
    len_r=entry[7]

    
    
    if c_dir==1:
        

        ins_len=0
        del_len=0

        prev_type='correct'
        prev_st=0
        prev_length=0

        for i in range(len(local_errors_list)):
            err=local_errors_list[i]
            

            err_type=err[2]
            err_st=err[0]
            err_end=err[1]
            err_len=err[3]
            err_source=err[4]

            if err_type.startswith('insertion') or err_type.startswith('wrong_gap') :
                ins_st_r=err_st-1-ins_len+del_len - c_st+start_r
                interv_list.append([ins_st_r,ins_st_r,err_type,err_len,err_st,err_end,err_source])
                ins_len+=err_len

            elif err_type.startswith('deletion'):
                del_st_r=err_st+1-ins_len+del_len -c_st+start_r
                del_end_r=err_st+err_len-ins_len+del_len -c_st+start_r
                interv_list.append([del_st_r,del_end_r,err_type,err_len,err_st,err_end,err_source])
                del_len+=err_len

            elif err_type.startswith('substitution') or err_type=='gap':
                if prev_type.startswith('deletion') and prev_st==err_st:
                    subst_st_r=err_st-ins_len+del_len -c_st+start_r-err_len-1
                    subst_end_r=err_st+err_len-1-ins_len+del_len - c_st+start_r-err_len-1
                else:
                    subst_st_r=err_st-ins_len+del_len -c_st+start_r
                    subst_end_r=err_st+err_len-1-ins_len+del_len - c_st+start_r
                interv_list.append([subst_st_r,subst_end_r,err_type,err_len,err_st,err_end,err_source])

            prev_type=err_type
            prev_st=err_st
            prev_length=err_len

            if err_source!='snps':
                interv_list[-1].append(err[5])
                
    else:
        ins_len=0
        del_len=0

        prev_type='correct'
        prev_st=0
        prev_length=0

        
        for i in range(len(local_errors_list)):
            err=local_errors_list[i]

            err_type=err[2]
            err_st=err[0]
            err_end=err[1]
            err_len=err[3]
            err_source=err[4]


            if err_type.startswith('insertion') or err_type.startswith('wrong_gap') :
               
                ins_st_r=end_r-(err_st-c_st+1)+1+ins_len-del_len
                interv_list.append([ins_st_r,ins_st_r,err_type,err_len,err_st,err_end,err_source])
                ins_len+=err_len

            
            elif err_type.startswith('deletion'):
                del_end_r=end_r-(err_st-c_st+1)+ins_len-del_len
                del_st_r=end_r-(err_st-c_st+1)+ins_len-del_len-err_len+1
                interv_list.append([del_st_r,del_end_r,err_type,err_len,err_st,err_end,err_source])
                del_len+=err_len

            elif err_type.startswith('substitution') or err_type=='gap':
                
                if prev_type.startswith('deletion') and prev_st==err_st:
                    subst_end_r=end_r-(err_st-c_st+1)+ins_len-del_len+prev_length+1
                    subst_st_r=end_r-(err_end-c_st+1)+ins_len-del_len+prev_length+1
                    interv_list.append([subst_st_r,subst_end_r,err_type,err_len,err_st,err_end,err_source])
                    
                else:
                   
                    subst_end_r=end_r-(err_st-c_st)+ins_len-del_len
                    subst_st_r=end_r-(err_end-c_st)+ins_len-del_len
                    interv_list.append([subst_st_r,subst_end_r,err_type,err_len,err_st,err_end,err_source])

            prev_type=err_type
            prev_st=err_st
            prev_length=err_len

            if err_source!='snps':
                interv_list[-1].append(err[5])


        


    return interv_list


def FIND_LOCAL_ERRORS(err_ref_cont_coord_errors_list,block_coord,cont_name):

    if block_coord[10]!=[]:
        ref_name=block_coord[8]
        c_dir=block_coord[4]
        interv_list=FIND_CORR_INTERVAL_AND_ERRORS(block_coord, cont_name)

        
        
        for entry in interv_list:
            r_st=entry[0]
            r_end=entry[1]
            err_type=entry[2]
            err_len=entry[3]
            c_st=entry[4]
            c_end=entry[5]
            err_source=entry[6]
            
            if r_st>r_end or c_st>c_end:
                print(ref_name,r_st,r_end, cont_name,c_st,c_end, err_type,err_len,'b')
                print('ERROR: wrong local difference')
            err_ref_cont_coord_errors_list.append([ref_name,r_st,r_end, cont_name,c_st,c_end, err_type,err_len,'b',c_dir,err_source])

            if err_source!='snps':
                err_ref_cont_coord_errors_list[-1].append(entry[7])
                
            
            

def FIND_TRANSP(err_ref_cont_coord_errors_list,misjoin_blocks,cont_name):
    if  len(list(misjoin_blocks.keys()))==1: 

                            block_name=list(misjoin_blocks.keys())[0]
                            block_coord=misjoin_blocks[block_name]['block']

                            FIND_LOCAL_ERRORS(err_ref_cont_coord_errors_list,block_coord, cont_name)

                            

                            
    else: 

                            for block_name in list(misjoin_blocks.keys()):
                                transp_block_coord=misjoin_blocks[block_name]['block']

                                FIND_LOCAL_ERRORS(err_ref_cont_coord_errors_list,transp_block_coord, cont_name)

                                transp_output_coord=misjoin_blocks[block_name]['transp_output']

                                if transp_output_coord!=[]:
                                    c_st=transp_output_coord[0][0]
                                    c_end=transp_output_coord[0][1]
                                    r_st=transp_output_coord[0][6]
                                    r_end=transp_output_coord[0][7]
                                    r_name=transp_output_coord[0][5]
                                    err_name=transp_output_coord[0][2]
                                    err_len=transp_output_coord[0][3]
                                    c_dir=transp_block_coord[4]
                                    
                                    

                                    if r_st>r_end or c_st>c_end:
                                            print(r_name,r_st,r_end,cont_name, c_st,c_end,err_name,err_len,'b')
                                            print(transp_output_coord)
                                            print('ERROR: wrong reshuffling output')

                                    err_ref_cont_coord_errors_list.append([r_name,r_st,r_end,cont_name, c_st,c_end,err_name,err_len,'b',c_dir])

                                between_output_coord=misjoin_blocks[block_name]['between_output']
                                if between_output_coord!=[]:
                                        
                                         

                                        for entry in between_output_coord:
                                           if entry[6]>entry[7] or entry[1]>entry[2]:
                                                print('ERROR: wrong difference between reshuffling blocks')
                                           err_ref_cont_coord_errors_list.append([entry[5],entry[6],entry[7], cont_name, entry[0],entry[1],entry[2],entry[3],'b']) 
                                        

                                invers_output_coord=misjoin_blocks[block_name]['invers_output']
                                if invers_output_coord!=[]:
                                    if invers_output_coord[0][2]=='inversion':
                                        c_st=invers_output_coord[0][0]
                                        c_end=invers_output_coord[0][1]
                                        r_st=invers_output_coord[0][6]
                                        r_end=invers_output_coord[0][7]
                                        r_name=invers_output_coord[0][5]
                                        err_name=invers_output_coord[0][2]
                                        err_len=invers_output_coord[0][3]
                                        c_dir=transp_block_coord[4]

                                        if r_st>r_end or c_st>c_end:
                                            print('ERROR: wrong inversion')

                                        err_ref_cont_coord_errors_list.append([r_name,r_st,r_end,cont_name, c_st,c_end,err_name,err_len,'b',c_dir])

                               

def FIND_MISJ_BLOCK_COORD(err_ref_cont_coord_errors_list,misjoin_blocks,cont_name):
    c_list=[]

    for block_name in list(misjoin_blocks.keys()):
                block_coord=misjoin_blocks[block_name]['block']

                c_st=block_coord[0]
                c_end=block_coord[1]
                r_st=block_coord[2]
                r_end=block_coord[3]
                c_dir=block_coord[4]
                r_name=block_coord[8]
                

                c_list.append([c_st,c_end,r_st,r_end,c_dir,r_name])
                

                block_between=misjoin_blocks[block_name]['between_output']
                
                if block_between!=[]:
                    c_st=block_between[0][0]
                    c_end=block_between[0][1]
                    r_st=block_between[0][6]
                    r_end=block_between[0][7]
                    r_name=block_between[0][5]
                    c_list.append([c_st,c_end,r_st,r_end,1,r_name])
                    
    sort_c_st_list=sorted(c_list, key=lambda inter:inter[0], reverse=False)
    c_st=sort_c_st_list[0][0]
    c_dir=sort_c_st_list[0][4]
    r_name=sort_c_st_list[0][5]

    if c_dir==1:
        r_st=sort_c_st_list[0][2]

        err_ref_cont_coord_errors_list.append([r_name, r_st, r_st, cont_name, c_st, c_st, 'misjoin_st',1,'r' ])
        
             
    else:
        r_st=sort_c_st_list[0][3]
        err_ref_cont_coord_errors_list.append([r_name, r_st, r_st, cont_name, c_st, c_st, 'misjoin_end',1,'r' ])
        

        
    

    sort_c_st_list=sorted(c_list, key=lambda inter:inter[1], reverse=False)
    c_end=sort_c_st_list[-1][1]
    c_dir=sort_c_st_list[-1][4]
    r_name=sort_c_st_list[-1][5]

    if c_dir==1:
        r_end=sort_c_st_list[-1][3]
        err_ref_cont_coord_errors_list.append([r_name, r_end, r_end, cont_name, c_end, c_end, 'misjoin_end',1,'r' ])
        
    
    else:
        r_end=sort_c_st_list[-1][2]
        err_ref_cont_coord_errors_list.append([r_name, r_end, r_end, cont_name, c_end, c_end, 'misjoin_st',1,'r' ])
        

            


def FIND_TRL_BLOCK_COORD(err_ref_cont_coord_errors_list,translocation_block,cont_name):

    c_list=[]

    for misj_gr in list(translocation_block.keys()):
        
        for block_name in list(translocation_block[misj_gr]['blocks'].keys()):
                block_coord=translocation_block[misj_gr]['blocks'][block_name]['block']

                c_st=block_coord[0]
                c_end=block_coord[1]
                r_st=block_coord[2]
                r_end=block_coord[3]
                c_dir=block_coord[4]
                r_name=block_coord[8]

                c_list.append([c_st,c_end,r_st,r_end,c_dir,r_name])
                

                block_between=translocation_block[misj_gr]['blocks'][block_name]['between_output']
                
                if block_between!=[]:
                    c_st=block_between[0][0]
                    c_end=block_between[0][1]
                    r_st=block_between[0][6]
                    r_end=block_between[0][7]
                    r_name=block_between[0][5]
                    c_list.append([c_st,c_end,r_st,r_end,1,r_name])
                    
                


    sort_c_st_list=sorted(c_list, key=lambda inter:inter[0], reverse=False)
    c_st=sort_c_st_list[0][0]
    c_dir=sort_c_st_list[0][4]
    r_name=sort_c_st_list[0][5]

    if c_dir==1:
        r_st=sort_c_st_list[0][2]
        err_ref_cont_coord_errors_list.append([r_name, r_st, r_st, cont_name, c_st, c_st, 'translocation_st',1,'r' ])
        
             
    else:
        r_st=sort_c_st_list[0][3]
        err_ref_cont_coord_errors_list.append([r_name, r_st, r_st, cont_name, c_st, c_st, 'translocation_end',1,'r' ])
        

        
    

    sort_c_st_list=sorted(c_list, key=lambda inter:inter[1], reverse=False)
    c_end=sort_c_st_list[-1][1]
    c_dir=sort_c_st_list[-1][4]
    r_name=sort_c_st_list[-1][5]

    if c_dir==1:
        r_end=sort_c_st_list[-1][3]
        err_ref_cont_coord_errors_list.append([r_name, r_end, r_end, cont_name, c_end, c_end, 'translocation_end',1,'r' ])
        
    
    else:
        r_end=sort_c_st_list[-1][2]
        err_ref_cont_coord_errors_list.append([r_name, r_end, r_end, cont_name, c_end, c_end, 'translocation_st',1,'r' ])
        

            


    
        


def FIND_REF_DIR_TRL(end_type,cur_trl_group,c_point):

    

    if end_type=='end':
        for msj_gr in list(cur_trl_group['blocks'].keys()):
            for block_name in list(cur_trl_group['blocks'][msj_gr]['blocks'].keys()):
                block_line=cur_trl_group['blocks'][msj_gr]['blocks'][block_name]['block']
                if block_line[1]==c_point:
                    r_st=block_line[2]
                    r_end=block_line[3]
                    c_dir=block_line[4]
                    r_name=block_line[8]
                   


                block_between=cur_trl_group['blocks'][msj_gr]['blocks'][block_name]['between_output']
                
                if block_between!=[]:
                    
                   for i in range(len(block_between)): 

                    if block_between[i][1]==c_point:
                        r_st=block_between[i][6]
                        r_end=block_between[i][7]
                        c_dir=1
                        r_name=block_between[i][5]
                    
                    



    elif end_type=='start':
         for msj_gr in list(cur_trl_group['blocks'].keys()):
            for block_name in list(cur_trl_group['blocks'][msj_gr]['blocks'].keys()):
                block_line=cur_trl_group['blocks'][msj_gr]['blocks'][block_name]['block']
                if block_line[0]==c_point:
                    r_st=block_line[2]
                    r_end=block_line[3]
                    c_dir=block_line[4]
                    r_name=block_line[8]


                block_between=cur_trl_group['blocks'][msj_gr]['blocks'][block_name]['between_output']
                
                if block_between!=[]:
                    
                   for i in range(len(block_between)): 
                    if block_between[i][0]==c_point:
                        r_st=block_between[i][6]
                        r_end=block_between[i][7]
                        c_dir=1
                        r_name=block_between[i][5]
            
            
            
        

    return r_st,r_end,c_dir,r_name   

def FIND_REF_DIR_MISJ(end_type,cur_misjoin_group,c_point):

    

    if end_type=='end':
        for block_name in list(cur_misjoin_group['blocks'].keys()):
            block_line=cur_misjoin_group['blocks'][block_name]['block']
            if block_line[1]==c_point:
                r_st=block_line[2]
                r_end=block_line[3]
                c_dir=block_line[4]
                r_name=block_line[8]
                r_len=block_line[7]
                

            block_between=cur_misjoin_group['blocks'][block_name]['between_output']
                
            if block_between!=[]:
                    
                for i in range(len(block_between)): 
                    if block_between[i][1]==c_point:
                        r_st=block_between[i][6]
                        r_end=block_between[i][7]
                        c_dir=1
                        r_name=block_between[i][5]
                        r_len=block_line[7]
                    



    elif end_type=='start':
        
        for block_name in list(cur_misjoin_group['blocks'].keys()):
            block_line=cur_misjoin_group['blocks'][block_name]['block']
            
            if block_line[0]==c_point:
                r_st=block_line[2]
                r_end=block_line[3]
                c_dir=block_line[4]
                r_name=block_line[8]
                r_len=block_line[7]
                
        
            block_between=cur_misjoin_group['blocks'][block_name]['between_output']
                
            if block_between!=[]:
                    
                for i in range(len(block_between)): 
                    if block_between[i][0]==c_point:
                        r_st=block_between[i][6]
                        r_end=block_between[i][7]
                        c_dir=1
                        r_name=block_between[i][5]
                        r_len=block_line[7]
            
            
            
      

    return r_st,r_end,c_dir,r_name, r_len

def TRANSLOCATION_REASON(err_ref_cont_coord_errors_list,cur_trl_group,next_trl_group,cont_name):
    coord_reason=cur_trl_group['reason'][0]
    err_ref_cont_coord_errors_list.append(['-','-','-',cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'a'])


    cur_output_line=cur_trl_group['output_line']
    c_end_cur=cur_output_line[1]
    r_st_cur, r_end_cur,c_dir_end_cur,r_name=FIND_REF_DIR_TRL('end',cur_trl_group,c_end_cur)


    next_output_line=next_trl_group['output_line']
    c_st_next=next_output_line[0]
    r_st_next, r_end_next,c_dir_st_next,r_name_next=FIND_REF_DIR_TRL('start',next_trl_group,c_st_next)


    translocation_reason=coord_reason[2]
    if translocation_reason=='translocation':
        if c_dir_end_cur==1:
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,c_end_cur,c_end_cur,translocation_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,c_end_cur,c_end_cur,translocation_reason+'_end',0,'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,c_end_cur,c_end_cur,translocation_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,c_end_cur,c_end_cur,translocation_reason+'_end',0,'a'])

        if c_dir_st_next==1:
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,c_st_next,c_st_next,translocation_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,c_st_next,c_st_next,translocation_reason+'_st',0,'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,c_st_next,c_st_next,translocation_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,c_st_next,c_st_next,translocation_reason+'_st',0,'a'])

                                

    elif translocation_reason=='translocation-insertion' or translocation_reason=='translocation-mixed_fragments' or translocation_reason=='translocation-wrong_scaffolding':
        if c_dir_end_cur==1:
            if  coord_reason[1]>coord_reason[2]:
                print('ERROR: wrong translocation reason')
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_end',coord_reason[3],'a'])
        else:
            if  coord_reason[1]>coord_reason[2]:
                print('ERROR: wrong translocation reason')
            err_ref_cont_coord_errors_list.append([r_name,max(1,r_st_cur-1),max(1,r_st_cur-1),cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,max(1,r_st_cur-1),max(1,r_st_cur-1),cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_end',coord_reason[3],'a'])

        if c_dir_st_next==1:
            if  coord_reason[1]>coord_reason[2]:
                print('ERROR: wrong translocation reason')
            err_ref_cont_coord_errors_list.append([r_name_next,max(r_st_next-1,1),max(1,r_st_next-1),cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,max(r_st_next-1,1),max(1,r_st_next-1),cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_st',coord_reason[3],'a'])
        else:
            if  coord_reason[1]>coord_reason[2]:
                print('ERROR: wrong translocation reason')
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_st',coord_reason[3],'a'])

    elif translocation_reason=='translocation-overlap':
        if c_dir_end_cur==1:
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_end',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_end',coord_reason[3],'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_st',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_end',coord_reason[3],'a'])

        if c_dir_st_next==1:
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_st',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_st',coord_reason[3],'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_end',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_st',coord_reason[3],'a'])

        

                            
def  MISJOIN_REASON(err_ref_cont_coord_errors_list,cur_misjoin_group,next_misjoin_group,cont_name):

    coord_reason=cur_misjoin_group['reason'][0]
    err_ref_cont_coord_errors_list.append(['-','-','-',cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'a'])


    cur_output_line=cur_misjoin_group['output_line']
    c_end_cur=cur_output_line[1]
    

    r_st_cur, r_end_cur,c_dir_end_cur,r_name, r_len=FIND_REF_DIR_MISJ('end',cur_misjoin_group,c_end_cur)


    next_output_line=next_misjoin_group['output_line']
    c_st_next=next_output_line[0]

    
    r_st_next, r_end_next,c_dir_st_next,r_name_next, r_len_next=FIND_REF_DIR_MISJ('start',next_misjoin_group,c_st_next)



    misjoin_reason=coord_reason[2]

    if misjoin_reason=='circular_genome_start':
        err_ref_cont_coord_errors_list.append([r_name, 1,1,cont_name,coord_reason[0],coord_reason[0],misjoin_reason,1,'r'])
        err_ref_cont_coord_errors_list.append([r_name, r_len_next,r_len_next,cont_name,coord_reason[0],coord_reason[0],misjoin_reason,1,'r'])


        if len(cur_misjoin_group['reason'])>1:
            for entry in cur_misjoin_group['reason'][1:]:
                if entry[2]=='deletion-collapsed_repeat':
                    if entry[4]=='12.1':
                        err_ref_cont_coord_errors_list.append([r_name, 1,entry[3],cont_name,entry[0],entry[1],entry[2],entry[3],'r'])
                        err_ref_cont_coord_errors_list.append([r_name, 1,entry[3],cont_name,entry[0],entry[1],entry[2],entry[3],'a'])
                    else:
                        err_ref_cont_coord_errors_list.append([r_name, r_len-entry[3]+1,r_len,cont_name,entry[0],entry[1],entry[2],entry[3],'r'])
                        err_ref_cont_coord_errors_list.append([r_name, r_len-entry[3]+1,r_len,cont_name,entry[0],entry[1],entry[2],entry[3],'a'])
                else:
                    err_ref_cont_coord_errors_list.append([r_name, r_len,r_len,cont_name,entry[0],entry[1],entry[2],entry[3],'r'])
                    err_ref_cont_coord_errors_list.append([r_name, r_len,r_len,cont_name,entry[0],entry[1],entry[2],entry[3],'a'])
            
                           
    elif misjoin_reason=='misjoin':
        if c_dir_end_cur==1:
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,c_end_cur,c_end_cur,misjoin_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,c_end_cur,c_end_cur,misjoin_reason+'_end',0,'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,c_end_cur,c_end_cur,misjoin_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,c_end_cur,c_end_cur,misjoin_reason+'_end',0,'a'])
        if c_dir_st_next==1:
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,c_st_next,c_st_next,misjoin_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,c_st_next,c_st_next,misjoin_reason+'_st',0,'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,c_st_next,c_st_next,misjoin_reason,0,'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,c_st_next,c_st_next,misjoin_reason+'_st',0,'a'])

                                

    elif misjoin_reason=='misjoin-insertion' or misjoin_reason=='misjoin-mixed_fragments' or misjoin_reason=='misjoin-wrong_scaffolding':
        if c_dir_end_cur==1:
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_end',coord_reason[3],'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name,max(1,r_st_cur-1),max(1,r_st_cur-1),cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,max(1,r_st_cur-1),max(1,r_st_cur-1),cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_end',coord_reason[3],'a'])

        if c_dir_st_next==1:
            err_ref_cont_coord_errors_list.append([r_name_next,max(r_st_next-1,1),max(1,r_st_next-1),cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,max(r_st_next-1,1),max(1,r_st_next-1),cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_st',coord_reason[3],'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[0],coord_reason[1],coord_reason[2],coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_st',coord_reason[3],'a'])

    elif misjoin_reason=='misjoin-overlap':
        if c_dir_end_cur==1:
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_end',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_end_cur,r_end_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_end',coord_reason[3],'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_st',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name,r_st_cur,r_st_cur,cont_name,coord_reason[1],coord_reason[1],coord_reason[2]+'_end',coord_reason[3],'a'])

        if c_dir_st_next==1:
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_st',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_st_next,r_st_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_st',coord_reason[3],'a'])
        else:
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_end',coord_reason[3],'r'])
            err_ref_cont_coord_errors_list.append([r_name_next,r_end_next,r_end_next,cont_name,coord_reason[0],coord_reason[0],coord_reason[2]+'_st',coord_reason[3],'a'])

        
                                    
                                
def FIND_MISJ(err_ref_cont_coord_errors_list,translocation_block,cont_name):
    temp_list=[]
    
    if len(list(translocation_block.keys()))==1: #no misjoins
        msj_gr=list(translocation_block.keys())[0]
        misjoin_blocks=translocation_block[msj_gr]['blocks']

        FIND_TRANSP(err_ref_cont_coord_errors_list,misjoin_blocks,cont_name)

        
       
    else:
                        
        for msj_gr in list(translocation_block.keys()):

            #print translocation_block[msj_gr]
            #raw_input('hgf')

            
            misjoin_blocks=translocation_block[msj_gr]['blocks']
            FIND_TRANSP(err_ref_cont_coord_errors_list,misjoin_blocks,cont_name)

            FIND_MISJ_BLOCK_COORD(temp_list,misjoin_blocks,cont_name)
            
            coord_reason=translocation_block[msj_gr]['reason']
            if not coord_reason==[]:
                cur_misjoin_group=translocation_block[msj_gr]
                next_misjoin_group=translocation_block[msj_gr+1]
                MISJOIN_REASON(temp_list,cur_misjoin_group,next_misjoin_group,cont_name)

            #for entry in temp_list:
            #    print entry
            #raw_input('jh')

            #for i in range(len(temp_list)):
            #    temp_list.pop(0)

#------------------------------------------------------------------------------------------------------------------------------------------------------
def FIND_CONTIGS_ENDS(err_ref_cont_coord_errors_list,struct_dict_contig, cont_name,end_err_dict):
    c_list=[]
    

    for trl_gr in list(struct_dict_contig.keys()):
        for misj_gr in list(struct_dict_contig[trl_gr]['blocks'].keys()):
            for block_name in list(struct_dict_contig[trl_gr]['blocks'][misj_gr]['blocks'].keys()):
                block_coord=struct_dict_contig[trl_gr]['blocks'][misj_gr]['blocks'][block_name]['block']
                c_st=block_coord[0]
                c_end=block_coord[1]
                r_st=block_coord[2]
                r_end=block_coord[3]
                c_dir=block_coord[4]
                r_name=block_coord[8]

                c_list.append([c_st,c_end,r_st,r_end,c_dir,r_name])
                

                block_between=struct_dict_contig[trl_gr]['blocks'][misj_gr]['blocks'][block_name]['between_output']
                
                if block_between!=[]:
                    c_st=block_between[0][0]
                    c_end=block_between[0][1]
                    r_st=block_between[0][6]
                    r_end=block_between[0][7]
                    r_name=block_between[0][5]
                    c_list.append([c_st,c_end,r_st,r_end,1,r_name])
                    
                


    sort_c_st_list=sorted(c_list, key=lambda inter:inter[0], reverse=False)

    c_st=sort_c_st_list[0][0]
    c_dir=sort_c_st_list[0][4]
    r_name=sort_c_st_list[0][5]

    
    if c_dir==1:
        r_st=sort_c_st_list[0][2]
        err_ref_cont_coord_errors_list.append([r_name, r_st, r_st, cont_name, c_st, c_st, 'contig_st',1,'r' ])

        if end_err_dict[cont_name]['wrong_beginning']!=[]:
            end_entry=end_err_dict[cont_name]['wrong_beginning'][0]
            
            err_ref_cont_coord_errors_list.append([r_name,max(1,r_st-1),max(1,r_st-1),cont_name,end_entry[0],end_entry[1],end_entry[2],end_entry[3],'struct'])

        if end_err_dict[cont_name]['duplication']!=[]:

            for dupl in end_err_dict[cont_name]['duplication']:
                if dupl[2]=='wrong_beginning':
                    err_ref_cont_coord_errors_list.append([r_name,max(1,r_st-1),max(1,r_st-1),cont_name,dupl[0],dupl[1],'insertion-multiple_copy',dupl[3],'struct',1,'ends'])

 
                  
    else:
        r_st=sort_c_st_list[0][3]
        err_ref_cont_coord_errors_list.append([r_name, r_st, r_st, cont_name, c_st, c_st, 'contig_end',1,'r' ])

        if end_err_dict[cont_name]['wrong_beginning']!=[]:
            end_entry=end_err_dict[cont_name]['wrong_beginning'][0]

            err_ref_cont_coord_errors_list.append([r_name,r_st,r_st,cont_name,end_entry[0],end_entry[1],end_entry[2],end_entry[3],'struct'])


        if end_err_dict[cont_name]['duplication']!=[]:

            for dupl in end_err_dict[cont_name]['duplication']:
                if dupl[2]=='wrong_beginning':
                    err_ref_cont_coord_errors_list.append([r_name,r_st,r_st,cont_name,dupl[0],dupl[1],'insertion-multiple_copy',dupl[3],'struct',1,'ends'])


                   

    

    

    sort_c_st_list=sorted(c_list, key=lambda inter:inter[1], reverse=False)
    c_end=sort_c_st_list[-1][1]
    c_dir=sort_c_st_list[-1][4]
    r_name=sort_c_st_list[-1][5]

    if c_dir==1:
        r_end=sort_c_st_list[-1][3]
        err_ref_cont_coord_errors_list.append([r_name, r_end, r_end, cont_name, c_end, c_end, 'contig_end',1,'r' ])

        if end_err_dict[cont_name]['wrong_end']!=[]:
            end_entry=end_err_dict[cont_name]['wrong_end'][0]
            err_ref_cont_coord_errors_list.append([r_name,r_end,r_end,cont_name,end_entry[0],end_entry[1],end_entry[2],end_entry[3],'struct'])

        if end_err_dict[cont_name]['duplication']!=[]:

            for dupl in end_err_dict[cont_name]['duplication']:
                if dupl[2]=='wrong_end':
                    err_ref_cont_coord_errors_list.append([r_name,r_end,r_end,cont_name,dupl[0],dupl[1],'insertion-multiple_copy',dupl[3],'struct',1,'ends'])
            
    else:
        r_end=sort_c_st_list[-1][2]
        err_ref_cont_coord_errors_list.append([r_name, r_end, r_end, cont_name, c_end, c_end, 'contig_st',1,'r' ])

        if end_err_dict[cont_name]['wrong_end']!=[]:
            end_entry=end_err_dict[cont_name]['wrong_end'][0]
            err_ref_cont_coord_errors_list.append([r_name,max(1,r_end-1),max(1,r_end-1),cont_name,end_entry[0],end_entry[1],end_entry[2],end_entry[3],'struct'])

        if end_err_dict[cont_name]['duplication']!=[]:

            

            for dupl in end_err_dict[cont_name]['duplication']:
                if dupl[2]=='wrong_end':
                    err_ref_cont_coord_errors_list.append([r_name,max(1,r_end-1),max(1,r_end-1),cont_name,dupl[0],dupl[1],'insertion-multiple_copy',dupl[3],'struct',1,'ends'])
            

           
                            
    
    
def FIND_ERRORS_ALL_COORD1(ref_dict,cont_dict, struct_dict, end_err_dict,unmapped_list):
    err_ref_cont_coord_errors_list=[]

    for cont_name in list(cont_dict.keys()):
        #wrong query sequences

        

        if cont_name not in struct_dict:
            if cont_name in unmapped_list:
                len_cont=len(cont_dict[cont_name])
                err_ref_cont_coord_errors_list.append(['-','-','-',cont_name, 1, len_cont,'wrong_sequence',len_cont,'a'])
            else:
                print()
                print('ERROR: not handled query sequence', cont_name)
                print()
                
                
            
        else:
            #query sequence ends and wrong ends and beginnings
            FIND_CONTIGS_ENDS(err_ref_cont_coord_errors_list,struct_dict[cont_name],cont_name,end_err_dict)
           
            
            if len(list(struct_dict[cont_name].keys()))==1: #no translocations
                    trl_group=list(struct_dict[cont_name].keys())[0]
                    translocation_block=struct_dict[cont_name][trl_group]['blocks']

                    FIND_MISJ(err_ref_cont_coord_errors_list,translocation_block,cont_name)
            else:
                for trl_group in list(struct_dict[cont_name].keys()):
                    translocation_block=struct_dict[cont_name][trl_group]['blocks']
                    FIND_MISJ(err_ref_cont_coord_errors_list,translocation_block,cont_name)

                    trl_block_coord=FIND_TRL_BLOCK_COORD(err_ref_cont_coord_errors_list,translocation_block,cont_name)
                            
                    coord_reason=struct_dict[cont_name][trl_group]['reason']
                    if not coord_reason==[]:
                        cur_trl_group=struct_dict[cont_name][trl_group]
                        next_trl_group=struct_dict[cont_name][trl_group+1]
                        TRANSLOCATION_REASON(err_ref_cont_coord_errors_list,cur_trl_group,next_trl_group,cont_name)
                               
     
            
    
    
    return err_ref_cont_coord_errors_list


                        
def FIND_INTERSECTION(interv_list):
    intervals_list=[]



    interv_list=sorted(interv_list, key=lambda inter:inter[0], reverse=False)

    

    if len(interv_list)>0:

        cur_interv_st=interv_list[0][0]
        cur_interv_end=interv_list[0][1]


        for i in range(1,len(interv_list)):
            if cur_interv_end<interv_list[i][0]:
                intervals_list.append([cur_interv_st,cur_interv_end])
                cur_interv_st=interv_list[i][0]
                cur_interv_end=interv_list[i][1]

            else:
                cur_interv_end=max(cur_interv_end,interv_list[i][1])

        intervals_list.append([cur_interv_st,cur_interv_end])

    
    return intervals_list

def FIND_UNCOV_REF(interv_list,len_ref):
    uncovered_list=[]

    if interv_list!=[]:

        if interv_list[0][0]!=1:
            uncovered_list.append([1,interv_list[0][0]-1])

        for i in range(len(interv_list)-1):
            cur=interv_list[i]
            nex=interv_list[i+1]

            if nex[0]-cur[1]>1:
                uncovered_list.append([cur[1]+1,nex[0]-1])

        if interv_list[-1][1]<len_ref:
            uncovered_list.append([interv_list[-1][1]+1,len_ref])

    return uncovered_list


    
def FIND_UNCOVERED_REF_FRAG(err_ref_cont_coord_errors_list,struct_dict, ref_dict,uncovered_dict):
    ref_interv_dict={}

    '''
    for ref_name in uncovered_dict.keys():
        if uncovered_dict[ref_name]!=[]:
            for entry in uncovered_dict[ref_name]:
                err_ref_cont_coord_errors_list.append([ref_name,entry[0],entry[1],'-','-','-','uncovered_region',entry[1]-entry[0]+1,'r'])
    '''
    #find clipped_repeated_regios
    for ref_name in list(ref_dict.keys()):
        ref_interv_dict[ref_name]=[]

    for cont_name in list(struct_dict.keys()):
        for trl_gr in list(struct_dict[cont_name].keys()):
            for msj_gr in list(struct_dict[cont_name][trl_gr]['blocks'].keys()):
                for bl_name in list(struct_dict[cont_name][trl_gr]['blocks'][msj_gr]['blocks'].keys()):
                    block_coord=struct_dict[cont_name][trl_gr]['blocks'][msj_gr]['blocks'][bl_name]['block']
                    ref_interv_dict[block_coord[8]].append([block_coord[2],block_coord[3]])

                    between_output=struct_dict[cont_name][trl_gr]['blocks'][msj_gr]['blocks'][bl_name]['between_output']
                    if between_output!=[]:
                        for entry in between_output:
                            
                            if entry[2].startswith('deletion'):
                                ref_interv_dict[entry[5]].append([entry[6],entry[7]])

    
    for ref_name in list(ref_interv_dict.keys()):
        if ref_interv_dict[ref_name]==[]:
            err_ref_cont_coord_errors_list.append([ref_name,1,len(ref_dict[ref_name]),'-','-','-','uncovered_region',len(ref_dict[ref_name]),'r'])
        else:
            cover_reg_list=FIND_INTERSECTION(ref_interv_dict[ref_name])
            uncovered_list=FIND_UNCOV_REF(cover_reg_list,len(ref_dict[ref_name]))

            

            for entry in uncovered_list:
                err_ref_cont_coord_errors_list.append([ref_name,entry[0],entry[1],'-','-','-','uncovered_region',entry[1]-entry[0]+1,'r'])
    
    '''
    for ref_name in ref_dict.keys():
        cover_reg_list=FIND_INTERSECTION(ref_interv_dict[ref_name])

        for entry in uncovered_dict[ref_name]:
            cover_reg_list.append(entry)

        cover_reg_list=FIND_INTERSECTION(cover_reg_list)

        uncovered_list=FIND_UNCOV_REF(cover_reg_list,len(ref_dict[ref_name]))

        for entry in uncovered_list:
                err_ref_cont_coord_errors_list.append([ref_name,entry[0],entry[1],'-','-','-','clipped_repeated_region',entry[1]-entry[0]+1,'r'])

    ''' 
        
def MERGE_BLOCK_TRL(reloc_list):
    #[c_st,c_end,r_st,r_end,ref_name,[c_st,c_st_Ref,c_st_Ref_type],[c_end,c_end_Ref,c_end_Tef_type]]

    if reloc_list==[]:
        return []
    min_st=reloc_list[0][0][0]
    min_st_Ref=reloc_list[0][0][5][1]
    min_ref_type=reloc_list[0][0][5][2]

    ref_name=reloc_list[0][0][4]
    

    max_end=reloc_list[0][0][1]
    max_end_Ref=reloc_list[0][0][6][1]
    max_ref_type=reloc_list[0][0][6][2]


    for entry in reloc_list[1:]:
        cur_st=entry[0][0]
        cur_end=entry[0][1]

        min_st=min(min_st,cur_st)

        if min_st==cur_st:
            min_st_Ref=entry[0][5][1]
            min_ref_type=entry[0][5][2]

                    
        max_end=max(max_end,cur_end)

        if max_end==cur_end:
            max_end_Ref=entry[0][6][1]
            max_ref_type=entry[0][6][2]


    min_st_r=reloc_list[0][0][2]
    max_end_r=reloc_list[0][0][3]

    
    for entry in reloc_list[1:]:
        cur_st=entry[0][2]
        cur_end=entry[0][3]

        min_st_r=min(min_st_r,cur_st)
        max_end_r=max(max_end_r,cur_end)
    
    

    return [min_st,max_end,min_st_r,max_end_r,ref_name,[min_st,min_st_Ref,min_ref_type],[max_end,max_end_Ref,max_ref_type]]

    

def MERGE_BLOCK_MSJ(block_list):
    
    min_st=block_list[0][0]

    if block_list[0][4]==1:
        min_st_Ref=block_list[0][2]
        min_ref_type='st'
    else:
        min_st_Ref=block_list[0][3]
        min_ref_type='end'

    ref_name=block_list[0][8]
        
    max_end=block_list[0][1]

    if block_list[0][4]==1:
        max_end_Ref=block_list[0][3]
        max_ref_type='end'
    else:
        max_end_Ref=block_list[0][2]
        max_ref_type='st'
    

    for entry in block_list[1:]:
        cur_st=entry[0]
        cur_end=entry[1]

        min_st=min(min_st,cur_st)

        if min_st==cur_st:
            if entry[4]==1:
                min_st_Ref=entry[2]
                min_ref_type='st'
            else:
                min_st_Ref=entry[3]
                min_ref_type='end'
                    
        max_end=max(max_end,cur_end)

        if max_end==cur_end:
            if entry[4]==1:
                max_end_Ref=entry[3]
                max_ref_type='end'
            else:
                max_end_Ref=entry[2]
                max_ref_type='st'


    min_st_r=block_list[0][2]

    if block_list[0][4]==1:
        min_st_Q=block_list[0][0]
    else:
        min_st_Q=block_list[0][1]
        
    max_end_r=block_list[0][3]

    if block_list[0][4]==1:
        max_end_Q=block_list[0][1]
    else:
        max_end_Q=block_list[0][0]
    

    for entry in block_list[1:]:
        cur_st=entry[2]
        cur_end=entry[3]

        min_st_r=min(min_st_r,cur_st)

        if min_st_r==cur_st:
            if entry[4]==1:
                min_st_Q=entry[0]
            else:
                min_st_Q=entry[1]
                    
        max_end_r=max(max_end_r,cur_end)

        if max_end_r==cur_end:
            if entry[4]==1:
                max_end_Q=entry[1]
            else:
                max_end_Q=entry[0]
    
    

    
    return [min_st,max_end,min_st_r,max_end_r,ref_name,[min_st,min_st_Ref,min_ref_type],[max_end,max_end_Ref,max_ref_type]]

def FIND_CONT_BLOCKS(structure_dict):
    cont_blocks_dict={}

    temp_dict={}
    result_dict={}
    
    for cont_name in list(structure_dict.keys()):
        temp_dict[cont_name]={'transloc':[], 'misj':[],'reshuf':[],'invers':[],'blocks':[]}
        cont_blocks_dict[cont_name]={'transloc':[], 'misj':[],'reshuf':[],'invers':[],'blocks':[]}
        result_dict[cont_name]=[]

       

        for transl_group in sorted(structure_dict[cont_name].keys()):
            for misj_group in sorted(structure_dict[cont_name][transl_group]['blocks'].keys()):
                for block_name in sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                    block=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block']
                    reshuf=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['transp_output']
                    invers=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['invers_output']

                    c_dir=block[4]
                    
                    if reshuf!=[]:
                        temp_dict[cont_name]['reshuf'].append(reshuf)

                    if invers!=[]:
                        temp_dict[cont_name]['invers'].append(invers)


                    temp_dict[cont_name]['blocks'].append(block)
                    
                
                if len(temp_dict[cont_name]['blocks'])>1:
                    merge_bl=MERGE_BLOCK(temp_dict[cont_name]['blocks'])
                    temp_dict[cont_name]['misj'].append(merge_bl)

                else:
                    temp_dict[cont_name]['misj'].append(temp_dict[cont_name]['blocks'][0])

                
                for entry in temp_dict[cont_name]['reshuf']:
                    entry=entry[0]
                    cont_blocks_dict[cont_name]['reshuf'].append([entry[0],entry[1],entry[2],entry[3],entry[5],entry[6],entry[7],c_dir])

                for i in range(len(temp_dict[cont_name]['reshuf'])):
                    temp_dict[cont_name]['reshuf'].pop(0)
                    

                for entry in temp_dict[cont_name]['invers']:
                    entry=entry[0]
                    if entry[2]=='inversion':
                        cont_blocks_dict[cont_name]['invers'].append([entry[0],entry[1],entry[2],entry[3],entry[5],entry[6],entry[7],c_dir])

                for i in range(len(temp_dict[cont_name]['invers'])):
                    temp_dict[cont_name]['invers'].pop(0)



                for entry in temp_dict[cont_name]['blocks']:
                    cont_blocks_dict[cont_name]['blocks'].append([entry[0],entry[1],'block',entry[1]-entry[0]+1,entry[8],entry[2],entry[3],c_dir])
                    
                for i in range(len(temp_dict[cont_name]['blocks'])):
                    temp_dict[cont_name]['blocks'].pop(0)

            
            if len(temp_dict[cont_name]['misj'])>1:
                merge_bl=MERGE_BLOCK(temp_dict[cont_name]['misj'])
                temp_dict[cont_name]['transloc'].append(merge_bl)

                for entry in temp_dict[cont_name]['misj']:
                    cont_blocks_dict[cont_name]['misj'].append([entry[0],entry[1],'relocation_block',entry[1]-entry[0]+1])

            else:
                temp_dict[cont_name]['transloc'].append(temp_dict[cont_name]['misj'][0])

            for i in range(len(temp_dict[cont_name]['misj'])):
                    temp_dict[cont_name]['misj'].pop(0)
            
        if len(temp_dict[cont_name]['transloc'])>1:
            for entry in temp_dict[cont_name]['transloc']:
                cont_blocks_dict[cont_name]['transloc'].append([entry[0],entry[1],'translocation_block',entry[1]-entry[0]+1])

            for i in range(len(temp_dict[cont_name]['transloc'])):
                temp_dict[cont_name]['transloc'].pop(0)

            

               
        
        for gr_name in ['blocks','invers','reshuf','misj','transloc']:
            for entry in cont_blocks_dict[cont_name][gr_name]:
                result_dict[cont_name].append(entry)

        result_dict[cont_name]=sorted(result_dict[cont_name], key=lambda inter:inter[0], reverse=False)

                        

        
    return result_dict
    
    

def FIND_REF_COORD_ERRORS1(struct_dict,end_err_dict,unmapped_list,file_ref, file_contigs):
    

    ref_dict=READ_FASTA_ENTRY(file_ref)
    cont_dict=READ_FASTA_ENTRY(file_contigs)

    
    
    err_ref_cont_coord_errors_list=FIND_ERRORS_ALL_COORD1(ref_dict,cont_dict, struct_dict, end_err_dict,unmapped_list)

    FIND_UNCOVERED_REF_FRAG1(err_ref_cont_coord_errors_list,struct_dict,ref_dict)

    cont_blocks_dict=FIND_CONT_BLOCKS1(struct_dict)

    return err_ref_cont_coord_errors_list,cont_blocks_dict



def FIND_ERRORS_ALL_COORD(ref_dict,cont_dict, structure_dict, end_err_dict,unmapped_list):
    err_ref_cont_coord_errors_list=[]


    cont_blocks_dict={}

    temp_dict={}
    result_dict={}
    temp_local_list=[]
    reloc_list1=[]
    reloc_list2=[]

    for cont_name in list(cont_dict.keys()):
        #wrong query sequences
        if cont_name not in structure_dict:
            if cont_name in unmapped_list:
                len_cont=len(cont_dict[cont_name])
                err_ref_cont_coord_errors_list.append(['-','-','-',cont_name, 1, len_cont,'wrong_sequence',len_cont,'a'])
            else:
                print()
                print('ERROR: not handled query sequence', cont_name)
                print()
                
                
            
        else:
            #query sequence ends and wrong ends and beginnings
            FIND_CONTIGS_ENDS(err_ref_cont_coord_errors_list,structure_dict[cont_name],cont_name,end_err_dict)
    
    
    for cont_name in list(structure_dict.keys()):
        temp_dict[cont_name]={'transloc':[], 'misj':[],'reshuf':[],'invers':[],'blocks':[],'between':[]}
        cont_blocks_dict[cont_name]={'transloc':[], 'misj':[],'reshuf':[],'invers':[],'blocks':[]}
        result_dict[cont_name]=[]


        for transl_group in sorted(structure_dict[cont_name].keys()):
            for misj_group in sorted(structure_dict[cont_name][transl_group]['blocks'].keys()):
                for block_name in sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                    block=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block']
                    reshuf=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['transp_output']
                    invers=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['invers_output']
                    between=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output']
                    local=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['local_output']


                    c_dir=block[4]

                    if reshuf!=[]:
                        temp_dict[cont_name]['reshuf'].append(reshuf)

                    if invers!=[]:
                        temp_dict[cont_name]['invers'].append(invers)

                    if between!=[]:
                        for entry in between:
                            temp_dict[cont_name]['between'].append(entry)


                    temp_dict[cont_name]['blocks'].append(block)


                msj_reason=structure_dict[cont_name][transl_group]['blocks'][misj_group]['reason']

                
                
                merge_bl=MERGE_BLOCK_MSJ(temp_dict[cont_name]['blocks'])
                temp_dict[cont_name]['misj'].append([merge_bl,msj_reason])
                

                for entry in temp_dict[cont_name]['reshuf']:
                    entry=entry[0]
                    blck_name='R'+entry[2][1:]
                    cont_blocks_dict[cont_name]['reshuf'].append([entry[0],entry[1],blck_name,entry[3],entry[5],entry[6],entry[7],c_dir])
                    err_ref_cont_coord_errors_list.append([entry[5],entry[6],entry[7],cont_name,entry[0],entry[1],entry[2],entry[3],'struct',c_dir,entry[4]])
                for i in range(len(temp_dict[cont_name]['reshuf'])):
                    temp_dict[cont_name]['reshuf'].pop(0)
                    

                for entry in temp_dict[cont_name]['invers']:
                    entry=entry[0]
                    if entry[2]=='inversion':
                        cont_blocks_dict[cont_name]['invers'].append([entry[0],entry[1],'Inversion',entry[3],entry[5],entry[6],entry[7],c_dir])
                        err_ref_cont_coord_errors_list.append([entry[5],entry[6],entry[7],cont_name,entry[0],entry[1],entry[2],entry[3],'struct',c_dir,entry[4]])
                for i in range(len(temp_dict[cont_name]['invers'])):
                    temp_dict[cont_name]['invers'].pop(0)


                for entry in temp_dict[cont_name]['between']:
                    err_ref_cont_coord_errors_list.append([entry[5],entry[6],entry[7],cont_name,entry[0],entry[1],entry[2],entry[3],'struct',1,entry[4]])
                for i in range(len(temp_dict[cont_name]['between'])):
                    temp_dict[cont_name]['between'].pop(0)



                for entry in temp_dict[cont_name]['blocks']:
                    
                    cont_blocks_dict[cont_name]['blocks'].append([entry[0],entry[1],'Block',entry[1]-entry[0]+1,entry[8],entry[2],entry[3],c_dir])
                    
                    if entry[10]!=[]:
                        FIND_LOCAL_ERRORS(temp_local_list,entry,cont_name)
                        for loc_err in temp_local_list:
                            if loc_err[10]=='snps':
                                loc_err[8]='snps'
                            else:
                                loc_err[8]='struct'
                            err_ref_cont_coord_errors_list.append(loc_err)

                        for i in range(len(temp_local_list)):
                            temp_local_list.pop(0)

                for i in range(len(temp_dict[cont_name]['blocks'])):
                    temp_dict[cont_name]['blocks'].pop(0)

            trl_reason=structure_dict[cont_name][transl_group]['reason']
            
            if len(temp_dict[cont_name]['misj'])>1:
                flag=0
                #reloc_list1
                for i in range(len(temp_dict[cont_name]['misj'])):
                    entry=temp_dict[cont_name]['misj'][i]
                    reloc_list1.append(entry)

                    
                    if entry[1]!=[] and entry[1][0][2]=='circular_genome_start':
                        flag=1
                        break

                merge_bl1=MERGE_BLOCK_TRL(reloc_list1)

                if flag==1:
                    for j in range(i+1,len(temp_dict[cont_name]['misj'])):
                        entry=temp_dict[cont_name]['misj'][j]
                        reloc_list2.append(entry)
                    
                merge_bl2=MERGE_BLOCK_TRL(reloc_list2)

                for i in range(len(reloc_list1)):
                    reloc_list1.pop(0)
                for i in range(len(reloc_list2)):
                    reloc_list2.pop(0)

                
                temp_dict[cont_name]['transloc'].append([merge_bl1,merge_bl2,trl_reason])

                

                entry=temp_dict[cont_name]['misj'][0]
                if entry[1]!=[] and entry[1][0][2]=='circular_genome_start':
                    cont_blocks_dict[cont_name]['misj'].append([entry[0][0],entry[0][1],'Circular_genome_block',entry[0][1]-entry[0][0]+1])
                else:
                    cont_blocks_dict[cont_name]['misj'].append([entry[0][0],entry[0][1],'Relocation_block',entry[0][1]-entry[0][0]+1])

                if entry[1][0][2]=='misjoin-overlap':
                    err_ref_cont_coord_errors_list.append([cont_name,entry[1][0][0],entry[1][0][1],entry[1][0][2],entry[1][0][3],entry[1][0][4],[entry[0]],[temp_dict[cont_name]['misj'][1][0]],'struct2',[entry[1][0][5],entry[1][0][6]]])
                else:
                    err_ref_cont_coord_errors_list.append([cont_name,entry[1][0][0],entry[1][0][1],entry[1][0][2],entry[1][0][3],entry[1][0][4],[entry[0]],[temp_dict[cont_name]['misj'][1][0]],'struct2'])

                
                for i in range(1,len(temp_dict[cont_name]['misj'])):
                    entry=temp_dict[cont_name]['misj'][i]
                    if entry[1]!=[] and entry[1][0][2]=='circular_genome_start':
                        cont_blocks_dict[cont_name]['misj'].append([entry[0][0],entry[0][1],'Circular_genome_block',entry[0][1]-entry[0][0]+1])
                    elif entry[1]!=[] and entry[1][0][2]!='circular_genome_start':
                        cont_blocks_dict[cont_name]['misj'].append([entry[0][0],entry[0][1],'Relocation_block',entry[0][1]-entry[0][0]+1])
                    else:
                        if temp_dict[cont_name]['misj'][i-1][1][0][2]=='circular_genome_start':
                            cont_blocks_dict[cont_name]['misj'].append([entry[0][0],entry[0][1],'Circular_genome_block',entry[0][1]-entry[0][0]+1])
                        else:
                            cont_blocks_dict[cont_name]['misj'].append([entry[0][0],entry[0][1],'Relocation_block',entry[0][1]-entry[0][0]+1])


                    if i!=len(temp_dict[cont_name]['misj'])-1:
                        if entry[1][0][2]=='misjoin-overlap':
                            err_ref_cont_coord_errors_list.append([cont_name,entry[1][0][0],entry[1][0][1],entry[1][0][2],entry[1][0][3],entry[1][0][4],[entry[0]],[temp_dict[cont_name]['misj'][i+1][0]],'struct2',[entry[1][0][5],entry[1][0][6]]])

                        else:
                            err_ref_cont_coord_errors_list.append([cont_name,entry[1][0][0],entry[1][0][1],entry[1][0][2],entry[1][0][3],entry[1][0][4],[entry[0]],[temp_dict[cont_name]['misj'][i+1][0]],'struct2'])

            else:
                merge_bl=MERGE_BLOCK_TRL(temp_dict[cont_name]['misj'])
                temp_dict[cont_name]['transloc'].append([merge_bl,[],trl_reason])
                

            for i in range(len(temp_dict[cont_name]['misj'])):
                    temp_dict[cont_name]['misj'].pop(0)

        if len(temp_dict[cont_name]['transloc'])>1:
            
            for entry in temp_dict[cont_name]['transloc']:
                if entry[1]==[]:
                    cont_blocks_dict[cont_name]['transloc'].append([entry[0][0],entry[0][1],'Translocation_block',entry[0][1]-entry[0][0]+1])
                else:
                    cont_blocks_dict[cont_name]['transloc'].append([entry[0][0],entry[1][1],'Translocation_block',entry[1][1]-entry[0][0]+1])
                    

            #print 'done'

            for i in range(len(temp_dict[cont_name]['transloc'])-1):
                entry=temp_dict[cont_name]['transloc'][i]
                next_entry=temp_dict[cont_name]['transloc'][i+1]

               

                if entry[1]==[] and next_entry[1]==[]:

                    if entry[2][0][2]=='translocation-overlap':
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[0]],[next_entry[0]],'struct2',[entry[2][0][5],entry[2][0][6]]])
                    else:
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[0]],[next_entry[0]],'struct2'])
                    
                    
                elif entry[1]==[] and next_entry[1]!=[]:
                    if entry[2][0][2]=='translocation-overlap':
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[0]],[next_entry[0]],'struct2',[entry[2][0][5],entry[2][0][6]]])

                    else:
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[0]],[next_entry[0]],'struct2'])

                    
                elif entry[1]!=[] and next_entry[1]==[]:
                    if entry[2][0][2]=='translocation-overlap':
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[1]],[next_entry[0]],'struct2',[entry[2][0][5],entry[2][0][6]]])
                    else:
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[1]],[next_entry[0]],'struct2'])

                   
            
                else:
                    if entry[2][0][2]=='translocation-overlap':
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[1]],[next_entry[0]],'struct2',[entry[2][0][5],entry[2][0][6]]])
                    else:
                        err_ref_cont_coord_errors_list.append([cont_name,entry[2][0][0],entry[2][0][1],entry[2][0][2],entry[2][0][3],entry[2][0][4],[entry[1]],[next_entry[0]],'struct2'])

                

            
        
            for i in range(len(temp_dict[cont_name]['transloc'])):
                temp_dict[cont_name]['transloc'].pop(0)

            

               
        
        for gr_name in ['blocks','invers','reshuf','misj','transloc']:
            for entry in cont_blocks_dict[cont_name][gr_name]:
                result_dict[cont_name].append(entry)

        result_dict[cont_name]=sorted(result_dict[cont_name], key=lambda inter:inter[0], reverse=False)
       
    return result_dict, err_ref_cont_coord_errors_list

    

    




def FIND_REF_COORD_ERRORS(struct_dict,end_err_dict,unmapped_list,file_ref, file_contigs, uncovered_dict):
    

    ref_dict=READ_FASTA_ENTRY(file_ref)
    cont_dict=READ_FASTA_ENTRY(file_contigs)

    
    
    cont_blocks_dict,err_ref_cont_coord_errors_list=FIND_ERRORS_ALL_COORD(ref_dict,cont_dict, struct_dict, end_err_dict,unmapped_list)

    FIND_UNCOVERED_REF_FRAG(err_ref_cont_coord_errors_list,struct_dict,ref_dict,uncovered_dict)

    

    return err_ref_cont_coord_errors_list,cont_blocks_dict

    
