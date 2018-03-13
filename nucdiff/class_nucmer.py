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
import subprocess
import os
import getopt
import multiprocessing

from . import general
from . import class_errors
from . import class_interv_coord
from . import ref_coord




    
class Nucmer:
    def __init__(self, nuc_prefix, reference_file, contigs_file, working_dir,delta_file,coord_file):
        self.prefix=nuc_prefix
        
        self.ref=reference_file
        self.cont=contigs_file
        self.working_dir=working_dir
        
        if delta_file=='':
            self.delta_file='unknown'
        else:
            self.delta_file=delta_file

        if coord_file=='':
            self.coord_file='unknown'
        else:
            self.coord_file=coord_file
        self.intervals={}



    def RUN(self,nucmer_opt,filter_opt):

       
        if self.delta_file=='unknown':

            if '--maxmatch' in nucmer_opt:
                call_line='nucmer '+nucmer_opt+' --prefix='+self.prefix+' '+self.ref+' '+self.cont
            else:
                call_line='nucmer '+nucmer_opt+' --maxmatch --prefix='+self.prefix+' '+self.ref+' '+self.cont

            try:
                subprocess.check_call(call_line, shell=True)
            except subprocess.CalledProcessError:
                sys.exit('\nCould not run nucmer. Please, check nucmer input parameters values.')
            

            self.delta_file=self.prefix+'.delta'

           

       
                       
        f=open(self.prefix+'.filter','w')
        call_line_list=[]
        call_line_list.append('delta-filter')
        for entry in filter_opt.split(' '):
            call_line_list.append(entry)
        if not '-q' in filter_opt:
            call_line_list.append('-q')
                
        call_line_list.append(self.prefix+'.delta')

        try:
            subprocess.check_call(call_line_list,stdout=f)
        except subprocess.CalledProcessError:
            sys.exit('\nCould not run delta-filter. Please, check delta-filter input parameters values.')        
        
        f.close()

                

        f=open(self.prefix+'.coords','w')
        subprocess.check_call(['show-coords', '-qcldH', self.prefix+'.filter'], stdout=f)
        f.close()

        self.coord_file=self.prefix+'.coords'

        
        
    
    def PARSE(self):

        
        intervals_temp_dict={}

        f=open(self.coord_file)
        lines=f.readlines()
        f.close()
        
    
        
        if os.stat(self.coord_file)[6]!=0:
            
         snps_errors_list=[]
         start_pos=0
         for i in range(len(lines)):
            temp=lines[i].split()

            ref_start=int(temp[0])
            ref_end=int(temp[1])
            cont_start=int(temp[3])
            cont_end=int(temp[4])
            len_ref=int(temp[11])
            len_cont=int(temp[12])
            ref_direction=int(temp[17])
            cont_direction=int(temp[18])
            ref_name=temp[19]
            cont_name=temp[20]

            proc_identity=float(temp[9])
            

            if cont_name not in intervals_temp_dict:
               intervals_temp_dict[cont_name]=[]

            if cont_direction==1 and ref_direction==1:
                intervals_temp_dict[cont_name].append([cont_start, cont_end,ref_start,ref_end,cont_direction,ref_direction, len_cont, len_ref, ref_name,0,[]])

                 
                
            elif cont_direction==-1 and ref_direction==1:
                intervals_temp_dict[cont_name].append([cont_end,cont_start,ref_start,ref_end,cont_direction,ref_direction, len_cont, len_ref, ref_name,0,[]])

            else:
                sys.exit('Error: unknown case during parsing coord_file')

            if proc_identity!=100:
                intervals_temp_dict[cont_name][-1][10].append(lines[i]) 


            
    


        for cont_name in list(intervals_temp_dict.keys()):
            sort_intervals=sorted(intervals_temp_dict[cont_name], key=lambda inter:inter[0], reverse=False)

            self.intervals[cont_name]=sort_intervals

    def FIND_ERR_INSIDE_FRAG(self,proc_num,file_contigs):

        frag_dict={}

        num_all=0
        for cont_name in list(self.intervals.keys()):
            frag_dict[cont_name]=[]
            
            for i in range(len(self.intervals[cont_name])):
                info=self.intervals[cont_name][i]
                frag_dict[cont_name].append([info[0], info[1], info[2], info[3], info[8], info[10],0,[]])
                num_all+=1
        
       
        FIND_SNPS(frag_dict, self.coord_file, self.delta_file,self.prefix, proc_num, file_contigs )

        contigs_dict, contig_seqs, contig_names, contig_full_names_dict=general.READ_FASTA_ENTRY(file_contigs)

        
        temp_dif_gap=[]            
    
        for cont_name in list(frag_dict.keys()):
            for i in range(len(frag_dict[cont_name])):
                if frag_dict[cont_name][i][:4]==self.intervals[cont_name][i][:4] and frag_dict[cont_name][i][4]==self.intervals[cont_name][i][8]:
                    start_seq=frag_dict[cont_name][i][0]
                    end_seq=frag_dict[cont_name][i][1]
                    seq=contigs_dict[cont_name][start_seq-1 :end_seq-1   +1]
    
                    gaps_list=FIND_GAPS(seq,start_seq,end_seq,frag_dict[cont_name][i][7])   

                    if self.intervals[cont_name][i][10]!=[]:
                        self.intervals[cont_name][i][10].pop(0)

                    for gap in gaps_list:
                       self.intervals[cont_name][i][10].append(gap)
                    

                    for err in frag_dict[cont_name][i][7]:
                        self.intervals[cont_name][i][10].append(err)

                    self.intervals[cont_name][i][10]=sorted(self.intervals[cont_name][i][10],key=lambda inter:inter[0], reverse=False)

                    self.intervals[cont_name][i][10]=MERGE_GAPS(self.intervals[cont_name][i][10])

                    

        print('The difference detection inside fragments step is complete')

    def REMOVE_REF_OVERLAP(self):

        temp_interv_list=[]
        temp_errors_list=[]
        remove_list=[]
        new_list=[]
        
        for cont_name in list(self.intervals.keys()):
            
            max_gr=0
            for entry in self.intervals[cont_name]:
                if entry[9][0]>max_gr:
                    max_gr=entry[9][0]

            for gr_num in range(max_gr+1):
                for i in range(len(self.intervals[cont_name])):
                    if self.intervals[cont_name][i][9][0]==gr_num:
                        temp_interv_list.append([self.intervals[cont_name][i][0],self.intervals[cont_name][i][1],self.intervals[cont_name][i][2],self.intervals[cont_name][i][3],self.intervals[cont_name][i][4],i,self.intervals[cont_name][i][10]])
                
                #print '\nnew_list'
                #for entry in temp_interv_list:
                #    print entry[:10]
                #raw_input('rfj')

                if len(temp_interv_list)>1:

                        temp_interv_list=sorted(temp_interv_list, key=lambda inter:inter[2], reverse=False)
                                   
                        for i in range(len(temp_interv_list)-1):

                            r1_end=temp_interv_list[i][3]
                            r2_st=temp_interv_list[i+1][2]

                            
                            if r1_end+1>r2_st: # reference fragments overlap
                                
                                c2_dir=temp_interv_list[i+1][4]
                                c2_st=temp_interv_list[i+1][0]
                                c2_end=temp_interv_list[i+1][1]

                                errors_2=temp_interv_list[i+1][6]

                                

                                r2_coord=r1_end

                                if c2_dir==1:
                                    c2_coord, last_err_end=class_interv_coord.FIND_CONT_COORD_FORWARD_START(r2_st, c2_st, r2_coord, errors_2,c2_end)

                                    if c2_coord<c2_end:

                                        temp_interv_list[i+1][0]=c2_coord+1
                                        temp_interv_list[i+1][2]=r2_coord+1

                                        self.intervals[cont_name][temp_interv_list[i+1][5]][0]=c2_coord+1
                                        self.intervals[cont_name][temp_interv_list[i+1][5]][2]=r2_coord+1

                                        for entry in errors_2:
                                            if entry[0]>last_err_end:
                                                temp_errors_list.append(entry)

                                        del temp_interv_list[i+1][6][:]
                                        del self.intervals[cont_name][temp_interv_list[i+1][5]][10][:]

                                        for entry in temp_errors_list:
                                            self.intervals[cont_name][temp_interv_list[i+1][5]][10].append(entry)
                                            
                                        del temp_errors_list[:]
                                        
                                    else:
                                        remove_list.append(temp_interv_list[i+1][5])
                                        
                                        
                                        
                                else: #c2_dir==-1
                                   
                                    
                                    c2_coord, last_err_end=class_interv_coord.FIND_CONT_COORD_REVERSE_END_SECOND(r2_st, c2_end, r2_coord, errors_2)

                                    if c2_coord>1:

                                        temp_interv_list[i+1][1]=c2_coord-1
                                        temp_interv_list[i+1][2]=r2_coord+1

                                        self.intervals[cont_name][temp_interv_list[i+1][5]][1]=c2_coord-1
                                        self.intervals[cont_name][temp_interv_list[i+1][5]][2]=r2_coord+1

                                        
                                        for entry in errors_2:
                                            if entry[0]<last_err_end:
                                                temp_errors_list.append(entry)

                                        
                                        del temp_interv_list[i+1][6][:]
                                        del self.intervals[cont_name][temp_interv_list[i+1][5]][10][:]

                                        for entry in temp_errors_list:
                                            self.intervals[cont_name][temp_interv_list[i+1][5]][10].append(entry)


                                        
                                        del temp_errors_list[:]

                                        
                                        
                                    else:
                                       remove_list.append(temp_interv_list[i+1][5]) 
                                        

                        
                        for j in range(len(self.intervals[cont_name])):
                            if not j in remove_list:
                                new_list.append(self.intervals[cont_name][j])

                        del remove_list[:]
                        

                        del self.intervals[cont_name][:]

                        for entry in new_list:
                             self.intervals[cont_name].append(entry)

                        del new_list[:]

                        
                        
                                
                        

            
                    
                    

                del temp_interv_list[:]

            
            
    
    def FIND_STRUCTURAL_ERRORS(self,contigs_dict, ref_dict,reloc_dist):
        
        errors_list=[]

        
        

        #-------1.find translocations
        structure_dict={}
        #end_err_dict={}

        #1.1. redefine num element values
        self.FIND_REF_ORDER_NUM()
        
        self.REMOVE_REF_OVERLAP()
        
        x=self.IMPROVE_PARSING_RESULTS(contigs_dict, ref_dict)

        
        for cont_name in list(self.intervals.keys()):

            

            #1.2 sort intervals by translocation groups (fill in 'blocks' field)
            structure_dict[cont_name]={}

            #print cont_name
            #for entry in self.intervals[cont_name]:
            #    print entry[:10]
            #raw_input('rfj')


            
            

            
            cur_transl_group_num=0
            cur_transl_group_name=0
            structure_dict[cont_name][cur_transl_group_name]={'blocks':[self.intervals[cont_name][0]], 'output_line':[],'reason':[],'temp':-1}
            

            
            
            for i in range(len(self.intervals[cont_name])-1):
                    if self.intervals[cont_name][i][9][0]==self.intervals[cont_name][i+1][9][0]: #from one ref seq
                        structure_dict[cont_name][cur_transl_group_name]['blocks'].append(self.intervals[cont_name][i+1])
                    else:
                        structure_dict[cont_name][cur_transl_group_name]['temp']=i
                        cur_transl_group_num+=1
                        cur_transl_group_name=cur_transl_group_num
                        structure_dict[cont_name][cur_transl_group_name]={'blocks':[self.intervals[cont_name][i+1]], 'output_line':[], 'reason':[], 'temp':-1}
                        
           
            if len(list(structure_dict[cont_name].keys()))>1: #if there are translocation differences
                for transl_group in list(structure_dict[cont_name].keys()):

                    #1.2 fill in the 'output_line' field
                    c_st=structure_dict[cont_name][transl_group]['blocks'][0][0]
                    r_st=structure_dict[cont_name][transl_group]['blocks'][0][2]
                    r_name=structure_dict[cont_name][transl_group]['blocks'][0][8]

                    c_end=structure_dict[cont_name][transl_group]['blocks'][-1][1]
                    r_end=structure_dict[cont_name][transl_group]['blocks'][-1][3]

                    structure_dict[cont_name][transl_group]['output_line']=[c_st,c_end,transl_group,c_end-c_st+1, r_name, r_st,r_end]

                    #1.3.  fill in the 'reason' field if it's necessary
                    if structure_dict[cont_name][transl_group]['temp']!=-1: 
                        contig_seq=contigs_dict[cont_name]
                        i=structure_dict[cont_name][transl_group]['temp']
                        c1_end=self.intervals[cont_name][i][1]
                        c2_st=self.intervals[cont_name][i+1][0]


                        if c1_end+1<c2_st:
                           type_gap=class_interv_coord.ANALYSE_SPACE_SIMB(contig_seq, c1_end+1, c2_st-1)
                           
                           if type_gap=='gap':
                                structure_dict[cont_name][transl_group]['reason'].append([c1_end+1, c2_st-1, 'translocation-wrong_scaffolding',c2_st-c1_end-1, 'transl' ])

                           elif type_gap=='mixed_fragment':
                               structure_dict[cont_name][transl_group]['reason'].append([c1_end+1, c2_st-1, 'translocation-mixed_fragments',c2_st-c1_end-1, 'transl' ])
                           else:                               
                                structure_dict[cont_name][transl_group]['reason'].append([c1_end+1, c2_st-1, 'translocation-insertion',c2_st-c1_end-1, 'transl' ])

                        elif c1_end+1==c2_st:
                            structure_dict[cont_name][transl_group]['reason'].append([c1_end,c2_st,'translocation',0, 'transl'])

                        else:
                            structure_dict[cont_name][transl_group]['reason'].append([c2_st, c1_end,'translocation-overlap', c1_end-c2_st+1, 'transl'])
                            
                            r_end_1=self.intervals[cont_name][i][3]
                            r_st_1=self.intervals[cont_name][i][2]
                                
                            if self.intervals[cont_name][i][4]==1:
                                ref_overl_1=class_interv_coord.FIND_REF_COORD_FORWARD_END_1(r_end_1, c1_end, c2_st, self.intervals[cont_name][i][10], r_st_1)
                                structure_dict[cont_name][transl_group]['reason'][-1].append([self.intervals[cont_name][i][8],ref_overl_1, r_end_1])
                            else:
                                ref_overl_1=class_interv_coord.FIND_REF_COORD_REVERSE_END_1(r_st_1, c1_end, c2_st, self.intervals[cont_name][i][10],r_end_1)
                                structure_dict[cont_name][transl_group]['reason'][-1].append([self.intervals[cont_name][i][8],r_st_1, ref_overl_1])

                            r_end_2=self.intervals[cont_name][i+1][3]
                            r_st_2=self.intervals[cont_name][i+1][2]

                            for err_ent in self.intervals[cont_name][i+1][10]:
                                errors_list.append([err_ent[0],err_ent[1],err_ent[2],err_ent[3],err_ent[4]])
                                
                            if self.intervals[cont_name][i+1][4]==1:
                                
                                ref_overl_2,x=class_interv_coord.FIND_REF_COORD_FORWARD_START(r_st_2, c2_st, c1_end, errors_list)
                                ref_overl_2=min(r_end_2,ref_overl_2)
                                structure_dict[cont_name][transl_group]['reason'][-1].append([self.intervals[cont_name][i+1][8],r_st_2, ref_overl_2])
                                
                            else:
                                ref_overl_2,x=class_interv_coord.FIND_REF_COORD_REVERSE_END(r_end_2, c2_st, c1_end, errors_list)
                                ref_overl_2=max(r_st_2,ref_overl_2)
                                structure_dict[cont_name][transl_group]['reason'][-1].append([self.intervals[cont_name][i+1][8],ref_overl_2, r_end_2])
                                
                            for i in range(len(errors_list)):
                                errors_list.pop(0)
                                
                            
                    #1.4. delete the 'temp' field
                    structure_dict[cont_name][transl_group].pop('temp', None)
            else:
                structure_dict[cont_name][0].pop('temp', None)


        

       
        
        #---------2.find misjoins and circular genome ends
                    
        #2. find misjoin_groups inside translocation_groups                
        misj_groups_dict={}
        temp_errors_misj=[]

        
        for cont_name in list(structure_dict.keys()):
            cur_misj_num=0
            
            
            
            for transl_group in  sorted(structure_dict[cont_name]):
                #2.1 find misjoin_groups. Delete temp variable temp_group num
                new_misj_groups=FIND_MISJOIN_GROUP(structure_dict[cont_name][transl_group]['blocks'],reloc_dist)

                

                el_num=-1
                misj_groups_dict[cur_misj_num]={'blocks':[], 'output_line':[],'reason':[],'temp':-1}

                for i in range(len(new_misj_groups[0])):
                        misj_groups_dict[cur_misj_num]['blocks'].append(new_misj_groups[0][i])
                        el_num+=1

                
                for j in range(1,len(new_misj_groups)):
                    misj_groups_dict[cur_misj_num]['temp']=el_num
                    cur_misj_num+=1
                    misj_groups_dict[cur_misj_num]={'blocks':[], 'output_line':[],'reason':[],'temp':-1}
                        
                    
                    for i in range(len(new_misj_groups[j])):
                        misj_groups_dict[cur_misj_num]['blocks'].append(new_misj_groups[j][i])
                        el_num+=1


                if len(list(misj_groups_dict.keys()))>1: #there are misjoin differences
                    for misj_group in list(misj_groups_dict.keys()):
                        #2.2 fill in the 'output_line' field

                        c_st=misj_groups_dict[misj_group]['blocks'][0][0]
                        r_st=misj_groups_dict[misj_group]['blocks'][0][2]
                        r_name=misj_groups_dict[misj_group]['blocks'][0][8]

                        c_end=misj_groups_dict[misj_group]['blocks'][-1][1]
                        r_end=misj_groups_dict[misj_group]['blocks'][-1][3]

                        misj_groups_dict[misj_group]['output_line']=[c_st,c_end,misj_group, c_end-c_st+1, r_name, r_st,r_end]

                        #2.3 fill in the 'reason' field if necessary
                        if  misj_groups_dict[misj_group]['temp']!=-1:
                            i=misj_groups_dict[misj_group]['temp']

                            
                            a=class_interv_coord.Interv_coord(structure_dict[cont_name][transl_group]['blocks'][i])
                            b=class_interv_coord.Interv_coord(structure_dict[cont_name][transl_group]['blocks'][i+1])

                            
                            interv_case, interv_type=a.FIND_CASE(b,reloc_dist)

                           
                            r1_st=structure_dict[cont_name][transl_group]['blocks'][i][2]
                            r2_end=structure_dict[cont_name][transl_group]['blocks'][i+1][3]
                            len_r2=structure_dict[cont_name][transl_group]['blocks'][i+1][7]

                            r2_st=structure_dict[cont_name][transl_group]['blocks'][i+1][2]
                            r1_end=structure_dict[cont_name][transl_group]['blocks'][i][3]
                            len_r1=structure_dict[cont_name][transl_group]['blocks'][i][7]

                            c1_end=structure_dict[cont_name][transl_group]['blocks'][i][1]
                            c2_st=structure_dict[cont_name][transl_group]['blocks'][i+1][0]
                            c2_end=structure_dict[cont_name][transl_group]['blocks'][i+1][1]

                                


                            if (interv_case in ['4.1','8.1','12.1']) and r1_st==1 and r2_end==len_r2: #circular genome start
                                if interv_case=='4.1':
                                    c_space_len=c2_st-c1_end-1
                                    contig_seq=contigs_dict[cont_name]

                                    errors_list_1=class_interv_coord.INSERTION_INSIDE_CONTIG(c1_end+1, c2_st-1,contig_seq, c_space_len,interv_case)
                                    for entry in errors_list_1:
                                        misj_groups_dict[misj_group]['reason'].append(entry)
                                    del errors_list_1[:]
                                    
                                    misj_groups_dict[misj_group]['reason'].append([c1_end, c1_end, 'circular_genome_start', 1, interv_case] )
               
                            
                                elif interv_case=='8.1':
                                     misj_groups_dict[misj_group]['reason'].append([c1_end, c1_end, 'circular_genome_start', 1, interv_case])
               
                            
                                elif interv_case=='12.1':
                                    misj_groups_dict[misj_group]['reason'].append([c1_end, c1_end, 'circular_genome_start', 1, interv_case])
               
                                    c_len=c1_end-c2_st+1

                                    
                                    misj_groups_dict[misj_group]['reason'].append([c1_end, c1_end, 'deletion-collapsed_repeat', c_len, interv_case])

                                    corresp_ref_coord, last_err_end=class_interv_coord.FIND_REF_COORD_REVERSE_END(r2_end, c2_st, c1_end+1, structure_dict[cont_name][transl_group]['blocks'][i+1][10])
                                    first_base=last_err_end+1

                                    for entry in structure_dict[cont_name][transl_group]['blocks'][i+1][10]:
                                        if entry[0]>=first_base:  
                                            temp_errors_misj.append(entry)

                                    if structure_dict[cont_name][transl_group]['blocks'][i+1][10]!=[]:
                                        if structure_dict[cont_name][transl_group]['blocks'][i+1][10]!=temp_errors_misj:
                                            del structure_dict[cont_name][transl_group]['blocks'][i+1][10][:]
                                            for entry in temp_errors_misj:
                                                structure_dict[cont_name][transl_group]['blocks'][i+1][10].append(entry)

                                            del temp_errors_misj[:]
                                        
                                    structure_dict[cont_name][transl_group]['blocks'][i+1][0]=structure_dict[cont_name][transl_group]['blocks'][i][1]+1
                                    structure_dict[cont_name][transl_group]['blocks'][i+1][3]=corresp_ref_coord

                                    
                                else:
                                    print('ERROR: need to describe this case (circular genome): '+interv_case)
                                    print()

                            elif (interv_case in ['1.11', '5.11', '9.11']) and r2_st==1 and r1_end ==len_r1: #circular genome start 


                                if interv_case=='1.11':
                                    c_space_len=c2_st-c1_end-1
                                    contig_seq=contigs_dict[cont_name]
                                    
                                    errors_list_1=class_interv_coord.INSERTION_INSIDE_CONTIG(c1_end+1, c2_st-1,contig_seq, c_space_len,interv_case)
                                    for entry in errors_list_1:
                                        misj_groups_dict[misj_group]['reason'].append(entry)
                                    del errors_list_1[:]
                                    
                                    misj_groups_dict[misj_group]['reason'].append([c2_st, c2_st, 'circular_genome_start', 1, interv_case])

                                elif interv_case=='5.11':
                                    misj_groups_dict[misj_group]['reason'].append([c2_st, c2_st, 'circular_genome_start', 1, interv_case])

                            
                                elif interv_case=='9.11':
                                    c_space_len=c1_end-c2_st+1
                                    misj_groups_dict[misj_group]['reason'].append([c2_st-1, c2_st-1,'deletion-collapsed_repeat',c_space_len,interv_case])

                                    misj_groups_dict[misj_group]['reason'].append([c2_st, c2_st, 'circular_genome_start', 1, interv_case])

                                    corresp_ref_coord, last_err_end=class_interv_coord.FIND_REF_COORD_FORWARD_END(r1_end, c1_end, c2_st-1,structure_dict[cont_name][transl_group]['blocks'][i][10])
                                    first_base=last_err_end-1

                                    for entry in structure_dict[cont_name][transl_group]['blocks'][i][10]:
                                        if entry[0]<=first_base: 
                                            temp_errors_misj.append(entry)

                                    if structure_dict[cont_name][transl_group]['blocks'][i][10]!=[]:
                                        if structure_dict[cont_name][transl_group]['blocks'][i][10]!=temp_errors_misj:
                                            del structure_dict[cont_name][transl_group]['blocks'][i][10][:]
                                            for entry in temp_errors_misj:
                                                structure_dict[cont_name][transl_group]['blocks'][i][10].append(entry)

                                            del temp_errors_misj[:]
                                        
                                    structure_dict[cont_name][transl_group]['blocks'][i][1]=structure_dict[cont_name][transl_group]['blocks'][i+1][0]-1
                                    structure_dict[cont_name][transl_group]['blocks'][i][3]=corresp_ref_coord

                                    misj_groups_dict[misj_group]['output_line'][1]=structure_dict[cont_name][transl_group]['blocks'][i+1][0]-1
                                    misj_groups_dict[misj_group]['output_line'][6]=corresp_ref_coord

                                    

                                else:
                                    print('ERROR: need to describe this case (circular genome): '+interv_case)
                                    print()

                            else: #real misjoin
                                if c1_end+1<c2_st:
                                   contig_seq=contigs_dict[cont_name]
                                   type_gap=class_interv_coord.ANALYSE_SPACE_SIMB(contig_seq, c1_end+1, c2_st-1)
                                   
                                   
                                   if type_gap=='gap':
                                        misj_groups_dict[misj_group]['reason'].append([c1_end+1, c2_st-1, 'misjoin-wrong_scaffolding',c2_st-c1_end-1, 'misj' ])
                                   elif type_gap=='mixed_fragment':
                                        misj_groups_dict[misj_group]['reason'].append([c1_end+1, c2_st-1, 'misjoin-mixed_fragments',c2_st-c1_end-1, 'misj' ])
                                 
                                   else:                               
                                         misj_groups_dict[misj_group]['reason'].append([c1_end+1, c2_st-1, 'misjoin-insertion',c2_st-c1_end-1, 'misj' ])

                                elif c1_end+1==c2_st:
                                     misj_groups_dict[misj_group]['reason'].append([c1_end,c2_st,'misjoin',0, 'misj'])

                                else:
                                    misj_groups_dict[misj_group]['reason'].append([c2_st, c1_end,'misjoin-overlap', c1_end-c2_st+1, 'misj'])

                                    r_end_1=structure_dict[cont_name][transl_group]['blocks'][i][3]
                                    r_st_1=structure_dict[cont_name][transl_group]['blocks'][i][2]
                                        
                                    if self.intervals[cont_name][i][4]==1:
                                        ref_overl_1=class_interv_coord.FIND_REF_COORD_FORWARD_END_1(r_end_1, c1_end, c2_st, self.intervals[cont_name][i][10], r_st_1)
                                        misj_groups_dict[misj_group]['reason'][-1].append([structure_dict[cont_name][transl_group]['blocks'][i][8],ref_overl_1, r_end_1])
                                    else:
                                        ref_overl_1=class_interv_coord.FIND_REF_COORD_REVERSE_END_1(r_st_1, c1_end, c2_st, self.intervals[cont_name][i][10],r_end_1)
                                        misj_groups_dict[misj_group]['reason'][-1].append([structure_dict[cont_name][transl_group]['blocks'][i][8],r_st_1, ref_overl_1])

                                    r_end_2=structure_dict[cont_name][transl_group]['blocks'][i+1][3]
                                    r_st_2=structure_dict[cont_name][transl_group]['blocks'][i+1][2]

                                    for err_ent in self.intervals[cont_name][i+1][10]:
                                        errors_list.append([err_ent[0],err_ent[1],err_ent[2],err_ent[3],err_ent[4]])
                                        
                                    if self.intervals[cont_name][i+1][4]==1:
                                        
                                        ref_overl_2,x=class_interv_coord.FIND_REF_COORD_FORWARD_START(r_st_2, c2_st, c1_end, errors_list)
                                        ref_overl_2=min(r_end_2,ref_overl_2)
                                        misj_groups_dict[misj_group]['reason'][-1].append([structure_dict[cont_name][transl_group]['blocks'][i+1][8],r_st_2, ref_overl_2])
                                        
                                    else:
                                        ref_overl_2,x=class_interv_coord.FIND_REF_COORD_REVERSE_END(r_end_2, c2_st, c1_end, errors_list)
                                        ref_overl_2=max(r_st_2,ref_overl_2)
                                        misj_groups_dict[misj_group]['reason'][-1].append([structure_dict[cont_name][transl_group]['blocks'][i+1][8],ref_overl_2, r_end_2])
                                        
                                    for i in range(len(errors_list)):
                                        errors_list.pop(0)
                                    

                        #2.4 delete the 'temp' field
                        misj_groups_dict[misj_group].pop('temp', None)

                       
                #2.5 substitute a list of intervals in transl_group['blocks'] by dict of misj_groups
                structure_dict[cont_name][transl_group].pop('blocks',None)
                structure_dict[cont_name][transl_group]['blocks']={}

                for misj_group in list(misj_groups_dict.keys()):
                    structure_dict[cont_name][transl_group]['blocks'][misj_group]={}

                    for field in list(misj_groups_dict[misj_group].keys()):
                        structure_dict[cont_name][transl_group]['blocks'][misj_group][field]=misj_groups_dict[misj_group][field]

                for key in list(misj_groups_dict.keys()):
                    misj_groups_dict.pop(key, None)


                 
       

        temp_block_list=[]

        for cont_name in list(structure_dict.keys()):
            for transl_group in list(structure_dict[cont_name].keys()):
                for misj_group in list(structure_dict[cont_name][transl_group]['blocks'].keys()):
                    if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])>1:
                        for block in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']:
                            temp_block_list.append(block)

                        tmp_res=ref_coord.MERGE_BLOCK_MSJ(temp_block_list)
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['bounds']=[tmp_res[0],tmp_res[1],tmp_res[2],tmp_res[3]]

                        for i in range(len(temp_block_list)):
                            temp_block_list.pop(0)


                    else:
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['bounds']=[]
                            

                    
                     
                        
                    

        
        #-------3. simplify blocks before detecting transposition and inversions
        temp_errors=[]
        remove_list=[]
        nested_frag_list=[]
        for cont_name in sorted(structure_dict.keys()):
            for transl_group in sorted(structure_dict[cont_name].keys()):
                for misj_group in sorted(structure_dict[cont_name][transl_group]['blocks'].keys()):

                    #3.1 take away overlaps between query fragments
                    if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])>1:
                        for i in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])-1):
                            c1_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][1]
                            c2_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][0]

                            if c1_end+1>c2_st: # query fragments overlap

                                c2_dir=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][4]
                                r2_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][2]
                                r2_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][3]

                                errors_2=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10]

                                c2_coord=c1_end

                                if c2_dir==1:
                                    r2_coord, last_err_end=class_interv_coord.FIND_REF_COORD_FORWARD_START(r2_st, c2_st, c2_coord, errors_2)

                                    

                                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][0]=c2_coord+1
                                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][2]=r2_coord+1

                                     
                                   

                                    
                                else: #c2_dir==-1
                                    r2_coord, last_err_end=class_interv_coord.FIND_REF_COORD_REVERSE_END(r2_end, c2_st, c2_coord, errors_2)

                                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][0]=c2_coord+1
                                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][3]=r2_coord-1

                                for entry in errors_2:
                                    if entry[0]>last_err_end:
                                        
                                            temp_errors.append(entry)
                                        
                                del structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10][:]

                                for entry in temp_errors:
                                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10].append(entry)

                                del temp_errors[:]


                                       
                    #3.2 delete nested reference fragments
                    if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])>1:
                        flag_rep=1
                        while flag_rep==1:
                            #redefine num values
                            num_entry=0
                            for entry in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']:
                                entry[9]=num_entry
                                num_entry+=1

                            temp_sorted=sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'], key=lambda inter:inter[2], reverse=False)
                            for i in range(len(temp_sorted)):
                                structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][temp_sorted[i][9]][9]=i


                            #detect nested reference fragments
                            flag_rep=0
                            for i in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])-1):
                                r1_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][2]
                                r1_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][3]

                                for j in range(i+1, len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])):
                                    r2_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][j][2]
                                    r2_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][j][3]


                                    #find nested ref fragment
                
                                    if r1_st<=r2_st and r1_end>=r2_end: # Ri contains or equal to Rj
                                         nested_frag_list.append([j,i])
                                         break
                                

                                    elif r1_st>=r2_st and r1_end<=r2_end: #Rj contains Ri
                                         nested_frag_list.append([i,j])
                                         break

                                if len(nested_frag_list)==1:
                                   
                                    flag_rep=1
                                    #delete all info about nested fragments
                                    j_what=nested_frag_list[0][0]
                                    
                                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].pop(j_what)
                                    nested_frag_list.pop(0)
                                    break

                         
                    #3.3 take away overlaps between reference fragments   
                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']=sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'], key=lambda inter:inter[2], reverse=False)
                    
                    if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])>1:
                        
                        for i in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])-1):
                            r1_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][3]
                            r2_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][2]

                            
                            if r1_end+1>r2_st: # reference fragments overlap
                                
                                c2_dir=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][4]
                                c2_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][0]
                                c2_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][1]

                                errors_2=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10]

                                

                                r2_coord=r1_end

                                if c2_dir==1:
                                    c2_coord, last_err_end=class_interv_coord.FIND_CONT_COORD_FORWARD_START(r2_st, c2_st, r2_coord, errors_2,c2_end)
                                    
                                    if c2_coord<c2_end:

                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][0]=c2_coord+1
                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][2]=r2_coord+1

                                        for entry in errors_2:
                                            if entry[0]>last_err_end:
                                                temp_errors.append(entry)

                                        del structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10][:]

                                        for entry in temp_errors:
                                            structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10].append(entry)

                                        del temp_errors[:]
                                        
                                    else:
                                        remove_list.append(i+1)
                                        
                                        
                                        
                                else: #c2_dir==-1
                                   
                                    
                                    c2_coord, last_err_end=class_interv_coord.FIND_CONT_COORD_REVERSE_END_SECOND(r2_st, c2_end, r2_coord, errors_2)

                                    
                                    if c2_coord>1:

                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][1]=c2_coord-1
                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][2]=r2_coord+1
                                        for entry in errors_2:
                                            if entry[0]<last_err_end:
                                                temp_errors.append(entry)

                                        del structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10][:]

                                        for entry in temp_errors:
                                            structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][10].append(entry)

                                        del temp_errors[:]
                                        
                                    else:
                                       remove_list.append(i+1) 
                                        

                                   

                                    
                    new_list=[]
                    for j in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])):
                        if not j in remove_list:
                            new_list.append(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][j])

                    del remove_list[:]
                    

                    del structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][:]

                    for entry in new_list:
                         structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].append(entry)

                    del new_list[:]

                    
                        
                                
                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']=sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'], key=lambda inter:inter[0], reverse=False)


        

        #--------4. detect transpositions and inversions
        blocks_info_dict={}
        
        for cont_name in sorted(structure_dict.keys()):
            transp_gr_num=0
            block_num=0
            
            for transl_group in sorted(structure_dict[cont_name].keys()):
                for misj_group in sorted(structure_dict[cont_name][transl_group]['blocks'].keys()):
                    #4.1. redefine num values
                    num_entry=0
                    for entry in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']:
                        entry[9]=num_entry
                        num_entry+=1

                    temp_sorted=sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'], key=lambda inter:inter[2], reverse=False)
                    for i in range(len(temp_sorted)):
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][temp_sorted[i][9]][9]=i

                    #4.2. check if query fragments have the same order as ref fragments
                    transp_flag=0
                    for i in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])):
                        if i!=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][9]:
                            transp_flag=1
                            break


                    #4.3. describe blocks info
                    if transp_flag==1:
                        gr_name='_gr_'+str(transp_gr_num)
                        transp_gr_num+=1 

                    if transp_flag==0 and len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])>1:
                        inv_flag=1
                    else:
                        inv_flag=0
                    for i in range(len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])):
                        
                        
                        block_name=block_num
                        block_num+=1

                        blocks_info_dict[block_name]={}
                        blocks_info_dict[block_name]['block']=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i]
                        blocks_info_dict[block_name]['transp_output']=[]
                        blocks_info_dict[block_name]['invers_output']=[]
                        blocks_info_dict[block_name]['between_output']=[]
                        blocks_info_dict[block_name]['local_output']=[]
                        


                        c_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][0]
                        c_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][1]
                        ref_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][8]
                        ref_st=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][2]
                        ref_end=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][3]
                        num=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][9]
                        c_dir=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][4]


                        

                        # find transposition differences
                        if transp_flag==1:
                            err_name='reshuffling-part_'+str(num+1)+gr_name
                            
                            blocks_info_dict[block_name]['transp_output'].append([c_st, c_end, err_name, c_end-c_st+1,'transp',ref_name, ref_st, ref_end])

                        
                        if len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])>1:
                                
                            #find inversion differences
                            if c_dir==-1:
                                blocks_info_dict[block_name]['invers_output'].append([c_st, c_end, 'inversion', c_end-c_st+1,'transp',ref_name, ref_st, ref_end])

                            else:
                                if inv_flag==1:
                                    blocks_info_dict[block_name]['invers_output'].append([c_st, c_end, 'forward', c_end-c_st+1,'transp',ref_name, ref_st, ref_end])

                            #find deletion differences based on gaps between reference fragments
                            for entry in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']:
                                if num+1==entry[9]:
                                    r_st_next=entry[2]

                                    if ref_end+1<r_st_next:
                                        #if c_dir==1:
                                            blocks_info_dict[block_name]['between_output'].append([c_end,c_end, 'deletion', r_st_next-ref_end-1,'transp',ref_name,ref_end+1,r_st_next-1 ])
                                        #else:
                                        #    blocks_info_dict[block_name]['between_output'].append([max(1,c_st-1),max(1,c_st-1), 'deletion', r_st_next-ref_end-1,'transp',ref_name,ref_end+1,r_st_next-1 ])


                                    

                         #find insertion differences based on gaps between query fragments
                        if i==0:
                                
                                if structure_dict[cont_name][transl_group]['blocks'][misj_group]['output_line']!=[]:
                                    st_misj=structure_dict[cont_name][transl_group]['blocks'][misj_group]['output_line'][0]
                                    if c_st>st_misj:
                                        
                                        #if c_dir==1:
                                        #    blocks_info_dict[block_name]['between_output'].append([st_misj,c_st-1, 'insertion-multiple_copy', c_st-1-st_misj+1,'transp',ref_name,max(1,ref_st-1),max(1,ref_st-1)])
                                        #else:
                                            blocks_info_dict[block_name]['between_output'].append([st_misj,c_st-1, 'insertion-multiple_copy', c_st-1-st_misj+1,'transp',ref_name,ref_end,ref_end])
                                            
                                else:
                                    if structure_dict[cont_name][transl_group]['output_line']!=[]:
                                        st_misj=structure_dict[cont_name][transl_group]['output_line'][0]

                                        if c_st>st_misj:
                                        
                                            #if c_dir==1:
                                            #    blocks_info_dict[block_name]['between_output'].append([st_misj,c_st-1, 'insertion-multiple_copy', c_st-1-st_misj+1,'transp',ref_name,max(1,ref_st-1),max(1,ref_st-1)])
                                            #else:
                                                blocks_info_dict[block_name]['between_output'].append([st_misj,c_st-1, 'insertion-multiple_copy', c_st-1-st_misj+1,'transp',ref_name,ref_end,ref_end])
                                                
                                    
   

                                    
                        if i!=len(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'])-1:
                                c_st_next=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i+1][0] 
                                
                                if c_end+1<c_st_next:
                                    
                                    seq=contigs_dict[cont_name]
                                    space_type=class_interv_coord.ANALYSE_SPACE_SIMB(seq, c_end+1, c_st_next-1)
                                    if space_type=='gap':
                                       #if c_dir==1:
                                            blocks_info_dict[block_name]['between_output'].append([c_end+1,c_st_next-1, 'wrong_gap', c_st_next-c_end-1,'transp',ref_name,ref_end,ref_end  ])
                                       #else:
                                       #     blocks_info_dict[block_name]['between_output'].append([c_end+1,c_st_next-1, 'wrong_gap', c_st_next-c_end-1,'transp',ref_name,max(1,ref_st-1),max(1,ref_st-1)])
                                       
                                   
                                    elif space_type=='nucleotides':
                                        #if c_dir==1:
                                            blocks_info_dict[block_name]['between_output'].append([c_end+1,c_st_next-1, 'insertion', c_st_next-c_end-1,'transp',ref_name,ref_end,ref_end  ])
                                        #else:
                                        #    blocks_info_dict[block_name]['between_output'].append([c_end+1,c_st_next-1, 'insertion', c_st_next-c_end-1,'transp',ref_name,max(1,ref_st-1),max(1,ref_st-1)])

                                        
                                    else: # space_type='mixed_fragment'
                                        
                                        errors_list=class_interv_coord.FIND_INSERTION_GAP_INTERVALS(seq, c_end+1, c_st_next-1,'transp')
                                        for entry in errors_list:
                                            #if c_dir==1:
                                                blocks_info_dict[block_name]['between_output'].append([entry[0],entry[1], entry[2], entry[1]-entry[0]+1,'transp',ref_name,ref_end,ref_end  ])
                                            #else:
                                            #    blocks_info_dict[block_name]['between_output'].append([entry[0],entry[1], entry[2], entry[1]-entry[0]+1,'transp',ref_name,max(1,ref_st-1),max(1,ref_st-1)])
                                            
                                        

                        else:
                                if structure_dict[cont_name][transl_group]['blocks'][misj_group]['output_line']!=[] :
                                    end_misj=structure_dict[cont_name][transl_group]['blocks'][misj_group]['output_line'][1]
                                    if c_end<end_misj:
                                        #if c_dir==1:
                                            blocks_info_dict[block_name]['between_output'].append([c_end+1,end_misj, 'insertion-multiple_copy', end_misj-c_end-1+1,'transp', ref_name, ref_end,ref_end ])
                                        #else:
                                        #    blocks_info_dict[block_name]['between_output'].append([c_end+1,end_misj, 'insertion-multiple_copy', end_misj-c_end-1+1,'transp', ref_name, max(1,ref_st-1),max(1,ref_st-1)])
                                        
                                else:
                                    if structure_dict[cont_name][transl_group]['output_line']!=[] :
                                        end_misj=structure_dict[cont_name][transl_group]['output_line'][1]
                                        if c_end<end_misj:
                                            #if c_dir==1:
                                                blocks_info_dict[block_name]['between_output'].append([c_end+1,end_misj, 'insertion-multiple_copy', end_misj-c_end-1+1,'transp', ref_name, ref_end,ref_end ])
                                            #else:
                                            #    blocks_info_dict[block_name]['between_output'].append([c_end+1,end_misj, 'insertion-multiple_copy', end_misj-c_end-1+1,'transp', ref_name, max(1,ref_st-1),max(1,ref_st-1)])
                                            

                            
                        
                                    
                            
                        #write down local differences in the local_output field
                        for entry in structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][i][10]:
                            blocks_info_dict[block_name]['local_output'].append(entry) 
                                
                            
                    #2.10. substitute a list of intervals in misj_group['blocks'] by dict of blocks_info_groups
                    structure_dict[cont_name][transl_group]['blocks'][misj_group].pop('blocks',None)
                    structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks']={}

                    for block_name in  list(blocks_info_dict.keys()):
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]={}

                        for field in list(blocks_info_dict[block_name].keys()):
                            structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name][field]=blocks_info_dict[block_name][field]

                    for block_name in  list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output']=sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output'], key=lambda inter:inter[0], reverse=False)
                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['local_output']=sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['local_output'], key=lambda inter:inter[0], reverse=False)

                    for key in list(blocks_info_dict.keys()):
                        blocks_info_dict.pop(key, None)

        for cont_name in list(structure_dict.keys()):
            
            for transl_group in list(structure_dict[cont_name].keys()):
                for misj_group in list(structure_dict[cont_name][transl_group]['blocks'].keys()):
                    if len(list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()))>1:
                        for block_name in list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                            block=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block']
                            temp_block_list.append(block)
                            
                        tmp_res=ref_coord.MERGE_BLOCK_MSJ(temp_block_list)
                        init_res=structure_dict[cont_name][transl_group]['blocks'][misj_group]['bounds']
                        if init_res!=[tmp_res[0],tmp_res[1],tmp_res[2],tmp_res[3]]:
                            if init_res[0]!=tmp_res[0]:
                                for block_name in list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                                    if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][0]==tmp_res[0]:
                                        if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][4]==1:
                                            ref_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][2]

                                        else:
                                            ref_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][3]
                                        ref_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][8]
                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output'].append([init_res[0], tmp_res[0]-1,'insertion-multiple_copy',tmp_res[0]-init_res[0],'transp2',ref_name,ref_pos-1,ref_pos-1])

                            if init_res[1]!=tmp_res[1]:
                                for block_name in list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                                    if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][1]==tmp_res[1]:
                                        if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][4]==1:
                                            ref_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][3]
                                        else:
                                            ref_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][2]
                                            
                                        ref_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][8]
                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output'].append([tmp_res[1]+1, init_res[1],'insertion-multiple_copy',init_res[1]-tmp_res[1],'transp2',ref_name,ref_pos-1,ref_pos-1])


                            if init_res[2]!=tmp_res[2]:
                                for block_name in list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                                    if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][2]==tmp_res[2]:
                                        if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][4]==1:
                                            cont_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][0]
                                            
                                        else:
                                            cont_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][1]
                                            
                                        ref_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][8]
                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output'].append([cont_pos-1, cont_pos-1,'deletion-collapsed_repeat',tmp_res[2]-init_res[2],'transp2',ref_name,init_res[2],tmp_res[2]-1])

                            if init_res[3]!=tmp_res[3]:
                                for block_name in list(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                                    if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][3]==tmp_res[3]:
                                        if structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][4]==1:
                                            cont_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][1]
                                            
                                        else:
                                            cont_pos=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][0]

                                        ref_name=structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block'][8]
                                        structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output'].append([cont_pos-1, cont_pos-1,'deletion-collapsed_repeat',init_res[3]-tmp_res[3],'transp2',ref_name,tmp_res[3]+1,init_res[3]])


                            
                        


                        for i in range(len(temp_block_list)):
                            temp_block_list.pop(0)

        
        '''        
        for cont_name in sorted(structure_dict.keys()):
            print_flag=0
            print '#-------------'
            print cont_name
            print

            for entry in self.intervals[cont_name]:
                print entry
            print

            for transl_group in sorted(structure_dict[cont_name].keys()):
                print transl_group
                print 'OUTPUT:', structure_dict[cont_name][transl_group]['output_line']
                print 'REASON:', structure_dict[cont_name][transl_group]['reason']
                print 'BLOCKS:'
                for misj_group in sorted(structure_dict[cont_name][transl_group]['blocks'].keys()):
                    print '\t', misj_group
                    print '\tOUTPUT:', structure_dict[cont_name][transl_group]['blocks'][misj_group]['output_line']
                    print '\tREASON:', structure_dict[cont_name][transl_group]['blocks'][misj_group]['reason']
                    print '\tBLOCKS:'
                    for block_name in sorted(structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'].keys()):
                        print '\t\t', block_name
                        print '\t\tTRANSP_OUTPUT',structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['transp_output']
                        print '\t\tINVERS_OUTPUT',structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['invers_output']
                        print '\t\tBETWEEN_OUTPUT',structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['between_output']
                        print '\t\tLOCAL_OUTPUT',structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['local_output']
                        print '\t\tBLOCK',structure_dict[cont_name][transl_group]['blocks'][misj_group]['blocks'][block_name]['block']
                        print
                    print
                print
           
        '''
        
        return structure_dict    
      
        
                   
        
    
    def FIND_NON_STRUCTURAL_ERRORS(self, contigs_dict, ref_dict,reloc_dist):
        

        local_case_dict={'1.2':0,'1.3':0,'1.4':0,'4.8':0,'4.9':0,'4.10':0,'5.2':0,'5.3':0,'5.4':0,'8.8':0,'8.9':0,'8.10':0,'9.2':0,'9.3':0,'9.4.1':0,'9.4.2':0,'9.4.3':0,'9.4.4':0,\
                     '12.8.1':0,'12.8.2':0,'12.8.3':0,'12.8.4':0,'12.9':0,'12.10':0}
    

    
        non_structural_errors=class_errors.Errors()
        new_intervals={}

        #1. assign a number to an interval that shows the order of this interval if entries were sorted by the reference start coordinate
        self.FIND_REF_ORDER_NUM()

        

        
        for cont_name in list(self.intervals.keys()):
            
            
            new_intervals[cont_name]=[]
            contig_sequence=contigs_dict[cont_name]

            

            if len(self.intervals[cont_name])==1:
                new_intervals[cont_name].append(self.intervals[cont_name][0])
                
            else:
                
                # 2. for each entry pair search for non-structural differences and resolve them 
                cur_interv=[]

                
                
                for entry in self.intervals[cont_name][0]:
                    cur_interv.append(entry)
                    
                
                for i in range(1,len(self.intervals[cont_name])):
                    next_interv=self.intervals[cont_name][i]

                    
                    a=class_interv_coord.Interv_coord(cur_interv)
                    b=class_interv_coord.Interv_coord(next_interv)
                                   
                    interv_case, interv_type=a.FIND_CASE(b,reloc_dist)

                    if interv_type=='non_structural':
                        reference_seq=ref_dict[cur_interv[8]]

                        cur_interv, new_interv_case=a.FIND_ERROR_MERGE(b,interv_case, contig_sequence, reference_seq,cont_name)
                        

                        local_case_dict[new_interv_case]+=1

                        if cur_interv==[]:
                            print(self.intervals[cont_name])
                            sys.exit('ERROR: '+interv_case+' unknown case, no merging is provided')
                        

                    else:
                        temp_interv=[]
                        for entry in cur_interv:
                            temp_interv.append(entry)

                        new_intervals[cont_name].append(temp_interv)

                        
                        
                        for j in range(len(next_interv)):
                            cur_interv[j]=next_interv[j]
                            

                new_intervals[cont_name].append(cur_interv)
                

        for cont_name in list(self.intervals.keys()):
            self.intervals[cont_name]=new_intervals[cont_name]

        

        
    def FIND_REF_ORDER_NUM(self):
        
        ref_interv_dict={}
        
        for cont_name in list(self.intervals.keys()):

            
            for i in range(len(self.intervals[cont_name])):

                #1. a temp value of query fragment order
                self.intervals[cont_name][i].append(i)

                #2. find reference seq names and their corresponding  intervals
                if self.intervals[cont_name][i][8] not in ref_interv_dict:
                    ref_interv_dict[self.intervals[cont_name][i][8]]=[self.intervals[cont_name][i]]
                else:
                    ref_interv_dict[self.intervals[cont_name][i][8]].append(self.intervals[cont_name][i])

                

            #3. find order of fragments on ref sequences
            ref_num=-1
            for ref_name in list(ref_interv_dict.keys()):
                ref_num+=1
                
                ref_interv_dict[ref_name]=sorted(ref_interv_dict[ref_name], key=lambda inter:inter[2], reverse=False)

                #4.make the correspondence between ref order and query order for self.intervals elements

                for i in range(len(ref_interv_dict[ref_name])):
                    ref_order_num=[ref_num,i]

                    
                    self_interv_num=ref_interv_dict[ref_name][i][11]

                    self.intervals[cont_name][self_interv_num][9]=ref_order_num


            #5. delete all elements from ref_interv_dict
            for ref_name in list(ref_interv_dict.keys()):
                ref_interv_dict.pop(ref_name, None)

        

        #6. delete temporal element from self.intervals elements
        for cont_name in list(self.intervals.keys()):
            for i in range(len(self.intervals[cont_name])):
                self.intervals[cont_name][i].pop(11)


        
        



        
    

    def IMPROVE_PARSING_RESULTS(self, contigs_dict, ref_dict):
        
        for cont_name in list(self.intervals.keys()):
           self.intervals[cont_name]=sorted(self.intervals[cont_name], key=lambda inter:inter[0], reverse=False)


        
        #1. assign a number to an interval that shows the order of this interval if the entries were sorted by a reference start coordinate
        
        self.FIND_REF_ORDER_NUM()

        #2. delete or modify elements with nested query coordinates

        nested_frag_list=[]
        nested_transloc_list=[]

        count_dict={}

        end_err_dict={}
        
        for cont_name in list(self.intervals.keys()):
            end_err_dict[cont_name]={'wrong_beginning':[], 'wrong_end':[],'duplication':[]}

            if self.intervals[cont_name][0][0]>1:
                    end_err_dict[cont_name]['wrong_beginning'].append([1,self.intervals[cont_name][0][0]-1,'wrong_beginning',self.intervals[cont_name][0][0]-1,'ends'])

            if self.intervals[cont_name][-1][1]<self.intervals[cont_name][-1][6]:
                    end_err_dict[cont_name]['wrong_end'].append([self.intervals[cont_name][-1][1]+1,self.intervals[cont_name][-1][6],'wrong_end',self.intervals[cont_name][-1][6]-self.intervals[cont_name][-1][1],'ends'])
  
            
            if len(self.intervals[cont_name])>1:
                    
                flag_rep=1

                # repeat these steps until all query overlaps are resolved
                while flag_rep==1:

                    #1. detect query fragment overlaps
                    flag_rep=0

                    
                    #2. assign a number to an interval that shows the order of this interval if the entries were sorted by a reference start coordinate
                    self.FIND_REF_ORDER_NUM()
                    

                    for i in range(len(self.intervals[cont_name])-1):
                        c_st1=self.intervals[cont_name][i][0]
                        c_end1=self.intervals[cont_name][i][1]
                        r_st1=self.intervals[cont_name][i][2]
                        r_end1=self.intervals[cont_name][i][3]
                        num1=self.intervals[cont_name][i][9]
                        
                        
                        for j in range(i+1,len(self.intervals[cont_name])):
                            c_st2=self.intervals[cont_name][j][0]
                            c_end2=self.intervals[cont_name][j][1]
                            r_st2=self.intervals[cont_name][j][2]
                            r_end2=self.intervals[cont_name][j][3]
                            num2=self.intervals[cont_name][j][9]



                            if c_st1<=c_st2 and c_end1>=c_end2 and self.intervals[cont_name][i][8]==self.intervals[cont_name][j][8]:#C1 contains C2  or C1=C2
                                nested_frag_list.append([j,i])# j in i
                                break

                            elif c_st1>=c_st2 and c_end1<=c_end2 and self.intervals[cont_name][i][8]==self.intervals[cont_name][j][8]:#C2 contains C1:
                                nested_frag_list.append([i,j])
                                break

                            elif c_st1<=c_st2 and c_end1>=c_end2 and self.intervals[cont_name][i][8]!=self.intervals[cont_name][j][8]:#C1 contains C2  or C1=C2 , different ref
                                nested_transloc_list.append([j,i])
                                break
                                    
                            elif c_st1>=c_st2 and c_end1<=c_end2 and self.intervals[cont_name][i][8]!=self.intervals[cont_name][j][8]:#C2 contains C1, different ref
                                nested_transloc_list.append([i,j])
                                break

                            elif c_st2>c_end1:
                                break

                        if len(nested_transloc_list)==1:
                            

                            for kj in range(len(self.intervals[cont_name])):
                                if self.intervals[cont_name][kj][8] not in count_dict:
                                    count_dict[self.intervals[cont_name][kj][8]]=0
                                count_dict[self.intervals[cont_name][kj][8]]+=self.intervals[cont_name][kj][3]-self.intervals[cont_name][kj][2]+1

                           
                            
                            flag_rep=1

                            first_ref=self.intervals[cont_name][nested_transloc_list[0][0]][8]
                            second_ref=self.intervals[cont_name][nested_transloc_list[0][1]][8]

                            if count_dict[first_ref]>count_dict[second_ref]:
                                self.intervals[cont_name].pop(nested_transloc_list[0][1])
                            elif count_dict[first_ref]<count_dict[second_ref]:
                                self.intervals[cont_name].pop(nested_transloc_list[0][0])
                            else:
                                if self.intervals[cont_name][nested_transloc_list[0][0]][1]-self.intervals[cont_name][nested_transloc_list[0][0]][0]+1>=self.intervals[cont_name][nested_transloc_list[0][1]][1]-self.intervals[cont_name][nested_transloc_list[0][1]][0]+1:
                                    self.intervals[cont_name].pop(nested_transloc_list[0][1])
                                else:
                                    self.intervals[cont_name].pop(nested_transloc_list[0][0])

                            count_dict.clear()

                           
                            nested_transloc_list.pop(0)
                            break
                        
                        elif len(nested_frag_list)==1:
                            flag_rep=1

                            
                            #4. delete contained fragment
                            self.intervals[cont_name].pop(nested_frag_list[0][0])

                            #5. empty nested_frag_list
                            nested_frag_list.pop(0)

                            

                            break
                            
                            
                flag_rep=1

                # repeat these steps until all reference overlaps are resolved
                while flag_rep==1:

                    #1. assign a number to an interval that shows the order of this interval if the entries were sorted by a reference start coordinate
                    self.FIND_REF_ORDER_NUM()
                    
                    #2. detect reference fragments overlap
                    flag_rep=0

                    for i in range(len(self.intervals[cont_name])-1):
                        c_st1=self.intervals[cont_name][i][0]
                        c_end1=self.intervals[cont_name][i][1]
                        r_st1=self.intervals[cont_name][i][2]
                        r_end1=self.intervals[cont_name][i][3]
                        cont_dir1=self.intervals[cont_name][i][4]
                        num1=self.intervals[cont_name][i][9]
                        
                        
                        for j in range(i+1,len(self.intervals[cont_name])):
                            c_st2=self.intervals[cont_name][j][0]
                            c_end2=self.intervals[cont_name][j][1]
                            r_st2=self.intervals[cont_name][j][2]
                            r_end2=self.intervals[cont_name][j][3]
                            cont_dir2=self.intervals[cont_name][j][4]
                            num2=self.intervals[cont_name][j][9]
                    
                            #find nested ref fragment
                
                            if r_st1<=r_st2 and r_end1>=r_end2 and self.intervals[cont_name][i][8]==self.intervals[cont_name][j][8]: # Ri contains or equal to Rj
                                 nested_frag_list.append(j)

                                 break
                                

                            elif r_st1>=r_st2 and r_end1<=r_end2 and self.intervals[cont_name][i][8]==self.intervals[cont_name][j][8]: #Rj contains Ri
                                 nested_frag_list.append(i)
                                 break

                        if len(nested_frag_list)==1:
                                flag_rep=1
    
                                #delete all info about nested fragments

                                self.intervals[cont_name].pop(nested_frag_list[0])
                                nested_frag_list.pop(0)
                                
                                break

        
        for cont_name in list(self.intervals.keys()):
            self.intervals[cont_name]=sorted(self.intervals[cont_name], key=lambda inter:inter[0], reverse=False)

            if self.intervals[cont_name][0][0]>1:
                if end_err_dict[cont_name]['wrong_beginning']==[]:
                    err_st=1
                    end_err_dict[cont_name]['duplication'].append([err_st,self.intervals[cont_name][0][0]-1,'wrong_beginning',self.intervals[cont_name][0][0]-err_st,'ends'])

                else:
                    if self.intervals[cont_name][0][0]>end_err_dict[cont_name]['wrong_beginning'][0][1]+1:
                        err_st=end_err_dict[cont_name]['wrong_beginning'][0][1]+1
                        end_err_dict[cont_name]['duplication'].append([err_st,self.intervals[cont_name][0][0]-1,'wrong_beginning',self.intervals[cont_name][0][0]-err_st,'ends'])


                    
                
                        

            if self.intervals[cont_name][-1][1]<self.intervals[cont_name][-1][6]:
                    if end_err_dict[cont_name]['wrong_end']==[]:
                        err_end=self.intervals[cont_name][-1][6]
                        end_err_dict[cont_name]['duplication'].append([self.intervals[cont_name][-1][1]+1,err_end,'wrong_end',err_end-self.intervals[cont_name][-1][1],'ends'])

                    else:
                        if self.intervals[cont_name][-1][1]<end_err_dict[cont_name]['wrong_end'][0][0]-1:
                            err_end=end_err_dict[cont_name]['wrong_end'][0][0]-1
                            end_err_dict[cont_name]['duplication'].append([self.intervals[cont_name][-1][1]+1,err_end,'wrong_end',err_end-self.intervals[cont_name][-1][1],'ends'])


                    


        for cont_name in list(self.intervals.keys()):
            for i in range(len(self.intervals[cont_name])):
                if self.intervals[cont_name][i][10]!=[]:
                    self.intervals[cont_name][i][10]=sorted(self.intervals[cont_name][i][10], key=lambda inter:inter[0], reverse=False)



        return end_err_dict


    def FIND_UNMAPPED_CONTIGS(self, contigs_dict):
        overlap_errors=class_errors.Errors()
        unmapped_list=[]
        
        #1. find unmapped contigs
        
        for cont_name in list(contigs_dict.keys()):
            if cont_name not in self.intervals:
                unmapped_list.append(cont_name)

        return unmapped_list

            
            
    def FIND_WRONG_END(self,structure_dict):

        end_err_dict={}

        for cont_name in list(structure_dict.keys()):
            end_err_dict[cont_name]={'wrong_beginning':[], 'wrong_end':[]}

            first_entry=structure_dict[cont_name][0]['blocks'][0]['blocks'][0]['block']

            last_trl=sorted(structure_dict[cont_name].keys())[-1]
            last_msj=sorted(structure_dict[cont_name][last_trl]['blocks'].keys())[-1]
            last_blk=sorted(structure_dict[cont_name][last_trl]['blocks'][last_msj]['blocks'].keys())[-1]

            last_entry=structure_dict[cont_name][last_trl]['blocks'][last_msj]['blocks'][last_blk]['block']

            if first_entry[0]>1:
                for entry in [1,first_entry[0]-1,'wrong_beginning',first_entry[0]-1,'nuc12']:
                    end_err_dict[cont_name]['wrong_beginning'].append(entry)

            if last_entry[1]<last_entry[6]:
                for entry in [last_entry[1]+1,last_entry[6],'wrong_end',last_entry[6]-last_entry[1],'nuc12']:
                    end_err_dict[cont_name]['wrong_end'].append(entry)


        return end_err_dict

    def FIND_UNCOVERED_REF_REGIONS(self, ref_dict):

        intervals_dict={}
        uncov_dict={}

        for ref_name in list(ref_dict.keys()):
            intervals_dict[ref_name]=[]
            uncov_dict[ref_name]=[]

        for cont_name in list(self.intervals.keys()):
            for entry in self.intervals[cont_name]:
                ref_name=entry[8]

                intervals_dict[ref_name].append([entry[2],entry[3]])


        for ref_name in list(intervals_dict.keys()):
            if intervals_dict[ref_name]!=[]:
                intervals_dict[ref_name]=sorted(intervals_dict[ref_name],key=lambda inter:inter[0], reverse=False)
                overlap_list=FIND_OVERLAP_INTERVALS(intervals_dict[ref_name])
                uncov_dict[ref_name]=FIND_ZERO_COV(overlap_list,len(ref_dict[ref_name]))
        return uncov_dict
                
    def FIND_ERRORS(self, file_contigs, file_ref,reloc_dist):

        
        contigs_dict, contig_seqs, contig_names, contigs_full_names_dict=general.READ_FASTA_ENTRY(file_contigs)
        ref_dict, ref_seqs, ref_names, ref_full_names_dict=general.READ_FASTA_ENTRY(file_ref)

        
 
        #1. find unmapped query sequences
        unmapped_list=self.FIND_UNMAPPED_CONTIGS(contigs_dict)

        uncovered_dict={}#self.FIND_UNCOVERED_REF_REGIONS(ref_dict )

           
        #2. delete nested query and references fragments 
        end_err_dict=self.IMPROVE_PARSING_RESULTS(contigs_dict, ref_dict)

        #3. find local differences
        self.FIND_NON_STRUCTURAL_ERRORS(contigs_dict, ref_dict,reloc_dist)

        #4. check difference correctness
        '''
        for cont_name in self.intervals.keys():
            first_simb=self.intervals[cont_name][0][0] 
            if first_simb>1:
                   self.intervals[cont_name][0][10].append([1,first_simb-1,'wrong_beginning', first_simb-1-1+1, 'nuc12'])

                  

            last_simb=self.intervals[cont_name][-1][1]
            len_cont=self.intervals[cont_name][-1][6]

           
                
            if last_simb<len_cont:
                self.intervals[cont_name][-1][10].append([last_simb+1,len_cont,'wrong_end',len_cont-last_simb-1+1, 'nuc13' ])

                
            self.intervals[cont_name][0][10]=sorted(self.intervals[cont_name][0][10], key=lambda inter:inter[0], reverse=False)
            self.intervals[cont_name][-1][10]=sorted(self.intervals[cont_name][-1][10], key=lambda inter:inter[0], reverse=False)


        '''
        #5. find structural differences
        struct_dict=self.FIND_STRUCTURAL_ERRORS(contigs_dict, ref_dict,reloc_dist)

        #end_err_dict=self.FIND_WRONG_END(struct_dict)
        

        
        
        return struct_dict,end_err_dict,unmapped_list,uncovered_dict
        

#------------------------------------
def FIND_ZERO_COV(overlap_list,len_seq):

    uncov_list=[]

    if overlap_list==[]:
        uncov_list.append([1,len_seq])
    else:
        if overlap_list[0][0]>1:
            uncov_list.append([1,overlap_list[0][0]-1])

        for i in range(len(overlap_list)-1):
            entry_1=overlap_list[i]
            entry_2=overlap_list[i+1]

            if entry_2[0]>entry_1[1]+1:
                uncov_list.append([entry_1[1]+1,entry_2[0]-1])


        if overlap_list[-1][1]<len_seq:
            uncov_list.append([overlap_list[-1][1]+1,len_seq])

    return uncov_list
def FIND_OVERLAP_INTERVALS(interv_list):

    
    if interv_list==[]:
        result_list=[]

    else:
        if len(interv_list)==1:
            result_list=interv_list

        else:
            result_list=[]

            st_interv=interv_list[0][0]
            end_interv=interv_list[0][1]
            

            for entry in interv_list[1:]:
                st_new=entry[0]
                end_new=entry[1]

                if end_interv<st_new:
                    result_list.append([st_interv, end_interv])

                    st_interv=st_new
                    end_interv=end_new

                   
                else:
                    end_interv=max(end_interv,end_new)

            
            result_list.append([st_interv, end_interv])
   

    return result_list


def FIND_GAP_POS(seq, start_seq, end_seq):

    gap_list=[]
    
    st=-1
    end=-1

    for i in range(start_seq, end_seq+1):
        if seq[i] in 'Nn':
            
            if st==-1:
                st=i
                end=i

            else:
                end=i

        else:
            if st==-1:
                a='do nothing'

            else:
                gap_list.append([st,end])
                st=-1
                end=-1
            

    if st!=-1:
        gap_list.append([st,end])
   
        

    return gap_list    


def FIND_GAPS(seq,start_seq,end_seq,snps_errors_list):

    gap_list=[]
    gap_interv=[]

    #before diff
    if snps_errors_list==[]:
        gap_interv.append([0,len(seq)-1])


    else:
        if snps_errors_list[0][0]>start_seq:
            gap_interv.append([0,snps_errors_list[0][0]-start_seq-1])
            
            
        #after diff
        if snps_errors_list[-1][1]<end_seq:
            gap_interv.append([snps_errors_list[-1][1]-start_seq+1, end_seq-start_seq])
            
        #between diff
        if len(snps_errors_list)>1:

            for i in range(len(snps_errors_list)-1):
                interv_st=snps_errors_list[i][1]-start_seq+1
                interv_end=snps_errors_list[i+1][0]-start_seq-1

                if interv_st<=interv_end:
                    gap_interv.append([interv_st,interv_end])
                    

        #inside diff
        for  i in range(len(snps_errors_list)):
            cur_err=snps_errors_list[i]

            if cur_err[2]=='deletion':

                if i==0:
                  if i<len(snps_errors_list)-1:
                    next_err=snps_errors_list[i+1]
                    if next_err[0]>cur_err[1]:
                        gap_interv.append([cur_err[0]-start_seq, cur_err[0]-start_seq])
                    
                  else:
                      gap_interv.append([cur_err[0]-start_seq, cur_err[0]-start_seq])
                      

                elif i==len(snps_errors_list)-1:
                    if i!=0:
                        prev_err=snps_errors_list[i-1]
                        if prev_err[1]<cur_err[0]:
                            gap_interv.append([cur_err[0]-start_seq, cur_err[0]-start_seq])
                            
                else:
                    next_err=snps_errors_list[i+1]
                    prev_err=snps_errors_list[i-1]
                    if next_err[0]>cur_err[1] and prev_err[1]<cur_err[0]:
                        gap_interv.append([cur_err[0]-start_seq, cur_err[0]-start_seq])
    
    for interv in gap_interv:
        gap_cur_list=FIND_GAP_POS(seq, interv[0], interv[1])
        for entry in gap_cur_list:
               gap_list.append([entry[0]+start_seq, entry[1]+start_seq, 'gap', entry[1]-entry[0]+1, 'snps'])

    
                        
    gap_list=sorted(gap_list, key=lambda inter:inter[0], reverse=False)

            

    return gap_list        
                    



def FIND_MISJOIN_GROUP(group_list,reloc_dist):
    #2.1 redefine num value
    num_entry=0
    for entry in group_list:
        entry[9]=num_entry
        num_entry+=1

    
    temp_sorted=sorted(group_list, key=lambda inter:inter[2], reverse=False)
    for i in range(len(temp_sorted)):
        group_list[temp_sorted[i][9]][9]=i

    
    #2.2 sorted intervals in transl_group by ref coordinate
    group_list=sorted(group_list, key=lambda inter:inter[9], reverse=False)


    
    #2.3 Examine the distance between ref fragments. If they are less than reloc_dist, these entries are placed in one group

                
                
    cur_temp_group=0
    temp_group_list=[[group_list[0][9]]]
    
    for i in range(len(group_list)-1):
        if group_list[i+1][2]-group_list[i][3]<=reloc_dist:
            temp_group_list[cur_temp_group].append(group_list[i+1][9])
                        

        else:
            temp_group_list.append([group_list[i+1][9]])
            cur_temp_group+=1
                        

   
                       
    #2.4 sorted intervals by a query coord
    group_list=sorted(group_list, key=lambda inter:inter[0], reverse=False)
        
    #2.5. temporary add a new variable to the intervals showing the number of temp_group
    for entry in group_list:
        for i in range(len(temp_group_list)):
            if entry[9] in temp_group_list[i]:
                entry.append(i)
                break

    del temp_group_list[:]

                
    #2.6 find misjoin_groups. Delete temp variable temp_group num

    if cur_temp_group==0:
        result_gr_list=[]

        result_gr_list.append(group_list)

    else:
        cur_gr=0
        temp_group_list.append([group_list[0]])
        
        for i in range(len(group_list)-1):
            if group_list[i][11]==group_list[i+1][11]:
                temp_group_list[cur_gr].append(group_list[i+1])
            else:
                cur_gr+=1
                temp_group_list.append([group_list[i+1]])

        for misj_group in temp_group_list:
            for entry in misj_group:
                    entry.pop(11)

        
        result_gr_list=[]

        for misj_group in temp_group_list:
            new_gr_list=FIND_MISJOIN_GROUP(misj_group,reloc_dist)

            for entry in new_gr_list:
                result_gr_list.append(entry)

    return result_gr_list
    
        
    


def PARSE_SNPS_FILE(snps_file,snps_raw_dict):


    f=open(snps_file,'r')
    lines=f.readlines()

   
    for line in lines:
        temp=line[:-1].split()

        ref_pos=int(temp[0])
        ref_simb=temp[1]
        cont_simb=temp[2]
        cont_pos=int(temp[3])
        cont_dir=int(temp[9])
        ref_name=temp[10]
        cont_name=temp[11]

        
        snps_raw_dict[cont_name].append([cont_pos, ref_pos, cont_simb, ref_simb, ref_name, cont_dir])
        
           

    f.close()

    
   
   
    return snps_raw_dict

def CHECK_OVERL_FRAG(coord_line,coord_lines_list):
    cont_name=coord_line[0]
    c_st=coord_line[1][0]
    c_end=coord_line[1][1]
    r_st=coord_line[1][2]
    r_end=coord_line[1][3]
    r_name=coord_line[1][4]

    flag=0
    for entry in  coord_lines_list:
        if cont_name==entry[0] and r_name==entry[1][4]:
            if not (c_end<entry[1][0] or c_st>entry[1][1]):
                if not (r_end<entry[1][2] or r_st>entry[1][3]):
                    flag=1
                    
                    break
        
    return flag
                
def FIND_SNPS_TYPES(lines_error_list):
    snps_errors_list=[]
    snps_dict={'insertion':[], 'wrong_gap':[], 'gap':[], 'substitution':[],'deletion':[] }

    #find differences for each base
    for line in lines_error_list: 
        cont_simb=line[2]
        ref_simb=line[3]
        cont_pos=line[0]
        

        if cont_simb in 'ATGCatgc':
                if ref_simb in 'ATGCatgc':
                    error='substitution'

                elif ref_simb in 'NnQWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm':
                    error='no error'
                    
                elif ref_simb=='.':
                    error='insertion'

        elif cont_simb in 'QWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm':
               if ref_simb in 'ATGCatgcNnQWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm':
                    error='no error'
                    

               elif ref_simb=='.':
                    error='insertion'

        elif cont_simb in 'Nn':
                if ref_simb in 'ATGCatgcNnQWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm':
                    error='no error'

                elif ref_simb=='.':
                    error='wrong_gap'

        elif cont_simb=='.':
                if ref_simb in 'ATGCatgcNnQWERYUIOPSDFHJKLZXVBMqweryuiopsdfhjklzxvbm':
                    error='deletion'
                elif ref_simb=='.':
                    error='no error'
        else:
                print(ref_simb, cont_simb)
                sys.exit('ERROR: unknown case in snps file during parsing')

        if error!='no error':
                snps_dict[error].append(cont_pos)
                


    #merge differences into intervals
    for err in list(snps_dict.keys()):
        snps_dict[err]=sorted(snps_dict[err])

    for err in ['insertion', 'substitution', 'wrong_gap', 'gap']:
        if len(snps_dict[err])>1:
            cur_st=-1
            cur_end=-1
            cur_len=-1

            for el in snps_dict[err]:
                if cur_len==-1:
                    cur_st=el
                    cur_end=el
                    cur_len=1

                else:
                    if el==cur_end+1:
                        cur_end+=1
                        cur_len+=1
                    else:
                        snps_errors_list.append([cur_st,cur_end, err, cur_len, 'snps'])
                        cur_st=el
                        cur_end=el
                        cur_len=1

            snps_errors_list.append([cur_st,cur_end, err, cur_len, 'snps'])


        elif  len(snps_dict[err])==1:
            snps_errors_list.append([snps_dict[err][0],snps_dict[err][0], err, 1, 'snps'])

    del_list=[]
    if len(snps_dict['deletion'])==1:
        cont_st=snps_dict['deletion'][0]
        del_list.append([cont_st,1])
    elif len(snps_dict['deletion'])>1:
        cur_st=-1
        cur_len=-1

        for el in snps_dict['deletion']:
            if cur_len==-1:
                cur_st=el
                cur_len=1
            else:
                if el==cur_st:
                   cur_len+=1
                else:
                   del_list.append([cur_st,cur_len])
                   cur_st=el
                   cur_len=1
        del_list.append([cur_st,cur_len])
                
        
    new_err=[]
    for entry in del_list:
        el=entry[0]

        for i in range(len(snps_errors_list)):
            err=snps_errors_list[i]
            if err[0]<=el and err[1]>el:
                new_err.append([el+1,snps_errors_list[i][1],snps_errors_list[i][2],snps_errors_list[i][1]-el,'snps'])
                snps_errors_list[i][1]=el
                snps_errors_list[i][3]=el-snps_errors_list[i][0]+1

        snps_errors_list.append([entry[0],entry[0],'deletion',entry[1],'snps'])

    for el in new_err:
        snps_errors_list.append(el)

    snps_errors_list=sorted(snps_errors_list,key=lambda inter:inter[0], reverse=False)

    
    return snps_errors_list
                            
            


   
    
def FIND_SNPS_SINGLE(input_list):
    
    coord_lines_list=input_list[0]
    delta_file=input_list[1]
    prefix=input_list[2]
    unique_name=input_list[3]

    snp_file=prefix+'_'+str(unique_name)+'.snps'
    coord_file=prefix+'_'+str(unique_name)+'.coord'

    f=open(coord_file,'w')

    frag_dict={}
    
    for entry in coord_lines_list:
        cont_name=entry[0]
        frag_line=entry[1]

        f.write(frag_line[5][0])
        

        if cont_name not in frag_dict:
            frag_dict[cont_name]=[]
        frag_dict[cont_name].append(frag_line)
                
    f.close()

    
    f=open(snp_file, 'w')
    f_in=open(coord_file,'r')

    try:
        subprocess.check_call(['show-snps', '-SqHT', delta_file], stdin=f_in, stdout=f)
    except subprocess.CalledProcessError:
        return 'failed'
    
    f.close()
    f_in.close()

    
    
    snps_raw_dict={}
    for cont_name in list(frag_dict.keys()):
        snps_raw_dict[cont_name]=[]

    PARSE_SNPS_FILE(snp_file,snps_raw_dict)

    frag_num_list=[]
    for cont_name in list(snps_raw_dict.keys()):
        for snp in snps_raw_dict[cont_name]:
            for i in range(len(frag_dict[cont_name])):
                frg=frag_dict[cont_name][i]
                if frg[4]==snp[4] and snp[0]>=frg[0] and snp[0]<=frg[1]and snp[1]>=frg[2] and snp[1]<=frg[3]:
                    frag_num_list.append(i)
            
            if len(frag_num_list)>1:
                
                answ=CHECK_OVERL_FRAG([cont_name,frag_dict[cont_name][frag_num_list[0]]],[[cont_name,frag_dict[cont_name][frag_num_list[0]]]])
                
            elif len(frag_num_list)==0:
                a='do_nothing'
            else:
                frag_dict[cont_name][frag_num_list[0]][7].append(snp)
                

            for i in range(len(frag_num_list)):
                frag_num_list.pop(0)
                
    
    for cont_name in list(frag_dict.keys()):
        for i in range(len(frag_dict[cont_name])):
            frag_dict[cont_name][i][7]=FIND_SNPS_TYPES(frag_dict[cont_name][i][7])

    return frag_dict
        

def FIND_SNPS(frag_dict,coord_file, delta_file,prefix,proc_num, file_contigs):

    snp_file=prefix+'_filtered.snps'

    f=open(snp_file, 'w')
    f_in=open(coord_file,'r')

    
    f.close()
    f_in.close()

    input_list=[]
    raund_num=0
    flag=0
    coord_lines_list=[]
    temp_list=[]
    while flag!=1:

        ind=0
        for cont_name in list(frag_dict.keys()):
            
            ind+=1
            for i in range(len(frag_dict[cont_name])):
                

                if frag_dict[cont_name][i][6]==0:
                    if frag_dict[cont_name][i][5]==[]:
                        frag_dict[cont_name][i][6]=1
                    else:
                        answ=CHECK_OVERL_FRAG([cont_name,frag_dict[cont_name][i]],temp_list)
                        if answ==0:
                                temp_list.append([cont_name,frag_dict[cont_name][i]])
                                frag_dict[cont_name][i][6]=1

                        
            for entry in temp_list:
                coord_lines_list.append(entry)

            for j in range(len(temp_list)):
                temp_list.pop(0)

        raund_num+=1
        if coord_lines_list==[]:
            flag=1
        else:
            input_list.append([[],delta_file,prefix,raund_num])
            for entry in coord_lines_list:
               input_list[-1][0].append(entry) 
            

        

        for j in range(len(coord_lines_list)):
            coord_lines_list.pop(0)

   
    
    pool=multiprocessing.Pool(processes=proc_num)                
    output_list=pool.map(FIND_SNPS_SINGLE,input_list)

    for frag_dict_part in output_list:
        if frag_dict_part=='failed':
            import sys
            sys.exit(0)
        for cont_name in list(frag_dict_part.keys()):
            for new_entry in frag_dict_part[cont_name]:

                for i in range(len(frag_dict[cont_name])):
                    if new_entry[:5]==frag_dict[cont_name][i][:5]:

                        for err in new_entry[7]:
                            frag_dict[cont_name][i][7].append(err)
                        break

    
    for cont_name in list(frag_dict.keys()):
        for entry in frag_dict[cont_name]:
            if entry[6]==0:
                print('ERROR: 0 instead 1')
                print(cont_name, entry)
                

            

def MERGE_GAPS(err_list):

    new_err_list=[]

    if len(err_list)>2:
        for i in range(len(err_list)-2):

            if err_list[i][0]!=-1:
                if err_list[i][2]!='gap':
                    new_err_list.append(err_list[i])

                else:
                    if err_list[i+1][2]=='gap' and err_list[i][1]+1==err_list[i+1][0] and err_list[i+1][3]==1:
                        if err_list[i+2][2]=='deletion' and err_list[i+1][0]==err_list[i+2][0]:
                            err_list[i][1]+=1
                            err_list[i][3]+=1
                            err_list[i+1][0]=-1
                            new_err_list.append(err_list[i])
                        else:
                            new_err_list.append(err_list[i])
                    
                    else:
                        new_err_list.append(err_list[i])
        
        for i in range(len(err_list)-2,len(err_list)):
            if err_list[i][0]!=-1:
                 new_err_list.append(err_list[i])

        new_err_list=sorted(new_err_list, key=lambda inter:inter[0], reverse=False)
    else:
       new_err_list=err_list

    
        
        
    return new_err_list
             
            

    
                
