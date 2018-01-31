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

class Interv_coord:
    def __init__(self, list_info):
        self.c_st=list_info[0]
        self.c_end=list_info[1]
        self.r_st=list_info[2]
        self.r_end=list_info[3]
        self.c_dir=list_info[4]
        self.r_dir=list_info[5]
        self.c_len=list_info[6]
        self.r_len=list_info[7]
        self.r_name=list_info[8]
        self.num=list_info[9]
        self.errors=list_info[10]
        
    def PRINT_FUNC(self):
        entry=[self.c_st,self.c_end,self.r_st,self.r_end,self.c_dir,self.r_dir,self.c_len, self.r_len,self.r_name,self.num]
        print(entry)

    def FIND_CASE(self, b, reloc_dist):
        interv_case=''

        # 1. find query sequence fragment directions
        if self.c_dir==b.c_dir and self.c_dir==1: #(--->  ---> c1c2)
            interv_case+='1.' 
        elif self.c_dir!=b.c_dir and self.c_dir==1:  #(---->   <-----c1c2)
            interv_case+='2.'
        
        elif self.c_dir!=b.c_dir and self.c_dir==-1: #(<---  ---> c1c2)
            interv_case+='3.'

        elif self.c_dir==b.c_dir and self.c_dir==-1: #(<---  <--- c1c2)
            interv_case+='4.'
        

        # 2. placement of query sequence fragments
        if self.c_end+1<b.c_st: #(-----...-----c1c2)
            interv_case+='1'

        elif self.c_end+1==b.c_st: #(------|-----c1c2)
            interv_case+='2'
        elif  self.c_end+1>b.c_st and self.c_end<b.c_end: #(-----*-*-****** c1c2)
            interv_case+='3'
        else:
            print('ERROR: unknown contig case in Interv_coord.find_case()')
            print('\t',self.c_st, self.c_end,self.r_st, self.r_end,self.c_dir, self.c_len,self.r_dir, self.r_len)
            print('\t',b.c_st, b.c_end,b.r_st, b.r_end,b.c_dir, b.c_len,b.r_dir, b.r_len)
            print()
            return '0.0', 'structural'
           
            


        case_dict={'1.1':'1.','1.2':'5.','1.3':'9.','2.1':'2.','2.2':'6.','2.3':'10.','3.1':'3.','3.2':'7.','3.3':'11.','4.1':'4.','4.2':'8.','4.3':'12.'}

        interv_case=case_dict[interv_case]    
        

        # 3. placement of reference fragments
        if self.r_end+1<b.r_st: #(-----...-----r1r2)
            if b.r_st-self.r_end>=reloc_dist:
                 interv_case+='1'
            else:
                interv_case+='2'

        elif self.r_end+1==b.r_st: #(------|------ r1r2)
            interv_case+='3'

        elif self.r_end+1>b.r_st and self.r_end<b.r_end: #(-----*-*-****** r1r2)
            interv_case+='4'

        elif self.r_st==b.r_st and self.r_end<b.r_end: # ( |----          r1
            interv_case+='5'                           #   |------------ r2)

        elif self.r_st==b.r_st and self.r_end==b.r_end: # ( |------------| r1
            interv_case+='6'                            #   |------------| r2)

        elif self.r_st==b.r_st and self.r_end>b.r_end: # ( |------------ r1
            interv_case+='7'                           #   |------       r2)

        elif self.r_st-1<b.r_end and self.r_end>b.r_end and self.r_st>b.r_st: # (       ------------ r1
            interv_case+='8'                                                  #   ------------       r2)

        elif self.r_st-1==b.r_end : # ( |------------|------------| r2r1)
            interv_case+='9'

        elif self.r_st-1>b.r_end: # ( |------------...------- r2r1)
            if self.r_st-b.r_end<reloc_dist:
                interv_case+='10'
            else:
                interv_case+='11'

        elif self.r_st>b.r_st and self.r_end<b.r_end:    # (    ------      r1
            interv_case+='12'                            #   ------------ r2)

        elif self.r_st<b.r_st and self.r_end>b.r_end:    # (   ------------ r1
            interv_case+='13'                            #      --------    r2)

        elif self.r_st<b.r_st and self.r_end==b.r_end:  # (  |------------| r1
            interv_case+='14'                            #         -------| r2)

        elif self.r_st>b.r_st and self.r_end==b.r_end:   # (       -------| r1
            interv_case+='15'                            #   |------------| r2)

        else:
            
           print('ERROR: unknown reference case in Interv_coord.find_case()')
           print('\t',self.c_st, self.c_end,self.r_st, self.r_end,self.c_dir, self.c_len,self.r_dir, self.r_len)
           print('\t',b.c_st, b.c_end,b.r_st, b.r_end,b.c_dir, b.c_len,b.r_dir, b.r_len)
           return '0.0', 'structural'
            
        non_structural_list=['1.2','1.3','1.4','1.5','1.6','1.7','1.12','1.13','1.14', '1.15',\
                             '2.5','2.6','2.7','2.12','2.13','2.14','2.15',\
                             '3.5','3.6','3.7','3.12','3.13','3.14','3.15',\
                             '4.5','4.6','4.7','4.8','4.9','4.10','4.12','4.13','4.14','4.15',\
                             '5.2','5.3','5.4','5.5','5.6','5.7','5.12','5.13','5.14','5.15',\
                             '6.5','6.6','6.7','6.12','6.13','6.14','6.15',\
                             '7.5','7.6','7.7','7.12','7.13','7.14','7.15',\
                             '8.5','8.6','8.7','8.8','8.9','8.10','8.12','8.13','8.14','8.15',\
                             '9.2','9.3','9.4','9.5','9.6','9.7','9.12','9.13','9.14','9.15',\
                             '10.5','10.6','10.7','10.12','10.13','10.14','10.15',\
                             '11.5','11.6','11.7','11.12','11.13','11.14','11.15',\
                             '12.5','12.6','12.7','12.8','12.9','12.10','12.12','12.13','12.14','12.15']

        self_first_list=['1.2','1.3','1.4','1.5','1.6','1.7','1.12','1.13','1.14', '1.15',\
                         '5.2','5.3','5.4','5.5','5.6','5.7','5.12','5.13','5.14','5.15',\
                         '9.2','9.3','9.4','9.5','9.6','9.7','9.12','9.13','9.14','9.15',\
                         '2.5','2.6','2.7','2.13','2.14',\
                         '6.5','6.6','6.7','6.13','6.14',\
                         '10.5','10.6','10.7','10.13','10.14',\
                         '3.5','3.6','3.7','3.13','3.14',\
                         '7.5','7.6','7.7','7.13','7.14',\
                         '11.5','11.6','11.7','11.13','11.14',\
                         '4.5','4.6','4.7','4.13','4.14',\
                         '8.5','8.6','8.7','8.13','8.14',\
                         '12.5','12.6','12.7','12.13','12.14']
        
        b_first_list=['2.5','2.6','2.7','2.12','2.15',\
                      '6.5','6.6','6.7','6.12','6.15',\
                      '10.5','10.6','10.7','10.12','10.15',\
                      '3.5','3.6','3.7','3.12','3.15',\
                      '7.5','7.6','7.7','7.12','7.15',\
                      '11.5','11.6','11.7','11.12','11.15',\
                      '4.5','4.6','4.7','4.8','4.9','4.10','4.12','4.15',\
                      '8.5','8.6','8.7','8.8','8.9','8.10','8.12','8.15',\
                      '12.5','12.6','12.7','12.8','12.9','12.10','12.12','12.15']


        
                         
        

        if interv_case in non_structural_list: # non-structural difference

            if self.r_name==b.r_name: # from same ref sequence

                        # check if placement order is correct
                        if self.num[1]+1==b.num[1]:
                            if interv_case in self_first_list:
                                interv_type='non_structural'
                            else:
                                interv_type='structural'

                        elif b.num[1]+1==self.num[1]:
                            if interv_case in b_first_list:
                                interv_type='non_structural'
                            else:
                                interv_type='structural'
                        else:
                            interv_type='structural'
           
            else:
                interv_type='structural'
            
            
        else:
            interv_type='structural'

        return interv_case, interv_type


    def FIND_ERROR_MERGE(self,b,interv_case, contig_seq, ref_seq,cont_name):
        overlap_interv=[]

        if interv_case=='1.2':
            c_space_len=b.c_st-self.c_end-1
            ref_space_len=b.r_st-self.r_end-1

            errors_list=INSERTION_CONTIG_REFERENCE(self.c_end+1, b.c_st-1,contig_seq, c_space_len, ref_space_len, interv_case)

            for entry in errors_list:
                entry.append([])

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)

            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num, temp_errors]
            

        elif interv_case=='1.3':
            c_space_len=b.c_st-self.c_end-1

            errors_list=INSERTION_INSIDE_CONTIG(self.c_end+1, b.c_st-1,contig_seq, c_space_len,interv_case)

            for entry in errors_list:
                entry.append([])

            #print errors_list
            #raw_input('grk')

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)
                        
            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num, temp_errors]
            

        

        elif interv_case=='1.4':
            c_space_len=b.c_st-self.c_end-1
            ref_space_len=self.r_end-b.r_st+1

             
            errors_list=INSERTION_INSIDE_CONTIG(self.c_end+1, b.c_st-1,contig_seq, c_space_len, interv_case)
            for entry in errors_list:
                entry.append([])
            
           
            corresp_cont_coord, last_err_end=FIND_CONT_COORD_FORWARD_START(b.r_st, b.c_st, self.r_end, b.errors, b.c_end)

            if min(corresp_cont_coord,b.c_end)-b.c_st+1!=0:
                overlap_st=FIND_CONT_COORD_FORWARD_END_1(self.r_end, self.c_end, b.r_st, self.errors, self.c_st)
                overlap_interv.append([self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end,1,[b.c_st, min(corresp_cont_coord,b.c_end),'insertion-multiple_copy',min(corresp_cont_coord,b.c_end)-b.c_st+1, interv_case]])

                errors_list.append([b.c_st, min(corresp_cont_coord,b.c_end),'insertion-multiple_copy',min(corresp_cont_coord,b.c_end)-b.c_st+1, interv_case,overlap_interv])

            first_base=last_err_end+1


            
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)

            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

    
            '''            
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print self.errors
            print
            print overlap_st
            print
            print self.r_end, self.c_end, b.r_st
            print self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end
            

            
            print
            print
            print ref_seq[b.r_st-1: self.r_end]
            print
            print contig_seq[overlap_st-1: self.c_end]
            raw_input('jdg1.4')    
            '''
        
        elif interv_case=='4.8':
            c_space_len=b.c_st-self.c_end-1
            ref_space_len=b.r_end-self.r_st+1

            
            errors_list=INSERTION_INSIDE_CONTIG(self.c_end+1, b.c_st-1,contig_seq, c_space_len, interv_case)
            for entry in errors_list:
                entry.append([])
            

            corresp_cont_coord,last_err_end=FIND_CONT_COORD_REVERSE_END(b.r_end, b.c_st, self.r_st, b.errors, b.c_end)

            if min(corresp_cont_coord,b.c_end)-b.c_st+1!=0:
                overlap_st=FIND_CONT_COORD_REVERSE_END_1(self.r_st, self.c_end, b.r_end, self.errors, self.c_st)
                overlap_interv.append([self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end,-1,[b.c_st, min(corresp_cont_coord,b.c_end),'insertion-multiple_copy',min(corresp_cont_coord,b.c_end)-b.c_st+1, interv_case]])

                errors_list.append([b.c_st, min(corresp_cont_coord,b.c_end),'insertion-multiple_copy',min(corresp_cont_coord,b.c_end)-b.c_st+1, interv_case,overlap_interv])

            first_base=last_err_end+1


            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)
                
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print overlap_st
            print
            print self.errors
            print
            print self.r_st, self.c_end, b.r_end
            print self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end
            

            
            print
            print
            import general
            print general.COMPL_STRING(ref_seq[self.r_st-1: b.r_end])
            print
            print contig_seq[overlap_st-1: self.c_end]
            print
            print contig_seq[b.c_st-1: min(corresp_cont_coord,b.c_end)]
            raw_input('jdg4.8')    
            '''
            
        elif interv_case=='4.9':
            c_space_len=b.c_st-self.c_end-1

            errors_list=INSERTION_INSIDE_CONTIG(self.c_end+1, b.c_st-1,contig_seq, c_space_len,interv_case)

            for entry in errors_list:
                entry.append([])

            #print errors_list
            #raw_input('flk')

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)

                        
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors] 
            
        elif interv_case=='4.10':
            c_space_len=b.c_st-self.c_end-1
            ref_space_len=self.r_st-b.r_end-1

            errors_list=INSERTION_CONTIG_REFERENCE(self.c_end+1, b.c_st-1,contig_seq, c_space_len, ref_space_len,interv_case)

            for entry in errors_list:
                entry.append([])

            #print errors_list
            #raw_input('flk')
                
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)


            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num, temp_errors]
 
        elif interv_case=='5.2':
            ref_space_len=b.r_st-self.r_end-1

            errors_list=[[self.c_end, self.c_end, 'deletion', ref_space_len, interv_case,[]]]

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)


            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]


        elif interv_case=='5.3':
            errors_list=[[]]

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)

            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

        elif interv_case=='5.4':
            ref_space_len=self.r_end-b.r_st+1

            errors_list=[]

            corresp_cont_coord, last_err_end=FIND_CONT_COORD_FORWARD_START(b.r_st, b.c_st, self.r_end, b.errors, b.c_end)

            if min(corresp_cont_coord,b.c_end)-self.c_end-1+1!=0:
                overlap_st=FIND_CONT_COORD_FORWARD_END_1(self.r_end, self.c_end, b.r_st, self.errors,self.c_st)
                overlap_interv.append([self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end,1,[self.c_end+1, min(corresp_cont_coord,b.c_end), 'insertion-tandem_multiple_copy', min(corresp_cont_coord,b.c_end)-self.c_end-1+1, interv_case]])

                errors_list.append([self.c_end+1, min(corresp_cont_coord,b.c_end), 'insertion-tandem_multiple_copy', min(corresp_cont_coord,b.c_end)-self.c_end-1+1, interv_case,overlap_interv])

           
            first_base=last_err_end+1


               
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)


            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print self.errors
            print
            print overlap_st
            print
            print self.r_end, self.c_end, b.r_st
            print self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end
            

            
            print
            print
            print ref_seq[b.r_st-1: self.r_end]
            print
            print contig_seq[overlap_st-1: self.c_end]
            print
            print contig_seq[self.c_end+1-1: min(corresp_cont_coord,b.c_end)]
            raw_input('jdg5.4')    
            '''


        elif interv_case=='8.8':
            ref_space_len=b.r_end-self.r_st+1
            
            errors_list=[]

            corresp_cont_coord,last_err_end=FIND_CONT_COORD_REVERSE_END(b.r_end, b.c_st, self.r_st, b.errors, b.c_end)

            
            
            if min(corresp_cont_coord,b.c_end)-self.c_end-1+1!=0:
                overlap_st=FIND_CONT_COORD_REVERSE_END_1(self.r_st, self.c_end, b.r_end, self.errors,self.c_st)
                overlap_interv.append([self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end,-1,[self.c_end+1,min(corresp_cont_coord,b.c_end) , 'insertion-tandem_multiple_copy', min(corresp_cont_coord,b.c_end)-self.c_end-1+1, interv_case]])

                errors_list.append([self.c_end+1,min(corresp_cont_coord,b.c_end) , 'insertion-tandem_multiple_copy', min(corresp_cont_coord,b.c_end)-self.c_end-1+1, interv_case,overlap_interv])
            first_base=last_err_end+1

            
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)

            
              
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num, temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print overlap_st
            print
            print self.errors
            print
            print self.r_st, self.c_end, b.r_end
            print self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end
            

            
            print
            print
            import general
            print general.COMPL_STRING(ref_seq[self.r_st-1:b.r_end])
            print
            print contig_seq[overlap_st-1: self.c_end]
            print
            print contig_seq[self.c_end+1-1:min(corresp_cont_coord,b.c_end)]
            raw_input('jdg8.8')    
            '''

        elif interv_case=='8.9':
            errors_list=[[]]


            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num, temp_errors]

        elif interv_case=='8.10':
            ref_space_len=self.r_st-b.r_end-1

            errors_list=[[self.c_end, self.c_end, 'deletion', ref_space_len, interv_case,[]]]

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                temp_errors.append(entry)


            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num, temp_errors]


           
        elif interv_case=='9.2':
            c_space_len=self.c_end-b.c_st+1
            ref_space_len=b.r_st-self.r_end-1

            
            errors_list=[[self.c_end, self.c_end, 'deletion', ref_space_len, interv_case,[]]]
            corresp_ref_coord, last_err_end=FIND_REF_COORD_FORWARD_START(b.r_st, b.c_st, self.c_end, b.errors)

            

            if min(corresp_ref_coord,b.r_end)-b.r_st+1 !=0:
                overlap_st=FIND_REF_COORD_FORWARD_END_1(self.r_end, self.c_end, b.c_st, self.errors,self.r_st)
                overlap_interv.append([self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end, 1,[self.c_end, self.c_end, 'deletion-collapsed_repeat',min(corresp_ref_coord,b.r_end)-b.r_st+1,interv_case]])


                errors_list.append([self.c_end, self.c_end, 'deletion-collapsed_repeat',min(corresp_ref_coord,b.r_end)-b.r_st+1,interv_case,overlap_interv])
            

            first_base=last_err_end+1

            
            
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)

                        
            
            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print overlap_st
            print
            print self.errors
            print
            print self.r_end, self.c_end, b.c_st
            print self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end
            

            
            print
            print
            print ref_seq[overlap_st-1:self.r_end]
            print
            print contig_seq[b.c_st-1: self.c_end]
            print
            print ref_seq[b.r_st-1: min(corresp_ref_coord,b.r_end)]
            
            raw_input('jdg9.2')    
            '''



            
        elif interv_case=='9.3':
            c_space_len=self.c_end-b.c_st+1

            errors_list=[]

            corresp_ref_coord, last_err_end=FIND_REF_COORD_FORWARD_START(b.r_st, b.c_st, self.c_end, b.errors)

            
            if min(corresp_ref_coord,b.r_end)-b.r_st+1!=0:
                overlap_st=FIND_REF_COORD_FORWARD_END_1(self.r_end, self.c_end, b.c_st, self.errors,self.r_st)
                overlap_interv.append([self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end,1, [self.c_end, self.c_end, 'deletion-collapsed_tandem_repeat',min(corresp_ref_coord,b.r_end)-b.r_st+1,interv_case]])

                errors_list.append([self.c_end, self.c_end, 'deletion-collapsed_tandem_repeat',min(corresp_ref_coord,b.r_end)-b.r_st+1,interv_case,overlap_interv])
            
            first_base=last_err_end+1    

            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)

           
            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print overlap_st
            print
            print self.errors
            print
            print self.r_end, self.c_end, b.c_st
            print self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end
            

            
            print
            print
            print ref_seq[overlap_st-1:self.r_end]
            print
            print contig_seq[b.c_st-1: self.c_end]
            print
            print ref_seq[b.r_st-1: min(corresp_ref_coord,b.r_end)]
            raw_input('jdg9.3')    
            '''

        elif interv_case=='9.4':
            errors_list=[]

            b_errors_del=[]
            b_errors_ins=[]
            for entry in b.errors:
            
                b_errors_del.append([entry[0],entry[1],entry[2],entry[3],entry[4]])
                b_errors_ins.append([entry[0],entry[1],entry[2],entry[3],entry[4]])
                

            ins_end_coord_2, last_err_end_2=FIND_REF_COORD_FORWARD_START(b.r_st, b.c_st, self.c_end, b_errors_del)


            ins_end_coord, last_err_end=FIND_CONT_COORD_FORWARD_START(b.r_st, b.c_st, self.r_end, b_errors_ins,b.c_end)
            
            first_base=last_err_end+1


            
        
            if min(ins_end_coord,b.c_end)-self.c_end-1+1>0:
                    overlap_st=FIND_REF_COORD_FORWARD_END_1(self.r_end, self.c_end, b.c_st, self.errors,self.r_st)
                    overlap_interv.append([self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end,1, [self.c_end+1, min(ins_end_coord,b.c_end),'insertion-tandem_multiple_copy',min(ins_end_coord,b.c_end)-self.c_end-1+1,interv_case+'.1']])

                    overlap_st=FIND_CONT_COORD_FORWARD_END_1(self.r_end, self.c_end, b.r_st, self.errors,self.c_st)
                    overlap_interv.append([self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end,1, [self.c_end+1, min(ins_end_coord,b.c_end),'insertion-tandem_multiple_copy',min(ins_end_coord,b.c_end)-self.c_end-1+1,interv_case+'.1']])

                
                    errors_list.append([self.c_end+1, min(ins_end_coord,b.c_end),'insertion-tandem_multiple_copy',min(ins_end_coord,b.c_end)-self.c_end-1+1,interv_case+'.1',overlap_interv])
                    interv_case='9.4.1'
                    b.errors=b_errors_ins

                    

                   
                

            else:
                ins_end_coord=ins_end_coord_2
                last_err_end=last_err_end_2
                b.errors=b_errors_del

                
                
                first_base=last_err_end+1
                if min(ins_end_coord,b.r_end)-self.r_end-1+1>0:

                        overlap_st=FIND_REF_COORD_FORWARD_END_1(self.r_end, self.c_end, b.c_st, self.errors, self.r_st)
                        overlap_interv.append([self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end, 1,[self.c_end, self.c_end,'deletion-collapsed_tandem_repeat',min(ins_end_coord,b.r_end)-self.r_end-1+1,interv_case+'.2']])

                        overlap_st=FIND_CONT_COORD_FORWARD_END_1(self.r_end, self.c_end, b.r_st, self.errors,self.c_st)
                        overlap_interv.append([self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end,1, [self.c_end, self.c_end,'deletion-collapsed_tandem_repeat',min(ins_end_coord,b.r_end)-self.r_end-1+1,interv_case+'.2']])
                    
                        errors_list.append([self.c_end, self.c_end,'deletion-collapsed_tandem_repeat',min(ins_end_coord,b.r_end)-self.r_end-1+1,interv_case+'.2',overlap_interv])

                       
                        
                
                interv_case='9.4.2'
                       
                

            temp_errors=[]
            for entry in self.errors: 
                    temp_errors.append(entry)
            for entry in errors_list:
                    temp_errors.append(entry)
            for entry in b.errors:
                    if entry[0]>=first_base:
                        temp_errors.append(entry)

            
       
            new_interv=[self.c_st,b.c_end,self.r_st,b.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            if interv_case=='9.4.1' or interv_case=='9.4.2' and errors_list!=[]:
                print errors_list
                overlap_st=FIND_REF_COORD_FORWARD_END_1(self.r_end, self.c_end, b.c_st, self.errors,self.r_st)
                overlap_interv.append([self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end,1, errors_list[-1]])

                print cont_name
                self.PRINT_FUNC()
                print
                b.PRINT_FUNC()
                print
                print errors_list[-1]
                print
                print overlap_st
                print
                print self.errors
                print
                print self.r_end, self.c_end, b.c_st
                print self.r_name, overlap_st,self.r_end, cont_name, b.c_st, self.c_end

                
                print
                print
                print ref_seq[overlap_st-1: self.r_end]
                print
                print contig_seq[b.c_st-1: self.c_end]
                   
                

                overlap_st=FIND_CONT_COORD_FORWARD_END_1(self.r_end, self.c_end, b.r_st, self.errors, self.c_st)
                overlap_interv.append([self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end,1, errors_list[-1]])

                
                print cont_name
                self.PRINT_FUNC()
                print
                b.PRINT_FUNC()
                print
                print errors_list[-1]
                print
                print overlap_st
                print
                print self.errors
                print
                print self.r_end, self.c_end, b.r_st
                print self.r_name, b.r_st,self.r_end, cont_name, overlap_st, self.c_end
                

                
                print
                print
                print ref_seq[b.r_st-1: self.r_end]
                print
                print contig_seq[overlap_st-1: self.c_end]
                raw_input('jdg9.4')    
                
            '''
    

            
        elif interv_case=='12.8':
           
            errors_list=[]

            old_list=[]
            for entry in b.errors:
                old_list.append(entry)


            b_errors_del=[]
            b_errors_ins=[]
            for entry in b.errors:
            
                b_errors_del.append([entry[0],entry[1],entry[2],entry[3],entry[4]])
                b_errors_ins.append([entry[0],entry[1],entry[2],entry[3],entry[4]])
                

           
            ins_end_coord_2, last_err_end_2=FIND_REF_COORD_REVERSE_END(b.r_end, b.c_st, self.c_end, b_errors_del)
            

            
            ins_end_coord,last_err_end=FIND_CONT_COORD_REVERSE_END(b.r_end, b.c_st, self.r_st, b_errors_ins, b.c_end)

            first_base=last_err_end+1

            if min(ins_end_coord,b.c_end)-self.c_end-1+1>0:
                overlap_end=FIND_REF_COORD_REVERSE_END_1(self.r_st, self.c_end, b.c_st, self.errors,self.r_end)
                overlap_interv.append([self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end,-1, [self.c_end+1, min(ins_end_coord,b.c_end),'insertion-tandem_multiple_copy',min(ins_end_coord,b.c_end)-self.c_end-1+1,interv_case+'.1']])

                overlap_st=FIND_CONT_COORD_REVERSE_END_1(self.r_st, self.c_end, b.r_end, self.errors, self.c_st)
                overlap_interv.append([self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end,-1,[self.c_end+1, min(ins_end_coord,b.c_end),'insertion-tandem_multiple_copy',min(ins_end_coord,b.c_end)-self.c_end-1+1,interv_case+'.1']])


                errors_list.append([self.c_end+1, min(ins_end_coord,b.c_end),'insertion-tandem_multiple_copy',min(ins_end_coord,b.c_end)-self.c_end-1+1,interv_case+'.1',overlap_interv])
                interv_case='12.8.1'
                b.errors=b_errors_ins

                

            else:
                ins_end_coord=ins_end_coord_2
                last_err_end=last_err_end_2
                b.errors=b_errors_del
                
                first_base=last_err_end+1
                if self.r_st-1-ins_end_coord+1>0:
                        overlap_end=FIND_REF_COORD_REVERSE_END_1(self.r_st, self.c_end, b.c_st, self.errors,self.r_end)
                        overlap_interv.append([self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end,-1, [self.c_end, self.c_end,'deletion-collapsed_tandem_repeat',self.r_st-1-ins_end_coord+1,interv_case+'.2']])

                        overlap_st=FIND_CONT_COORD_REVERSE_END_1(self.r_st, self.c_end, b.r_end, self.errors,self.c_st)
                        overlap_interv.append([self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end,-1,[self.c_end, self.c_end,'deletion-collapsed_tandem_repeat',self.r_st-1-ins_end_coord+1,interv_case+'.2']])

                        
                        errors_list.append([self.c_end, self.c_end,'deletion-collapsed_tandem_repeat',self.r_st-1-ins_end_coord+1,interv_case+'.2', overlap_interv])

                        
                interv_case='12.8.2'
                       

             
                

            temp_errors=[]
            for entry in self.errors: 
                    temp_errors.append(entry)
            for entry in errors_list:
                    temp_errors.append(entry)
            for entry in b.errors:
                    if entry[0]>=first_base:
                        temp_errors.append(entry)


             
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            '''
            if interv_case=='12.8.1' or interv_case=='12.8.2' and errors_list!=[]:
                print errors_list
                overlap_end=FIND_REF_COORD_REVERSE_END_1(self.r_st, self.c_end, b.c_st, self.errors,self.r_end)
                overlap_interv.append([self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end,-1, errors_list[-1]])

                print cont_name
                self.PRINT_FUNC()
                print
                b.PRINT_FUNC()
                print
                print errors_list[-1]
                print
                print overlap_st
                print
                print self.errors
                print
                print self.r_st, self.c_end, b.c_st
                print self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end
                

                
                print
                print
                import general
                print general.COMPL_STRING(ref_seq[self.r_st-1:overlap_end])
                print
                print contig_seq[b.c_st-1: self.c_end]
                   
                

                overlap_st=FIND_CONT_COORD_REVERSE_END_1(self.r_st, self.c_end, b.r_end, self.errors, self.c_st)
                overlap_interv.append([self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end,-1,errors_list[-1]])

                
                print cont_name
                self.PRINT_FUNC()
                print
                b.PRINT_FUNC()
                print
                print errors_list[-1]
                print
                print overlap_st
                print
                print self.errors
                print
                print self.r_st, self.c_end, b.r_end
                print self.r_name, self.r_st,b.r_end, cont_name, overlap_st, self.c_end
                

                
                print
                print
                import general
                print general.COMPL_STRING(ref_seq[self.r_st-1:b.r_end])
                print
                print contig_seq[overlap_st-1: self.c_end]
                raw_input('jdg12.8')    
                
            '''
    

        
 
        elif interv_case=='12.9':
            c_space_len=self.c_end-b.c_st+1

            errors_list=[]

           

           
            corresp_ref_coord, last_err_end=FIND_REF_COORD_REVERSE_END(b.r_end, b.c_st, self.c_end, b.errors)
            if b.r_end-corresp_ref_coord+1!=0:
                overlap_end=FIND_REF_COORD_REVERSE_END_1(self.r_st, self.c_end, b.c_st, self.errors,self.r_end)
                overlap_interv.append([self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end,-1, [self.c_end, self.c_end, 'deletion-collapsed_tandem_repeat',b.r_end-corresp_ref_coord+1,interv_case]])

                errors_list.append([self.c_end, self.c_end, 'deletion-collapsed_tandem_repeat',b.r_end-corresp_ref_coord+1,interv_case,overlap_interv])

                
            first_base=last_err_end+1

                
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base:
                    temp_errors.append(entry)
            
                
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print overlap_end
            print
            print self.errors
            print
            print self.r_st, self.c_end, b.c_st
            print self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end
            

            
            print
            print
            import general
            print general.COMPL_STRING(ref_seq[self.r_st-1:overlap_end])
            print
            print contig_seq[b.c_st-1: self.c_end]
            print general.COMPL_STRING(ref_seq[corresp_ref_coord-1:b.r_end])
            raw_input('jdg12.9')    
            '''
    

        elif interv_case=='12.10':
            c_space_len=self.c_end-b.c_st+1
            ref_space_len=self.r_st-b.r_end-1

            errors_list=[[self.c_end, self.c_end, 'deletion', ref_space_len,interv_case,[]]]

           
            
           
            
            corresp_ref_coord, last_err_end=FIND_REF_COORD_REVERSE_END(b.r_end, b.c_st, self.c_end, b.errors)
            
            if b.r_end-corresp_ref_coord+1!=0:
                overlap_end=FIND_REF_COORD_REVERSE_END_1(self.r_st, self.c_end, b.c_st, self.errors, self.r_end)
                overlap_interv.append([self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end,-1, [self.c_end, self.c_end, 'deletion-collapsed_repeat',b.r_end-corresp_ref_coord+1,interv_case]])

                errors_list.append([self.c_end, self.c_end, 'deletion-collapsed_repeat',b.r_end-corresp_ref_coord+1,interv_case,overlap_interv])

                
            first_base=last_err_end+1
            
            temp_errors=[]
            for entry in self.errors: 
                temp_errors.append(entry)
            for entry in errors_list:
                temp_errors.append(entry)
            for entry in b.errors:
                if entry[0]>=first_base: 
                    temp_errors.append(entry)
            

            
                       
            new_interv=[self.c_st,b.c_end,b.r_st,self.r_end,self.c_dir,self.r_dir, self.c_len, self.r_len, self.r_name,b.num,temp_errors]

            
            '''
            print cont_name
            self.PRINT_FUNC()
            print
            b.PRINT_FUNC()
            print
            print errors_list[-1]
            print
            print overlap_end
            print
            print self.errors
            print
            print self.r_st, self.c_end, b.c_st
            print self.r_name, self.r_st,overlap_end, cont_name, b.c_st, self.c_end
            

            
            print
            print
            import general
            print general.COMPL_STRING(ref_seq[self.r_st-1:overlap_end])
            print
            print contig_seq[b.c_st-1: self.c_end]
            print general.COMPL_STRING(ref_seq[corresp_ref_coord-1:b.r_end])
            raw_input('jdg12.10')    
            '''
    
    
        else:
            errors_list=[[]]
            new_interv=[]

     
        return new_interv, interv_case







#----------------------------------------------------------------------------------
def ANALYSE_SPACE_SIMB(sequence, start, end):

    N_num=0
    
    
    for i in range(start, end+1):
        if sequence[i-1]=='N' or sequence[i-1]=='n':
            N_num+=1

    space_len=end-start+1

    if N_num==space_len:
        space_type='gap'
    elif N_num==0:
        space_type='nucleotides'
    else:
        space_type='mixed_fragment'

    
        
            
    return space_type


def FIND_INSERTION_GAP_INTERVALS(sequence, st, end,interv_case):

    insertion_gap_intervals=[]

   
    if not sequence[st-1]in 'Nn':
        current_interval=[st,st,'insertion',1,interv_case]
    else:
        current_interval=[st,st,'wrong_gap',1,interv_case]

    
    for i in range(st+1, end+1):
        if sequence[i-1] not in 'Nn':
            type_simb='insertion'
        else:
            type_simb='wrong_gap'


        if current_interval[2]==type_simb:
            current_interval[1]+=1
            current_interval[3]+=1

        else:
            insertion_gap_intervals.append(current_interval)
            current_interval=[i,i,type_simb,1,interv_case]
    
    insertion_gap_intervals.append(current_interval)

    return insertion_gap_intervals



def INSERTION_INSIDE_CONTIG(start, end, seq, space_len, interv_case):

    space_type=ANALYSE_SPACE_SIMB(seq, start, end)

    if space_type=='gap':

        errors_list=[[start,end,'wrong_gap',space_len,interv_case]]
        
    elif space_type=='nucleotides':

        errors_list=[[start,end,'insertion',space_len,interv_case]]
        
    else: # space_type='mixed_fragment'
        errors_list=FIND_INSERTION_GAP_INTERVALS(seq, start, end,interv_case)


    return errors_list


def INSERTION_CONTIG_REFERENCE(start, end,seq, c_space_len, ref_space_len,interv_case):
    space_type=ANALYSE_SPACE_SIMB(seq, start, end)

    if space_type=='gap':
        if c_space_len<ref_space_len:
            errors_list=[[start,end,'gap',c_space_len, interv_case]]
            errors_list.append([end,end,'deletion-gap_underestimated',ref_space_len-c_space_len,interv_case])
            
        elif c_space_len>ref_space_len:
            errors_list=[[start,start+ref_space_len-1,'gap',ref_space_len,interv_case]]
            errors_list.append([start+ref_space_len, end, 'wrong_gap-overestimated',c_space_len-ref_space_len ,interv_case])

        else:
            errors_list=[[start, end, 'gap', c_space_len,interv_case]]    
        
    elif space_type=='nucleotides':
        if c_space_len<ref_space_len:
            errors_list=[[start,end,'substitution',c_space_len, interv_case]]
            errors_list.append([end,end,'deletion',ref_space_len-c_space_len,interv_case])
            

        elif c_space_len>ref_space_len:
            errors_list=[[start,start+ref_space_len-1,'substitution',ref_space_len,interv_case]]
            errors_list.append([start+ref_space_len, end, 'insertion',c_space_len-ref_space_len ,interv_case])

        else:
            errors_list=[[start, end, 'substitution', c_space_len,interv_case]]    
        
    else: # space_type='mixed_fragment'
        

        
        if c_space_len<ref_space_len: 
            errors_list_temp=FIND_INSERTION_GAP_INTERVALS(seq, start, end,interv_case)

            for entry in errors_list_temp:
                if entry[2]=='insertion':
                    entry[2]='substitution'
                elif entry[2]=='wrong_gap':
                    entry[2]='gap'
                else:
                    sys.exit('unknown case in INSERTION_CONTIG_REFERENCE()')

            errors_list=errors_list_temp
            errors_list.append([end, end, 'deletion',ref_space_len-c_space_len,interv_case ])

        elif c_space_len>ref_space_len:
            errors_list_temp=FIND_INSERTION_GAP_INTERVALS(seq, start, start+ref_space_len-1,interv_case)

            for entry in errors_list_temp:
                if entry[2]=='insertion':
                    entry[2]='substitution'
                elif entry[2]=='wrong_gap':
                    entry[2]='gap'
                else:
                    sys.exit('unknown case in INSERTION_CONTIG_REFERENCE()')

            errors_list=errors_list_temp

            errors_list_temp=FIND_INSERTION_GAP_INTERVALS(seq, start+ref_space_len, end,interv_case)

            for entry in errors_list_temp:
                errors_list.append(entry)

        else: #c_space_len==ref_space_len:
            errors_list_temp=FIND_INSERTION_GAP_INTERVALS(seq, start, start+ref_space_len-1,interv_case)

            for entry in errors_list_temp:
                if entry[2]=='insertion':
                    entry[2]='substitution'
                elif entry[2]=='wrong_gap':
                    entry[2]='gap'
                else:
                    sys.exit('unknown case in INSERTION_CONTIG_REFERENCE()')

            errors_list=errors_list_temp
                                
        

    return errors_list

def FIND_CONT_COORD_FORWARD_START(r_st, c_st, ref_pos, errors_list,c_end):
    c_pos=c_st+(ref_pos-r_st+1)-1
    used_i_list=[]

    save_list=[]

    flag_del=0
    for i in range(len(errors_list)):
        if errors_list[i][0]<=c_pos:

            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    c_pos+=errors_list[i][3]
                    used_i_list.append(i)
                else:
                    c_pos=errors_list[i][1]+1+c_pos-errors_list[i][0]
                    used_i_list.append(i)
                    
                    

            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                 if errors_list[i][0]<c_pos and errors_list[i][1]<=c_pos:
                     used_i_list.append(i)
                 elif errors_list[i][0]==c_pos and errors_list[i][1]==c_pos:
                     used_i_list.append(i)
                     
                    
                 else:
                     errors_list[i][0]=c_pos+1
                     errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    c_pos-=errors_list[i][3]

                    if c_pos >errors_list[i][0]:
                        used_i_list.append(i)
                    elif c_pos==errors_list[i][0]:
                        errors_list[i][3]=1
                        used_i_list.append(i)
                    else:
                        
                        del_len=errors_list[i][0]-c_pos
                        if c_end >=errors_list[i][0]+1:
                            c_pos=errors_list[i][0]
                            
                            if del_len==0:
                               used_i_list.append(i)
                            else:
                                errors_list[i][3]=del_len
                                save_list.append(i)
                                
                        else:
                            
                            c_pos=c_end
                            errors_list[i][3]=del_len-(c_end-errors_list[i][0]+1)
                            save_list.append(i)
                        
                        flag_del=1
                        
                        
                else:
                    flag_del=1
                    save_list.append(i)
            

    

    if used_i_list==[]:
        
        if flag_del==1:
            last_err_end=-3
        else:
            last_err_end=-3
    else:
        max_i=max(used_i_list)
        last_err_end=errors_list[max_i][1]

        if flag_del==1:
            list_temp=[]
            for i in range(len( errors_list)):
                if errors_list[i][0]<last_err_end:
                    list_temp.append(errors_list[i])
                elif errors_list[i][0]==last_err_end:
                    if errors_list[i][2].startswith('deletion'):
                        if i in save_list:
                            list_temp.append(errors_list[i])
                else:
                   list_temp.append(errors_list[i])

            for i in range (len(errors_list)):
                errors_list.pop(0)

            for entry in list_temp:
                errors_list.append(entry)

            last_err_end-=1

            

       
    return c_pos, last_err_end
                    
def FIND_CONT_COORD_REVERSE_END_1(r_st, c_end, ref_pos, errors_list,c_st):
    c_pos=c_end-(ref_pos-r_st+1)+1

    
    flag_del=0
    for i in range(len(errors_list)-1,-1,-1):
        if errors_list[i][1]>=c_pos:
            
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                c_pos-=errors_list[i][3]
                    
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                 a='do_nothing'    
                 
            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]>=c_pos and errors_list[i][1]>=c_pos:
                    c_pos+=errors_list[i][3]


    if c_pos>c_end:
        c_pos=c_end
    elif c_pos<c_st:
        c_pos=c_st
   
    return c_pos
def FIND_CONT_COORD_FORWARD_END_1(r_end, c_end, ref_pos, errors_list,c_st):

    c_pos=c_end-(r_end-ref_pos+1)+1

    
    flag_del=0
    for i in range(len(errors_list)-1,-1,-1):
        if errors_list[i][1]>=c_pos:
            
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                c_pos-=errors_list[i][3]
                    
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                 a='do_nothing'    
                 
            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]>=c_pos and errors_list[i][1]>=c_pos:
                    c_pos+=errors_list[i][3]

                
    if c_pos>c_end:
        c_pos=c_end
    elif c_pos<c_st:
        c_pos=c_st
   
    
    return c_pos


def FIND_CONT_COORD_REVERSE_END_SECOND(r_st, c_end, r_coord, errors_list):
    c_pos=c_end-(r_coord-r_st+1)+1

    used_i_list=[]
    save_list=[]

    
    flag_del=0
    for i in range(len(errors_list)-1,-1,-1):
        if errors_list[i][1]>=c_pos:

            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]>c_pos and errors_list[i][1]>c_pos:
                    c_pos-=errors_list[i][3]
                    used_i_list.append(i)
                else:
                    c_pos-=errors_list[i][3]
                    used_i_list.append(i)
                    
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                 if errors_list[i][0]>=c_pos and errors_list[i][1]>c_pos:
                     used_i_list.append(i)
                     
                 elif errors_list[i][0]==c_pos and errors_list[i][1]==c_pos:
                     used_i_list.append(i)
                     
                 else:
                     errors_list[i][1]=c_pos-1
                     errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]>=c_pos and errors_list[i][1]>=c_pos:
                    c_pos+=errors_list[i][3]

                    if c_pos <errors_list[i][0]:
                        used_i_list.append(i)
                    elif c_pos==errors_list[i][0]:
                        used_i_list.append(i)
                    else:
                        del_len=c_pos-errors_list[i][0]-1
                        if 1<=errors_list[i][0]-1:
                            c_pos=errors_list[i][0]+1
                            if del_len==0:
                                used_i_list.append(i)
                            else:
                                errors_list[i][3]=del_len
                                save_list.append(i)
                            
                        else:
                            
                            c_pos=1
                            errors_list[i][3]=del_len-(errors_list[i][0]-1+1)
                            save_list.append(i)
                        
                        flag_del=1
                        
                        
                
            

    

    if used_i_list==[]:
        
        if flag_del==1:
            last_err_end=c_end+1
        else:
            last_err_end=c_end+1
    else:
        min_i=min(used_i_list)
        last_err_end=errors_list[min_i][0]

        if flag_del==1:
            list_temp=[]
            for i in range(len( errors_list)):
                if errors_list[i][0]<last_err_end:
                    list_temp.append(errors_list[i])
                elif errors_list[i][0]==last_err_end:
                    if errors_list[i][2].startswith('deletion'):
                        
                        if i in save_list:
                            list_temp.append(errors_list[i])
                else:
                   list_temp.append(errors_list[i])

            for i in range (len(errors_list)):
                errors_list.pop(0)

            for entry in list_temp:
                errors_list.append(entry)

            last_err_end+=1

            

    
    
    
       
    return c_pos, last_err_end
                    
 

def FIND_CONT_COORD_REVERSE_END(r_end, c_st, ref_pos, errors_list, c_end):

    c_pos=c_st+(r_end-ref_pos+1)-1

    used_i_list=[]
    save_list=[]

    
    flag_del=0
    for i in range(len(errors_list)):
        if errors_list[i][0]<=c_pos:

            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    c_pos+=errors_list[i][3]
                    used_i_list.append(i)
                else:
                    c_pos=errors_list[i][1]+1+c_pos-errors_list[i][0]
                    used_i_list.append(i)

            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                 if errors_list[i][0]<c_pos and errors_list[i][1]<=c_pos:
                     used_i_list.append(i)
                 elif errors_list[i][0]==c_pos and errors_list[i][1]==c_pos:
                     used_i_list.append(i)
                     
                    
                 else:
                     errors_list[i][0]=c_pos+1
                     errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    c_pos-=errors_list[i][3]

                    if c_pos >errors_list[i][0]:
                        used_i_list.append(i)
                    elif c_pos==errors_list[i][0]:
                        errors_list[i][3]=1
                        used_i_list.append(i)
                    else:
                        
                        del_len=errors_list[i][0]-c_pos
                        if c_end >=errors_list[i][0]+1:
                            c_pos=errors_list[i][0]
                            if del_len==0:
                                used_i_list.append(i)
                            else:
                                errors_list[i][3]=del_len
                                save_list.append(i)
                        else:
                            
                            c_pos=c_end
                            errors_list[i][3]=del_len-(c_end-errors_list[i][0]+1)
                            save_list.append(i)
                        
                        flag_del=1
                        
                        
                else:
                    flag_del=1
                    save_list.append(i)
            

    

    if used_i_list==[]:
        
        if flag_del==1:
            last_err_end=-3
        else:
            last_err_end=-3
    else:
        max_i=max(used_i_list)
        last_err_end=errors_list[max_i][1]

        if flag_del==1:
            list_temp=[]
            for i in range(len( errors_list)):
                if errors_list[i][0]<last_err_end:
                    list_temp.append(errors_list[i])
                elif errors_list[i][0]==last_err_end:
                    if errors_list[i][2].startswith('deletion'):
                        if i in save_list:
                            list_temp.append(errors_list[i])
                else:
                   list_temp.append(errors_list[i])

            for i in range (len(errors_list)):
                errors_list.pop(0)

            for entry in list_temp:
                errors_list.append(entry)

            last_err_end-=1

            

       
    return c_pos, last_err_end
                    
 
##use in interv
def FIND_REF_COORD_FORWARD_START(r_st, c_st, c_pos, errors_list):
    
    r_pos=r_st+(c_pos-c_st+1)-1
    
   
    used_i_list=[]
    save_list=[]

    
    
    for i in range(len(errors_list)):
        if errors_list[i][0]<=c_pos:
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    r_pos-=errors_list[i][3]
                    used_i_list.append(i)                    

                elif (errors_list[i][0]<c_pos and errors_list[i][1]==c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]==c_pos): 
                
                    r_pos-=errors_list[i][3]
                    used_i_list.append(i)

                

                elif (errors_list[i][0]<c_pos and errors_list[i][1]>c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]>c_pos):
                    r_pos-=c_pos-errors_list[i][0]+1
                    errors_list[i][0]=c_pos+1
                    errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

                    if r_pos<r_st:
                        r_pos=r_st-1

                        

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    r_pos+=errors_list[i][3]
                    used_i_list.append(i)
                else:
                    save_list.append(i)

            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    used_i_list.append(i)                    

                elif (errors_list[i][0]<c_pos and errors_list[i][1]==c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]==c_pos):
                    used_i_list.append(i)

                elif (errors_list[i][0]<c_pos and errors_list[i][1]>c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]>c_pos):
                    errors_list[i][0]=c_pos+1
                    errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1
           

    if used_i_list==[]:
        last_err_end=-3
    else:
        max_i=max(used_i_list)
        last_err_end=errors_list[max_i][1]

        if save_list!=[]:
            list_temp=[]
            for i in range(len( errors_list)):
                if errors_list[i][0]<last_err_end:
                    list_temp.append(errors_list[i])
                elif errors_list[i][0]==last_err_end:
                    if errors_list[i][2].startswith('deletion'):
                        
                        if i in save_list:
                            list_temp.append(errors_list[i])
                else:
                   list_temp.append(errors_list[i])

            for i in range (len(errors_list)):
                errors_list.pop(0)

            for entry in list_temp:
                errors_list.append(entry)

            last_err_end-=1
    
    

    return r_pos, last_err_end

def FIND_REF_COORD_REVERSE_END_1(r_st, c_end, c_pos, errors_list,r_end):

    

    r_pos=r_st+(c_end-c_pos+1)-1

    

    
    for i in range(len(errors_list)-1,-1,-1):
        if errors_list[i][1]>=c_pos:
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]>c_pos and errors_list[i][1]>c_pos:
                    r_pos-=errors_list[i][3]
                   
                elif (errors_list[i][0]<c_pos and errors_list[i][1]==c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]==c_pos): 
                
                    r_pos-=errors_list[i][1]-c_pos+1
                   
                elif (errors_list[i][0]<c_pos and errors_list[i][1]>c_pos):
                    r_pos-=errors_list[i][1]-c_pos+1

                elif (errors_list[i][0]==c_pos and errors_list[i][1]>c_pos):
                    r_pos-=errors_list[i][3]
                   
                        

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]>=c_pos and errors_list[i][1]>=c_pos:
                    r_pos+=errors_list[i][3]
                    
                   
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                  a='do_nothing'
    if r_pos>r_end:
        r_pos=r_end
    elif r_pos<r_st:
        r_pos=r_st
   
    

    return r_pos


            
def FIND_REF_COORD_FORWARD_END_1(r_end, c_end, c_pos, errors_list,r_st):

    

    r_pos=r_end-(c_end-c_pos+1)+1

    

    
    for i in range(len(errors_list)-1,-1,-1):
        if errors_list[i][1]>=c_pos:
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]>c_pos and errors_list[i][1]>c_pos:
                    r_pos+=errors_list[i][3]
                   
                elif (errors_list[i][0]<c_pos and errors_list[i][1]==c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]==c_pos): 
                
                    r_pos+=errors_list[i][1]-c_pos+1
                   
                elif (errors_list[i][0]<c_pos and errors_list[i][1]>c_pos):
                    r_pos+=errors_list[i][1]-c_pos+1

                elif (errors_list[i][0]==c_pos and errors_list[i][1]>c_pos):
                    r_pos+=errors_list[i][3]
                   
                        

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]>=c_pos and errors_list[i][1]>=c_pos:
                    r_pos-=errors_list[i][3]
                    
                   
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                  a='do_nothing'

    if r_pos>r_end:
        r_pos=r_end
    elif r_pos<r_st:
        r_pos=r_st       

    return r_pos


##use in nuc
def FIND_REF_COORD_FORWARD_END(r_end, c_end, c_pos, errors_list):

     

    r_pos=r_end-(c_end-c_pos+1)+1

    

    used_i_list=[] 
    for i in range(len(errors_list)):
        
        if errors_list[i][0]>c_pos:
            
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                r_pos+=errors_list[i][3]
            elif errors_list[i][2].startswith('deletion'):
                r_pos-=errors_list[i][3]
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                a='do nothing'
            elif errors_list[i][2]=='wrong_end' or errors_list[i][2]=='wrong_beginning':
                r_pos+=errors_list[i][3]

            used_i_list.append(i)

        elif errors_list[i][0]==c_pos:
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                
                r_pos+=1
                errors_list[i][3]-=1
                errors_list[i][0]+=1
            elif errors_list[i][2].startswith('deletion'):
                
                r_pos-=errors_list[i][3]
            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                a='do nothing'
            elif errors_list[i][2]=='wrong_end' or errors_list[i][2]=='wrong_beginning':
                r_pos+=errors_list[i][3]

            used_i_list.append(i)
        else:
            break

    
    if used_i_list==[]:
        last_err_end=c_pos
    else:
        max_i=min(used_i_list)
        last_err_end=errors_list[max_i][1]
    

    return r_pos, last_err_end

#use in interv
def FIND_REF_COORD_REVERSE_END(r_end, c_st, c_pos, errors_list):

    
    r_pos=r_end-(c_pos-c_st+1)+1

    used_i_list=[]

    used_i_list=[]
    save_list=[]
    add_errors=[]

    for i in range(len(errors_list)):
        if errors_list[i][0]<=c_pos:
            if errors_list[i][2].startswith('insertion') or errors_list[i][2].startswith('wrong_gap'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    r_pos+=errors_list[i][3]
                    used_i_list.append(i)                    

                elif (errors_list[i][0]<c_pos and errors_list[i][1]==c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]==c_pos):
                    r_pos+=errors_list[i][3]
                    used_i_list.append(i)

                elif (errors_list[i][0]<c_pos and errors_list[i][1]>c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]>c_pos):
                    r_pos+=c_pos-errors_list[i][0]+1
                    
                    if r_pos>r_end:
                        
                        errors_list[i][0]=c_pos+1
                        errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

                        r_pos=r_end+1

                        if errors_list[i][3]==0:
                            used_i_list.append(i)
                           
                        
                    else:
                        errors_list[i][0]=c_pos+1
                        errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

                        

            elif errors_list[i][2].startswith('deletion'):
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    r_pos-=errors_list[i][3]
                    used_i_list.append(i)
                else:
                    save_list.append(i)

            elif errors_list[i][2]=='substitution' or errors_list[i][2]=='gap' :
                if errors_list[i][0]<c_pos and errors_list[i][1]<c_pos:
                    used_i_list.append(i)                    

                elif (errors_list[i][0]<c_pos and errors_list[i][1]==c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]==c_pos):
                    used_i_list.append(i)

                elif (errors_list[i][0]<c_pos and errors_list[i][1]>c_pos) or (errors_list[i][0]==c_pos and errors_list[i][1]>c_pos):
                    errors_list[i][0]=c_pos+1
                    errors_list[i][3]=errors_list[i][1]-errors_list[i][0]+1

            

    if used_i_list==[]:
        last_err_end=-3
    else:
        max_i=max(used_i_list)
        last_err_end=errors_list[max_i][1]

        if save_list!=[]:
            list_temp=[]
            for i in range(len( errors_list)):
                if errors_list[i][0]<last_err_end:
                    list_temp.append(errors_list[i])
                elif errors_list[i][0]==last_err_end:
                    if errors_list[i][2].startswith('deletion'):
                        if i in save_list:
                            list_temp.append(errors_list[i])
                else:
                   list_temp.append(errors_list[i])

            for i in range (len(errors_list)):
                errors_list.pop(0)

            for entry in list_temp:
                errors_list.append(entry)

            last_err_end-=1

    
    return r_pos, last_err_end

