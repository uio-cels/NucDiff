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

from . import class_nucmer

def FIND_ERRORS_ASSEMBLY(file_ref,file_contigs, working_dir, nucmer_opt, prefix, proc_num,filter_opt,delta_file,reloc_dist,coord_file ):
    mummer_prefix=working_dir+prefix

    
    #1. find fragments for each query sequence
    print()
    print('Run NUCmer...')
    print()
    nuc=class_nucmer.Nucmer(mummer_prefix, file_ref, file_contigs, working_dir,delta_file,coord_file)    
    nuc.RUN(nucmer_opt,filter_opt)
    nuc.PARSE()

    print()
    print('Find differences...')
    print()

    #2. find differences inside fragments
    nuc.FIND_ERR_INSIDE_FRAG(proc_num,file_contigs)
    
   
    #3. find differences  between fragments
    struct_dict,end_err_dict,unmapped_list,uncovered_dict=nuc.FIND_ERRORS(file_contigs, file_ref,reloc_dist)

    print()
    print('The difference detection between fragments step is complete')

    return struct_dict,end_err_dict,unmapped_list,uncovered_dict
    
