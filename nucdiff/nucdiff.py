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
import argparse
import os

from .find_errors import FIND_ERRORS_ASSEMBLY
from .generate_output import GENERATE_OUTPUT
from .initial_preparation import CHECK_INPUT_DATA
from .ref_coord import FIND_REF_COORD_ERRORS
from .statistics import FIND_STATISTICS

import subprocess




def START(args):
    
    file_ref=args['file_ref']
    file_contigs=args['file_contig']
    working_dir=args['working_dir']
    prefix=args['prefix']

    reloc_dist=int(args['reloc_dist'])
    
    nucmer_opt=args['nucmer_opt']
    filter_opt=args['filter_opt']
    proc_num=args['proc']

    delta_file=args['delta_file']
    coord_file=''


    asmb_name_full=args['query_name_full']
    ref_name_full=args['ref_name_full']

    vcf_flag=args['vcf']
    
    
          
    #1. check input data correctness
    file_ref, file_contigs, working_dir, delta_file=CHECK_INPUT_DATA(file_ref, file_contigs, working_dir, prefix, nucmer_opt,filter_opt,delta_file)

    

    #2. find differences and gaps in assembly
    struct_dict,end_err_dict,unmapped_list,uncovered_dict=FIND_ERRORS_ASSEMBLY(file_ref,file_contigs, working_dir, nucmer_opt, prefix,proc_num, filter_opt,delta_file,reloc_dist,coord_file ) #class Errors


    #3. find reference-oriented difference coordinates
    err_ref_cont_coord_errors_list,cont_blocks_dict=FIND_REF_COORD_ERRORS(struct_dict,end_err_dict,unmapped_list,file_ref, file_contigs,uncovered_dict)

    #4. find statistics
    statistics_output_lines=FIND_STATISTICS(err_ref_cont_coord_errors_list)


    #5. generate output
    GENERATE_OUTPUT(struct_dict,end_err_dict,unmapped_list, file_ref, file_contigs, working_dir, prefix,err_ref_cont_coord_errors_list,statistics_output_lines,asmb_name_full,ref_name_full,cont_blocks_dict,vcf_flag)


   

    


def main():
    
    argv=sys.argv
    parser = argparse.ArgumentParser()

    parser.add_argument('file_ref',metavar='Reference.fasta', type=str, help='- Fasta file with the reference sequences')
    parser.add_argument('file_contig',metavar='Query.fasta', type=str, help='- Fasta file with the query sequences')
    parser.add_argument('working_dir',metavar='Output_dir', type=str, help='- Path to the directory where all intermediate and final results will be stored')
    parser.add_argument('prefix',metavar='Prefix', type=str, help='- Name that will be added to all generated files including the ones created by NUCmer')
    parser.add_argument('--reloc_dist',metavar='int', type=int,nargs='?',default="10000", help='- Minimum distance between two relocated blocks [10000]')
    parser.add_argument('--nucmer_opt',  type=str, nargs='?',default="", help='- NUCmer run options. By default, NUCmer will be run with its default parameters values, except the --maxmatch parameter. --maxmatch is hard coded and cannot be changed. To change any other parameter values, type parameter names and new values inside single or double quotation marks.')
    parser.add_argument('--filter_opt',  type=str, nargs='?',default="-q", help='- Delta-filter run options. By default, it will be run with -q parameter only. -q is hard coded and cannot be changed. To add any other parameter values, type parameter names and their values inside single or double quotation marks.')
    parser.add_argument('--delta_file',  type=str, nargs='?',default="", help='- Path to the already existing delta file (NUCmer output file) ')
    parser.add_argument('--proc',  metavar='int', type=int, nargs='?',default=1, help='- Number of processes to be used [1] ')
    parser.add_argument('--ref_name_full',  type=str, nargs='?',default='no',choices=['yes','no'], help="- Print full reference names in output files ('yes' value). In case of 'no', everything after the first space will be ignored. ['no']")
    parser.add_argument('--query_name_full',  type=str, nargs='?',default='no',choices=['yes','no'], help="- Print full query names in output files ('yes' value). In case of 'no', everything after the first space will be ignored.['no'] ")
    parser.add_argument('--vcf',  type=str, nargs='?',default='no',choices=['yes','no'], help="- Output small and medium local differences in the VCF format")
    parser.add_argument('--version', action='version', version='NucDiff version 2.0.2')
    
    
    
    
    args=vars(parser.parse_args())

    START(args)

main()
