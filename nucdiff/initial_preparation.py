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

import os
import subprocess
import sys



def CHECK_INPUT_DATA(file_ref, file_contigs, working_dir, prefix, nucmer_opt,filter_opt, delta_file):
    # check if a file with the reference exists
    cur_dir=os.getcwd()

    
    if not os.path.exists(file_ref):
        print()
        print('ERROR: the provided file with a reference genome does not exist')
        print()
        sys.exit(0)
    else:
        file_ref=os.path.abspath(file_ref)
        
    
    if not os.path.exists(file_contigs):
        print()
        print('ERROR: The provided file with an assembly does not exist')
        print()
        sys.exit(0)
    else:
        file_contigs=os.path.abspath(file_contigs)
        

    if delta_file!='':
        if not os.path.exists(delta_file):
            print()
            print('ERROR: the provided delta file does not exist')
            print()
            sys.exit(0)
        else:
            delta_file=os.path.abspath(delta_file)

    working_dir=os.path.abspath(working_dir)

    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

        if not os.path.exists(working_dir):
            print()
            print('ERROR: it is not possible to create working directory')
            print()
            sys.exit(0)


    if not working_dir.endswith('/'):
        working_dir+='/'

    if not os.path.exists(working_dir+'/results'):
        os.makedirs(working_dir+'/results')


    #check nucmer run options
    nucmer_opt_list=['--maxmatch', '--breaklen', '-b', '-c', '--mincluster', '--delta', '--nodelta',
                     '--depend', '-d', '--diagfactor', '--extend', '--noextend', '-f', '--forward', '--maxgap','-l', '--minmatch', '-o',
                     '--coords','--optimize','--nooptimize',  '--reverse', '--simplify','--nosimplify','-g','-r']

    temp=nucmer_opt.split(' ')
    for opt in temp:
        if opt.startswith('-'):
            if not opt in nucmer_opt_list:
                if opt in ['--prefix','-p', '--mum', '--mumreference']:
                    if opt in ['--prefix','-p']:
                        print()
                        print("ERROR: It is not possible to use {0} option in --nucmer_opt".format(opt))
                        print()
                        sys.exit(0)
                    elif opt in ['--mum', '--mumreference']:
                        print()
                        print("ERROR: It is not possible to use {0} option. Use --maxmatch instead.".format(opt))
                        print()
                        sys.exit(0)
                   
                else:
                    print()
                    print("ERROR: the wrong nucmer option '{0}'".format(opt))
                    print() 
                    sys.exit(0)
            
    #check filter run options
    filter_opt_list=['-i', '-l', '-q', '-u','-o']
    



    temp=filter_opt.split(' ')
    for opt in temp:
        if opt.startswith('-'):
            if not opt in filter_opt_list:
                if opt in ['-g','-r']:
                        print()
                        print("ERROR: It is not possible to use {0} option. Use -q instead.".format(opt))
                        print()
                        sys.exit(0)
                   
                else:
                    print()
                    print("ERROR: the wrong delta-filter option '{0}'".format(opt))
                    print() 
                    sys.exit(0)
            
    


    return file_ref, file_contigs, working_dir, delta_file

    
    
   

    
    

        
