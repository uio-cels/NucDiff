# copyright (c) 2014-2016 Ksenia Khelik
#
# This Source Code Form is subject to the 
# terms of the Mozilla Public License, v. 2.0. 
# If a copy of the MPL was not distributed 
# with this file, You can obtain one at 
# https://mozilla.org/MPL/2.0/.
#
#-------------------------------------------------------------------------------



import os
import subprocess
import sys



def CHECK_INPUT_DATA(file_ref, file_contigs, working_dir, prefix, nucmer_opt):
    # check if a file with the reference exists
    cur_dir=os.getcwd()

    
    if not os.path.exists(file_ref):
        print
        print 'ERROR: the provided file with a reference genome does not exist'
        print
        sys.exit(0)

    
    if not os.path.exists(file_contigs):
        print
        print 'ERROR: The provided file with an assembly does not exist'
        print
        sys.exit(0)

    
    if not os.path.exists(working_dir):
        if working_dir.startswith('/'):
            os.makedirs(working_dir)

        else:
            cur_dir=os.getcwd()
            
            os.makedirs(working_dir)
            if working_dir.startswith('./'):
                working_dir=cur_dir+'/'+working_dir.split('./')[1]
            else:
                working_dir=cur_dir+'/'+working_dir

        
       

            
        if not os.path.exists(working_dir):
            print
            print 'ERROR: it is not possible to create working directory'
            print
            sys.exit(0)


    if not working_dir.endswith('/'):
        working_dir+='/'

    if not os.path.exists(working_dir+'/results'):
        os.makedirs(working_dir+'/results')


    #check nucmer run options
    nucmer_opt_list=['--mum', '--mumreference', '--maxmatch', '--breaklen', '-b', '-c', '--mincluster', '--delta', '--nodelta',
                     '--depend', '-d', '--diagfactor', '--extend', '--noextend', '-f', '--forward', '-g','--maxgap','-l', '--minmatch', '-o',
                     '--coords','--optimize','--nooptimize', '-r', '--reverse', '--simplify','--nosimplify']

    temp=nucmer_opt.split(' ')
    for opt in temp:
        if opt.startswith('-'):
            if not opt in nucmer_opt_list:
                if opt=='--prefix' or opt=='-p':
                    print
                    print "ERROR: It is not possible to use {0} option in --nucmer_opt".format(opt)
                    print
                    sys.exit(0)
                else:
                    print
                    print "ERROR: the wrong nucmer option '{0}'".format(opt)
                    print 
                    sys.exit(0)
            

    


    return working_dir

    
    
   

    
    

        
