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

class Errors:
    def __init__(self):
        self.data={}

    def ADD_LIST(self, errors_list, cont_name):
        if cont_name not in self.data:
            self.data[cont_name]=[]

        for entry in errors_list:
            if entry!=[[]] and entry!=[]  : 
                self.data[cont_name].append(entry)
            


    def MERGE(self, new_errors):
        

        for cont_name in list(new_errors.data.keys()):

            if cont_name not in self.data:
                self.data[cont_name]=[]

            for entry in new_errors.data[cont_name]:
                self.data[cont_name].append(entry)



    def SORT(self):

        for cont_name in list(self.data.keys()):
            
            if len(self.data[cont_name])>1:
                temp_sorted=sorted(self.data[cont_name], key=lambda inter:inter[0], reverse=False)
        
                self.data[cont_name]=temp_sorted

    def PRINT(self):
        for cont_name in list(self.data.keys()):
            print(cont_name)

            for entry in self.data[cont_name]:
                print(entry)
            print()

    def CLASS_TO_LIST(self):
        error_list=[]
        for cont_name in list(self.data.keys()):

            for entry in self.data[cont_name]:    
                error_list.append(entry)
        
        
        return error_list





