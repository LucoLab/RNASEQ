'''
Created on Sep 6, 2016

:platform: Unix
:synopsis: Load files via Json
:author: Villemin Jean-Philippe
'''

import json
import datetime

class Configuration(object):
    """Configuration Class.
    
        Configuration is a class that let us load parameters from a file.
    
        Args:
            path_to_dir (str): Pathway to the config file.
            type_of_file (str): Type of data to parse.
        Returns:
            instance of Configuration (obj) : A config with parameters and time associated.
        Note:
            I should do a Class with inheritance for specific files
        Todo:
            * Class with inheritance
            
    """
 
    def __init__(self,path_to_dir,type_of_file):
       
        self.path_to_dir = path_to_dir
        timeline = datetime.datetime.now()
        self.chrono = str(timeline.day)+"_"+str(timeline.month)+"_"+str(timeline.year)+"__"+str(timeline.hour)+"_"+str(timeline.minute)+"_"+str(timeline.second)

        if type_of_file =="json": 
            with open(path_to_dir,'r') as stream:
                try:
                    self.parameters = json.load(stream)
                    stream.close()
                
                except json.JSONDecodeError as exc:
                    print("JSON is corrupted :")    
                    print(exc)                        
