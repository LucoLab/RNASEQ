
'''

:date: Sep 6, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Align the samples based on config set up in json format file.

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import datetime
from pandas.io.excel import ExcelWriter
import sys
import os
import pandas as pd
from pathlib import Path


###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################
def write_subprocess_log(completedProcess,logger):
    """
    Write in log the stdout or stderr of subprocess.
    Tcheck if everything was ok.
  
    Args:
        completedProcess (obj): Instance of CompletedProcess send by subprocess.run().
        logger (obj): Instance of logging().
  
    """
    try :
        completedProcess.check_returncode()
        logger.info(completedProcess.stdout)
    except subprocess.CalledProcessError as exc:
                logger.error("===> Exception Caught : ")
                logger.error(exc)  
                logger.error("====> Standard Error : ")
                logger.error(completedProcess.stderr) 
                

def create_logger(config,id):
    """
    Define a logger instance to write in a file and stream.
  
    Args:
        config (obj): Configuration instance.

    Returns:
        logger (obj): Logger instance to log messages in file and output mainstream.
    
        
    """
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+"/"+id+"."+'supperWrapper_activity.log', 'a', 1000000, 1)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    return logger


###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
  
    
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to config file.",required=True,type=str,dest='file_config')
    parser.add_argument("-ct","--control",action="store_true",help="AS or Control",default=False,dest='iscontrol')
    parser.add_argument("-i","--init",action="store",help="Path to a json file.",required=False,default=str(Path(os.path.dirname(os.path.realpath(__file__))).parents[0])+"/config/init.json",type=str,dest='init')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")
    init = custom_parser.Configuration(parameters.init,"json")

    logger = create_logger(config,config.parameters['project'])

    iscontrol      = parameters.iscontrol
    read_cutoff    = config.parameters['reads_junction']
    path_to_output =   config.parameters['path_to_output'] 
    
    applicatifDir   = init.parameters['scriptDir']
    
    print(os.path.basename(__file__))
    
    for event in ["SE","A3SS","A5SS","RI"] : #,"SE" "MXE","SE","A3SS","A5SS","RI""SE"
        
        for analyse in config.parameters['tabs'] : 
            
            subprocess.run(("mkdir -p "+config.parameters['path_to_output']+"/"+analyse+"/"+event),shell=True)

            logger.info(analyse)
            bashcommand =applicatifDir+"bash/merge.and.clean.splicing.sh "+event+" "+analyse+" "+read_cutoff+" "+parameters.file_config+" "+path_to_output+" "+str(iscontrol)+" "+applicatifDir+" "+init.parameters['python3']
            
            logger.info(bashcommand)
            
            bash = subprocess.run((bashcommand),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            logger.info(bash.stdout)
            logger.info(bash.stderr)

            write_subprocess_log(bash,logger)

    