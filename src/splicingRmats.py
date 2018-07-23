##!/usr/bin/python3.5
##!/home/jean-philippe.villemin/bin/python-3.5-2/bin/python3

'''

:date: September 9, 2017
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Launch Splicing with Rmats

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import re
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
                
def create_logger(path,LEVEL,id):
    """
    Define a logger instance to write in a file and stream.
  
    Args:
        config (obj): Configuration instance.

    Returns:
        logger (obj): Logger instance to log messages in file and output mainstream.
    
        
    """
    
    logger = logging.getLogger()
    logger.setLevel(LEVEL)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    script = "rmats_activity"
    print(path+"/"+id+"."+script+'.log')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(path+"/"+id+"."+script+'.log', 'a', 1000000, 1)
    file_handler.setLevel(LEVEL)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(LEVEL)
    logger.addHandler(stream_handler)

    return logger


###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will launch rmats on fastq.  
    Example : 
    python ./splicingRmats.py -c /home/jp/workspace/RNA-SEQ/bash/rmats.json -a T5_vs_unT

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")
    
    #ugly trick
    path_above = Path(config.parameters['path_to_output']).parents[2]

    logger = create_logger(str(path_above),"INFO",config.parameters['project'])

    print('\n=========> START:')
    
    print("file_config : "+parameters.file_config)
    print("path_to_output : "+config.parameters['path_to_output'])
    print("path_to_log : "+str(path_above))
    print("gtf : "+config.parameters["gtf"])
    print("softPath : "+config.parameters["softPath"])
    print("index : "+config.parameters["index"])
    print("    ")

    '''
    ######################################################################
    DataFrame Object Construction : First Part
    ######################################################################
    '''
    
    for analysis in config.parameters["analysis"] :
        
        print("Analysis : "+analysis)
        #pathlib.Path(  config.parameters['path_to_output']+analysis).mkdir(parents=True, exist_ok=True) 
        print("len : "+ config.parameters["analysis"][analysis]["type_read"])
        type_read = config.parameters["analysis"][analysis]["type_read"]
        print("len : "+config.parameters["analysis"][analysis]["libType"])
        libType = config.parameters["analysis"][analysis]["libType"]
        print("diff : "+config.parameters["analysis"][analysis]["diff"])
        diff = config.parameters["analysis"][analysis]["diff"]
        print("len : "+config.parameters["analysis"][analysis]["len"])
        leng = config.parameters["analysis"][analysis]["len"]
        #print("test : "+config.parameters["analysis"][analysis]["test"])
        #test = config.parameters["analysis"][analysis]["test"]
        print("Nb samples in Test & Control :  ")
        print( int(len(config.parameters["analysis"][analysis].get("fastq")["s1"])) )
        print( int(len(config.parameters["analysis"][analysis].get("fastq")["s2"])) )
        
        s1=",".join(config.parameters["analysis"][analysis].get("fastq")["s1"])
        s2=",".join(config.parameters["analysis"][analysis].get("fastq")["s2"])
        
        #A: Yes, paired-analysis (-analysis P) option requires at least 3 replicates per sample..
        analysisType = ""
        if( int(len(config.parameters["analysis"][analysis].get("fastq")["s1"]))>=3 and int( len(config.parameters["analysis"][analysis].get("fastq")["s2"]))>=3 ) :
            analysisType = " -analysis P "
        #if( int(len(config.parameters["analysis"][analysis].get("fastq")["s1"]))>=1 and int( len(config.parameters["analysis"][analysis].get("fastq")["s2"]))>=1 ) :

        print( "old version called" +analysisType )

        #splicingRmats = config.parameters["softPath"]+" -s1 "+s1+" -s2 "+s2+" -gtf "+config.parameters["gtf"]+" -t "+type_read+" -libType "+libType+" -len "+leng +" -bi "+ config.parameters["index"]+" -c "+diff+" -o "+config.parameters['path_to_output']+analysis
        splicingRmats = config.parameters["softPath"]+" -s1 "+s1+" -s2 "+s2+analysisType+" -gtf "+config.parameters["gtf"]+" -t "+type_read+" -libType "+libType+" -len "+leng +" -bi "+ config.parameters["index"]+" -c "+diff+" -o "+config.parameters['path_to_output']+analysis

        #else : 
            
           # print( "new version called" )

           # s1_path = config.parameters['path_to_output']+"s1.txt"   
           # s1_file = open(s1_path, 'w')
           # s1_file.write(s1)
           # s1_file.close()
            
           # s2_path = config.parameters['path_to_output']+"s2.txt"   
           # s2_file = open(s2_path, 'w')
           # s2_file.write(s2)
           # s2_file.close()
            
           # splicingRmats = config.parameters["softPath"]+" --s1 "+s1_path+" --s2 "+s2_path+" --nthread 28 --gtf "+config.parameters["gtf"]+" -t "+type_read+" --libType "+libType+" --readLength "+leng +" --bi "+ config.parameters["index"]+" --cstat "+diff+" --tmpPath "+config.parameters['path_to_output']+analysis+"/tmp --od "+config.parameters['path_to_output']+analysis

        logger.info(splicingRmats)
        splicingRmats_execLine = subprocess.run((splicingRmats),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        write_subprocess_log(splicingRmats_execLine,logger)
        
        proc_find_files_project = subprocess.getoutput(("find "+config.parameters['path_to_output']+analysis+" -not -name *edgeCount -name "+"*.bam*"))
        proc_find_files_project = proc_find_files_project.split("\n")
        for aFile in proc_find_files_project:
            logger.info("Delete bam file from rmats analysis :")
            subprocess.run("rm "+aFile,shell=True)
            logger.info(aFile)
            
        proc_find_files_project = subprocess.getoutput(("find "+config.parameters['path_to_output']+analysis+" -not -name *edgeCount -name "+"*.sam*"))
        proc_find_files_project = proc_find_files_project.split("\n")
        for aFile in proc_find_files_project:
            logger.info("Delete sam file from rmats analysis :")
            subprocess.run("rm "+aFile,shell=True)
            logger.info(aFile)   
            
            
    logger.info("End")

