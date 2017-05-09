
'''

:date: Apr 04, 2017
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Use Salmon on Fastq to compute TPM values for genes

'''
import argparse,textwrap

from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import datetime

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
                

def core_sample_name(config,group_pair):
    '''
    Clean up the names of the paired sample.
    Remove fastq.gz extension.
    Remove '_R1_' to have the same core syntax to use with '_R2_'.
  
    Args:
        config (obj): Configuration instance.
        group_pair (str): The pair identifier of fastq sequence to trait.

    Returns:
        clean_name (string) : Clean core name of sample without R1 / R2 disctinction.
    
    
    '''
    clean_name =  config.parameters["files"][group_pair]["R1"].replace(".fastq.gz","",1)
    clean_name =  clean_name.replace("_R1_","",1)
    
    return clean_name


def create_logger(config):
    """
    Define a logger instance to write in a file and stream.
  
    Args:
        config (obj): Configuration instance.

    Returns:
        logger (obj): Logger instance to log messages in file and output mainstream.
    
        
    """
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+'activity.log', 'a', 1000000, 1)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)

    return logger


###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    1- Concat all fastqs  per R1 and R2 group separetely.  
    2- Make salmon count TPM values by genes
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")


    subprocess.run(("mkdir -p "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    print("mkdir -p "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono)
    
    logger = create_logger(config)

    logger.info('=========> START:PARSE EACH FASTQ PAIRS FOR AN EXPERIMENT')
   
    list_core_sample_names_for_R1  = []
    list_core_sample_names_for_R2  = []   
    
    # Reads fastq of your analysis
    for group_pair in config.parameters.get("files").keys():
       
        logger.info(group_pair)
        logger.info (config.parameters["files"][group_pair]["R1"])
        logger.info (config.parameters["files"][group_pair]["R2"])
        logger.info ("")
        
        list_core_sample_names_for_R1.append(config.parameters['path_to_input']+'/'+config.parameters["files"][group_pair]["R1"])
        list_core_sample_names_for_R2.append(config.parameters['path_to_input']+'/'+config.parameters["files"][group_pair]["R2"])
        
    ####################### concat R1 fastq #######################
    tr_files_1 = ' '.join(list_core_sample_names_for_R1)
    logger.info("cat "+tr_files_1+" > "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_1.fastq.gz")
    
    proc_cat1= subprocess.run(("cat "+tr_files_1+" > "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_1.fastq.gz")
    ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
    write_subprocess_log(proc_cat1,logger)
    
    
    ####################### concat R2 fastq #######################
    tr_files_2 = ' '.join(list_core_sample_names_for_R2)
    logger.info("cat "+tr_files_2+" > "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_2.fastq.gz")
    
    proc_cat2= subprocess.run(("cat "+tr_files_2+" > "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_2.fastq.gz"
    ),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        
    write_subprocess_log(proc_cat2,logger)
    
    
    ####################### Salmonize ##########################
    logger.info("salmon quant -i "+config.parameters['transcriptome_index'] \
        +" -p " +config.parameters['salmon']['cpu'] \
        +" -l A " \
        +" -g "+config.parameters['path_to_gtf'] \
        +" -1 "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_1.fastq.gz" \
        +" -2 "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_2.fastq.gz" \
        +" -o "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.parameters['final_bam_name']+"_quant")
     
    salmon= subprocess.run("salmon quant -i "+config.parameters['transcriptome_index'] \
        +" -p " +config.parameters['salmon']['cpu'] \
        +" -l A " \
        +" -g "+config.parameters['path_to_gtf']  \
        +" -1 "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_1.fastq.gz" \
        +" -2 "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_2.fastq.gz" 
        +" -o "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.parameters['final_bam_name']+"_quant"
    ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
           
    write_subprocess_log(salmon,logger)
    timeline = datetime.datetime.now()
    endtime = str(timeline.hour)+"_"+str(timeline.minute)+"_"+str(timeline.second)

  
    logger.info('=========> END '+endtime)
    logger.info('=========> Finish !')
