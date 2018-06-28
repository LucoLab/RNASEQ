
'''

:date: Jan 23, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Merge Expression and Splicing

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import glob
import glob
import pandas as pd
from collections import OrderedDict
import os
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
                
def create_logger(path,LEVEL):
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
    script = "diffExp_activity"
    print(path+"/"+script+'.log')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(path+"/"+script+'.log', 'a', 1000000, 1)
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
    Thi script will merge all intermediate outputs.  
    Example : 
    python3 ./diffGeneExp.py -c test.json -p Cond.TEST_vs_Cond.CONTROL

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')
    parser.add_argument("-p","--conditionPrevalenceList",action="store",help="Order to make the comparison",required=True,type=str,dest='conditionPrevalenceList')
    parser.add_argument("-i","--init",action="store",help="Path to a json file.",required=False,default=str(Path(os.path.dirname(os.path.realpath(__file__))).parents[0])+"/config/init.json",type=str,dest='init')

    parameters = parser.parse_args()
    
    init = custom_parser.Configuration(parameters.init,"json")

    config = custom_parser.Configuration(parameters.file_config,"json")
    #if not os.path.exists(config.parameters['path_to_output']):
        #os.makedirs(config.parameters['path_to_output'])
        
    orderedConditionList = parameters.conditionPrevalenceList.split("_")
    logger = create_logger(os.path.dirname(parameters.file_config ),"INFO")

    logger.info('=========> START:')
    logger.info("orderedConditionList : "+str(orderedConditionList))

    logger.info("file_config : "+parameters.file_config)
    logger.info("path_to_output : "+config.parameters['path_to_output'])
    logger.info("lib-type : "+config.parameters.get("lib-type"))
    logger.info("project : "+config.parameters.get("project"))


    # SEE CONF HTSEQ-COUNT OR STAR
    index=""
    if(config.parameters.get("lib-type")=="Reverse") :
        index=2
    if(config.parameters.get("lib-type")=="Forward") :#Yes
        index=1
    if(config.parameters.get("lib-type")=="Unstranded") :
        index=0

    
    matrix=OrderedDict()

    result = config.parameters.get("path_to_output")+"output/"+config.parameters.get("project")+"/Design.csv"
    
    mydesign  = open(result,"w")
    mydesign.write("sample_id,condition,emt\n")
    
    conditions = []
    for sample in config.parameters.get("files"):
       
        logger.info(sample.get("sample_id"))
        logger.info(sample.get("condition"))
        logger.info(sample.get("emt"))

        conditions.append(sample.get("condition"))
        logger.info(config.parameters.get("path_to_input")+"output/"+config.parameters.get("project")+"/"+sample.get("sample_id")+"/"+"star_output/"+sample.get("sample_id")+"_ReadsPerGene.tab")
        file = glob.glob(config.parameters.get("path_to_input")+"output/"+config.parameters.get("project")+"/"+sample.get("sample_id")+"/"+"star_output/"+sample.get("sample_id")+"_ReadsPerGene.tab")
        logger.info(file)
        if len(file) > 1 : 
            logger.info("error")
            break
        logger.info(file[0])
        data = pd.read_csv(file[0],sep='\t', index_col='Gene' )
        mydesign.write(sample.get("sample_id")+","+sample.get("condition")+","+str(sample.get("emt"))+"\n")

        data.sort_index(inplace=True)


        matrix[sample.get("sample_id")]=data.iloc[:,index]
    
    mydesign.close()
    pd.DataFrame(matrix)
    pd.DataFrame(matrix).to_csv(config.parameters.get("path_to_output")+"output/"+config.parameters.get("project")+"/Raw_read_counts.csv")
    
    #list_of_pairs = [(p1, p2) for p1 in list(set(conditions)) for p2 in list(set(conditions)) if p1 != p2]
    list_of_pairs = [(conditions[p1], conditions[p2]) for p1 in range(len(conditions)) for p2 in range(p1+1,len(conditions))]
    list_of_pairs = list(set(list_of_pairs))
    
    logger.info(list_of_pairs)

    logger.info( config.parameters.get("organism") )
    script=""
    if(config.parameters.get("organism")=="human" or config.parameters.get("organism")=="homo_sapiens") :
        script="diff_exp.R"
    if(config.parameters.get("organism")=="mouse" or config.parameters.get("organism")=="mus_musculus")   :
        script="diff_exp_mouse.R"
    
    logger.info("-> DIFF EXPRESSION")

    for i in list_of_pairs :
        
        logger.info(i)

        if (i[0] != i[1]) :
            
            test    = i[0] if "TEST" in i[0] else i[1]
            control = i[1] if "CONTROL" in i[1] else i[0]
  
            logger.info("-> "+ test+"  VS "+control)

            logger.info(("\n Rscript "+init.parameters['scriptDir']+"Rscript/"+script+" " \
                +"  --dir "+config.parameters.get("path_to_output")+""+test+"_vs_"+control    \
                +" "+config.parameters.get("path_to_output")+"output/"+config.parameters.get("project")+"/Design.csv"     \
                +" --cond1 "+test   \
                +" --cond2 "+control    \
                +" "+config.parameters.get("path_to_output")+"output/"+config.parameters.get("project")+"/"+"Raw_read_counts.csv \n"))
                        
            proc = subprocess.run(("Rscript "+init.parameters['scriptDir']+"Rscript/"+script+" " \
                +"  --dir "+config.parameters.get("path_to_output")+""+test+"_vs_"+control    \
                +" "+config.parameters.get("path_to_output")+"output/"+config.parameters.get("project")+"/Design.csv"     \
                +" --cond1 "+test    \
                +" --cond2 "+control    \
                +" "+config.parameters.get("path_to_output")+"output/"+config.parameters.get("project")+"/Raw_read_counts.csv") \
               ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            write_subprocess_log(proc,logger)
    
    logger.info("- > End diffexp")
   
