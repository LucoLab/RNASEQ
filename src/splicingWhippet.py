
'''

:date: September 9, 2017
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Launch Splicing with Whippet

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import os 
from pathlib import Path


###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################

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
    script = "whippet_activity"
    print("Log is : "+path+"/"+script+'.log')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(path+"/"+script+".log", 'a', 1000000, 1)
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
    python ./splicingRmats.py -c /home/jp/workspace/RNA-SEQ/bash/rmats.json

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')
    parser.add_argument("-a","--analyse",action="store",help="Type of analyse you want to execute.",default="dpsi",required=False,type=str,dest='analyse')
    parser.add_argument("-k","--kmer",action="store",help="KmerSize Index to use.",default="normal",required=False,type=str,dest='kmer')
    parser.add_argument("-i","--init",action="store",help="Path to a json file.",required=False,default=str(Path(os.path.dirname(os.path.realpath(__file__))).parents[0])+"/config/init.json",type=str,dest='init')

    #parser.add_argument("-a","--analysis",action="store",help="Name of analysis",required=True,type=str,dest='analysis')

    parameters = parser.parse_args()
    print(parameters.file_config)
    
    config = custom_parser.Configuration(parameters.file_config,"json")
    
    init = custom_parser.Configuration(parameters.init,"json")

    
    path_above = Path(os.path.dirname(config.parameters['path_to_output'])).parents[2]

    logger = create_logger(str(path_above),"INFO")

    print('=========> START:')
    print("file_config : "+parameters.file_config)
    print("Kmer size : "+parameters.kmer)

    #print("analysis : "+parameters.analysis)

    print("path_to_output : "+config.parameters['path_to_output'])
    print("    ")

    
    index_path="None"
    if (config.parameters['organism']=="human") :
        if (parameters.kmer=="short") :
            logger.info("Kmer Index used is short ")
            index_path=init.parameters['index_whippet_short'] #/human_index_whippet
        if (parameters.kmer=="normal") :
            logger.info("Kmer Index used is default ")
            index_path=init.parameters['index_whippet_normal'] #/human_index_whippet
    #if (config.parameters['organism']=="mouse") :
        #index_path="/home/jean-philippe.villemin/data/data/index_whippet/mouse/mouse_index_whippet"
 
    logger.info(config.parameters['organism'])    

    logger.info(index_path)    
    '''
    ######################################################################
    DataFrame Object Construction : First Part
    ######################################################################
    '''
    for analyse_key,analyse_values in config.parameters["analysis"].items(): 
       
        #if(analyse_key==parameters.analysis) :
        logger.info(analyse_key)
        for element_key,element_value in analyse_values.items() : 
            logger.info("-> "+element_key+" :  "+element_value)
            for sample_group in config.parameters["files"][element_value] : 
                logger.info("====>")
                logger.info(element_value)
        
                for hash_key,hash_values in sample_group.items() : 
                        
                        # you don't to redo something already done.
                        if(os.path.isfile(config.parameters['path_to_output']+hash_key+".psi.gz") or os.path.isfile(config.parameters['path_to_output']+hash_key+".psi")) : 
                                logger.info("No need to process again. File "+hash_key+".psi.gz already exist")
                                continue
                                
                        logger.info(hash_key)
                        logger.info(hash_values.get("R1"))
                        logger.info(hash_values.get("R2"))
                        if(hash_values.get("R2")=="None" ) : 
                            if(config.parameters['path_to_input']=="None") :
                                
                                execLine = init.parameters['pathToJulia']+" "+init.parameters['pathToWhippet']+"whippet-quant.jl -x "+index_path+" --url http://"+hash_values.get("R1")[6:]+" -o "+config.parameters['path_to_output']+hash_key
                            else :
                                execLine = init.parameters['pathToJulia']+" "+init.parameters['pathToWhippet']+"whippet-quant.jl -x "+index_path+" "+config.parameters['path_to_input']+hash_values.get("R1")+" -o "+config.parameters['path_to_output']+hash_key
                        else :
                            if(config.parameters['path_to_input']=="None") : 
                                execLine = init.parameters['pathToJulia']+" "+init.parameters['pathToWhippet']+"whippet-quant.jl -x "+index_path+"  --url http://"+hash_values.get("R1")[6:]+" http://"+hash_values.get("R2")[6:]+" -o "+config.parameters['path_to_output']+hash_key
                            else :
                                execLine = init.parameters['pathToJulia']+" "+init.parameters['pathToWhippet']+"whippet-quant.jl -x "+index_path+"  "+config.parameters['path_to_input']+hash_values.get("R1")+" "+config.parameters['path_to_input']+hash_values.get("R2")+" -o "+config.parameters['path_to_output']+hash_key

                        logger.info(execLine)
        
                        whippetquant = subprocess.run((execLine),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)                
                        write_subprocess_log(whippetquant,logger)
                            
    
    # Never used need to dev something here...
    if  (parameters.analyse =="psi" ) :   
    
        logger.info("START TO PARSE ANALYSIS FOR PSI ====>")

        for analyse_key,analyse_values in config.parameters["analysis"].items(): 
                #if(analyse_key==parameters.analysis) :
            logger.info(analyse_key) 
            logger.info("TEST :    "+analyse_values.get("TEST"))
            for sample_groupTest in config.parameters["files"][analyse_values.get("TEST")] : 
                for hash_keyTest,hash_valuesTest in sample_groupTest.items() : 
                    logger.info(hash_keyTest) 
                    
            logger.info("CONTROL :    "+analyse_values.get("CONTROL"))
            for sample_groupCon in config.parameters["files"][analyse_values.get("CONTROL")] : 
                for hash_keyCon,hash_valuesCon in sample_groupCon.items() : 
                    logger.info(hash_keyCon)         
                    
    if  (parameters.analyse =="dpsi" ) :   

        logger.info("START TO PARSE ANALYSIS FOR DPSI ====>")

        for analyse_key,analyse_values in config.parameters["analysis"].items(): 
                #if(analyse_key==parameters.analysis) :
            logger.info(analyse_key) 
            logger.info("TEST :    "+analyse_values.get("TEST"))
            test = [] 
            trickForOneSample =""
            extension = ".psi.gz"
            
            ### Sorry mum , it's so ugly, but that for correction of the bug that need a comma & no extension after sample name if only one repicate
            tcheckBug=0
            for sample_groupTest in config.parameters["files"][analyse_values.get("TEST")] : 
                for hash_keyTest,hash_valuesTest in sample_groupTest.items() :
                    tcheckBug+=1
            print(tcheckBug)
            if(tcheckBug==1):
                #extension = ""
                trickForOneSample = ","
    
            for sample_groupTest in config.parameters["files"][analyse_values.get("TEST")] : 
                for hash_keyTest,hash_valuesTest in sample_groupTest.items() : 
                    logger.info(hash_keyTest) 
                    test.append(config.parameters['path_to_output']+hash_keyTest+extension) #.gz
            logger.info("CONTROL :     "+analyse_values.get("CONTROL"))
    
            contr = [] 
            for sample_groupCon in config.parameters["files"][analyse_values.get("CONTROL")] : 
                for hash_keyCon,hash_valuesCon in sample_groupCon.items() : 
                    logger.info(hash_keyCon) 
                    contr.append(config.parameters['path_to_output']+hash_keyCon+extension)   #.gz
              
            skip=0    
            if(os.path.isfile(config.parameters['path_to_output']+analyse_key+".diff.gz") or os.path.isfile(config.parameters['path_to_output']+analyse_key+".diff")) : 
                logger.info("No need to process again. File "+hash_key+".diff.gz already exist")
                skip=1
           
            if(skip==0) :
                logger.info("DELTA ====>")
                execLine_delta = init.parameters['pathToJulia']+" "+init.parameters['pathToWhippet']+"whippet-delta.jl -a "+",".join(test)+trickForOneSample+" -b "+",".join(contr)+trickForOneSample+" -o "+config.parameters['path_to_output']+analyse_key
                logger.info(execLine_delta)
                whippetdelta = subprocess.run((execLine_delta),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                write_subprocess_log(whippetdelta,logger)
                logger.info(" ")
        
                logger.info("GUNZIP ====>")
        
                gunzip_command = "gunzip  -k "+config.parameters['path_to_output']+analyse_key+".diff.gz"
                logger.info(gunzip_command)
        
                gunzip = subprocess.run((gunzip_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                write_subprocess_log(gunzip,logger)
            
            logger.info(" ")
            logger.info("FILTER ALL EVENTS ====>")
    
            for event in ["CE","AA","AD","RI"] :
                
                if(os.path.isfile(config.parameters['path_to_output']+analyse_key+".clean."+event+".diffannoted.csv") and os.path.isfile(config.parameters['path_to_output']+analyse_key+".bad."+event+".diffannoted.csv") ) : 
                    logger.info("Already annotated and filtered -> "+config.parameters['path_to_output']+analyse_key+".clean."+event+".diffannoted.csv")
                    logger.info("Already annotated and filtered -> "+config.parameters['path_to_output']+analyse_key+".bad."+event+".diffannoted.csv")
    
                    continue
    
                filter_command = config.parameters["path_to_cleaner"]+" "+ event +" "+analyse_key+" "+config.parameters['path_to_output']
                logger.info("EVENT : "+event)
                logger.info(filter_command)
                filter = subprocess.run((filter_command),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                write_subprocess_log(filter,logger)
                
                Rcommand ="Rscript "+init.parameters['scriptDir']+"Rscript/annotSyymbol.R --organism="+config.parameters['organism']+" --file="+config.parameters['path_to_output']+analyse_key+".clean."+event+".diff"
                logger.info(Rcommand)
                R = subprocess.run((Rcommand),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                write_subprocess_log(R,logger)
                
                Rcommand2 ="Rscript "+init.parameters['scriptDir']+"Rscript/annotSyymbol.R --organism="+config.parameters['organism']+" --file="+config.parameters['path_to_output']+analyse_key+".bad."+event+".diff"
                logger.info(Rcommand2)
                R2 = subprocess.run((Rcommand2),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                write_subprocess_log(R2,logger)
    
        logger.info("End")


