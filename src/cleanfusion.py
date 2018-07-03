'''

:date: Jan 23, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Merge Expression and Splicing

'''
import argparse,textwrap
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import numpy as np
import os
import statistics
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
                
def create_logger(path,type,LEVEL):
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
    script = "cleanFusion_activity"
    print(path+"/"+type+"."+script+'.log')
    ''' 1 st handler in file'''
    file_handler = RotatingFileHandler(path+"/"+type+"."+script+'.log', 'a', 1000000, 1)
    file_handler.setLevel(LEVEL)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(LEVEL)
    logger.addHandler(stream_handler)

    return logger


def selectTheHighestScore(elements) :
    
    selectedEvent = ""
    #chr1    77947452    77947552    77947452,77947452    77947547,77947552    0.223_FUBP1_ENSG00000162613_0.49_22,0.225_FUBP1_ENSG00000162613_0.49_22    0.223,0.225    -,-
    arrayOfSignal = elements[6].strip().split(",")

    #logger.info(arrayOfSignal)

    countSameSignal = 0
    sens = ""
    dpsi_retained = 0
    
    for signal in arrayOfSignal :
        if (float(signal) > 0 ) :
            countSameSignal+=1
            sens="pos"
        if (float(signal) < 0 ) :
            countSameSignal+=1
            sens="neg"
    #logger.info(countSameSignal)

    if(countSameSignal!=len(arrayOfSignal)  ):
        logger.info("SIGNAL DIFFERENT SIGN"+(",".join(elements)))
        return "NA"

    elif ( sens == "neg" ) :
        dpsi_retained =  min(map(float, arrayOfSignal))
    elif  ( sens == "pos" ):
        dpsi_retained =  max(map(float, arrayOfSignal))
 
    index = arrayOfSignal.index(str(dpsi_retained))
    #logger.info(index)

    selectedEvent = "\t".join([elements[0],elements[3].split(",")[index] ,elements[4].split(",")[index],elements[5].split(",")[index],elements[6].split(",")[index],elements[7].split(",")[index]])+"\n"
    #logger.info(selectedEvent)

    return selectedEvent 
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will merge all intermediate outputs.  
    Example : 
    python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/cleanfusion.py  -s1 TMP -p ${FINAL_OUT} -o ${ID1} -ct ${ISCONTROL}

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s1","--suffix1",action="store", help="suffix first file to read",required=True,type=str,dest='suffix1')
    parser.add_argument("-p","--path",action="store",help="path to file",required=True,type=str,dest='path')
    parser.add_argument("-o","--out",action="store",help="path to output file",required=True,type=str,dest='out')
    parser.add_argument("-ct","--control",action="store",help="control or not",required=True,type=str,dest='iscontrol')

    parameters = parser.parse_args()
   

    path_above = Path(os.path.dirname(parameters.path)).parents[1]
   
    type = os.path.basename(os.path.normpath(parameters.path+"/"))

    logger = create_logger(str(path_above),type,"INFO")
    logger.info("type : "+str(type))

    logger.info("Log is here : "+str(path_above))
    logger.info(os.path.dirname(parameters.path))

    logger.info("find "+parameters.path+" -name '*"+parameters.suffix1+"*'")
    proc_find = subprocess.getoutput(("find "+parameters.path+" -name '*"+parameters.suffix1+"*'"))
    
    logger.info(proc_find)
    logger.info(parameters.iscontrol)

    proc_find = proc_find.split("\n")

    final = proc_find[0]
    #"-s "
    proc_sort = subprocess.run(("sort -k1,1V -k2,2n  "+final+ " > "+parameters.path+parameters.out+".final.sorted.bed"),shell=True, check=True)
    write_subprocess_log(proc_sort,logger)
    
    proc_count = subprocess.getoutput(("wc -l  "+parameters.path+parameters.out+".final.sorted.bed"))
    logger.info("Lines before merged "+str(proc_count.split(" ")[0]))

    proc_merge = subprocess.run(("bedtools merge -i "+parameters.path+parameters.out+".final.sorted.bed "+ "  -c 2,3,4,5,6 -o collapse,collapse,collapse,collapse,collapse > "+parameters.path+parameters.out+"_TMP_final.sorted.merged.bed"),shell=True, check=True) 
    write_subprocess_log(proc_merge,logger)

    result = parameters.path+parameters.out+".bed"
    myresultfile  = open(result,"w")
    mylistgene= []

    mylastfile  = parameters.path+parameters.out+"_TMP_final.sorted.merged.bed"
    countAllLinesMerged = 0
    countLinesOverlaped = 0
    
    complexEvent =  parameters.path+parameters.out+".complex.txt"
    #garbage= open(complexEvent,"w") 

    with open(mylastfile) as lines:
        for line in lines: 
            countAllLinesMerged+=1

            elements = line.split("\t")
            #logger.info(elements)
            #chr1    6886239    6886346    6886239    6886346    0.282_CAMTA1_ENSG00000171735_0.09_8    0.282    +
            if ',' not in line:

                myresultfile.write("\t".join([elements[0],elements[1],elements[2],elements[5],elements[6],elements[7]]))
                ensembl=elements[5].split("_")[2]
                name=elements[5].split("_")[1]

            else :
            #chr1    77947452    77947552    77947452,77947452    77947547,77947552    0.223_FUBP1_ENSG00000162613_0.49_22,0.225_FUBP1_ENSG00000162613_0.49_22    0.223,0.225    -,-
                selectedEvent = selectTheHighestScore(elements)
                if (selectedEvent=="NA") : 
                    logger.info("NA")

                    continue
                #logger.info("YourAreIn")

                myresultfile.write(selectedEvent)
                elementsBis = selectedEvent.split("\t")
                ensembl=elementsBis[3].split("_")[2]
                name=elementsBis[3].split("_")[1]
                #chr1    6886239   6886346    0.282_CAMTA1_ENSG00000171735_0.09_8    0.282    +
                countLinesOverlaped+=1
                #garbage.write(line)
 
            mylistgene.append (name+"\t"+ensembl+"\n")
    #garbage.close()
    lines.close()       
    myresultfile.close()
    logger.info("Lines after merged,total not overlapping "+str(countAllLinesMerged))
    logger.info("Lines arranged because of merged "+str(countLinesOverlaped))
    logger.info("Number of genes seen : "+str(len(np.unique(mylistgene))))
    # If Control , dont rewrite gene list.
    if (parameters.iscontrol!="True") :

        finalgeneSymbol = parameters.path+parameters.out+".gene.symbol.txt"
        finalgeneEnsembl = parameters.path+parameters.out+".gene.ensembl.txt"

        myfilegeneSymbol  = open(finalgeneSymbol,"w")
        myfilegeneEnsembl  = open(finalgeneEnsembl,"w")

        for mygene in  np.unique(mylistgene) : 
            myfilegeneSymbol.write (mygene.split("\t")[0]+"\n")
            myfilegeneEnsembl.write (mygene.split("\t")[1])     
        myfilegeneEnsembl.close()
        myfilegeneSymbol.close()
    logger.info("Clean Fusion finished !")

