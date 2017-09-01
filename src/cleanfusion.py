##!/usr/bin/python3.5
##!/home/jean-philippe.villemin/bin/python-3.5-2/bin/python3

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
import re
import pprint
import string
import statistics
from math import log2
import numpy
import scipy.stats as stats
import pandas as pd
import sys
from collections import OrderedDict
from pyliftover import LiftOver
import math
import numpy as np
from  more_itertools import unique_everseen
#import pybedtools
#from pybedtools import bedtool

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
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+"/"+config.chrono+"/"+'activity.log', 'a', 1000000, 1)
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
    Thi script will merge all intermediate outputs.  
    Example : 
    python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/cleanfusion.py  -s1 Early_PolyA.bed -s2 Early_Ribo0.bed -p /home/jean-philippe.villemin/RNASEQ_2017_RESULTS/RMATS/28_8_2017__10_14_4/test/

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    #-c /home/jean-philippe.villemin/code/RNA-SEQ/configs/GHRC38/RMATS/MERGED/SE.splicing.MERGED.json
    #parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')
    parser.add_argument("-s1","--suffix1",action="store", help="suffix first file to read",required=True,type=str,dest='suffix1')
    parser.add_argument("-s2","--suffix2",action="store",help="suffix second file to read",required=True,type=str,dest='suffix2')
    parser.add_argument("-p","--path",action="store",help="path to file",required=True,type=str,dest='path')
    parser.add_argument("-o","--out",action="store",help="path to output file",required=True,type=str,dest='out')


    parameters = parser.parse_args()


    proc_find = subprocess.getoutput(("find "+parameters.path+" -name '*"+parameters.suffix1+"' -o  -name '*"+parameters.suffix2+"'")) # --ignoreDuplicates
    proc_find = proc_find.split("\n")
    print(proc_find)
    # report interact of a with b, and a with nothing
    proc_bedtools_1 = subprocess.run(("bedtools intersect -a "+ proc_find[0]+" -b "+ proc_find[1]+" -loj  > "+parameters.path+parameters.out+".interA.bed"),shell=True) # --ignoreDuplicates
    # report b that interact with anything
    proc_bedtools_2 = subprocess.run(("bedtools intersect -b "+proc_find[0]+" -a "+proc_find[1]+" -loj  > "+parameters.path+parameters.out+".interB.bed"),shell=True) # --ignoreDuplicates
    proc_bed_cat = subprocess.run(("cat "+parameters.path+parameters.out+".interA.bed "+parameters.path+parameters.out+".interB.bed > "+parameters.path+parameters.out+".inter.bed "),shell=True) # --ignoreDuplicates

    output=parameters.path+parameters.out+".inter.bed"
    
    unique_hash_event = {}
    all_events = []
    with open(output) as exon_line:
        for line in exon_line: 
            elements = line.strip().split("\t")
            # name, signal ,other start,other end, other name, other signal, strand
            if (elements[10] != "-1" ) :
                if (float(elements[4]) > 0 and float(elements[10]) < 0 ) or (float(elements[4]) < 0 and float(elements[10]) > 0 )   :
                    continue# Sorry soooo ugly...
                unique_hash_event[elements[0]+"_"+elements[1]+"_"+elements[2]]= [elements[0],elements[1],elements[2],elements[3],elements[4],elements[5],"MULTI",elements[7],elements[8],elements[9],elements[10]]
            else :
                unique_hash_event[elements[0]+"_"+elements[1]+"_"+elements[2]]= [elements[0],elements[1],elements[2],elements[3],elements[4],elements[5],"UNIQ"]

        # indice - value
        # 0 chr
        # 1 start
        # 2 end
        # 3 name
        # 4 signal
        # 5 strand
        # 6 start
        # 7 end
        # 8 name
        # 9 signal
        all_events.append( line );

    exon_line.close()    
    
    final = parameters.path+parameters.out+".final.bed"
    myfile  = open(final,"w")

    mylistgene= []
    for key, value in unique_hash_event.items():

        
        chrom    = value[0]
        start    = value[1]
        end      = value[2]
        name     = value[3]
        signal   = value[4]
        strand   = value[5]
        

        
        # Print directly case with no overlap
        if (value[6] == "UNIQ") :
            
            myfile.write (chrom+"\t"+start+"\t"+end+"\t"+"_".join(name.split("_")[0:2])+"\t"+signal+"\t"+strand+"\n")
            mylistgene.append (name.split("_")[1]+"\t"+name.split("_")[2]+"\n")
            continue
        
        # When A & B have intersection
        for line in all_events :
            elements = line.split("\t")
            test_id_line= "_".join(elements[0:3])

            if (test_id_line == key) :
                if (value[2] == '-1'):
                    break

                if (abs(float(elements[9])) > abs(float(signal) )) :
                    start   = elements[7]
                    end     = elements[8]
                    name    = elements[9]
                    signal  = elements[10]

        myfile.write (chrom+"\t"+start+"\t"+end+"\t"+"_".join(name.split("_")[0:2])+"\t"+signal+"\t"+strand+"\n")
        mylistgene.append (name.split("_")[1]+"\t"+name.split("_")[2]+"\n")


    myfile.close()
    
    finalgene = parameters.path+parameters.out+".final.gene.txt"
    myfilegene  = open(finalgene,"w")
    
    for mygene in  np.unique(mylistgene) : 
        myfilegene.write (mygene)   
    myfilegene.close()
    
    proc_sort = subprocess.run(("sort -k1,1V -k2,2n  "+final+ " > "+parameters.path+parameters.out+".final.sorted.bed"),shell=True, check=True) # --ignoreDuplicates
    proc_merge = subprocess.run(("bedtools merge -i "+parameters.path+parameters.out+".final.sorted.bed "+ "-c 4,5,6 -o collapse,mean,collapse > "+parameters.path+parameters.out+".final.sorted.merged.bed"),shell=True, check=True) # --ignoreDuplicates
 
    result = parameters.path+parameters.out+".bed"
    myresultfile  = open(result,"w")
 
    mylastfile  = parameters.path+parameters.out+".final.sorted.merged.bed"
    with open(mylastfile) as lines:
        for line in lines: 
            elements = line.split("\t")
            #chr6    31972621    31972919    -0.294_STK19,-0.301_STK19    -0.2975    +,+
            subelement=elements[3].split(",")
            if (len(subelement)>1) :
                symbol=subelement[0].split("_")[1]
                myresultfile.write(elements[0]+"\t"+elements[1]+"\t"+elements[2]+"\t"+elements[4]+"_"+symbol+"\t"+elements[4]+"\t"+elements[5].split(",")[0]+"\n")
            else : 
                myresultfile.write(line)

    lines.close()       
    myresultfile.close()
    
    #proc = subprocess.run(("sort -k1,1n -k2,2n "+parameters.path+"final.bed > "+parameters.path+"final.sorted.bed"),shell=True) # --ignoreDuplicates
    #proc = subprocess.run(("bedtools sort -i "+parameters.path+"final.bed > "+parameters.path+"final.sorted.bed"),shell=True) # --ignoreDuplicates


    #b'./SE_reads_5_Early_PolyA.bed\n./SE_reads_5_Early_Ribo0.bed\n'
    #bedtools intersect -b SE_reads_5_Early_PolyA.bed -a SE_reads_5_Early_Ribo0.bed -loj > loj_inv.bed
    #output_bedgraph=path_with_time+"/"+group_pair+"."+region_for_file+".bedgraph"
    #proc_bam_compare = subprocess.run(("bamCompare --bamfile1 "+config.parameters['files'][group_pair]['test']+" --bamfile2 "+config.parameters['files'][group_pair]['control']+" --normalizeUsingRPKM --ratio subtract --binSize 1 --outFileFormat bedgraph -o "+output_bedgraph+" --region "+region+" --numberOfProcessors 12 "),shell=True) # --ignoreDuplicates
        
    
