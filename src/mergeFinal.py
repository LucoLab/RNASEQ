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
#logging.getLogger('requests').setLevel(logging.ERROR)

import re
import statistics
import numpy
import scipy.stats as stats
import pandas as pd
import math
import numpy as np
import itertools
import os
import requests
import json
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
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+"/"+id+"."+'mergeFinal_activity.log', 'a', 1000000, 1)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logger.addHandler(stream_handler)

    return logger


def parse_rmats(lineElements,chrom,strand,as_type,prefix):
    
    exon_size =  "."
    trickMXE     =  ""

    gonext=""
    if (as_type == "A3SS") : 
        if (strand == "-") :
            id_ucsc_event = chrom+":"+lineElements[5]+"-"+lineElements[10]+":"+strand  #flankingES    flankingEE
            highlight = prefix+"."+chrom+":"+str(int(lineElements[8])+1)+"-"+lineElements[6]
            if ( int(lineElements[8]) - int(lineElements[6]) >= 200 ) : # Correct bug A3SS , A5SS. Sometimes it uses intron. Correted by length selection.
                return "next",id_ucsc_event,highlight,exon_size,trickMXE
  
        if (strand == "+") :
            id_ucsc_event = chrom+":"+lineElements[9]+"-"+lineElements[6]+":"+strand  #flankingES    flankingEE
            highlight = prefix+"."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[7]
            if ( int(lineElements[7]) - int(lineElements[5]) >= 200 ) : # Correct bug A3SS , A5SS. Sometimes it uses intron. Correted by length selection.
                return "next",id_ucsc_event,highlight,exon_size,trickMXE
    
    if (as_type == "A5SS") :  
        if (strand == "-") :
            id_ucsc_event = chrom+":"+lineElements[9]+"-"+lineElements[6]+":"+strand  #flankingES    flankingEE
            highlight = prefix+"."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[7]
            if ( int(lineElements[7]) - int(lineElements[5]) >= 200 ) : # Correct bug A3SS , A5SS. Sometimes it uses intron. Correted by length selection.
                return "next",id_ucsc_event,highlight,exon_size,trickMXE
  
        if (strand =="+") :
            id_ucsc_event = chrom+":"+lineElements[5]+"-"+lineElements[10]+":"+strand  #flankingES    flankingEE
            highlight = prefix+"."+chrom+":"+str(int(lineElements[8])+1)+"-"+lineElements[6]  
            if ( int(lineElements[6]) - int(lineElements[8]) >= 200 ) : # Correct bug A3SS , A5SS. Sometimes it uses intron. Correted by length selection.
                return "next",id_ucsc_event,highlight,exon_size,trickMXE
 
                
    if (as_type == "SE" ) :   
        id_ucsc_event = chrom+":"+lineElements[7]+"-"+lineElements[10]+":"+strand  #upstreamES downstreamEE
        exon_size = int(lineElements[6]) - int(lineElements[5]) # + 1 0-Based pr Rmats   example 0 99 
        highlight = prefix+"."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[6]
    
    if (as_type == "RI" ) :
        id_ucsc_event = chrom+":"+lineElements[5]+"-"+lineElements[6]+":"+strand  #upstreamES downstreamEE
        highlight = prefix+"."+chrom+":"+str(int(lineElements[8])+1)+"-"+lineElements[9]  
    
    if (as_type == "MXE" ) :
        
        id_ucsc_event = chrom+":"+lineElements[9]+"-"+lineElements[12]+":"+strand #upstreamES downstreamEE
        highlight = prefix+"."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[6]+"|"+prefix+"."+chrom+":"+str(int(lineElements[7])+1)+"-"+lineElements[8]
        trickMXE=lineElements[5]+"\t"+lineElements[7]
        

    return gonext,id_ucsc_event,highlight,exon_size,trickMXE


def filter_rmats(incLevel1,incLevel2,diffinc,ic_sample_1,sc_sample_1,ic_sample_2,sc_sample_2,readsNumber,isControl) : 
    
    ''' LOG RATIO // FOLD CHANGE '''
    incLevel1 = remove_values_from_list(incLevel1.split(","),"NA")
    incLevel2 = remove_values_from_list(incLevel2.split(","),"NA")
    
    average_incLevel1 =  statistics.mean(map(float, incLevel1))
    average_incLevel2 =  statistics.mean(map(float, incLevel2))
    
    max_incLevel1 =  max(map(float, incLevel1))
    max_incLevel2 =  max(map(float, incLevel2))
    
    max_incLevel   = max_incLevel1 if max_incLevel1 > max_incLevel2 else max_incLevel2
    if(max_incLevel==0):
        psi_classifier = "NaN"
    else :
        psi_classifier = abs(float(diffinc)/max_incLevel)
    
    logratio     = "NaN" 
    retval       = "NaN"
    #print('average'+str(average_incLevel1)+' '+str(average_incLevel2))
    if(average_incLevel1 != 0 and average_incLevel2 != 0): 
        retval   = foldchange (average_incLevel1,average_incLevel2) 
        logratio = foldchange2logratio(retval)
        #print('logratio'+str(logratio))
    
    # List of reads per group
    list_int_ic_sample1 = list(map(int,ic_sample_1.split(",")))
    list_int_sc_sample1 = list(map(int,sc_sample_1.split(",")))

    list_int_ic_sample2 = list(map(int,ic_sample_2.split(",")))
    list_int_sc_sample2 = list(map(int,sc_sample_2.split(",")))
    
    index_read_count_ic_sample1  = 0
    index_read_count_sc_sample2  = 0
    
    # TEST 1
    for read_count_ic1 in  list_int_ic_sample1  :
        if read_count_ic1 >= readsNumber :
            index_read_count_ic_sample1 = index_read_count_ic_sample1 + 1
    
    
    for read_count_sc2 in  list_int_sc_sample2  :
        if read_count_sc2 >= readsNumber :
            index_read_count_sc_sample2 = index_read_count_sc_sample2 + 1
        
    # More than a half has read count > X for INCLUSION in TEST and more than a half has read count > X for EXLUCSION in CONTROL
    test1 = 0 
    if ( (round(index_read_count_ic_sample1/len(list_int_ic_sample1),1) >=  0.5 ) and (round(index_read_count_sc_sample2/len(list_int_sc_sample2),1) >=  0.5 ) ):
        test1 = 1
    
    # TEST 2        
    index_read_count_sc_sample1  = 0  
    index_read_count_ic_sample2  = 0
    
    for read_count_sc1 in  list_int_sc_sample1  :
        if read_count_sc1 >= readsNumber :
            index_read_count_sc_sample1 = index_read_count_sc_sample1 + 1
    
    
    for read_count_ic2 in  list_int_ic_sample2  :
        if read_count_ic2 >= readsNumber :
            index_read_count_ic_sample2 = index_read_count_ic_sample2 + 1
            
    # More than a half has read count > X for EXCLUSION in TEST and more than a half has read count > X for INCLUSION in CONTROL
    test2 = 0 
    if ( (round(index_read_count_ic_sample2/len(list_int_ic_sample2),1) >=  0.5 ) and  (round(index_read_count_sc_sample1/len(list_int_sc_sample1),1) >=  0.5 ) ) :
        test2 = 1
    
    # IF TEST1 OR TEST2 failed , go next
    test = ""
    if(not isControl) :
        if ( test1 + test2 == 0 ) :
            test = "next"  
         
    ic_sample_1_sum = sum(list_int_ic_sample1)/len(ic_sample_1)
    sc_sample_1_sum = sum(list_int_sc_sample1)/len(sc_sample_1)
    ic_sample_2_sum = sum(list_int_ic_sample2)/len(ic_sample_2)
    sc_sample_2_sum = sum(list_int_sc_sample2)/len(sc_sample_2)
    
    ''' FISHER TEST '''
    pvalue_fisher = 0
    oddsratio, pvalue_fisher = stats.fisher_exact([[ic_sample_1_sum, sc_sample_1_sum], [ic_sample_2_sum, sc_sample_2_sum]])
    fisher = "-"
    if (pvalue_fisher < 0.05) : fisher = "PASS"

    return  test,average_incLevel1, average_incLevel2,logratio,psi_classifier,retval,pvalue_fisher,fisher 

def get_position_exons(ext) : 
    
    server = "https://rest.ensembl.org"

    nb_transcripts      = 0
    nb_exon_transcripts = 0
    position_list       = []
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
        r.raise_for_status() 
        print("WTF")                        

    params = json.loads(r.text)
    for item in params:
        for k, v in item.items():
            #logger.info( k+" "+str(v))
            if(k=="feature_type" and v=="transcript") : 
                nb_transcripts=nb_transcripts+1

            if(k=="rank") : 
                position_list.append(str(v))
                nb_exon_transcripts=nb_exon_transcripts+1
    
    # Remove if first or last
    skip = 0
    #decodedObjects = r.json()
    #for object in decodedObjects:
        #print(object)
    #    if(skip==1) : break
    #    for name,value in object.items() :
    #        if(name=="rank") : 
    #            if ( int(value) <= 0) :
    #                skip = 1
    #            break 
                        

        #if(skip==1) : continue
    
                    
    return nb_transcripts,nb_exon_transcripts,position_list,skip

def parse_all_splicing_files(path_to_dir,list_as_type,dict_for_analysis, readsNumber,gene_length,soft,list_replicates_test,list_replicates_control,organism,genes,isControl):
    """
    Parse all files produced by RMATS.
    All files have fixed names and in the same directories.
  
    Args:
        path_to_dir (str): Path to the directory.
        list_as_type (list): List of the as event to  retrieve the filenames in the directory
        readsNumber : Reads Number on junction for cutoff
        gene_length :  Dict Gene Name -> Length
    Returns:
        splicing_dict (obj): Return all the informations in a structured object.
    
        
    """
    prefix=""
    if(organism=="homo_sapiens") :
        prefix="hg38"
    if(organism=="mus_musculus") :
        prefix="mm10"
    
                    ############################################################# #############################################################
                    ###############           RMATS            ###############
                    ###############                               ###############
                    ############################################################# #############################################################
                    
    if (soft=="RMATS")  :
        for as_type in list_as_type :
            logger.info("LOOK SPLICING FOR AS_TYPE: "+as_type)
            array_splicing_data = { as_type :{} }
            namefile = ""
            if (parameters.modeCount == "2" ) :
                if(config.parameters["soft_version_rmats"]=="old"):
                    namefile = "MATS.JunctionCountOnly.txt"
                if(config.parameters["soft_version_rmats"]=="new"):
                    namefile = "MATS.JC.txt"
            else  :
                if(config.parameters["soft_version_rmats"]=="new"):
                    namefile = "MATS.JCEC.txt"
                if(config.parameters["soft_version_rmats"]=="old"):
                    namefile = "MATS.ReadsOnTargetAndJunctionCounts.txt"

                   
            with open(path_to_dir+"/"+as_type+"."+namefile) as f:
                count = 0
                for line in f:
                    #print(line)
                    if count == 0 : 
                        count = count +1
                        continue

                    id_ucsc_event = ""
                    lineElements = line.replace("\"", "").strip().split("\t")
                    m            = re.search('^(\w+)\.(\d+)$', lineElements[1])
                    ensembl_id   = m.group(1)
                    gene_size = gene_length[ensembl_id]
                    gene         = lineElements[2]
                    chrom        = lineElements[3]
                    strand       = lineElements[4]
                    ic_sample_1 = lineElements[12]
                    sc_sample_1 = lineElements[13]
                    ic_sample_2 = lineElements[14]
                    sc_sample_2 = lineElements[15]
                    incLevel1   = lineElements[20]
                    incLevel2   = lineElements[21]
                    diffinc     = lineElements[22]
                    pvalue      = lineElements[18]
                    fdr         = lineElements[19]
                    highlight   = ""
                
                    #############################################################
                    ###############           SPLICING            ###############
                    ###############                               ###############
                    #############################################################
                    if(not isControl):
                        if( ( (as_type != "MXE") and (abs(float(diffinc)) >= 0.2 ) ) or ( (as_type == "MXE") and (abs(float(lineElements[24])) >= 0.2) ) ) : 
                            
                            key_id_ucsc_event = chrom+":"+(":".join(lineElements[5:11]))+":"+strand  #flankingES    flankingEE   
                            gonext,id_ucsc_event,highlight,exon_size,trickMXE = parse_rmats(lineElements,chrom,strand,as_type,prefix)
                            if(gonext=="next"): continue

                            # Need to rewrite some stuffs because format is different for MXE !!
                            if (as_type == "MXE" ) :
        
                                key_id_ucsc_event = chrom+":"+(":".join(lineElements[5:13]))+":"+strand  #flankingES    flankingEE    
    
                                ic_sample_1 = lineElements[14]
                                sc_sample_1 = lineElements[15]
                                ic_sample_2 = lineElements[16]
                                sc_sample_2 = lineElements[17]
                                incLevel1   = lineElements[22]
                                incLevel2   = lineElements[23]
                                diffinc     = lineElements[24]
                                pvalue      = lineElements[20]
                                fdr         = lineElements[21]
                                
                            #logger.info(highlight[highlight.index(":"):][1:])
                            
                            # Avant j'avais mis and   avec le float 
                            if (float(fdr)   < 0.05 ) : 
                                    test,average_incLevel1, average_incLevel2,logratio,psi_classifier,retval,pvalue_fisher,fisher  = filter_rmats(incLevel1,incLevel2,diffinc,ic_sample_1,sc_sample_1,ic_sample_2,sc_sample_2,readsNumber,isControl)
                                    
                                    if(test=="next"): continue
                                    
                                    ext = "/overlap/region/"+config.parameters['organism']+"/"+chrom.replace("chr","")+":"+highlight[highlight.index(":"):][1:]+"?feature=exon;feature=transcript;content-type=application/json"

                                    # Too much time when a lot of exons
                                    if (as_type != "MXE" ) :

                                        nb_transcripts,nb_exon_transcripts,position_list,skip = get_position_exons(ext)
                                    else :
                                        nb_transcripts,nb_exon_transcripts,position_list,skip = 0,0,["0"],0
                                        
                                    #if(skip==1): continue
 
                                    array_splicing_data[as_type][key_id_ucsc_event] = { "Ensembl"  :ensembl_id, 
                                                                                    "Symbol"       :gene, 
                                                                                    "Chromosome"   :chrom, 
                                                                                    "Strand"       :strand,
                                                                                    "ic_sample_1"  : ic_sample_1,
                                                                                    "sc_sample_1"  : sc_sample_1,
                                                                                    "ic_sample_2"  : ic_sample_2,
                                                                                    "sc_sample_2"  : sc_sample_2,
                                                                                    "incLevel1"    : incLevel1,
                                                                                    "psiLevel1"    : average_incLevel1,
                                                                                    "psiLevel2"    : average_incLevel2,
                                                                                    "highlight"    : highlight,
                                                                                     "trickMXE"    : trickMXE,
                                                                                    "incLevel2"    : incLevel2,
                                                                                    "diffinc"      : diffinc ,
                                                                                    "logRatioIncLevel"     : logratio,
                                                                                    'fdr'          : fdr,
                                                                                    'log10fdr'     :  int(abs(math.log10(float(fdr)))) if float(fdr) > 0 else -math.inf,
                                                                                    "psi_classifier"     : psi_classifier,
                                                                                    "FCIncLevel"   : retval,
                                                                                    "pvalueFisher" : pvalue_fisher,
                                                                                    "pval"         : "NA",
                                                                                    "fisher"       : fisher,
                                                                                    "id_ucsc"      : id_ucsc_event[:-2],
                                                                                    "gene_size"    : gene_size,
                                                                                    "exon_size"    : exon_size,
                                                                                    "cleanBed"     : highlight[highlight.index(":"):][1:].replace("-","\t"),
                                                                                    "NB_total_transcripts": nb_transcripts,
                                                                                    "NB_exon_transcripts_": nb_exon_transcripts,
                                                                                    "pos_exon_in_transcripts": "::".join(position_list)                                                                       
                            
                                                                                  }
                    #############################################################
                    ###############           CONTROL            ###############
                    ###############                               ###############
                    #############################################################
                    if(isControl):
                        
                        if (ensembl_id not in genes.keys() ) : continue
                        
                        if( ( (as_type != "MXE") and (abs(float(diffinc)) <= 0.01 ) ) or ( (as_type == "MXE") and (abs(float(lineElements[24])) <= 0.01 ) ) ) : 
                            
                            key_id_ucsc_event = chrom+":"+(":".join(lineElements[5:11]))+":"+strand  #flankingES    flankingEE   
                            gonext,id_ucsc_event,highlight,exon_size,trickMXE = parse_rmats(lineElements,chrom,strand,as_type,prefix)

                            if(gonext=="next"): continue

                            # Need to rewrite some stuffs because format is different for MXE !!
                            if (as_type == "MXE" ) :
        
                                key_id_ucsc_event = chrom+":"+(":".join(lineElements[5:13]))+":"+strand  #flankingES    flankingEE    
    
                                ic_sample_1 = lineElements[14]
                                sc_sample_1 = lineElements[15]
                                ic_sample_2 = lineElements[16]
                                sc_sample_2 = lineElements[17]
                                incLevel1   = lineElements[22]
                                incLevel2   = lineElements[23]
                                diffinc     = lineElements[24]
                                pvalue      = lineElements[20]
                                fdr         = lineElements[21]
                            
                            #if (float(pvalue)   > 0.05  and float(fdr)   > 0.05 ) : 
                              
                            test,average_incLevel1, average_incLevel2,logratio,psi_classifier,retval,pvalue_fisher,fisher  = filter_rmats(incLevel1,incLevel2,diffinc,ic_sample_1,sc_sample_1,ic_sample_2,sc_sample_2,readsNumber,isControl)
                            
                            if(average_incLevel1 < 0.9 and average_incLevel2 < 0.9) : continue
                            
                            if(test=="next"): continue
                            
                            ext = "/overlap/region/"+config.parameters['organism']+"/"+chrom.replace("chr","")+":"+highlight[highlight.index(":"):][1:]+"?feature=exon;feature=transcript;content-type=application/json"

                            # Too much time when a lot of exons
                            nb_transcripts,nb_exon_transcripts,position_list,skip = get_position_exons(ext)
                            if(skip==1): continue

                            array_splicing_data[as_type][key_id_ucsc_event] = { "Ensembl"  :ensembl_id, 
                                                                            "Symbol"       :gene, 
                                                                            "Chromosome"   :chrom, 
                                                                            "Strand"       :strand,
                                                                            "ic_sample_1"  : ic_sample_1,
                                                                            "sc_sample_1"  : sc_sample_1,
                                                                            "ic_sample_2"  : ic_sample_2,
                                                                            "sc_sample_2"  : sc_sample_2,
                                                                            "incLevel1"    : incLevel1,
                                                                            "psiLevel1"    : average_incLevel1,
                                                                            "psiLevel2"    : average_incLevel2,
                                                                            "highlight"    : highlight,
                                                                             "trickMXE"    : trickMXE,
                                                                            "incLevel2"    : incLevel2,
                                                                            "diffinc"      : diffinc ,
                                                                            "logRatioIncLevel"     : logratio,
                                                                            'fdr'          : fdr,
                                                                            'log10fdr'     :  int(abs(math.log10(float(fdr)))) if float(fdr) > 0 else -math.inf,
                                                                            "psi_classifier"     : psi_classifier,
                                                                            "FCIncLevel"   : retval,
                                                                            "pvalueFisher" : pvalue_fisher,
                                                                            "pval"         : "NA",
                                                                            "fisher"       : fisher,
                                                                            "id_ucsc"      : id_ucsc_event[:-2],
                                                                            "gene_size"    : gene_size,
                                                                            "exon_size"    : exon_size,
                                                                            "cleanBed"     : highlight[highlight.index(":"):][1:].replace("-","\t"),
                                                                            "NB_total_transcripts": nb_transcripts,
                                                                            "NB_exon_transcripts_": nb_exon_transcripts,
                                                                            "pos_exon_in_transcripts": "::".join(position_list)                                                               
            
                                                                          }                
                                        
                f.close()
                dict_for_analysis.update(array_splicing_data)
    
        return dict_for_analysis
    
                    ############################################################# #############################################################
                    ###############           WHIPPET            ###############
                    ###############                               ###############
                    ############################################################# #############################################################
    
    if (soft=="WHIPPET"):
        for as_type in list_as_type :
            
            logger.info("LOOK SPLICING FOR AS_TYPE: "+as_type)
            array_splicing_data = { as_type :{} }
            namefile = ""
         
            with open(path_to_dir) as f:
                count = 0
                for line in f:
                    #print(line)
                    if count == 0 : 
                        count = count +1
                        continue

                    id_ucsc_event   = ""
                    trickMXE        = ""
                    lineElements    = line.replace("\"", "").strip().split(",")
                    type_event_SE   = lineElements[5]
  
                    if(type_event_SE in ["RI","A5SS","AD","AA","A3SS","SE","CE"]) :
                        ensembl   = lineElements[0]
                       
                        if (isControl) :   
                            if (ensembl not in genes.keys() ) : continue
                            
                        #logger.info(ensembl)  
                        strand           = lineElements[4]
                        ic_sample_1      = "NA"
                        sc_sample_1      = "NA"
                        ic_sample_2      = "NA"
                        sc_sample_2      = "NA"
                        psi_test_mean    = lineElements[6]
                        psi_control_mean = lineElements[7]
                        dpsi             = lineElements[8]
                        fisher          = "-"
                        pvalue_fisher   = '-'
                        highlight       = ""
                        psi_classifier = "NaN" 
                        logratio       = "NaN" 
                        retval         = "NaN"
                        firstelement   = lineElements[3].split(":")
                        chrom          = firstelement[0]
                        section        = lineElements[3].split(":")[1]
                        start_end      = section.split("-")
                        exon_size      = int(start_end[1])-int(start_end[0]) +1 
                        id_ucsc_event  = lineElements[3]+":"+strand  
                        highlight      = prefix+"."+lineElements[3]
                        gene_size      = gene_length[ensembl]

                        key_id_ucsc_event   = lineElements[3]+":"+strand

                        
                        ext = "/overlap/region/"+config.parameters['organism']+"/"+chrom.replace("chr","")+":"+highlight[highlight.index(":"):][1:]+"?feature=exon;feature=transcript;content-type=application/json"
                        nb_transcripts,nb_exon_transcripts,position_list,skip = get_position_exons(ext)
                        if(isControl):
                            if(skip==1): continue

                        array_splicing_data[as_type][key_id_ucsc_event] = { "Ensembl":ensembl, 
                                                                                "Symbol":lineElements[1], 
                                                                                "Chromosome": chrom, 
                                                                                "NB_total_transcripts": nb_transcripts,
                                                                                "NB_exon_transcripts_": nb_exon_transcripts,
                                                                                "pos_exon_in_transcripts": "::".join(position_list),
                                                                                "Strand": strand,
                                                                                "ic_sample_1" : "NA" ,
                                                                                "sc_sample_1" : "NA" ,
                                                                                "ic_sample_2" : "NA" ,
                                                                                "sc_sample_2" : "NA" ,
                                                                                "incLevel1"   : "NA" ,
                                                                                "incLevel2"   : "NA" ,
                                                                                "psiLevel1"   : float(psi_test_mean),
                                                                                "psiLevel2"   : float(psi_control_mean),
                                                                                "highlight"   : highlight,
                                                                                 "trickMXE"   : "NA",
                                                                                "diffinc"     : dpsi,
                                                                                "logRatioIncLevel"   : logratio,
                                                                                'fdr'                : -math.inf,
                                                                                'log10fdr'           :   -math.inf,
                                                                                "psi_classifier"     : psi_classifier,
                                                                                "FCIncLevel"         : retval,
                                                                                "pvalueFisher"       : pvalue_fisher,
                                                                                "pval"      : "NA",
                                                                                "fisher"    : fisher,
                                                                                "id_ucsc"      : id_ucsc_event[:-2],
                                                                                "gene_size" : gene_size,
                                                                                "exon_size" : exon_size,
                                                                                "cleanBed"     : highlight[highlight.index(":"):][1:].replace("-","\t"),
                                                                                 "NB_total_transcripts": nb_transcripts,
                                                                                "NB_exon_transcripts_": nb_exon_transcripts,
                                                                                 "pos_exon_in_transcripts": "::".join(position_list)                       

                                                                              }
                f.close()
                dict_for_analysis.update(array_splicing_data)
                logger.info("    Whippet parsing output done :: ")
    
        return dict_for_analysis

    
def complete_with_expression(path_to_expression_file,dict_for_analysis):
    """
    Add to dic with splicing data, the gene expression Fold Change
  
    Args:
        path_to_expression_file (str): Path to the file generated with DESEQ2.
        dict_for_analysis (dic): dict object to upgrade

    Returns:
        splicing_dict (obj): Return all the informations in a structured object.(gene FC added)
    
        
    """
    geneExpression = {}
    if (os.path.isfile(path_to_expression_file)) : 
        with open(path_to_expression_file) as f:
            logger.info("LOOK EXPRESSION FOR: "+path_to_expression_file)
    
            count = 1
            for line in f:
                #print(line)
                if count == 1 : 
                    count = 2
                    continue
                id_ucsc_event = ""
                lineElements  = line.replace("\"", "").strip().split(",")
                if (lineElements[12] != "NA") : 
                    if (float(lineElements[12]) > 0.05) : 
                        lineElements[12] = "NaN" 
                geneExpression[lineElements[0]] = {"Symbol_biomart":lineElements[1],"gene_biotype":lineElements[2],"log2fc" :lineElements[8],"fc":lineElements[13],"padj": ( "NaN" if (lineElements[12] == "NA") else lineElements[12].replace("e","E"))} 
                #print(line)
        f.close()
    #pp = pprint.PrettyPrinter()

    for event in dict_for_analysis:
        
        logger.info("UPGRADE DICTIONNARY FOR EVENT : "+event)

        for id_ucsc_event in dict_for_analysis[event].keys():

                #print("id_ucsc_event : "+id_ucsc_event)

                if  dict_for_analysis[event][id_ucsc_event]["Ensembl"] in geneExpression :                    
                    dict_for_analysis[event][id_ucsc_event].update(geneExpression[dict_for_analysis[event][id_ucsc_event]["Ensembl"]])
                else : 
                    dict_for_analysis[event][id_ucsc_event].update( {"gene_biotype":"-","log2fc" :"NaN","fc":"NaN","padj":"NaN"} )

            #pp.pprint(analysis[event][id_ucsc_event])

    return dict_for_analysis

def complete_with_quantif(list_all_replicates,dict_for_analysis):
    """
    Add to dictionary  TPM values fro each replicates in the considered analysis
  
    Args:
        list_replicates (list): List of dictionaries with id, path etc....
        dict_for_analysis (dic): dict object to upgrade
        tpm_cutoff (int): tpm cutoff
    Returns:
        splicing_dict (obj): Return all the informations in a structured object.(TPM added)
    
        
    """
    quantification            = {}
    #list_of_replicates       = []
    tpm_cutoff = 0
    #for replicate_object in list_all_replicates:
        #list_of_replicates.append([replicate_object["replicat_id"],replicate_object["file_path"]])
    cutOffDict = {}
    for replicat in list_all_replicates:
        #print(replicat)
        logger.info("LOOK QUANTIFICATION IN : "+replicat["replicat_id"] )
        tpms = []
        numReads = []
        
        if (os.path.isfile(replicat["file_path"])) : 
            with open(replicat["file_path"]) as f:
                count = -1
                for line in f:
                   
                    if count == -1 : 
                        count = 0
                        continue
        
                    lineElements = line.strip().split("\t")
                    
                    m            = re.search('^(\w+)\.(\d+)$', lineElements[0])
                    #NB : ENST00000638165.1|ENSG00000147862.15|OTTHUMG00000021027.6|OTTHUMT00000488972.1|NFIB-018|NFIB|1783|processed_transcript|    1783    1619.23    0.0452832    7.15216
                    if(m):
                        ensembl_id = m.group(1)
                        if ensembl_id not in quantification : quantification[ensembl_id] = {}
                        quantification[ensembl_id].update( { replicat["replicat_id"]  : {"TPM":lineElements[3],"NumReads" :lineElements[4]} })
                        if(float(lineElements[3]) > 0 ) : 
                            tpms.append( float(lineElements[3]))
                            numReads.append( float(lineElements[4]))
                    else : 
                        if ensembl_id not in quantification : quantification[ensembl_id] = {}
                        quantification[ensembl_id].update( {  replicat["replicat_id"]  : {"TPM":"NaN","NumReads" :"NaN"} })
    
                    #print(line)
            f.close()
        
        # USE DEFINED USER CUT OFF SET FOR ALL SAMPLES . DEFAULT 0
        cutOffDict[replicat["replicat_id"]]=tpm_cutoff

    #pp = pprint.PrettyPrinter()
    
    for event in dict_for_analysis.keys():
     
            logger.info("UPGRADE DICTIONARY IN "+event)

            for id_ucsc_event in dict_for_analysis[event].keys():
        
                #print("id_ucsc_event : "+id_ucsc_event)

                if  dict_for_analysis[event][id_ucsc_event]["Ensembl"] in quantification :
                    
                    dict_for_analysis[event][id_ucsc_event].update(quantification[dict_for_analysis[event][id_ucsc_event]["Ensembl"]])
                else : 
                    logger.info("Error in complete_with_quantif")

            #pp.pprint(dict_for_analysis[event][id_ucsc_event])

    return dict_for_analysis,cutOffDict

def return_all_uniq_replicates_object(dict_samples) :
    """
    Return list of unique name of replicates from your json config
  
    Args:
        dict_samples(dict(list(dict))): Dict of samples in json config

    Returns:
       all_replicates_object (list): Return list of replicates names
   
    """
    dict_rep = {}
    for sample in sorted(dict_samples.keys()) :
        # ordered because it's a
        for replicat in dict_samples[sample]  :
            dict_rep[ replicat["replicat_id"]] = 0
    
    all_replicates_object   = []
    
    for sample in sorted(dict_samples.keys()) :
        # list so it's sorted
        for replicat in  dict_samples[sample] :
            
            if dict_rep[replicat["replicat_id"]] == 0 : 
                dict_rep[replicat["replicat_id"]] =  1
                all_replicates_object.append(replicat) 
            
    return all_replicates_object

def replicates_per_condition(dict_samples,condition) :
    """
    Return list of unique name of replicates from your json config
  
    Args:
        dict_samples(dict(list(dict))): Dict of samples in json config

    Returns:
       all_replicates_object (list): Return list of replicates names
   
    """
    reps = []
    for sample in sorted(dict_samples.keys()) :
        if sample == condition :
            # ordered because it's a
            for replicat in dict_samples[sample]  :
                reps.append(replicat["replicat_id"])

         
    return reps

def remove_values_from_list(the_list, val):
    """
    Remove specific value from list
  
    Args:
        the_list(list): List of values
        val(str|int): Value to remove

    Returns:
        list (list): Return list
   
    """
    return [value for value in the_list if value != val]

def median(lst):
    """
    Compute median from a list
  
    Args:
        list(list): List of values
   
    Returns:
        meian (float): Return median float
   
    """
    return numpy.median(numpy.array(lst))

def _convert_to_number(cell):
    if cell.isnumeric():
        return int(cell)
    try:
        return float(cell)
    except ValueError:
        return cell
    

            
            
def create_header(dict_samples, list_analysis_to_check):
    """
    Create header of the file
  
    Args:
        dict_samples(dict(list(dict))): Dict of condition in json config
        list_analysis_to_check:  List of analysis/condition you want in tab

    Returns:
        headerListOfFields (list): Return list of fields for the header
   
    """
    conditions = list_analysis_to_check[0].split("_vs_")
    logger.info("  conditions  :: "+" - ".join(conditions))

    headerListOfFields                                    = []
    list_all_TPM_replicates_in_samples_for_quantification = []
    header_variabe = ["Coordinates","Test:inclusion.reads|exclusion.reads","Control:inclusion.reads|exclusion.reads","Test.PSI|Control.PSI","Test.Average.PSI","Control.Average.PSI","dPSI","Fisher","LOG10-PSI-FDR","PSI-FoldChange","Gene-Log2FoldChange","Gene-Padj","Gene.Raw.Reads","Test:Gene.Raw.Reads","Control:Gene.Raw.Reads","TotalRawReadsNormalisedPerGeneSize"] 
    coreHeaderFields =  ['ID-UCSC','Epissage','Event','Symbol','Ensembl','Track Bed Style','Strand','Gene_biotype','Gene_size','Exon_size',"TranscriptsWithExon","AllTranscripts","Pos_exon_in_transcripts"]               
                                                                                  
    
    #for cond in sorted(dict_samples.keys()) :
        # ordered because it's a
        # if cond in conditions :

            #for replicat in dict_samples[cond]  :
                #list_all_TPM_replicates_in_samples_for_quantification.append(replicat["replicat_id"]+"-TPM")
                #list_all_TPM_replicates_in_samples_for_quantification.append(replicat["replicat_id"]+"-NumReads")

    #keep unique only
    # tpm_header  = list(OrderedDict.fromkeys(list_all_TPM_replicates_in_samples_for_quantification))
    tpm_header= []
    tpm_header.append("Test:Gene.Estimated.Reads")
    tpm_header.append("Test:Gene.TPMs")

    tpm_header.append("Control:Gene.Estimated.Reads")
    tpm_header.append("Control:Gene.TPMs")
    
    tpm_header.append("Test:Gene.Mean.Estimated.Reads")
    tpm_header.append("Test:Gene.Mean.TPM")
    
    tpm_header.append("Control:Gene.Mean.Estimated.Reads")
    tpm_header.append("Control:Gene.Mean.TPM")

    
    #logger.info("  TPM_HEADER  :: "+" - ".join(tpm_header))

    for condition in (list_analysis_to_check) :
            
                for x in range(0,len(header_variabe)) :
                    
                    headerListOfFields.append(condition+"-"+header_variabe[x])

    
    headerListOfFields = sum([coreHeaderFields,headerListOfFields,tpm_header],[]) 
    

    return headerListOfFields,header_variabe

   
def return_formated(n):
   
    return ("NaN" if  n=="NaN" else ("%0.3f" % float(n)))    #.replace('.',',')
    """
    Compute fold change from two values num and denum
  
    Args:
        num (float|int): numerator of ratio
        denum (float|int):  denumerator of ratio

    Returns:
        ratio (float): Return fold change
   
    """
def foldchange (num,denom) :
  
    if (num >= denom ) : 
        return num/denom
    else :
        return -denom/num

    """
    Compute fold change from logratio
  
    Args:
        logratio (float|int): logratio   

    Returns:
        fc (float): Return fold change
   
    """  
def logratio2foldchange (logratio, base) :

    retval = base^(logratio)
    if retval < 1 :
        return -1/retval
    else :
        return  retval

    """
    Compute logratio from fold change  
  
    Args:
        fold change (float|int): fold change 

    Returns:
              (float): Return logratio
   
    """  
def foldchange2logratio (foldchange) :

    if foldchange < 0 :
        retval = 1/-foldchange
    else : 
        retval =  foldchange 
        
    retval = math.log2(retval)
    return retval


def check_distribution(list_values) :
    
        a = np.array(list_values)
        cutOffPercentile = np.percentile(a, 50) # return 50th percentile, e.g median.
        logger.info("GENES : "+str(len(list_values)))
        logger.info("    MIN - MAX : ["+str(np.amin(a))+" - "+str(np.amax(a))+"]")   
        logger.info("    MEAN : "+str(np.mean(a)))   
        logger.info("    MEDIAN : "+str(np.median(a)))   
        
        return cutOffPercentile

def blacklist_gene_under_median_at_least_in_one_sample(matrice_path,names_test_replicates,names_control_replicates,gene_length):
    """
    Return a list of ensembl genes where distribution of their RPKM value is under median of its distribution in at leat one replicate
    Args:
        names_test_replicates (list): List of column to parse in matrice for test....
        names_control_replicates (list): List of column to parse in matrice for control....
        matrice_path (str) : path to file with raw reads count from STAR
        gene_length (dic) : Hash gene_ensembl to size all exons
    Returns:
        blacklist_gene (array): Return all the genes name that sould be skipped for splicing analysis
    
    """
    
    logger.info("BlackList Genes from Raw Count : "+matrice_path )
    data = pd.read_csv(matrice_path, sep=",")
    indexed_data = data.set_index(['Gene'])
    
    all_replicates = list(itertools.chain(names_test_replicates,names_control_replicates))
    rpkm           = {}
    rpkm_gene_dic      = {}
    blacklist_gene     = []
     
    for replicate in all_replicates :
        logger.info("Replicate-id from DESEQ2 FC matrice : "+replicate)
        total_mapped_reads = indexed_data[replicate].sum()
        scaling_factor = total_mapped_reads / 1000000
        logger.info("Total Uniquely Mapped Reads "+str(total_mapped_reads))

        rpkm[replicate] = []
        rpkm_gene_dic[replicate] = []
        for gene in  indexed_data.index.values :
            m            = re.search('^(\w+)\.(\d+)$', gene) # Ok you miss somes of them ENSEMBL_PARY that you count in the sum but that nothing compared to the whole thing
            if(m):
                ensembl_id     = m.group(1)    
                gene_length_kb = gene_length[ensembl_id]/1000
                if(indexed_data.get_value(gene,replicate) > 0 ) :
                    rpkm_gene      = (indexed_data.get_value(gene,replicate)/scaling_factor)/gene_length_kb 
                    rpkm[replicate].append(rpkm_gene)
                    rpkm_gene_dic[replicate].append(ensembl_id)
                if(indexed_data.get_value(gene,replicate) == 0 ) :    
                    blacklist_gene.append(ensembl_id)
   

    for rep in all_replicates :
            logger.info("Replicate: "+rep)

            cut_off = check_distribution(rpkm[rep])
            
            for index,rpkm_value in enumerate(rpkm[rep]):

                if (rpkm_value < cut_off ) : 
                    blacklist_gene.append(rpkm_gene_dic[rep][index])
                    
                    #logger.info(rpkm_gene_dic[rep][index])
                    #logger.info(rpkm_value)
    logger.info("")
                
    return list(set(blacklist_gene))
                    
def complete_with_raw_read_count(matrice_path,dict_for_analysis,names_test_replicates,names_control_replicates):
    """
    Add to dictionary  raw count given by star, this will be use for expression filtering
  
    Args:
        names_test_replicates (list): List of column to parse in matrice for test....
        names_control_replicates (list): List of column to parse in matrice for control....
        dict_for_analysis (dic): dict object to upgrade
    Returns:
        splicing_dict (obj): Return all the informations in a structured object.(gene FC added)
    
    """
    raw_read_count = {}

    logger.info("Check Raw counts : "+matrice_path )
    data = pd.read_csv(matrice_path, sep=",")
    indexed_data = data.set_index(['Gene'])

    for gene in  indexed_data.index.values :

        m            = re.search('^(\w+)\.(\d+)$', gene)
                        #NB : ENST00000638165.1|ENSG00000147862.15|OTTHUMG00000021027.6|OTTHUMT00000488972.1|NFIB-018|NFIB|1783|processed_transcript|    1783    1619.23    0.0452832    7.15216
        if(m):
            ensembl_id = m.group(1)
            if ensembl_id not in raw_read_count : raw_read_count[ensembl_id] = {}
            for control in names_control_replicates :
                raw_read_count[ensembl_id].update( {  control+"_rc"  : {"rawReads":indexed_data.get_value(gene,control)} })
            for test in names_test_replicates :
                raw_read_count[ensembl_id].update( { test+"_rc"  : {"rawReads":indexed_data.get_value(gene,test)} } )
        else : 
            if ensembl_id not in raw_read_count : raw_read_count[ensembl_id] = {}
            for control in names_control_replicates :
                raw_read_count[ensembl_id].update( { control+"_rc"  : {"rawReads":"NaN"} })
            for test in names_test_replicates :
                raw_read_count[ensembl_id].update( { test+"_rc"  : {"rawReads":"NaN"} })            


    #pp = pprint.PrettyPrinter()
    for event in dict_for_analysis.keys():
     
            logger.info("UPGRADE DICTIONARY IN "+event)

            for id_ucsc_event in dict_for_analysis[event].keys():
        

                if  dict_for_analysis[event][id_ucsc_event]["Ensembl"] in raw_read_count :
                    
                    dict_for_analysis[event][id_ucsc_event].update(raw_read_count[dict_for_analysis[event][id_ucsc_event]["Ensembl"]])
                else : 
                    logger.info("Error in complete_with_raw_read_count")

            #pp.pprint(dict_for_analysis[event][id_ucsc_event])

    return dict_for_analysis

def createBackgroundGeneList(matrice_path,names_test_replicates,names_control_replicates,blacklist_gene,output,analyse_name):
    """
    Return a list of ensembl genes for Go Background purpose
    Args:
        names_test_replicates (list): List of column to parse in matrice for test....
        names_control_replicates (list): List of column to parse in matrice for control....
        matrice_path (str) : path to file with raw reads count from STAR
        blacklist_gene (array) : Blacklisted genes
    Returns:
        null
    """
    
    logger.info("Use : "+matrice_path )
    data = pd.read_csv(matrice_path, sep=",")
    indexed_data = data.set_index(['Gene'])
    black_list = []
    for gene in  indexed_data.index.values :
        m            = re.search('^(\w+)\.(\d+)$', gene) #
        if(m):
            ensembl_id     = m.group(1)    
            if(ensembl_id not in blacklist_gene) :
                black_list.append(ensembl_id+"\n")

   
    black_file = open(output+analyse_name+'.txt', 'w')

    for  line in  list(set( black_list ))   :
        black_file.write(line)
    black_file.close()

def complete_with_gene_length(gene_path):
    """
    Add to dictionary gene count
  
    Args:
        gene_path (string): path to gene length file
    Returns:
        dict_length (obj): Return hash nameGene=> Length
    
    """
    gene_length = {}
    
    with open(gene_path) as f:
        logger.info("    File used for gene sizes : "+gene_path)

        for line in f:

            lineElements  = line.strip().split(",")
            m            = re.search('^(\w+)\.(\d+)$', lineElements[0])
            if(m):
                ensembl_id = m.group(1)
                gene_length[ensembl_id]= int(lineElements[1])
           
    f.close()
   

    return gene_length

def sublist(path):
    """
    Add to dictionary gene count
  
    Args:
        gene_path (string): path to gene  file
    Returns:
        dict_length (obj): Return hash ensembl=> ensembl
    
    """
    gene = {}
    logger.info("Gene list to test is here : : "+path)

    with open(path) as f:

        for line in f:

            lineElements  = line.strip().split("\t")
            gene[lineElements[0]]= lineElements[0]
           
    f.close()

    return gene

def rewriteBed(output_crosslink,rewritedBed) :
    
    newfile = open(rewritedBed,"w")
    with open(output_crosslink) as f:

        for line in f:

            lineElements  = line.strip().split("\t")
            
            if(lineElements[7]=="-1") :
                logger.info ("==========> EXON WAS NOT FOUND WTF !!!!")
                continue
            
            chr_init = lineElements[0]
            start_init = lineElements[1]
            end_init = lineElements[2]
           
            name_init = lineElements[3]
            
            match_init = re.search('^(.*)_(.*)_(.*)_(.*)_(.*)$',name_init)        
            
            dpsi_init     = match_init.group(1)
            symbol_init   = match_init.group(2)
            ensembl_init  = match_init.group(3)
            psi_init      = match_init.group(4)
            pos_init      = match_init.group(5)

            signal_init = lineElements[4]
            strand_init = lineElements[5]
            
            name_new = lineElements[9]
            match_new = re.search('^(.*)_(.*)$',name_new)        

            ensembl_new  = match_new.group(1)
            pos_new      = match_new.group(2)
            
            if(ensembl_new != ensembl_init) :
                continue
            
            name_rewrited = "_".join([dpsi_init,symbol_init,ensembl_init,psi_init,pos_new])
            #

            newfile.write("\t".join([chr_init,start_init,end_init,name_rewrited,signal_init,strand_init])+"\n")
           
    f.close()
    newfile.close()
        #chrX    100627108    100629986    -0.708_CD44_ENSG00000026508_0.96_3    0    -    chrX    100627108    100629986    ENSG00000000003_10    0    -
        #chr11    35209964    35210054    -0.708_CD44_ENSG00000026508_0.96_3    -0.708    +
        #chrX    100007108    100009986    ENSG00000000003_10    0    -    .    -1    -1    .    -1    .
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will merge all intermediate outputs.  
    Example : 
    python3 ./mergeFinal.py -c /home/jp/Desktop/MANTrep1.json 

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')
    parser.add_argument("-mc","--modeCount",action="store", default="2" ,help="with or without ReadsOnTarget",required=False,type=str,dest='modeCount')
    parser.add_argument("-r","--readsJunction",action="store",help="Number of reads on junction used to filter out events.",required=False,default=10,type=int,dest='readsNumber')
    parser.add_argument("-e","--event",action="store",help="Type of events",required=True,type=str,dest='event')
    parser.add_argument("-ct","--control",action="store_true",help="AS or Control",default=False,dest='iscontrol')
    parser.add_argument("-i","--init",action="store",help="Path to a json file.",required=False,default=str(Path(os.path.dirname(os.path.realpath(__file__))).parents[0])+"/config/init.json",type=str,dest='init')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")
    init = custom_parser.Configuration(parameters.init,"json")


    logger = create_logger(config,config.parameters['project'])

    logger.info('\n =========> START MERGING AND CLEANING : ')
    
    logger.info("Iscontrol ? : "+str(parameters.iscontrol))
    logger.info("file_config : "+parameters.file_config)
    logger.info("path_to_output : "+config.parameters['path_to_output'])
    logger.info("organism : "+config.parameters['organism'])
    logger.info("event : "+parameters.event)
    logger.info("soft_version_rmats : "+config.parameters['soft_version_rmats'])
    logger.info("postitionGenomicExon : "+config.parameters['postitionGenomicExon'])

    events = []
    events.append(parameters.event)

    sessionName=config.parameters['sessionName']
    logger.info("sessionName :  "+sessionName)

    namefile = ""
    if (parameters.modeCount == "2" ) :
        if(config.parameters["soft_version_rmats"]=="old"):
            namefile = "MATS.JunctionCountOnly.txt"
        if(config.parameters["soft_version_rmats"]=="new"):
            namefile = "MATS.JC.txt"
    else  :
        if(config.parameters["soft_version_rmats"]=="new"):
            namefile = "MATS.JCEC.txt"
        if(config.parameters["soft_version_rmats"]=="old"):
            namefile = "MATS.ReadsOnTargetAndJunctionCounts.txt"


    logger.info ("########################################################################################################################")
    logger.info (" DataFrame Object Construction : First Part)")
    logger.info ("########################################################################################################################")
   
    catalog         = {}
    cutOffDict      = {}
    gene_length=complete_with_gene_length(config.parameters['gene_length'])
    
    mode="clean"
    if(parameters.iscontrol) :
        mode="bad"
        
    logger.info("    mode : "+mode)

 
    for analyse_name_loop in (config.parameters.get("analysis").keys() ):
            

        genes = {}
        if(parameters.iscontrol) :
            genes = sublist(config.parameters['path_to_output']+"/"+analyse_name_loop+"/"+parameters.event+"/"+analyse_name_loop+"."+'gene.ensembl.txt')

        catalog[analyse_name_loop] = {}
    
        logger.info("    Analyse_name : "+analyse_name_loop)
        
        if (config.parameters["analysis"][analyse_name_loop]["soft"] == "WHIPPET") :
            if (parameters.event == "SE") :
                config.parameters["analysis"][analyse_name_loop]["splicing_dir"] = config.parameters["analysis"][analyse_name_loop]["splicing_dir"]+analyse_name_loop+"."+mode+".CE.diffannoted.csv"
            if (parameters.event == "A5SS") :
                config.parameters["analysis"][analyse_name_loop]["splicing_dir"] = config.parameters["analysis"][analyse_name_loop]["splicing_dir"]+analyse_name_loop+"."+mode+".AD.diffannoted.csv"
            if (parameters.event == "A3SS") :
                config.parameters["analysis"][analyse_name_loop]["splicing_dir"] = config.parameters["analysis"][analyse_name_loop]["splicing_dir"]+analyse_name_loop+"."+mode+".AA.diffannoted.csv"
            if (parameters.event == "RI") :
                config.parameters["analysis"][analyse_name_loop]["splicing_dir"] = config.parameters["analysis"][analyse_name_loop]["splicing_dir"]+analyse_name_loop+"."+mode+".RI.diffannoted.csv"
        
        logger.info ("    splicing : "+config.parameters["analysis"][analyse_name_loop]["splicing_dir"])
        logger.info ("    expression : "+config.parameters["analysis"][analyse_name_loop]["expression_file_path"])
        logger.info ("    sample_control_for_quantification : "+config.parameters["analysis"][analyse_name_loop]["sample_control_for_quantification"])
        logger.info ("    sample_test_for_quantification : "+config.parameters["analysis"][analyse_name_loop]["sample_test_for_quantification"])
        logger.info ("    list_splicing_files : "+" ".join(events))
        logger.info ("")
        
        # If there is an  inversion in sample_control_for_quantification & sample_test_for_quantification no impact...don't worry Here it's just print stuffs
        #########################

        test = config.parameters["analysis"][analyse_name_loop]["sample_test_for_quantification"]

        for replicat_id in config.parameters["samples_for_quantification"][test] :
            logger.info ("    Replicat_id - test : "+replicat_id.get("replicat_id"))
        logger.info ("")
        
        control = config.parameters["analysis"][analyse_name_loop]["sample_control_for_quantification"]

        for replicat_id in config.parameters["samples_for_quantification"][control] :
            logger.info ("    Replicat_id - control : "+replicat_id.get("replicat_id"))
       
        logger.info ("GO...")
        logger.info ("")
        logger.info ("Soft : "+config.parameters["analysis"][analyse_name_loop]["soft"])
        #########################

        catalog[analyse_name_loop] = parse_all_splicing_files(config.parameters["analysis"][analyse_name_loop]["splicing_dir"],events, catalog[analyse_name_loop],parameters.readsNumber,gene_length,config.parameters["analysis"][analyse_name_loop]["soft"],config.parameters["analysis"][analyse_name_loop]["replicates_test"],config.parameters["analysis"][analyse_name_loop]["replicates_control"],config.parameters["organism"],genes,parameters.iscontrol)
        
        logger.info ("")

        catalog[analyse_name_loop] = complete_with_expression(config.parameters["analysis"][analyse_name_loop]["expression_file_path"], catalog[analyse_name_loop])
        
        logger.info ("")
        
        names_test_replicates_for_raw_count       =  config.parameters["analysis"][analyse_name_loop]["replicates_test"]
        names_control_replicates_for_raw_count    =  config.parameters["analysis"][analyse_name_loop]["replicates_control"]
       
        replicates    =  list(itertools.chain(names_test_replicates_for_raw_count,names_control_replicates_for_raw_count))
        logger.info ("REPLICATES FOR "+analyse_name_loop)
        logger.info ("CHECK FOR BLACKLISTED GENES FROM (low expression under median distribution) "+analyse_name_loop)

        genes_blacklisted = blacklist_gene_under_median_at_least_in_one_sample(config.parameters["read_count_matrice"],names_test_replicates_for_raw_count,names_control_replicates_for_raw_count,gene_length)

        logger.info ("CREATE BACKGOUND FOR FUTURE GO ANALYSES "+analyse_name_loop)

        #createBackgroundGeneList(config.parameters["read_count_matrice"],names_test_replicates_for_raw_count,names_control_replicates_for_raw_count,genes_blacklisted,config.parameters['path_to_output']+parameters.event+"/",analyse_name_loop)

        logger.info (replicates)
        
        catalog[analyse_name_loop] = complete_with_raw_read_count(config.parameters["read_count_matrice"], catalog[analyse_name_loop],names_test_replicates_for_raw_count,names_control_replicates_for_raw_count)
        
        logger.info ("")
        
        # If inversion in sample_control_for_quantification & sample_test_for_quantification no impact...don't worry
        #control_replicates   =  config.parameters["samples_for_quantification"][config.parameters["analysis"][analyse_name_loop]["sample_control_for_quantification"]]
        #test_replicates      =  config.parameters["samples_for_quantification"][config.parameters["analysis"][analyse_name_loop]["sample_test_for_quantification"]]
  

        catalog[analyse_name_loop],cutOffDict = complete_with_quantif(return_all_uniq_replicates_object(config.parameters["samples_for_quantification"]), catalog[analyse_name_loop])
        #print(cutOffDict)
        logger.info ("")
        logger.info ("Next...")
        logger.info ("")
   
    logger.info ("Finish...")
    

    logger.info ("########################################################################################################################")
    logger.info ("That the tricky part to structure data to fit your envy (biologist want an evil format : excel...)")
    logger.info ("########################################################################################################################")
    namefile=namefile.replace(".txt","")
    namefile=namefile.replace(".MATS","")
    
   # 2018-06-20 16:05:30,326 :: INFO :: Excel : /home/jean-philippe.villemin/data/data/PROJECT/WHOLE_EMT//paired.HMLE.TAMOXIFEN.SNAIL.TEST.0.5_vs_paired.HMLE.TAMOXIFEN.SNAIL.CONTROL.0/RI/RI.paired.HMLE.TAMOXIFEN.SNAIL.TEST.0.5_vs_paired.HMLE.TAMOXIFEN.SNAIL.CONTROL.0.clean.SPLICING.xlsx

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    file = config.parameters['path_to_output']+"/"+analyse_name_loop+"/"+parameters.event+"/"+parameters.event+"."+analyse_name_loop+"."+mode+'.SPLICING.xlsx'
    logger.info ("Excel : "+file)

    writer = pd.ExcelWriter(file, engine='xlsxwriter')
    logger.info (namefile)

    dict_analysis               = config.parameters.get("analysis")
    dict_samples                = config.parameters["samples_for_quantification"]
    list_analysis               = list(sorted(config.parameters.get("analysis").keys()))

    
    ''' That the tricky part to structure data to fit your envy'''
    #Foreach tab you want to create
    for tab in config.parameters.get("tabs").keys():
    
        subprocess.run(("mkdir -p "+config.parameters['path_to_output']+"/"+tab+"/"+parameters.event),shell=True)

        logger.info("")
        logger.info("  TAB  :: "+tab)
        logger.info("  analysisTotcheck  :: "+" - "+config.parameters["tabs"][tab]["analysisTocheck"][0])
        
        conditions = []
        conditions = config.parameters["tabs"][tab]["analysisTocheck"][0].split("_vs_")
        logger.info("    Cond1 : "+conditions[0])
        logger.info("    Cond2 : "+conditions[1])

        lines_for_my_tab = []
        list_analysis_with_event    = set(config.parameters["tabs"][tab]["analysisTocheck"])
        
        list_analysisNotTotcheck = [ x for x in list_analysis if x not in list_analysis_with_event ]
        logger.info("  Analysis Not To check  :: "+" - ".join(list_analysisNotTotcheck))
        logger.info("  GO...")
        
        replicates_filter = replicates_per_condition(config.parameters["samples_for_quantification"],conditions[0]) + replicates_per_condition(config.parameters["samples_for_quantification"],conditions[1])
        logger.info("  Replicates used for TPM  :: "+" - ".join(replicates_filter))

        headerListOfFields,header_variabe = create_header(dict_samples,config.parameters["tabs"][tab]["analysisTocheck"])

        length_fields_sometimes_empty = len(header_variabe)
        
        
        #print("HEADER COLUMNS : "+"\n" +"\n".join(headerListOfFields)) 
        # You are going to check the first analysis in list vs all others for you tab definition 
        for analyse_name in config.parameters["tabs"][tab]["analysisTocheck"]:
          
            logger.info("  ") 
            logger.info("    ANALYSE TO CHECK :: "+analyse_name)

            # By event
          
            for event in events :

                logger.info("          EVENT  :: "+event)
                logger.info("          NB EVENTS  :: "+str(len(catalog[analyse_name][event])))
                logger.info("     Run trought the list of keys coordinates in current analysis...")
                counter = 0
                
                
                for id_ucsc in catalog[analyse_name][event]:
                   
                    counter=counter + 1

                    if (float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > 0 ) :  epissage = 'INC' 
                    else : 
                        epissage = 'EXCL'
   
                    #  Check For Presence in 
                    analysis_to_ucscKey = {}
                    analysis_to_ucscKey[analyse_name] = id_ucsc

                    # USE BLACK LIST GENE UNDER MEDIAN  BASED ON EXPRESSION
                    if(catalog[analyse_name][event][id_ucsc]["Ensembl"] in genes_blacklisted) : continue

                    index = 0

                    for another_analyse1 in config.parameters["tabs"][tab]["analysisTocheck"]:
                        
                        if (analyse_name == another_analyse1   ) : continue

                        #logger.info(another_analyse1)
                        index =  config.parameters["tabs"][tab]["analysisTocheck"].index(another_analyse1)

                        if id_ucsc in catalog[another_analyse1][event] :
                           
                            if (math.isnan(float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]))== True or math.isnan(float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]))==True) :
                                analysis_to_ucscKey[another_analyse1] = id_ucsc
                                break
                        

                    '''
                    We find the event in the others analysis
                    '''
                    if (len(analysis_to_ucscKey) == len(config.parameters["tabs"][tab]["analysisTocheck"]) ) :            
                       
                
                        features         = []
                        chro             = catalog[analyse_name][event][id_ucsc]["Chromosome"]
                        id_ucsc_clean    = catalog[analyse_name][event][id_ucsc]["id_ucsc"]
                        cleanBed         = catalog[analyse_name][event][id_ucsc]["cleanBed"]
                        
                        if (event == 'MXE') : cleanBed = catalog[analyse_name][event][id_ucsc]["trickMXE"]
                        
                        features.extend([ 
                                         "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jp&hgS_otherUserSessionName="+sessionName+"&position="+id_ucsc_clean+"&highlight="+catalog[analyse_name][event][id_ucsc]["highlight"],
                                         epissage,
                                         event,
                                         catalog[analyse_name][event][id_ucsc]["Symbol"],
                                         catalog[analyse_name][event][id_ucsc]["Ensembl"],
                                         chro+"\t"+cleanBed+"\t"+str(catalog[analyse_name][event][id_ucsc]["diffinc"])+"_"+catalog[analyse_name][event][id_ucsc]["Symbol"]+"\t"+str(catalog[analyse_name][event][id_ucsc]["diffinc"])+"\t"+catalog[analyse_name][event][id_ucsc]["Strand"] ,
                                         catalog[analyse_name][event][id_ucsc]["Strand"],
                                         catalog[analyse_name][event][id_ucsc]["gene_biotype"],
                                         catalog[analyse_name][event][id_ucsc]["gene_size"],
                                         catalog[analyse_name][event][id_ucsc]["exon_size"],
                                         str(catalog[analyse_name][event][id_ucsc]["NB_exon_transcripts_"]),
                                         str(catalog[analyse_name][event][id_ucsc]["NB_total_transcripts"]),
                                         
                                         catalog[analyse_name][event][id_ucsc]["pos_exon_in_transcripts"]


                                        ]
                                        )
                        # Add Gnene FC , Inclusion Level in each other analysis present in the lab
                        for analyse_name_bis in config.parameters["tabs"][tab]["analysisTocheck"] :
                                
                                if analyse_name_bis in analysis_to_ucscKey :
                                    
                                    if ( config.parameters["tabs"][tab][str(index)] != "notIn" ) :
                                        raw_read = []
                                        raw_read_test = []
                                        raw_read_control = []

                                        names_test_replicates_for_raw_count       =  config.parameters["analysis"][analyse_name]["replicates_test"]
                                        names_control_replicates_for_raw_count    =  config.parameters["analysis"][analyse_name]["replicates_control"]
                                       
                                        myreplicates    =  list(itertools.chain(names_test_replicates_for_raw_count,names_control_replicates_for_raw_count))

                                        for rep in myreplicates :
                                            #logger.info(catalog[analyse_name][event][analysis_to_ucscKey[analyse_name_bis]])
                                            raw_read.append(catalog[analyse_name][event][analysis_to_ucscKey[analyse_name]][rep+"_rc"]["rawReads"])
                                            if( rep in names_test_replicates_for_raw_count) :
                                                raw_read_test.append(catalog[analyse_name][event][analysis_to_ucscKey[analyse_name]][rep+"_rc"]["rawReads"])
                                            if( rep in names_control_replicates_for_raw_count) :
                                                raw_read_control.append(catalog[analyse_name][event][analysis_to_ucscKey[analyse_name]][rep+"_rc"]["rawReads"])

                                        features.extend([
                                                          analysis_to_ucscKey[analyse_name_bis],

                                                          str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["ic_sample_1"]).replace(",","::")+"|"+str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["sc_sample_1"]).replace(",","::")
                                                         
                                                          ,str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["ic_sample_2"]).replace(",","::")+"|"+str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["sc_sample_2"]).replace(",","::")
                                                           ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["incLevel1"].replace(",","::")+"|"+catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["incLevel2"].replace(",","::")

                                                          #,"::".join(str(x) for x in catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["incLevel1"])+"|"+"::".join(str(x) for x in catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["incLevel2"])
                                                          ,"{0:.2f}".format(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["psiLevel1"])
                                                          ,"{0:.2f}".format(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["psiLevel2"])
                                                          ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["diffinc"]
                                                          ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["fisher"] 
                                                          ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["log10fdr"]
                                                          #  ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["psi_classifier"]  
                                                          #,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["logRatioIncLevel"])
                                                          ,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["FCIncLevel"])
                                                          #,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["pval"].replace("e","E")                                           
                                                            ,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["log2fc"] )
                                                          #,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["fc"] )
                                                           ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["padj"]
                                                            ,"::".join(str(x) for x in raw_read)
                                                            ,"{0:.2f}".format(statistics.mean(raw_read_test))
                                                            ,"{0:.2f}".format(statistics.mean(raw_read_control))
                                                           ,(sum(raw_read)/catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["gene_size"])
                                                          

                                                         ])
                                    else :   
                                        for i in range(0,length_fields_sometimes_empty) : features.extend(["."])
    
    
                        tpm_test  = []
                        read_test = []

                        for rep_test in replicates_per_condition(config.parameters["samples_for_quantification"],conditions[0]) :  #config.parameters["samples_for_quantification"][sample]
                                tpm_test.append(catalog[analyse_name][event][id_ucsc][rep_test]["TPM"])
                                read_test.append(catalog[analyse_name][event][id_ucsc][rep_test]["NumReads"])

                        tpm_control  = []                        
                        read_control = []

                        for rep_control in replicates_per_condition(config.parameters["samples_for_quantification"],conditions[1]) :  #config.parameters["samples_for_quantification"][sample]
                                tpm_control.append(catalog[analyse_name][event][id_ucsc][rep_control]["TPM"])
                                read_control.append(catalog[analyse_name][event][id_ucsc][rep_control]["NumReads"])
                        
                        features.extend( ["::".join(str(x) for x in read_test)]) 
                        features.extend( ["::".join(str(x) for x in tpm_test)]) 
                        
                        features.extend( ["::".join(str(x) for x in read_control)]) 
                        features.extend( ["::".join(str(x) for x in tpm_control)]) 
                        
                        features.extend( [round(statistics.mean(map(float,read_test)))]) 
                        features.extend( [round(statistics.mean(map(float,tpm_test)))]) 

                        features.extend( [round(statistics.mean(map(float,read_control)))]) 
                        features.extend( [round(statistics.mean(map(float,tpm_control)))]) 

                        lines_for_my_tab.append(features)
            break
        
        logger.info("  Dataframe with pandas : ")

        df = pd.DataFrame(lines_for_my_tab,columns=headerListOfFields,dtype=float)
        df.apply(lambda x: pd.to_numeric(x, errors='coerce') )
        
        onglet= "-"
        if( len(tab) >= 31 ) :
            onglet = tab[:30]
            print(tab)
        else : onglet = tab
        
        logger.info("write tab :"+onglet)
        
        # Convert the dataframe to an XlsxWriter Excel object.
        df.to_excel(writer, sheet_name=onglet,index_label=False, index=False,header=True)
       
        clean_tab_name = tab.replace(" ", "_")

        bed_output_gene=config.parameters['path_to_output']+"/"+tab+"/"+parameters.event+"/"+parameters.event+"_TMP_"+clean_tab_name+".bed"
        logger.info("     Bed output : "+bed_output_gene)
        bed_list = open(bed_output_gene, 'w')

        for line in lines_for_my_tab :
            #logger.info(line[5]) # ajout 4
            elements = line[5].split("\t")
            elements[1]= str(int(elements[1])-1) # Correction for the bed to be 0-based (map(float,read_control))
            #print(line[12].split("::"))
            elements[3]=  elements[3]+"_"+line[4]+"_"+line[18]+"_"+str(min(map(int, line[12].split("::"))))
            bed_list.write("\t".join(elements)+"\n")

        bed_list.close()
        
        postitionGenomiceExon=init.parameters['postitionGenomicExon']
        output_crosslink=config.parameters['path_to_output']+"/"+tab+"/"+parameters.event+"/"+parameters.event+"_"+clean_tab_name+"_posRewrited.bed"
        proc_intersect = subprocess.run(("bedtools intersect -loj -a "+bed_output_gene+" -b "+postitionGenomiceExon+ " > "+output_crosslink),shell=True, check=True) 
        write_subprocess_log(proc_intersect,logger)
        subprocess.run(("rm "+bed_output_gene),shell=True)
        rewriteBed(output_crosslink,bed_output_gene)
        subprocess.run(("rm "+output_crosslink),shell=True)
        
        workbook = writer.book

        worksheet = writer.sheets[onglet]
        #headerListOfFields
        # Format the first column
        # Add the standard url link format.
        url_format = workbook.add_format({
            'font_color': 'blue',
            'underline':  1
        })

        worksheet.set_column('A:C', 5)
        worksheet.set_column('D:E', 20)
        worksheet.set_column('F:F', 30)
        worksheet.set_column('G:G', 10)
        worksheet.set_column('H:H', 10)
        worksheet.set_column('I:BZ', 50)

 

    #Close the Pandas Excel writer and output the Excel file.
    writer.save()
    writer.close()
    logger.info("Output XLS : " +file)
    logger.info ("Output BED : " +bed_output_gene)


    
    logger.info('=========> Finish !')

    
    
