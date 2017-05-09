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

def parse_all_splicing_files(path_to_dir,list_as_type,dict_for_analysis, readsNumber):
    """
    Parse all files produced by RMATS.
    All files have fixed names and in the same directories.
  
    Args:
        path_to_dir (str): Path to the directory.
        list_as_type (list): List of the as eventt retrieve the filenames in the directory

    Returns:
        splicing_dict (obj): Return all the informations in a structured object.
    
        
    """

    for as_type in list_as_type :
        logger.info("LOOK SPLICING FOR AS_TYPE: "+as_type)
        array_splicing_data = { as_type :{} }
        namefile = ""
        if (parameters.modeCount == "2" ) :
            namefile = "MATS.JunctionCountOnly.txt"
        else  :
            namefile = "MATS.ReadsOnTargetAndJunctionCounts.txt"
               
        with open(path_to_dir+"/"+as_type+"."+namefile) as f:
            count = 0
            #counter_for_test = 0
            for line in f:
                #print(line)
                if count == 0 : 
                    count = count +1
                    continue
                

                #print(count)

                #counter_for_test=counter_for_test+1
                #if counter_for_test == 1000 : break
                id_ucsc_event = ""
                lineElements = line.replace("\"", "").strip().split("\t")
                m            = re.search('^(\w+)\.(\d+)$', lineElements[1])
                ensembl_id   = m.group(1)
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
                key_id_ucsc_event = chrom+":"+(":".join(lineElements[5:11]))+":"+strand  #flankingES    flankingEE   
                #print(as_type+key_id_ucsc_event) 
                #highlight=<DB>.<CHROM>:<START>-<END>#<COLOR>|q
                if (as_type == "A3SS") : 
                    if (strand == "-") :
                        id_ucsc_event = chrom+":"+lineElements[5]+"-"+lineElements[10]+":"+strand  #flankingES    flankingEE
                        highlight = "hg38."+chrom+":"+str(int(lineElements[8])+1)+"-"+lineElements[6]
                    if (strand == "+") :
                        id_ucsc_event = chrom+":"+lineElements[9]+"-"+lineElements[6]+":"+strand  #flankingES    flankingEE
                        highlight = "hg38."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[7]
                        
                if (as_type == "A5SS") :  
                    if (strand == "-") :
                        id_ucsc_event = chrom+":"+lineElements[9]+"-"+lineElements[6]+":"+strand  #flankingES    flankingEE
                        highlight = "hg38."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[7]
                    if (strand =="+") :
                        id_ucsc_event = chrom+":"+lineElements[5]+"-"+lineElements[10]+":"+strand  #flankingES    flankingEE
                        highlight = "hg38."+chrom+":"+str(int(lineElements[8])+1)+"-"+lineElements[6]  
                        
                if (as_type == "SE" ) :   
                    id_ucsc_event = chrom+":"+lineElements[7]+"-"+lineElements[10]+":"+strand  #upstreamES downstreamEE
                    highlight = "hg38."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[6]
                
                if (as_type == "RI" ) :
                    id_ucsc_event = chrom+":"+lineElements[5]+"-"+lineElements[6]+":"+strand  #upstreamES downstreamEE
                    highlight = "hg38."+chrom+":"+str(int(lineElements[8])+1)+"-"+lineElements[9]  

                if (as_type == "MXE" ) :  
                    key_id_ucsc_event = chrom+":"+(":".join(lineElements[5:13]))+":"+strand  #flankingES    flankingEE    
                    #print(as_type+key_id_ucsc_event) 

                    id_ucsc_event = chrom+":"+lineElements[9]+"-"+lineElements[12]+":"+strand #upstreamES downstreamEE
                    ic_sample_1 = lineElements[14]
                    sc_sample_1 = lineElements[15]
                    ic_sample_2 = lineElements[16]
                    sc_sample_2 = lineElements[17]
                    incLevel1   = lineElements[22]
                    incLevel2   = lineElements[23]
                    diffinc     = lineElements[24]
                    pvalue      = lineElements[20]
                    fdr         = lineElements[21]
     
                    highlight = "hg38."+chrom+":"+str(int(lineElements[5])+1)+"-"+lineElements[6]+"|"+"hg38."+chrom+":"+str(int(lineElements[7])+1)+"-"+lineElements[8]

                # FILTER OUT BAD EVENTS
                if (float(pvalue)   < 0.05  and float(fdr)   < 0.01 ) : 
                    count = count+1
                    ''' LOG RATIO // FOLD CHANGE '''
                    incLevel1 = remove_values_from_list(incLevel1.split(","),"NA")
                    incLevel2 = remove_values_from_list(incLevel2.split(","),"NA")
                   
                    average_incLevel1 =  statistics.mean(map(float, incLevel1))
                    average_incLevel2 =  statistics.mean(map(float, incLevel2))
                    
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
                    
                    # List of reads per Inclusion (Control+Test) & Exclusion(Control + Test)
                    list_reads_ic = list_int_ic_sample1 + list_int_ic_sample2
                    list_reads_sc = list_int_sc_sample1 + list_int_sc_sample2
                
                    nbr_zero_sc = list_reads_sc.count(0)
                    
                    nbr_zero_ic = list_reads_ic.count(0)
                    
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

                    
                    if ( test1 + test2 == 0 ) :
                        continue  

                    '''
                    low_read_count = "Under20"
                    
                    #NB READS AT LEAST
                    for read_count in  list_reads_ic  :
                        if read_count >= 20 :
                            low_read_count = "-"
    
                    for read_count in  list_reads_sc   :
                        if read_count >= 20 :
                            low_read_count = "-"
                    
                    if (low_read_count == "Under20") : continue        
                    
                    count_sc = 0
                    count_ic = 0

                    if (nbr_zero_sc  != 0 ) :      
                        if (round(nbr_zero_sc/len(list_reads_sc),1) >=  0.5 ) :
                            count_sc = 1
                            
                    if (nbr_zero_ic  != 0 ) :       
                        if (round(nbr_zero_ic/len(list_reads_ic),1) >=  0.5 ) :    
                            count_ic = 1
                    
                    if ( count_sc + count_ic == 2 ) :
                        low_read_count = "LowInIncAndExc"    
                        continue  
                    '''
                            
                    ic_sample_1_sum = sum(list_int_ic_sample1)/len(ic_sample_1)
                    sc_sample_1_sum = sum(list_int_sc_sample1)/len(sc_sample_1)
                    ic_sample_2_sum = sum(list_int_ic_sample2)/len(ic_sample_2)
                    sc_sample_2_sum = sum(list_int_sc_sample2)/len(sc_sample_2)
                    
                    #if ((ic_sample_1_sum  <=  10 or ic_sample_2_sum <= 10) or (sc_sample_1_sum  <=  10 or sc_sample_2_sum <= 10 ) ) :
                        #low_read_count = "LowSumReadsPerGroup"    
                

                           
                    ''' FISHER TEST '''
                    #print(str(ic_sample_1_sum)+" ->  "+str(ic_sample_2_sum)+ " -> "+str(sc_sample_1_sum)+" -> "+str(sc_sample_2_sum))
                    pvalue_fisher = 0
                    oddsratio, pvalue_fisher = stats.fisher_exact([[ic_sample_1_sum, sc_sample_1_sum], [ic_sample_2_sum, sc_sample_2_sum]])
                    fisher = "-"
                    #print(pvalue)
                    if (pvalue_fisher < 0.05) : fisher = "PASS"
                    
                    array_splicing_data[as_type][key_id_ucsc_event] = { "Ensembl":ensembl_id, 
                                                                    "Symbol":gene, 
                                                                    "Chromosome":chrom, 
                                                                    "Strand":strand,
                                                                    "ic_sample_1" : ic_sample_1,
                                                                    "sc_sample_1" : sc_sample_1,
                                                                    "ic_sample_2" : ic_sample_2,
                                                                    "sc_sample_2" : sc_sample_2,
                                                                    "incLevel1"  : incLevel1,
                                                                    "highlight" : highlight,
                                                                    "incLevel2"   : incLevel2,
                                                                    "diffinc"     : diffinc,
                                                                    "logRatioIncLevel"     : logratio,
                                                                    'fdr' : fdr,
                                                                    'log10fdr' :  int(abs(math.log10(float(fdr)))) if float(fdr) > 0 else -math.inf,

                                                                    "FCIncLevel"     : retval,
                                                                    "pvalueFisher" : pvalue_fisher,
                                                                    "pval" : pvalue,
                                                                    "fisher"     : fisher,
                                                                    "id_ucsc"     : id_ucsc_event[:-2]
                                                                    #'low_read_count' : low_read_count
                                                                   
    
                                                                  }
        #logger.info("NUMBER OF EVENTS         : "+str(count))
        f.close()
        dict_for_analysis.update(array_splicing_data)

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
        logger.info("    UPGRADE DICTIONNARY FOR EVENT : "+event)

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

    Returns:
        splicing_dict (obj): Return all the informations in a structured object.(gene FC added)
    
        
    """
    quantification            = {}
    #list_of_replicates       = []

    #for replicate_object in list_all_replicates:
        #list_of_replicates.append([replicate_object["replicat_id"],replicate_object["file_path"]])
    cutOffDict = {}
    for replicat in list_all_replicates:
        #print(replicat)
        logger.info("LOOK QUANTIFICATION IN : "+replicat["replicat_id"] )
        tpms = []
        numReads = []
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
        
        #print(tpms)
        a = np.array(tpms)
        cutOffPercentile = np.percentile(a, 50) # return 50th percentile, e.g median.
        logger.info("Genes : "+str(len(tpms)))
        logger.info("MIN : "+str(np.amin(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MAX : "+str(np.amax(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MEAN : "+str(np.mean(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MEDIAN : "+str(np.median(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("STD : "+str(np.std(a)))   #  Compute the standard deviation along the specified axis.
        logger.info("CUTOFF 50 : "+str(cutOffPercentile))
        cutOffDict[replicat["replicat_id"]]=cutOffPercentile
        
        b = np.array(numReads)
        logger.info("Genes : "+str(len(numReads)))
        logger.info("MIN : "+str(np.amin(b)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MAX : "+str(np.amax(b)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MEAN : "+str(np.mean(b)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MEDIAN : "+str(np.median(b)))   # Compute the arithmetic mean along the specified axis.
        
        
    for event in dict_for_analysis.keys():
     
            logger.info("    UPGRADE DICTIONARY IN "+event)

            for id_ucsc_event in dict_for_analysis[event].keys():
        
                #print("id_ucsc_event : "+id_ucsc_event)

                if  dict_for_analysis[event][id_ucsc_event]["Ensembl"] in quantification :
                    
                    dict_for_analysis[event][id_ucsc_event].update(quantification[dict_for_analysis[event][id_ucsc_event]["Ensembl"]])
                else : 
                    logger.info("WTF")

            #pp.pprint(analysis[event][id_ucsc_event])

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
    

            
            
def create_header(dict_samples, list_analysisTocheck):
    """
    Create header of the file
  
    Args:
        dict_samples(dict(list(dict))): Dict of samples in json config
        list_analysis_to_check:  List of _analysis you want in tab


    Returns:
        headerListOfFields (list): Return list of fields for the header
   
    """
    headerListOfFields                                    = []
    list_all_TPM_replicates_in_samples_for_quantification = []
    header_variabe = ["Test:inclusion.reads|exclusion.reads","Control:inclusion.reads|exclusion.reads","Test.PSI|Control.PSI","PSI-Diff","PSI-Log2Ratio","PSI-FoldChange","PSI-Pval","LOG10-PSI-FDR","Fisher","Gene-Log2FoldChange","Gene-FoldChange","GenePadj","Match"] #,""
    coreHeaderFields =  ['ID-UCSC-HG38','ID-UCSC-HG19','Epissage','Event','Symbol','Ensembl','Window Edges HG38|HG19',"Inc/Excl Edges HG38",'Strand','Gene_biotype']               

    
    for sample in sorted(dict_samples.keys()) :
        # ordered because it's a
        for replicat in dict_samples[sample]  :
            list_all_TPM_replicates_in_samples_for_quantification.append('TPM-'+replicat["replicat_id"])
            list_all_TPM_replicates_in_samples_for_quantification.append("Is-"+replicat["replicat_id"]+"-Expressed")
            list_all_TPM_replicates_in_samples_for_quantification.append("Is-"+replicat["replicat_id"]+"-EstimatedNumReads")

    #keep unique only
    tpm_header  = list(OrderedDict.fromkeys(list_all_TPM_replicates_in_samples_for_quantification))
    
    logger.info("  TPM_HEADER  :: "+" - ".join(tpm_header))

    for analyse_name1 in (list_analysisTocheck) :
               
            for x in range(0,len(header_variabe)) :
                    
                headerListOfFields.append(analyse_name1+"-"+header_variabe[x])

    
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
    parser.add_argument("-mc","--modeCount",action="store", default=1 ,help="with or without ReadsOnTarget",required=False,type=str,dest='modeCount')
    parser.add_argument("-r","--readsNumber",action="store",help="Number of reads used to filter out events.",required=False,default=10,type=int,dest='readsNumber')


    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")
    
    subprocess.run(("mkdir -p "+config.parameters['path_to_output']+"/"+config.chrono),shell=True)

    logger = create_logger(config)

    logger.info('=========> START:')
    
    logger.info("file_config : "+parameters.file_config)
    logger.info("path_to_output : "+config.parameters['path_to_output'])
    logger.info("    ")

    namefile = ""
    if (parameters.modeCount == "2" ) :
        namefile = "MATS.JunctionCountOnly.txt"
    else  :
        namefile = "MATS.ReadsOnTargetAndJunctionCounts.txt"
    
    '''
    ######################################################################
    DataFrame Object Construction : First Part
    ######################################################################
    '''
    lo = LiftOver('hg38', 'hg19')
    catalog         = {}
    cutOffDict      = {}
    #read_count         = {"str1":{},"str2":{}}
    for analyse_name_loop in (config.parameters.get("analysis").keys() ):
             
        catalog[analyse_name_loop] = {}
        #logger.info(" ".join(list(range(0,6))))
    
        logger.info("    Analyse_name : "+analyse_name_loop)
        logger.info ("    splicing : "+config.parameters["analysis"][analyse_name_loop]["splicing_dir"])
        logger.info ("    expression : "+config.parameters["analysis"][analyse_name_loop]["expression_file_path"])
        logger.info ("    sample_control_for_quantification : "+config.parameters["analysis"][analyse_name_loop]["sample_control_for_quantification"])
        logger.info ("    sample_test_for_quantification : "+config.parameters["analysis"][analyse_name_loop]["sample_test_for_quantification"])
        logger.info ("    list_splicing_files : "+" ".join(config.parameters["list_files_splicing"]))
        logger.info ("")

        test = config.parameters["analysis"][analyse_name_loop]["sample_test_for_quantification"]

        for replicat_id in config.parameters["samples_for_quantification"][test] :
            logger.info ("    Replicat_id - test : "+replicat_id.get("replicat_id"))
        logger.info ("")
        
        control = config.parameters["analysis"][analyse_name_loop]["sample_control_for_quantification"]

        for replicat_id in config.parameters["samples_for_quantification"][control] :
            logger.info ("    Replicat_id - control : "+replicat_id.get("replicat_id"))
       
        logger.info ("GO...")
        logger.info ("")

        catalog[analyse_name_loop] = parse_all_splicing_files(config.parameters["analysis"][analyse_name_loop]["splicing_dir"],config.parameters["list_files_splicing"], catalog[analyse_name_loop],    parameters.readsNumber)
        
        logger.info ("")

        catalog[analyse_name_loop] = complete_with_expression(config.parameters["analysis"][analyse_name_loop]["expression_file_path"], catalog[analyse_name_loop])
        
        logger.info ("")
        
        control_replicates   =  config.parameters["samples_for_quantification"][config.parameters["analysis"][analyse_name_loop]["sample_control_for_quantification"]]
        test_replicates      =  config.parameters["samples_for_quantification"][config.parameters["analysis"][analyse_name_loop]["sample_test_for_quantification"]]
  
        #catalog[analyse_name] = complete_with_quantif(control_replicates, catalog[analyse_name])
        #catalog[analyse_name] = complete_with_quantif(test_replicates, catalog[analyse_name])
 
        catalog[analyse_name_loop],cutOffDict = complete_with_quantif(return_all_uniq_replicates_object(config.parameters["samples_for_quantification"]), catalog[analyse_name_loop])
        #print(cutOffDict)
        logger.info ("")
        logger.info ("Next...")
        logger.info ("")
   
    logger.info ("Finish...")
    
    print(cutOffDict)
    
    logger.info ("######################################################################")
    logger.info ("That the tricky part to structure data to fit your envy")
    logger.info ("######################################################################")

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(config.parameters['path_to_output']+"/"+config.chrono+'/Events.'+config.parameters["list_files_splicing"][0]+"."+namefile+'.xlsx', engine='xlsxwriter')
    logger.info (namefile)

    dict_analysis               = config.parameters.get("analysis")
    dict_samples                = config.parameters["samples_for_quantification"]
    list_analysis               = list(sorted(config.parameters.get("analysis").keys()))

    
    ''' That the tricky part to structure data to fit your envy'''
    #Foreach tab you want to create
    for tab in config.parameters.get("tabs").keys():
    
        #output = open(config.parameters["path_to_output"]+"/"+analyse_name+".csv",'w')
        #list_output_csv.append(config.parameters["path_to_output"]+"/"+tab+".csv")
        logger.info("")
        logger.info("  TAB  :: "+tab)
        logger.info("  analysisTotcheck  :: "+" - ".join(config.parameters["tabs"][tab]["analysisTocheck"]))
        
        lines_for_my_tab = []
        list_analysis_with_event    = set(config.parameters["tabs"][tab]["analysisTocheck"])
        
        list_analysisNotTotcheck = [ x for x in list_analysis if x not in list_analysis_with_event ]
        logger.info("  analysisNotTocheck  :: "+" - ".join(list_analysisNotTotcheck))
        logger.info("  GO...")
            
        headerListOfFields,header_variabe = create_header(dict_samples,config.parameters["tabs"][tab]["analysisTocheck"])
        length_fields_sometimes_empty = len(header_variabe)
        
        # Column fields list are specific of what you compare in each tab...
        #columnNames_to_remove_per_tab = []
        #for analysis_to_ban in list_analysisNotTotcheck : 
        #    for indice in range(0,len(header_variabe)) :
        #        columnNames_to_remove_per_tab.append(analysis_to_ban+"-"+header_variabe[indice])
        #indexes_to_remove = [i for i, item in enumerate(headerListOfFields) if item in set(columnNames_to_remove_per_tab)] 
       
        #for index in sorted(indexes_to_remove, reverse=True):
        #    del headerListOfFields[index]
        
        print("Header : " +"\n".join(headerListOfFields)) 
        # You are going to check the first analysis in list vs all others for you tab definition 
        for analyse_name in config.parameters["tabs"][tab]["analysisTocheck"]:
          
            logger.info("  ") 
            logger.info("    ANALYSE TO CHECK :: "+analyse_name)

            # By event
          
            for event in config.parameters["list_files_splicing"] :

                logger.info("          EVENT  :: "+event)
                logger.info("          NB EVENTS  :: "+str(len(catalog[analyse_name][event])))
                logger.info("     Run trought the list of keys coordinates in current analysis...")
                counter = 0
                
                
                for id_ucsc in catalog[analyse_name][event]:
                   
                    #logger.info(id_ucsc+" "+analyse_name+" "+event+str(counter))
                    counter=counter + 1
                    coordinates = id_ucsc.split(":")
                    chrom       = coordinates[0]
                    strand      = coordinates[-1]
                    #print(id_ucsc)
                    regions = ""
                    for i in range(1,len(coordinates)-2,2) :
                        regions+=chrom+" "+coordinates[i]+" "+coordinates[i+1]+" "+"Region"+str(i)+" "+"."+" "+strand+" \n"
                    #print("        id_ucsc: "+id_ucsc)
                    #print("        reg:\n"+regions)

                    #regionsBedTools = pybedtools.BedTool(regions,from_string=True)

                   
                    if (float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > 0 ) :  epissage = 'INC' 
                    else : 
                        epissage = 'EXCL'
   
                    #  Check For Presence in 
                    analysis_to_ucscKey = {}
                    analysis_to_ucscKey[analyse_name] = id_ucsc
                    
                    #logger.info("            Run trought the different analysis to check it the event exist ...")


                    '''
                    Remove the envents by TPM filter for the one you test
                    
                    '''
                    break_because_of_tpm = "no" 
                    testAgain_test1 = config.parameters["analysis"][analyse_name]["sample_test_for_quantification"]

                    for replicat_id in config.parameters["samples_for_quantification"][testAgain_test1] :
                            
                        if (float(catalog[analyse_name][event][id_ucsc][replicat_id.get("replicat_id")]["TPM"]) < cutOffDict[replicat_id.get("replicat_id")] ) : 
                            break_because_of_tpm = "yes"
                            break
                        
                    testAgain_control1 = config.parameters["analysis"][analyse_name]["sample_control_for_quantification"]

                    for replicat_id in config.parameters["samples_for_quantification"][testAgain_control1] :
                        
                        if (float(catalog[analyse_name][event][id_ucsc][replicat_id.get("replicat_id")]["TPM"]) < cutOffDict[replicat_id.get("replicat_id")] ) : 
                            break_because_of_tpm = "yes"
                            break
                        
                    if(break_because_of_tpm == "yes" ) : continue


                    index = 0

                    for another_analyse1 in config.parameters["tabs"][tab]["analysisTocheck"]:
                        
                        if (analyse_name == another_analyse1   ) : continue

                        #logger.info(another_analyse1)
                        index =  config.parameters["tabs"][tab]["analysisTocheck"].index(another_analyse1)

                        if id_ucsc in catalog[another_analyse1][event] :
                            
                                '''
                                Remove the envents by TPM filter for the other you test
                        
                                '''
                            
                                break_because_of_tpm = "no" 
                                testAgain_test2 = config.parameters["analysis"][another_analyse1]["sample_test_for_quantification"]
            
                                for replicat_id in config.parameters["samples_for_quantification"][testAgain_test2] :
                                        
                                    if (float(catalog[another_analyse1][event][id_ucsc][replicat_id.get("replicat_id")]["TPM"]) < cutOffDict[replicat_id.get("replicat_id")] ) :                        
                                        break_because_of_tpm = "yes" 
                                        break
        
                                testAgain_control2 = config.parameters["analysis"][another_analyse1]["sample_control_for_quantification"]
            
                                for replicat_id in config.parameters["samples_for_quantification"][testAgain_control2] :
                                        
                                    if (float(catalog[another_analyse1][event][id_ucsc][replicat_id.get("replicat_id")]["TPM"]) < cutOffDict[replicat_id.get("replicat_id")] ) :                         
                                        break_because_of_tpm = "yes" 
                                        break
        
                                if(break_because_of_tpm == "yes" ) : break
                                
                                if (math.isnan(float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]))== True or math.isnan(float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]))==True) :
                                    analysis_to_ucscKey[another_analyse1] = id_ucsc
                                    break
                            
                                if ( config.parameters["tabs"][tab][str(index)] == "lower" ) :
                                        
                                            if (abs(float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"])) >= abs(float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"])) ) :
                                            
                                                if (( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) < 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) < 0)  or ( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) > 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) > 0) ): 
                                                    analysis_to_ucscKey[another_analyse1] = id_ucsc
                                                    break
                                            
                                if ( config.parameters["tabs"][tab][str(index)] == "upper" ) :
                                        
                                            if (abs(float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"])) <= abs(float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"])) ) :
                                                
                                                if (( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) < 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) < 0)  or ( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) > 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) > 0) ): 
                                                    analysis_to_ucscKey[another_analyse1] = id_ucsc
                                                    break
                                       
                                if ( config.parameters["tabs"][tab][str(index)] == "inverse" ) :
                                        
                                            if (( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) > 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) < 0)  or ( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) < 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) > 0) ): 
                                                    analysis_to_ucscKey[another_analyse1] = id_ucsc
                                                    break            
                                                
                                if ( config.parameters["tabs"][tab][str(index)] == "idem" ) :
                                       
                                        if (( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) < 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) < 0)  or ( float(catalog[analyse_name][event][id_ucsc]["FCIncLevel"]) > 0 and float(catalog[another_analyse1][event][id_ucsc]["FCIncLevel"]) > 0) ): 

                                            analysis_to_ucscKey[another_analyse1] = id_ucsc
                                            break  
                         
                        else : 
                            
                            if ( config.parameters["tabs"][tab][str(index)] == "notIn" ) :
                                       
                                analysis_to_ucscKey[another_analyse1] = id_ucsc
                                          
                            
                            continue
                            '''
                            #print(another_analyse1)
                            # Tricky change id_ucsc of id_uccsc2
                            for id_ucsc2 in catalog[another_analyse1][event] :
    
                                coordinates2 = id_ucsc2.split(":")
                                chrom2       = coordinates2[0]
                                strand2      = coordinates2[-1]
                                
                                regions2 = ""
                                for i2 in range(1,len(coordinates2)-2,2) :
                                    regions2+=chrom2+" "+coordinates2[i2]+" "+coordinates2[i2+1]+" "+"RegionBis"+str(i2)+" "+"."+" "+strand2+" \n"
                                    
                                #print("        id_ucsc2 :"+id_ucsc2)
                                #print("        reg2 :\n"+regions2)
    
                                regionsBedTools2 = pybedtools.BedTool(regions2,from_string=True)
    
                                intersection = regionsBedTools.intersect(regionsBedTools2,u=True)
                                
                                if( (len(intersection) >= 4 and event == "MXE") or (len(intersection) >= 3  and event != "MXE" ) ): 
                                    
                                    # Check they are going in the same way...
                                    if( float(catalog[analyse_name][event][id_ucsc]["diffinc"]) < 0 and float(catalog[another_analyse1][event][id_ucsc2]["diffinc"]) > 0  ):     
                                        continue      
                                    if( float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > 0 and float(catalog[another_analyse1][event][id_ucsc2]["diffinc"]) < 0):     
                                        continue
                                    
                                    # Check how comparison should behave
                                    index =  config.parameters["tabs"][tab]["analysisTocheck"].index(another_analyse1)
                                    
                                    if ( config.parameters["tabs"][tab][str(index)] == "upper" ) :
                                        
                                        if( float(catalog[analyse_name][event][id_ucsc]["diffinc"]) < 0 ) :
                                            
                                            if (float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > float(catalog[another_analyse1][event][id_ucsc2]["diffinc"]) ) :
                                            
                                                analysis_to_ucscKey[another_analyse1] = id_ucsc2
                                                break
                                            
                                        if( float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > 0 ) :
                                            
                                            if (float(catalog[analyse_name][event][id_ucsc]["diffinc"]) < float(catalog[another_analyse1][event][id_ucsc2]["diffinc"]) ) :
    
                                                analysis_to_ucscKey[another_analyse1] = id_ucsc2
                                                break
                                            
                                    if ( config.parameters["tabs"][tab][str(index)] == "lower" ) :
                                        
                                        if( float(catalog[analyse_name][event][id_ucsc]["diffinc"]) < 0 ) :
                                            
                                            if (float(catalog[analyse_name][event][id_ucsc]["diffinc"]) < float(catalog[another_analyse1][event][id_ucsc2]["diffinc"]) ) :
                                                
                                                analysis_to_ucscKey[another_analyse1] = id_ucsc2
                                                break
                                            
                                        if( float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > 0 ) :
                                    
                                            if (float(catalog[analyse_name][event][id_ucsc]["diffinc"]) > float(catalog[another_analyse1][event][id_ucsc2]["diffinc"]) ) :
    
                                                analysis_to_ucscKey[another_analyse1] = id_ucsc2
                                                break
                                            
                                    if ( config.parameters["tabs"][tab][str(index)] == "any" ) :
                                       
                                        analysis_to_ucscKey[another_analyse1] = id_ucsc2
                                        break                                          


                            '''
                    '''
                    We find the event in the others analysis
                    '''
                    if (len(analysis_to_ucscKey) == len(config.parameters["tabs"][tab]["analysisTocheck"]) ) :            
                       
                
                        features         = []
                        chro             = catalog[analyse_name][event][id_ucsc]["Chromosome"]
                        id_ucsc_clean    = catalog[analyse_name][event][id_ucsc]["id_ucsc"]
                        #logger.info(event+" "+id_ucsc_clean)
                        ######UGLY Conversion UGLY#####
                        m            = re.search('^(.*):(\d+)-(\d+)$', id_ucsc_clean)
    
                        starthg19 = lo.convert_coordinate(chro,int(m.group(2)),catalog[analyse_name][event][id_ucsc]["Strand"])
                        endhg19   = lo.convert_coordinate(chro, int(m.group(3)),catalog[analyse_name][event][id_ucsc]["Strand"])
                        
                        highlight_hg19 = "" 
                        if(event == "MXE") :
                            
                            m1           = re.search('^hg38.(.*):(\d+)-(\d+)\|hg38.(.*):(\d+)-(\d+)$',catalog[analyse_name][event][id_ucsc]["highlight"])
                            
                            starthg19_e1 = lo.convert_coordinate(chro,int(m1.group(2)),catalog[analyse_name][event][id_ucsc]["Strand"])
                            endhg19_e1   = lo.convert_coordinate(chro, int(m1.group(3)),catalog[analyse_name][event][id_ucsc]["Strand"])
                            #print("e1"+str(starthg19_e1[0][1])+"-"+str(endhg19_e1[0][1]))
                            starthg19_e2 = lo.convert_coordinate(chro,int(m1.group(5)),catalog[analyse_name][event][id_ucsc]["Strand"])
                            endhg19_e2   = lo.convert_coordinate(chro, int(m1.group(6)),catalog[analyse_name][event][id_ucsc]["Strand"])
                            #print("e1"+str(starthg19_e2[0][1])+"-"+str(endhg19_e2[0][1]))
                            # ADDED +1 TO CORRECT  O-BASED ERROR DISPLAY
                            highlight_hg19 = "&highlight=hg19."+chro+":"+str(starthg19_e1[0][1]+1 if  starthg19_e1 else "None" )+"-"+str(endhg19_e1[0][1] if  endhg19_e1 else "None")+"|"+chro+":"+str(starthg19_e2[0][1]+1 if  starthg19_e2 else "None" )+"-"+str(endhg19_e2[0][1] if  endhg19_e2 else "None" )
                            
                        else : 
                            
                            m2            = re.search('^hg38.(.*):(\d+)-(\d+)$',catalog[analyse_name][event][id_ucsc]["highlight"])
                            starthg19_simple = lo.convert_coordinate(chro,int(m2.group(2)),catalog[analyse_name][event][id_ucsc]["Strand"])
                            endhg19_simple   = lo.convert_coordinate(chro, int(m2.group(3)),catalog[analyse_name][event][id_ucsc]["Strand"])
                            #print(starthg19_simple)
                            #print(endhg19_simple)
                            # ADDED +1 TO CORRECT  O-BASED ERROR DISPLAY
                            highlight_hg19 = "&highlight=hg19."+chro+":"+str(starthg19_simple[0][1]+1 if  starthg19_simple else "None")+"-"+str(endhg19_simple[0][1] if  endhg19_simple else "None")
    
                        features.extend([ 
                                         "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jp&hgS_otherUserSessionName=EMT_RNASEQ_hg38&position="+id_ucsc_clean+"&highlight="+catalog[analyse_name][event][id_ucsc]["highlight"],
                                         "http://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=jp&hgS_otherUserSessionName=Eneritz_EMT&position="+chro+":"+str(starthg19[0][1] if  starthg19 else "None" )+"-"+str(endhg19[0][1] if  endhg19 else "None" )+highlight_hg19,
                                         epissage,
                                         event,
                                         catalog[analyse_name][event][id_ucsc]["Symbol"],
                                         catalog[analyse_name][event][id_ucsc]["Ensembl"],
                                         id_ucsc_clean+"|"+chro+":"+str(starthg19[0][1]  if  starthg19 else "None" )+"-"+str(endhg19[0][1]  if  endhg19 else "None"),
                                         chro+" "+catalog[analyse_name][event][id_ucsc]["highlight"][catalog[analyse_name][event][id_ucsc]["highlight"].index(":"):][1:].replace("-"," ")+" "+catalog[analyse_name][event][id_ucsc]["diffinc"]+"_"+catalog[analyse_name][event][id_ucsc]["Symbol"]+" "+str(catalog[analyse_name][event][id_ucsc]["FCIncLevel"])+" "+catalog[analyse_name][event][id_ucsc]["Strand"] ,

                                         catalog[analyse_name][event][id_ucsc]["Strand"],
                                         catalog[analyse_name][event][id_ucsc]["gene_biotype"]
                                        ]
                                        )
                        
                        # Add Gnene FC , Inclusion Level in each other analysis present in the lab
                        for analyse_name_bis in config.parameters["tabs"][tab]["analysisTocheck"] :
                                
                                if analyse_name_bis in analysis_to_ucscKey :
                                    
                                    if ( config.parameters["tabs"][tab][str(index)] != "notIn" ) :
                                      
                                        features.extend([
                                                          str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["ic_sample_1"]).replace(",","::")+"|"+str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["sc_sample_1"]).replace(",","::")
                                                         
                                                          ,str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["ic_sample_2"]).replace(",","::")+"|"+str(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["sc_sample_2"]).replace(",","::")
                                                          
                                                          ,"::".join(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["incLevel1"])+"|"+"::".join(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["incLevel2"])
                                                          ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["diffinc"]
                                                          ,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["logRatioIncLevel"])
                                                          ,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["FCIncLevel"])
                                                          ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["pval"].replace("e","E")
                                                         ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["log10fdr"]
                                                          ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["fisher"]
                                                          ,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["log2fc"] )
                                                          ,return_formated(catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["fc"] )
                                                           ,catalog[analyse_name_bis][event][analysis_to_ucscKey[analyse_name_bis]]["padj"]
                                                           ,analysis_to_ucscKey[analyse_name_bis]
                                                         ])
                                    else :   
                                        for i in range(0,length_fields_sometimes_empty) : features.extend(["."])
    
    
                        for replicat in return_all_uniq_replicates_object(config.parameters["samples_for_quantification"]) :  #config.parameters["samples_for_quantification"][sample]
                                
                            express = "-"
                            if (float(catalog[analyse_name][event][id_ucsc][replicat["replicat_id"]]["TPM"]) < cutOffDict[replicat["replicat_id"]] ) :
                                express = "No"
                            features.extend([catalog[analyse_name][event][id_ucsc][replicat["replicat_id"]]["TPM"],express,catalog[analyse_name][event][id_ucsc][replicat["replicat_id"]]["NumReads"]]) 
                        
                        lines_for_my_tab.append(features)
            break

        df = pd.DataFrame(lines_for_my_tab,columns=headerListOfFields,dtype=float)
        df.apply(lambda x: pd.to_numeric(x, errors='coerce') )
        # Convert the dataframe to an XlsxWriter Excel object.
        df.to_excel(writer, sheet_name=tab,index_label=False, index=False,header=True)
        #output.close()              

    workbook = writer.book
    for tab in config.parameters.get("tabs").keys():

        worksheet = writer.sheets[tab]
       
        #headerListOfFields
        # Format the first column
        # Add the standard url link format.
        url_format = workbook.add_format({
            'font_color': 'blue',
            'underline':  1
        })

        worksheet.set_column('A:D', 5)
        worksheet.set_column('E:E', 12)
        worksheet.set_column('F:F', 20)
        worksheet.set_column('G:G', 25)
        worksheet.set_column('H:H', 25)
        worksheet.set_column('I:I', 5)
        worksheet.set_column('J:BZ', 20)



    #Close the Pandas Excel writer and output the Excel file.
    writer.save()
    #writer.close()
    logger.info('=========> Finish !')

    
    
