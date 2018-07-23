
'''

:date: Jan 23, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Coverage

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import re
import pprint
import statistics
import numpy
import pandas as pd
import math
import numpy as np
import itertools

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
    script = "coverage_activity"
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
        #logger.info(a)
        cutOffPercentile = np.percentile(a, 50) # return 50th percentile, e.g median.
        
        logger.info("Genes : "+str(len(list_values)))
        logger.info("MIN : "+str(np.amin(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MAX : "+str(np.amax(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MEAN : "+str(np.mean(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("MEDIAN : "+str(np.median(a)))   # Compute the arithmetic mean along the specified axis.
        logger.info("STD : "+str(np.std(a)))   #  Compute the standard deviation along the specified axis.
        logger.info("CUTOFF 50 : "+str(cutOffPercentile))
        
        return cutOffPercentile

def blacklist_gene_under_median_at_least_in_one_sample(catalog,matrice_path,names_test_replicates,names_control_replicates,gene_length,project_name):
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
    
    logger.info("Return all the genes name that sould be skipped for splicing analysis: "+matrice_path )
    data = pd.read_csv(matrice_path, sep=",")
    indexed_data = data.set_index(['Gene'])
    
    logger.info(names_test_replicates)
    logger.info(names_control_replicates)

    all_replicates = list(itertools.chain(names_test_replicates,names_control_replicates))
    logger.info(all_replicates)

    rpkm               = {}
    rpkm_gene_dic      = {}
    blacklist_gene     = []
     
    for replicate in all_replicates :
        
        logger.info("Replicate: "+replicate)
        
        total_mapped_reads       = indexed_data[replicate].sum()
        
        logger.info(total_mapped_reads)
        
        scaling_factor           = total_mapped_reads / 1000000
        rpkm[replicate]          = []
        rpkm_gene_dic[replicate] = []
        
        for gene in  indexed_data.index.values :
            
            m            = re.search('^(\w+)\.(\d+)$', gene) # Ok you miss somes of them ENSEMBL_PARY that you count in the sum but that nothing compared to the whole thing
            if(m):
                ensembl_id     = m.group(1)
                
                if ensembl_id in catalog :

                    catalog[ensembl_id].update(  { replicate  : {"rpkm-like-val":"NaN","expression":"undef"} }  )
                    gene_length_kb = gene_length[ensembl_id]/1000
                    
                    if(indexed_data.get_value(gene,replicate) > 0 ) :
                        rpkm_gene      = (indexed_data.get_value(gene,replicate)/scaling_factor)/gene_length_kb 
                        rpkm[replicate].append(rpkm_gene)
                        rpkm_gene_dic[replicate].append(ensembl_id)
                        continue
                    
                    if(indexed_data.get_value(gene,replicate) == 0 ) :    
                        blacklist_gene.append(ensembl_id)
                        catalog[ensembl_id].update(  {"filter":"YES"}   )
                        catalog[ensembl_id].update(  { replicate  : {"rpkm-like-val":"0","expression":"0"} }  )
                        continue
   
    distrib_path= config.parameters['path_to_output']+"output/"+project_name+"/"+parameters.analyse+".distribution.tsv"
    
    distrib_file = open(distrib_path, 'w')

    distrib_file.write("sample"+"\t"+"cutoff"+"\n")

    for rep in all_replicates :
        logger.info("Replicate : "+rep)

        cut_off = check_distribution(rpkm[rep])
        distrib_file.write(rep+"\t"+str(cut_off)+"\n")

        for index,rpkm_value in enumerate(rpkm[rep]):
            
            
            if rpkm_gene_dic[rep][index] in catalog :
                if (rpkm_value < cut_off ) : 
                    blacklist_gene.append(rpkm_gene_dic[rep][index])
                    catalog[rpkm_gene_dic[rep][index]].update(  { rep  : {"rpkm-like-val":rpkm_value,"expression":"LOW"} }  )
                    catalog[rpkm_gene_dic[rep][index]].update(  {"filter":"YES"}   )

                    continue 
                
                catalog[rpkm_gene_dic[rep][index]].update(  { rep  : {"rpkm-like-val":rpkm_value,"expression":"OK"} }  )
                

                #logger.info(rpkm_gene_dic[rep][index])
                #logger.info(rpkm_value)
    distrib_file.close()  
                  
    return list(set(blacklist_gene)),catalog
                    
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

    logger.info("Check Raw counts : "+matrice_path )
    data = pd.read_csv(matrice_path, sep=",")
    indexed_data = data.set_index(['Gene'])

    for gene in  indexed_data.index.values :
        
        #logger.info(gene)

        m            = re.search('^(\w+)\.(\d+)$', gene)
        if(m):
            ensembl_id = m.group(1)
            if ensembl_id in dict_for_analysis :
                    
                    for control in names_control_replicates :
                        dict_for_analysis[ensembl_id][control].update( {"rawReads":indexed_data.get_value(gene,control)} )
        
                    for test in names_test_replicates :
                        dict_for_analysis[ensembl_id][test].update( {"rawReads":indexed_data.get_value(gene,test)} )
        
            #else : 
                    #for control in names_control_replicates :
                    #    dict_for_analysis[ensembl_id][control].update( {"rawReads":"NaN"} )
        
                    #for test in names_test_replicates :
                        #dict_for_analysis[ensembl_id][test].update(   {"rawReads":"NaN"})


    logger.info("UPGRADE DICTIONARY ")


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
    
    logger.info("createBackgroundGeneList : "+matrice_path )
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
        logger.info("LOOK Gene Size IN: "+gene_path)

        for line in f:

            lineElements  = line.strip().split(",")
            m            = re.search('^(\w+)\.(\d+)$', lineElements[0])
            if(m):
                ensembl_id = m.group(1)
                gene_length[ensembl_id]= int(lineElements[1])
           
    f.close()
   

    return gene_length


def complete_with_quantif(list_all_replicates,dict_for_analysis,tpm_cutoff):
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
                    if(ensembl_id in dict_for_analysis) : 
                        dict_for_analysis[ensembl_id][replicat["replicat_id"]].update(  {"TPM":lineElements[3],"NumReads" :lineElements[4]})
                        if(float(lineElements[3]) > 0 ) : 
                            tpms.append( float(lineElements[3]))
                            numReads.append( float(lineElements[4]))
                #else : 
                    #dict_for_analysis[ensembl_id][replicat["replicat_id"]].update(    {"TPM":"NaN","NumReads" :"NaN"} )

                #print(line)
        f.close()
        # USE DEFINED USER CUT OFF SET FOR ALL SAMPLES . DEFAULT 0
        cutOffDict[replicat["replicat_id"]]=tpm_cutoff


    return dict_for_analysis,cutOffDict

def read_file (path): 
    
    hashContentFile = {}
    with open(path) as file_handle:
        for line in file_handle: 
            if(line==""): continue
            lineElements = line.strip().split("\t")
            hashContentFile[lineElements[0]]=lineElements[1]

   
    file_handle.close() 
    return hashContentFile

def initialise_subset(gene_subset,gene_length) :
    
    dict_for_analysis = {}

    for gene_to_test in gene_subset : 
        
        if (gene_to_test  not in gene_length ) :continue
        dict_for_analysis[gene_to_test] = { "symbol" : gene_subset[gene_to_test] , "length" : gene_length[gene_to_test],"filter": "-" }

        
    return  dict_for_analysis


###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
        '''This is just a main '''
        
        parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
        Thi script will merge all intermediate outputs.  
        Example : 
        python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/coverage.py --config=/home/jean-philippe.villemin/code/configs/MERGE_SPLICING/Reik/WHIPPET/WT.Heart_vs_KO.Heart/SE.splicing_KO2.json --analyse=WT.Heart_vs_KO.Heart --filter=/home/jean-philippe.villemin/data/data/motif_databases/AllMotifsSymbolSearched.tsv 

        '''),formatter_class=argparse.RawDescriptionHelpFormatter)
        
        parser.add_argument("-c","--config",action="store",help="Path to a json file.",required=True,type=str,dest='file_config')
        parser.add_argument("-a","--analyse",action="store",help="Name of analyse",required=True,type=str,dest='analyse')
        parser.add_argument("-f","--filter",action="store",help="File with genes to filter out",required=True,type=str,dest='filter')
        parser.add_argument("-tpm","--tpmNumber",action="store",help="Number of reads per million used to filter out events by expression.",required=False,default=0,type=int,dest='tpmNumber')

        parameters = parser.parse_args()
        
        config  = custom_parser.Configuration(parameters.file_config,"json")
        catalog = {}
        catalog[parameters.analyse] = {}

        logger = create_logger(config.parameters['path_to_output'],"INFO",config.parameters['project'])
        
        logger.info('\n=========> START:')
        logger.info("Analysis : "+parameters.analyse)

        gene_subset = read_file(parameters.filter)

        gene_length=complete_with_gene_length(config.parameters['gene_length'])
        logger.info(config.parameters['gene_length'])
        
        catalog = initialise_subset(gene_subset,gene_length)
        #logger.info(catalog)

        pp = pprint.PrettyPrinter()

        
        names_test_replicates_for_raw_count       =  config.parameters["analysis"][parameters.analyse]["replicates_test"]
        names_control_replicates_for_raw_count    =  config.parameters["analysis"][parameters.analyse]["replicates_control"]
        
        replicates    =  list(itertools.chain(names_test_replicates_for_raw_count,names_control_replicates_for_raw_count))
        
        logger.info ("REPLICATES FOR "+parameters.analyse)
        logger.info (replicates)
        logger.info ("")
        
        control_replicates   =  config.parameters["samples_for_quantification"][config.parameters["analysis"][parameters.analyse]["sample_control_for_quantification"]]
        test_replicates      =  config.parameters["samples_for_quantification"][config.parameters["analysis"][parameters.analyse]["sample_test_for_quantification"]]
        
        genes_blacklisted,catalog = blacklist_gene_under_median_at_least_in_one_sample(catalog,config.parameters["read_count_matrice"],names_test_replicates_for_raw_count,names_control_replicates_for_raw_count,gene_length,config.parameters["project"])
        
        #createBackgroundGeneList(config.parameters["read_count_matrice"],names_test_replicates_for_raw_count,names_control_replicates_for_raw_count,genes_blacklisted,config.parameters['path_to_output']+config.parameters['list_files_splicing'][0]+"/",parameters.analyse)
        
        catalog = complete_with_raw_read_count(config.parameters["read_count_matrice"], catalog,names_test_replicates_for_raw_count,names_control_replicates_for_raw_count)
        
        #logger.info(catalog)
        catalog,cutOffDict = complete_with_quantif(return_all_uniq_replicates_object(config.parameters["samples_for_quantification"]), catalog,parameters.tpmNumber)
       
        #pp.pprint(catalog)
        lines = []
        logger.info("Print : ")
        
        filter_path = config.parameters['path_to_output']+"output/"+ config.parameters['project']+"/"+parameters.analyse+".subset.tsv"
        filter_file = open(filter_path, 'w')
        
        for gene in catalog : 
            header = []
            for key in sorted(catalog[gene].keys()):
                if (key == "length" ) : continue

                header.append(key)
            break
        filter_file.write("Gene"+"\t"+"\t".join(header)+"\t"+"mean"+"\n")
    
        for gene in catalog : 
            
            elements_lines = []
            elements_numeric = []

            for key in sorted(catalog[gene].keys()):
                if (key == "length" ) : continue
                
                if ( key !=  "symbol" and key != "filter" ) :
                    elements_lines.append(str(catalog[gene][key]["rpkm-like-val"]))
                    elements_numeric.append(catalog[gene][key]["rpkm-like-val"])
                    continue 
                
                elements_lines.append(str(catalog[gene][key]))
            allvalues = remove_values_from_list(elements_numeric,"NaN")
            mean =  "NaN"
            if len(allvalues) != 0 : 
                mean =  statistics.mean(map(float, allvalues))

            lines.append(gene+"\t"+"\t".join(elements_lines)+"\t"+str(mean))
            
        for line in lines:
            filter_file.write(line+"\n")
           
        filter_file.close()
   
        
        logger.info ("Finish...")
    
