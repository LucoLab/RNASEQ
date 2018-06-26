
'''

:date: June 13, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Whole Pipeline/Wrapper for RNA-SEQ

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import glob
import pandas as pd
from collections import OrderedDict
import os
import json
import re
import sys
import numpy
import gzip
import pathlib
import shutil
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
    script = "signature_activity"
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


def write_align_conf(project,id_condition,projet_id,uniq_id_sample,path,path,repNumber):
    """
    Write Json configuration file for alignment.
  
    Args:
        project (dict): project.
        id_condition (str): id_condition.
        projet_id (str): projet_id.
        uniq_id_sample (dict): uniq_id_sample.
        path (str): path to files.
        path (str): path to files.
        repNumber (str): number of the replicate.

    """  

    logger.info("\n-> Write Alignment & TPM values generation JSON Configuration File \n")  

    with open(path +projet_id+"."+uniq_id_sample+"_align.json", 'w') as outfile:  

        outfile.write("{\n")
        outfile.write('"project" : "'+projet_id+'" , \n')
        outfile.write('"star":\n')
        outfile.write('{"runThreadN":"32",\n'+
        '"genomeDir":"/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/indexes/GRCh38_PRIM_GENCODE_R25/",\n'+
        '"outSAMtype":"BAM SortedByCoordinate"},\n')
        outfile.write('"fastqc":\n'+
                    '{  "threads":"24"},\n'+
                    '"salmon":\n'+
                    '{ "cpu":"24"},\n'+
                    '"samtools":\n'+
                    '{ "cpu":"24"},\n'+
                    '"softs":\n'+
                    '{ "samtools":"samtools","star":"STAR",  "wigToBigWig":"wigToBigWig"},\n')
        outfile.write('"type": "'+project[projet_id][uniq_id_sample]["pair"].lower()+'End"\n, "gzip": "yes"\n, ')\
        #"strand": "'+project[projet_id][uniq_id_sample]["strand"].capitalize()+'",\n')
        logger.info(project[projet_id][uniq_id_sample]["pair"].lower())
        #logger.info(project[projet_id][uniq_id_sample]["strand"].capitalize())

        outfile.write('"files":{\n')
        url_fastq = project[projet_id][uniq_id_sample]["url"].split(";")
        filePath =""
        if(project[projet_id][uniq_id_sample]["pair"]=="PAIRED") :
            
            url_fastq = project[projet_id][uniq_id_sample]["url"].split(";")
            url_1 = url_fastq[0].replace("\"","")
            url_2 = "None" if len(url_fastq)==1 else url_fastq[1].replace("\"","")
            
            # if ftp url
            if("ftp" in url_1) :
                filePath = path+""+os.path.basename(url_1.replace("ftp://ftp.sra.ebi.ac.uk",""))
                # You check file is already DL else you download
                if not os.path.isfile(filePath):
                    logger.info("Download")
                    #download_if_not_there_already()

            if("ftp" in url_2) :
                filePath = path+""+os.path.basename(url_1.replace("ftp://ftp.sra.ebi.ac.uk",""))
                # You check file is already DL else you download
                if not os.path.isfile(filePath):
                    logger.info("Download")
                #download_if_not_there_already()
          
 
            outfile.write(' "GroupeA": { "R1": "'+os.path.basename(url_1.replace("ftp://ftp.sra.ebi.ac.uk",""))+'","R2": "'+os.path.basename(url_2.replace("ftp://ftp.sra.ebi.ac.uk",""))+'"}\n') 
        else : 
            
            url_fastq = project[projet_id][uniq_id_sample]["url"].split(";")
            url_1    = url_fastq[0].replace("\"","")
            outfile.write(' "GroupeA": { "R1": "'+os.path.basename(url_1.replace("ftp://ftp.sra.ebi.ac.uk",""))+'","R2": "None"}\n') 
            
            filePath = path+""+os.path.basename(url_1.replace("ftp://ftp.sra.ebi.ac.uk",""))

            if("ftp" in url_1) :
                # You check file is already DL else you download
                if not os.path.isfile(filePath):
                    logger.info("Download")
                #download_if_not_there_already()

        #json.dump(dataformatedforexport2Json, outfile)

        #{ "GroupeA": { "R1": "SRR3578179.fastq.gz","R2": "None"}  }
        outfile.write("},\n")
        outfile.write('"path_to_input" : "'+path+'",\n'+
              '"path_to_output": "'+path+'",\n'+
              '"path_to_chrom_length_ref":"/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genome/GRCh38_PRIM_GENCODE_R25/GRCh38.primary_assembly.genome.sizes",\n'+
              '"final_bam_name": "'+projet_id+"."+uniq_id_sample+"."+repNumber+'",\n'
              '"path_to_gtf":"/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/gencode.v25.primary_assembly.annotation.gtf",\n'+
              '"transcriptome_index":"/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/transcriptome/GRCh38_PRIM_GENCODE_R25/gencode.v25.transcripts.index"\n')

        outfile.write("}")
    outfile.close()
  
def write_diff_exp_conf(project,projet_id_current,path,hash_fc,analysis,samples_by_project_conditions,libtype):
    """
    Write Json configuration file for differential expression.
  
    Args:
        project (dict): project.
        projet_id_current (str): projet_id_current.
        path_to_out (str): path.
        hash_fc (dict): hash_fc.
        analysis (str): analysis.
        samples_by_project_conditions (dict): samples_by_project_conditions.
        libtype (str): libtype.

    """  
    
    logger.info("\n-> Write Differential Expression JSON Configuration File \n")  
    
    for projet_id in sorted(project.keys()) :

        if(projet_id==projet_id_current) :
            logger.info(path+projet_id)
            outfile = open(path+projet_id+"_"+"diff_exp.json", 'w')  
            outfile.write(" {\n")
            outfile.write('"project" : "'+projet_id+'" , \n')

            outfile.write('   "files" : [ \n')
            tot_lines=0
            for condition in samples_by_project_conditions[projet_id] : 
                for sample_id in samples_by_project_conditions[projet_id][condition] : 
                    tot_lines+=1
            a=0         
            for condition in samples_by_project_conditions[projet_id] : 
                logger.info(condition)
               
                for sample_id in samples_by_project_conditions[projet_id][condition] : 
                    outfile.write(  '{"sample_id": "'+projet_id+"."+sample_id+'" ,"condition": "'+condition+'"}')
                    logger.info(sample_id)
                    if(a!=tot_lines-1) :  outfile.write(",\n")
                    a+=1

            outfile.write('   ], \n')
            outfile.write('"organism": "human", \n')
            outfile.write('     "lib-type": "'+ libtype+'", \n')
            outfile.write('     "path_to_output":"'+ path+'", \n')
            outfile.write('     "path_to_input":"'+ path+'" \n')
            outfile.write(' } \n')
            outfile.close()

   # if(config.parameters.get("lib-type")=="Reverse") :
   #     index=2
   # if(config.parameters.get("lib-type")=="Yes") :
   #      index=1
   #  if(config.parameters.get("lib-type")=="Unstranded") :


def  write_rmats_conf(project,projet_id_current,path,hash_fc,dataformatedforexport2Json_2,analysis,samples_by_project_test_or_control,sizeToclip,libType,single_or_paired,reads_for_rmats) :
    """
    Write Json configuration file for launching Rmats command for splicing.
  
    Args:
        project (dict): project.
        projet_id_current (str): projet_id_current.
        path (str): path.
        hash_fc (dict): hash_fc.
        dataformatedforexport2Json_2 (dict): dataformatedforexport2Json_2.
        samples_by_project_test_or_control (dict): samples_by_project_test_or_control.
        sizeToclip (str): sizeToclip.
        libType (str): libType.
        single_or_paired (str): single_or_paired.
        analysis (str): analysis.
        reads_for_rmats (str): reads_for_rmats.

    """  
    
    logger.info("\n-> Write RMATS JSON Configuration File \n")  
    #logger.info(hash_fc)  
    #logger.info(dataformatedforexport2Json_2)  
    #logger.info(analysis)  
    #logger.info(samples_by_project_test_or_control)  


    for projet_id in sorted(project.keys()) :


        if(projet_id==projet_id_current) :
            
            # Create RMats Ouput

            pathlib.Path(path+"output/"+projet_id+"/RMATS/").mkdir(parents=True, exist_ok=True) 

            logger.info("Config file is : ")
            outfile = open(path+projet_id+"_"+"rmats.json", 'w')  
            
            outfile.write(" {\n")
            outfile.write('"gtf" : "/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/gencode.v25.primary_assembly.annotation.gtf", \n')
            outfile.write('"project" : "'+projet_id+'", \n')
            outfile.write(' "path_to_output":"'+path+"output/"+projet_id+"/RMATS/"+'",\n')
            logger.info(path+"output/"+projet_id+"/RMATS/")
            outfile.write('"index" : "/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/indexes/GRCh38_PRIM_GENCODE_R25/" , \n')
            
            if ( (len(samples_by_project_test_or_control[projet_id_current]["CONTROL"]) > 1 ) and  (len(samples_by_project_test_or_control[projet_id_current]["TEST"]) > 1 ) ) : 
                outfile.write('"softPath" : "/home/jean-philippe.villemin/bin/anaconda3/envs/python2/bin/python2 /home/jean-philippe.villemin/bin/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py", \n')
            else : 
                outfile.write('"softPath" : "/home/jean-philippe.villemin/bin/anaconda3/envs/python2/bin/python2 /home/jean-philippe.villemin/bin/rMATS.3.2.5/RNASeq-MATS.py", \n')

            a=0
            outfile.write('"analysis" : {  \n')
            logger.info(hash_fc[projet_id_current])

            for comparision in hash_fc[projet_id_current] :

               outfile.write('"'+comparision+'" : {  \n')
               outfile.write('    "fastq" :  { \n')
               outfile.write('        "s1" : \n')
               outfile.write( str(reads_for_rmats[projet_id_current][0]).replace("'","\"")+' ,\n')
               outfile.write('        "s2" :  \n')
               outfile.write( str(reads_for_rmats[projet_id_current][1]).replace("'","\"")+'\n')
               outfile.write('}, \n')
               #logger.info(reads_for_rmats[projet_id_current][0])
               #logger.info(reads_for_rmats[projet_id_current][1])

               outfile.write('"len" : "'+str(sizeToclip)+'", \n')
               if(libType=="Forward") :
                    outfile.write('"libType" : "fr-firststrand", \n')
               if(libType=="Reverse") :
                   outfile.write('"libType" : "fr-secondstrand", \n')
               if(libType=="Unstranded") :
                    outfile.write('"libType" : "fr-unstranded", \n')
                    
               if(single_or_paired=="SINGLE") :
                   outfile.write('"type_read" : "single", \n')
               if(single_or_paired=="PAIRED") :
                   outfile.write('"type_read" : "paired", \n')
    
               outfile.write('"diff" : "0.1" \n')
               outfile.write('} \n')
               if(a!=len(hash_fc)-1) :  outfile.write(",")
               a+=1
            outfile.write('} \n')
            outfile.write('} \n')

           #"test" : "P", by default Unpaired and in the last version no different test to select
           #outfile.write('"diff" : "0.1", \n')

           #for  uniq_id_sample in project[projet_id] :
    outfile.close()
           
             
def write_main_conf(project,projet_id_current,path,hash_fc,dataformatedforexport2Json_2,analysis,samples_by_project_test_or_control,splicing_soft):
    """
    Write Json configuration file for applying Rmats / Whippet merge and filters actions.
  
    Args:
        project (dict): project.
        projet_id_current (str): projet_id_current.
        path_to_out (str): path.
        hash_fc (dict): hash_fc.
        dataformatedforexport2Json_2 (dict): dataformatedforexport2Json_2.
        samples_by_project_test_or_control (dict): samples_by_project_test_or_control.
        splicing_soft (str): splicing_soft.
  
    """  
    
    logger.info("\n-> Write Splicing "+splicing_soft+" JSON Configuration File \n")  
    
    for projet_id in sorted(project.keys()) :


        if(projet_id==projet_id_current) :
           logger.info("Config file is : ")
           outfile = open(path+projet_id+"_"+splicing_soft+"_"+"main.json", 'w')  

           outfile.write(" {\n")
           outfile.write('"organism" : "homo_sapiens", \n')
           outfile.write('"project" : "'+projet_id+'", \n')
           outfile.write(' "reads_junction":"10",\n')
           outfile.write('"project" : "'+projet_id+'" , \n')
           outfile.write('"sessionName" : "'+projet_id+'" , \n')
           outfile.write('"path_to_output" : "'+path+'", \n')
           outfile.write('"list_files_splicing" : ["SE"], \n')
           outfile.write('"gene_length": "/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_PRIM_GENCODE_R25/gencode.v25.gene.length.exon.sums.csv", \n')
           outfile.write('"read_count_matrice": "'+path+"output/"+projet_id+'/Raw_read_counts.csv",\n')
           outfile.write('"analysis": { \n')
           a=0
           for comparision in hash_fc[projet_id_current] :
               
               conds=comparision.split("_vs_")
               outfile.write('"'+comparision+'": { \n' )
               #"splicing_dir":  "/home/jean-philippe.villemin/data/data/PROJECT/WHOLE_EMT/paired", 
               if(splicing_soft=="RMATS") :
                   outfile.write('"soft":  "'+"RMATS"+'", \n')
                   outfile.write('"splicing_dir":  "'+path+"output/"+projet_id+"/RMATS/"+comparision+'/MATS_output/", \n')
               if (splicing_soft=="WHIPPET") :
                   outfile.write('"soft":  "'+"WHIPPET"+'", \n')
                   outfile.write('"splicing_dir":  "'+path+"output/"+projet_id+'/WHIPPET/", \n')

               outfile.write('"expression_file_path": "'+path+"output/"+projet_id+'", \n')
               outfile.write('"sample_test_for_quantification": "'+conds[0]+'",\n')
               outfile.write('"sample_control_for_quantification": "'+conds[1]+'",\n')
               outfile.write('"replicates_test": '+str(samples_by_project_test_or_control[projet_id]["TEST"]).replace("'",'"')+',\n')#.replace("{'","").replace("}'",""+',\n')
               outfile.write('"replicates_control":  '+str(samples_by_project_test_or_control[projet_id]["CONTROL"]).replace("'",'"')+'\n')#.replace("{'","").replace("}'","")
               outfile.write(' }')

               if(a!=len(hash_fc)-1) :  outfile.write(",")
               a+=1
           outfile.write(' }, \n')
           outfile.write(' "samples_for_quantification": \n')
           outfile.write(' { \n')

           count_condition = 0
           for condition in analysis[projet_id] : 
                count = 0
                outfile.write(' "'+condition+'": [  \n')
                logger.info(condition)
                for replicatLineQuantif in dataformatedforexport2Json_2[condition] :
                    logger.info(str(replicatLineQuantif))
                    outfile.write(str(replicatLineQuantif).replace("{'","").replace("}'",""+'\n'))
                    if(count!=len(dataformatedforexport2Json_2[condition])-1) : outfile.write(",")
                    count+=1
                outfile.write('  ] \n')
                if(count_condition!=len(analysis[projet_id])-1) : outfile.write(",\n")
                count_condition+=1
           outfile.write(' } \n')
           outfile.write(', \n')
           outfile.write('"tabs": \n')
           outfile.write(' { \n')
           i = 0
           for comparision in  hash_fc[projet_id_current] :
               outfile.write('"'+ comparision+'": { "analysisTocheck": ["'+comparision+'"],"0":"nevermind" }')
               if(i!=len(hash_fc)-1) :  outfile.write(",")
               i+=1
           outfile.write('  } \n')
           outfile.write('  } \n')

           #for  uniq_id_sample in project[projet_id] :
    outfile.close()

def  read_salmon_output_for_libtype(salmon_output_meta_file,typeEnd):
    """
    Read salmon output meta infos to get libtype
  
    Args:
        salmon_output_meta_file (str): filepath to STR.
        typeEnd (typeEnd): typeEnd.

    Returns:
        libType (str): Unstranded,Reverse,Forward
    
        
    """
    
    
    
    json1_file = open(salmon_output_meta_file)
    json1_str = json1_file.read()
    json1_data = json.loads(json1_str)
    
    logger.info(json1_data["expected_format"])
    
    if (typeEnd=="singleEnd") :
        if (json1_data["expected_format"]=="SR") :
            return "Reverse"
        if (json1_data["expected_format"]=="SF") :
            return "Forward"
        if (json1_data["expected_format"]=="U" or json1_data["expected_format"]=="IU" ) :
            return "Unstranded"
        
    if (typeEnd=="pairedEnd") :

        if (json1_data["expected_format"]=="ISR") :
            return "Reverse"
        if (json1_data["expected_format"]=="ISF") :
            return "Forward"
        if (json1_data["expected_format"]=="U" or json1_data["expected_format"]=="IU" ) :
            return "Unstranded"
    '''
                   Paired-end    Single-end
    -fr-unstranded    -l IU    -l U
    -fr-firststrand    -l ISR    -l SR
    -fr-secondstrand    -l ISF    -l SF  
    
    F = read 1 (or single-end read) comes from the forward strand
    R = read 1 (or single-end read) comes from the reverse strand
    '''   
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Thi script will launch all the steps for splicing and expression analysis for RNA-SEQ.  
    
    Example : 
    python3 ./signature.py -l listofRnaSeqProject.tsv

    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-l","--listing",action="store",help="File with all files to treat, read docs.",required=True,type=str,dest='path2Refs')

    parameters = parser.parse_args()
    
    path    = "/home/jean-philippe.villemin/data/data/PROJECT/WHOLE_EMT/"
    
    logger = create_logger(path,"INFO")
    
    logger.info("\n-> Log is written here : "+path)
    
    project = {}
    project_global = {}
    count   = 0
    
    with open(parameters.path2Refs) as f:
        for line in f:
            if(count==0) : 
                count=1 
                lineElements  = line.strip().split("\t")
                logger.info(lineElements)
                
                index_study                 = lineElements.index("STUDY")
                
                index_runaccession          = lineElements.index("RUNACCESSION")
                index_librarylayout         = lineElements.index("LIBRARY_LAYOUT")
                index_fastq                 = lineElements.index("FASTQ")
                index_kmer                  = lineElements.index("KMER")
    
                index_kmer                   = lineElements.index("KMER")
                
                index_condition              = lineElements.index("CONDITION")
                index_repNumber              = lineElements.index("REP_NUMBER")
                
                index_cellLine               = lineElements.index("CELL_LINE")
    
                index_treatment              = lineElements.index("TREATMENT")
                index_treatmentDay           = lineElements.index("TREATMENT_DAY")

                continue

            lineElements  = line.strip().split("\t")

            if lineElements[index_study] not in project : project[lineElements[index_study]]               = {}
            if lineElements[index_study] not in project_global : project_global[lineElements[index_study]] = {}
                
            project_global[lineElements[index_study]]["pair"] =  lineElements[index_librarylayout]    

            project[lineElements[index_study]][lineElements[index_runaccession]] = {
                     "pair"           : lineElements[index_librarylayout],    
                     "url"            : lineElements[index_fastq],
                     
                    # "read_count"              : lineElements[4],    
                    # "average_read_length"     : lineElements[7],
                    
                     "inducer"        : lineElements[index_treatment], # 
                     "cell_line"      : lineElements[index_cellLine],  # 
                     "condition"      : lineElements[index_condition], # CONTROL / TEST
                     "rep_number"     : lineElements[index_repNumber], # Can be Rep1, Rep2, Rep3, or just 1 ,2 ,3
                     "treatment_day"  : lineElements[index_treatmentDay], # Number 0.5, 1, 3
                     
                     "kmer"           : lineElements[index_kmer]

            }

        f.close()
    

    analysis                               = {}
    analysis_2_files                       = {}
    hash_fc                                = {}
    samples_by_project_test_or_control     = {}
    samples_by_project_conditions          = {}
    replicat_by_project                    = {}
    fastqNameAfterTrimGalore               = {}
    
    for projet_id in sorted(project.keys()) :
        
        dataformatedforexport2Json   = {}
        dataformatedforexport2Json_2 = {}
        # default
        kmer= "normal"
        
        logger.info("\n====> PROJECT : "+projet_id)
        # if only 2 samples...but you can have several comparison 3vs3, 1vs3 . 
        replicat_by_project[projet_id]   = len(project[projet_id])

        if projet_id not in samples_by_project_conditions                  :   samples_by_project_conditions[projet_id]                     = {}
        if projet_id not in samples_by_project_test_or_control             :   samples_by_project_test_or_control[projet_id]                = {}
        if "CONTROL" not  in samples_by_project_test_or_control[projet_id] : samples_by_project_test_or_control[projet_id]["CONTROL"]       = []  
        if "TEST" not  in samples_by_project_test_or_control[projet_id]    : samples_by_project_test_or_control[projet_id]["TEST"]          = []  
        if projet_id not in fastqNameAfterTrimGalore : fastqNameAfterTrimGalore[projet_id]         = {}

        for  uniq_id_sample in sorted(project[projet_id]) :
            samples_by_project_test_or_control[projet_id][project[projet_id][uniq_id_sample]["condition"]].append(projet_id+"."+uniq_id_sample+"."+project[projet_id][uniq_id_sample]["rep_number"])

            id_condition = ".".join( [projet_id,project[projet_id][uniq_id_sample]["cell_line"],project[projet_id][uniq_id_sample]["inducer"],project[projet_id][uniq_id_sample]["condition"],project[projet_id][uniq_id_sample]["treatment_day"]])
            logger.info(id_condition)
            kmer = project[projet_id][uniq_id_sample]["kmer"] 
             
            if id_condition not in samples_by_project_conditions[projet_id] : samples_by_project_conditions[projet_id][id_condition]  = [] 
            samples_by_project_conditions[projet_id][id_condition].append( uniq_id_sample+"."+project[projet_id][uniq_id_sample]["rep_number"])

            if id_condition not in dataformatedforexport2Json : dataformatedforexport2Json[id_condition]     = []
            if id_condition not in dataformatedforexport2Json_2 : dataformatedforexport2Json_2[id_condition] = []
            
            if id_condition not in fastqNameAfterTrimGalore[projet_id]   : fastqNameAfterTrimGalore[projet_id]  [id_condition]         = []

            url_fastq = project[projet_id][uniq_id_sample]["url"].split(";")

            url_1 = url_fastq[0].replace("\"","")
            url_2 = "None" if len(url_fastq)==1 else url_fastq[1].replace("\"","")
            
            if projet_id not in analysis_2_files : analysis_2_files[projet_id] = [] 
            analysis_2_files[projet_id].append(url_1)
            
            dataformatedforexport2Json_2[id_condition].append ( { '{"replicat_id": "'+  projet_id+"."+uniq_id_sample+"."+project[projet_id][uniq_id_sample]["rep_number"]+'","file_path":"'+path+"output/"+projet_id+"/"+projet_id+"."+uniq_id_sample+"."+project[projet_id][uniq_id_sample]["rep_number"]+"/"+'salmon_output/quant.genes.sf"}' } )
            fastqNameAfterTrimGalore[projet_id][id_condition].append (projet_id+"."+uniq_id_sample+"."+project[projet_id][uniq_id_sample]["rep_number"])
            
            if url_2 !=  "None" : analysis_2_files[projet_id].append(url_2)
            # set to none
            if "ftp" in url_fastq[0]:
                #path =  "None"
                # will use ftp DL so paste complete url
                dataformatedforexport2Json[id_condition].append ( {id_condition+"."+project[projet_id][uniq_id_sample]["rep_number"] :  { "R1" : url_1, "R2":url_2 } }  )

            else : 
                # remove path , path is then concatened from input
                dataformatedforexport2Json[id_condition].append ( {id_condition+"."+project[projet_id][uniq_id_sample]["rep_number"] :  { "R1" : os.path.basename(url_1), "R2":os.path.basename(url_2) } }  )

            write_align_conf(project,id_condition,projet_id,uniq_id_sample,path,path,project[projet_id][uniq_id_sample]["rep_number"])

            
            if projet_id not in analysis : 
                analysis[projet_id] =[]
                
            analysis[projet_id] .append(id_condition.strip())
        
        logger.info("\n-> Write Whippet JSON Configuration File \n")  

        with open(path+projet_id+"_"+"whippet.json", 'w') as outfile:  
            outfile.write("{")
            outfile.write('"project":"'+projet_id+'",\n')
            outfile.write('"organism":"human",\n')
            outfile.write('"path_to_output" : "'+path+"output/"+projet_id+"/WHIPPET/"'",\n')

            if "ftp" in url_fastq[0]:
                outfile.write('"path_to_input" : "'+"None"+'",\n')
            else :    
                outfile.write('"path_to_input" : "'+path+'",\n')

            outfile.write('"path_to_cleaner" : "/home/jean-philippe.villemin/code/RNA-SEQ/bash/whippet_filter_eventType_wrapped.sh",\n')
            outfile.write('"files": \n')
            json.dump(dataformatedforexport2Json, outfile)
            outfile.write(",\n")
            outfile.write('"analysis":{\n')                                  
            
           
            for projet_id_embed in analysis :
                if projet_id_embed == projet_id : 
                    #logger.info(projet_id)
                    
                    if projet_id not in hash_fc :  hash_fc[projet_id] = []
                    conditions = list(set(analysis[projet_id]))
                    list_of_pairs = [(conditions[p1], conditions[p2]) for p1 in range(len(conditions)) for p2 in range(p1+1,len(conditions))]
                    list_of_pairs = list(set(list_of_pairs))
                    logger.info("Pair of conditions :")
                    logger.info(list_of_pairs)
                    count_comparison = 0
                    for pair in list_of_pairs :
                        if ".TEST." in  pair[0] and ".TEST." not in  pair[1]: count_comparison+=1
                        if ".TEST." not in  pair[0] and ".TEST." in  pair[1]: count_comparison+=1

                    i = 0
                    #logger.info( count_comparison)
                    for pair in list_of_pairs :
    
                        if ".TEST." in  pair[0] and ".TEST." in  pair[1]: continue

                        if ".TEST." in  pair[0] :
                            #logger.info("TEST IS FIRST")
                            #logger.info(pair)
                            outfile.write('"'+pair[0]+"_vs_"+pair[1]+'": {\n')
                            outfile.write('"TEST":"'+pair[0]+'","CONTROL":"'+pair[1]+'"}\n')
                            hash_fc[projet_id].append(pair[0]+"_vs_"+pair[1])

                            
                        if ".TEST." in  pair[1] :
                            #logger.info("TEST IS SECOND")
                            #logger.info(pair)
                            
                            outfile.write('"'+pair[1]+"_vs_"+pair[0]+'": {\n')
                            hash_fc[projet_id].append(pair[1]+"_vs_"+pair[0])
                            outfile.write('"TEST":"'+pair[1]+'","CONTROL":"'+pair[0]+'"}\n')
                        if(i!=count_comparison-1) :  outfile.write(",\n")
                        i+=1
                   
                    outfile.write("}\n")
                    outfile.write("}")
        outfile.close()
        
        write_main_conf(project,projet_id,path,hash_fc,dataformatedforexport2Json_2,analysis,samples_by_project_test_or_control,"RMATS")
        write_main_conf(project,projet_id,path,hash_fc,dataformatedforexport2Json_2,analysis,samples_by_project_test_or_control,"WHIPPET")

        subprocess.run(("mkdir -p "+path+"logs/"+projet_id+"/trimmed"),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    
    
    reads_for_rmats = {}
    # Dowload and launch align     
    for projet_id in sorted(project.keys()) :
        logger.info("\nProject : "+projet_id+"\n")
        
        logger.info("\n-> Alternative Splicing Analyse with WHIPPET \n")

        subprocess.run(("mkdir -p "+path+"output/"+projet_id+"/WHIPPET/"),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

        whippet= "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/splicingWhippet.py -c "+path+projet_id+"_whippet.json -k "+kmer
        whippet_exec = subprocess.run((whippet),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        write_subprocess_log(whippet_exec,logger)
  
        proc_find = subprocess.getoutput(("find "+path+" -maxdepth 1 -name '"+projet_id+"*_align.json"+"'"))
        proc_find = proc_find.split("\n")
        logger.info(("find "+path+" -name -maxdepth 1 '"+projet_id+"*_align.json"+"'"))
        logger.info(proc_find)
        
        logger.info("Fastq files : ")
        for file in analysis_2_files[projet_id] :
            logger.info(os.path.basename(file) )
            
            filePath = path+os.path.basename(file)
            if not os.path.isfile(filePath):
            #prod_dl = subprocess.getoutput(("rsync -av rsync://"+file))s
                prod_dl= "wget "+file+" -P "+path
                prod_dl_exec = subprocess.run((prod_dl),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                write_subprocess_log(prod_dl_exec,logger)
                
        #outfile_average = open(path +projet_id+"_average_read_size.tsv", 'w') 
        logger.info("\n-> Alignment + Compute TPM values with Salmon \n")
        for oneSample in proc_find:

            m           = re.search('^\.', os.path.basename(oneSample))
            if m  : continue
            logger.info("Config file is : ")
            logger.info(os.path.basename(oneSample))

            name=os.path.basename(oneSample).replace("_align.json",".fastq.gz").split(".")
            nameFastq=name[1]


            mapping= "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/alignment.py -c "+oneSample#+"_align.json"
            logger.info("\n"+mapping+" \n")
            mapping_exec = subprocess.run((mapping),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            write_subprocess_log(mapping_exec,logger)

        end_type =""
        if (   project_global[projet_id]["pair"]  == "SINGLE") :
            end_type="singleEnd"
        if (   project_global[projet_id]["pair"]  == "PAIRED") : 
            end_type="pairedEnd"

        # GO READ SOMEWHERE read library in at least one file randomly choosen because they all should be the same for a same project 
        proc_find = subprocess.getoutput(("find "+path+"output/"+projet_id+" -maxdepth 3 -name lib_format_counts.json"))
        proc_find = proc_find.split("\n")
        logger.info(("\nfind "+path+"output/"+projet_id+" -maxdepth 3 -name lib_format_counts.json"))
        logger.info(proc_find[0])
        
        libType = read_salmon_output_for_libtype(proc_find[0],end_type)

        write_diff_exp_conf(project,projet_id,path,hash_fc,analysis,samples_by_project_conditions,libType)
    
        logger.info("\n-> Differential Expression \n")
        for comparison in hash_fc[projet_id] : 
            # Test vs Control
            logger.info("Comparison : "+comparison+"\n")
            
            conditions = comparison.split("_vs_")

            expressionAnalyse= "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/diffGeneExp.py -c "+path+projet_id+"_diff_exp.json -p "+conditions[0]+"_"+conditions[1]
            expressionAnalyse_exec = subprocess.run((expressionAnalyse),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            write_subprocess_log(expressionAnalyse_exec,logger)
            
            
            subprocess.run(("mkdir -p "+path+"output/"+projet_id+"/trimmed"),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

            coverageAnalyse= "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/coverage.py -c "+path+projet_id+"_RMATS_main.json -a "+conditions[0]+"_vs_"+conditions[1]+" "+"--filter=/home/jean-philippe.villemin/mount_archive2/commun.luco/ref/genes/GRCh38_ENSEMBL_85_HUMAN/genes.biomart.ensembl85.human.txt"
            logger.info("\n"+coverageAnalyse)
            coverageAnalyse_exec = subprocess.run((coverageAnalyse),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            write_subprocess_log(coverageAnalyse_exec,logger)

        list_size = []
        for file in analysis_2_files[projet_id] :
            logger.info(os.path.basename(file) )
            res = subprocess.getoutput(("/home/jean-philippe.villemin/code/RNA-SEQ/bash/read_average_size_fastq.sh "+file))
            logger.info("/home/jean-philippe.villemin/code/RNA-SEQ/bash/read_average_size_fastq.sh "+file)
            logger.info(res.split(" ")[0])
            list_size.append(int(res.split(" ")[0]))
            
        logger.info(list_size)
        sizeToclip = round(int(numpy.mean(list_size)))
        logger.info('Size to clip is : '+str(sizeToclip))

        cwd = os.getcwd()
       
        if projet_id not in reads_for_rmats : reads_for_rmats[projet_id] = {}
        
        avoid_more_add_trick =  {}
        
        logger.info("Nb files :  "+str(len(analysis_2_files[projet_id])))


        for indice in [0,1] :
            
            if indice not in reads_for_rmats[projet_id] : reads_for_rmats[projet_id][indice] = []
            
            logger.info("conditions[indice] : "+conditions[indice])
            
            # set to 0 
            if conditions[indice] not in avoid_more_add_trick : avoid_more_add_trick[conditions[indice]] = 0
                
            avoid_more_add_trick[conditions[indice]] +=1
            
            for file_trimgalored in fastqNameAfterTrimGalore[projet_id][conditions[indice]] :
                
                logger.info("Group files : "+file_trimgalored)
                # ATTTENTION ICI jp.data.Rep1
                if (project[projet_id][file_trimgalored.split(".")[1]]["pair"]=="PAIRED") :
                    
                    fileTrimmed_1=path+"output/"+projet_id+"/"+file_trimgalored+"/trim_galore_output/"+file_trimgalored.split(".")[1]+"_1_val_1.fq"
                    logger.info("\n-> Unzip "+fileTrimmed_1+".gz\n")
                    with gzip.open(fileTrimmed_1+".gz", 'rb') as f_in:
                        with open(fileTrimmed_1, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    logger.info("\n-> Trimming for Rmats R1 \n")
                    trimingForRmats = "/home/jean-philippe.villemin/bin/anaconda3/envs/python2/bin/python2 /home/jean-philippe.villemin/bin/rMATS.3.2.5/bin/trimFastq.py "+fileTrimmed_1 +" "+fileTrimmed_1.replace("_1_val_1.fq","_1.trimmed.clipped.fastq")+" "+str(sizeToclip)
                    logger.info(trimingForRmats)
                    trimingForRmats_exec = subprocess.run((trimingForRmats),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(trimingForRmats_exec,logger)
                    subprocess.run(("mv "+cwd+"/logtrimFastq*.txt "+path+"logs/"+projet_id+"/trimmed" ),shell=True)
                    # removed unziped
                    subprocess.run(("rm "+path+file_trimgalored+"_1.fastq"   ),shell=True) 

                    # Jdevrais créer un parseur pour ce log qui renvoie la mediane/moyenne

                    fileTrimmed_2=path+"output/"+projet_id+"/"+file_trimgalored+"/trim_galore_output/"+file_trimgalored.split(".")[1]+"_2_val_2.fq"
                    logger.info("\n-> Unzip  "+fileTrimmed_2+".gz\n")
                    with gzip.open(fileTrimmed_2+".gz", 'rb') as f_in:
                        with open(fileTrimmed_2, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    logger.info("\n-> Trimming for Rmats R2 \n")
                    trimingForRmats = "/home/jean-philippe.villemin/bin/anaconda3/envs/python2/bin/python2 /home/jean-philippe.villemin/bin/rMATS.3.2.5/bin/trimFastq.py "+fileTrimmed_2 +" "+fileTrimmed_2.replace("_2_val_2.fq","_2.trimmed.clipped.fastq")+" "+str(sizeToclip)
                    logger.info(trimingForRmats)
                    trimingForRmats_exec = subprocess.run((trimingForRmats),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(trimingForRmats_exec,logger)
                    subprocess.run(("mv "+cwd+"/logtrimFastq*.txt "+path+"logs/"+projet_id+"/trimmed"),shell=True)

                    if(avoid_more_add_trick[conditions[indice]] ==1) :  
                        reads_for_rmats[projet_id][indice].append(fileTrimmed_1.replace("_1_val_1.fq","_1.trimmed.clipped.fastq")+":"+fileTrimmed_2.replace("_2_val_2.fq","_2.trimmed.clipped.fastq"))
                    
                else : 
                    
                    fileTrimmed_1=path+"output/"+projet_id+"/"+file_trimgalored+"/trim_galore_output/"+file_trimgalored.split(".")[1]+"_1_trimmed.fq"
                    logger.info("unzip "+fileTrimmed_1+".gz\n")
                    with gzip.open(fileTrimmed_1+".gz", 'rb') as f_in:
                        with open(fileTrimmed_1, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    logger.info("\n-> Trimming for Rmats R1 \n")
                    trimingForRmats = "/home/jean-philippe.villemin/bin/anaconda3/envs/python2/bin/python2 /home/jean-philippe.villemin/bin/rMATS.3.2.5/bin/trimFastq.py "+fileTrimmed_1 +" "+fileTrimmed_1.replace("_1_trimmed.fq","_1.trimmed.clipped.fastq")+" "+str(sizeToclip)
                    logger.info(trimingForRmats)
                    trimingForRmats_exec = subprocess.run((trimingForRmats),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
                    write_subprocess_log(trimingForRmats_exec,logger)
                    subprocess.run(("mv "+cwd+"/logtrimFastq*.txt "+path+"logs/"+projet_id+"/trimmed" ),shell=True)
                    
                    if(avoid_more_add_trick[conditions[indice]] ==1) :  

                        reads_for_rmats[projet_id][indice].append(fileTrimmed_1.replace("_1_trimmed.fq","_1.trimmed.clipped.fastq"))

        logger.info("\n-> Alternative Splicing Analyse with Rmats \n")
         
        write_rmats_conf(project,projet_id,path,hash_fc,dataformatedforexport2Json_2,analysis,samples_by_project_test_or_control,sizeToclip,libType,project_global[projet_id]["pair"],reads_for_rmats)

        splicingRmats = "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/splicingRmats.py  -c "+path+projet_id+"_rmats.json"
        logger.info(splicingRmats)
        splicingRmats_exec = subprocess.run((splicingRmats),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        write_subprocess_log(splicingRmats_exec,logger)
        
        quit()
        
        logger.info("\n-> Wrap, Merge and Filter Splicing Analyse with RMATS \n")

        splicingSupperWrapper = "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/supperWrapperForMerge.py -c "+path+projet_id+"_RMATS_main.json"
        logger.info(splicingSupperWrapper)
        splicingSupperWrapper_exec = subprocess.run((splicingSupperWrapper),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        write_subprocess_log(splicingSupperWrapper_exec,logger)

        move ="mv "+path+"/"+comparison+"/ "+path+"output/"+projet_id+"/MERGE.RMATS."+comparison
        subprocess.run(move,shell=True)
        
        logger.info("\n-> Wrap, Merge and Filter Splicing Analyse with WHIPPET \n")

        splicingSupperWrapper = "python3 /home/jean-philippe.villemin/code/RNA-SEQ/src/supperWrapperForMerge.py -c "+path+projet_id+"_WHIPPET_main.json"
        logger.info(splicingSupperWrapper)
        splicingSupperWrapper_exec = subprocess.run((splicingSupperWrapper),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        write_subprocess_log(splicingSupperWrapper_exec,logger)
        
        move ="mv "+path+"/"+comparison+"/ "+path+"output/"+projet_id+"/MERGE.WHIPPET."+comparison
        subprocess.run(move,shell=True)
        
        subprocess.run(("mkdir -p "+path+"logs/"+projet_id),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        subprocess.run(("mkdir -p "+path+"configs/"+projet_id),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        
        subprocess.run(("mv "+path+"*.json "+path+"/configs/"+projet_id ),shell=True)
        subprocess.run(("mv "+path+"*.log "+path+"/logs/"+projet_id ),shell=True)


        logger.info(move)
