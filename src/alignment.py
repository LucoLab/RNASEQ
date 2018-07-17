
'''

:date: Oct, 2017
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Align the samples based on config set up in json format file.

'''
import argparse,textwrap
from utility import custom_parser
import subprocess
import logging
from logging.handlers import RotatingFileHandler
import datetime
import json

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
                

def core_sample_name(config,group_pair,opt_fastqtype):
    '''
    Clean up the names of the paired sample.
    Remove fastq.gz extension.
    Remove '_R1_' to have the same core syntax to use with '_R2_'.
  
    Args:
        config (obj): Configuration instance.
        group_pair (str): The pair identifier of fastq sequence to trait.
        opt_fastqtype (str) : say if remove fastq.gz ou juste fastq
    Returns:
        clean_name (string) : Clean core name of sample without R1 / R2 disctinction.
    
    
    '''
    clean_name =  config.parameters["files"][group_pair]["R1"].replace(opt_fastqtype,"",1)
    clean_name =  clean_name.replace("_R1","",1)
    clean_name =  clean_name.replace("_1","",1)

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
    file_handler = RotatingFileHandler(config.parameters['path_to_output']+'alignment_activity.log', 'a', 1000000, 1)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    ''' 2 nd handler in stream'''
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    logger.addHandler(stream_handler)

    return logger

def  read_salmon_output_for_strandness(salmon_output_meta_file):
    
    json1_file = open(salmon_output_meta_file)
    json1_str = json1_file.read()
    json1_data = json.loads(json1_str)
    
    #logger.info(json1_data["expected_format"])
    
    if (json1_data["expected_format"]=="U" or json1_data["expected_format"]=="IU" ) :
        return "Unstranded"
    else :  return "Stranded"

def  read_salmon_output_for_libtype(salmon_output_meta_file,typeEnd):
    
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

def print_count_reads_per_gene(list_core_sample_names,config,stranded_or_not,type):
    '''
    Pull all the intermediates files _ReadsPerGene.
    Compute a final reads count per gene.
  
    Args:
        list_core_sample_names (list): List of the sample names.
        config (object): Config parameters from json file.
    '''
    read_count         = {"str0":{},"str1":{},"str2":{}}
    read_stat         = {"N_unmapped":{"str0":0,"str1":0,"str2":0},"N_multimapping":{"str0":0,"str1":0,"str2":0},"N_noFeature":{"str0":0,"str1":0,"str2":0},"N_ambiguous":{"str0":0,"str1":0,"str2":0}}
    
    '''For stranded=no, a read is considered overlapping with a feature regardless of whether
     it is mapped to the same or the opposite strand as the feature. For stranded=yes and single-end reads, 
     the read has to be mapped to the same strand as the feature. 
     For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. 
     For stranded=reverse, these rules are reversed
    '''
    
    logger.info("Parse STAR counts output.\n")

    for core_sample_name in list_core_sample_names :
        logger.info(core_sample_name)
        with open(config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+core_sample_name+"_ReadsPerGene.out.tab", 'r') as f:
            for line in f:
                cleanedLine = line.strip()
                if cleanedLine: # is not empty
                    lineElements = cleanedLine.split()
                    gene_id  = lineElements[0]
                    if gene_id in ("N_unmapped","N_ambiguous","N_noFeature","N_multimapping") : 
                        read_stat[gene_id]['str0']+= int(lineElements[1])
                        read_stat[gene_id]['str1']+= int(lineElements[2])
                        read_stat[gene_id]['str2']+= int(lineElements[3])
                        continue
                        
                    count_r0 = lineElements[1]
                    count_r1 = lineElements[2]
                    count_r2 = lineElements[3]
                    #if gene_id == "ENSG00000227232.5" : print(cleanedLine)
                    if not read_count['str0'].get(gene_id) : 
                        read_count['str0'][gene_id] = 0
                    if not read_count['str1'].get(gene_id) : 
                        read_count['str1'][gene_id] = 0
                    if not read_count['str2'].get(gene_id) : 
                        read_count['str2'][gene_id] = 0
                    read_count['str1'][gene_id] += int(count_r1)
                    read_count['str2'][gene_id] += int(count_r2)
                    read_count['str0'][gene_id] += int(count_r0)

        f.close()
    count = {}
    
    logger.info("Write .tab used after for diff Exp : "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+core_sample_name+"_ReadsPerGene.tab")
    
    file_Reads_PerGene = open(config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+core_sample_name+"_ReadsPerGene.tab", "w")

    #file_Reads_PerGene.write("Gene\tUnstranded\tBadPairs(F1R2+/F2R1-)\tGoodPairs(F2R1+/F1R2-)"+"\n")
    file_Reads_PerGene.write("Gene\tUnstranded\tForward(F1R2+/F2R1-)\tReverse(F2R1+/F1R2-)"+"\n")

    for key1,value1 in read_count['str1'].items() :
        file_Reads_PerGene.write(str(key1)+"\t"+str(read_count['str0'][key1])+"\t"+str(value1)+"\t"+str(read_count['str2'][key1])+"\n")
        
        if 'Unstranded' not in count : count['Unstranded']              = 0 
        if  'Stranded_forward' not in count : count['Stranded_forward'] = 0 
        if 'Stranded_reverse'  not in count : count['Stranded_reverse'] = 0
        
        count['Unstranded']+=read_count['str0'][key1]
        count['Stranded_forward']+=value1
        count['Stranded_reverse']+=read_count['str2'][key1]
        
    file_Reads_PerGene.close()

    logger.info(count)
    #strand_orientation="Unstranded"
#    if (count['Stranded_reverse'] > count['Stranded_forward'] ) :
    #   strand_orientation = "Reverse"
    #else : strand_orientation = "Forward"
        
    #return strand_orientation

def read_log_final(list_core_sample_names,config,stranded_or_not):
    '''
    Pull all the intermediates files _ReadsPerGene.
    Compute a final reads count per gene.
  
    Args:
        list_core_sample_names (list): List of the sample names.
        config (object): Config parameters from json file.
    '''

    stats ={}
    for core_sample_name in list_core_sample_names :
        stats[core_sample_name] = {}
        with open(config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+core_sample_name+"_Log.final.out", 'r') as f:
            for line in f:
                cleanedLine = line.strip().replace("\t","")
                if cleanedLine: # is not empty
                    if '|' in cleanedLine : 
                        #logger.info(cleanedLine)
                        lineElements = cleanedLine.split('|')
    
                        stats[core_sample_name][lineElements[0]]= lineElements[1]
                        if(lineElements[0]=="Uniquely mapped reads %") :
                            if(float(lineElements[1]) < 50 ) : 
                                logger.info (">>>>>>>>>>>>>>>>>>>>>> BAM CONTAINS LESS THAN 50 % READS MAPPED <<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
                        if(lineElements[0]=="Average input read length") :
                                logger.info ("Average read size is "+lineElements[1])
                                                        
        f.close()
        
    return stats   
'''
                                 Started job on |    Mar 15 11:56:27
                             Started mapping on |    Mar 15 12:02:30
                                    Finished on |    Mar 15 12:06:19
       Mapping speed, Million of reads per hour |    316.72

                          Number of input reads |    20146908
                      Average input read length |    250
                                    UNIQUE READS:
                   Uniquely mapped reads number |    16728164
                        Uniquely mapped reads % |    83.03%
                          Average mapped length |    244.91
                       Number of splices: Total |    17272633
            Number of splices: Annotated (sjdb) |    17127586
                       Number of splices: GT/AG |    17105007
                       Number of splices: GC/AG |    126140
                       Number of splices: AT/AC |    13125
               Number of splices: Non-canonical |    28361
                      Mismatch rate per base, % |    0.22%
                         Deletion rate per base |    0.01%
                        Deletion average length |    1.87
                        Insertion rate per base |    0.01%
                       Insertion average length |    1.46
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |    616055
             % of reads mapped to multiple loci |    3.06%
        Number of reads mapped to too many loci |    6838
             % of reads mapped to too many loci |    0.03%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |    0.00%
                 % of reads unmapped: too short |    13.85%
                     % of reads unmapped: other |    0.02%
                                  CHIMERIC READS:
                       Number of chimeric reads |    0
                            % of chimeric reads |    0.00%
'''
###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    -- Align the fastq.gz files per directory/condition. -- 
    
    Path , tools, memory and number of threads are given in JSON file.
    It describes where files are installed.
    One JSON file is needed per condition/directory.
    You treat each condition separately.
    
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-c","--config",action="store",help="Path to a json file for the condition.",required=True,type=str,dest='file_config')
    parser.add_argument("-a","--aligner",action="store",help="Star is default setting.",default="star",type=str,dest='aligner')

    parameters = parser.parse_args()

    config = custom_parser.Configuration(parameters.file_config,"json")

    #print(parameters.file_config)
    outputdirname="output/"+config.parameters['project']

    trim_galore_output=config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"trim_galore_output"
    salmon_output=config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"salmon_output"

    subprocess.run(("mkdir -p "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R1"),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    subprocess.run(("mkdir -p "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R2"),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    subprocess.run(("mkdir -p "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"trim_galore_output"),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

    print("mkdir -p "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R1")
    print("mkdir -p "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R2")
    print("mkdir -p "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"trim_galore_output")

    logger = create_logger(config)
    #logger.info(config.chrono)

    logger.info('=========> START:PARSE EACH FASTQ PAIRS FOR AN EXPERIMENT')
   
    list_core_sample_names_for_star  = []
    list_core_sample_names           = []
    
    opt_fastqtype=".fastq.gz"
    opt_gzippedFile=" --readFilesCommand zcat "

    if (config.parameters["gzip"] =='no') :
        opt_gzippedFile=" "
        opt_fastqtype=".fastq"
        
    ##############################################
    #### BE AWARE NOW ONLY R1 and R2 ARE HANDLED
    ##############################################
    

    
    for group_pair in config.parameters.get("files").keys():
       

        logger.info(group_pair)
        logger.info (config.parameters["files"][group_pair]["R1"])
        logger.info (config.parameters["files"][group_pair]["R2"])
        logger.info ("")

        clean_sample_name = core_sample_name(config,group_pair,opt_fastqtype)
        #list_core_sample_names.append(clean_sample_name)
        logger.info(clean_sample_name)
        logger.info(config.parameters['final_bam_name'])

        list_core_sample_names.append(config.parameters['final_bam_name'])
        
       
        ###########################################################################################################
        ###########################################################################################################
        #default
        logger.info('=========> FASTQC PER FASTQ')
        if (config.parameters["type"] =='singleEnd') :
            logger.info('=========>Fastqc SingleEND')
       
            logger.info((
                "fastqc --extract --threads " \
                +config.parameters["fastqc"]["threads"] \
                +" --outdir "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R1 " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']))
            
            proc_fastqc_single = subprocess.run((
                "fastqc --extract --threads " \
                +config.parameters["fastqc"]["threads"] \
                +" --outdir "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R1  " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1'])
                 ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            write_subprocess_log(proc_fastqc_single,logger)

            logger.info((
                "trim_galore --fastqc --output_dir "+trim_galore_output+" " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']))
            
            proc_fastqc_single_trim = subprocess.run((
               "trim_galore --fastqc --output_dir "+trim_galore_output+" " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1'])
                 ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            write_subprocess_log(proc_fastqc_single_trim,logger)


        if (config.parameters["type"] =='pairedEnd') :
            logger.info('=========>Fastqc PairedEND')
       
            logger.info((
                "fastqc --extract --threads " \
                +config.parameters["fastqc"]["threads"] \
                +" --outdir "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R1 " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']))
            
            proc_fastqc_paired_R1 = subprocess.run((
                "fastqc --extract --threads " \
                +config.parameters["fastqc"]["threads"] \
                +" --outdir "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R1 " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1'])
                 ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            write_subprocess_log(proc_fastqc_paired_R1,logger)

            logger.info((
                "fastqc --extract --threads " \
                +config.parameters["fastqc"]["threads"] \
                +" --outdir "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R2 " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2']))
            
            proc_fastqc_paired_R2= subprocess.run((
                "fastqc --extract --threads " \
                +config.parameters["fastqc"]["threads"] \
                +" --outdir "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/R2 " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2'])
                 ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            write_subprocess_log(proc_fastqc_paired_R2,logger)
            
            logger.info((
                "trim_galore --paired --fastqc --output_dir "+trim_galore_output+" " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']+' '+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2']))
            
            proc_fastqc_paired_trim = subprocess.run((
               "trim_galore --paired --fastqc --output_dir "+trim_galore_output+" " \
                +config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']+' '+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2'])
                 ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            write_subprocess_log(proc_fastqc_paired_trim,logger)

        ####################################################################################################
        ##########     USE SALMON TO GET LIBRARY TYPE UNSTRAND / STRANDED #########################################
        ###########################################################################################################
        if (config.parameters["type"] =='singleEnd') :

            logger.info('=========> SALMON QUANT')

            logger.info("salmon quant -i "+config.parameters['transcriptome_index'] \
                +" -p " +config.parameters['salmon']['cpu'] \
                +" -l A " \
                +" -g "+config.parameters['path_to_gtf'] \
                +" -r "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R1'].replace('.fastq','_trimmed.fq')+" " \
                +" -o "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"salmon_output")
             
            salmon= subprocess.run("salmon quant -i "+config.parameters['transcriptome_index'] \
                +" -p " +config.parameters['salmon']['cpu'] \
                +" -l A " \
                +" -g "+config.parameters['path_to_gtf']  \
                +" -r "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R1'].replace('.fastq','_trimmed.fq')+" " \
                +" -o "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"salmon_output"
            ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)            
                   
            write_subprocess_log(salmon,logger)

        if (config.parameters["type"] == "pairedEnd") :

            logger.info('=========> SALMON QUANT')
    
        ####################### Salmonize ##########################
            logger.info("salmon quant -i "+config.parameters['transcriptome_index'] \
                +" -p " +config.parameters['salmon']['cpu'] \
                +" -l A " \
                +" -g "+config.parameters['path_to_gtf'] \
                +" -1 "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R1'].replace('.fastq','_val_1.fq') \
                +" -2 "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R2'].replace('.fastq','_val_2.fq') \
                +" -o "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"salmon_output")
             
            salmon= subprocess.run("salmon quant -i "+config.parameters['transcriptome_index'] \
                +" -p " +config.parameters['salmon']['cpu'] \
                +" -l A " \
                +" -g "+config.parameters['path_to_gtf']  \
                +" -1 "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R1'].replace('.fastq','_val_1.fq') \
                +" -2 "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R2'].replace('.fastq','_val_2.fq') \
                +" -o "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"salmon_output"
            ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            
            write_subprocess_log(salmon,logger)


        if parameters.aligner == "star" :
            logger.info('=========> ALIGN WITH STAR')
            
            #https://github.com/alexdobin/STAR/issues/143
            #Default pairend' abc'.translate(str.maketrans('ac','xy'))
            opt_input=" --readFilesIn "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R1'].replace('.fastq','_val_1.fq')+" "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R2'].replace('.fastq','_val_2.fq')+" "
            #opt_input=" --readFilesIn <(zcat "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']+") <(zcat  "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2']+") "

            if (config.parameters["type"] == "singleEnd") :
                logger.info('=========> SINGLE END')

                opt_input=" --readFilesIn "+trim_galore_output+'/'+config.parameters['files'][group_pair]['R1'].replace('.fastq','_trimmed.fq')+" "
                #opt_input=" --readFilesIn <(zcat "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']+") "

            # --sjdbFileChrStartEnd For the second pass
            logger.info(config.parameters['softs']['star'] \
                +" --runThreadN "+config.parameters['star']['runThreadN']         \
                +" --genomeDir "+config.parameters['star']['genomeDir']         \
                +opt_input         \
                +opt_gzippedFile      \
                +" --outFileNamePrefix "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_"        \
                +" --outSAMtype "+config.parameters['star']['outSAMtype']      \
                +" --quantMode GeneCounts "        \
                #+" --outSAMattrIHstart "+config.parameters['star']['outSAMattrIHstart']        \
                #+" --outWigType bedGraph "        \
                #+" --outWigStrand Stranded"  \
                #+" --outWigNorm None"
                )         
            #TESTAligned.sortedByCoord.out.bam
            proc_star = subprocess.run((config.parameters['softs']['star'] \
                +" --runThreadN "+config.parameters['star']['runThreadN']         \
                +" --genomeDir "+config.parameters['star']['genomeDir']         \
                +opt_input         \
                +opt_gzippedFile      \
                +" --outFileNamePrefix "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_"        \
                +" --outSAMtype "+config.parameters['star']['outSAMtype']      \
                +" --quantMode GeneCounts"        \
               # +" --outSAMattrIHstart "+config.parameters['star']['outSAMattrIHstart']        \
                #+" --outWigType bedGraph "       \
                #+" --outWigStrand Stranded" \
                #+" --outWigNorm None"
                ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            #nohup python3 ./align.py -a=star -nfqc -c ../configs/GHRC38/MANTrep1.json &
            list_core_sample_names_for_star.append( config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Aligned.sortedByCoord.out.bam")

            write_subprocess_log(proc_star,logger)    
    '''               
    logger.info('=========> MERGE DIFFERENTS BAM')
    
  
    list_core_sample_names_final = []
    if parameters.aligner == "star" :
        list_core_sample_names_final = list_core_sample_names_for_star

    logger.info(
        config.parameters['softs']['samtools']+" merge -r -@ " \
        +config.parameters['samtools']['cpu']+" " \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+".bam"+" " \
        +(" ".join(list_core_sample_names_final))  )

    proc_sam2 = subprocess.run((
        config.parameters['softs']['samtools']+" merge -r -@ " \
        +config.parameters['samtools']['cpu'] \
        +" "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+".bam"+" " \
        +(" ".join(list_core_sample_names_final))  ),shell=True)
    
    write_subprocess_log(proc_sam2,logger)
    
    logger.info('=========> SORT BAM')


    proc_sam4 =  subprocess.run((
        config.parameters['softs']['samtools']+" sort -@ 12 -o "
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/SORTED."+config.parameters['final_bam_name']+".bam"+" " \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+".bam") ,shell=True)
    
    write_subprocess_log(proc_sam4,logger)
    '''
    
     #list_core_sample_names_for_star.append( config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Aligned.sortedByCoord.out.bam")


    subprocess.run(("mv "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Aligned.sortedByCoord.out.bam "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/SORTED."+config.parameters['final_bam_name']+".bam"),shell=True)

    logger.info('=========> INDEX BAM')

    proc_sam3 =  subprocess.run((
        config.parameters['softs']['samtools']+" index " \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/SORTED."+config.parameters['final_bam_name']+".bam") ,shell=True)   
    
    write_subprocess_log(proc_sam3,logger)

    ###############################################################################################
    ###################################       FASTQC      ######################################### 
    ###############################################################################################

    logger.info("=========> COUNT READS FOR GENE") 

    # You want to know strandness and libtype (forward/reverse)

    stranded_or_not    = read_salmon_output_for_strandness(salmon_output+"/lib_format_counts.json")
    reverse_or_forward_or_unstrand = read_salmon_output_for_libtype(salmon_output+"/lib_format_counts.json",config.parameters["type"])
    logger.info(">>>>>>>>>>>>>>>>>>> From Salmon, Library is "+reverse_or_forward_or_unstrand)

   
     # You want to know if it's forward ou reverse when Stranded using count from STAR
    print_count_reads_per_gene(list_core_sample_names,config,stranded_or_not,config.parameters["type"])
    #strand_orientation = 
    
    logger.info(stranded_or_not)
    logger.info(reverse_or_forward_or_unstrand)
    
    #if(reverse_or_forward_or_unstrand != strand_orientation) :
        #logger.info("Salmon and STAR define different library types.")

    stats = read_log_final(list_core_sample_names,config,stranded_or_not)
    #logger.info(stats)
   
    name=config.parameters['final_bam_name'].split(".")#Projet.ID_SAMPLE.REP#
    logger.info(name)
    nameFastq=name[1]
  
    
    logger.info("=========> STAR GIVE WIG FILES ") 

    logger.info(config.parameters['softs']['star'] \
    +" --runMode inputAlignmentsFromBAM"    \
    +" --runThreadN "+config.parameters['star']['runThreadN']         \
    +" --inputBAMfile "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/SORTED."+config.parameters['final_bam_name']+".bam"        \
    +" --outFileNamePrefix "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_"        \
    #+" --quantMode GeneCounts "        \
    +" --outWigType wiggle "        \
    +" --outWigStrand "+stranded_or_not  \
    +" --outWigNorm RPM"
    )         
    proc_star_for_wig = subprocess.run((config.parameters['softs']['star'] \
    +" --runMode inputAlignmentsFromBAM" \
    +" --runThreadN "+config.parameters['star']['runThreadN']         \
    +" --inputBAMfile "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/SORTED."+config.parameters['final_bam_name']+".bam"         \
    +" --outFileNamePrefix "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_"        \
    #+" --quantMode GeneCounts "        \
    +" --outWigType wiggle "        \
    +" --outWigStrand "+stranded_or_not  \
    +" --outWigNorm RPM"
    ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    
    write_subprocess_log(proc_star_for_wig,logger)

    if (stranded_or_not=="Stranded") :

        logger.info("=========> WIG GOES BIGWIG") 


        proc_bigwiggle_1 = subprocess.run((
        config.parameters['softs']['wigToBigWig'] \
        +" -clip "    \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.Unique.str1.out.wig "        \
        +config.parameters['path_to_chrom_length_ref']+" "      \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.Unique.negStrand.out.bw "  
        ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        
        write_subprocess_log(proc_bigwiggle_1,logger)
    
        
        proc_bigwiggle_2 = subprocess.run((
        config.parameters['softs']['wigToBigWig'] \
        +" -clip "    \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.Unique.str2.out.wig "        \
        +config.parameters['path_to_chrom_length_ref']+" "      \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.Unique.posStrand.out.bw "  
        ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        
        write_subprocess_log(proc_bigwiggle_2,logger)
    
        
        proc_bigwiggle_3 = subprocess.run((
        config.parameters['softs']['wigToBigWig'] \
        +" -clip "    \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.str1.out.wig "        \
        +config.parameters['path_to_chrom_length_ref']+" "      \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.negStrand.out.bw "  
        ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    
        write_subprocess_log(proc_bigwiggle_3,logger)
        
        proc_bigwiggle_4 = subprocess.run((
        config.parameters['softs']['wigToBigWig'] \
        +" -clip "    \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.str2.out.wig "        \
        +config.parameters['path_to_chrom_length_ref']+" "      \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.posStrand.out.bw "  
        ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    
        write_subprocess_log(proc_bigwiggle_4,logger)

    if (stranded_or_not=="Unstranded") :

        proc_bigwiggle_5 = subprocess.run((
        config.parameters['softs']['wigToBigWig'] \
        +" -clip "    \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.Unique.str1.out.wig "        \
        +config.parameters['path_to_chrom_length_ref']+" "      \
        +config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+"_Signal.Unique.unStrand.out.bw "  
        ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
        
        write_subprocess_log(proc_bigwiggle_5,logger)

    
    #logger.info('=========> REMOVE INTERMEDIATE FILES')
    #subprocess.run(("mv "+config.parameters['path_to_input']+"*_trimming_report.txt "+config.parameters['path_to_output']+outputdirname+"/" ),shell=True)
    #subprocess.run(("mv "+config.parameters['path_to_input']+"*.html "+config.parameters['path_to_output']+outputdirname+"/" ),shell=True)
    #subprocess.run(("mv "+config.parameters['path_to_input']+"*.zip "+config.parameters['path_to_output']+outputdirname+"/" ),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+".bam" ),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*_Aligned.sortedByCoord.out.bam" ),shell=True)
    subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*.wig" ),shell=True)
    subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*.bam" ),shell=True)
    subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*.bam.bai" ),shell=True)

    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+config.chrono+"/*_ReadsPerGene.out.tab"),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Log.out" ),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Log.final.out" ),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Log.progress.out" ),shell=True)
    #subprocess.run(("cat "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*SJ.out.tab > "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/"+config.parameters['final_bam_name']+".All.SplicingJunctions.out.tab"),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*SJ.out.tab" ),shell=True)
    #subprocess.run(("rm  "+config.parameters['path_to_output']+outputdirname+"/"+config.parameters['final_bam_name']+"/"+"star_output"+"/*.zip"),shell=True)

    timeline = datetime.datetime.now()
    endtime = str(timeline.hour)+"_"+str(timeline.minute)+"_"+str(timeline.second)

  
    logger.info('=========> END '+endtime)
    
    
    

