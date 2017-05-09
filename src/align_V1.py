
'''

:date: Sep 6, 2016
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

def print_count_reads_per_gene(list_core_sample_names,config):
    '''
    Pull all the intermediates files _ReadsPerGene.
    Compute a final reads count per gene.
  
    Args:
        list_core_sample_names (list): List of the sample names.
        config (object): Config parameters from json file.
    '''
    read_count         = {"str1":{},"str2":{}}
    file_Reads_PerGene = open(config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_ReadsPerGene.tab", "w")
   
    for core_sample_name in list_core_sample_names :
       
        with open(config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+core_sample_name+"_ReadsPerGene.out.tab", 'r') as f:
            for line in f:
                cleanedLine = line.strip()
                if cleanedLine: # is not empty
                    lineElements = cleanedLine.split()
                    gene_id  = lineElements[0]
                    if gene_id in ("N_unmapped","N_ambiguous","N_noFeature","N_multimapping") : continue
                    count_r1 = lineElements[2]
                    count_r2 = lineElements[3]
                    #if gene_id == "ENSG00000227232.5" : print(cleanedLine)
                    if not read_count['str1'].get(gene_id) : 
                        read_count['str1'][gene_id] = 0
                    if not read_count['str2'].get(gene_id) : 
                        read_count['str2'][gene_id] = 0
                    read_count['str1'][gene_id] += int(count_r1)
                    read_count['str2'][gene_id] += int(count_r2)
   
        f.closed
        
    file_Reads_PerGene.write("Gene\tBadPairs(F1R2+/F2R1-)\tGoodPairs(F2R1+/F1R2-)"+"\n")
    for key1,value1 in read_count['str1'].items() :
        file_Reads_PerGene.write(str(key1)+"\t"+str(value1)+"\t"+str(read_count['str2'][key1])+"\n")
    file_Reads_PerGene.close()

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
    #print(parameters.file_config)

    config = custom_parser.Configuration(parameters.file_config,"json")

    subprocess.run(("mkdir -p "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    print("mkdir -p "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono)
    
    logger = create_logger(config)
    logger.info(config.chrono)

    logger.info('=========> START:PARSE EACH FASTQ PAIRS FOR AN EXPERIMENT')
   
    list_core_sample_names_for_star  = []
    list_core_sample_names           = []
    

    for group_pair in config.parameters.get("files").keys():
       
        logger.info(group_pair)
        logger.info (config.parameters["files"][group_pair]["R1"])
        logger.info (config.parameters["files"][group_pair]["R2"])
        logger.info ("")

        clean_sample_name = core_sample_name(config,group_pair)
        list_core_sample_names.append(clean_sample_name)
  
        ###########################################################################################################
        ###########################################################################################################

        if parameters.aligner == "star" :
            logger.info('=========> ALIGN WITH STAR')
            
            # --sjdbFileChrStartEnd For the second pass
            logger.info(config.parameters['softs']['star'] \
                +" --runThreadN "+config.parameters['star']['runThreadN']         \
                +" --genomeDir "+config.parameters['star']['genomeDir']         \
                +" --readFilesIn "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']+" "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2']         \
                +" --readFilesCommand zcat "      \
                +" --outFileNamePrefix "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+clean_sample_name+"_"        \
                +" --outSAMtype "+config.parameters['star']['outSAMtype']      \
                +" --quantMode GeneCounts "        \
                +" --outSAMattrIHstart "+config.parameters['star']['outSAMattrIHstart']        \
                #+" --outWigType bedGraph "        \
                #+" --outWigStrand Stranded"  \
                #+" --outWigNorm None"
                )         
            #TESTAligned.sortedByCoord.out.bam
            proc_star = subprocess.run((config.parameters['softs']['star'] \
                +" --runThreadN "+config.parameters['star']['runThreadN']         \
                +" --genomeDir "+config.parameters['star']['genomeDir']         \
                +" --readFilesIn "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R1']+" "+config.parameters['path_to_input']+''+config.parameters['files'][group_pair]['R2']         \
                +" --readFilesCommand zcat"      \
                +" --outFileNamePrefix "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+clean_sample_name+"_"        \
                +" --outSAMtype "+config.parameters['star']['outSAMtype']      \
                +" --quantMode GeneCounts"        \
                +" --outSAMattrIHstart "+config.parameters['star']['outSAMattrIHstart']        \
                #+" --outWigType bedGraph "       \
                #+" --outWigStrand Stranded" \
                #+" --outWigNorm None"
                ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
            #nohup python3 ./align.py -a=star -nfqc -c ../configs/GHRC38/MANTrep1.json &
            list_core_sample_names_for_star.append( config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+clean_sample_name+"_Aligned.sortedByCoord.out.bam")

            write_subprocess_log(proc_star,logger)    
                   
    logger.info('=========> MERGE DIFFERENTS BAM')
    
    list_core_sample_names_final = []
    if parameters.aligner == "star" :
        list_core_sample_names_final = list_core_sample_names_for_star

    logger.info(
        config.parameters['softs']['samtools']+" merge -r -@ " \
        +config.parameters['samtools']['cpu']+" " \
        +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam"+" " \
        +(" ".join(list_core_sample_names_final))  )

    proc_sam2 = subprocess.run((
        config.parameters['softs']['samtools']+" merge -r -@ " \
        +config.parameters['samtools']['cpu'] \
        +" "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam"+" " \
        +(" ".join(list_core_sample_names_final))  ),shell=True)
    
    write_subprocess_log(proc_sam2,logger)
    
    proc_sam3 =  subprocess.run((
        config.parameters['softs']['samtools']+" index " \
        +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam") ,shell=True)
    
    write_subprocess_log(proc_sam3,logger)

    ###############################################################################################
    ###################################       FASTQC      ######################################### 
    ###############################################################################################
    
    logger.info('=========> FASTQC FINAL BAM')

    logger.info((
        "fastqc --extract --threads " \
        +config.parameters["fastqc"]["threads"] \
        +" --outdir "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+" " \
        +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam"))
    
    proc_fastqc_final = subprocess.run((
        "fastqc --extract --threads " \
        +config.parameters["fastqc"]["threads"] \
        +" --outdir "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+" " \
        +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam")
        ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

    write_subprocess_log(proc_fastqc_final,logger)
    
    logger.info("=========> STAR GIVE WIG FILES STRANDED") 

    logger.info(config.parameters['softs']['star'] \
    +" --runMode inputAlignmentsFromBAM"    \
    +" --runThreadN "+config.parameters['star']['runThreadN']         \
    +" --inputBAMfile "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam"        \
    +" --outFileNamePrefix "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_"        \
    #+" --quantMode GeneCounts "        \
    +" --outWigType wiggle "        \
    +" --outWigStrand Stranded"  \
    +" --outWigNorm RPM"
    )         
    proc_star_for_wig = subprocess.run((config.parameters['softs']['star'] \
    +" --runMode inputAlignmentsFromBAM" \
    +" --runThreadN "+config.parameters['star']['runThreadN']         \
    +" --inputBAMfile "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".bam"         \
    +" --outFileNamePrefix "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_"        \
    #+" --quantMode GeneCounts "        \
    +" --outWigType wiggle "        \
    +" --outWigStrand Stranded"  \
    +" --outWigNorm RPM"
    ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

    write_subprocess_log(proc_star_for_wig,logger)

    logger.info("=========> WIG GOES BIGWIG") 

    proc_bigwiggle_1 = subprocess.run((
    config.parameters['softs']['wigToBigWig'] \
    +" -clip "    \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.Unique.str1.out.wig "        \
    +config.parameters['path_to_chrom_length_ref']+" "      \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.Unique.negStrand.out.bw "  
    ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    
    proc_bigwiggle_2 = subprocess.run((
    config.parameters['softs']['wigToBigWig'] \
    +" -clip "    \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.Unique.str2.out.wig "        \
    +config.parameters['path_to_chrom_length_ref']+" "      \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.Unique.posStrand.out.bw "  
    ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)
    
    proc_bigwiggle_3 = subprocess.run((
    config.parameters['softs']['wigToBigWig'] \
    +" -clip "    \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.str1.out.wig "        \
    +config.parameters['path_to_chrom_length_ref']+" "      \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.negStrand.out.bw "  
    ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

    proc_bigwiggle_4 = subprocess.run((
    config.parameters['softs']['wigToBigWig'] \
    +" -clip "    \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.str2.out.wig "        \
    +config.parameters['path_to_chrom_length_ref']+" "      \
    +config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+"_Signal.UniqueMultiple.posStrand.out.bw "  
    ) ,stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True,shell=True)

    write_subprocess_log(proc_bigwiggle_1,logger)
    write_subprocess_log(proc_bigwiggle_2,logger)
    write_subprocess_log(proc_bigwiggle_3,logger)
    write_subprocess_log(proc_bigwiggle_4,logger)

    logger.info("=========> COUNT READS FOR GENE") 

    print_count_reads_per_gene(list_core_sample_names,config)

    logger.info('=========> REMOVE INTERMEDIATE FILES')

    subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Aligned.sortedByCoord.out.bam" ),shell=True)
    subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*.wig" ),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*_ReadsPerGene.out.tab"),shell=True)
    subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Log.out" ),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Log.final.out" ),shell=True)
    subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*_Log.progress.out" ),shell=True)
    subprocess.run(("cat "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*SJ.out.tab > "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/"+config.parameters['final_bam_name']+".All.SplicingJunctions.out.tab"),shell=True)
    #subprocess.run(("rm "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*SJ.out.tab" ),shell=True)
    subprocess.run(("rm  "+config.parameters['path_to_output']+config.parameters['final_bam_name']+"/"+config.chrono+"/*.zip"),shell=True)

    timeline = datetime.datetime.now()
    endtime = str(timeline.hour)+"_"+str(timeline.minute)+"_"+str(timeline.second)

  
    logger.info('=========> END '+endtime)
    logger.info('=========> Files has to be Moved...')
    logger.info('=========> TODO : Check the path for moving the files...')
    
    
    
<<<<<<< HEAD

=======
>>>>>>> 642f291badb64b3ee04f91773db884d09052cc60
