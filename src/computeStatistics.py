##!/usr/bin/python3.5
##!/home/jean-philippe.villemin/bin/python-3.5-2/bin/python3

'''

:date: Oct 25, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Extract Theme

'''
import argparse,textwrap
import subprocess
import os
import re
###########################################################################################################
########################################   Functions   ####################################################
###########################################################################################################




###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Merge statistics from STAR output
    Example : 
    python3 /home/jp/workspace/RNA-SEQ/src/computeStatistics.py -d /home/jp/Desktop/garfield/mount_archive2/commun.luco/EMT_2015_SEQ_data/EMT_RNASEQ_FILES/RNASEQ_2016_RESULTS/T0rep1_2015/10_11_2016__7_56_38/
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-d","--dir",action="store",help="Directory to input.",required=True,type=str,dest='dir')
    parser.add_argument("-n","--name",action="store",help="Name of the files.",required=True,type=str,dest='name')

    parameters = parser.parse_args()
    
  
    d         = {} 
    d["input"]        = 0
    d["uniq"]         = 0
    d["uniq_pc"]      = 0
    d["multiple"]     = 0
    d["multiple_pc"]  = 0
    d["many"]         = 0
    d["many_pc"]      = 0
    total_file = 0
    print(parameters.name)
    for filename in os.listdir(parameters.dir):
      
            m           = re.search('.*_Log\.final\.out',filename)
            if m:
                print(filename)
                total_file +=1
                with open(parameters.dir+filename) as f:
                    for line in f:
                        lineElements = line.strip().split("|")
                        theme = lineElements[0].strip()
                        
                        if (theme == "Number of input reads"):
                            d["input"]        += int(lineElements[1])
                        if (theme == "Uniquely mapped reads number"):
                            d["uniq"]        += int(lineElements[1])
                            
                        if (theme == "Uniquely mapped reads %"):
                            percent           = re.search('(\d+.\d+)%',lineElements[1])
                            if percent :   d["uniq_pc"] += float(percent.group(1))
                            
                        if (theme == "Number of reads mapped to multiple loci"):
                                d["multiple"] += int(lineElements[1])
                        if (theme == "Number of reads mapped to too many loci"):
                                d["multiple"] += int(lineElements[1])
                                    
                        if (theme == "% of reads mapped to multiple loci"):
                                percent           = re.search('(\d+.\d+)%',lineElements[1])
                                if percent :    d["multiple_pc"] += float(percent.group(1))
                        if (theme == "% of reads mapped to too many loci"):
                                percent           = re.search('(\d+.\d+)%',lineElements[1])
                                if percent :    d["multiple_pc"] += float(percent.group(1))                                                           
                  
                f.closed
    f = open("statistics.csv",'a')
    f.write(parameters.name+","+str(d["input"])+","+   str(d["uniq"]) +"("+str("%.2f" % (d["uniq_pc"] /total_file))+"),"+str( d["multiple"])+"("+str("%.2f" %(d["multiple_pc"] /total_file))+")\n")
    f.close
    