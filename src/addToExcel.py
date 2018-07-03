
'''

:date: Sep 6, 2016
:platform: Ubuntu 16.04

:author: Villemin Jean-Philippe
:team: Epigenetic Component of Alternative Splicing - IGH

:synopsis: Convert to excel.

'''
import argparse,textwrap
import os
import pandas as pd
from openpyxl import load_workbook

###########################################################################################################
########################################   Main   #########################################################
###########################################################################################################

if __name__ == '__main__':
    '''This is just a main '''
    
    parser = argparse.ArgumentParser(description=textwrap.dedent ('''\
    Add regular file as new tab in an Excel File.
    
    '''),formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-f","--file",action="store",help="Tsv/csv file added to excel.",required=True,type=str,dest='file')
    parser.add_argument("-e","--excel",action="store",help="Excel file to upgrade.",required=True,type=str,dest='excel')
    parameters = parser.parse_args()

    print(parameters.file)
    print(parameters.excel)
 
    book = load_workbook(parameters.excel)
    writer = pd.ExcelWriter(parameters.excel, engine='openpyxl')
    writer.book = book
    #Maximum 31 characters allowed in sheet title
    print(os.path.basename(parameters.file))#[:-18])
    sheetTile=""
    if( len(os.path.basename(parameters.file)) >= 31 ) :
        sheetTile = os.path.basename(parameters.file)[:30]
        print(len(os.path.basename(parameters.file)))
        print(sheetTile)

    else : sheetTile = os.path.basename(parameters.file)
    
    pd.read_csv(parameters.file,sep='\t',header=None).to_excel(writer,sheetTile,startrow=0,startcol=0,header=False,index=False)
    writer.save()
    
    
        
