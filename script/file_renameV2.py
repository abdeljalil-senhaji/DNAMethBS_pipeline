# file_rename.py
# Created: 24/06/2021

"""
This will rename a group of files in a given directory, 
once you have provided the list of forward and reverse reads
"""

# just checking
__author__ = 'Abdeljalil Senhaji Rachik'
__version__ = '1.0'

import argparse
import os

def get_parser():
    parser = argparse.ArgumentParser(description='change extension of files in a working directory')
    parser.add_argument('work_dir', metavar='WORK_DIR', type=str, nargs=1, help='the directory where to change file')
    parser.add_argument('list_forword', metavar='LIST_FORWORD', type=str, nargs=1, help='forword list R1')
    parser.add_argument('list_reverse', metavar='LIST_REVERSE', type=str, nargs=1, help='reverse list R2')
    return parser

def main():
    """
    This will be called if the script is directly invoked.
    """
    # adding command line argument
    parser = get_parser()
    args = vars(parser.parse_args())
    # adding command line argument
    originalefiles=os.open(args['work_dir'])
    for filename in os.open(args['list_forword']):
        if filename in originalefiles :
            stripfile=filename.replace('.fastq.gz','').replace('.fq.gz', '').replace('_R1', '').replace('_R', '').replace('.R1', '').replace('.R', '')
            os.rename(filename ,f"{stripfile}._R1.fastq.gz")
       
    for filename in os.open(args['list_reverce']):     
        if filename in originalefiles:
            stripfile=filename.replace('.fastq.gz','').replace('.fq.gz', '').replace('_R1', '').replace('_R', '').replace('.R1', '').replace('.R', '')
            os.rename(filename ,f"{stripfile}._R2.fastq.gz")
             
if __name__ == '__main__':
    main()



