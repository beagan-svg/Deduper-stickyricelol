'''
Beagan Nguy
Bi_624
11/03/2021

This program will remove all PCR duplicates and retains only a single copy of each read from a sam file.

Usage: python3 final_deduper.py -s <sorted_sam_file> -u <umi.txt>

Input: 1) a sorted sam file, 2) a text file containing a list of umi

Output: 3 Files
1) duplicates.sam, File containing all duplicate reads from the input sam file
2) unique.sam, File containing all non-duplicated reads from the input sam file
3) stat.txt, File containing the chromosome and number of non-duplicated reads per chromosome

'''

import sys
import argparse
import re

# Arg_parser
def getArgs():
    parser = argparse.ArgumentParser(
        description='Optional arguments: -s, -p, -u', usage='python3 final_deduper.py -s <sorted_sam_file> -u <umi.txt>'
    )

    parser.add_argument(
        '-s', '-samfile', help='Input sorted samfile', type=str, default='/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam', required=True
    )
 
    parser.add_argument(
        '-p', '-paired_end', help='Specify paired_end option', type=str, required=False
    )

    parser.add_argument(
        '-u', '-umi', help='Specify umi option', type=str, required=False
    )

    return parser.parse_args()

# Read umi list and return a list of umi
def readUmi(umi_read):
    umi_list = list()
    for line in umi_read:
        umi_list.append(line.rstrip())
    return umi_list

# Check strandness
def checkStrand(flag):
    strand = ''
    if ((int(flag) & 16) == 16):
        strand = 'reverse'
    else:
        strand = 'forward'
    return strand

# Function to softclip on the forward strand 
def softClipPlus(pos, cigar):
    adjusted_pos = int(pos)
    try:
        clip = re.match('^([0-9]+)S', cigar) 
        adjusted_pos = int(pos) - int(clip[1])
    except:
        pass
    return str(adjusted_pos)

# Function to softclip on the reverse strand
def softClipMinus(pos, cigar):
    soft_clip = 0
    cigar_list = list()
    try:
        char_s = re.search('([0-9]+)S$', cigar)
        soft_clip += int(char_s[1])
    except:
        pass

    char_list = re.findall('[0-9]+[MDN]', cigar)
    for char in char_list:
        cigar_list.append(int(char[:-1]))
    
    adjusted_pos = sum(cigar_list) + soft_clip + int(pos)
    return str(adjusted_pos)

# Creates a unique identifier for each read (umi, chrom, strand, adjusted_pos)
# Duplicated and unique reads will be outputed to duplicates.txt and unique.txt respectively.
def deduplicate(sam_open, umi_list, dup_out, uniq_out, stat_out):
    dedup_dict = dict()     # Key: Unique Identifier, Value: ''
    chrom_count_dict = dict()
    prev_chrom = 1
    invalid_umi_count = 0   # Hold counts for each umi in the sam file not in the umi list

    # Read sam file
    for line in sam_open:
        
        # Write headers to output files
        if line[0] == "@":
            dup_out.write(line) 
            uniq_out.write(line)
        else:
            sam_col = line.split('\t')
            qname = sam_col[0].split(':')

            # Define column fields
            flag = sam_col[1]
            chrom = sam_col[2]
            pos = sam_col[3]
            cigar = sam_col[5]
            umi = qname[7]
            
            # Next chromsome condition
            # Clears running dedup dictionary for memory limits and write out stats file
            if chrom != prev_chrom:
                prev_chrom = chrom
                dedup_dict = dict()

            # Create a dictionary counts for each chromosome
            if chrom not in chrom_count_dict:
                chrom_count_dict[chrom] = 0

            # Check to see if sam umi is in umi list
            if umi_list != None:
                if umi in umi_list:
                    pass
                else:
                    invalid_umi_count += 1
                    continue
            else:
                pass
            
            # Check strand of the read based on bitwise flag
            strand = checkStrand(flag)
            adjusted_pos = pos 

            # Adjust start position base on strand
            if strand == 'forward':
                adjusted_pos = softClipPlus(pos, cigar)
                
            if strand == 'reverse':
                adjusted_pos = softClipMinus(pos, cigar)

            # Create a unique identifier for each read
            identifier = (chrom, strand, adjusted_pos, umi)
            
            # Check for duplicates, if a read has already been recorded output to duplicate file, otherwise output to unique output file
            if identifier in dedup_dict:
                dup_out.write(line)
            else:
                uniq_out.write(line)
                chrom_count_dict[chrom] += 1
                dedup_dict[identifier] = ''

    # Create a stat output file containing the chromosome and number of non-duplicated reads per chrom
    for chrom, count in chrom_count_dict.items():
        stat_out.write('{} : {}\n'.format(chrom, count))
    
    # Write number of invalid umis to stat.txt file
    stat_out.write('Invalid Umis: {}'.format(invalid_umi_count))

def main():
    # Parse through given files
    args = getArgs()
    sam_read = args.s
    
    # Handles cases where no umi list are inputed
    umi_list = list()
    if args.u:
        try:
            umi_read = args.u
            umi_read = open(umi_read, 'r')  # Open inputed umi file
            umi_list = readUmi(umi_read)    # Read umi file
        except:
            print('Incorrect umi list file format, will continue assuming no umi list is inputed')
    else:
        umi_list = None

    # Adds error when paired-end option given.
    if args.p != None:
        print("No paired-end functionality yet")
        sys.exit()
    
    # Open Sam and Umi file
    sam_open = open(sam_read, 'r') 

    # Name output file
    out_dup_str = 'duplicates.sam'
    out_unique_str = 'unique.sam'
    out_stat_str = 'stat.txt'

    # Open output file
    dup_out = open('{}'.format(out_dup_str),'w')
    uniq_out = open('{}'.format(out_unique_str),'w')
    stat_out = open('{}'.format(out_stat_str), 'w')

    # Main deduplicate method
    deduplicate(sam_open, umi_list, dup_out, uniq_out, stat_out)

    # Close files
    dup_out.close()
    uniq_out.close()
    
if __name__ == "__main__":
    main()


