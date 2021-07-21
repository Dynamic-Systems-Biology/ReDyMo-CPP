#! /bin/python3
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

def CDS_to_SQL():
    with open(sys.argv[1]) as ann_fasta_file:
        output_file = open(sys.argv[2], 'w')
        for section in SimpleFastaParser(ann_fasta_file):
            header = section[0]
            properties = header.split("|")

            # It is a bit ugly, but this part extracts the location from the encoded string
            chrm = properties[3].split("=")[1].split(":")[0]
            start = properties[3].split("=")[1].split(":")[1].split("-")[0]
            end = properties[3].split("=")[1].split( ":")[1].split("-")[1].split("(")[0]
            direction = properties[3].split("=")[1].split("(")[1][:-2]

            # Revert direction
            if direction == '-':
                tmp = start
                start = end
                end = tmp

            sql_statement = f"INSERT INTO TranscriptionRegion (start, end, chromosome_code) VALUES ({start}, {end}, \"{chrm}\");\n"
            output_file.write(sql_statement)
        output_file.close()

def print_help():
    print(f"""
        Usage: {sys.argv[0]} <input_file.fasta> <output_file.sql>

        Options:
         -h, --help \t Show this message
    """)

if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] == "--help" or sys.argv[1] == "-h":
        print_help()
    else:
        CDS_to_SQL()
