import pandas as pd
from Bio import SeqIO, pairwise2
import os
import sys
from optparse import OptionParser, OptionGroup


if __name__ == "__main__":
    parser = OptionParser("Cut the fastq sequences into single size.")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "")
    group_reqd.add_option("-f", "--fastq_path",
                      type="string", default="/home/vered/EMIRGE/data/s26_mock/dada2/SampRS26_rev_fasta_mod.fastq",
                      help="path to output new database location")
    group_reqd.add_option("-l", "--length",
                          type="int", default=145,
                          help="length of sequence")

    parser.add_option_group(group_reqd)
    (options, args) = parser.parse_args(sys.argv[1:])

    new_fastq_path = "{}_{}{}".format(options.fastq_path.split('.fastq')[0], options.length, ".fastq")

    with open(options.fastq_path, 'r') as f_i:
        with open(new_fastq_path, 'w') as f_o:
            for line in f_i:
                line = "".join(line.strip()[:options.length+1])
                f_o.write(line + "\n")