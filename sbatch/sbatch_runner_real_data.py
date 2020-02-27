import subprocess
import string
import os
import sys
import pandas as pd
from optparse import OptionParser, OptionGroup
import shutil

# Change this consts if needed

params_dict = {}

# Don't touch this consts

# Template of sbatch file
SBATCH_TEMPLATE = """#!/bin/bash
#module load anaconda2
#source activate env2
#SBATCH --partition=work
#SBATCH --output="OUTPUT"
#SBATCH --error="ERROR"
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=4
srun "/home/veredhi/EMIRGE_SMURF/run_emirge_smurf_real_data.sh" "CODE_PATH" "FASTQ1" "FASTQ2" "WORKDIR" "REF" "READ_LEN"
"""


def replace_in_string(s, args_dict):
    """
    Replace list of strings in another list of strings, in the same size
    :param s: original string
    :param current_strings: list of strings to replace
    :param replace_strings: replacement strings
    :return: original string, with replacements
    """
    for key, value in args_dict.items():
        s = s.replace(key, value)
    return s


class Experiment():
    def __init__(self, workdir, fastq1, fastq2, ref, read_len):
        self.data = {}
        self.data['FASTQ1'] = fastq1
        self.data['FASTQ2'] = fastq2
        self.data['REF'] = ref
        self.tmp = os.path.join(workdir, "tmp")
        self.data['WORKDIR'] = os.path.join(self.tmp, "workdir")
        self.data['OUTPUT'] = os.path.join(self.tmp, "output")
        self.data['ERROR'] = os.path.join(self.tmp, "error")
        self.data['CODE_PATH'] = "/home/veredhi/EMIRGE_SMURF/Code/"
        self.data["READ_LEN"] = str(read_len)
        self.name = workdir.split("/")[-2]
        print(self.name)

    def run(self):
        os.mkdir(self.tmp, 0755)
        os.mkdir(self.data['WORKDIR'], 0755)

        sbatch_text = replace_in_string(SBATCH_TEMPLATE, self.data)

        sbatch_file_path = os.path.join(self.data["WORKDIR"], self.name)
        print sbatch_file_path
        with open(sbatch_file_path, 'w') as f:
            f.write(sbatch_text)

        command = 'sbatch "{}"'.format(sbatch_file_path)
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        status = process.wait()
        print("Run {} status = {}".format(self.name, status))


def main(argv = sys.argv[1:]):
    """
    Main entry point, create sbatch and run it, for each cloud.
    """

    parser = OptionParser()
    # OPTIONAL
    group_opt = OptionGroup(parser, "Optional parameters",
                            "optional")
    group_opt.add_option("--clean",
                         action="store_true", default=False,
                         help="clean all experiment data")
    (options, args) = parser.parse_args(argv)

    experiments = []
    experiments.append(Experiment("/home/veredhi/results/04-15-R02N1_S34/",
                                  "/home/veredhi/results/04-15-R02N1_S34/04-15-R02N1_S34_L001_R1_001.fastq",
                                  "/home/veredhi/results/04-15-R02N1_S34/04-15-R02N1_S34_L001_R2_001.fastq",
                                  "/home/veredhi/EMIRGE_SMURF/ReferenceFullDB_151",
                                  151))
    experiments.append(Experiment("/home/veredhi/results/06-15-T03N1_S7/",
                                  "/home/veredhi/results/06-15-T03N1_S7/06-15-T03N1_S7_L001_R1_001.fastq",
                                  "/home/veredhi/results/06-15-T03N1_S7/06-15-T03N1_S7_L001_R2_001.fastq",
                                  "/home/veredhi/EMIRGE_SMURF/ReferenceFullDB_151",
                                  151))
    experiments.append(Experiment("/home/veredhi/results/06-15-T15S1_S46/",
                                  "/home/veredhi/results/06-15-T15S1_S46/06-15-T15S1_S46_L001_R1_001.fastq",
                                  "/home/veredhi/results/06-15-T15S1_S46/06-15-T15S1_S46_L001_R2_001.fastq",
                                  "/home/veredhi/EMIRGE_SMURF/ReferenceFullDB_151",
                                  151))
    experiments.append(Experiment("/home/veredhi/results/Batch184_RDB92_GCGTATCA/",
                                  "/home/veredhi/results/Batch184_RDB92_GCGTATCA/RDB92_GCGTATCA_L001_R1_001.fastq",
                                  "/home/veredhi/results/Batch184_RDB92_GCGTATCA/RDB92_GCGTATCA_L001_R2_001.fastq",
                                  "/home/veredhi/EMIRGE_SMURF/ReferenceFullDB_151",
                                  151))
    experiments.append(Experiment("/home/veredhi/results/Batch380_RDB100_TAAGTGGC/",
                                  "/home/veredhi/results/Batch380_RDB100_TAAGTGGC/RDB100_TAAGTGGC_L008_R1_001.fastq",
                                  "/home/veredhi/results/Batch380_RDB100_TAAGTGGC/RDB100_TAAGTGGC_L008_R2_001.fastq",
                                  "/home/veredhi/EMIRGE_SMURF/ReferenceFullDB",
                                  126))

    if options.clean:
        for exp in experiments:
            shutil.rmtree(exp.tmp)
    else:
        for exp in experiments:
            try:
                exp.run()
            except Exception as ex:
                print("Failed running exp = {}: ex = {}".format(exp.name, ex))





if __name__ == '__main__':
    main()
