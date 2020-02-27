import subprocess
import string
import os
import pandas as pd

# Change this consts if needed
ARGS_PATH = "/home/veredhi/EMIRGE_SMURF/cluster_summary.csv"
REPEATS = 50

# Don't touch this consts
HOME_PATH = '/RG/compbio/groupData/results/vered/cluster_mutated_mix2/CONFIG_INDEX'
CODE_PATH = '/home/veredhi/EMIRGE_SMURF/Code'
REFERENCE_PATH = '/home/veredhi/EMIRGE_SMURF/ReferenceFullDB'

SLURM_OUTPUT_FOLDER = 'HOME_PATH/slurm_output'.replace('HOME_PATH', HOME_PATH)
SBATCH_PATH = 'HOME_PATH/sbatch_files'.replace('HOME_PATH', HOME_PATH)
PYTHON_OUTPUT_FOLDER = 'HOME_PATH/python_results'.replace('HOME_PATH', HOME_PATH)
TMP_FOLDER = 'HOME_PATH/tmp'.replace('HOME_PATH', HOME_PATH)
TAXA_PATH = '/home/veredhi/EMIRGE_SMURF/taxa_and_header.csv'

# Template of sbatch file
SBATCH_TEMPLATE = """#!/bin/bash
#module load anaconda2
#source activate env2
#SBATCH --partition=work
#SBATCH --output="SLURM_OUTPUT_FOLDER/out_file"
#SBATCH --error="SLURM_OUTPUT_FOLDER/err_file"
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=12GB
#srun /home/veredhi/EMIRGE_SMURF/run_emirge_smurf.sh c_bases c_references mix_size unique "CODE_PATH" "PYTHON_OUTPUT_FOLDER" "REFERENCE_PATH" "TMP_FOLDER" index
#srun /home/veredhi/EMIRGE_SMURF/run_emirge_smurf_mutated.sh c_bases c_references mix_size unique "CODE_PATH" "PYTHON_OUTPUT_FOLDER" "REFERENCE_PATH" "TMP_FOLDER" index
srun /home/veredhi/EMIRGE_SMURF/run_emirge_smurf_mutated_mix2.sh c_bases c_references mix_size 97 "CODE_PATH" "PYTHON_OUTPUT_FOLDER" "REFERENCE_PATH" "TMP_FOLDER" index
#srun /home/veredhi/EMIRGE_SMURF/run_emirge_smurf_family.sh c_bases c_references mix_size unique "CODE_PATH" "PYTHON_OUTPUT_FOLDER" "REFERENCE_PATH" "TMP_FOLDER" index "TAXA_PATH"
"""

def replace_in_string(s, current_strings, replace_strings):
    """
    Replace list of strings in another list of strings, in the same size
    :param s: original string
    :param current_strings: list of strings to replace
    :param replace_strings: replacement strings
    :return: original string, with replacements
    """
    for i in range(len(current_strings)):
        s = s.replace(current_strings[i], replace_strings[i])
    return s


def main():
    """
    Main entry point, create sbatch and run it, for each cloud.
    """
    args_df = pd.read_csv(ARGS_PATH)
    #for i in range(6, args_df.shape[0]):
    for i in range(6,10):
        CONFIG_INDEX=str(i)

        os.mkdir(HOME_PATH.replace("CONFIG_INDEX", CONFIG_INDEX), 0755)
        os.mkdir(SLURM_OUTPUT_FOLDER.replace("CONFIG_INDEX", CONFIG_INDEX), 0755)
        os.mkdir(SBATCH_PATH.replace("CONFIG_INDEX", CONFIG_INDEX), 0755)
        os.mkdir(PYTHON_OUTPUT_FOLDER.replace("CONFIG_INDEX", CONFIG_INDEX), 0755)
        os.mkdir(TMP_FOLDER.replace("CONFIG_INDEX", CONFIG_INDEX), 0755)

        row = args_df.iloc[i]
        c_bases = str(row['number of changed bases'])
        c_references = str(row['number of changed bacterias'])
        mix_size = str(row['mock_mixure_size'])
        unique = str(row['unique_bacteria_in_mixture'])

        for j in range(REPEATS):
            index=str(j)
            out_file = index + '_out'
            err_file = index + '_err'

            test_name="test_" + index

            sbatch_text = replace_in_string(SBATCH_TEMPLATE,
                                            ["HOME_PATH", "CODE_PATH","REFERENCE_PATH","SLURM_OUTPUT_FOLDER","SBATCH_PATH", "PYTHON_OUTPUT_FOLDER","TMP_FOLDER", "c_bases", "c_references", "mix_size", "unique", "index", "out_file", "err_file", "CONFIG_INDEX", "TAXA_PATH"],
                                            [HOME_PATH, CODE_PATH,REFERENCE_PATH,SLURM_OUTPUT_FOLDER,SBATCH_PATH, PYTHON_OUTPUT_FOLDER,TMP_FOLDER, c_bases, c_references, mix_size, unique, index, out_file, err_file, CONFIG_INDEX, TAXA_PATH])

            sbatch_file_path = os.path.join(SBATCH_PATH.replace("CONFIG_INDEX", CONFIG_INDEX), test_name)
            with open(sbatch_file_path, 'w') as f:
                f.write(sbatch_text)

            command = 'sbatch "{}"'.format(sbatch_file_path)
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            process.wait()


if __name__ == '__main__':
    main()
