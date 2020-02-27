import subprocess
import string
import os
import sys
import pandas as pd
from optparse import OptionParser, OptionGroup

# Change this consts if needed

params_dict = {}

# Don't touch this consts

params_dict['CODE_PATH'] = '/home/veredhi/EMIRGE_SMURF/Code'
params_dict['REFERENCE_PATH'] = '/home/veredhi/EMIRGE_SMURF/ReferenceFullDB'

params_dict['SLURM_OUTPUT_FOLDER'] = 'HOME_PATH/slurm_output'
params_dict['SBATCH_PATH'] = 'HOME_PATH/sbatch_files'
params_dict['PYTHON_OUTPUT_FOLDER'] = 'HOME_PATH/python_results'
params_dict['TMP_FOLDER'] = 'HOME_PATH/tmp'
params_dict['TAXA_PATH'] = '/home/veredhi/EMIRGE_SMURF/taxa_and_header.csv'

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
srun "SCRIPT_PATH" c_bases c_references mix_size unique "CODE_PATH" "PYTHON_OUTPUT_FOLDER" "REFERENCE_PATH" "TMP_FOLDER" index "TAXA_PATH"
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
    for key, value in args_dict.items():
        s = s.replace(key, value)
    for key, value in args_dict.items():
        s = s.replace(key, value)
    return s


def mkdir_if_not_exists(dir_path):
    if not os.path.isdir(dir_path):
        os.mkdir(replace_in_string(dir_path, params_dict), 0755)


def main(argv = sys.argv[1:]):
    """
    Main entry point, create sbatch and run it, for each cloud.
    """

    parser = OptionParser("""usage: %prog  RESULT_DIR RUNNIGN_SCRIPT <optional parameters>""")

    # OPTIONAL
    group_opt = OptionGroup(parser, "Optional parameters",
                            "Defaults should normally be fine for these options in order to run batch")
    group_opt.add_option("-r", "--repeats",
                         type="int", default=50,
                         help="""number of repeats of each configuration""")
    group_opt.add_option("-a", "--args",
                         type="string", default="//home//veredhi//EMIRGE_SMURF//cluster_summary.csv",
                         help="arguments csv path")

    parser.add_option_group(group_opt)

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    # minimal sanity checking of input
    if len(args) != 2:
        parser.error(
            "RESULT_DIR ans RUNNIGN_SCRIPT are required, and all options except should have a flag associated with them (options without flags: %s)" % args)

    params_dict['HOME_PATH'] = os.path.join(os.path.abspath(args[0]), "CONFIG_INDEX")
    params_dict['REPEATS'] = str(options.repeats)

    args_path = os.path.abspath(options.args)
    params_dict['SCRIPT_PATH'] = os.path.abspath(args[1])

    args_df = pd.read_csv(args_path)

    for i in range(args_df.shape[0]):
        params_dict['CONFIG_INDEX'] = str(i)

        mkdir_if_not_exists(replace_in_string(params_dict['HOME_PATH'], params_dict))
        mkdir_if_not_exists(replace_in_string(params_dict['SLURM_OUTPUT_FOLDER'], params_dict))
        mkdir_if_not_exists(replace_in_string(params_dict['SBATCH_PATH'], params_dict))
        mkdir_if_not_exists(replace_in_string(params_dict['PYTHON_OUTPUT_FOLDER'], params_dict))
        mkdir_if_not_exists(replace_in_string(params_dict['TMP_FOLDER'], params_dict))

        row = args_df.iloc[i]
        params_dict["c_bases"] = str(row['number of changed bases'])
        params_dict["c_references"] = str(row['number of changed bacterias'])
        params_dict["mix_size"] = str(row['mock_mixure_size'])
        params_dict["unique"] = str(row['unique_bacteria_in_mixture'])

        for j in range(int(params_dict['REPEATS'])):

            index=str(j)
            test_name = "test_" + index

            params_dict["index"] = index

            if os.path.isfile(os.path.join(replace_in_string(params_dict["PYTHON_OUTPUT_FOLDER"], params_dict),
                                           test_name, "emirge_smurf_WFalseSTrue.csv")):
                continue

            params_dict["out_file"] = index + '_out'
            params_dict["err_file"] = index + '_err'

            sbatch_text = replace_in_string(SBATCH_TEMPLATE, params_dict)

            sbatch_file_path = os.path.join(replace_in_string(params_dict['SBATCH_PATH'], params_dict), test_name)
            with open(sbatch_file_path, 'w') as f:
                f.write(sbatch_text)

            command = 'sbatch "{}"'.format(sbatch_file_path)
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            status = process.wait()
            while status is not 0:
                status = process.wait()
            print("Add {} {}".format(i, test_name))



if __name__ == '__main__':
    main()
