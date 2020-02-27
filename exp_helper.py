import os
from optparse import OptionParser, OptionGroup
import sys

ACTUAL_RES_FILE_NAME = "emirge_smurf_WFalseSTrue.csv"

def move_test(test_name, source_path, target_path, suffix=''):
    dirs = ['python_results', 'sbatch_files']
    for dir in dirs:
        dest_dir = os.path.join(target_path, dir)
        if not os.path.isdir(dest_dir):
            os.mkdir(dest_dir)
        copy_from = os.path.join(source_path, dir, test_name)
        copy_to =  os.path.join(dest_dir, test_name+suffix)
        # print("copy {} _to {}".format(copy_from, copy_to))
        os.rename(copy_from,copy_to)

    # slurm_output_dest = os.path.join(target_path, 'slurm_output')
    # slurm_output_source = os.path.join(source_path, 'slurm_output')
    # if not os.path.isdir(slurm_output_dest):
    #     os.mkdir(slurm_output_dest)
    # index = test_name.split('_')[1]
    # os.rename(os.path.join(slurm_output_source, index + '_err' + suffix), os.path.join(slurm_output_dest, index + '_err' + suffix))
    # os.rename(os.path.join(slurm_output_source, index + '_out' + suffix), os.path.join(slurm_output_dest, index + '_out'+ suffix))

    tmp_dest = os.path.join(target_path, 'tmp')
    tmp_source = os.path.join(source_path, 'tmp')
    if not os.path.isdir(tmp_dest):
        os.mkdir(tmp_dest)
    index = test_name.split('_')[1]
    try:
        os.rename(os.path.join(tmp_source, 'mock_ix_' + index), os.path.join(tmp_dest, 'mock_ix_' + index + suffix))
    except Exception as ex:
        print (ex)


def combine_experiment(from_path, to_path):
    errors_path = os.path.join(to_path, "Errors")
    results_path = os.path.join(to_path, "python_results")

    if not os.path.isdir(to_path):
        os.mkdir(to_path)
        os.mkdir(results_path)

    if not os.path.isdir(errors_path):
        os.mkdir(errors_path)

    test_dirs = os.listdir(results_path)
    # move errors to error directory
    for test_dir in test_dirs:
        test_path = os.path.join(results_path, test_dir)
        actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME)
        if not os.path.isfile(actual_path):
            move_test(test_dir, to_path, errors_path, '_to')

    # move errors to error directory and success to to_path
    from_results_path = os.path.join(from_path, "python_results")
    from_test_dirs = os.listdir(from_results_path)
    if os.listdir(from_path):
        for test_dir in from_test_dirs:
            test_path = os.path.join(from_results_path, test_dir)
            actual_path = os.path.join(test_path, ACTUAL_RES_FILE_NAME)
            if not os.path.isfile(actual_path):
                # print("path not exist = {}".format(actual_path))
                move_test(test_dir, from_path, errors_path, '_a')
            else:
                move_test(test_dir, from_path, to_path, '_a')

    print("Counter = {} in path = {}".format(len(os.listdir(results_path)), to_path))
    errors_path = os.path.join(errors_path, "python_results")
    if os.path.isdir(errors_path):
        print("Errors counter = {} in path = {}, errors = {}".format(len(os.listdir(errors_path)), errors_path, os.listdir(errors_path)))
        print("----")

def get_cmd_arguments(argv = sys.argv[1:]):
    USAGE = \
        """usage: %prog FROM_PAH TO_PATH [options]

        Experiments helper
        """

    parser = OptionParser(USAGE)

    # group_opt = OptionGroup(parser, "Optional flags",
    #                          "These flags are all optional.")
    #
    # group_opt.add_option("-c", "--combine",
    #                      action="store_true", default=False,
    #                      help="combine two expriments into one")
    # parser.add_option_group(group_opt)

    (options, args) = parser.parse_args(argv)

    if len(args) != 2:
        parser.error(
            "FROM_PATH and TO_PATH is required, and all options except DIR should have a flag associated with them (options without flags: %s)" % args)
        return

    return options, args

def combine_experiments(from_path, to_path):
    experiments = os.listdir(from_path)
    for experiment in experiments:
        if os.path.isdir(os.path.join(from_path, experiment)):
            from_experiment_path = os.path.join(from_path, experiment)
            combine_experiment(from_experiment_path, os.path.join(to_path, experiment))
            # print("remove path = {}".format(from_experiment_path))
            #os.remove(from_experiment_path)


def main():

    # ACTUALLY PARSE ARGS
    (options, args) = get_cmd_arguments()

    from_path = os.path.abspath(args[0])
    to_path = os.path.abspath(args[1])
    combine_experiments(from_path, to_path)


if __name__ == "__main__":
    main()