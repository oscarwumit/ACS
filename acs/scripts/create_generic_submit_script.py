import argparse
from copy import deepcopy
import os

from acs.common import ACS_PATH
from acs.script import generic_submit_script


def parse_command_line_arguments(command_line_args=None):
    """
    Parse command-line arguments.
    Args:
        command_line_args: The command line arguments.
    Returns:
        The parsed command-line arguments by key words.
    """

    parser = argparse.ArgumentParser(description='Automatic Conformer Search (ACS)')
    parser.add_argument('--step', type=int, nargs=1, help='which step in the workflow is being run')
    args = parser.parse_args(command_line_args)
    args.step = args.step[0]

    return args


def main():
	dic  = {1: {'job_name': 'gen_conf',
				'script': os.path.join(ACS_PATH, 'acs', 'screening', 'gen_conf.py'),
				'input': 'input.yml',
				'sub_script_name': 'gen_conf.sh'
				},
			3: {'job_name': 'analyze_screen',
				'script': os.path.join(ACS_PATH, 'acs', 'screening', 'analyze_screen.py'),
				'input': 'initial_conf_screening_result.yml',
				'sub_script_name': 'analyze_screen.sh'
				},
			5: {'job_name': 'analyze_opt_freq_result',
				'script': os.path.join(ACS_PATH, 'acs', 'analysis', 'analyze_opt_freq_result.py'),
				'input': 'opt_project_info.yml',
				'sub_script_name': 'opt_project_info.sh'
				},
			7: {'job_name': 'analyze_fine_opt_freq_result',
				'script': os.path.join(ACS_PATH, 'acs', 'analysis', 'analyze_opt_freq_result.py'),
				'input': 'fine_opt_project_info.yml',
				'sub_script_name': 'fine_opt_project_info.sh'
				},
			9: {'job_name': 'analyze_sp_result',
				'script': os.path.join(ACS_PATH, 'acs', 'analysis', 'analyze_sp_result.py'),
				'input': 'sp_project_info.yml',
				'sub_script_name': 'sp_project_info.sh'
				}
			}
	args = parse_command_line_arguments()

	sub_script = generic_submit_script.format(job_name=dic[args.step]['job_name'],
								   			  stdout=dic[args.step]['job_name'],
								   			  stderr=dic[args.step]['job_name'],
								   			  script=dic[args.step]['script'],
								   			  input=dic[args.step]['input']
								   			   )

	# write the script to the current directory of run_ACS.sh
	with open(dic[args.step]['sub_script_name'], 'w') as f:
	    f.write(sub_script)


if __name__ == '__main__':
    main()