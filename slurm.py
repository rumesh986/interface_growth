from datetime import timedelta

import setup.cmdline as cmd
from setup.slurm import RunSlurm, SlurmParams

if __name__ == '__main__':
	params, system, args = cmd.parse_cmdline_args()

	prog = args.pop('prog')
	source_dir = args.pop('source')

	run = RunSlurm(
		prog=prog,
		params=params,
		source_dir=source_dir,
		system=system,
		workspacesf='workspaces_slurm',
		**args
	)

	if args.pop('build') is True:
		run.build_all()
	
	slurmparams = SlurmParams(
		1,
		1,
		1,
		3882,
		timedelta(hours=24)
	)

	print(slurmparams._time_str)

	# run.prepare_slurmfile(f'{run.wd}/{run.name}.slurm', slurmparams)
	run.prepare_array_job(f'{run.name}.slurm', slurmparams)
	run.test()

	