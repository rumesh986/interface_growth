import setup.cmdline as cmd
from setup.local import RunLocal

from visualisation.visualisation import Visualiser

if __name__ == '__main__':
	params, system, args = cmd.parse_cmdline_args()

	prog = args.pop('prog')
	source_dir = args.pop('source')

	run = RunLocal(
		prog=prog,
		params=params,
		source_dir=source_dir,
		system=system,
		**args
	)

	if args.pop('build') is True:
		run.build_all()

	if args.pop('run') is True:
		try:
			run.run_all()
		except Exception as e:
			if not args.pop('force_vis'):
				raise e

	if args.pop('vis') is True:
		vis = Visualiser(
			prefix=prog,
			analysis_type=params.analysis_type,
			**args
		)

		vis.default_run(**args)