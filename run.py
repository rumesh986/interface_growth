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

	run.build_all()

	if args.pop('no_run'):
		exit(0)

	run.run_all(**args)

	if args.pop('no_vis'):
		exit(0)

	vis = Visualiser(
		prefix=prog,
		analysis_type=params.analysis_type,
		**args
	)

	vis.default_run(**args)
	# vis.anim_profile(vis.runs[0])
