import setup.cmdline as cmd
from setup.local import RunLocal

from visualisation.visualisation import Visualiser

if __name__ == '__main__':
	params, system, args = cmd.parse_cmdline_args()

	run = RunLocal(
		prog=args.prog,
		params=params,
		source_dir=args.source,
		system=system,
		dname=args.dname,
		nthreads=args.nthreads
	)

	run.build_all()

	run.run_all()

	print(args)

	vis = Visualiser(
		prefix=args.prog,
		**vars(args)
	)

	vis.default_run(**vars(args))
