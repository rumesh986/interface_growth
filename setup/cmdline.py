import os
from itertools import product
from argparse import ArgumentParser, ArgumentError, BooleanOptionalAction

from systems.params import SimParams, AnalysisType
from systems.materials import MaterialSystem
from visualisation.visualisation import PltStyle

def parse_cmdline_args():
	parser = ArgumentParser()

	params_group = parser.add_argument_group('Simulation parameters')
	materials_group = parser.add_argument_group('Material properties')
	post_group = parser.add_argument_group('Post processing arguments')


	parser.add_argument(
		'prog',
		type=str,
		help="Name of program to run"
	)

	parser.add_argument(
		'source',
		type=str,
		help='Name of project'
	)

	params_group.add_argument('--name',
		type=str,
		dest='dname',
		default=None,
		help='Name of directory to store data'
	)

	params_group.add_argument('--xs', '--x',
		nargs='+', 
		type=int,
		dest='xs',
		default=[2],
		help="NNode per elements in one dimension to simulate with, accepts multiple arguments to be in a list"
	)

	params_group.add_argument('--ts', '--t',
		nargs='+', 
		type=int,
		dest='ts',
		default=[1],
		help="Order of timestepper to simulate with, accepts multiple in list"
	)

	nx_group = params_group.add_mutually_exclusive_group(required=True)

	nx_group.add_argument('--nx', '--nxs',
		nargs='*',
		dest='nxs',
		help='Number of elements to simulate with in x-direction. Expects a list of tuples. eg: 10,100 will have 10 elements in solid phase and 100 in liquid phase'
	)

	nx_group.add_argument('--nx1',
		nargs='*', 
		type=int,
		dest='nx1',
		help='Number of elements to simulate with in x-direction for solid phase, accepts multiple in list'
	)

	params_group.add_argument('--nx2',
		nargs='*',
		type=int,
		dest='nx2',
		help='Number of elements to simulate with in x-direction for liquid phase, accepts multiple in list'
	)

	params_group.add_argument('--ny',
		nargs='+',
		type=int,
		default=[1],
		help='Number of elements in y-direction of domain'
	)

	params_group.add_argument('--k',
		type=int,
		default=1,
		help='wavenumber of perturbation - 2pi * k'
	)

	params_group.add_argument('--A',
		type=float,
		default=0.001,
		help='Amplitude of perturbation'
	)

	params_group.add_argument('--dts', '--dt',
		nargs='+', 
		type=float,
		dest='dts',
		default=[1e-4],
		help='timesteps to simulate with, accepts multiple in list'
	)

	params_group.add_argument('--tstart', '--start',
		type=float,
		dest='tstart',
		default=0.5,
		help='Start time of simulation'
	)

	params_group.add_argument('--tend', '--tend',
		type=float,
		dest='tend',
		default=1.0,
		help='End time of simulation'
	)

	params_group.add_argument('--write-freq', '--wf', 
		type=int,
		dest='wf',
		default=1,
		help='Number of timesteps to wait for before saving results'
	)

	params_group.add_argument('--nthreads',
		type=int,
		choices=range(1,os.cpu_count()),
		default=os.cpu_count(),
		help='Number of processes that can be started'
	)

	materials_group.add_argument('--material',
		type=str,
		choices=MaterialSystem.list_default_materials(),
		default='water_ice_und',
		help='Material system to use'
	)

	materials_group.add_argument('--Ts',
		type=float,
		default=None,
		help='Boundary condition at x = 0.0 (solid)'
	)

	materials_group.add_argument('--Tm',
		type=float,
		default=None,
		help='Boundary condition at interface'
	)

	materials_group.add_argument('--Tl',
		type=float,
		default=None,
		help='Boundary condition at x = 1.0 (liquid)'
	)

	post_group.add_argument('--type',
		type=AnalysisType,
		dest='analysis_type',
		default=AnalysisType.standard,
		choices=list(AnalysisType),
		help="Type of analysis to perform a: standard run\ndx: vary dx\ndt: vary dt\ns: sensitivity analysis"
	)

	post_group.add_argument('--style',
		type=PltStyle,
		dest='style',
		default=PltStyle.standard,
		choices=list(PltStyle),
		help='Stylesheet for matplotlib'
	)

	post_group.add_argument('--anim-length',
		type=int,
		help='Length of animation when using standard processing, in seconds'
	)

	parser.add_argument('--clear-results',
		action='store_true',
		help='Delete existing results folders if they exist'
	)

	parser.add_argument('--build',
		action=BooleanOptionalAction,
		default=True,
		help='Build program'
	)

	parser.add_argument('--run',
		action=BooleanOptionalAction,
		default=True,
		help='run program'
	)

	parser.add_argument('--vis',
		action=BooleanOptionalAction,
		default=True,
		help='Perform post-processing'
	)

	parser.add_argument('--force-vis',
		action='store_true',
		help='Perform post-processing even if main run fails to converge'
	)

	args = parser.parse_args()

	if args.nxs is not None and args.nx2 is not None:
		parser.error('argument --nx2 not allowed with argument --nx/--nxs')

	args = vars(args)

	nxs = []
	if args.get('nxs') is not None:
		for pair in args.pop('nxs'):
			nx1, nx2 = pair.split(',')
			nxs.append((int(nx1), int(nx2)))
		args.pop('nx1')
		args.pop('nx2')
	else:
		if args.get('nx2') is not None:
			for nx1, nx2 in product(args.pop('nx1'), args.pop('nx2')):
				nxs.append((nx1, nx2))
		else:
			for nx1 in args.pop('nx1'):
				nxs.append((nx1, nx1))
		args.pop('nx')

	nys = []
	for ny in args.pop('ny'):
		nys.append(ny)

	params = SimParams(
		args.pop('xs'),
		args.pop('ts'),
		nxs,
		nys,
		args.pop('dts'),
		args.pop('tstart'),
		args.pop('tend'),
		args.pop('wf'),
		args.pop('analysis_type')
	)

	system = MaterialSystem.default_material(args.pop('material'))

	# update temperatures if needed
	for T in ['Ts', 'Tm', 'Tl']:
		temp = args.pop(T)
		if temp is not None:
			system.__setattr__(T, temp)

	return params, system, args