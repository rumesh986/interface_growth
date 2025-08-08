#!/usr/bin/env python3

import os
import subprocess

from datetime import datetime
from argparse import ArgumentParser
from concurrent.futures import wait, ProcessPoolExecutor

from visualizer import Visualizer

class Run:
	def __init__(self, prog, xs, ts, nxs, dts, tsteps, dxdt, nthreads, **kwargs):
		self.prog = prog
		self.xs = [xs] if isinstance(xs, int) else xs
		self.ts = [ts] if isinstance(ts, int) else ts
		self.nxs = [nxs] if isinstance(nxs, int) else nxs
		self.dts = [dts] if isinstance(dts, float) else dts
		self.tsteps = tsteps
		self.dxdt = dxdt
		self.nthreads = nthreads
		self.kwargs = kwargs
		self._exes = {}

		self.od = os.getcwd()
		self.wd = f'{self.od}/runs/{datetime.today().strftime("%d-%m-%Y_%H:%M:%S")}'
		# self._wd = f'/home/u5665436/repo/oomph-lib/user_drivers/rumesh/erf/test_runs/25-06-2025_15:11:05'

		os.mkdir(f'{self.wd}')
		os.mkdir(f'{self.wd}/RESLT')
		os.mkdir(f'{self.wd}/imgs')

		print(f'__init__ {prog=} {self.xs=} {ts=} {nxs=} {dts=} {tsteps=} {kwargs=}')
		# print(f'Running in directorys: {self.wd}')

		# self.build()

		# os.chdir(self.wd)

		# self.run_all()

		# _dxdt = ''
		# if self.dxdt == 'bxt':
		# 	if len(self.xs) > 1:
		# 		_dxdt = f'{_dxdt}x'
		# 	if len(self.ts) > 1:
		# 		_dxdt = f'{_dxdt}t'

		# for mode in self.dxdt:
		# 	match mode:
		# 		case x if x in 'abxt': Visualizer(self.prog, x)
		# 		case 'n': pass
		# Visualizer(self.prog, _dxdt)
		# os.chdir(self._od)

	# not parallelized as all processes will be writing to the same file in main directory
	# unsure how to deal with this race condition yet
	def build(self):
		with open("Makefile", "r") as f:
			while line := f.readline():
				if line.startswith("CXXFLAGS"):
					flags_base = line.strip()

		for x in self.xs:
			for t in self.ts:
				flags = f'{flags_base} -DX_ORDER={x} -DT_ORDER={t}'
				subprocess.run(["make", "mostlyclean-compile"])

				with open(f"{self.wd}/build_{x}_{t}.stdout", "w") as f:
					subprocess.run(
						["make", self.prog, flags], 
						check=True, 
						stderr=subprocess.STDOUT,
						stdout=f
					)

				prog_name = f'{self.prog}_{x}_{t}'
				os.rename(self.prog, f'{self.wd}/{prog_name}')
				self._exes[(x, t)] = prog_name

	def _run_base(self, x, t, nx, dt, tsteps, tshift=None, write_freq=None):
		nxargs = ('--nx', str(nx))
		dtargs = ('--dt', str(dt)) if dt is not None else ('--vardt')
		tstepsargs = ('--tsteps', str(tsteps))
		tshiftargs = ('--tshift', str(tshift) if tshift is not None else "0.0")
		write_freqargs = ('--write-freq', str(write_freq) if write_freq is not None else "1")
		
		print(f'running config ({x}, {t}) with {nx=} and {dt=}')

		with open(f'RESLT/{x}n{nx}_{t}t{dt}.stdout', 'w') as f:
			subprocess.run(
				[f'./{self._exes[(x, t)]}', *nxargs, *dtargs, *tstepsargs, *tshiftargs, *write_freqargs], 
				check=True,
				stderr=subprocess.STDOUT,
				stdout=f
			)
	
	def run_all(self):
		futures = []
		with ProcessPoolExecutor(max_workers=self.nthreads) as executor:
			for x in self.xs:
				for t in self.ts:
					for nx in self.nxs:
						for dt in self.dts:
							futures.append(executor.submit(self._run_base, x, t, nx, dt, self.tsteps))
		print("Finished Running for all configs")
	
	def run_dt(self):
		max_dt = max(self.dts)
		futures = []
		with ProcessPoolExecutor(max_workers=self.nthreads) as executor:
			for x in self.xs:
				for t in self.ts:
					tshift = (t+1) * max_dt
					for nx in self.nxs:
						for dt in self.dts:
							tsteps = self.tsteps / dt 
							wf = max_dt / dt
							futures.append(executor.submit(self._run_base, x, t, nx, dt, tsteps, tshift=tshift, write_freq=wf))
		print("Finished Running for all configs")

def trial(xs, ts, **kwargs):
	print('in trial')
	print(xs)
	print(ts)
	print(kwargs)

if __name__ == '__main__':
	parser = ArgumentParser()

	parser.add_argument(
		'prog',
		type=str,
		help="Name of program to run"
	)

	parser.add_argument('--xs', '--x',
		nargs='+', 
		type=int,
		dest='xs',
		default=2,
		help="NNode per elements to simulate with, accepts multiple arguments to be in a list"
	)

	parser.add_argument('--ts', '--t',
		nargs='+', 
		type=int,
		dest='ts',
		default=1,
		help="Order of timestepper to simulate with, accepts multiple in list"
	)

	parser.add_argument('--nxs', '--nx',
		nargs='+', 
		type=int,
		dest='nxs',
		default=100,
		help='Number of elements to simulate with, accepts multiple in list'
	)

	parser.add_argument('--dts', '--dt',
		nargs='+', 
		type=float,
		dest='dts',
		default=1e-4,
		help='timesteps to simulate with, accepts multiple in list'
	)

	parser.add_argument('--tsteps', 
		type=int,
		default=100,
		help='Number of timesteps to simulate for, only single number'
	)

	parser.add_argument('--type',
		type=str,
		choices=['a', 'x', 't', 'b', 'ax', 'at', 'ab'],
		dest='dxdt',
		default='bxt',
		help="Type of analysis to perform a: standard run\nx: vary dx\nt: vary dt"
	)

	parser.add_argument('--nthreads',
		type=int,
		choices=range(1,os.cpu_count()),
		default=os.cpu_count(),
		help='Number of processes that can be started'
	)

	args = parser.parse_args()
	run = Run(**vars(args))