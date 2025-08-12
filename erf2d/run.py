#!/usr/bin/env python3

import os
import shutil
import subprocess

from datetime import datetime
from concurrent.futures import wait, ProcessPoolExecutor

from visualizer import Visualizer

class Run:
	def __init__(self, prog, xs, ts, nxs, dts, tsteps, dxdt, nthreads, debug=False, interactive=False, **kwargs):
		self.prog = prog
		self.xs = [xs] if isinstance(xs, int) else xs
		self.ts = [ts] if isinstance(ts, int) else ts
		self.nxs = [nxs] if isinstance(nxs, int) else nxs
		self.dts = [dts] if isinstance(dts, float) else dts
		self.tsteps = tsteps
		self.dxdt = dxdt
		self.nthreads = nthreads
		self.debug = debug
		self.interactive= interactive
		self.kwargs = kwargs
		self._exes = {}

		self.od = os.getcwd()

		# self._wd = f'/home/u5665436/repo/oomph-lib/user_drivers/rumesh/erf/test_runs/25-06-2025_15:11:05'

		if self.debug:
			print("Working in debug directory")
			self.wd = f'{self.od}/runs/debug/'
			shutil.rmtree(f'{self.wd}/RESLT/')
			os.mkdir(f'{self.wd}/RESLT')
		else:
			self.wd = f'{self.od}/runs/{datetime.today().strftime("%d-%m-%Y_%H:%M:%S")}'
			os.mkdir(f'{self.wd}')
			os.mkdir(f'{self.wd}/RESLT')
			os.mkdir(f'{self.wd}/imgs')

		print(f'__init__ {prog=} {self.xs=} {ts=} {nxs=} {dts=} {tsteps=} {kwargs=}')
		# print(f'Running in directorys: {self.wd}')

		if not self.interactive:
			self.build()

			os.chdir(self.wd)

			if 't' in self.dxdt:
				self.run_dt()
			else:
				self.run_all()

			print(self.dxdt)
			Visualizer(self.prog, self.dxdt)
			os.chdir(self.od)

	# not parallelized as all processes will be writing to the same file in main directory
	# unsure how to deal with this race condition yet
	def build(self):
		with open("Makefile", "r") as f:
			while line := f.readline():
				if line.startswith("CXXFLAGS"):
					flags_base = line.strip()

		for x in self.xs:
			for t in self.ts:
				if self.debug:
					flags = f'-Wall -O0 -g -DRUN_SCRIPT -DX_ORDER={x} -DT_ORDER={t}'
				else:
					flags = f'{flags_base} -DRUN_SCRIPT -DX_ORDER={x} -DT_ORDER={t}'

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