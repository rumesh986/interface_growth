#!/usr/bin/env python3

import os
import time
import shutil
import subprocess

from datetime import datetime
from concurrent.futures import wait, ProcessPoolExecutor

import numpy as np

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

			self.run(**kwargs)

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
					flags = f'CXXFLAGS=-Wall -O0 -g -DRUN_SCRIPT -DX_ORDER={x} -DT_ORDER={t}'
				else:
					flags = f'{flags_base} -DRUN_SCRIPT -DX_ORDER={x} -DT_ORDER={t}'

				subprocess.run(["make", "mostlyclean-compile"])

				try:
					with open(f"{self.wd}/build_{x}_{t}.stdout", "w") as f:
						subprocess.run(
							["make", self.prog, flags], 
							check=True, 
							stderr=subprocess.STDOUT,
							stdout=f
						)
				except subprocess.CalledProcessError:
					print(f"Build failed with X={x} T={t}")
					with open(f"{self.wd}/build_{x}_{t}.stdout", "r") as f:
						print(f.read())
					raise

				prog_name = f'{self.prog}_{x}_{t}'
				os.rename(self.prog, f'{self.wd}/{prog_name}')
				self._exes[(x, t)] = prog_name

	def _run_base(self, x, t, nx, dt, tsteps, tshift=None, write_freq=None, **kwargs):
		nxargs = ('--nx', str(nx))
		dtargs = ('--dt', str(dt)) if dt is not None else ('--vardt')
		tstepsargs = ('--tsteps', str(tsteps))
		tshiftargs = ('--tshift', str(tshift) if tshift is not None else "0.0")
		write_freqargs = ('--write-freq', str(write_freq) if write_freq is not None else "1")

		translated_kwargs = []
		for key, value in kwargs.items():
			if key == 'tshift' or value is None or key == 'eps':
				continue
			translated_kwargs.extend([f'--{key}', str(value)])# if value is not str else value])
		
		print(f'Running config {x=}, {t=}, {nx=} and {dt=}')

		try:
			if 'dname' not in kwargs.keys():
				kwargs['dname'] = f'RESLT/{x}n{nx}_{t}t{dt}'
			with open(f'{kwargs["dname"]}.stdout', 'w') as f:
				print(f'{kwargs["dname"]}: {translated_kwargs}')
				subprocess.run(
					[f'./{self._exes[(x, t)]}', *nxargs, *dtargs, *tstepsargs, *tshiftargs, *write_freqargs, *translated_kwargs], 
					check=True,
					stderr=subprocess.STDOUT,
					stdout=f
				)
		except subprocess.CalledProcessError as e:
			print(f"Run failed with {x=} {t=}")
			with open(f'{kwargs["dname"]}.stdout', 'r') as f:
				print(f.read())
			raise
	
	def run(self, tshift=None, **kwargs):
		num_workers = len(self.xs) * len(self.ts) * len(self.nxs) * len(self.dts)
		if (num_workers > self.nthreads):
			num_workers = self.nthreads
		
		with open('config', 'w') as f:
			f.write(f'xs: {self.xs}\n')
			f.write(f'ts: {self.ts}\n')
			f.write(f'nxs: {self.nxs}\n')
			f.write(f'dts: {self.dts}\n')
			f.write(f'tsteps: {self.tsteps}\n')
			f.write(f'dxdt: {self.dxdt}\n')
			f.write(f'kwargs: {self.kwargs}')

		futures = {}
		wf = None
		with ProcessPoolExecutor(max_workers=num_workers) as executor:
			if 't' in self.dxdt:
				max_dt = max(self.dts)
				tshift = (max(self.ts) + 1) * max_dt
				for x in self.xs:
					for t in self.ts:
						for nx in self.nxs:
							for dt in self.dts:
								tsteps = self.tsteps / dt
								wf = max_dt / dt
								futures[(x, t, nx, dt)] = executor.submit(self._run_base, x, t, nx, dt, tsteps, tshift=tshift, write_freq=wf, **kwargs)
			elif 's' in self.dxdt:
				for x in self.xs:
					for t in self.ts:
						for nx in self.nxs:
							for dt in self.dts:
								kwargs['dname'] = f'RESLT/{x}n{nx}_{t}t{dt:.2e}_{kwargs["k1"]:9.7f}_{kwargs["k2"]:9.7f}_{kwargs["rho1"]:9.7f}_{kwargs["rho2"]:9.7f}_{kwargs["cp1"]:9.7f}_{kwargs["cp2"]:9.7f}'
								futures[(x, t, nx, dt, np.random.rand(1)[0])] = executor.submit(self._run_base, x, t, nx, dt, self.tsteps, tshift=tshift, write_freq=wf, **kwargs)
								for v in ['k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2']:
									kwargs_cp = kwargs.copy()
									kwargs_cp[v] *= (1.0 + kwargs_cp['eps'])
									kwargs_cp['dname'] = f'RESLT/{x}n{nx}_{t}t{dt:.2e}_{kwargs_cp["k1"]:9.7f}_{kwargs_cp["k2"]:9.7f}_{kwargs_cp["rho1"]:9.7f}_{kwargs_cp["rho2"]:9.7f}_{kwargs_cp["cp1"]:9.7f}_{kwargs_cp["cp2"]:9.7f}'
									futures[(x, t, nx, dt, np.random.rand(1)[0])] = executor.submit(self._run_base, x, t, nx, dt, self.tsteps, tshift=tshift, write_freq=wf, **kwargs_cp)
			else:
				for x in self.xs:
					for t in self.ts:
						for nx in self.nxs:
							for dt in self.dts:
								futures[(x, t, nx, dt)] = executor.submit(self._run_base, x, t, nx, dt, self.tsteps, tshift=tshift, write_freq=wf, **kwargs)
			
			while True:
				if len(futures.items()) == 0:
					break
				
				dellist = []
				for key, future in futures.items():
					x, t, nx, dt = key[:4]
					if future.done():
						print(f"Completed config {x=}, {t=}, {nx=} and {dt=}")
						if future.exception() is None:
							dellist.append(key)
						else:
							raise future.exception()

				for key in dellist:
					del futures[key]

				time.sleep(1)

def trial(xs, ts, **kwargs):
	print('in trial')
	print(xs)
	print(ts)
	print(kwargs)