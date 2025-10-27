#!/usr/bin/env python3

import os
import time
import shutil
import subprocess

from datetime import datetime
from concurrent.futures import wait, ProcessPoolExecutor

import numpy as np
from scipy.optimize import newton
from scipy.special import erf

from visualizer import Visualizer

class Run:
	def __init__(self, prog, xs, ts, nxs, dts, tsteps, dxdt, nthreads, dname=None, debug=False, interactive=False, **kwargs):
		self.prog = prog
		self.xs = [xs] if isinstance(xs, int) else xs
		self.ts = [ts] if isinstance(ts, int) else ts
		self.nxs = [nxs] if isinstance(nxs, int) else nxs
		self.dts = [dts] if isinstance(dts, float) else dts
		self.tsteps = tsteps
		self.dxdt = dxdt
		self.nthreads = nthreads
		self.debug = debug
		self.interactive = interactive # flag to avoid extra computations when running in jupyter notebooks
		self.kwargs = kwargs
		self._exes = {}

		self._cleanup = False

		# save current directory for reference
		self.od = os.getcwd()

		# prepare required directories for runs
		if not os.path.exists(f'{self.od}/runs'):
			os.mkdir(f'{self.od}/runs')
		
		if debug:
			self.wd = f'{self.od}/runs/debug/'
		else:
			if dname is None:
				self.wd = f'{self.od}/runs/{datetime.today().strftime("%Y-%m-%d_%H:%M:%S")}'
			else:
				self.wd = f'{self.od}/runs/{dname}'
		
		print(f"Working in {self.wd}")

		# make run directory and set cleanup flag
		if not os.path.exists(self.wd):
			self._cleanup = True
			os.mkdir(self.wd)

		# make remaining directories
		for path in [f'{self.wd}/imgs']:
			if not os.path.exists(path):
				os.mkdir(path)

		if os.path.exists(f'{self.wd}/RESLT'):
			shutil.rmtree(f'{self.wd}/RESLT')
		os.mkdir(f'{self.wd}/RESLT')

		print(f'Run.__init__ {prog=} {self.xs=} {ts=} {nxs=} {dts=} {tsteps=} {kwargs=}')

		stylef = kwargs.pop('stylef', None)
		if not self.interactive:
			# build required executables
			self.build()

			os.chdir(self.wd)

			# perform all computations
			try:
				self.run(**kwargs)
			except:
				os.chdir(self.od)
				if self._cleanup:
					shutil.rmtree(self.wd)
				raise

			# post-process results
			Visualizer(self.prog, self.dxdt, stylef=stylef, reference=kwargs['reference'])
			os.chdir(self.od)

	# not parallelized as all processes will be writing to the same file in main directory
	# unsure how to deal with this race condition yet
	#
	# stdio output from build commands are saved to a file for later reference
	def build(self):
		# we use precompiler flags to set certain options in the code
		# this code block gets the current flags from the Makefile
		with open("Makefile", "r") as f:
			while line := f.readline():
				if line.startswith("CXXFLAGS"):
					flags_base = line.strip()

		for x in self.xs:
			for t in self.ts:
				# set custom pre-compiler flags
				if self.debug:
					flags = f'CXXFLAGS=-Wall -O0 -g -DRUN_SCRIPT -DX_ORDER={x} -DT_ORDER={t}'
				else:
					flags = f'{flags_base} -DRUN_SCRIPT -DX_ORDER={x} -DT_ORDER={t}'

				# delete existing build artefacts to ensure proper compilation
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
					# print stdio output in case of error
					with open(f"{self.wd}/build_{x}_{t}.stdout", "r") as f:
						print(f.read())
					raise

				# save all executables with different names to avoid clashing filenames
				prog_name = f'{self.prog}_{x}_{t}'
				os.rename(self.prog, f'{self.wd}/{prog_name}')
				self._exes[(x, t)] = prog_name

	def get_De(self, k1, k2, rho1, rho2, cp1, cp2, L, **kwargs):
		def func(x):
			D1 = k1/(rho1*cp1)
			D2 = k2/(rho2*cp2)
			Ts = -1.0
			Tm = 0.0
			Tl = 1.0

			e1 = np.emath.sqrt(k1*rho1*cp1)
			e2 = np.emath.sqrt(k2*rho2*cp2)

			match self.prog:
				case 'erf2d11':
					return cp2*(Tm-Ts)/L - 0.5*np.emath.sqrt(np.pi*x/D2) * np.exp(0.25*x/D2) * ((e2/e1) + erf(0.5*x/D2))
				case 'erf2d12':
					alpha = L / (cp1 * (Tm - Ts))
					return 1 - alpha * np.emath.sqrt(0.25 * x * np.pi) * np.exp(0.25 * x) * erf(np.emath.sqrt(0.25 * x))
				case 'erf2d13':
					St = L / (cp1 * (Tl - Ts))
					beta = x / D1
					Tp = (Tl - Tm)/(Tl - Ts)
					er = np.emath.sqrt((k2*rho2*cp2)/(k1*rho1*cp1))
					return 0.5 * St * np.emath.sqrt(np.pi * beta) + (Tp-1.0)*np.exp(-0.25*beta)/erf(0.5*np.emath.sqrt(beta)) + er*Tp*np.exp(-0.25 * x/D2)/(1.0-erf(0.5*np.emath.sqrt(x/D2)))
				case 'erf2d14':
					St = L / (cp1 * (Tm - Ts))
					T_l = (Tm - Tl) / (Tm - Ts)
					D = D1 / D2
					k = k1 / k2
					return 0.5 * St * np.emath.sqrt(np.pi * x) - (k * T_l * np.exp(-0.25 * x/D))/(np.emath.sqrt(D)*(1 - erf(0.5 * np.emath.sqrt(x/D)))) - np.exp(-0.25 * x)/erf(0.5*np.emath.sqrt(x))
				case _:
					return 0.5*rho1*L*np.emath.sqrt(x*np.pi) - k1*(Tm-Ts)*np.exp(-0.25*x/D1)/(np.sqrt(D1)*(1+erf(0.5*np.emath.sqrt(x/D1)))) + k2*(Tl-Tm)*np.exp(-0.25*x/D2)/(np.sqrt(D2)*(1-erf(0.5*np.emath.sqrt(x/D2))))

		
		De = newton(func, k1/(rho1*cp1))

		if (not np.isclose(func(De), 0.0)):
			print("Trying with D2")
			De = newton(func, k2/(rho2*cp2))
		
		De = np.real(De)
			
		print(f'f({De}) = {func(De)}')

		if func(De) is None:
			raise Exception("De is problem")
		return De

	# wrapper for actual command that gets run
	# this method also handles the command line arguments that need to be passed to the executable
	def _run_base(self, x, t, nx, dt, tsteps, tshift=None, write_freq=None, **kwargs):
		nxargs = ('--nx', str(nx))
		dtargs = ('--dt', str(dt)) if dt is not None else ('--vardt')
		tstepsargs = ('--tsteps', str(tsteps))
		tshiftargs = ('--tshift', str(tshift) if tshift is not None else "0.0")
		write_freqargs = ('--write-freq', str(write_freq) if write_freq is not None else "1")

		if 'dname' not in kwargs.keys():
			kwargs['dname'] = f'RESLT/{x}n{nx}_{t}t{dt}'

		# ke = self.get_ke(kwargs['k1'], kwargs['k2'], kwargs['rho1'], kwargs['rho2'], kwargs['cp1'], kwargs['cp2'], kwargs['L'])
		# De = self.get_De(**kwargs)
		# kwargs['De'] = De

		translated_kwargs = []
		for key, value in kwargs.items():
			if key in ['eps'] or value is None:
				continue

			translated_kwargs.extend([f'--{key}', str(value)])# if value is not str else value])
		
		# print(translated_kwargs)
		print(f'Running config {x=}, {t=}, {nx=} and {dt=}')

		# call executable with command line argumets
		try:
			with open(f'{kwargs["dname"]}.stdout', 'w') as f:
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
	
	# method that organises and performs required runs in parallel
	# will run all permutations of options passed to the program, with small adjustments for type of computation
	def run(self, tshift=None, wf=None, **kwargs):
		# calculate max number of CPUs usable/available
		num_workers = len(self.xs) * len(self.ts) * len(self.nxs) * len(self.dts)
		if (num_workers > self.nthreads):
			num_workers = self.nthreads

		kwargs['De'] = self.get_De(**kwargs)
		
		with open('config', 'w') as f:
			f.write(f'xs: {self.xs}\n')
			f.write(f'ts: {self.ts}\n')
			f.write(f'nxs: {self.nxs}\n')
			f.write(f'dts: {self.dts}\n')
			f.write(f'tsteps: {self.tsteps}\n')
			f.write(f'dxdt: {self.dxdt}\n')
			f.write(f'kwargs: {self.kwargs}')

		# run processes in parallel
		futures = {}
		with ProcessPoolExecutor(max_workers=num_workers) as executor:
			# for dt error analysis
			if 't' in self.dxdt:
				max_dt = max(self.dts)
				if tshift is None:
					tshift = (max(self.ts) + 1) * max_dt
				else:
					tshift += (max(self.ts) + 1) * max_dt
				for x in self.xs:
					for t in self.ts:
						for nx in self.nxs:
							for dt in self.dts:
								tsteps = self.tsteps / dt
								wf = max_dt / dt
								futures[(x, t, nx, dt)] = executor.submit(self._run_base, x, t, nx, dt, tsteps, tshift=tshift, write_freq=wf, **kwargs)
			# for sensitivity analysis
			elif 's' in self.dxdt:
				for x in self.xs:
					for t in self.ts:
						for nx in self.nxs:
							for dt in self.dts:
								kwargs['dname'] = f'RESLT/{x}n{nx}_{t}t{dt:.2e}_{kwargs["k1"]:9.7f}_{kwargs["k2"]:9.7f}_{kwargs["rho1"]:9.7f}_{kwargs["rho2"]:9.7f}_{kwargs["cp1"]:9.7f}_{kwargs["cp2"]:9.7f}_{kwargs["L"]:9.7f}'
								futures[(x, t, nx, dt, np.random.rand(1)[0])] = executor.submit(self._run_base, x, t, nx, dt, self.tsteps, tshift=tshift, write_freq=wf, **kwargs)
								for v in ['k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2', 'L']:
									kwargs_cp = kwargs.copy()
									kwargs_cp[v] *= (1.0 + kwargs_cp['eps'])
									kwargs_cp['dname'] = f'RESLT/{x}n{nx}_{t}t{dt:.2e}_{kwargs_cp["k1"]:9.7f}_{kwargs_cp["k2"]:9.7f}_{kwargs_cp["rho1"]:9.7f}_{kwargs_cp["rho2"]:9.7f}_{kwargs_cp["cp1"]:9.7f}_{kwargs_cp["cp2"]:9.7f}'
									kwargs_cp['De'] = self.get_De(**kwargs_cp)
									futures[(x, t, nx, dt, np.random.rand(1)[0])] = executor.submit(self._run_base, x, t, nx, dt, self.tsteps, tshift=tshift, write_freq=wf, **kwargs_cp)
			# for all other run types
			else:
				for x in self.xs:
					for t in self.ts:
						for nx in self.nxs:
							for dt in self.dts:
								futures[(x, t, nx, dt)] = executor.submit(self._run_base, x, t, nx, dt, self.tsteps, tshift=tshift, write_freq=wf, **kwargs)
			
			# Check for process completion
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