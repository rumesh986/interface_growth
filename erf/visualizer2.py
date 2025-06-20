import os
import sys
import asyncio

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import numpy as np
import pandas as pd

import timeit
from visualizer import Visualizer

class _Config:
	def __init__(self, file):
		self.nx = 0
		self.x_order = 0
		self.dt = 0.0
		self.t_shift = 0.0
		self.t_order = 0
		self.t_steps = 0
		self.write_freq = 1

		with open(file) as f:
			while line := f.readline():
				k, v = line.split('=')

				match k:
					case 'nx': self.nx = int(v)
					# case 'nx1': self.nx1 = int(v)
					# case 'nx2': self.nx2 = int(v)
					case 'x_order': self.x_order = int(v)
					case 'dt': self.dt = float(v)
					case 't_shift': self.t_shift = float(v)
					case 't_order': self.t_order = int(v)
					case 't_steps': self.t_steps = int(v)
					case 'write_freq': self.write_freq = int(v)
					# case _: raise Exception(f"Unexpected key in config {k}")

		self._title = f'nx={self.nx}, dt={self.dt}, X({self.x_order}), T({self.t_order})'

	def __str__(self):
		return f'{self.nx = }\n{self.x_order = }\n{self.dt = }\n{self.t_shift = }\n{self.t_order = }\n{self.t_steps = }\n{self.write_freq = }'
	
	@property
	def title(self):
		return self._title
 
class _RunData:
	def __init__(
		self, 
		inp_folder, 
		solnf='soln', 
		exact_solnf='exact_soln', 
		errorf='error',
		configf='config',
		ext='dat'
	):
		self.inp_folder = inp_folder
		self.solnf = solnf
		self.exact_solnf = exact_solnf
		self.errorf = errorf
		self.ext = ext

		self._exact_solns = None
		self._solns = None
		self._errors = None

		self.config = _Config(f'{self.inp_folder}/{configf}')
		self.data = pd.DataFrame()

		try:
			self.times = np.fromfile(f'{self.inp_folder}/times.{self.ext}', sep='\n')
		except FileNotFoundError:
			raise Exception("Time data not availble")
			
	@property
	def title(self):
		return self.config.title
	
	@property
	def errors(self):
		if not self._errors:
			self._errors = self._read_series(self.errorf)
		return self._errors
	
	@property
	def solns(self):
		if not self._solns:
			self._solns = self._read_series(self.solnf)
		return self._solns
	
	@property
	def exact_solns(self):
		if not self._exact_solns:
			self._exact_solns = self._read_series(self.exact_solnf)
		return self._exact_solns
	
	@property
	def error_norms(self):
		abs_errors = np.array([np.linalg.norm(x['error']) for x in self.errors])
		abs_errors /= len(self.errors[0]['error'])

		return abs_errors
	
	@property
	def error_max(self):
		node = np.argmax(np.fabs(self.errors[1]['error'][1:])) + 1
		errors = np.array([df['error'][node] for df in self.errors])
		exact_vals = np.array([df['u'][node] for df in self.exact_solns])
		mask = exact_vals != 0.0

		abs_errors = np.fabs(errors)
		rel_errors = abs_errors[mask] / np.fabs(exact_vals)[mask]

		return abs_errors, rel_errors, self.errors[0]['x'][node]
	
	@property
	def total_error_norm(self):
		errors self.error_norms

		return np.linalg.norm(errors) / len(errors)

	def _read_file(self, ftype, index):
		file = f'{self.inp_folder}/{ftype}{index}.{self.ext}'

		if ftype == self.errorf:
			cols = ['x', 'norm', 'error']
		else:
			cols = ['x', 'u']

		data = []
		with open(file, 'r') as f:
			while line := f.readline():
				if not line.startswith('ZONE'):
					data.append([float(x) for x in line.split()])

		return pd.DataFrame(data, columns=cols)

	def _read_series(self, ftype):
		data = []
		for i, t in enumerate(self.times):
			_data = self._read_file(ftype, i)
			_data['time'] = t

			data.append(_data)
		
		return data

class Visualizer:
	def __init__(self, resd, outd, prefix, dxdt):
		self.resd = resd
		self.outd = outd
		self.prefix = prefix
		self.dxdt = dxdt

		self._data = []
		self._markers = ['+', 'O', '.', '*', 'x', 'D', '^', 'v']
		self._linestyles = ['-', '--', ':']

		plt.style.use('pltstyle.mplstyle')

		for rundir in os.listdir(self.resd):
			rund = f'{self.resd}/{rundir}'
			if not os.path.isdir(rund):
				continue
			
			if not os.path.exists(f'{self.outd}/{rundir}'):
				os.mkdir(f'{self.outd}/{rundir}')

			print(f'Processing results in {rund}')

			self._data.append(_RunData(rund))
		
	
	async def standard(self):
		async def _make(run):
			self.make_anims(run)
			self.plot_errors(run)

		async with asyncio.TaskGroup() as tg:
			for run in self._data:
				tg.create_task(_make(run))

	def plot_errors(self, run):
		plt.cla()

		norm_err = run.error_norms
		abs_max_err, rel_max_err, pos = run.error_max

		fig, (norm_ax, max_ax) = plt.subplots(1,2)
		fig.suptitle(f'Errors ({run.title})')

		norm_ax.semilogy(norm_err)
		
		norm_ax.set_xlabel('Timestep')
		norm_ax.set_ylabel('Norm of Error at timestep')
		norm_ax.set_title('L2 norm of errors')

		max_ax.semilogy(abs_max_err, label='Absolute')
		max_ax.semilogy(rel_max_err, label='Relative')
		
		max_ax.legend()
		max_ax.set_xlabel('Timestep')
		max_ax.set_ylabel('Error magnitude')
		max_ax.set_title(f'Magnitude of errors at x={pos}')

		plt.savefig('trial.png')

	def make_anims(self, run):
		def _update(n):
			time_text.set_text(f't={run.times[n]}')

			line_soln.set_data(solns[n]['x'], solns[n]['u'])
			line_exact_soln.set_data(exact_solns[n]['x'], exact_solns[n]['u'])
			line_diff.set_data(errors[n]['x'], errors[n]['error'])

			diff_ax.set_ylim(errors[n]['error'].min() * 1.1, errors[n]['error'].max() * 1.1)

		plt.cla()

		fig, (diff_ax, prof_ax) = plt.subplots(2, sharex=True, figsize=(8,10))
		fig.suptitle(f'Profiles ({run.title})')

		time_text = prof_ax.annotate(
			f't={run.times[0]}', 
			xy=(0.8, 0.9), 
			xycoords='axes fraction'
		)

		solns = run.solns
		exact_solns = run.exact_solns
		errors = run.errors

		line_soln = prof_ax.plot(
			solns[0]['x'], 
			solns[0]['u'], 
			label='Numerical', 
			lw=5.0
		)[0]
		line_exact_soln = prof_ax.plot(
			exact_solns[0]['x'], 
			exact_solns[0]['u'], 
			label='Analytical', 
			linestyle='--', 
			lw=5.0
		)[0]
		
		line_diff = diff_ax.plot(
			errors[0]['x'], 
			errors[0]['error'], 
		)[0]

		prof_ax.legend(loc='upper left')
		prof_ax.set_ylim([-1.5, 1.5])
		prof_ax.set_ylabel('Temperature-ish')
		prof_ax.set_title('Temperature profile')

		diff_ax.set_xlabel('x')
		diff_ax.set_ylabel('Error')
		diff_ax.set_title('Error in profile')
		diff_ax.yaxis.set_major_formatter('{x:3.1e}')

		anim = FuncAnimation(fig, _update, len(run.times))
		anim.save('trial.mp4')

if __name__ =='__main__':
	vis = Visualizer('TRIALRESLT', 'TRIALIMG', 'erf2', 'a')

	asyncio.run(vis.standard())
	

	# run = _RunData('RESLT/2n1000_4t1.00e-01', 'erf2')

	# print(run.config)
	# rd = run.plot_errors()

	# print(len(rd))

	# # timer1 = timeit.Timer(stmt='for i in range(100): run.read_file("error", i)', setup='run = _RunData("RESLT/2n1000_4t1.00e-01", "erf2")', globals=globals())
	# # print(timer1.timeit(100))

	# vis = Visualizer("RESLT/2n1000_4t1.00e-01", "", "erf2", times="times.dat")
	# vd  = vis.plot_errors()

	# print(rd)
	# print(vd)

	# if np.array_equal(rd, vd):
	# 	print("Were still golden")
	# else:
	# 	print("Damn got issues now")


	# for i in range(100):
	# 	rd = run.read_file('error', i)['error']
	# 	vd = vis.read_file('error', i)[2][:, 2]

	# 	print(f'[{i}] {rd.shape} {vd.shape}')

	# 	if not np.array_equal(rd, vd):
	# 		raise Exception(f'whoops {i}')
			
	# print("Its all good so far")

	# timer2 = timeit.Timer(stmt='for i in range(100): run.read_file("error", i)', setup='run = Visualizer("RESLT/2n1000_4t1.00e-01", "", "erf2", times="times.dat")', globals=globals())
	# print(timer2.timeit(100))