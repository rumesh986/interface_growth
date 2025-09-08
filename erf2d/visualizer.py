#!/usr/bin/env python3

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, ArtistAnimation

import numpy as np
import pandas as pd
import scipy.optimize as scopt

# import cProfile

class _Config:
	def __init__(self, file):
		self.nx = 0
		self.x_order = 0
		self.dt = 0.0
		self.t_shift = 0.0
		self.t_order = 0
		self.t_steps = 0
		self.write_freq = 1

		self.params = {}
		self.fixed_pos = []

		with open(file) as f:
			while line := f.readline():
				k, v = line.split('=')

				match k:
					case 'nx': self.nx = int(v)
					case 'x_order': self.x_order = int(v)
					case 'dt': self.dt = float(v)
					case 't_shift': self.t_shift = float(v)
					case 't_order': self.t_order = int(v)
					case 't_steps': self.t_steps = int(v)
					case 'write_freq': self.write_freq = int(v)
					case x if x.startswith('x'): self.fixed_pos.append(float(v))
					case _: self.params[k] = float(v)

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
		out_folder,
		solnf='solns/soln', 
		exact_solnf='exact_solns/exact_soln', 
		errorf='errors/error',
		configf='config',
		ext='dat'
	):
		self.inp_folder = inp_folder
		self.out_folder = out_folder
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
		
		self.interface = None
		if os.path.exists(f'{self.inp_folder}/interface.{self.ext}'):
			self.interface = np.fromfile(f'{self.inp_folder}/interface.{self.ext}', sep='\n')
			
	@property
	def title(self):
		return self.config.title
	
	@property
	def errors(self):
		if not self._errors:
			self._errors = self._read_series(self.errorf)#, ignore_history=True)
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
		# skip_steps = int(self.config.t_shift / self.config.dt)+1
		skip_steps = self.config.t_order + 1

		abs_errors = np.fabs(errors)
		# rel_errors = abs_errors[mask[skip_steps:]] / np.fabs(exact_vals)[mask][skip_steps:]
		rel_errors = abs_errors[mask] / np.fabs(exact_vals[mask])

		return abs_errors, rel_errors, self.errors[0]['x'][node]
	
	@property
	def total_error_norm(self):
		errors = self.error_norms

		return np.linalg.norm(errors) / len(errors)

	def _read_file(self, ftype, index):
		file = f'{self.inp_folder}/{ftype}{index}.{self.ext}'

		if ftype == self.errorf:
			cols = ['x', 'y', 'norm', 'error']
		else:
			cols = ['x', 'y', 'u']

		data = []
		with open(file, 'r') as f:
			while line := f.readline():
				if not line.startswith('ZONE'):
					data.append([float(x) for x in line.split()])

		return pd.DataFrame(data, columns=cols)

	def _read_series(self, ftype, ignore_history=False):
		data = []
		for i, t in enumerate(self.times):
			# ignore data set in initial condition (including history values)
			if ignore_history and t <= self.config.t_shift:
				continue
			_data = self._read_file(ftype, i)
			_data['time'] = t

			data.append(_data)
		
		return data

	def create_outdir(self):
		if not os.path.exists(self.out_folder):
			os.mkdir(self.out_folder)

class _RunData2:
	_STEP_HEADERS = ['x', 'y', 'exact', 'u', 'error']
	_OVERALL_HEADERS = ['time', 'error', 'interface']

	def __init__(
		self,
		inp_folder,
		out_folder,
		load_data=True,
		stepsf='steps/step',
		configf='config',
		resultsf='results',
		ext='dat'
	):
		self.inp_folder = inp_folder
		self.out_folder = out_folder
		self.stepsf = stepsf
		self.resultsf = resultsf
		self.ext = ext
		
		self.config = _Config(f'{self.inp_folder}/{configf}')

		try:
			self.results = pd.read_csv(f'{self.inp_folder}/{self.resultsf}.{self.ext}', sep=' ', names=self._OVERALL_HEADERS)	
			self.times = self.results['time']
			self.interface = self.results['interface']
		except FileNotFoundError as e:
			print("Results file not found")
			raise e
		
		self.data = []
		self.errors = []
		self.exact_solns = []
		self.solns = []

		if load_data:
			for i, t in enumerate(self.results['time']):
				data = pd.read_csv(f'{inp_folder}/{stepsf}{i}.{self.ext}', sep=' ', names=self._STEP_HEADERS)
				data['time'] = t

				self.errors.append(data[['time', 'x', 'y', 'error']])
				self.exact_solns.append(data[['time', 'x', 'y', 'exact']].rename(columns={'exact': 'u'}))
				self.solns.append(data[['time', 'x', 'y', 'u']])

				self.data.append(data)

	@property
	def title(self):
		return self.config.title
	
	@property
	def error_norms(self):
		return self.results['error']
		# abs_errors = np.array([np.linalg.norm(x['error']) for x in self.errors])
		# abs_errors /= len(self.errors[0]['error'])

		# if np.all(np.equal(abs_errors, self.results['error'])):
		# 	print("Yup its all equal")
		# else:
		# 	print("There appears to be some difference")
		# 	print(f"Mag: {np.linalg.norm(abs_errors - self.results['error'])}")
		# 	print(np.c_[self.results['error'], abs_errors])
		# return abs_errors
	
	@property
	def error_max(self):
		node = np.argmax(np.fabs(self.errors[10]['error'][1:])) + 1
		errors = np.array([df['error'][node] for df in self.errors])
		exact = np.array([df['u'][node] for df in self.exact_solns])
		mask = exact != 0.0

		abs_errors = np.fabs(errors)
		rel_errors = abs_errors[mask] / np.fabs(exact[mask])

		return abs_errors, rel_errors, self.errors[0]['x'][node]
	
	@property
	def total_error_norm(self):
		return np.linalg.norm(self.results['error'])/len(self.results['error'])

	def create_outdir(self):
		if not os.path.exists(self.out_folder):
			os.mkdir(self.out_folder)

class Visualizer:
	def __init__(self, prefix, dxdt, resd='RESLT', outd='imgs', stylef='pltstyle.mplstyle', interactive=False):
		self.prefix = prefix
		self.dxdt = dxdt
		self.resd = resd
		self.outd = outd
		self.interactive = interactive

		self.runs = []
		self._markers = ['+', 'o', 'x', '*', '.', 'D', '^', 'v']
		self._linestyles = ['-', '--', ':']

		plt.style.use(f'../../{stylef}')

		for rundir in os.listdir(self.resd):
			rund = f'{self.resd}/{rundir}'
			if not os.path.isdir(rund):
				continue

			print(f'Processing results in {rund}')

			self.runs.append(_RunData2(rund, f'{self.outd}/{rundir}', self.dxdt != 's'))
		
		if 'a' in self.dxdt:
			for run in self.runs:
				run.create_outdir()
		
		if not self.interactive:
			self.prepare_pd()
			for rtype in self.dxdt:
				match rtype:
					case 'a': self.standard()
					case x if x in 'bxt': self.plot_analysis(x)
					case 's': self.sensitivity_analysis()

	def prepare_pd(self):
		self._pd = pd.DataFrame(columns=['x_order', 't_order', 'dx', 'dt', 'error'])
		config = lambda inp: (inp.x_order, inp.t_order, 1/inp.nx, inp.dt)

		for i, run in enumerate(self.runs):
			x_order, t_order, dx, dt = config(run.config)
			self._pd.loc[i] = [x_order, t_order, dx, dt, run.total_error_norm]

	def _make(self, run):
		try:
			# self.make_anims(run)
			# self.plot_errors(run)
			self.plot_error_norms(run)
			# cProfile.runctx("self.make_surf_anims(run)", {"self": self}, {"run": run}, sort='cumtime')
			# self.make_surf_anims(run)
			# self.make_results_surf(run)
			# self.make_surf_profile_anim(run)
			self.subplot_surf_prof(run)
			if (run.interface is not None):
				self.plot_interface(run)

		except:
			raise Exception(f"Crashing for some reason {run.config}")

	def standard(self):
		with ProcessPoolExecutor(max_workers=10) as executor:
			for run, thread in zip(self.runs, executor.map(self._make, self.runs)):
				print(f'Standard processing for \n{run.config}')
		self.plot_interfaces()

	def plot_errors(self, run):
		plt.cla()

		norm_err = run.error_norms
		abs_max_err, rel_max_err, pos = run.error_max

		fig, (norm_ax, max_ax) = plt.subplots(1,2, figsize=(10, 5))
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

		plt.savefig(f'{run.out_folder}/{self.prefix}-errors.png')

		if self.interactive:
			plt.show()
			
		plt.close(fig)
	
	def plot_error_norms(self, run):
		plt.cla()

		norm_err = run.error_norms
		# abs_max_err, rel_max_err, pos = run.error_max

		fig, ax = plt.subplots(1)
		# fig.suptitle(f'Errors ({run.title})')

		ax.semilogy(norm_err)
		print(f'{run.title=} max error={norm_err.max()} at {norm_err.argmax()}')
		
		ax.set_xlabel('Timestep')
		ax.set_ylabel('Norm of Error at timestep')
		ax.set_title('L2 norm of errors')

		# max_ax.semilogy(abs_max_err, label='Absolute')
		# max_ax.semilogy(rel_max_err, label='Relative')
		
		# max_ax.legend()
		# max_ax.set_xlabel('Timestep')
		# max_ax.set_ylabel('Error magnitude')
		# max_ax.set_title(f'Magnitude of errors at x={pos}')

		plt.savefig(f'{run.out_folder}/{self.prefix}-error_norms.png')

		if self.interactive:
			plt.show()
			
		plt.close(fig)
	
	def plot_interface(self, run):
		def f(x, a, b):
			return np.sqrt(x * a) + b

		plt.cla()

		fig, ax = plt.subplots(1)

		ax.plot(run.times, run.interface, label='Interface position')

		try:
			optimal, _ = scopt.curve_fit(f, run.times, run.interface, p0=[1.0, 0.0])
			x_ref = np.linspace(run.times.iloc[0], run.times.iloc[-1])
			y_ref = f(x_ref, optimal[0], optimal[1])

			plt.plot(x_ref, y_ref, ls='--', label=fr'$y=\sqrt{{{optimal[0]:.2f} * t}} + {optimal[1]:.2f}$')
		except:
			print("Failed to fit h-curve")

		ax.legend()
		ax.set_xlabel('time')
		ax.set_ylabel('interface posiion h(t)')
		ax.set_title('Interface position with time')

		fig.savefig(f'{run.out_folder}/{self.prefix}-interface.png')

		if self.interactive:
			plt.show()
		
		plt.close()
	
	def plot_interfaces(self):
		def f(x, a, b):
			return np.sqrt(x*a) + b
		
		plt.cla()

		fig, ax = plt.subplots(1,1)


		for run in self.runs:
			ax.plot(run.times, run.interface, label=run.title)
		
		try:
			optimal, _ = scopt.curve_fit(f, self.runs[0].times, self.runs[0].interface, p0=[1.0, 0.0])
			x_ref = np.linspace(self.runs[0].times.iloc[0], self.runs[0].times.iloc[-1])
			y_ref = f(x_ref, optimal[0], optimal[1])

			ax.plot(x_ref[::3], y_ref[::3], marker=self._markers[0], ls='', ms=10, label=fr'$y=\sqrt{{{optimal[0]:.4f} * t}} + {optimal[1]:.4f}$')
		except:
			print("Failed to fit interface curve")

		ax.set_xlabel('time')
		ax.set_ylabel('Interface position $h(t)$')
		ax.set_title('Interface position with time')

		box=ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
		ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), ncol=3, prop={'size': 8})

		plt.savefig(f'{self.outd}/{self.prefix}_interfaces.png')

		if self.interactive:
			plt.show()
		
		plt.close()
	
	@mpl.rc_context({'font.size': 20})
	def subplot_surf_prof(self, run):
		def _extract_data(data):
			ys = data['y'].unique()
			ref_y = ys[ys.shape[0] //2]

			mask = np.isclose(data['y'].values, ref_y)

			return pd.DataFrame(data.values[mask], data.index[mask], data.columns)

		plt.cla()

		fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(14,10), layout="compressed")

		num_frames = len(run.solns) // 4

		for i, ax in enumerate(axs.flatten()):
			data = _extract_data(run.solns[i*num_frames])
			exact_data = _extract_data(run.exact_solns[i*num_frames])

			ax.plot(data['x'], data['u'])
			ax.scatter(exact_data['x'][::10], exact_data['u'][::10], marker='X', c='C1', zorder=2)
			ax.axvline(run.interface[i*num_frames], c='black', ls='--')

			ax.set_title(f't={run.times[i*num_frames]}')

			if i == 3:
				ax.legend(['Numerical', 'Analytical', 'Interface'], loc='lower right')
		
		fig.supxlabel('            X')
		fig.supylabel('Temperature')
		
		fig.suptitle('Temperature Profiles')
		plt.tight_layout()
		fig.savefig(f'{run.out_folder}/{self.prefix}-surf_prof.png')
		plt.close(fig)

	def _reshape_2D_data(self, data, key):
		xs = data['x'].unique()
		ys = data['y'].unique()
		u = np.zeros((len(ys), len(xs)))
		
		# the loop below has some odd pandas operations, this was done to optimize it
		# pandas seems to be pretty slow and you can get a significant speedup by changing how you interact with it
		for j, y in enumerate(ys):
			# optimized way of filtering dataframe
			mask1 = data['y'].values == y
			d2 = pd.DataFrame(data.values[mask1], data.index[mask1], data.columns)
			for i, x in enumerate(xs):
				# access by index to avoid creating another dataframe
				mask2 = np.isclose(d2['x'].values, x)
				# print(f'{j} {i} {mask2}')
				try:
					index = mask2.nonzero()[0][0]
				except IndexError:
					print(f"{j=} {i=}")
					raise Exception("Its failing again")

				u[j,i] = d2[key].iloc[index]
		
		return xs, ys, u
	
	def make_surf_anims(self, run):
		def _plot_surface(ax, data, key):
			xs, ys, us = self._reshape_2D_data(data, key)
			X, Y = np.meshgrid(xs, ys)
			return ax.plot_surface(X, Y, us, color='b')

		plt.cla()

		fig, (soln_ax, exact_soln_ax, error_ax) = plt.subplots(3, figsize=(8, 15), subplot_kw={'projection': '3d'})

		soln_ax.view_init(elev=0, azim=-90, roll=0)
		exact_soln_ax.view_init(elev=0, azim=-90, roll=0)
		error_ax.view_init(elev=0, azim=-90, roll=0)

		artists = []

		for t, time in enumerate(run.times):
			print(t)
			soln_artist = _plot_surface(soln_ax, run.solns[t], 'u')
			exact_soln_artist = _plot_surface(exact_soln_ax, run.exact_solns[t], 'u')
			error_artist = _plot_surface(error_ax, run.errors[t], 'error')

			time_text = soln_ax.annotate(
				f't={time}',
				xy=(0.8, 0.9),
				xycoords='axes fraction'
			)

			artists.append([soln_artist, exact_soln_artist, error_artist, time_text])
		
		soln_ax.set_xlabel('x')
		soln_ax.set_ylabel('y')
		soln_ax.set_zlabel('Temperature-ish')
		soln_ax.set_title('Simulation')

		exact_soln_ax.set_xlabel('x')
		exact_soln_ax.set_ylabel('y')
		exact_soln_ax.set_zlabel('Temperature-ish')
		exact_soln_ax.set_title('Analytical')

		error_ax.set_xlabel('x')
		error_ax.set_ylabel('y')
		error_ax.set_zlabel('Error')
		error_ax.set_title('Errors')

		anim = ArtistAnimation(fig, artists)
		anim.save(f'{run.out_folder}/{self.prefix}-results.mp4')

	def make_results_surf(self, run):
		plt.cla()

		fig, ax = plt.subplots(1, subplot_kw={'projection': '3d'})	
		ax.view_init(elev=0, azim=-90, roll=0)
		artists = []
		for t, time in enumerate(run.times):
			xs, ys, us = self._reshape_2D_data(run.solns[t], 'u')
			X, Y = np.meshgrid(xs, ys)
			soln_artist = ax.plot_surface(X, Y, us, color='b')

			time_text = ax.annotate(
				f't={time}',
				xy=(0.8, 0.9),
				xycoords='axes fraction'
			)

			artists.append([soln_artist, time_text])
		
		ax.set_xlabel('x')
		ax.set_ylabel('y')
		ax.set_zlabel('Temperature-ish')
		ax.set_title('simulation')

		anim = ArtistAnimation(fig, artists)
		anim.save(f'{run.out_folder}/{self.prefix}-soln.mp4')

	def make_surf_profile_anim(self, run):
		def _extract_data(data):
			ys = data['y'].unique()
			ref_y = ys[ys.shape[0] //2]

			mask = np.isclose(data['y'].values, ref_y)

			return pd.DataFrame(data.values[mask], data.index[mask], data.columns)

		def _update(n):
			time_text.set_text(f't={run.times[n]}')

			data = _extract_data(run.solns[n])
			exact_data = _extract_data(run.exact_solns[n])
			errors = _extract_data(run.errors[n])

			# account for repeated points 
			# 	avoids a line crossing the plot unnecessarily
			num_points = len(data['x'])
			if len(data['x']) != len(data['x'].unique()):
				num_points //= 2

			line_soln.set_data(data['x'][:num_points], data['u'][:num_points])
			line_exact.set_data(exact_data['x'][:num_points], exact_data['u'][:num_points])
			line_diff.set_data(errors['x'][:num_points], errors['error'][:num_points])

			if run.interface is not None:
				line_interface.set_xdata([run.interface[n]])
				diff_interface.set_xdata([run.interface[n]])

			diff_ax.set_ylim(errors['error'].min() * 1.1, errors['error'].max() * 1.1)
		
		plt.cla()
		
		fig, (prof_ax, diff_ax) = plt.subplots(2, sharex=True, figsize=(8,10))

		time_text = prof_ax.annotate(
			f't={run.times[0]}',
			xy=(0.8,0.9),
			xycoords='axes fraction'
		)

		data = _extract_data(run.solns[0])
		exact_data = _extract_data(run.exact_solns[0])
		errors = _extract_data(run.errors[0])

		line_soln = prof_ax.plot(data['x'], data['u'], label='Numerical')[0]
		line_exact = prof_ax.plot(exact_data['x'], exact_data['u'], label='Analytical', ls=':')[0]
		line_diff = diff_ax.plot(errors['x'], errors['error'], label='error')[0]

		if run.interface is not None:
			line_interface = prof_ax.axvline(run.interface[0], c='black', ls='--', label='interface')
			diff_interface = diff_ax.axvline(run.interface[0], c='black', ls='--', label='interface')
		
		for pos in run.config.fixed_pos:
			prof_ax.axvline(pos, c='black', ls='--', label='fixed interface')

		prof_ax.axhline(0.0, c='blue', ls=':', label='Expected interface temp')

		prof_ax.legend(loc='upper left')
		prof_ax.set_ylim([-1.5, 1.5])
		prof_ax.set_ylabel('Temperature-ish')
		prof_ax.set_title('Temperature profile')

		diff_ax.set_xlabel('x')
		diff_ax.set_ylabel('Error')
		diff_ax.set_title('Error in profile')
		diff_ax.yaxis.set_major_formatter('{x:3.1e}')

		anim = FuncAnimation(fig, _update, len(run.times))
		anim.save(f'{run.out_folder}/{self.prefix}-surf_prof.mp4')	
		plt.close(fig)

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
		)[0]
		line_exact_soln = prof_ax.plot(
			exact_solns[0]['x'], 
			exact_solns[0]['u'], 
			label='Analytical', 
			linestyle='--', 
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
		anim.save(f'{run.out_folder}/{self.prefix}-results.mp4')	
		plt.close(fig)

	def plot_analysis(self, rtype):
		match rtype:
			case 'x':
				xlabel = 'dx'
				title = f'{self.prefix} dx error analysis'
				fname = 'dx_errors'
				label = 'order'
				ref_order = 'x_order'
			case 'b':
				xlabel = 'dx'
				title = f'{self.prefix} dxdt error analysis'
				fname = 'dxdt_errors'
				label = 'order'
				ref_order = 'x_order'
			case 't':
				xlabel = 'dt'
				title = f'{self.prefix} dt error analysis'
				fname = 'dt_errors'
				label = 'BDF'
				ref_order = 't_order'
			case _:
				raise Exception("Unknown analysis type provided")

		xs = np.logspace(-1, -4)
		x_orders = self._pd['x_order'].unique()
		x_orders.sort()

		t_orders = self._pd['t_order'].unique()
		t_orders.sort()

		fig, ax = plt.subplots(1, subplot_kw={
			'xlabel': xlabel,
			'ylabel': 'Normalized Error',
			'xscale': 'log',
			'yscale': 'log',
			'ylim': [self._pd['error'].min() * 1e-1, self._pd['error'].max() * 1e1],
			'title': title
		})

		marker_i = 0
		for x in x_orders:
			for t in t_orders:
				df = self._pd.loc[self._pd['x_order'] == x]
				df = df.loc[df['t_order'] == t]

				ax.scatter(df[xlabel], df['error'], label=f'X({x}) T({t})', marker=self._markers[marker_i % len(self._markers)])
				marker_i += 1

		ref_orders = self._pd[ref_order].unique()
		ref_orders.sort()

		for x in ref_orders:
			df = self._pd.loc[self._pd[ref_order] == x]

			max_point = df.loc[df['error'] == df['error'].max()]
			ax.plot(xs, (xs / max_point[xlabel].iloc[0]) ** x * max_point['error'].iloc[0], label=f'Reference (order={x})', linestyle='--', marker='')

		ax.legend(ncols=ref_orders.shape[0])

		fig.savefig(f'{self.outd}/{self.prefix}_{fname}.png')

		if self.interactive:
			plt.show()

		plt.close(fig)
	
	def sensitivity_analysis(self):
		def f(x, a, b):
			return np.sqrt(x*a) + b
		
		eps = 1e-6
		table = pd.DataFrame(columns=['key', 'k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2', 'deff', 'sensitivity'])
		for i, run in enumerate(self.runs):
			try:
				optimal, _ = scopt.curve_fit(f, run.times, run.interface, p0=[1.0, 0.0])
				params = run.config.params
				table.loc[i] = ['', params['k1'], params['k2'], params['rho1'], params['rho2'], params['cp1'], params['cp2'], optimal[0], 0.0]
			except:
				print("Something went wrong")
				raise
		pd.set_option("display.precision", 10)

		keys = ['k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2']
		orig = {}
		for k in keys:
			orig[k] = table[k].min()

		orig_row = table[
			(table['k1'] == orig['k1']) & 
			(table['k2'] == orig['k2']) & 
			(table['rho1'] == orig['rho1']) & 
			(table['rho2'] == orig['rho2']) & 
			(table['cp1'] == orig['cp1']) &
			(table['cp2'] == orig['cp2'])
		]

		key_names = {
			'k1': r'$\lambda_s$',
			'k2': r'$\lambda_l$',
			'rho1': r'$\rho_s$',
			'rho2': r'$\rho_l$',
			'cp1': r'$C_{p_s}$',
			'cp2': r'$C_{p_l}$',
		}

		sensitivity = {}
		for k in keys:
			index = table.index[table[k] == table[k].max()]
			h = table.loc[index, k] - orig_row[k].iloc[0]

			table.loc[index, 'key'] = key_names[k]
			table.loc[index, 'sensitivity'] = (table.loc[index, 'deff'] - orig_row['deff'].iloc[0]) / h
			table.loc[index, 'scaled_sensitivity'] = table.loc[index, k] * table.loc[index, 'sensitivity']

		print(table)

		filtered_table = table[table['key'] != '']

		plt.plot(filtered_table['key'], filtered_table['scaled_sensitivity'].abs())
		plt.xlabel('Parameters')
		plt.ylabel('Scaled sensitivities')
		plt.title(r"Scaled sensitivities of $D_{\mathrm{eff}}$")
		plt.savefig(f'{self.outd}/{self.prefix}-sensitivity.png')

if __name__ =='__main__':
	assert len(sys.argv) == 3
	assert sum(c1 == c2 for c1 in "abxts" for c2 in sys.argv[2]) > 0
	assert sum(c1 == c2 for c1 in "bxt" for c2 in sys.argv[2]) < 2

	vis = Visualizer(sys.argv[1], sys.argv[2])
