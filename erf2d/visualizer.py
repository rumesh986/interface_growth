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

# Class to store parameters for each run
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

# class to hold results from each run
class _RunData:
	_STEP_HEADERS = ['x', 'y', 'exact', 'u', 'error']
	_OVERALL_HEADERS = ['time', 'error', 'interface', 'expected_interface', 'dx']
	# _OVERALL_HEADERS = ['time', 'error', 'interface']

	def __init__(
		self,
		inp_folder,
		out_folder,
		load_data=True,
		reference=None,
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
			self.interface_errors = np.fabs(self.results['interface'] - self.results['expected_interface'])
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
		else:
			if reference is None:
				i = len(self.results['time']) - 1

				data = pd.read_csv(f'{inp_folder}/{stepsf}{i}.{self.ext}', sep=' ', names=self._STEP_HEADERS)
				data['time'] = self.results['time'].values[-1]
				self.solns.append(data[['time', 'x', 'y', 'u']])
			else:
				step = reference.config.nx // self.config.nx
				count = 0
				print(f"Reference nx = {reference.config.nx} self nx = {self.config.nx} step = {step}")
				
				count = 0
				for t, time in enumerate(reference.times):
					if time not in self.times.values:
						continue
					
					# print(f'[{self.config.nx}] {t =}')
					data = pd.read_csv(f'{inp_folder}/{stepsf}{t}.{self.ext}', sep=' ', names=self._STEP_HEADERS)
					data['time'] = t

					_, _, ref_2d = self._reshape_2D_data(reference.solns[t], 'u')
					_, _, data_2d = self._reshape_2D_data(data, 'u')

					_, _, ref_2dx = self._reshape_2D_data(reference.solns[t], 'x')
					_, _, data_2dx = self._reshape_2D_data(data, 'x')

					errors = np.zeros(data_2d.shape)
					xerrors = np.zeros(data_2d.shape)
					for i in range(data_2d.shape[1]):
						errors[:, i] = data_2d[:, i] - ref_2d[:, i*step]
						xerrors[:, i] = data_2dx[:, i] - ref_2dx[:, i*step]

					print(f'[{self.config.nx}][{t}] errors={np.linalg.norm(errors)} xerrors={np.linalg.norm(xerrors)}')
					self.results['error'][count] = np.linalg.norm(errors)
					count += 1
				
				print(self.results)

	@property
	def title(self):
		return self.config.title
	
	@property
	def error_norms(self):
		return self.results['error']
	
	# get errors at one spatial node
	# node at which errors is max at a certain step is chosen
	@property
	def error_max(self):
		try:
			node = np.argmax(np.fabs(self.errors[10]['error'][1:])) + 1 # node at highest error some steps after start of simulation
		except:
			print("Not enough steps, choosing first available step")
			node = np.argmax(np.fabs(self.errors[0]['error'][1:])) + 1 # node at highest error at start of simulation

		errors = np.array([df['error'][node] for df in self.errors])
		exact = np.array([df['u'][node] for df in self.exact_solns])
		mask = exact != 0.0

		abs_errors = np.fabs(errors)
		rel_errors = abs_errors[mask] / np.fabs(exact[mask])

		return abs_errors, rel_errors, self.errors[0]['x'][node]
	
	# normalized error for whole run
	@property
	def total_error_norm(self):
		return np.linalg.norm(self.results['error'])/len(self.results['error'])

	@property
	def interface_error_norm(self):
		return np.linalg.norm(self.interface_errors) / len(self.interface_errors)

	# create folder to store plots
	def create_outdir(self):
		if not os.path.exists(self.out_folder):
			os.mkdir(self.out_folder)

	# helper method to reshape the data into 2D for easy access
	# function has been optimised to reduce time with pandas operations
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

# Main class to handle post-processing
class Visualizer:
	def __init__(self, prefix, dxdt, resd='RESLT', outd='imgs', stylef='pltstyle.mplstyle', interactive=False, reference=None):
		self.prefix = prefix
		self.dxdt = dxdt
		self.resd = resd
		self.outd = outd
		self.interactive = interactive # flag to avoid extra computations when running in jupyter notebooks
		self.report = stylef == 'report.mplstyle' # check if figures should be optimized for reports

		self.runs = []
		self._markers = ['+', 'o', 'x', '*', '.', 'D', '^', 'v']
		self._linestyles = ['-', '--', ':']

		# set stylesheet for plots
		plt.style.use(f'../../{stylef}')

		self.reference = None
		if reference is not None:
			for rundir in os.listdir(f'../../{reference}/{self.resd}'):
				rund = f'../../{reference}/{self.resd}/{rundir}'
				if os.path.isdir(rund):
					self.reference = _RunData(rund, f'{self.outd}/reference', True)

		# prepare run data
		for rundir in os.listdir(self.resd):
			rund = f'{self.resd}/{rundir}'
			if not os.path.isdir(rund):
				continue

			print(f'Processing results in {rund}')

			self.runs.append(_RunData(rund, f'{self.outd}/{rundir}', 'a' in self.dxdt, reference=self.reference))
		
		# create output directories if needed
		if 'a' in self.dxdt:
			for run in self.runs:
				run.create_outdir()
		
		# perform all post-processing if run in script mode
		if not self.interactive:
			self.prepare_pd()

			if 'a' in self.dxdt:
				self.standard()
			elif 's' in self.dxdt:
				self.sensitivity_analysis()
			else:
				self.plot_interfaces()
				self.plot_analysis()

			# for rtype in self.dxdt:
			# 	match rtype:
			# 		case 'a': self.standard()
			# 		case x if x in 'ibxt': 
			# 			self.plot_interfaces()
			# 			self.plot_analysis(x)
			# 		case 's': self.sensitivity_analysis()

	# collect information about all runs, needed for analysis processing
	def prepare_pd(self):
		self._pd = pd.DataFrame(columns=['x_order', 't_order', 'dx', 'dt', 'error', 'interface_error'])
		config = lambda inp: (inp.x_order, inp.t_order, 1/inp.nx, inp.dt)

		for i, run in enumerate(self.runs):
			x_order, t_order, dx, dt = config(run.config)
			if 'x' in self.dxdt:
				# get largest dx for dx analysis
				# phase 1 of the final step is assumed to have the largest dx of the simulation
				# phase 1 is the growing phase (solid phase)
				# dx = run.solns[-1]['x'][x_order-1] - run.solns[-1]['x'][0]

				# In newer versions, the maximum element size (calculated by oomph-lib) is stored in the results file
				dx = run.results['dx'].values.max()
			self._pd.loc[i] = [x_order, t_order, dx, dt, run.total_error_norm, run.interface_error_norm]

	# post-processing for individual runs (part of 'a' type runs)
	def _make(self, run):
		try:
			# self.make_anims(run)
			self.plot_errors(run)
			# self.plot_error_norms(run)
			# cProfile.runctx("self.make_surf_anims(run)", {"self": self}, {"run": run}, sort='cumtime')
			# self.make_surf_anims(run)
			# self.make_results_surf(run)
			# self.subplot_surf_prof(run)
			if run.config.nx == 10:
				self.subplot_mesh(run)
			if (run.interface is not None):
				self.plot_interface(run)
				self.plot_de(run)
			
			self.make_surf_profile_anim(run)

		except:
			raise Exception(f"Crashing while processing {run.config}")

	# post-processing for individual runs (part of 'a' type runs)
	# runs post-processing for each run in parallel
	def standard(self):
		with ProcessPoolExecutor(max_workers=10) as executor:
			for run, thread in zip(self.runs, executor.map(self._make, self.runs)):
				print(f'Standard processing for \n{run.config}')
		self.plot_interfaces()

	# creates a figure with two subplots
	# one subplot shows how the norm of the total error changes with each timestep
	# second subplot shows how the error changes at a particular spatial point
	def plot_errors(self, run):
		plt.cla()

		norm_err = run.error_norms
		abs_max_err, rel_max_err, pos = run.error_max

		fig, (norm_ax, max_ax) = plt.subplots(1,2, figsize=(10, 5))
		if not self.report:
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
	
	# similar to previous function but only plots how the norm of error changes with each timestep (for reports)
	def plot_error_norms(self, run):
		plt.cla()

		norm_err = run.error_norms

		fig, ax = plt.subplots(1)
		ax.semilogy(norm_err)
		ax.set_xlabel('Timestep')
		ax.set_ylabel('Norm of Error at timestep')

		if not self.report:
			ax.set_title('L2 norm of errors')

		plt.savefig(f'{run.out_folder}/{self.prefix}-error_norms.png')

		if self.interactive:
			plt.show()
			
		plt.close(fig)
	
	# plots the position of the interface against time
	# also does some curve fitting to estimate the value of D_eff
	def plot_interface(self, run):
		# expected functional form for curve fitting
		def f(x, a, b):
			return np.sqrt(x * a) + b

		plt.cla()

		fig, (iface_ax, err_ax) = plt.subplots(2, figsize=(7,10))
		err_ax2 = err_ax.twinx()
		
		iface_ax.plot(run.times, run.interface, label='Interface position')

		try:
			# perform curve-fitting and plot results
			optimal, _ = scopt.curve_fit(f, run.times, run.interface, p0=[1.0, 0.0])
			x_ref = np.linspace(run.times.iloc[0], run.times.iloc[-1])
			y_ref = f(x_ref, optimal[0], optimal[1])

			print(f'{run.title} D_eff = {optimal[0]}')
			iface_ax.plot(x_ref, y_ref, ls='--', label=fr'$y=\sqrt{{{optimal[0]:.6e} * t}} + {optimal[1]:.2f}$')
		except:
			print("Failed to fit h-curve")

		if 'expected_interface' in run.results.columns:
			if 'De' in run.config.params.keys():
				label = fr'analytical $y=\sqrt{{{run.config.params["De"]:.6e} t}}$'
			else:
				label = 'analytical'
			iface_ax.plot(run.results['time'], run.results['expected_interface'], label=label, ls=':')
			abs_plot, = err_ax.plot(run.results['time'], run.interface_errors, label='Absolute Error')
			rel_plot, = err_ax2.plot(run.results['time'], run.interface_errors / run.results['expected_interface'], label='Relative Error', c='C1')
			
			err_ax.legend(handles=[abs_plot, rel_plot])
			
		iface_ax.legend()
		iface_ax.set_xlabel('time')
		iface_ax.set_ylabel('interface posiion h(t)')
		
		# err_ax.set_yscale('log')
		err_ax.set_xlabel('Time')
		err_ax.set_ylabel('Absolute Error')
		err_ax2.set_ylabel('Relative error')

		if not self.report:
			iface_ax.set_title('Interface position with time')
			err_ax.set_title('Error in interface position')

		fig.savefig(f'{run.out_folder}/{self.prefix}-interface.png')

		if self.interactive:
			plt.show()
		
		plt.close()
	
	# plot interfaces from all runs in one plot
	# mostly just test code
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
		
		if 'expected_interface' in self.runs[0].results.columns:
			plt.plot(self.runs[0].results['time'], self.runs[0].results['expected_interface'], label='analytical', ls=':')


		ax.set_xlabel('time')
		ax.set_ylabel('Interface position $h(t)$')

		if not self.report:
			ax.set_title('Interface position with time')

		# move legend outside plot
		box=ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.05, box.width, box.height * 0.95])
		ax.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), ncol=3, prop={'size': 8})

		plt.savefig(f'{self.outd}/{self.prefix}_interfaces.png')

		if self.interactive:
			plt.show()
		
		plt.close()

	def plot_de(self, run):
		def f(x, a, b):
			return np.sqrt(x * a) + b

		optimal, _ = scopt.curve_fit(f, run.times, run.interface, p0=[1.0, 0.0])
		
		De = np.square(run.interface) / run.times
		print(De)
		print(run.config.params['De'])

		plt.plot(De[4:], label='Instantaneous De')
		plt.axhline(run.config.params['De'], c='C1', linestyle='--', label='Analytical')
		plt.axhline(optimal[0], c='C2', linestyle='--', label='Numerical')
		
		plt.xlabel('time step')
		plt.ylabel('$D_e$')
		plt.legend()
		plt.savefig(f'{run.out_folder}/{self.prefix}-de_conv.png')

	# plot temperature profile at different steps in one figure (for report)
	@mpl.rc_context({'font.size': 20})
	def subplot_surf_prof(self, run):
		# function to extract data from middle of simulation (get 1D profile from the 2D data)
		def _extract_data(data):
			ys = data['y'].unique()
			ref_y = ys[ys.shape[0] //2]

			mask = np.isclose(data['y'].values, ref_y)

			return pd.DataFrame(data.values[mask], data.index[mask], data.columns)

		plt.cla()

		fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(14,10), layout="compressed")

		num_frames = len(run.solns) // 4 # stride for timesteps to choose which frames to plot

		for i, ax in enumerate(axs.flatten()):
			data = _extract_data(run.solns[i*num_frames])
			exact_data = _extract_data(run.exact_solns[i*num_frames])

			ax.plot(data['x'], data['u'])
			ax.scatter(exact_data['x'][::10], exact_data['u'][::10], marker='X', c='C1', zorder=2)
			ax.axvline(run.interface[i*num_frames], c='black', ls='--')

			ax.set_title(f't={run.times[i*num_frames]}')

			if i == 0:
				ax.legend(['Numerical', 'Analytical', 'Interface'], loc='upper left')
		
		fig.supxlabel('            X')
		fig.supylabel('Temperature')
		
		if not self.report:
			fig.suptitle('Temperature Profiles')
		
		plt.tight_layout()
		fig.savefig(f'{run.out_folder}/{self.prefix}-surf_prof.png')
		plt.close(fig)
	
	# plot node positions at different steps in one figure (for report)
	@mpl.rc_context({'font.size': 20})
	def subplot_mesh(self, run):
		plt.cla()
		
		fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(14,10), layout="compressed")

		num_frames = len(run.exact_solns) // 4

		mesh_pos = np.zeros((run.config.nx+1, 4))

		for i, ax in enumerate(axs.flatten()):
			# data = run.exact_solns[i]
			data = run.exact_solns[i*num_frames]
			mesh_pos[:, i] = data[data['y'] == 0.0]['x'].values
			print(f'time in col {i}={run.times[i]}')

			ax.scatter(data['x'], data['y'])
			ax.axvline(run.interface[i*num_frames], c='black', ls='--')
			ax.set_title(f't={run.times[i*num_frames]}')
		
		print(mesh_pos)
		
		fig.supxlabel('            X')
		fig.supylabel('Y')
		
		fig.savefig(f'{run.out_folder}/{self.prefix}-meshes.png')
		plt.close(fig)

	# helper method to reshape the data into 2D for easy access
	# function has been optimised to reduce time with pandas operations
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
	
	# plot surface animation of numerical, analytical results and errors
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

	# plot surface animation of only numerical results
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

	# plot cross-section of the 2D surface (easier to visualize for the 1D problem)
	def make_surf_profile_anim(self, run):
		def _extract_data(data):
			ys = data['y'].unique()
			ref_y = ys[ys.shape[0] //2]

			mask = np.isclose(data['y'].values, ref_y)

			return pd.DataFrame(data.values[mask], data.index[mask], data.columns)

		def _update(n):
			time_text.set_text(f't={run.times[n]:6.4f}')

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
			
			if 'expected_interface' in run.results.columns:
				line_interface_exp.set_xdata([run.results['expected_interface'][n]])
				diff_interface_exp.set_xdata([run.results['expected_interface'][n]])

			diff_ax.set_ylim(errors['error'].min() * 1.1, errors['error'].max() * 1.1)
		
		plt.cla()
		
		fig, (prof_ax, diff_ax) = plt.subplots(2, sharex=True, figsize=(8,10))
		# fig, prof_ax= plt.subplots(1, figsize=(10,5))

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
		
		if 'expected_interface' in run.results.columns:
			line_interface_exp = prof_ax.axvline(run.results['expected_interface'][0], c='blue', ls=':', label='Expected interface')
			diff_interface_exp = diff_ax.axvline(run.results['expected_interface'][0], c='blue', ls='--', label='expected interface')
		
		for pos in run.config.fixed_pos:
			prof_ax.axvline(pos, c='black', ls='--', label='fixed interface')

		prof_ax.axhline(0.0, c='blue', ls=':', label='Expected interface temp')

		prof_ax.set_ylim([-1.5, 1.5])
		prof_ax.set_xlim([data['x'].min()-0.1, data['x'].max()])
		prof_ax.set_ylabel('Temperature-ish')

		if self.report:
			box = prof_ax.get_position()
			prof_ax.legend(loc='center left', bbox_to_anchor=(1.0,0.5), ncol=1)
			plt.tight_layout()
		else:
			prof_ax.legend(loc='upper left')
			prof_ax.set_title('Temperature profile')

		diff_ax.set_xlabel('x')
		diff_ax.set_ylabel('Error')
		diff_ax.set_title('Error in profile')
		diff_ax.yaxis.set_major_formatter('{x:3.1e}')

		anim = FuncAnimation(fig, _update, len(run.times))
		anim.save(f'{run.out_folder}/{self.prefix}-surf_prof.mp4')	
		plt.close(fig)

	# plot 1D data (for older erf code that was run with 1D elements)
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

	# perform and plot error analysis
	def plot_analysis(self):
		# self.plot_interfaces()
		# helper function to get legend entries
		def _legend(x, t):
			match x:
				case 2: return f"Linear BDF{t}"
				case 3: return f"Quadratic BDF{t}"
				case 4: return f"Cubic BDF{t}"
				case _: raise Exception("Unknown element order")
		
		if 'x' in self.dxdt:
			xlabel = 'dx'
			ylabel = 'error'
			title = f'{self.prefix} dx error analysis'
			fname = 'dx_errors'
			ref_order = 'x_order'
			ref_variable = 'x'
			xs = np.logspace(0, -3.5)
		elif 't' in self.dxdt:
			xlabel = 'dt'
			ylabel = 'error'
			title = f'{self.prefix} dt error analysis'
			fname = 'dt_errors'
			ref_order = 't_order'
			ref_variable = 't'
			xs = np.logspace(0, -4)
		
		if 'i' in self.dxdt:
			ylabel = 'interface_error'
			title += ' (interface)'
			fname += '_iface'

		x_orders = self._pd['x_order'].unique()
		x_orders.sort()

		t_orders = self._pd['t_order'].unique()
		t_orders.sort()

		fig, ax = plt.subplots(1, subplot_kw={
			'xlabel': xlabel,
			'ylabel': 'Normalized Error',
			'xscale': 'log',
			'yscale': 'log',
			'ylim': [self._pd[ylabel].min() * 1e-1, self._pd[ylabel].max() * 1e1],
			'xlim': [self._pd[xlabel].min() * 0.8, self._pd[xlabel].max() * 1.2],
		})

		marker_i = 0
		for x in x_orders:
			for t in t_orders:
				df = self._pd.loc[self._pd['x_order'] == x]
				df = df.loc[df['t_order'] == t]

				ax.scatter(df[xlabel], df[ylabel], label=_legend(int(x), int(t)), marker=self._markers[marker_i % len(self._markers)])
				marker_i += 1

		# plot reference lines
		# lines are made to intersect the point at the largest dx/dt
		ref_orders = self._pd[ref_order].unique()
		ref_orders.sort()
		for x in ref_orders:
			df = self._pd.loc[self._pd[ref_order] == x]

			if 'i' in self.dxdt:
				ref_power = x - 0.8
			else:
				ref_power = x

			if 't' in self.dxdt:
				ref_power = 1.2

			# max_point = df.loc[df[ylabel] == df[ylabel].min()]
			max_point = df.iloc[(df[xlabel] - 0.02).abs().argsort()[:1]]
			ax.plot(xs, (xs / max_point[xlabel].iloc[0]) ** ref_power * max_point[ylabel].iloc[0], label=fr'$\mathcal{{O}}({ref_variable}^{{{ref_power:.2f}}})$', linestyle='--', marker='', alpha=0.8)

		if not self.report:
			fig.suptitle(title)

		# put legend outside axis
		box=ax.get_position()
		ax.legend(loc='center left', bbox_to_anchor=(1.0,0.5), ncol=1)
		plt.tight_layout()

		fig.savefig(f'{self.outd}/{self.prefix}_{fname}.png')

		if self.interactive:
			plt.show()

		plt.close(fig)

		print(self._pd)
	
	# plot results of sensitivity analysis
	def sensitivity_analysis(self):
		def f(x, a, b):
			return np.sqrt(x*a) + b
		
		eps = 1e-6
		table = pd.DataFrame(columns=['key', 'k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2', 'L', 'deff', 'sensitivity'])
		for i, run in enumerate(self.runs):
			try:
				# calculate D_eff (quantity of interest in sensitivity analysis)
				optimal, _ = scopt.curve_fit(f, run.times, run.interface, p0=[1.0, 0.0])
				params = run.config.params
				table.loc[i] = ['', params['k1'], params['k2'], params['rho1'], params['rho2'], params['cp1'], params['cp2'], params['L'], optimal[0], 0.0]
			except:
				print("Something went wrong")
				raise
		pd.set_option("display.precision", 10)

		keys = ['k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2', 'L']
		orig = {}
		for k in keys:
			orig[k] = table[k].min()

		# find run with mean parameter values (not shifted to calculate sensitivity)
		orig_row = table[
			(table['k1'] == orig['k1']) & 
			(table['k2'] == orig['k2']) & 
			(table['rho1'] == orig['rho1']) & 
			(table['rho2'] == orig['rho2']) & 
			(table['cp1'] == orig['cp1']) &
			(table['cp2'] == orig['cp2']) &
			(table['L'] == orig['L'])
		]

		# proper display names for each parameter
		key_names = {
			'k1': r'$\lambda_s$',
			'k2': r'$\lambda_l$',
			'rho1': r'$\rho_s$',
			'rho2': r'$\rho_l$',
			'cp1': r'$C_{p_s}$',
			'cp2': r'$C_{p_l}$',
			'L': r'$\mathcal{L}_f$'
		}

		# calculate sensitivities
		for k in keys:
			index = table.index[table[k] == table[k].max()]
			h = table.loc[index, k] - orig_row[k].iloc[0]

			table.loc[index, 'key'] = key_names[k]
			table.loc[index, 'sensitivity'] = (table.loc[index, 'deff'] - orig_row['deff'].iloc[0]) / h
			table.loc[index, 'scaled_sensitivity'] = table.loc[index, k] * table.loc[index, 'sensitivity']

		# filter and sort data for plotting
		filtered_table = pd.DataFrame(columns=['key', 'k1', 'k2', 'rho1', 'rho2', 'cp1', 'cp2', 'L', 'deff', 'sensitivity', 'scaled_sensitivity'])
		for i, k in enumerate(key_names.keys()):
			filtered_table.loc[i] = table[table['key'] == key_names[k]].iloc[0]

		plt.plot(filtered_table['key'], filtered_table['scaled_sensitivity'].abs())
		plt.xlabel('Parameters')
		plt.ylabel('Scaled sensitivities')

		if not self.report:
			plt.title(r"Scaled sensitivities of $D_{\mathrm{eff}}$")
	
		plt.savefig(f'{self.outd}/{self.prefix}-sensitivity.png')

# simple way to run visualizer for debug purposes outside of proper runs
if __name__ =='__main__':
	# assert len(sys.argv) == 3
	# assert sum(c1 == c2 for c1 in "abxts" for c2 in sys.argv[2]) > 0
	# assert sum(c1 == c2 for c1 in "bxt" for c2 in sys.argv[2]) < 2

	os.chdir(sys.argv[1])
	# vis = Visualizer(sys.argv[2], sys.argv[3], stylef='report.mplstyle')
	vis = Visualizer(sys.argv[2], sys.argv[3])
	# vis = Visualizer(sys.argv[2], 'x', reference='runs/reference')
