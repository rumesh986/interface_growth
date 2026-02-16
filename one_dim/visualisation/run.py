#!/usr/bin/env python3

import os
from typing import Self

import numpy as np
import polars as pl
from scipy.optimize import curve_fit
from matplotlib.animation import FuncAnimation

from systems.systems import Params, AnalysisType, System

class Run:
	def __init__(
		self,
		inp_dir: str,
		out_dir: str,
		analysis_type: AnalysisType = AnalysisType.standard,
		reference: Self = None,
		stepsf: str = 'steps/step',
		resultsf: str = 'results',
		ext: str = 'dat'
	):
		self.inp_dir = inp_dir
		self.out_dir = out_dir
		self.params = Params.from_file(f'{self.inp_dir}/params')
		self.system = System.from_file('system')
		self.reference = reference

		self._De = None

		self.results = pl.read_csv(f'{self.inp_dir}/{resultsf}.{ext}', separator=' ').with_columns(
				abs_error_h=(pl.col('exact_h') - pl.col('h')).abs(),
				rel_error_h=((pl.col('exact_h') - pl.col('h')) / pl.col('h') ).abs()
			).with_row_index('step')
		self.times = self.results['time'] 

		self.start_index = self.params.t+1

		size = self.results[self.start_index:, 'error'].shape[0]
		self.total_error_norm = np.linalg.norm(self.results[self.start_index:, 'error']) / size
		self.interface_error_norm = np.linalg.norm(self.results[self.start_index:, 'abs_error_h']) / size

		self.steps = pl.DataFrame()
		if analysis_type == AnalysisType.standard:
			for i, t in enumerate(self.times):
				df = pl.read_csv(f'{self.inp_dir}/{stepsf}{i}.{ext}', separator=' ').with_columns(
					time=self.times[i],
					step=i
				).with_row_index('node')
				self.steps.vstack(df, in_place=True)
			self.steps.rechunk()
			
			self.errors = self.steps['time', 'x', 'y', 'error']
			self.solns = self.steps['time', 'x', 'y', 'u']
			self.exact = self.steps['time', 'x', 'y', 'exact_u']

			self.interface_node = self.step(0).filter(pl.col('x').is_close(self.results[0, 'h']))[0, 'node']

			if not self.steps.filter(pl.col('node').eq(self.interface_node) & ~pl.col('u').is_close(0.0)).is_empty():
				print("\n\nWARNING: non-zero temperature detected at interface, check results!")
				print(f"system: {self.params}\n\n")

			self.create_outdir()	

	def step(self, n: int) -> pl.DataFrame:
		return self.steps.filter(step=n)

	def errors_at_node(self, node: int | None = None) -> pl.DataFrame:
		if node is None:
			node =  self.step(self.start_index)['error'].arg_max()
		
		return self.steps.filter(
			pl.col('node') == node
		).with_columns(
			abs_error=pl.col('error').abs(),
			rel_error=(pl.col('error') / pl.col('exact_u')).abs()
		)
	
	def reshape_data(self, step: int, key: str) -> pl.DataFrame:
		df = self.step(step)
		return df.pivot(on='y', index='x', values=key)
	
	@property
	def max_dA(self):
		return self.results['dA'].max()
	
	def create_outdir(self) -> None:
		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)

	@property
	def fitted_De(self) -> float | None:
		if self._De is None:
			def f(x, a, b):
				return np.sqrt(x*a) + b

			try:
				optimal, _ = curve_fit(f, self.results['time'], self.results['h'], p0=[1.0, 0.0])
				self._De = optimal[0]
			except Exception as e:
				print(f"Failed to fit De for {self.params.title}, error: {e}")
		
		return self._De
	
	def plot_error_norm(self, ax) -> None:
		ax.semilogy(self.results['time'], self.results['error'])
		ax.set_xlabel('Time')
		ax.set_ylabel('$L_2$ norm of errors')

	def plot_error_at_node(self, ax, data=None, node=None) -> pl.DataFrame | None:
		if data is None:
			df = self.errors_at_node(node)
		else:
			df = data

		ax.semilogy(df['time'], df['abs_error'], label='Absolute')
		ax.semilogy(df['time'], df['rel_error'], label='Relative')

		ax.legend()
		ax.set_xlabel('Time')
		ax.set_ylabel('Magnitude of errors')

		if data is None:
			return df

	def plot_interface(self, ax, de=True, label=None) -> None:
		if label is None:
			label = 'Interface location'
		ax.plot(self.results['time'], self.results['h'], label=label)

		if de: 
			if self.fitted_De is not None:
				x = np.linspace(self.results[0, 'time'], self.results[-1, 'time'])
				ax.plot(x, np.sqrt(x * self.fitted_De), ls='--', label=fr'Fitted $D_e={{{self.fitted_De:.4e}}}$')
			
			ax.plot(self.results['time'], self.results['exact_h'], ls=':', label=fr'Analytical $D_e={{{self.system.De:4e}}}$')

		ax.legend()
		ax.set_xlabel('time')
		ax.set_ylabel(r'Interface location $h(t)$')
	
	def plot_de(self, ax) -> None:
		ax.plot(np.square(self.results['h']) / self.results['time'], label='Instantaneous $D_e$')

		ax.axhline(self.system.De, c='C1', ls='--', label='Analytical')
		ax.axhline(self.fitted_De, c='C2', ls='--', label='Numerical')


		ax.set_xlabel('time step')
		ax.set_ylabel('$D_e$')
		ax.legend()
	
	def plot_interface_error(self, abs_ax) -> None:
		rel_ax = abs_ax.twinx()

		abs_plot, = abs_ax.plot(self.results['time'], self.results['abs_error_h'], label='Absolute Error')
		rel_plot, = rel_ax.plot(self.results['time'], self.results['rel_error_h'], label='Relative Error', c='C1')

		abs_ax.legend(handles=[abs_plot, rel_plot])

		abs_ax.set_xlabel('Time')
		abs_ax.set_ylabel('Absolute Error')
		rel_ax.set_ylabel('Relative Error')
	
	def plot_mesh(self, ax, step: int, cmap: str = 'jet') -> None:
		frame = self.steps.filter(step=step)

		plot = ax.scatter(frame['x'], frame['y'], c=frame['u'], cmap=cmap)

		ax.axvline(self.results[step, 'h'], c='black', ls='--', label='Interface')
		ax.set_title(f't={self.results[step, "time"]:.4f}')

		return plot

	def plot_profile(self, ax, step: int, y_idx: int = 0, markerstep: int=10) -> None:
		frame = self.reshape_data(step, 'u')
		exact_frame = self.reshape_data(step, 'exact_u')

		ax.plot(frame['x'], frame[:, y_idx+1], label='Numerical')
		ax.scatter(exact_frame[::markerstep, 'x'], exact_frame[::markerstep, y_idx+1], c='C1', label='Analytical')

		ax.axvline(self.results[step, 'h'], c='C2', ls='--', label='Numerical Interface')
		ax.axvline(self.results[step, 'exact_h'], c='C3', ls='--', label='Analytical Interface')

		ax.set_title(f't={self.results[step, "time"]:.4f}')

	def anim_profile(self, fig, prof, error, nodes: int | None = 15, zoom: bool = False):
		def _update(step):
			
			for time in times:
				time.set_text(f't={self.results[step, "time"]:.4f}')

			if prof is not None:
				frame_num = self.reshape_data(step, 'u')
				frame_exact = self.reshape_data(step, 'exact_u')

				prof_lines['num'].set_data(frame_num['x'], frame_num[:, 1])
				prof_lines['ana'].set_data(frame_exact['x'], frame_exact[:, 1])

				prof_lines['iface_num'].set_xdata([self.results[step, 'h']])
				prof_lines['iface_ana'].set_xdata([self.results[step, 'exact_h']])

				if nodes is not None:
					prof_lines['nodes'].set_offsets(frame_num[node_range, 'x':2])
				
				if zoom:
					zoom_lines['num'].set_data(frame_num['x'], frame_num[:, 1])
					zoom_lines['ana'].set_data(frame_exact['x'], frame_exact[:, 1])

					h = self.results[step, 'h']
					exact_h = self.results[step, 'exact_h']

					zoom_lines['iface_num'].set_xdata([h])
					zoom_lines['iface_ana'].set_xdata([exact_h])

					x_halfrange = np.abs(h - exact_h)
					zoom_ax.set_xlim(h-2*x_halfrange, h+2*x_halfrange)

				
			if error is not None:
				frame = self.reshape_data(step, 'error')

				error_lines['error'].set_data(frame['x'], frame[:, 1])
				error.set_ylim(frame[:, 1].min() * 1.1, frame[:, 1].max() * 1.1)

				error_lines['iface_num'].set_xdata([self.results[step, 'h']])
				error_lines['iface_ana'].set_xdata([self.results[step, 'exact_h']])

				if nodes is not None:
					error_lines['nodes'].set_offsets(frame[node_range, 'x':2])
	
		times = []
		node_range = list(range((self.params.nx1*(self.params.x-1))-nodes, (self.params.nx1*(self.params.x-1))+nodes))
		
		if prof is not None:
			prof_lines = {}
			times.append(prof.annotate(
				f't={self.results[0, "time"]:.4f}',
				xy=(0.6, 0.9),
				xycoords='axes fraction'
			))

			frame_num = self.reshape_data(0, 'u')
			frame_exact = self.reshape_data(0, 'exact_u')

			prof_lines['num'], = prof.plot(frame_num['x'], frame_num[:, 1], label='Numerical')
			prof_lines['ana'], = prof.plot(frame_exact['x'], frame_exact[:, 1], ls='--', label='Analytical')

			prof_lines['iface_num'] = prof.axvline(self.results[0, 'h'], ls='--', label='Numerical Interface', c='C2')
			prof_lines['iface_ana'] = prof.axvline(self.results[0, 'exact_h'], ls='--', label='Analytical Interface', c='C3')

			if nodes is not None:
				prof_lines['nodes'] = prof.scatter(frame_num[node_range, 'x'], frame_num[node_range, 1], label='nodes')
			
			if zoom:
				zoom_lines = {}
				zoom_ax = prof.inset_axes(
					[0.4, 0.1, 0.2, 0.7],
					xlim=(0.0, 0.2),
					ylim=(-5e-5, 5e-5)
				)
				prof.indicate_inset_zoom(zoom_ax)
				fig.add_axes(zoom_ax)

				zoom_lines['num'], = zoom_ax.plot(frame_num['x'], frame_num[:, 1], label='Numerical')
				zoom_lines['ana'], = zoom_ax.plot(frame_exact['x'], frame_exact[:, 1], ls='--', label='Analytical')

				zoom_lines['iface_num'] = zoom_ax.axvline(self.results[0, 'h'], ls='--', label='Numerical Interface', c='C2')
				zoom_lines['iface_ana'] = zoom_ax.axvline(self.results[0, 'exact_h'], ls='--', label='Analytical Interface', c='C3')


			prof.legend()

			prof.set_ylabel('Temperature')
			prof.set_xlabel('x')
		
		if error is not None:
			error_lines = {}
			times.append(error.annotate(
				f't={self.results[0, "time"]:.4f}',
				xy=(0.6, 0.9),
				xycoords='axes fraction'
			))

			frame = self.reshape_data(0, 'error')
			error_lines['error'], = error.plot(frame['x'], frame[:, 1], label='Error')
			error_lines['iface_num'] = error.axvline(self.results[0, 'h'], ls='--', label='Numerical Interface', c='C2')
			error_lines['iface_ana'] = error.axvline(self.results[0, 'exact_h'], ls='--', label='Analytical Interface', c='C3')

			if nodes is not None:
				error_lines['nodes'] = error.scatter(frame[node_range, 'x'], frame[node_range, 1], label='nodes')

			error.legend()

			error.set_ylabel('Error')
			error.set_xlabel('x')

		return FuncAnimation(
			fig, 
			_update,
			self.results.shape[0]
		)