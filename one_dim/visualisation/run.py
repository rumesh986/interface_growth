#!/usr/bin/env python3

import os
from typing import Self

import numpy as np
import polars as pl
from scipy.optimize import curve_fit

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
		# self.results = pd.read_csv(f'{self.inp_dir}/{resultsf}.{self.ext}', sep=' ')
		self.times = self.results['time'] #.values
		# self.interface = self.results['h']
		# self.error_norms = self.results['error']

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
			self.start_index = self.params.t+1

			size = self.results[self.start_index:, 'error'].shape[0]
			self.total_error_norm = np.linalg.norm(self.results[self.start_index:, 'error']) / size
			self.interface_error_norm = np.linalg.norm(self.results[self.start_index:, 'abs_error_h']) / size

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
	def dx(self):
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
	
	def plot_interface_error(self, abs_ax) -> None:
		rel_ax = abs_ax.twinx()

		abs_plot, = abs_ax.plot(self.results['time'], self.results['abs_error_h'], label='Absolute Error')
		rel_plot, = rel_ax.plot(self.results['time'], self.results['rel_error_h'], label='Relative Error', c='C1')

		abs_ax.legend(handles=[abs_plot, rel_plot])

		abs_ax.set_xlabel('Time')
		abs_ax.set_ylabel('Absolute Error')
		rel_ax.set_ylabel('Relative Error')
	
	def plot_mesh(self, ax, step: int, cmap: str = 'jet') :
		frame = self.steps.filter(step=step)

		plot = ax.scatter(frame['x'], frame['y'], c=frame['u'], cmap=cmap)

		ax.axvline(self.results[step, 'h'], c='black', ls='--', label='Interface')
		ax.set_title(f't={self.results[step, "time"]:.4f}')

		return plot

	def plot_profile(self, ax, step: int, y_idx: int = 0, markerstep: int=10):
		frame = self.reshape_data(step, 'u')
		exact_frame = self.reshape_data(step, 'exact_u')

		ax.plot(frame['x'], frame[:, y_idx+1], label='Numerical')
		# ax.plot(exact_frame[::markerstep, 'x'], exact_frame[::markerstep, y_idx+1], marker='x', ls=' ')
		ax.scatter(exact_frame[::markerstep, 'x'], exact_frame[::markerstep, y_idx+1], c='C1', label='Analytical')
		# ax.scatter(frame[::10, 'x'], frame[::10, y_idx+1], c='C1', zorder=2)
		ax.axvline(self.results[step, 'h'], c='C2', ls='--', label='Numerical Interface')
		ax.axvline(self.results[step, 'exact_h'], c='C3', ls='--', label='Analytical Interface')

		ax.set_title(f't={self.results[step, "time"]:.4f}')