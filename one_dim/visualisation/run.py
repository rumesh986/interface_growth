#!/usr/bin/env python3

import os
from typing import Self

import numpy as np
import polars as pl

from systems.systems import Params, AnalysisType

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
		self.params = Params.from_file('config')
		self.reference = reference

		print(os.getcwd())

		self.results = pl.read_csv(f'{self.inp_dir}/{resultsf}.{ext}', separator=' ').with_columns(
				interface_errors=(pl.col('exact_h') - pl.col('h')).abs()
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
			self.interface_error_norm = np.linalg.norm(self.results[self.start_index:, 'interface_errors']) / size

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
