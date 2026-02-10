#!/usr/bin/env python3

import os
from enum import StrEnum
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from visualisation.run import Run
from systems.systems import AnalysisType

class PltStyle(StrEnum):
	report = 'styles/report.mplstyle'
	standard = 'styles/standard.mplstyle'

class Visualiser:
	def __init__(
		self,
		prefix: str,
		analysis_type: AnalysisType,
		resd: str = 'RESLT',
		outd: str = 'imgs',
		interactive: bool = False,
		style: PltStyle = PltStyle.standard,
		reference: Run = None,
		*args,
		**kwargs
	):
		self.prefix = prefix
		self.analysis_type = analysis_type
		self.outd = outd
		self.interactive = interactive
		self.style = style
		self.reference = reference
		self.report = self.style == PltStyle.report

		plt.style.use(f'{os.path.dirname(__file__)}/{self.style}')

		if not os.path.exists(self.outd):
			os.mkdir(self.outd)

		self.runs = []
		for rundir in os.listdir(resd):
			rund = f'{resd}/{rundir}'
			if not os.path.isdir(rund):
				continue

			print(f'Processing results in {rund}')

			self.runs.append(Run(
				rund,
				f'{self.outd}/{rundir}',
				analysis_type,
				reference = self.reference
			))
		
		if not self.interactive:
			self.default_run(*args, **kwargs)
	
	def default_run(self, *args, **kwargs) -> None:
		self.aggregate_data = pl.DataFrame()
		for run in self.runs:
			df = pl.DataFrame(data={
				'x': [run.params.x],
				't': [run.params.t],
				'dx': [run.dx],
				'dt': [run.params.dt],
				'error': [run.total_error_norm],
				'interface_error': [run.interface_error_norm]
			})
			self.aggregate_data.vstack(df, in_place=True)
		
		self.aggregate_data.rechunk()

		match self.analysis_type:
			case AnalysisType.standard:	self.standard(*args, **kwargs)
			case AnalysisType.dx | AnalysisType.dt:	self.plot_analysis(*args, **kwargs)
			case AnalysisType.sensitivity: 	self.sensitivity_analysis(*args, **kwargs)
		self.plot_interfaces()
	
	def _make(self, run) -> None:
		try:
			self.plot_errors(run)
			self.plot_interface(run)
			self.plot_de(run)
			self.plot_interface_temp(run)
			self.anim_surf_prof(run)

			# if self.report:
				# self.subplot_surf_prof(run)
			
			if run.params.nx == 20:
				self.subplot_mesh(run)
		except:
			raise Exception("SMTH WENT WRONG")
	
	def standard(self, nworkers=10, *args, **kwargs) -> None:
		for run in self.runs:
			print(f'Standard processing for: \n{run.params}')
			self._make(run)
		# with ProcessPoolExecutor(max_workers=nworkers) as executor:
		# 	for run, threads in zip(self.runs, executor.map(self._make, self.runs)):
		# 		print(f'Standard processing for: \n{run.params}')
	
	def trial(self):
		for run in self.runs:
			self.plot_errors(run)

	def savefig(self, fig, outdir, title):
		if self.report:
			title = f'{outdir}/{self.prefix}-report-{title}.png'
		else:
			title = f'{outdir}/{self.prefix}-{title}.png'

		fig.savefig(title)

	def plot_errors(self, run) -> None:
		if self.report:
			axs = {}
			fig1, axs['norm'] = plt.subplots(1)
			fig2, axs['max'] = plt.subplots(1)
		else:
			fig, axs = plt.subplot_mosaic(
				[['norm', 'max']],
				figsize=(10, 5)
			)
		
		run.plot_error_norm(axs['norm'])
		df = run.plot_error_at_node(axs['max'], node=None)

		if self.report:
			self.savefig(fig1, run.out_dir, 'error-norms')
			self.savefig(fig2, run.out_dir, 'error-max')
		else:
			axs['max'].set_title(f'Magnitude of errors at x={df[0, "x"]}')
			fig.suptitle(fr'Errors ({run.params.title})')

			self.savefig(fig, run.out_dir, 'errors')	

	def plot_interface(self, run) -> None:
		if self.report:
			axs = {}
			fig, axs['iface'] = plt.subplots(1)
		else:
			fig, axs = plt.subplot_mosaic(
				[
					['iface'],
					['error']
				],
				figsize=(7, 10)
			)

		run.plot_interface(axs['iface'])

		if not self.report:
			run.plot_interface_error(axs['error'])

			axs['iface'].set_title('Interface position with time')
			axs['error'].set_title('Error in interface position')
		
		self.savefig(fig, run.out_dir, 'interface')



	def plot_de(self, run) -> None:
		pass
		# raise NotImplementedError()

	def plot_interface_temp(self, run) -> None:
		pass
		# raise NotImplementedError()

	def anim_surf_prof(self, run) -> None:
		pass
		# raise NotImplementedError()


	def plot_interfaces(self) -> None:
		pass

	def plot_analysis(self) -> None:
		pass

	def sensitivity_analysis(self) -> None:
		pass