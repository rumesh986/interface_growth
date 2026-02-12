#!/usr/bin/env python3

import os
from enum import StrEnum, auto
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from visualisation.run import Run
from systems.systems import AnalysisType

class PltStyle(StrEnum):
	report = auto()
	standard = auto()


class Visualiser:
	def __init__(
		self,
		prefix: str,
		analysis_type: AnalysisType,
		results_dir: str = 'RESLT',
		out_dir: str = 'imgs',
		style: PltStyle = PltStyle.standard,
		reference: Run = None,
		*args,
		**kwargs
	):
		self.prefix = prefix
		self.analysis_type = analysis_type
		self.out_dir = out_dir
		self.style = style
		self.reference = reference
		self.report = self.style == PltStyle.report

		plt.style.use(f"visualisation.styles.{self.style}")

		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)

		self.runs = []
		for rundir in os.listdir(results_dir):
			rund = f'{results_dir}/{rundir}'
			if not os.path.isdir(rund):
				continue

			print(f'Processing results in {rund}')

			self.runs.append(Run(
				rund,
				f'{self.out_dir}/{rundir}',
				analysis_type,
				reference = self.reference
			))
	
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
		
	def default_run(self, *args, **kwargs) -> None:
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

	def savefig(self, fig, outdir, title, ext='png'):
		if self.report:
			title = f'{outdir}/{self.prefix}-report-{title}.{ext}'
		else:
			title = f'{outdir}/{self.prefix}-{title}.{ext}'

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

	def plot_interfaces(self) -> None:
		fig, ax = plt.subplots(1,1)

		self.runs[0].plot_interface(ax, label=self.runs[0].params.title)
		for run in self.runs[1:]:
			run.plot_interface(ax, False, label=run.params.title)

		box = ax.get_position()
		ax.set_position([box.x0, box.y0 + box.height * 0.005, box.width, box.height * 0.995])
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, prop={'size': 8})

		self.savefig(fig, self.out_dir, 'interfaces')
	
	def plot_mesh(self, run, num_ax=4):
		nrows = int(np.ceil(np.sqrt(num_ax)))
		ncols = int(np.ceil(num_ax / nrows))

		fig, axs = plt.subplots(
			nrows, 
			ncols,
			sharex=True,
			sharey=True,
			figsize=(14, 10),
			layout="compressed"
		)

		frame_step = int(np.floor(run.results.shape[0] / num_ax))
		for i, ax in enumerate(axs.flat):
			if (i*frame_step > run.results.shape[0]):
				run.plot_mesh(ax, run.results.shape[0]-1)
				break
			cbar_plot = run.plot_mesh(ax, i * frame_step)
		
		fig.colorbar(cbar_plot, ax=axs.flatten())
		
		if not self.report:
			fig.suptitle('Nodal positions')
		fig.supxlabel('X')
		fig.supylabel('Y')

		self.savefig(fig, run.out_dir, 'mesh')

	def plot_profiles(self, run, num_ax=4) -> None:
		nrows = int(np.ceil(np.sqrt(num_ax)))
		ncols = int(np.ceil(num_ax / nrows))

		fig, axs = plt.subplots(
			nrows,
			ncols,
			sharex=True,
			sharey=True,
			figsize=(14, 10),
			layout="compressed"
		)

		frame_step = int(np.floor(run.results.shape[0]) / num_ax)

		for i, ax in enumerate(axs.flat):
			frame = i * frame_step
			if frame > run.results.shape[0]:
				run.plot_profile(ax, run.results.shape[0] - 1)
				break

			run.plot_profile(ax, frame)

		axs.flat[0].legend()

		if not self.report:
			fig.suptitle('Temperature profiles')
		fig.supxlabel('X')
		fig.supylabel('Temperature')

		self.savefig(fig, run.out_dir, 'profiles')

	def plot_de(self, run) -> None:
		pass
		# raise NotImplementedError()

	def anim_surf_prof(self, run) -> None:
		pass
		# raise NotImplementedError()

	def plot_analysis(self) -> None:
		pass

	def sensitivity_analysis(self) -> None:
		pass