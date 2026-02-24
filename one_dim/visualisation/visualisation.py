import os
from enum import StrEnum, auto
# from concurrent.futures import ProcessPoolExecutor

import numpy as np
import polars as pl
import matplotlib.pyplot as plt

from visualisation.run import Run
from systems.systems import AnalysisType, System

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

		system = System.from_file(f'{results_dir}/system')

		self.runs = []
		for rundir in os.listdir(results_dir):
			rund = f'{results_dir}/{rundir}'
			if not os.path.isdir(rund):
				continue

			print(f'Processing results in {rund}')

			self.runs.append(Run(
				rund,
				f'{self.out_dir}/{rundir}',
				system,
				analysis_type == AnalysisType.standard,
				reference = self.reference
			))
	
		self.aggregate_data = pl.DataFrame()
		for run in self.runs:
			df = pl.DataFrame(data={
				'x': [run.params.x],
				't': [run.params.t],
				'dA': [run.max_elem_size],
				'dt': [run.params.dt],
				'bulk_error': [run.total_error_norm],
				'iface_error': [run.interface_error_norm]
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
		self.plot_errors(run)
		self.plot_interface(run)
		self.plot_de(run)
		self.anim_profile(run, not self.report)
		
		self.plot_mesh(run)
		self.plot_profiles(run)

		if run.params.nx == 20:
			self.plot_mesh(run)
	
	def standard(self, nworkers=10, *args, **kwargs) -> None:
		for run in self.runs:
			print(f'Standard processing for: \n{run.params}')
			self._make(run)
		# with ProcessPoolExecutor(max_workers=nworkers) as executor:
		# 	for run, threads in zip(self.runs, executor.map(self._make, self.runs)):
		# 		print(f'Standard processing for: \n{run.params}')

	def savefig(self, fig, outdir, title, ext='png'):
		if self.report:
			title = f'{outdir}/{self.prefix}-report-{title}.{ext}'
		else:
			title = f'{outdir}/{self.prefix}-{title}.{ext}'

		fig.savefig(title)
		plt.close()

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
		fig, ax = plt.subplots(1)

		run.plot_de(ax)

		if not self.report:
			ax.set_title('$D_e$ analytical vs numerical')

		self.savefig(fig, run.out_dir, 'de')

	def anim_profile(self, run, zoom: bool = False) -> None:
		if self.report:
			fig, ax = plt.subplots(1)
			anim = run.anim_profile(fig, prof=ax, zoom=zoom)
			title = f'{run.out_dir}/{self.prefix}-report-profile.mp4'
		else:
			fig, axs = plt.subplot_mosaic(
				[
					['prof'],
					['error']
				],
				sharex=True,
				figsize=(8, 10),
			)
				
			anim = run.anim_profile(fig, prof=axs['prof'], error=axs['error'], zoom=zoom)

			axs['prof'].set_xlabel('')
			axs['prof'].set_title('Temperature profile')
			axs['error'].set_title('Error in profile')
			title = f'{run.out_dir}/{self.prefix}-profile.mp4'

		anim.save(title)
		plt.close()

	def plot_analysis(self, ref: int | str = 0.09) -> None:
		def _label(col):
			x, t = col.strip('{}').split(',')

			match int(x):
				case 2: return f'Linear BDF{t}'
				case 3: return f'Quadratic BDF{t}'
				case 4: return f'Cubic BDF{t}'
				case _: raise ValueError("Invalid element order")

		match self.analysis_type:
			case AnalysisType.dx:
				xlabel = 'dA'
				fname = 'dx_errors'
				ref_variable = 'x'
			case AnalysisType.dt:
				xlabel = 'dt'
				fname = 'dt_errors'
				ref_variable = 't'

		if self.report:
			axs = {}
			fig1, axs['bulk'] = plt.subplots(
				1, 
				figsize=(9, 6),
				subplot_kw={
					'xscale': 'log',
					'yscale': 'log'
				}
			)
			fig2, axs['iface'] = plt.subplots(
				1, 
				figsize=(9, 6),
				subplot_kw={
					'xscale': 'log',
					'yscale': 'log'
				}
			)
		else:
			fig, axs = plt.subplot_mosaic(
				[
					['bulk'],
					['iface']
				],
				sharex=True,
				figsize=(8, 10),
				subplot_kw={
					'xlabel': xlabel,
					'ylabel': 'Normalized $L_2$ norm of error',
					'xscale': 'log',
					'yscale': 'log'
				},
				per_subplot_kw={
					'bulk': {
						'title': 'Bulk error norms'
					},
					'iface': {
						'title': 'Interface error norms'
					}
				}
			)

		iface_df = self.aggregate_data.pivot(index=xlabel, on=['x', 't'], values='iface_error')
		bulk_df = self.aggregate_data.pivot(index=xlabel, on=['x', 't'], values='bulk_error')

		for i in range(1, bulk_df.shape[1]):
			axs['bulk'].scatter(bulk_df[:, 0], bulk_df[:, i], label=_label(bulk_df.columns[i]))
			axs['iface'].scatter(iface_df[:, 0], iface_df[:, i], label=_label(iface_df.columns[i]))

		xs = np.linspace(self.aggregate_data[xlabel].min(), self.aggregate_data[xlabel].max())
		for i in self.aggregate_data[ref_variable].unique().sort():
			if type(ref) == float:
				ref_point = self.aggregate_data.filter(pl.col(ref_variable) == i).sort((pl.col(xlabel) - ref).abs())[0]
			elif ref == 'max':
				ref_point = self.aggregate_data.filter(pl.col(ref_variable) == i).sort(xlabel, descending=True)[0]
			elif ref == 'min':
				ref_point = self.aggregate_data.filter(pl.col(ref_variable) == i).sort(xlabel, descending=False)[0]
			else:
				raise ValueError('Invalid reference point argument')
			
			axs['bulk'].plot(
				xs, 
				(xs/ref_point[0, xlabel]) ** i * ref_point[0, 'bulk_error'], 
				label=f'$\mathcal{{O}}({ref_variable}^{i})$', 
				ls='--', 
				alpha=0.8
			)
			axs['iface'].plot(
				xs, 
				(xs/ref_point[0, xlabel]) ** i * ref_point[0, 'iface_error'], 
				label=f'$\mathcal{{O}}({ref_variable}^{i})$', 
				ls='--', 
				alpha=0.8
			)

		if self.report:
			for ax in axs.values():
				box = ax.get_position()
				ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
				plt.tight_layout()
			self.savefig(fig1, self.out_dir, f'{fname}_bulk')
			self.savefig(fig2, self.out_dir, f'{fname}_iface')
		else:
			axs['bulk'].legend()
			self.savefig(fig, self.out_dir, fname)

	def sensitivity_analysis(self) -> None:
		raise NotImplementedError()