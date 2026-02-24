import os
from typing import Self
from enum import StrEnum

import polars as pl
from matplotlib.pyplot import Axes, Figure
from matplotlib.animation import Animation

from systems.systems import Params, System

class OneDCols(StrEnum):
	time = 'time'

	h = 'h'
	exact_h = 'exact_h'
	
	x = 'x'
	y = 'y'
	
	temp = 'u'
	exact_temp = 'exact_u'
	
	error = 'error'
	abs_error_h = 'abs_error_h'
	rel_error_h = 'rel_error_h'
	abs_error = 'abs_error'
	rel_error = 'rel_error'
	
	elem_size = 'dA'
	
	node = 'node'
	step = 'step'

class RunBase:
	COLS: StrEnum = OneDCols
	
	_PL_SEPARATOR = ' '
	_RESULTS_EXPRS = [
		((pl.col(COLS.exact_h) - pl.col(COLS.h)).abs()).alias(COLS.abs_error_h),
		(((pl.col(COLS.exact_h) - pl.col(COLS.h)) / pl.col(COLS.h)).abs()).alias(COLS.rel_error_h)
	]

	def __init__(
		self,
		inp_dir: str,
		out_dir: str,
		system: System,
		load_data: bool = True,
		reference: Self | None = None,
		stepsf: str = 'steps/step',
		resultsf: str = 'results',
		paramsf: str = 'params',
		systemf: str = 'system',
		ext: str = 'dat',
	) -> None:
		self._inp_dir = inp_dir
		self.out_dir = out_dir
		self.system = system
		self.params = Params.from_file(f'{inp_dir}/{paramsf}')
		# self.system = System.from_file(f'{inp_dir}/{systemf}')
		self.reference = reference
		self._resultsf = resultsf
		self._stepsf = stepsf
		self._ext = ext

		self._De: float = None
		self.results: pl.DataFrame = None
		self.steps: pl.DataFrame = None
		self.interface_node: int = -1

		self._total_error_norm: float | None = None
		self._interface_error_norm: float | None = None

		self.load_results()
		if load_data:
			self.load_steps()

			self.interface_node = self.step(0).filter(pl.col(self.COLS.x).is_close(self.results[0, self.COLS.h]))[0, self.COLS.node]
			
			if not self.steps.filter(pl.col(self.COLS.node).eq(self.interface_node) & ~pl.col(self.COLS.temp).is_close(0.0)).is_empty():
				print("\n\nWARNING: non-zero temperature detected at interface, check results!")
				print(f"system: {self.params}\n\n")
			
			self.create_outdir()
	
	def load_results(self) -> None:
		self.results = pl.read_csv(
			f'{self._inp_dir}/{self._resultsf}.{self._ext}', 
			separator=self._PL_SEPARATOR
		).with_row_index(self.COLS.step).with_columns(self._RESULTS_EXPRS)
	
	def load_steps(self) -> None:
		if self.steps is not None:
			return

		if self.results is None:
			self.load_results()

		self.steps = pl.DataFrame()
		for i, t in enumerate(self.results[self.COLS.time]):
			df = pl.read_csv(
				f'{self._inp_dir}/{self._stepsf}{i}.{self._ext}', 
				separator=self._PL_SEPARATOR
			).with_columns([
				pl.lit(self.results[i, self.COLS.time]).alias(self.COLS.time),
				pl.lit(i).alias(self.COLS.step)
			]).with_row_index(self.COLS.node)

			self.steps.vstack(df, in_place=True)
		
		self.steps.rechunk()

	def reshape_data(self, step: int, key: str) -> pl.DataFrame:
		return self.step(step).pivot(on=self.COLS.y, index=self.COLS.x, values=key)

	def step(self, n: int) -> pl.DataFrame:
		return self.steps.filter(pl.col(self.COLS.step).eq(n))
	
	def create_outdir(self) -> None:
		if not os.path.exists(self.out_dir):
			os.mkdir(self.out_dir)
	
	@property
	def total_error_norm(self) -> float:
		raise NotImplementedError()
	
	@property
	def interface_error_norm(self) -> float:
		raise NotImplementedError()
	
	@property
	def max_elem_size(self) -> float:
		raise NotImplementedError()
	
	@property
	def fitted_de(self) -> float | None:
		raise NotImplementedError()
	
	def errors_at_node(self, node: int | None = None) -> pl.DataFrame:
		if node is None:
			node = self.step(self.params.t+1)[self.COLS.error].arg_max()
		
		return self.steps.filter(
			pl.col(self.COLS.node).eq(node)
		).with_columns([
			pl.col(self.COLS.error).abs().alias(self.COLS.abs_error),
			(pl.col(self.COLS.error) / pl.col(self.COLS.exact_temp)).abs().alias(self.COLS.rel_error),
		])
	
	def plot_error_norm(self, ax: Axes, *args, **kwargs) -> None:
		raise NotImplementedError()

	def plot_error_at_node(self, ax: Axes, data: pl.DataFrame | None = None, node: int | None = None, *args, **kwargs) -> pl.DataFrame | None:
		raise NotImplementedError()

	def plot_interface(self, ax: Axes, de: bool = True, label: str = 'Interface location', *args, **kwargs) -> None:
		raise NotImplementedError()

	def plot_interface_error(self, abs_ax: Axes, *args, **kwargs) -> None:
		raise NotImplementedError()

	def plot_de(self, ax: Axes) -> None:
		raise NotImplementedError()

	def plot_mesh(self, ax: Axes, step: int, cmap: str = 'jet', *args, **kwargs) -> None:
		raise NotImplementedError()

	def plot_profile(self, ax: Axes, step: int, y_idx: int = 0, markerstep: int = 10, *args, **kwargs) -> None:
		raise NotImplementedError()

	def anim_profile(self, fig: Figure, prof: Axes | None = None, error: Axes | None = None, nodes: int | None = 15, zoom: bool = False, *args, **kwargs) -> Animation:
		raise NotImplementedError()