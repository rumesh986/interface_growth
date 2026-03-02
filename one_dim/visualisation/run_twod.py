from typing import override

import numpy as np
import polars as pl
from matplotlib.pyplot import Axes, Figure
from matplotlib.animation import Animation, FuncAnimation

from visualisation.run_base import RunBase

class Run_TwoD(RunBase):
	@override
	def load_additional_data(self) -> None:
		self.interface = pl.read_csv(
			f'{self._inp_dir}/{self._ifacef}.{self._ext}',
			separator=self._PL_SEPARATOR
		).with_row_index(self.COLS.step)
	
	@property
	@override
	def total_error_norm(self) -> float:
		pass

	@property
	@override
	def interface_error_norm(self) -> float:
		pass

	@property
	@override
	def fitted_De(self) -> float | None:
		pass

	@override
	def plot_error_norm(self, ax: Axes, *args, **kwargs) -> None:
		pass

	@override
	def plot_error_at_node(self, ax: Axes, data: pl.DataFrame | None = None, node: int | None = None, *args, **kwargs) -> pl.DataFrame | None:
		return None

	@override
	def plot_interface(self, ax: Axes, de: bool = True, label: str = 'Interface location', *args, **kwargs) -> None:
		pass

	@override
	def plot_interface_error(self, abs_ax: Axes, *args, **kwargs) -> None:
		pass

	@override
	def plot_de(self, ax: Axes) -> None:
		pass

	@override
	def plot_mesh(self, ax: Axes, step: int, cmap: str = 'jet', *args, **kwargs) -> None:
		pass

	@override
	def plot_profile(self, ax: Axes, step: int, y_idx: inst = 0, markerstep: int = 10, *args, **kwargs) -> None:
		frame = self.steps.filter(pl.col(self.COLS.step).eq(step))
		iface = self.interface.filter(pl.col(self.COLS.step).eq(step))

		xs = iface.row(0)[2:]
		ys = tuple([float(i) for i in iface.columns[2:]])

		ax.tricontourf(frame[self.COLS.x], frame[self.COLS.y], frame[self.COLS.temp])
		ax.plot(xs, ys, label='Interface', c='black')



	@override
	def anim_profile(self, fig: Figure, prof: Axes | None = None, error: Axes | None = None, nodes: int | None = 15, zoom: bool = False, *args, **kwargs) -> Animation:
		pass