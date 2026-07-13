from typing import override
from warnings import warn

import numpy as np
import polars as pl
from scipy.optimize import curve_fit
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

	@override
	def fitted_De(self, yi: str) -> float | None:
		def f(x, a, b):
			return np.sqrt(x * a) + b

		try:
			optimal, _ = curve_fit(f, self.interface[self.COLS.time], self.interface[yi], p0=[1.0, 0.0])
			de = optimal[0]
			return de
		except Exception as e:
			print(f"Failed to fit De for {self.params.title}, error: {e}")
		
	@override
	def plot_error_norm(self, ax: Axes, *args, **kwargs) -> None:
		pass

	@override
	def plot_error_at_node(self, ax: Axes, data: pl.DataFrame | None = None, node: int | None = None, *args, **kwargs) -> pl.DataFrame | None:
		return None

	@override
	def plot_interface(self, ax: Axes, de: bool = True, label: str = 'Interface location', *args, **kwargs) -> None:
		x = [float(i) for i in self.interface.columns[2:]]

		start = 0
		end = self.interface.shape[0]
		step = 10

		for t in range(start, end, step):
			ax.plot(self.interface[t, 2:].transpose().to_series(), x, label=f'{self.interface[t, self.COLS.time]}')
		
		ax.legend()

	@override
	def plot_interface_error(self, abs_ax: Axes, *args, **kwargs) -> None:
		pass

	@override
	def plot_de(self, ax: Axes) -> None:
		for yi in self.interface.columns[2:]:
			de = self.fitted_De(yi)
			ax.plot(self.interface[self.COLS.time], self.interface[yi], label=f'y={float(yi):.3f} de={de:e}')

			# print(f'y={yi} de={de}')
		
		ax.plot(self.results[self.COLS.time], self.results[self.COLS.exact_h], ls=':', label=fr'Analytical $D_e={{{self.system.De:.5f}}}$')
		
		# ax.legend()

	@override
	def plot_mesh(self, ax: Axes, step: int, cmap: str = 'jet', *args, **kwargs) -> None:
		frame = self.step(step)
		iface = self.interface.filter(step=step)

		plot = ax.scatter(frame[self.COLS.x], frame[self.COLS.y], c=frame[self.COLS.temp], cmap=cmap)
		ax.plot(iface.row(0)[2:], [float(i) for i in iface.columns[2:]], label='Interface', c='black', ls='--')

		ax.set_title(f't={self.results[step, "time"]:.4f}')

		return plot

	@override
	def plot_profile(self, ax: Axes, step: int, y_idx: int = 0, markerstep: int = 10, *args, **kwargs) -> None:
		frame = self.step(step)
		iface = self.interface.filter(step=step)

		xs = iface.row(0)[2:]
		ys = tuple([float(i) for i in iface.columns[2:]])

		ax.tricontourf(frame[self.COLS.x], frame[self.COLS.y], frame[self.COLS.temp])
		ax.plot(xs, ys, label='Interface', c='black')

		ax.set_title(f't={self.results[step, self.COLS.time]:.4f}')


	@override
	def anim_profile(self, fig: Figure, prof: Axes | None = None, error: Axes | None = None, nodes: int | None = 15, zoom: bool = False, *args, **kwargs) -> Animation:
		def _update(step):
			for time in times:
				time.set_text(f't={self.results[step, self.COLS.time]:.4f}')

			if prof is not None:
				frame = self.step(step)
				iface = self.interface.filter(step=step)

				prof_lines['contour'].remove()

				prof_lines['contour'] = prof.tricontourf(frame[self.COLS.x], frame[self.COLS.y], frame[self.COLS.temp], cmap='jet', levels=40)
				prof_lines['iface'].set_data(iface.row(0)[2:], [float(i) for i in iface.columns[2:]])

		times = []

		if prof is not None:
			prof_lines = {}
			times.append(prof.annotate(
				f't={self.results[0, self.COLS.time]:.4f}',
				xy=(0.6, 0.9),
				xycoords='axes fraction'
			))

			frame = self.step(0)
			iface = self.interface.filter(step=0)
			prof_lines['contour'] = prof.tricontourf(frame[self.COLS.x], frame[self.COLS.y], frame[self.COLS.temp], cmap='jet', levels=40)
			prof_lines['iface'], = prof.plot(iface.row(0)[2:], [float(i) for i in iface.columns[2:]], ls='--', c='black')

			prof.set_ylabel('y')
			prof.set_xlabel('x')
		
		fig.colorbar(prof_lines['contour'])

		if 'anim_length' in kwargs:
			warn("Animation length control has not been implemented yet")

		return FuncAnimation(
			fig, 
			_update,
			self.results.shape[0]
		)