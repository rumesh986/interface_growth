from warnings import warn

import numpy as np
import polars as pl
from scipy.optimize import curve_fit
from matplotlib.pyplot import Axes, Figure
from matplotlib.animation import Animation, FuncAnimation

from visualisation.run_base import RunBase

class Run(RunBase):
	@property
	def total_error_norm(self) -> float:
		if self._total_error_norm is None:
			self._total_error_norm = np.linalg.norm(self.results[self.params.t+1:, self.COLS.error]) / self.results[self.params.t+1:].shape[0]
		
		return self._total_error_norm

	@property
	def interface_error_norm(self) -> float:
		if self._interface_error_norm is None:
			self._interface_error_norm = np.linalg.norm(self.results[self.params.t+1:, self.COLS.abs_error_h]) / self.results[self.params.t+1:].shape[0]
		
		return self._interface_error_norm

	@property
	def max_elem_size(self) -> float:
		return self.results[self.COLS.elem_size].max()

	@property
	def fitted_De(self) -> float | None:
		if self._De is None:
			def f(x, a, b):
				return np.sqrt(x*a) + b

			try:
				optimal, _ = curve_fit(f, self.results[self.COLS.time], self.results[self.COLS.h], p0=[1.0, 0.0])
				self._De = optimal[0]
			except Exception as e:
				print(f"Failed to fit De for {self.params.title}, error: {e}")
		
		return self._De
	
	def plot_error_norm(self, ax: Axes, *args, **kwargs) -> None:
		ax.semilogy(self.results[self.COLS.time], self.results[self.COLS.error])
		ax.set_xlabel('Time')
		ax.set_ylabel('$L_2$ norm of errors')

	def plot_error_at_node(self, ax: Axes, data: pl.DataFrame | None = None, node: int | None = None, *args, **kwargs) -> pl.DataFrame | None:
		if data is None:
			df = self.errors_at_node(node)
		else:
			df = data

		ax.semilogy(df[self.COLS.time], df[self.COLS.abs_error], label='Absolute')
		ax.semilogy(df[self.COLS.time], df[self.COLS.rel_error], label='Relative')

		ax.legend()
		ax.set_xlabel('Time')
		ax.set_ylabel('Magnitude of errors')

		if data is None:
			return df

	def plot_interface(self, ax: Axes, de: bool = True, label: str = 'Interface location', *args, **kwargs) -> None:
		ax.plot(self.results[self.COLS.time], self.results[self.COLS.h], label=label)

		if de: 
			if self.fitted_De is not None:
				x = np.linspace(self.results[0, self.COLS.time], self.results[-1, self.COLS.time])
				ax.plot(x, np.sqrt(x * self.fitted_De), ls='--', label=fr'Fitted $D_e={{{self.fitted_De:.4e}}}$')
			
			ax.plot(self.results[self.COLS.time], self.results[self.COLS.exact_h], ls=':', label=fr'Analytical $D_e={{{self.system.De:4e}}}$')

		ax.legend()
		ax.set_xlabel('time')
		ax.set_ylabel(r'Interface location $h(t)$')
	
	def plot_de(self, ax: Axes) -> None:
		ax.plot(np.square(self.results[self.COLS.h]) / self.results[self.COLS.time], label='Instantaneous $D_e$')

		ax.axhline(self.system.De, c='C1', ls='--', label='Analytical')
		ax.axhline(self.fitted_De, c='C2', ls='--', label='Numerical')


		ax.set_xlabel('time step')
		ax.set_ylabel('$D_e$')
		ax.legend()
	
	def plot_interface_error(self, abs_ax: Axes, *args, **kwargs) -> None:
		rel_ax = abs_ax.twinx()

		abs_plot, = abs_ax.plot(self.results[self.COLS.time], self.results[self.COLS.abs_error_h], label='Absolute Error')
		rel_plot, = rel_ax.plot(self.results[self.COLS.time], self.results[self.COLS.rel_error_h], label='Relative Error', c='C1')

		abs_ax.legend(handles=[abs_plot, rel_plot])

		abs_ax.set_xlabel('Time')
		abs_ax.set_ylabel('Absolute Error')
		rel_ax.set_ylabel('Relative Error')
	
	def plot_mesh(self, ax: Axes, step: int, cmap: str = 'jet', *args, **kwargs) -> None:
		frame = self.steps.filter(step=step)

		plot = ax.scatter(frame[self.COLS.x], frame[self.COLS.y], c=frame[self.COLS.temp], cmap=cmap)

		ax.axvline(self.results[step, self.COLS.h], c='black', ls='--', label='Interface')
		ax.set_title(f't={self.results[step, "time"]:.4f}')

		return plot

	def plot_profile(self, ax: Axes, step: int, y_idx: int = 0, markerstep: int = 10, *args, **kwargs) -> None:
		frame = self.reshape_data(step, self.COLS.temp)
		exact_frame = self.reshape_data(step, self.COLS.exact_temp)

		ax.plot(frame[self.COLS.x], frame[:, y_idx+1], label='Numerical')
		ax.scatter(exact_frame[::markerstep, self.COLS.x], exact_frame[::markerstep, y_idx+1], c='C1', label='Analytical')

		ax.axvline(self.results[step, self.COLS.h], c='C2', ls='--', label='Numerical Interface')
		ax.axvline(self.results[step, self.COLS.exact_h], c='C3', ls='--', label='Analytical Interface')

		ax.set_title(f't={self.results[step, self.COLS.time]:.4f}')

	def anim_profile(
		self, 
		fig: Figure, 
		prof: Axes | None = None, 
		error: Axes | None = None, 
		nodes: int = 15, 
		zoom: bool = False,
		*args, 
		**kwargs
	) -> Animation:

		def _update(step):
			for time in times:
				time.set_text(f't={self.results[step, self.COLS.time]:.4f}')

			if prof is not None:
				frame_num = self.reshape_data(step, self.COLS.temp)
				frame_exact = self.reshape_data(step, self.COLS.exact_temp)

				prof_lines['num'].set_data(frame_num[self.COLS.x], frame_num[:, 1])
				prof_lines['ana'].set_data(frame_exact[self.COLS.x], frame_exact[:, 1])

				prof_lines['iface_num'].set_xdata([self.results[step, self.COLS.h]])
				prof_lines['iface_ana'].set_xdata([self.results[step, self.COLS.exact_h]])

				if nodes > 0:
					prof_lines['nodes'].set_offsets(frame_num[node_range, self.COLS.x:2])
				
				if zoom:
					zoom_lines['num'].set_data(frame_num[self.COLS.x], frame_num[:, 1])
					zoom_lines['ana'].set_data(frame_exact[self.COLS.x], frame_exact[:, 1])

					h = self.results[step, self.COLS.h]
					exact_h = self.results[step, self.COLS.exact_h]

					zoom_lines['iface_num'].set_xdata([h])
					zoom_lines['iface_ana'].set_xdata([exact_h])

					x_halfrange = np.abs(h - exact_h)
					zoom_ax.set_xlim(h-2*x_halfrange, h+2*x_halfrange)

				
			if error is not None:
				frame = self.reshape_data(step, self.COLS.error)

				error_lines['error'].set_data(frame[self.COLS.x], frame[:, 1])
				error.set_ylim(frame[:, 1].min() * 1.1, frame[:, 1].max() * 1.1)

				error_lines['iface_num'].set_xdata([self.results[step, self.COLS.h]])
				error_lines['iface_ana'].set_xdata([self.results[step, self.COLS.exact_h]])

				if nodes > 0:
					error_lines['nodes'].set_offsets(frame[node_range, self.COLS.x:2])
	
		times = []
		nodes = min([self.params.nx1, self.params.nx2, nodes])
		node_range = list(range((self.params.nx1*(self.params.x-1))-nodes, (self.params.nx1*(self.params.x-1))+nodes))
		
		if prof is not None:
			prof_lines = {}
			times.append(prof.annotate(
				f't={self.results[0, "time"]:.4f}',
				xy=(0.6, 0.9),
				xycoords='axes fraction'
			))

			frame_num = self.reshape_data(0, self.COLS.temp)
			frame_exact = self.reshape_data(0, self.COLS.exact_temp)

			prof_lines['num'], = prof.plot(frame_num[self.COLS.x], frame_num[:, 1], label='Numerical')
			prof_lines['ana'], = prof.plot(frame_exact[self.COLS.x], frame_exact[:, 1], ls='--', label='Analytical')

			prof_lines['iface_num'] = prof.axvline(self.results[0, self.COLS.h], ls='--', label='Numerical Interface', c='C2')
			prof_lines['iface_ana'] = prof.axvline(self.results[0, self.COLS.exact_h], ls='--', label='Analytical Interface', c='C3')

			if nodes > 0:
				prof_lines['nodes'] = prof.scatter(frame_num[node_range, self.COLS.x], frame_num[node_range, 1], label='nodes')
			
			if zoom:
				zoom_lines = {}
				zoom_ax = prof.inset_axes(
					[0.4, 0.1, 0.2, 0.7],
					xlim=(0.0, 0.2),
					ylim=(-5e-5, 5e-5)
				)
				prof.indicate_inset_zoom(zoom_ax)
				fig.add_axes(zoom_ax)

				zoom_lines['num'], = zoom_ax.plot(frame_num[self.COLS.x], frame_num[:, 1], label='Numerical')
				zoom_lines['ana'], = zoom_ax.plot(frame_exact[self.COLS.x], frame_exact[:, 1], ls='--', label='Analytical')

				zoom_lines['iface_num'] = zoom_ax.axvline(self.results[0, self.COLS.h], ls='--', label='Numerical Interface', c='C2')
				zoom_lines['iface_ana'] = zoom_ax.axvline(self.results[0, self.COLS.exact_h], ls='--', label='Analytical Interface', c='C3')


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

			frame = self.reshape_data(0, self.COLS.error)
			error_lines['error'], = error.plot(frame[self.COLS.x], frame[:, 1], label='Error')
			error_lines['iface_num'] = error.axvline(self.results[0, self.COLS.h], ls='--', label='Numerical Interface', c='C2')
			error_lines['iface_ana'] = error.axvline(self.results[0, self.COLS.exact_h], ls='--', label='Analytical Interface', c='C3')

			if nodes > 0:
				error_lines['nodes'] = error.scatter(frame[node_range, self.COLS.x], frame[node_range, 1], label='nodes')

			error.legend()

			error.set_ylabel('Error')
			error.set_xlabel('x')

		if 'anim_length' in kwargs:
			warn("Animation length control has not been implemented yet")

		return FuncAnimation(
			fig, 
			_update,
			self.results.shape[0]
		)