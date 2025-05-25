import matplotlib.pyplot as plt
import matplotlib.animation as anim
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import numpy as np

class Visualizer:
	def __init__(self, inp_folder, out_folder, prefix, n, times=None, solnf='soln', exact_solnf='exact_soln', errorf='error', ext="dat"):
		self.inp_folder = inp_folder
		self.out_folder = out_folder
		self.n = n
		self.times = times
		self.solnf = solnf
		self.exact_solnf = exact_solnf
		self.errorf = errorf
		self.ext = ext
		self.prefix = prefix

		if times:
			self.times = []
			with open(f'{inp_folder}/{times}', 'r') as f:
				while line := f.readline():
					self.times.append(float(line))
			self.times = np.array(self.times)
		else:
			self.times = range(n)

	def read_file(self, prefix, index):
		file = f"{self.inp_folder}/{prefix}{index}.{self.ext}"

		with open(file, 'r') as f:
			lines = f.readlines()

		data_start = False
		data = []

		for line in lines:
			if line.startswith("ZONE"):
				data_start = True
				continue
			
			if data_start:
				data.append(list(map(float, line.split())))
		
		data = np.array(data)
		
		return data[:, 0], data[:, 1], data

	def plot_line(self, ax, prefix, n, options):
		x, y, _ = self.read_file(prefix, n)

		if options:
			return ax.plot(x, y, **options)
		else:
			return ax.plot(x, y)

	def make_anim(self, fname, opts, ax_options):
		print(f"Making animation of {opts.keys()}")

		plt.cla()

		fig, ax = plt.subplots(1,1, subplot_kw=ax_options)
		lines = {}
		time_text = ax.annotate(f't={self.times[0]}', xy=(0.8, 0.9), xycoords='axes fraction')

		for k, v in opts.items():
			line = self.plot_line(ax, k, 0, v)
			lines[k] = line[0]
		
		def update(n: int):
			time_text.set_text(f't={self.times[n]}')
			for k, _ in opts.items():
				x, y, _ = self.read_file(k, n)

				lines[k].set_data(x, y)

		plt.legend(loc='upper left')
		ani = anim.FuncAnimation(fig, update, self.n)
		ani.save(f'{self.out_folder}/{self.prefix}-{fname}')

	def make_diffs_anim(self, fname, prefixs, ax_options):
		plt.cla()

		fig, ax = plt.subplots(1,1, subplot_kw=ax_options)

		time_text = ax.annotate(f't={self.times[0]}', xy=(0.8, 0.9), xycoords='axes fraction')
		x, y1, _ = self.read_file(prefixs[0], 0)
		_, y2, _ = self.read_file(prefixs[1], 0)

		line = plt.plot(x, y1-y2)[0]

		def update(n):
			time_text.set_text(f't={self.times[n]}')
			x, y1, _ = self.read_file(prefixs[0], n)
			_, y2, _ = self.read_file(prefixs[1], n)
			line.set_data(x, y1-y2)

		ani = anim.FuncAnimation(fig, update, self.n)
		ani.save(f'{self.out_folder}/{self.prefix}-{fname}')

	
	def plot_errors(self):
		plt.cla()
		x, _, data = self.read_file('error', 1)
		error = data[:, 2]
		node = np.argmax(np.fabs(error))

		abs_errors = np.array([np.fabs(self.read_file('error', i)[2][node, 2]) for i in range(self.n)])
		exact_vals = np.array([np.fabs(self.read_file('exact_soln', i)[1][node]) for i in range(self.n)])
		rel_errors = abs_errors / exact_vals

		plt.semilogy(abs_errors, label='absolute')
		plt.semilogy(rel_errors, label='relative')

		plt.title(f'Errors at x={x[node]}')
		plt.xlabel('timestep')
		plt.ylabel('Error magnitude')
		plt.legend()

		plt.savefig(f'{self.out_folder}/{self.prefix}-errors.png')

	def plot_error_norms(self):
		plt.cla()
		
		abs_error = np.array([np.linalg.norm(self.read_file('error', i)[2][:,2]) for i in range(self.n)])
		# exact_vals = np.array([np.fabs(self.read_file('exact_soln', i)[1]) for i in range(self.n)])
		# rel_error = np.array([np.linalg.norm(self.read_file('error', i)[2][:,2]/self.read_file('exact_soln', i)[1]) for i in range(self.n)])

		plt.semilogy(abs_error, label='Absolute error')
		# plt.semilogy(rel_error, label='Relative error')

		plt.title('L2 norm of errors')
		plt.xlabel('timestep')
		plt.ylabel('Error noem')
		plt.legend()

		plt.savefig(f'{self.out_folder}/{self.prefix}-error-norms.png')
	
	def get_error_norm(self):
		plt.cla()

		abs_error = np.zeros(self.n)

		for i in range(self.n):
			errors = self.read_file('error', i)[2][:, 2]
			abs_error[i] = np.linalg.norm(errors)/len(errors)
		
		# abs_error = np.array([np.linalg.norm(np.sqrt(np.abs(self.read_file('error', i)[2][:,2])))/len(self.read_file('error', i)[2][:, 2]) for i in range(self.n)])
		return np.linalg.norm(abs_error) / len(abs_error)

	def plot_triangle(self, ax, start, gradient):
		vertices = [start, start * (1e2, 1), start * (1e2, 10**(2+gradient)) + (0, start[1]), (0,0)]
		codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]

		path = Path(vertices, codes)
		patch = PathPatch(path, facecolor=None, fill=False, edgecolor='b')

		ax.add_patch(patch)