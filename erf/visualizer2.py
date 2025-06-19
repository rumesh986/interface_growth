import matplotlib.pyplot as plt
import matplotlib.animation as anim

import numpy as np
import pandas as pd

import timeit
from visualizer import Visualizer

class _Config:
	def __init__(self, file):
		self.nx = 0
		self.x_order = 0
		self.dt = 0.0
		self.t_shift = 0.0
		self.t_order = 0
		self.t_steps = 0
		self.write_freq = 1

		with open(file) as f:
			while line := f.readline():
				k, v = line.split('=')

				match k:
					case 'nx': self.nx = int(v)
					# case 'nx1': self.nx1 = int(v)
					# case 'nx2': self.nx2 = int(v)
					case 'x_order': self.x_order = int(v)
					case 'dt': self.dt = float(v)
					case 't_shift': self.t_shift = float(v)
					case 't_order': self.t_order = int(v)
					case 't_steps': self.t_steps = int(v)
					case 'write_freq': self.write_freq = int(v)
					# case _: raise Exception(f"Unexpected key in config {k}")

	def __str__(self):
		return f'{self.nx = }\n{self.x_order = }\n{self.dt = }\n{self.t_shift = }\n{self.t_order = }\n{self.t_steps = }\n{self.write_freq = }'
	
	def title(self):
		return f'nx={self.nx}, dt={self.dt}, X({self.x_order}), T({self.t_order})'
 
class _RunData:
	def __init__(
		self, 
		inp_folder, 
		prefix, 
		solnf='soln', 
		exact_solnf='exact_soln', 
		errorf='error',
		configf='config',
		ext='dat'
	):
		self.inp_folder = inp_folder
		self.prefix = prefix
		self.solnf = solnf
		self.exact_solnf = exact_solnf
		self.errorf = errorf
		self.ext = ext

		self.config = _Config(f'{self.inp_folder}/{configf}')
		self.data = pd.DataFrame()

		try:
			self.times = np.fromfile(f'{self.inp_folder}/times.{self.ext}', sep='\n')
		except FileNotFoundError:
			raise Exception("Time information not availble")
			


	def _read_file(self, ftype, index):
		file = f'{self.inp_folder}/{ftype}{index}.{self.ext}'

		if ftype == self.errorf:
			cols = ['x', 'norm', 'error']
		else:
			cols = ['x', 'u']

		data = []
		with open(file, 'r') as f:
			while line := f.readline():
				if not line.startswith('ZONE'):
					data.append([float(x) for x in line.split()])

		return pd.DataFrame(data, columns=cols)

	def read_series(self, ftype):
		data = []
		for i, t in enumerate(self.times):
			_data = self._read_file(ftype, i)
			_data['time'] = t

			data.append(_data)
		
		return data
	
	def get_errors(self, error_norms=True, error_max=True):
		data = self.read_series(self.errorf)

		norms_ret = None
		max_ret = None
		if error_norms:
			norms_ret = self._get_error_norms(data)

		if error_max:
			max_ret = self._get_error_max(data)

		return norms_ret, max_ret


	def _get_error_norms(self, data):
		plt.cla()

		abs_errors = np.array([np.linalg.norm(x['error']) for x in data])
		abs_errors /= len(data[0]['error'])

		return abs_errors

	def _get_error_max(self, data):
		plt.cla()

		node = np.argmax(np.fabs(data[1]['error'][1:])) + 1
		
		errors = np.array([df['error'][node] for df in data])

		exact_solns = self.read_series(self.exact_solnf)
		exact_vals = np.array([df['u'][node] for df in exact_solns])
		mask = exact_vals != 0.0

		abs_errors = np.fabs(errors)
		rel_errors = abs_errors[mask] / np.fabs(exact_vals)[mask]

		return abs_errors, rel_errors

if __name__ =='__main__':
	run = _RunData('RESLT/2n1000_4t1.00e-01', 'erf2')

	print(run.config)
	rd = run.plot_errors()

	print(len(rd))

	# timer1 = timeit.Timer(stmt='for i in range(100): run.read_file("error", i)', setup='run = _RunData("RESLT/2n1000_4t1.00e-01", "erf2")', globals=globals())
	# print(timer1.timeit(100))

	vis = Visualizer("RESLT/2n1000_4t1.00e-01", "", "erf2", times="times.dat")
	vd  = vis.plot_errors()

	print(rd)
	print(vd)

	if np.array_equal(rd, vd):
		print("Were still golden")
	else:
		print("Damn got issues now")


	# for i in range(100):
	# 	rd = run.read_file('error', i)['error']
	# 	vd = vis.read_file('error', i)[2][:, 2]

	# 	print(f'[{i}] {rd.shape} {vd.shape}')

	# 	if not np.array_equal(rd, vd):
	# 		raise Exception(f'whoops {i}')
			
	# print("Its all good so far")

	# timer2 = timeit.Timer(stmt='for i in range(100): run.read_file("error", i)', setup='run = Visualizer("RESLT/2n1000_4t1.00e-01", "", "erf2", times="times.dat")', globals=globals())
	# print(timer2.timeit(100))