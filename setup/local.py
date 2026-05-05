import os
import subprocess
from datetime import datetime
from time import sleep

from concurrent.futures import ProcessPoolExecutor

from setup.build import Build
from systems.params import SimParams, Params, AnalysisType

class RunLocal(Build):
	def run(self, params: Params, *args, **kwargs) -> None:
		try:
			dname = f'{self.wd}/{self.results_dir}/{params.directory}'

			if kwargs['clear_results'] and os.path.exists(dname):
				shutil.rmtree(dname)
			os.mkdir(dname)

			params.to_file(f'{dname}/params')
			self.system.to_file(f'{dname}/system')

			kwargs['dname'] = dname

			cmd_args = self.process_args(params, *args, **kwargs)
			outfile = f'{dname}/stdout'

			with open(outfile, 'w') as f:
				cmd = [f'./{self.bins[(params.x, params.t)]}', *cmd_args]
				f.write(f'COMMAND: {cmd}\n\n')
				f.flush()
				subprocess.run(
					cmd,
					check=True,
					stderr=subprocess.STDOUT,
					stdout=f
				)
		except subprocess.CalledProcessError as e:
			print(f'run failed with params: {params}')
			with open(outfile, 'r') as f:
				print(f.read())
			raise

	def run_all(self, *args, **kwargs) -> None:
		num_workers = self.params.num_jobs
		if (num_workers > self.nthreads):
			num_workers = self.nthreads
		
		os.chdir(self.wd)
		self.params.to_file(f'{self.results_dir}/params')
		self.system.to_file(f'{self.results_dir}/system')

		futures = {}
		with ProcessPoolExecutor(max_workers=num_workers) as executor:
			match self.params.analysis_type:
				case AnalysisType.standard | AnalysisType.dx | AnalysisType.dt:
					for params in self.params.jobs:
						print(params)

						futures[params] = executor.submit(self.run, params, *args, **kwargs)
					
				case _:
					raise NotImplementedError()

			while len(futures) > 0:
				dellist = []
				for params, future in futures.items():
					if future.done():
						print(f'Completed {params}')
						if future.exception() is None:
							dellist.append(params)
						else:
							raise future.exception()
				
				for params in dellist:
					del futures[params]

				sleep(5)