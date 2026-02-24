import os
import subprocess
from datetime import datetime
from time import sleep

from concurrent.futures import ProcessPoolExecutor

from systems.params import SimParams, Params, AnalysisType
from systems.materials import MaterialSystem

class Build:
	def __init__(self, 
		prog: str,
		params: SimParams,
		source_dir: str,
		system: MaterialSystem,
		dname: str | None = None,
		workspacesf: str = 'workspaces',
		results_dir: str = 'RESLT',
		nthreads: int = os.cpu_count()
	) -> None:
		self.prog = prog
		self.params = params
		self.source_dir = source_dir
		self.system = system
		self.results_dir = results_dir
		self.nthreads = nthreads

		if dname is None:
			dname = f'{self.prog}_{datetime.today().strftime("%Y-%m-%d_%H:%M:%S")}'

		self.od = os.getcwd()
		self.wd = f'{self.od}/{workspacesf}/{prog}/{dname}'

		if not os.path.exists(f'{self.od}/{workspacesf}'):
			os.mkdir(f'{self.od}/{workspacesf}')
		
		if not os.path.exists(f'{self.od}/{workspacesf}/{prog}'):
			os.mkdir(f'{self.od}/{workspacesf}/{prog}')

	def build_all(self) -> None:
		os.chdir(self.source_dir)

		self.bins = {}

		command = ['cmake', '--build', 'build', '--target']
		for x in self.params.xs:
			for t in self.params.ts:
				bin_name = f'{self.prog}_{x}_{t}'
				command.append(bin_name)
				self.bins[(x, t)] = bin_name
		
		result = subprocess.run(command)
		if result.returncode != 0:
			raise Exception(f"Build failed with error code {result.returncode}")

		if not os.path.exists(self.wd):
			os.mkdir(self.wd)

		for item in self.bins.values():
			os.rename(f'build/{item}', f'{self.wd}/{item}')
	
	def process_args(self, params: Params, *args, **kwargs) -> list[str]:
		ret = [
			*params.args,
			*self.system.args
		]

		if len(args) > 0:
			print(args)
			print(len(args))
			ret.extend(args)
		
		for k, v in kwargs.items():
			print(f'{k}={v}')
			ret.extend([f'--{k}', str(v)])

		return ret