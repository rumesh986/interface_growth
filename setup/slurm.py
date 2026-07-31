import os
from datetime import timedelta
from typing import TextIO
from dataclasses import dataclass, fields, field

from setup.build import Build

@dataclass(frozen=True)
class SlurmParams():
	nodes: int
	ntasks_per_node: int
	ncpus_per_task: int
	memory_per_cpu: int
	max_time: timedelta

	@property
	def _time_str(self):
		h = self.max_time.seconds // 3600
		m = (self.max_time.seconds - 3600*h) // 60
		s = (self.max_time.seconds - 3600*h - 60*m) // 60

		# account for number of days which is stored separately from seconds for some reason
		h += self.max_time.days * 24

		return f'{h:02}:{m:02}:{s:02}'

	@property
	def slurm_args(self):
		return [
			f'#SBATCH --nodes={self.nodes}\n',
			f'#SBATCH --ntasks-per-node={self.ntasks_per_node}\n',
			f'#SBATCH --cpus-per-task={self.ncpus_per_task}\n',
			f'#SBATCH --mem-per-cpu={self.memory_per_cpu}\n',
			f'#SBATCH --time={self._time_str}\n',
		]

class RunSlurm(Build):
	def prepare_slurmfile(self, fname: str, slurm_params: SlurmParams):
		with open(fname, 'w') as f:
			f.write("#!/bin/bash\n\n")
			f.writelines(slurm_params.slurm_args)
			# f.write(f'#SBATCH --output={self.wd}/stdout-%A-%a.out\n')
			if len(self.params) > 1:
				f.write(f'#SBATCH --array=0-{len(self.params)-1}\n')
			f.write("\nsource scrtp_modules\n\n")
			f.write(f"cd {self.wd}\n\n")

			for param in self.params.jobs:
				dname =f'{self.wd}/{self.results_dir}/{param.directory}'
				os.mkdir(dname)
				param.to_file(f'{dname}/params')
				self.system.to_file(f'{dname}/system')
				f.write(f"srun ./{self.bins[(param.x, param.t)]} {param.args_cmd} --dname {dname}\n\n")			


	def prepare_array_job(self, fname: str, slurm_params: SlurmParams):
		if len(self.bins) > 1:
			raise NotImplementedError("Array jobs with multiple executables have not been implemented yet")

		with open('scrtp_modules', 'r') as f:
			module_cmd = f.read()

		os.chdir(self.wd)

		with open(fname, 'w') as f:
			f.write('#!/bin/bash\n\n')
			f.writelines(slurm_params.slurm_args)
			f.write(f'#SBATCH --array=0-{len(self.params)-1}\n\n')
			f.write(f'{module_cmd}\n\n')
			# f.write('source scrtp_modules\n\n')
			# f.write("IFS=$'\\n' ARGS=$(cat args.inp)\n\n")
			f.write(f"IFS=$'\\n' read -a ARGS -d '' <<< `cat args.inp`\n\n")
			f.write(f'srun ./{list(self.bins.values())[0]} ${{ARGS[$SLURM_ARRAY_TASK_ID]}}')


		with open('args.inp', 'w') as f:
			for param in self.params.jobs:
				dname = f'{self.wd}/{self.results_dir}/{param.directory}'
				os.mkdir(dname)
				param.to_file(f'{dname}/params')
				self.system.to_file(f'{dname}/system')
				f.write(f'{param.args_cmd} --dname {dname}\n')




	def prepare_directory(self):
		os.chdir(self.wd)
		for param in self.params.jobs:
			dname =f'{self.wd}/{self.results_dir}/{param.directory}'
			param.to_file(f'{dname}/params')
			self.system.to_file(f'{dname}/system')

	def test(self):
		# os.chdir(self.wd)
		self.params.to_file(f'{self.results_dir}/params')
		self.system.to_file(f'{self.results_dir}/system')