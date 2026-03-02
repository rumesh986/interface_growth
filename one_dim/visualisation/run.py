from visualisation.run_oned import Run_OneD
from visualisation.run_twod import Run_TwoD

class Run:
	def __new__(cls, *args, **kwargs):
		with open(f'{args[0]}/dim') as f:
			dim = int(f.readline())

		match dim:
			case 1: return Run_OneD(*args, **kwargs)
			case 2: return Run_TwoD(*args, **kwargs)
			case _: raise NotImplementedError()
