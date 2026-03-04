import os
import sys

from systems.params import AnalysisType
from visualisation.visualisation import Visualiser, PltStyle

if __name__ == '__main__':
	os.chdir(sys.argv[2])

	vis = Visualiser(
		prefix=sys.argv[1],
		analysis_type=AnalysisType.standard,
		style=PltStyle.report
	)

	# vis.default_run()
	# vis.plot_mesh(vis.runs[0], 1)
	# vis.plot_interface(vis.runs[0])
	vis.anim_profile(vis.runs[0], False)
