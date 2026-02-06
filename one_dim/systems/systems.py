#!/usr/bin/env python3

import typing

from dataclasses import dataclass
from enum import StrEnum, auto

@dataclass(frozen=True)
class Material:
	name: str
	k: float
	rho: float
	cp: float

	def to_file(
		self,
		file: typing.TextIO,
		tag: str
	) -> None:
		file.write(f'{tag}.k={self.k}\n')
		file.write(f'{tag}.rho={self.rho}\n')
		file.write(f'{tag}.cp={self.cp}\n')

@dataclass
class System:
	solid: Material
	liquid: Material
	L: float
	De: float

	def to_file(
		self,
		fname: str
	) -> None:
		with open(fname, 'w') as f:
			self.solid.to_file(f, 'solid')
			self.liquid.to_file(f, 'liquid')
			f.write(f'L={self.L}\n')
			f.write(f'De={self.De}\n')

class AnalysisType(StrEnum):
	standard = auto()
	dx = auto()
	dt = auto()
	sensitivity = auto()

@dataclass
class SimParams:
	xs: int | list[int]
	ts: int | list[int]
	nxs: int | list[int]
	dts: float | list[float]
	tstart: float
	tend: float
	write_freq: int
	type: AnalysisType
