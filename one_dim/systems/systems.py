#!/usr/bin/env python3

from itertools import product
from typing import TextIO, Self, Iterator

from dataclasses import dataclass, fields, field
from enum import StrEnum, auto

@dataclass(frozen=True)
class Base:
	def to_file(self, fname: str) -> None:
		with open(fname, 'w') as f:
			for k in fields(self):
				f.write(f'{k.name}={self.__getattribute__(k.name)}\n')

	@classmethod
	def from_file(cls, fname: str) -> Self:
		clsinfo = fields(cls)
		obj = {}
		types = {k.name: k.type for k in clsinfo}

		with open(fname, 'r') as f:
			while line := f.readline():
				k, v = line.split('=')
				obj[k] = types[k](v)
		
		return Params(**obj)

	def __str__(self) -> str:
		string = ""
		for k in fields(self):
			string += f"{k.name}={self.__getattribute__(k.name)} "
		
		return string

@dataclass(frozen=True)
class Material:
	name: str
	k: float
	rho: float
	cp: float

	def to_file(
		self,
		file: TextIO,
		tag: str
	) -> None:
		file.write(f'{tag}.k={self.k}\n')
		file.write(f'{tag}.rho={self.rho}\n')
		file.write(f'{tag}.cp={self.cp}\n')
	
	def args(self, index: int) -> list[str]:
		return [
			f'--k{index}', str(self.k),
			f'--rho{index}', str(self.rho),
			f'--cp{index}', str(self.cp),
		]

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
	
	def from_file(self) -> Self:
		water = Material(
			name="water",
			k=0.55575,
			rho=999.89,
			cp=4220.0
		)

		ice = Material(
			name="ice",
			k=2.2,
			rho=916.2,
			cp=2050.0
		)

		system = System(
			ice,
			water,
			334000.0,
			0.011587153186771986
		)

		return system
	
	@property
	def args(self) -> list[str]:
		return [
			*self.solid.args(1),
			*self.liquid.args(2),
			'--L', str(self.L),
			'--De', str(self.De)
		]

class AnalysisType(StrEnum):
	standard = auto()
	dx = auto()
	dt = auto()
	sensitivity = auto()

@dataclass(frozen=True)
class Params(Base):
	x: int
	t: int
	nx1: int
	nx2: int
	dt: float
	tstart: float
	tend: float
	write_freq: int
	nx: int = field(init=False)

	def __post_init__(self):
		object.__setattr__(self, 'nx', self.nx1 + self.nx2)

	@property
	def args(self) -> list[str]:
		return [
			'--nx1', str(self.nx1),
			'--nx2', str(self.nx2),
			'--dt', str(self.dt),
			'--tstart', str(self.tstart),
			'--tend', str(self.tend),
			'--write-freq', str(self.write_freq)
		]

	@property
	def title(self) -> str:
		print(self.nx1)
		return f'x={self.x} t={self.t} nx=({self.nx1}, {self.nx2}) dt={self.dt}'
	
	@property
	def short_title(self) -> str:
		return f'{self.x}n{self.nx1}+{self.nx2}_{self.t}t{self.dt}'

@dataclass(frozen=True)
class SimParams(Base):
	xs: list[int]
	ts: list[int]
	nxs: list[tuple[int, int]]
	dts: list[float]
	tstart: float
	tend: float
	write_freq: int
	analysis_type: AnalysisType
	num_jobs: int = field(init=False)

	def __post_init__(self):
		if type(self.xs) is int:
			object.__setattr__(self, 'xs', [self.xs])
		
		if type(self.ts) is int:
			object.__setattr__(self, 'ts', [self.ts])
		
		if type(self.nxs) is int:
			object.__setattr__(self, 'nxs', [(self.nxs, self.nxs)])
		elif type(self.nxs) is tuple:
			object.__setattr__(self, 'nxs', [self.nxs])

		if type(self.dts) is float:
			object.__setattr__(self, 'dts', [self.dts])
		
		object.__setattr__(self, 'num_jobs', len(self.xs) * len(self.ts) * len(self.nxs) * len(self.dts))

	@property
	def jobs(self) -> Iterator[Params]:
		for x, t, (nx1, nx2), dt in product(self.xs, self.ts, self.nxs, self.dts):
			if self.analysis_type == AnalysisType.dt:
				wf = int(max(self.dts) / dt)
			else:
				wf = self.write_freq

			yield Params(
				x,
				t,
				nx1,
				nx2,
				dt,
				self.tstart,
				self.tend,
				wf
			)