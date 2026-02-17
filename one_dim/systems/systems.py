#!/usr/bin/env python3

from typing import TextIO, Self

from dataclasses import dataclass, fields
from enum import StrEnum, auto

@dataclass
class Base:
	def to_file(self, fname: str) -> None:
		with open(fname, 'w') as f:
			for k in fields(self):
				f.write(f'{k.name}={self.__getattribute__(k.name)}\n')

	@classmethod
	def from_file(cls, fname: str) -> Self:
		clsinfo = fields(cls)
		obj = {k.name: None for k in clsinfo}
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
			water,
			ice,
			334000.0,
			0.011587153186771986
		)

		return system

class AnalysisType(StrEnum):
	standard = auto()
	dx = auto()
	dt = auto()
	sensitivity = auto()

@dataclass
class SimParams(Base):
	xs: int | list[int]
	ts: int | list[int]
	nxs: tuple[int, int] | list[tuple[int, int]]
	dts: float | list[float]
	tstart: float
	tend: float
	write_freq: int
	analysis_type: AnalysisType

@dataclass
class Params(Base):
	x: int
	t: int
	nx1: int
	nx2: int
	dt: float
	tstart: float
	tend: float
	write_freq: int

	@property
	def title(self) -> str:
		print(self.nx1)
		return f'x={self.x} t={self.t} nx=({self.nx1}, {self.nx2}) dt={self.dt}'
	
	@property
	def nx(self) -> int:
		return self.nx1 + self.nx2