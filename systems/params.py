from itertools import product
from enum import StrEnum, auto
from typing import Self, Iterator
from dataclasses import dataclass, fields, field

class AnalysisType(StrEnum):
	standard = auto()
	dx = auto()
	dt = auto()
	sensitivity = auto()

@dataclass(frozen=True)
class _Base:
	def to_file(self, fname: str) -> None:
		with open(fname, 'w') as f:
			for k in fields(self):
				f.write(f'{k.name}={self.__getattribute__(k.name)}\n')

	@classmethod
	def from_file(cls, fname: str) -> Self:
		clsinfo = fields(cls)
		obj = {}
		types = {k.name: k.type for k in clsinfo if k.init}

		with open(fname, 'r') as f:
			while line := f.readline():
				k, v = line.split('=')
				if k in types.keys():
					obj[k] = types[k](v)
		
		return cls(**obj)

	def __str__(self) -> str:
		string = ""
		for k in fields(self):
			string += f"{k.name}={self.__getattribute__(k.name)} "
		
		return string


@dataclass(frozen=True)
class Params(_Base):
	x: int
	t: int
	nx1: int
	nx2: int
	ny: int
	dt: float
	tstart: float
	tend: float
	write_freq: int
	k: int
	A: float
	nx: int = field(init=False)

	def __post_init__(self):
		object.__setattr__(self, 'nx', self.nx1 + self.nx2)

	@property
	def args(self) -> list[str]:
		return [
			'--nx1', str(self.nx1),
			'--nx2', str(self.nx2),
			'--ny', str(self.ny),
			'--dt', str(self.dt),
			'--tstart', str(self.tstart),
			'--tend', str(self.tend),
			'--write-freq', str(self.write_freq),
			'--k', str(self.k),
			'--A', str(self.A)
		]

	@property
	def args_cmd(self) -> str:
		return (
			f" --nx1 {self.nx1}"
			f" --nx2 {self.nx2}"
			f" --ny {self.ny}"
			f" --dt {self.dt}"
			f" --tstart {self.tstart}"
			f" --tend {self.tend}"
			f" --write-freq {self.write_freq}"
			f" --k {self.k}"
			f" --A {self.A}"
		)

	@property
	def title(self) -> str:
		print(self.nx1)
		return f'x={self.x} t={self.t} nx=({self.nx1}, {self.nx2}) ny={self.ny} dt={self.dt}'
	
	@property
	def short_title(self) -> str:
		return f'{self.x}n{self.nx1}+{self.nx2}_{self.ny}_{self.t}t{self.dt}'
	
	@property
	def directory(self) -> str:
		return f'{self.x}n{self.nx1}+{self.nx2}_{self.ny}_{self.t}t{self.dt:.2e}_{self.k}k{self.A}'

@dataclass(frozen=True)
class SimParams(_Base):
	xs: list[int]
	ts: list[int]
	nxs: list[tuple[int, int]]
	nys: list[int]
	dts: list[float]
	ks: list[int]
	As: list[float]
	tstart: float
	tend: float
	write_freq: int
	analysis_type: AnalysisType

	def __len__(self):
		return len(self.xs) * len(self.ts) * len(self.nxs) * len(self.nys) * len(self.dts) * len(self.ks) * len(self.As)

	@property
	def jobs(self) -> Iterator[Params]:
		for x, t, (nx1, nx2), ny, dt, k, A in product(self.xs, self.ts, self.nxs, self.nys, self.dts, self.ks, self.As):
			if self.analysis_type == AnalysisType.dt:
				wf = int(max(self.dts) / dt)
			else:
				wf = self.write_freq

			yield Params(
				x=x,
				t=t,
				nx1=nx1,
				nx2=nx2,
				ny=ny,
				dt=dt,
				tstart=self.tstart,
				tend=self.tend,
				write_freq=wf,
				k=k,
				A=A
			)