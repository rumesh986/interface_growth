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
class SimDomain(_Base):
	x0: float
	x1: float
	x2: float
	y0: float
	y1: float

	def to_file(self, file: TextIO) -> None:
		for k in fields(self):
			file.write(f'domain.{k.name}={self.__getattribute__(k.name)}\n')
	
	@property
	def args(self) -> list[str]:
		return [
			'--x0', str(self.x0),
			'--x1', str(self.x1),
			'--x2', str(self.x2),
			'--y0', str(self.y0),
			'--y1', str(self.y1),
		]

@dataclass(frozen=True)
class Params(_Base):
	x: int
	t: int
	nx1: int
	nx2: int
	ny: int
	domain: SimDomain
	dt: float
	tstart: float
	tend: float
	write_freq: int
	k: int
	A: float
	v: float
	nx: int = field(init=False)

	def __post_init__(self):
		object.__setattr__(self, 'nx', self.nx1 + self.nx2)
	
	def to_file(self, fname: str) -> None:
		with open(fname, 'w') as f:
			for k in fields(self):
				if k.type is SimDomain:
					self.domain.to_file(f)
					continue
				f.write(f'{k.name}={self.__getattribute__(k.name)}\n')
	
	@classmethod
	def from_file(cls, fname: str) -> Self:
		clsinfo = fields(cls)
		obj = {}
		cls_types = {k.name: k.type for k in clsinfo if k.init}
		domain_types = {k.name: k.type for k in fields(SimDomain)}
		domain = {}

		with open(fname, 'r') as f:
			while line := f.readline():
				k, v = line.split('=')
				match k:
					case tag if tag.startswith('domain'):
						_, k2 = tag.split('.')
						domain[k2] = domain_types[k2](v.strip())
					case key if key in cls_types.keys():
						obj[k] = cls_types[k](v)
		
		obj['domain'] = SimDomain(**domain)
		
		return cls(**obj)

	@property
	def args(self) -> list[str]:
		return [
			'--nx1', str(self.nx1),
			'--nx2', str(self.nx2),
			'--ny', str(self.ny),
			'--dt', str(self.dt),
			*self.domain.args,
			'--tstart', str(self.tstart),
			'--tend', str(self.tend),
			'--write-freq', str(self.write_freq),
			'--k', str(self.k),
			'--A', str(self.A),
			'--v', str(self.v)
		]

	@property
	def args_cmd(self) -> str:
		return ' '.join(self.args)

	@property
	def title(self) -> str:
		return f'x={self.x} t={self.t} nx=({self.nx1}, {self.nx2}) ny={self.ny} dt={self.dt}'
	
	@property
	def short_title(self) -> str:
		return f'{self.x}n{self.nx1}+{self.nx2}_{self.ny}_{self.t}t{self.dt}'
	
	@property
	def directory(self) -> str:
		return f'{self.x}n{self.nx1}+{self.nx2}_{self.ny}_{self.t}t{self.dt:.2e}_{self.k}k{self.A}_v{self.v}'#_{self.__hash__()}'

@dataclass(frozen=True)
class SimParams(_Base):
	xs: list[int]
	ts: list[int]
	nxs: list[tuple[int, int]]
	nys: list[int]
	domains: list[SimDomain]
	dts: list[float]
	ks: list[int]
	As: list[float]
	vs: list[float]
	tstart: float
	tend: float
	write_freq: int
	analysis_type: AnalysisType

	def __len__(self):
		ret = 1
		for k in fields(self):
			attr = self.__getattribute__(k.name)
			if isinstance(attr, list):
				ret *= len(attr)
		return ret

	@property
	def jobs(self) -> Iterator[Params]:
		for x, t, (nx1, nx2), ny, dt, k, A, domain, v in product(self.xs, self.ts, self.nxs, self.nys, self.dts, self.ks, self.As, self.domains, self.vs):
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
				domain=domain,
				dt=dt,
				tstart=self.tstart,
				tend=self.tend,
				write_freq=wf,
				k=k,
				A=A,
				v=v
			)