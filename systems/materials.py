import os
from importlib import resources
from typing import TextIO, Self, ClassVar
from dataclasses import dataclass, fields, InitVar, field

import numpy as np
from scipy.special import erf
from scipy.optimize import newton

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
		for k in fields(self):
			file.write(f'{tag}.{k.name}={self.__getattribute__(k.name)}\n')
	
	def args(self, index: int) -> list[str]:
		return [
			f'--k{index}', str(self.k),
			f'--rho{index}', str(self.rho),
			f'--cp{index}', str(self.cp),
		]
	
	@property
	def alpha(self):
		return self.rho * self.cp
	@property
	def beta(self):
		return self.k

	@property
	def D(self):
		return self.beta / self.alpha

@dataclass
class MaterialSystem:
	solid: Material
	liquid: Material
	L: float
	gamma: float

	Ts: float
	Tm: float
	Tl: float

	De: float = field(default=None)

	def __post_init__(self):
		if self.De is None:
			self.De = self._calculate_De()

	def to_file(
		self,
		fname: str
	) -> None:
		with open(fname, 'w') as f:
			self.solid.to_file(f, 'solid')
			self.liquid.to_file(f, 'liquid')
			for k in fields(self):
				if k.type is float:
					f.write(f'{k.name}={self.__getattribute__(k.name)}\n')
	
	@classmethod
	def from_file(cls, fname: str) -> Self:
		solid = {}
		liquid = {}
		L = None
		De = None
		gamma = None

		material_types = {k.name: k.type for k in fields(Material)}
		cls_types = {k.name: k.type for k in fields(cls)}
		obj = {}

		with open(fname, 'r') as f:
			while line := f.readline():
				k,v = line.split('=')

				match k:
					case tag if tag.startswith('solid'):
						_, k2 = tag.split('.')
						solid[k2] = material_types[k2](v.strip())
					case tag if tag.startswith('liquid'):
						_, k2 = tag.split('.')
						liquid[k2] = material_types[k2](v.strip())
					# case key if key in ['L', 'De', 'Ts', 'Tm', 'Tl']:
					case key if cls_types[key] is float:
						obj[key] = float(v)
					case _: raise ValueError('Unknown key in system file')

		obj['solid'] = Material(**solid)
		obj['liquid'] = Material(**liquid)

		return cls(**obj)

	@property
	def args(self) -> list[str]:
		return [
			*self.solid.args(1),
			*self.liquid.args(2),
			'--L', str(self.L),
			'--De', str(self.De),
			'--gamma', str(self.gamma)
		]
	
	@classmethod
	def list_default_materials(cls) -> List[str]:
		material_dir = resources.files('systems').joinpath('default_materials')
		return os.listdir(material_dir)
	
	@classmethod
	def default_material(cls, name: str) -> Self:
		material_dir = resources.files('systems').joinpath('default_materials')
		return cls.from_file(f'{material_dir}/{name}')
	
	def _calculate_De(self):
		alpha = self.liquid.alpha / self.solid.alpha
		beta = self.liquid.beta / self.solid.beta
		D = beta / alpha

		def func(x):
			if self.Tl > self.Tm:
				# standard case 
				delT = self.Tm - self.Ts
				Tl = (self.Tl - self.Tm) / delT
				St = self.L / (self.solid.cp * delT)

				return 0.5 * St * np.sqrt(np.pi * x) - np.exp(-0.25 * x) / erf(0.5*np.sqrt(x)) + Tl * np.sqrt(alpha * beta) * np.exp(-0.25*x/D) / (1.0 - erf(0.5 * np.sqrt(x/D)))
			else:
				# undercooled case
				delT = self.Tm - self.Tl
				Ts = (self.Ts - self.Tm) / delT
				St = self.L / (self.solid.cp * delT)

				return 0.5 * St * np.sqrt(np.pi * x) + Ts * np.exp(-0.25 * x) / erf(0.5*np.sqrt(x)) - np.sqrt(alpha * beta) * np.exp(-0.25*x/D) / (1.0 - erf(0.5 * np.sqrt(x/D)))
		
		De  = newton(func, self.solid.D)

		if not np.isclose(func(De), 0.0):
			raise Exception("Unable to get analytical De for this system, check parameters")
		
		return De