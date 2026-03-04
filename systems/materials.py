from typing import TextIO, Self
from dataclasses import dataclass, fields

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

@dataclass
class MaterialSystem:
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
	
	@classmethod
	def from_file(cls, fname: str) -> Self:
		solid = {}
		liquid = {}
		L = None
		De = None

		material_types = {k.name: k.type for k in fields(Material)}

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
					case 'L': 
						L = float(v)
					case 'De': 
						De = float(v)
					case _: raise ValueError('Unknown key in system file')

		if L is None or De is None:
			raise Exception("Could not find all required information in system file")

		return cls(Material(**solid), Material(**liquid), L, De)

	@property
	def args(self) -> list[str]:
		return [
			*self.solid.args(1),
			*self.liquid.args(2),
			'--L', str(self.L),
			'--De', str(self.De)
		]

class DefaultMaterials:
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

class DefaultMaterialSystems:
	water_ice = MaterialSystem(
		DefaultMaterials.ice,
		DefaultMaterials.water,
		334000.0,
		0.011587153186771986
	)