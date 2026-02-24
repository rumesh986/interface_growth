from typing import TextIO, Self
from dataclasses import dataclass

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
	water_ice = System(
		DefaultMaterials.ice,
		DefaultMaterials.water,
		334000.0,
		0.011587153186771986
	)