import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pybamm

model = pybamm.lithium_ion.SPM()
param = pybamm.ParameterValues("Prada2013") # preciso alterar somente os que servem para SPM.

print(30*'=')
print('PARAMETROS PRADA 2013')
print(30*'=')
for chave, valor in param.items():
    print(f'{chave} ----- {valor}')

print(f'Quantidade de parâmetros: {len(param)}')


referencia = pybamm.ParameterValues("OKane2022")

print(30*'=')
print('PARAMETROS OKane2022')
print(30*'=')
print(f'Quantidade de parâmetros: {len(referencia)}')
for chave, valor in referencia.items():
    if chave not in param:
        print(f'{chave} ---- {valor}')
chaves_degradacao = [

    # SEI Growth
    "Ratio of lithium moles to SEI moles",
    "SEI partial molar volume [m3.mol-1]",
    "SEI resistivity [Ohm.m]",
    "SEI solvent diffusivity [m2.s-1]",
    "Bulk solvent concentration [mol.m-3]",
    "Initial SEI thickness [m]",
    "Initial SEI on cracks thickness [m]",
    "SEI growth activation energy [J.mol-1]",

    # Lithium Plating
    "Lithium metal partial molar volume [m3.mol-1]",
    "Lithium plating kinetic rate constant [m.s-1]",
    "Exchange-current density for plating [A.m-2]",
    "Exchange-current density for stripping [A.m-2]",
    "Initial plated lithium concentration [mol.m-3]",
    "Typical plated lithium concentration [mol.m-3]",
    "Lithium plating transfer coefficient",
    "Dead lithium decay constant [s-1]",
    "Dead lithium decay rate [s-1]",

    # Particle Mechanics + Cracking
    "Negative electrode Poisson's ratio",
    "Negative electrode Young's modulus [Pa]",
    "Negative electrode reference concentration for free of deformation [mol.m-3]",
    "Negative electrode partial molar volume [m3.mol-1]",
    "Negative electrode volume change",
    "Negative electrode initial crack length [m]",
    "Negative electrode initial crack width [m]",
    "Negative electrode number of cracks per unit area [m-2]",
    "Negative electrode Paris' law constant b",
    "Negative electrode Paris' law constant m",
    "Negative electrode cracking rate",

    # Loss of Active Material
    "Negative electrode LAM constant proportional term [s-1]",
    "Negative electrode LAM constant exponential term",
    "Negative electrode critical stress [Pa]",
]





