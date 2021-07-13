import numpy as np
from Function import Function

class Vehicle:
    def __init__(self, Tire, vehicle_mass, CG_height, load_mass, load_height, ab, Ixx, Iyy, Izz, a, b, l, Cphi, Kphi):
        self.Tire = Tire
        self.vehicleMass = vehicle_mass
        self.CG_height = CG_height
        self.load_mass = load_mass
        self.load_height = load_height
        self.ab = ab
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.a = a
        self.b = b
        self.l = l
        self.Cphi = Cphi
        self.Kphi = Kphi