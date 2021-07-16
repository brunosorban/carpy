import numpy as np
from Function import Function
from TMEasy import Tire

class Vehicle:
    def __init__(self, tire, vehicle_mass, Izz, lf, lr, wf, wr):
        self.Tire = tire
        self.n_tires = 4
        self.vehicle_mass = vehicle_mass
        self.vehicle_equivalent_mass = self.vehicle_mass + self.n_tires * tire.Jzz / tire.tire_radius**2
        self.lf = lf
        self.lr = lr
        self.wf = wf
        self.wr = wr
        # self.CG_height = CG_height
        # self.load_mass = load_mass
        # self.load_height = load_height
        # self.ab = ab
        # self.Ixx = Ixx
        # self.Iyy = Iyy
        self.Izz = Izz
        # self.Cphi = Cphi
        # self.Kphi = Kphi

        self.Fz = self.vehicle_mass * 9.81