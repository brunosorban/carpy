import numpy as np
from Function import Function

class Vehicle:
    def __init__(self, Tire, vehicle_mass, Ixx, Iyy, Izz, lf, lr, wf, wr, af, cd, CG_height, rho=1.225):
        self.Tire = Tire
        self.n_tires = 4
        self.vehicle_mass = vehicle_mass
        self.vehicle_equivalent_mass = self.vehicle_mass + self.n_tires * Tire.Jzz / Tire.tire_radius**2
        self.lf = lf
        self.lr = lr
        self.wf = wf
        self.wr = wr
        self.af = af
        self.cd = cd
        self.rho = rho
        self.CG_height = CG_height
        self.h = CG_height - Tire.tire_radius
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.K_arf = 0
        self.K_arr = 0
        self.Ff0 = self.lr / (self.lf + self.lr) * self.vehicle_mass * 9.81
        self.Fr0 = self.lf / (self.lf + self.lr) * self.vehicle_mass * 9.81

    def set_suspension(self, K_sf, K_sr, C_sf, C_sr):
        self.K_sf = K_sf
        self.K_sr = K_sr
        self.C_sf = C_sf
        self.C_sr = C_sr

    def set_anti_roll_bar(self, d, a, b, G, position):
        K_ar = G * (np.pi * d**4 / 32) * b / a**2
        if position == 'f': self.K_arf = K_ar
        elif position == 'r': self.K_arr = K_ar
        else: return print('Please insert a valid position. Position must be "f" or "r"')
        return self.set_suspension(self.K_sf, self.K_sr, self.C_sf, self.C_sr)

    def get_vertical_load(self, phi, theta, phi_dot, theta_dot):
        F10 = F20 = self.Ff0
        F30 = F40 = self.Fr0
        F1 = F10 + self.K_sf * (-self.wf * np.sin(phi) + self.lf * np.sin(theta)) + self.C_sf * (-self.wf * np.sin(phi_dot) + self.lf * np.sin(theta_dot))
        F2 = F20 + self.K_sf * (self.wf * np.sin(phi) + self.lf * np.sin(theta)) + self.C_sf * (self.wf * np.sin(phi_dot) + self.lf * np.sin(theta_dot))
        F3 = F30 + self.K_sf * (-self.wf * np.sin(phi) - self.lf * np.sin(theta)) + self.C_sf * (-self.wf * np.sin(phi_dot) - self.lf * np.sin(theta_dot))
        F4 = F40 + self.K_sf * (self.wf * np.sin(phi) - self.lf * np.sin(theta)) + self.C_sf * (self.wf * np.sin(phi_dot) - self.lf * np.sin(theta_dot))
        return F1, F2, F3, F4

    def get_steer_wbw(self, delta_sw):
        # return the ackerman geometry wheel by wheel
        if delta_sw == 0: 
            return 0, 0, 0, 0
        else:
            a = self.lf + self.lr
            s = self.wf
            r = a / np.tan( abs(delta_sw) )
            d1 = np.arctan(a / r)
            d2 = np.arctan(a / (r + s))
            if delta_sw > 0:
                return d1, d2, 0, 0
            else:
                return -d2, -d1, 0, 0