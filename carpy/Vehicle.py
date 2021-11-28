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
        self.Ff0 = self.lr / (2 * (self.lf + self.lr) ) * self.vehicle_mass * 9.81
        self.Fr0 = self.lf / (2 * (self.lf + self.lr) ) * self.vehicle_mass * 9.81

    def set_suspension(self, K_sf, K_sr, C_sf, C_sr):
        self.K_sf = K_sf
        self.K_sr = K_sr
        self.C_sf = C_sf
        self.C_sr = C_sr

    def set_anti_roll_bar(self, d, a, b, G, position):
        K_arz = G * (np.pi * d**4 / 32) * b / a**2
        if position == 'f': self.K_sf += K_arz
        elif position == 'r': self.K_sr += K_arz
        else: return print('Please insert a valid position. Position must be "f" or "r"')
        return None

    def get_vertical_load(self, z, vz, phi, theta, phi_dot, theta_dot, z1, z2, z3, z4, vz1, vz2, vz3, vz4, zc1, zc2, zc3, zc4, vzc1, vzc2, vzc3, vzc4):
        # Retrieve suspension properties
        Kt = self.Tire.cz
        Ct = 0
        K_sf = self.K_sf
        K_sr = self.K_sr
        C_sf = self.C_sf
        C_sr = self.C_sr
        
        # Get static forces distribution
        F10 = F20 = self.Ff0
        F30 = F40 = self.Fr0
        
        # Calculate sines used
        sphi = np.sin(phi)
        stheta = np.sin(theta)
        sphidot = np.sin(phi_dot)
        sthetadot = np.sin(theta_dot)

        # Calculate tire forces        
        Ft1 = F10 + Kt * (zc1 - z1) + Ct * (vzc1 - vz1)
        Ft2 = F20 + Kt * (zc2 - z2) + Ct * (vzc2 - vz2)
        Ft3 = F30 + Kt * (zc3 - z3) + Ct * (vzc3 - vz3)
        Ft4 = F40 + Kt * (zc4 - z4) + Ct * (vzc4 - vz4)
        
        # Calculate suspension forces
        Fs1 = F10 + K_sf * (z1 - z - self.wf/2 * sphi + self.lf * stheta) + C_sf * (vz1 - vz - self.wf/2 * sphidot + self.lf * sthetadot)
        Fs2 = F20 + K_sf * (z2 - z + self.wf/2 * sphi + self.lf * stheta) + C_sf * (vz2 - vz + self.wf/2 * sphidot + self.lf * sthetadot)
        Fs3 = F30 + K_sr * (z3 - z - self.wr/2 * sphi - self.lr * stheta) + C_sr * (vz3 - vz - self.wr/2 * sphidot - self.lr * sthetadot)
        Fs4 = F40 + K_sr * (z4 - z + self.wr/2 * sphi - self.lr * stheta) + C_sr * (vz4 - vz + self.wr/2 * sphidot - self.lr * sthetadot)
        return Ft1, Ft2, Ft3, Ft4, Fs1, Fs2, Fs3, Fs4

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