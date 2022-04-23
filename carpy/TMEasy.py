from math import gamma
import numpy as np
from Function import Function

''' Here lies an inplementation of TMEasy tire model from George Rill and Abel
Castro. The implmentation is based on the Matlab code presented at Road Vehicle 
Dynamics, second edition.'''

class Tire:
    def __init__(self, radius, mass, Jzz, cz, dfx0, dfy0, Fxm, Fym, Sxm, Sym, Fxs, Fys, Sxs, Sys, Sy0, SyE, lamb, n2L0, frr):
        self.tire_radius = radius
        self.mass = mass
        self.Jzz = Jzz
        self.cz = cz
        self.dfx0 = dfx0
        self.dfy0 = dfy0
        self.Fxm = Fxm
        self.Fym = Fym
        self.Sxm = Sxm
        self.Sym = Sym
        self.Fxs = Fxs
        self.Fys = Fys
        self.Sxs = Sxs
        self.Sys = Sys
        self.Sy0 = Sy0
        self.SyE = SyE
        self.lamb = lamb
        self.n2L0 = n2L0
        self.frr = frr

        # Calculate tire parameters
        self.hsxn = Sxm / (Sxm + Sym) + Fxm / (dfx0 * (Fxm / dfx0 + Fym / dfy0))
        self.hsyn = Sym / (Sxm + Sym) + Fym / (dfy0 * (Fxm / dfx0 + Fym / dfy0))

    def __call__(self, *args):
        Sx, Sy, Sg, Sb = self.calculateSlip(*args)
        return self.get_tire_forces(Sx, Sy, Sg, Sb, args[6])
        
    def __tire_offset__(self, Sy):
        Sy = abs(Sy) # absolute value because n(ny) = n(-Sy)
        if Sy < self.SyE: # 4th - order polynomial 0 <= Sy < SyE
            a = -(2 * self.Sy0 ** 2 + (self.Sy0 + self.SyE) ** 2)
            b = 2 * (self.Sy0 + self.SyE) ** 2 / self.SyE
            c = -(2 * self.Sy0 + self.SyE) / self.SyE
            n2L = self.n2L0 * (1 + (Sy / (self.Sy0 * self.SyE)) ** 2 * (a + Sy * (b + Sy * c)))
        else:
            n2L = 0 # Full sliding
        return n2L
    
    def __tire_combined_forces__(self, S, df0, Fm, Sm, Fs, Ss):
        # Curve adaptation for too small df0
        if df0 > 0: 
            Smloc = max([2 * Fm / df0, Sm])
            Ssloc = Ss + Smloc - Sm
        else: return print("df0 must be greater than 0. Please insert a valid df0.")
        
        # Normal operating conditions
        if S > 0 and Smloc > 0:
            if S > Ssloc: # Full sliding
                F = Fs
                dF = F / S
            
            elif S < Smloc and Fm > 0: # Full adhesion. The force varies linearly with the slip
                p = df0 * Smloc / Fm - 2
                Sn = S / Smloc
                dF = df0 / (1 + (Sn + p) * Sn)
                F = dF * S
            else: # Adhesion and sliding. The function is composed of two parabolas or a cubic function
                a = (Fm / Smloc) ** 2 / (df0 * Smloc)
                Sstar = Smloc + (Fm - Fs) / (a * (Ssloc - Smloc))

                if Sstar <= Ss: 
                    if S <= Sstar: # Sm < S < Sstar
                        F = Fm - a * (S - Smloc)**2 
                    else: # Sstar < S < Ss
                        b = a * (Sstar - Smloc) / (Ssloc - Sstar)
                        F = Fs + b * (Ssloc - S)**2
                else: # A cubic function if the curve doesn's fit in a coherent way (Sstar must be < than Ss)
                    Sn = (S - Smloc) / (Ssloc - Smloc)
                    F = Fm - (Fm - Fs) * Sn**2 * (3 - 2 * Sn)
                
                dF = F / S
        else: 
            F = 0
            dF = df0
        return F, dF

    def calculateSlip(self, Vx, Vy, delta, gamma, omega_wheel, Wn, Fz):
        # Calculate slip parameters
        Vn = 0.01
        Rd = self.lamb * self.tire_radius + (1 - self.lamb) * (self.tire_radius - Fz / self.cz)
        L = 2 * np.sqrt(self.tire_radius * Fz / self.cz)
        Vxt = Vx * np.cos(delta) + Vy * np.sin(delta)
        Vyt = (-Vx * np.sin(delta)) + Vy * np.cos(delta)
        Sx = -(Vxt - Rd * omega_wheel) / (Rd * abs(omega_wheel) + Vn)
        Sy = -Vyt / (Rd * abs(omega_wheel) + Vn)
        Sg = -omega_wheel * np.sin(gamma) / (Rd * abs(omega_wheel) + Vn) * (L/2)
        Sb = -(L/3) * Wn / (Rd * abs(omega_wheel) + Vn)
        return Sx, Sy, Sg, Sb

    def get_tire_forces(self, Sx, Sy, Sg, Sb, Fz):
        # Normalize slip
        Sxn = Sx / self.hsxn
        Syn = Sy / self.hsyn
        Sn = np.sqrt(Sxn**2 + Syn**2)

        # Get tire offset
        L = 2 * np.sqrt(self.tire_radius * Fz / self.cz)
        N = L * self.__tire_offset__(Sy) 

        if Sn > 0: # Get slip angle sine and cosine
            cphi = Sxn / Sn
            sphi = Syn / Sn
        else:
            cphi = sphi = np.sqrt(2) / 2
        
        df0 = np.sqrt((self.dfx0 * self.hsxn * cphi)**2 + (self.dfy0 * self.hsyn * sphi)**2)
        Fm = np.sqrt((self.Fxm * cphi)**2 + (self.Fym * sphi)**2)
        Sm = np.sqrt((self.Sxm / self.hsxn * cphi)**2 + (self.Sym / self.hsyn * sphi)**2)
        Fs = np.sqrt((self.Fxs * cphi)**2 + (self.Fys * sphi)**2)
        Ss = np.sqrt((self.Sxs / self.hsxn * cphi)**2 + (self.Sys / self.hsyn * sphi)**2)
        F, dF = self.__tire_combined_forces__(Sn, df0, Fm, Sm, Fs, Ss)
        Fx = F * cphi
        Fy = F * sphi + df0 * Sg / 3
        Mz = -N * Fy + (L/3) * self.__tire_combined_forces__(Sb, df0, Fm, Sm, Fs, Ss)[0] # Rb = L/3
        return Fx, Fy, Mz, Sx, Sy

    def get_rolling_resistance(self, omega, Fz):
        Rd = self.lamb * self.tire_radius + (1 - self.lamb) * (self.tire_radius - Fz / self.cz)
        Ty = - Rd * omega / (abs(Rd * omega) + 0.001) * Fz * self.tire_radius * self.frr
        return Ty

    def all_info(self, Fz=3500, camber_slip=0, bore_slip=0, matplotlib=True):
        '''Prints information about the tires. Fz is in Newtons and camber_angle in radians'''
        print("Vertical force considered: {:.0F}N".format(Fz))
        print("Camber slip = {:.2f}".format(camber_slip))
        print("Bore slip = {:.2f}".format(bore_slip))
        # Longitudinal slip
        Sx = np.linspace(-0.5, 0.5, 1001)
        sy = [0, 0.11, 0.22, 0.33, 0.44, 0.55]
        forces = []
        for slipy in sy:
            Fx = [self.get_tire_forces(Sx[i], slipy, camber_slip, bore_slip, Fz)[0] for i in range(len(Sx))]
            forces.append(Function(Sx, Fx, 'Slip longitudinal (sx)', 'Força longitudinal (N)', name='Sy={:.2F}, Sg={:.2F}'.format(slipy, camber_slip)))
        forces[0].comparaNPlots(forces[1:], title='Slip Longitudinal TMEasy')

        # Longitudinal slip
        Sy = np.linspace(-0.5, 0.5, 1001)
        sx = [0, 0.088, 0.176, 0.264, 0.352, 0.44]
        forces = []
        for slipx in sx:
            Fy = [self.get_tire_forces(slipx, Sy[i], camber_slip, bore_slip, Fz)[1] for i in range(len(Sy))]
            forces.append(Function(Sy, Fy, 'Slip Lateral (sy)', 'Força Lateral (N)', name='Sx={:.2F}, Sg={:.2F}'.format(slipx, camber_slip)))
        print('Coenering stiffness = {:.1f}N/৹'.format(forces[0].derivate(0, np.deg2rad(1)) * np.deg2rad(1)))
        forces[0].comparaNPlots(forces[1:], title='Slip Lateral TMEasy')

        # Longitudinal slip
        Sy = np.linspace(-0.5, 0.5, 1001)
        sx = [0, 0.088, 0.176, 0.264, 0.352, 0.44]
        forces = []
        for slipx in sx:
            Mz = [self.get_tire_forces(slipx, Sy[i], camber_slip, bore_slip, Fz)[2] for i in range(len(Sy))]
            forces.append(Function(Sy, Mz, 'Slip Lateral (sy)', 'Momento restaurador (Nm)', name='Sx={:.2F}, Sg={:.2F}'.format(slipx, camber_slip)))
        forces[0].comparaNPlots(forces[1:], title='Momento restaurador TMEasy')
        
        # # Heatmap
        # x = y = np.linspace(-0.5, 0.5, 1001)
        # [SX, SY] = np.meshgrid(x, y)
        # FX = np.zeros((len(x), len(y)))
        # FY = np.zeros((len(x), len(y)))
        # MZ = np.zeros((len(x), len(y)))
        # for i in range(len(x)):
        #     for j in range(len(y)):
        #         FX[i, j], FY[i, j], MZ[i, j], temp1, temp2 = self.get_tire_forces(SX[i, j], SY[i, j], Fz)
        
        # fig_fx = go.Figure(data=[go.Surface(x=SX, y=SY, z=FX)])

        # fig_fx.update_layout(title='Longitudinal force', autosize=False,
        #                 width=500, height=500,
        #                 margin=dict(l=65, r=50, b=65, t=90))

        # fig_fx.show()

        # fig_fy = go.Figure(data=[go.Surface(x=SX, y=SY, z=FY)])

        # fig_fy.update_layout(title='Lateral force', autosize=False,
        #                 width=500, height=500,
        #                 margin=dict(l=65, r=50, b=65, t=90))

        # fig_fy.show()