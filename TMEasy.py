import numpy as np
from Function import Function

''' Here lies an inplementation of TMEasy tire model from George Rill and Abel
Castro. The implmentation is based on the Matlab code presented at Road Vehicle 
Dynamics, second edition.'''

class Tire:
    def __init__(self, radius, cz, dfx0, dfy0, Fxm, Fym, Sxm, Sym, Fxs, Fys, Sxs, Sys, Sy0, SyE, lamb, n2L0):
        self.tire_radius = radius
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

        # Calculate tire parameters
        self.hsxn = Sxm / (Sxm + Sym) + Fxm / (dfx0 * (Fxm / dfx0 + Fym / dfy0))
        self.hsyn = Sym / (Sxm + Sym) + Fym / (dfy0 * (Fxm / dfx0 + Fym / dfy0))
        
    
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

    def calculateSlip(self, Vx, Vy, delta, omega_wheel, Fz):
        # Calculate slip parameters
        Vn = 0.01
        Rd = self.lamb * self.tire_radius + (1 - self.lamb) * (self.tire_radius - Fz / cz)
        Vxt = Vx * np.cos(delta) + Vy * np.sin(delta)
        Vyt = (-Vx * np.sin(delta)) + Vy * np.cos(delta)
        Sx = -(Vxt - Rd * omega_wheel) / (Rd * abs(omega_wheel) + Vn)
        Sy = -(Vxt - Rd * omega_wheel) / (Rd * abs(omega_wheel) + Vn)
        return Sx, Sy

    def get_tire_forces(self, Sx, Sy, Fz):
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
        Fm = np.sqrt((self.Fxm * cphi)**2 + (self.Fym * sphi)**2);
        Sm = np.sqrt((self.Sxm / self.hsxn * cphi)**2 + (self.Sym / self.hsyn * sphi)**2);
        Fs = np.sqrt((self.Fxs * cphi)**2 + (self.Fys * sphi)**2);
        Ss = np.sqrt((self.Sxs / self.hsxn * cphi)**2 + (self.Sys / self.hsyn * sphi)**2);
        F, dF = self.__tire_combined_forces__(Sn, df0, Fm, Sm, Fs, Ss);
        Fx = F * cphi;
        Fy = F * sphi;
        Mz = -N * Fy;
        return Fx, Fy, Mz