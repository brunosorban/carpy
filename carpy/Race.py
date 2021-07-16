import numpy as np
from Function import Function
from scipy import integrate

class Race:
    def __init__(self, 
        Vehicle, 
        Driver, 
        initialSolution=None, 
        maxTime=60, 
        rtol=1e-6, 
        atol=1e-6,
        maxStep=1e-3
    ):

        self.Vehicle = Vehicle
        self.Driver = Driver
        self.maxTime = maxTime
        self.initialSolution = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0] if initialSolution is None else initialSolution
        self.solution = [[self.initialSolution[i]] for i in range(len(self.initialSolution))]
        self.atol = atol
        self.rtol = rtol
        self.maxStep = maxStep

        self.T1 = [0]
        self.T2 = [0]
        self.T3 = [0]
        self.T4 = [0]
        self.delta_1 = [0]
        self.delta_2 = [0]
        self.delta_3 = [0]
        self.delta_4 = [0]
        self.ax = [0]
        self.ay = [0]
        self.avx = [0]
        self.avy = [0]
        self.omega1= [0]
        self.omega2= [0]
        self.omega3 = [0]
        self.omega4 = [0]
        self.omegaDot1 = [0]
        self.omegaDot2 = [0]
        self.omegaDot3 = [0]
        self.omegaDot4 = [0]

    
    def uDot(self, t, sol):
        '''Calculate the state space derivatives'''
        # Retrieve integration data
        x, y, vx, vy, psi, psiDot, omega1, omega2, omega3, omega4 = sol

        # Get velocities in vehicle coordinate system
        vvx = vx * np.cos(psi) + vy * np.sin(psi)
        vvy = -vx * np.sin(psi) + vy * np.cos(psi)
        v = np.sqrt(vx**2 + vy**2)

        # Get control variables
        # T1 = self.Driver.T1_control(v)
        # T2 = self.Driver.T2_control(v)
        # T3 = self.Driver.T3_control(v)
        # T4 = self.Driver.T4_control(v)
        T1 = T2 = (60-v) * 200
        T3 = T4 = 0

        # delta_1 = self.Driver.delta_1_control(t)
        # delta_2 = self.Driver.delta_2_control(t)
        # delta_3 = self.Driver.delta_3_control(t)
        # delta_4 = self.Driver.delta_4_control(t)
        delta_1 = delta_2 = delta_3 = delta_4 = 0

        # Calculate tire forces
        ftx_1, fty_1, Mz_1 = self.Vehicle.Tire(vvx - psiDot * self.Vehicle.wf/2, vvy + psiDot * self.Vehicle.lf, delta_1, omega1, self.Vehicle.Fz / self.Vehicle.n_tires)
        ftx_2, fty_2, Mz_2 = self.Vehicle.Tire(vvx + psiDot * self.Vehicle.wf/2, vvy - psiDot * self.Vehicle.lf, delta_2, omega2, self.Vehicle.Fz / self.Vehicle.n_tires)
        ftx_3, fty_3, Mz_3 = self.Vehicle.Tire(vvx - psiDot * self.Vehicle.wr/2, vvy + psiDot * self.Vehicle.lr, delta_3, omega3, self.Vehicle.Fz / self.Vehicle.n_tires)
        ftx_4, fty_4, Mz_4 = self.Vehicle.Tire(vvx + psiDot * self.Vehicle.wr/2, vvy - psiDot * self.Vehicle.lr, delta_4, omega4, self.Vehicle.Fz / self.Vehicle.n_tires)

        # Convert tire forces to vehicle coordinate system
        fvx_1 = ftx_1 * np.cos(delta_1) - fty_1 * np.sin(delta_1) 
        fvy_1 = ftx_1 * np.sin(delta_1) + fty_1 * np.cos(delta_1)

        fvx_2 = ftx_2 * np.cos(delta_2) - fty_2 * np.sin(delta_2)
        fvy_2 = ftx_2 * np.sin(delta_2) + fty_2 * np.cos(delta_2)

        fvx_3 = ftx_3 * np.cos(delta_3) - fty_3 * np.sin(delta_3) 
        fvy_3 = ftx_3 * np.sin(delta_3) + fty_3 * np.cos(delta_3)

        fvx_4 = ftx_4 * np.cos(delta_4) - fty_4 * np.sin(delta_4)
        fvy_4 = ftx_4 * np.sin(delta_4) + fty_4 * np.cos(delta_4)

        # print(fvx_1, fvx_2, fvx_3, fvx_4)

        # Calculate other forces
        Fd = 0
        Frr = 0

        # Calculate accelerations in vehivle coordinate system
        avx = (fvx_1 + fvx_2 + fvx_3 + fvx_4 + Frr + Fd) / self.Vehicle.vehicle_equivalent_mass
        avy = (fvy_1 + fvy_2 + fvy_3 + fvy_4) / self.Vehicle.vehicle_equivalent_mass

        # Calculete state space accelerations
        ax = avx * np.cos(psi) - avy * np.sin(psi)
        ay = avx * np.sin(psi) + avy * np.cos(psi)
        tau_z = (fvy_1 + fvy_2) * self.Vehicle.lf - (fvy_3 + fvy_4) * self.Vehicle.lr + (fvx_2 - fvx_1) * self.Vehicle.wf / 2 + (fvx_4 - fvx_3) * self.Vehicle.wr / 2
        psiDotDot = tau_z / self.Vehicle.Izz

        # Calculate accelerations in the wheels
        omegaDot1 = (T1 - ftx_1 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz
        omegaDot2 = (T2 - ftx_2 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz
        omegaDot3 = (T3 - ftx_3 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz
        omegaDot4 = (T4 - ftx_4 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz

        self.ax.append(ax)
        self.ay.append(ay)
        self.avx.append(avx)
        self.avy.append(avy)
        self.T1.append(T1)
        self.T2.append(T2)
        self.T3.append(T3)
        self.T4.append(T4)
        self.delta_1.append(delta_1)
        self.delta_2.append(delta_2)
        self.delta_3.append(delta_3)
        self.delta_4.append(delta_4)
        self.omega1.append(omega1)
        self.omega2.append(omega2)
        self.omega3.append(omega3)
        self.omega4.append(omega4)
        self.omegaDot1.append(omegaDot1)
        self.omegaDot2.append(omegaDot2)
        self.omegaDot3.append(omegaDot3)
        self.omegaDot4.append(omegaDot4)
        
        uDot = [
            vx,
            vy,
            ax,
            ay,
            psiDot,
            psiDotDot,
            omegaDot1,
            omegaDot2,
            omegaDot3,
            omegaDot4
        ]
        return uDot        

    def solver(self):
        time = [0]
        solver = integrate.LSODA(self.uDot, t0=0, y0=self.initialSolution, max_step=self.maxStep, t_bound=self.maxTime, rtol=self.rtol, atol=self.atol)
        while solver.status != 'finished':
            solver.step()
            time.append(solver.t)
            for i in range(len(self.solution)):
                self.solution[i].append(solver.y[i])

        self.XY = Function(self.solution[0], self.solution[1])
        self.X = Function(time, self.solution[0])
        self.Y = Function(time, self.solution[1])
        self.R = Function(np.sqrt(np.array( self.solution[0])**2 + np.array( self.solution[1])**2))
        self.Vx = Function(time, self.solution[2])
        self.Vy = Function(time, self.solution[3])
        self.V = Function(time, np.sqrt(np.array( self.solution[2])**2 + np.array( self.solution[3])**2))
        self.psi = Function(time, self.solution[4])
        self.psiDot = Function(time, self.solution[5])
        self.omega = Function(time, self.solution[6])
        self.omegaDot = Function(time, self.solution[7])
        self.T1 = Function(time, self.T1)
        self.T2 = Function(time, self.T2)
        self.T3 = Function(time, self.T3)
        self.T4 = Function(time, self.T4)
        self.delta_1 = Function(time, 180 / np.pi * np.array(self.delta_1))
        self.delta_2 = Function(time, 180 / np.pi * np.array(self.delta_2))
        self.delta_3 = Function(time, 180 / np.pi * np.array(self.delta_3))
        self.delta_4 = Function(time, 180 / np.pi * np.array(self.delta_4))
        self.ax = Function(time, self.ax)
        self.ay = Function(time, self.ay)
        self.avx = Function(time, self.avx)
        self.avy = Function(time, self.avy)
        self.omega1= Function(time, self.omega1)
        self.omega2= Function(time, self.omega2)
        self.omega3 = Function(time, self.omega3)
        self.omega4 = Function(time, self.omega4)
        self.omegaDot1 = Function(time, self.omegaDot1)
        self.omegaDot2 = Function(time, self.omegaDot2)
        self.omegaDot3 = Function(time, self.omegaDot3)
        self.omegaDot4 = Function(time, self.omegaDot4)
