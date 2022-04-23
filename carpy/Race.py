import numpy as np
from Function import Function
from scipy import integrate
import matplotlib.pyplot as plt

# For the animation class
import sys, pygame
import imageio # Gif da animação

class Race:
    def __init__(self, 
        Vehicle, 
        Driver,
        Track, 
        initialSolution=None, 
        maxTime=60, 
        rtol=1e-6, 
        atol=1e-6,
        maxStep=1e-3
    ):
        self.Vehicle = Vehicle
        self.Driver = Driver
        self.Track = Track
        self.maxTime = maxTime
        self.initialSolution = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] if initialSolution is None else initialSolution
        self.t0 = 0
        self.atol = atol
        self.rtol = rtol
        self.maxStep = maxStep
        self.time_count = 0
        self.solver()

    def rotation_matrix(self, phi, theta, psi):
        c11 = np.cos(theta) * np.cos(psi)
        c12 = np.cos(theta) * np.sin(psi)
        c13 = -np.sin(theta)
        c21 = np.sin(phi) * np.sin(theta) * np.cos(psi) - np.cos(phi) * np.sin(psi)
        c22 = np.sin(phi) * np.sin(theta) * np.sin(psi) + np.cos(phi) * np.cos(psi)
        c23 = np.sin(phi) * np.cos(theta)
        c31 = np.cos(phi) * np.sin(theta) * np.cos(psi) + np.sin(phi) * np.sin(psi)
        c32 = np.cos(phi) * np.sin(theta) * np.sin(psi) - np.sin(phi) * np.cos(psi)
        c33 = np.cos(phi) * np.cos(theta)
        return np.array([[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]])
    
    def uDot(self, t, sol, post_process=False):
        '''Calculate the state space derivatives'''
        # Retrieve integration data
        x, y, z, vx, vy, vz, phi, theta, psi, phi_dot, theta_dot, psi_dot, z1, z2, z3, z4, vz1, vz2, vz3, vz4, omega1, omega2, omega3, omega4 = sol

        # Mount rotation matrix
        C = self.rotation_matrix(phi, theta, psi)

        # Get velocities in vehicle coordinate system
        [vvx, vvy, vvz] = C @ [vx, vy, vz]
        g = C @ np.array([0, 0, -9.81])

        # Get track variables; A base rotation is applied to consider distance from CG to each wheel
        Ct = self.rotation_matrix(0, 0, psi)
        [xc1, yc1, zc1] = np.transpose(Ct) @ [self.Vehicle.lf, self.Vehicle.wf, 0] + [x, y, 0]
        [xc2, yc2, zc2] = np.transpose(Ct) @ [self.Vehicle.lf, -self.Vehicle.wf, 0] + [x, y, 0]
        [xc3, yc3, zc2] = np.transpose(Ct) @ [-self.Vehicle.lr, self.Vehicle.wr, 0] + [x, y, 0]
        [xc4, yc4, zc2] = np.transpose(Ct) @ [-self.Vehicle.lr, -self.Vehicle.wr, 0] + [x, y, 0]
        
        zc1 = self.Track.get_track_height(xc1, yc1)
        zc2 = self.Track.get_track_height(xc2, yc2)
        zc3 = self.Track.get_track_height(xc3, yc3)
        zc4 = self.Track.get_track_height(xc4, yc4)

        vzc1 = self.Track.get_track_speed(xc1, yc1, vx, vy)
        vzc2 = self.Track.get_track_speed(xc2, yc2, vx, vy)
        vzc3 = self.Track.get_track_speed(xc3, yc3, vx, vy)
        vzc4 = self.Track.get_track_speed(xc4, yc4, vx, vy)
        
        # Get control variables
        T1, T2, T3, T4 = self.Driver.get_torque(vvx, t)
        delta_sw = self.Driver.get_steering(t) # Get driver steering wheel angle
        delta_1, delta_2, delta_3, delta_4 = self.Vehicle.get_steer_wbw(delta_sw) # Get steering angle wheel-by-wheel
        Wn_1, Wn_2, Wn_3, Wn_4 = 0, 0, 0, 0 # Wheel velocity
        # gamma1, gamma2, gamma3, gamma4 = np.deg2rad(0), np.deg2rad(0), 0, 0 # Camber angle

        # Get tire weight distribution
        ftz_1, ftz_2, ftz_3, ftz_4, fs_1, fs_2, fs_3, fs_4 = self.Vehicle.get_vertical_load(z, vvz, phi, theta, phi_dot, theta_dot, z1, z2, z3, z4, vz1, vz2, vz3, vz4, zc1, zc2, zc3, zc4, vzc1, vzc2, vzc3, vzc4)
        if min(ftz_1, ftz_2, ftz_3, ftz_4, fs_1, fs_2, fs_3, fs_4) < 0: 
            self.solver.status = 'finished'
            return [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        # Calculate tire forces
        ftx_1, fty_1, Mz_1, sx_1, sy_1 = self.Vehicle.Tire(vvx - psi_dot * self.Vehicle.wf/2, vvy + psi_dot * self.Vehicle.lf, delta_1, self.Vehicle.gamma1, omega1, Wn_1, ftz_1)
        ftx_2, fty_2, Mz_2, sx_2, sy_2 = self.Vehicle.Tire(vvx + psi_dot * self.Vehicle.wf/2, vvy + psi_dot * self.Vehicle.lf, delta_2, self.Vehicle.gamma2, omega2, Wn_2, ftz_2)
        ftx_3, fty_3, Mz_3, sx_3, sy_3 = self.Vehicle.Tire(vvx - psi_dot * self.Vehicle.wr/2, vvy - psi_dot * self.Vehicle.lr, delta_3, self.Vehicle.gamma3, omega3, Wn_3, ftz_3)
        ftx_4, fty_4, Mz_4, sx_4, sy_4 = self.Vehicle.Tire(vvx + psi_dot * self.Vehicle.wr/2, vvy - psi_dot * self.Vehicle.lr, delta_4, self.Vehicle.gamma4, omega4, Wn_4, ftz_4)

        # Get rolling resistance torques
        Trr1 = self.Vehicle.Tire.get_rolling_resistance(omega1, ftz_1)
        Trr2 = self.Vehicle.Tire.get_rolling_resistance(omega2, ftz_2)
        Trr3 = self.Vehicle.Tire.get_rolling_resistance(omega3, ftz_3)
        Trr4 = self.Vehicle.Tire.get_rolling_resistance(omega4, ftz_4)

        # Convert tire forces to vehicle coordinate system
        [fvx_1, fvy_1, fvz_1] = np.transpose(self.rotation_matrix(phi, theta, delta_1)) @ np.array([ftx_1, fty_1, fs_1])
        [fvx_2, fvy_2, fvz_2] = np.transpose(self.rotation_matrix(phi, theta, delta_2)) @ np.array([ftx_2, fty_2, fs_2])
        [fvx_3, fvy_3, fvz_3] = np.transpose(self.rotation_matrix(phi, theta, delta_3)) @ np.array([ftx_3, fty_3, fs_3])
        [fvx_4, fvy_4, fvz_4] = np.transpose(self.rotation_matrix(phi, theta, delta_4)) @ np.array([ftx_4, fty_4, fs_4])

        # Convert tire moments to vehicle coordinate system
        [mvx_1, mvy_1, mvz_1] = np.transpose(self.rotation_matrix(phi, theta, delta_1)) @ np.array([0, 0, Mz_1])
        [mvx_2, mvy_2, mvz_2] = np.transpose(self.rotation_matrix(phi, theta, delta_2)) @ np.array([0, 0, Mz_2])
        [mvx_3, mvy_3, mvz_3] = np.transpose(self.rotation_matrix(phi, theta, delta_3)) @ np.array([0, 0, Mz_3])
        [mvx_4, mvy_4, mvz_4] = np.transpose(self.rotation_matrix(phi, theta, delta_4)) @ np.array([0, 0, Mz_4])

        # Calculate drag force
        Fd = 1 / 2 * self.Vehicle.rho * self.Vehicle.cd * vvx**2 * self.Vehicle.af * np.sign(vvx)

        # Calculate accelerations in vehivle coordinate system
        avx = (fvx_1 + fvx_2 + fvx_3 + fvx_4 - Fd) / self.Vehicle.vehicle_equivalent_mass + g[0]
        avy = (fvy_1 + fvy_2 + fvy_3 + fvy_4) / self.Vehicle.vehicle_mass + g[1]
        avz = (fvz_1 + fvz_2 + fvz_3 + fvz_4) / self.Vehicle.vehicle_mass + g[2]

        # Calculete state space accelerations
        [ax, ay, az] = np.transpose(C) @ [avx, avy, avz]
        
        # Calculate tires vertical acceleration
        az1 = (ftz_1 - fs_1) / self.Vehicle.Tire.mass
        az2 = (ftz_2 - fs_2) / self.Vehicle.Tire.mass
        az3 = (ftz_3 - fs_3) / self.Vehicle.Tire.mass
        az4 = (ftz_4 - fs_4) / self.Vehicle.Tire.mass

        # phi angular acceleration
        tau_x = (fvz_1 - fvz_2) * self.Vehicle.wf / 2 + (fvx_3 - fvx_4) * self.Vehicle.wr / 2 + (fvy_1 + fvy_2 + fvy_3 + fvy_4) * self.Vehicle.CG_height + mvx_1 + mvx_2 + mvx_3 + mvx_4
        phi_dot_dot = (tau_x + (self.Vehicle.Iyy - self.Vehicle.Izz) * (psi_dot * theta_dot)) / self.Vehicle.Ixx 
        
        # theta angular acceleration
        tau_y = -(fvz_1 + fvz_2) * self.Vehicle.lf + (fvz_3 + fvz_4) * self.Vehicle.lr - (fvx_1 + fvx_2 + fvx_3 + fvx_4) * self.Vehicle.CG_height + mvy_1 + mvy_2 + mvy_3 + mvy_4
        theta_dot_dot = (tau_y + (self.Vehicle.Izz - self.Vehicle.Ixx) * (psi_dot * phi_dot)) / self.Vehicle.Iyy

        # psi angular acceleraion
        tau_z = (fvy_1 + fvy_2) * self.Vehicle.lf - (fvy_3 + fvy_4) * self.Vehicle.lr + (fvx_2 - fvx_1) * self.Vehicle.wf / 2 + (fvx_4 - fvx_3) * self.Vehicle.wr / 2 + mvz_1 + mvz_2 + mvz_3 + mvz_4
        psi_dot_dot = (tau_z + (self.Vehicle.Ixx - self.Vehicle.Iyy) * (phi_dot * theta_dot)) / self.Vehicle.Izz
        

        # Calculate accelerations in the wheels
        # Deveria aqui ser considerado o raio dinâmico?
        omega1_dot = (T1 - ftx_1 * self.Vehicle.Tire.tire_radius + Trr1) / self.Vehicle.Tire.Jzz
        omega2_dot = (T2 - ftx_2 * self.Vehicle.Tire.tire_radius + Trr2) / self.Vehicle.Tire.Jzz
        omega3_dot = (T3 - ftx_3 * self.Vehicle.Tire.tire_radius + Trr3) / self.Vehicle.Tire.Jzz
        omega4_dot = (T4 - ftx_4 * self.Vehicle.Tire.tire_radius + Trr4) / self.Vehicle.Tire.Jzz

        if not post_process:
            uDot = [
                vx,
                vy,
                vz,
                ax,
                ay,
                az,
                phi_dot,
                theta_dot,
                psi_dot,
                phi_dot_dot,
                theta_dot_dot,
                psi_dot_dot,
                vz1,
                vz2,
                vz3,
                vz4,
                az1,
                az2,
                az3,
                az4,
                omega1_dot,
                omega2_dot,
                omega3_dot,
                omega4_dot
            ]
            return uDot   

        else:
            # If post_process, save the data
            self.T1.append(T1)
            self.T2.append(T2)
            self.T3.append(T3)
            self.T4.append(T4)
            self.delta_sw.append(delta_sw)
            self.delta_1.append(delta_1)
            self.delta_2.append(delta_2)
            self.delta_3.append(delta_3)
            self.delta_4.append(delta_4)
            self.ftx_1.append(ftx_1)
            self.fty_1.append(fty_1)
            self.Mz_1.append(Mz_1)
            self.sx_1.append(sx_1)
            self.sy_1.append(sy_1)
            self.ftx_2.append(ftx_2)
            self.fty_2.append(fty_2)
            self.sx_2.append(sx_2)
            self.sy_2.append(sy_2)
            self.Mz_2.append(Mz_2)
            self.ftx_3.append(ftx_3)
            self.fty_3.append(fty_3)
            self.sx_3.append(sx_3)
            self.sy_3.append(sy_3)
            self.Mz_3.append(Mz_3)
            self.ftx_4.append(ftx_4)
            self.fty_4.append(fty_4)
            self.Mz_4.append(Mz_4)
            self.sx_4.append(sx_4)
            self.sy_4.append(sy_4)
            self.Trr.append(Trr1 + Trr2 + Trr3 + Trr4)
            self.Fd.append(Fd)
            self.ftz_1.append(ftz_1)
            self.ftz_2.append(ftz_2)
            self.ftz_3.append(ftz_3)
            self.ftz_4.append(ftz_4)
            self.fs_1.append(fs_1)
            self.fs_2.append(fs_2)
            self.fs_3.append(fs_3)
            self.fs_4.append(fs_4)

            self.x.append(x)
            self.y.append(y)
            self.z.append(z)
            self.vx.append(vx)
            self.vz.append(vz)
            self.vy.append(vy)
            self.vvx.append(vvx)
            self.vvy.append(vvy)
            self.vvz.append(vvz)
            self.ax.append(ax)
            self.ay.append(ay)
            self.az.append(az)
            self.avx.append(avx)
            self.avy.append(avy)
            self.avz.append(avz)
            self.phi.append(phi)
            self.phi_dot.append(phi_dot)
            self.phi_dot_dot.append(phi_dot_dot)
            self.theta.append(theta)
            self.theta_dot.append(theta_dot)
            self.theta_dot_dot.append(theta_dot_dot)
            self.psi.append(psi)
            self.psi_dot.append(psi_dot)
            self.psi_dot_dot.append(psi_dot_dot)
            self.z1.append(z1)
            self.z2.append(z2)
            self.z3.append(z3)
            self.z4.append(z4)
            self.vz1.append(vz1)
            self.vz2.append(vz2)
            self.vz3.append(vz3)
            self.vz4.append(vz4)
            self.omega1.append(omega1)
            self.omega2.append(omega2)
            self.omega3.append(omega3)
            self.omega4.append(omega4)
            self.omega1_dot.append(omega1_dot)
            self.omega2_dot.append(omega2_dot)
            self.omega3_dot.append(omega3_dot)
            self.omega4_dot.append(omega4_dot)

    def solver(self):
        # Initialize the data
        self.time = [self.t0]
        self.solution = [self.initialSolution]
        self.last_avy = 0

        # Create the solver
        self.solver = integrate.LSODA(self.uDot, t0=0, y0=self.initialSolution, t_bound=self.maxTime, max_step=self.maxStep, rtol=self.rtol, atol=self.atol)
        
        # Iterate until max_time is reached
        while self.solver.status != 'finished':
            self.solver.step()
            self.time.append(self.solver.t)
            self.solution.append(self.solver.y)
        print("Solution Finished")

    def post_process(self):
        # Initialize the data
        self.T1 = [0]
        self.T2 = [0]
        self.T3 = [0]
        self.T4 = [0]
        self.delta_sw = [0]
        self.delta_1 = [0]
        self.delta_2 = [0]
        self.delta_3 = [0]
        self.delta_4 = [0]
        self.ftx_1 = [0]
        self.fty_1 = [0]
        self.Mz_1 = [0]
        self.sx_1 = [0]
        self.sy_1 = [0]
        self.ftx_2 = [0]
        self.fty_2 = [0]
        self.Mz_2 = [0]
        self.sx_2 = [0]
        self.sy_2 = [0]
        self.ftx_3 = [0]
        self.fty_3 = [0]
        self.sx_3 = [0]
        self.sy_3 = [0]
        self.Mz_3 = [0]
        self.ftx_4 = [0]
        self.fty_4 = [0]
        self.Mz_4 = [0]
        self.sx_4 = [0]
        self.sy_4 = [0]
        self.Trr = [0]
        self.Fd = [0]
        self.ftz_1 = [0]
        self.ftz_2 = [0]
        self.ftz_3 = [0]
        self.ftz_4 = [0]
        self.fs_1 = [0]
        self.fs_2 = [0]
        self.fs_3 = [0]
        self.fs_4 = [0]

        self.x = [0]
        self.y = [0]
        self.z = [0]
        self.vx = [0]
        self.vy = [0]
        self.vz = [0]
        self.vvx = [0]
        self.vvy = [0]
        self.vvz = [0]
        self.ax = [0]
        self.ay = [0]
        self.az = [0]
        self.avx = [0]
        self.avy = [0]
        self.avz = [0]
        self.phi = [0]
        self.phi_dot = [0]
        self.phi_dot_dot = [0]
        self.theta = [0]
        self.theta_dot = [0]
        self.theta_dot_dot = [0]
        self.psi = [0]
        self.psi_dot = [0]
        self.psi_dot_dot = [0]
        self.z1 = [0]
        self.z2 = [0]
        self.z3 = [0]
        self.z4 = [0]
        self.vz1 = [0]
        self.vz2 = [0]
        self.vz3 = [0]
        self.vz4 = [0]
        self.omega1 = [0]
        self.omega2 = [0]
        self.omega3 = [0]
        self.omega4 = [0]
        self.omega1_dot = [0]
        self.omega2_dot = [0]
        self.omega3_dot = [0]
        self.omega4_dot = [0]

        # Iterate with udot function to generate data from solution vector
        for i in range(1, len(self.solution)):
            self.uDot(self.time[i], self.solution[i], post_process=True)

        # Create function objects from the retrieved data
        self.T1 = Function(self.time, self.T1, xlabel='Time (s)', ylabel='Front Left tire torque (Nm)', name='T1')
        self.T2 = Function(self.time, self.T2, xlabel='Time (s)', ylabel='Front Right tire torque (Nm)', name='T2')
        self.T3 = Function(self.time, self.T3, xlabel='Time (s)', ylabel='Rear Left tire torque (Nm)', name='T3')
        self.T4 = Function(self.time, self.T4, xlabel='Time (s)', ylabel='Rear Right tire torque (Nm)', name='T4')
        self.delta_sw = Function(self.time, 180 / np.pi * np.array(self.delta_sw), xlabel='Time (s)', ylabel='Front Left tire steering angle (degrees)', name='δsw')
        self.delta_1 = Function(self.time, 180 / np.pi * np.array(self.delta_1), xlabel='Time (s)', ylabel='Front Left tire steering angle (degrees)', name='δ1')
        self.delta_2 = Function(self.time, 180 / np.pi * np.array(self.delta_2), xlabel='Time (s)', ylabel='Front Right tire steering angle (degrees)', name='δ2')
        self.delta_3 = Function(self.time, 180 / np.pi * np.array(self.delta_3), xlabel='Time (s)', ylabel='Rear Left tire steering angle (degrees)', name='δ3')
        self.delta_4 = Function(self.time, 180 / np.pi * np.array(self.delta_4), xlabel='Time (s)', ylabel='Rear Right tire steering angle (degrees)', name='δ4')
        self.ftx_1 = Function(self.time, self.ftx_1, xlabel='Time (s)', ylabel="Front Left tire force in tire's x axis (N)", name='ftx_1')
        self.fty_1 = Function(self.time, self.fty_1, xlabel='Time (s)', ylabel="Front Left tire force in tire's y axis (N)", name='fty_1')
        self.Mz_1 = Function(self.time, self.Mz_1, xlabel='Time (s)', ylabel='Front Left tire moment in tire z axis (Nm)', name='Mz_1')
        self.sx_1 = Function(self.time, self.sx_1, xlabel='Time (s)', ylabel='Longitudinal slip', name='sx_1')
        self.sy_1 = Function(self.time, self.sy_1, xlabel='Time (s)', ylabel='Lateral slip', name='sy_1')
        self.ftx_2 = Function(self.time, self.ftx_2, xlabel='Time (s)', ylabel='Front Right tire force in tire x axis (N)', name='ftx_2')
        self.fty_2 = Function(self.time, self.fty_2, xlabel='Time (s)', ylabel='Front Right tire force in tire y axis (N)', name='fty_2')
        self.Mz_2 = Function(self.time, self.Mz_2, xlabel='Time (s)', ylabel='Front Right tire moment in tire z axis (Nm)', name='Mz_2')
        self.sx_2 = Function(self.time, self.sx_2, xlabel='Time (s)', ylabel='Longitudinal slip', name='sx_2')
        self.sy_2 = Function(self.time, self.sy_2, xlabel='Time (s)', ylabel='Lateral slip', name='sy_2')
        self.ftx_3 = Function(self.time, self.ftx_3, xlabel='Time (s)', ylabel='Rear Left tire force in tire x axis (N)', name='ftx_3')
        self.fty_3 = Function(self.time, self.fty_3, xlabel='Time (s)', ylabel='Rear Left tire force in tire y axis (N)', name='fty_3')
        self.Mz_3 = Function(self.time, self.Mz_3, xlabel='Time (s)', ylabel='Rear Left tire moment in tire z axis (Nm)', name='Mz_3')
        self.sx_3 = Function(self.time, self.sx_3, xlabel='Time (s)', ylabel='Longitudinal slip', name='sx_3')
        self.sy_3 = Function(self.time, self.sy_3, xlabel='Time (s)', ylabel='Lateral slip', name='sy_3')
        self.ftx_4 = Function(self.time, self.ftx_4, xlabel='Time (s)', ylabel='Rear Right tire force in tire x axis (N)', name='ftx_4')
        self.fty_4 = Function(self.time, self.fty_4, xlabel='Time (s)', ylabel='Rear Right tire force in tire y axis (N)', name='fty_4')
        self.Mz_4 = Function(self.time, self.Mz_4, xlabel='Time (s)', ylabel='Rear Right tire moment in tire z axis (Nm)', name='Mz_4')
        self.sx_4 = Function(self.time, self.sx_4, xlabel='Time (s)', ylabel='Longitudinal slip', name='sx_4')
        self.sy_4 = Function(self.time, self.sy_4, xlabel='Time (s)', ylabel='Lateral slip', name='sy_4')
        self.Trr = Function(self.time, self.Trr, xlabel='Time (s)', ylabel='Total roling resistance (Nm)', name='Trr')
        self.Fd = Function(self.time, self.Fd, xlabel='Time (s)', ylabel='Drag force (N)', name='Fd')
        self.ftz_1 = Function(self.time, self.ftz_1, xlabel='Time (s)', ylabel='Front left tire vertical force (N)', name='ftz_1')
        self.ftz_2 = Function(self.time, self.ftz_2, xlabel='Time (s)', ylabel='Front right tire vertical force (N)', name='ftz_2')
        self.ftz_3 = Function(self.time, self.ftz_3, xlabel='Time (s)', ylabel='Rear left tire vertical force (N)', name='ftz_3')
        self.ftz_4 = Function(self.time, self.ftz_4, xlabel='Time (s)', ylabel='Rear right tire vertical force (N)', name='ftz_4')
        self.fs_1 = Function(self.time, self.fs_1, xlabel='Time (s)', ylabel='Front left suspension vertical force (N)', name='fs_1')
        self.fs_2 = Function(self.time, self.fs_2, xlabel='Time (s)', ylabel='Front right suspension vertical force (N)', name='fs_2')
        self.fs_3 = Function(self.time, self.fs_3, xlabel='Time (s)', ylabel='Rear left suspension vertical force (N)', name='fs_3')
        self.fs_4 = Function(self.time, self.fs_4, xlabel='Time (s)', ylabel='Rear right suspension vertical force (N)', name='fs_4')

        self.r = Function(self.time, np.sqrt(np.array(self.x)**2 + np.array(self.y)**2), xlabel='Time (s)', ylabel='Distance from the origin (m)', name='r')
        self.v = Function(self.time, np.sqrt(np.array(self.vx)**2 + np.array(self.vy)**2), xlabel='Time (s)', ylabel='Velocity (m/s)', name='v')
        self.xy = Function(self.x, self.y, xlabel='Position in inertial frame x axis (m)', ylabel='Position in inertial frame y axis (m)', name='xy')
        self.x = Function(self.time, self.x, xlabel='Time (s)', ylabel='Position in inertial frame x axis (m)', name='x')
        self.y = Function(self.time, self.y, xlabel='Time (s)', ylabel='Position in inertial frame y axis (m)', name='y')
        self.z = Function(self.time, self.z, xlabel='Time (s)', ylabel='Position in inertial frame z axis (m)', name='z')
        self.vx = Function(self.time, self.vx, xlabel='Time (s)', ylabel='Velocity in inertial frame x axis (m/s)', name='vx')
        self.vy = Function(self.time, self.vy, xlabel='Time (s)', ylabel='Velocity in inertial frame y axis (m/s)', name='vy')
        self.vz = Function(self.time, self.vz, xlabel='Time (s)', ylabel='Velocity in inertial frame z axis (m/s)', name='vz')
        self.vvx = Function(self.time, self.vvx, xlabel='Time (s)', ylabel='Velocity in vehicle frame x axis (m/s)', name='vvx')
        self.vvy = Function(self.time, self.vvy, xlabel='Time (s)', ylabel='Velocity in vehicle frame y axis (m/s)', name='vvy')
        self.vvz = Function(self.time, self.vvz, xlabel='Time (s)', ylabel='Velocity in vehicle frame z axis (m/s)', name='vvz')
        self.ax = Function(self.time, self.ax, xlabel='Time (s)', ylabel='Acceleration in inertial frame x axis (m/s²)', name='ax')
        self.ay = Function(self.time, self.ay, xlabel='Time (s)', ylabel='Acceleration in inertial frame y axis (m/s²)', name='ay')
        self.az = Function(self.time, self.az, xlabel='Time (s)', ylabel='Acceleration in inertial frame z axis (m/s²)', name='az')
        self.avx = Function(self.time, self.avx, xlabel='Time (s)', ylabel='Acceleration in vehicle frame x axis (m/s²)', name='avx')
        self.avy = Function(self.time, self.avy, xlabel='Time (s)', ylabel='Acceleration in vehicle frame y axis (m/s²)', name='avy')
        self.avz = Function(self.time, self.avz, xlabel='Time (s)', ylabel='Acceleration in vehicle frame z axis (m/s²)', name='avz')
        self.phi = Function(self.time,  180 / np.pi * np.array(self.phi), xlabel='Time (s)', ylabel='Roll angle (degrees)', name='phi')
        self.phi_dot = Function(self.time, self.phi_dot, xlabel='Time (s)', ylabel='Roll angular velocity (rad/s)', name='phi_dot')
        self.phi_dot_dot = Function(self.time, self.phi_dot_dot, xlabel='Time (s)', ylabel='Roll angular acceleration (rad/s²)', name='phi_dot_dot')
        self.theta = Function(self.time,  180 / np.pi * np.array(self.theta), xlabel='Time (s)', ylabel='Pitch angle (degrees)', name='theta')
        self.theta_dot = Function(self.time, self.theta_dot, xlabel='Time (s)', ylabel='Pitch angular velocity (rad/s)', name='theta_dot')
        self.theta_dot_dot = Function(self.time, self.theta_dot_dot, xlabel='Time (s)', ylabel='Pitch angular acceleration (rad/s²)', name='theta_dot_dot')
        self.psi = Function(self.time,  180 / np.pi * np.array(self.psi), xlabel='Time (s)', ylabel='Yaw angle (degrees)', name='psi')
        self.psi_dot = Function(self.time, self.psi_dot, xlabel='Time (s)', ylabel='Yaw angular velocity (rad/s)', name='psi_dot')
        self.psi_dot_dot = Function(self.time, self.psi_dot_dot, xlabel='Time (s)', ylabel='Yaw angular acceleration (rad/s²)', name='psi_dot_dot')
        self.z1 = Function(self.time, self.z1, xlabel='Time (s)', ylabel='Tire 1 position in inertial frame z axis (m)', name='z1')
        self.z2 = Function(self.time, self.z2, xlabel='Time (s)', ylabel='Tire 2 position in inertial frame z axis (m)', name='z2')
        self.z3 = Function(self.time, self.z3, xlabel='Time (s)', ylabel='Tire 3 position in inertial frame z axis (m)', name='z3')
        self.z4 = Function(self.time, self.z4, xlabel='Time (s)', ylabel='Tire 4 position in inertial frame z axis (m)', name='z4')
        self.vz1 = Function(self.time, self.vz1, xlabel='Time (s)', ylabel='Tire 1 velocity in inertial frame z axis (m/s)', name='vz1')
        self.vz2 = Function(self.time, self.vz2, xlabel='Time (s)', ylabel='Tire 2 velocity in inertial frame z axis (m/s)', name='vz2')
        self.vz3 = Function(self.time, self.vz3, xlabel='Time (s)', ylabel='Tire 3 velocity in inertial frame z axis (m/s)', name='vz3')
        self.vz4 = Function(self.time, self.vz4, xlabel='Time (s)', ylabel='Tire 4 velocity in inertial frame z axis (m/s)', name='vz4')
        self.omega1= Function(self.time, self.omega1, xlabel='Time (s)', ylabel='Front Left tire angular velocity (rad/s)', name='ω1')
        self.omega2= Function(self.time, self.omega2, xlabel='Time (s)', ylabel='Front Right tire angular velocity (rad/s)', name='ω2')
        self.omega3 = Function(self.time, self.omega3, xlabel='Time (s)', ylabel='Rear Left tire angular velocity (rad/s)', name='ω3')
        self.omega4 = Function(self.time, self.omega4, xlabel='Time (s)', ylabel='Rear Right tire angular velocity (rad/s)', name='ω4')
        self.omega1_dot = Function(self.time, self.omega1_dot, xlabel='Time (s)', ylabel='Front Left tire angular acceleration (rad/s²)', name='ω1_dot')
        self.omega2_dot = Function(self.time, self.omega2_dot, xlabel='Time (s)', ylabel='Front Right tire angular acceleration (rad/s²)', name='ω2_dot')
        self.omega3_dot = Function(self.time, self.omega3_dot, xlabel='Time (s)', ylabel='Rear Left tire angular acceleration (rad/s²)', name='ω3_dot')
        self.omega4_dot = Function(self.time, self.omega4_dot, xlabel='Time (s)', ylabel='Rear Right tire angular acceleration (rad/s²)', name='ω4_dot')

    # Animation 
    def blitRotateCenter(self, surf, image, topleft, angle):
        # Auxiliar function to blit the animation image rotated
        rotated_image = pygame.transform.rotate(image, angle)
        new_rect = rotated_image.get_rect(center = image.get_rect(topleft = topleft).center)
        surf.blit(rotated_image, new_rect)
        
    def mapi(self, x, y, xlim, ylim, xlim_inf, ylim_inf, size, realScale=True):
        # Auxiliar function to map the trajectory to fit in view
        x_factor = size[0] / xlim
        y_factor = size[1] / ylim
        if realScale:
            if x_factor > y_factor:
                x_factor = y_factor
            else:
                y_factor = x_factor
        return np.array([(x + abs(xlim_inf)) * x_factor, (y - abs(ylim_inf)) * y_factor])

    def animate(self, timeMax, save=False):
        # Initialization
        xlim_inf = abs(min(self.x.__Y_source__))
        xlim = max(self.x.__Y_source__) + 1e-3 + xlim_inf
        ylim_inf = abs(min(self.y.__Y_source__))
        ylim = max(self.y.__Y_source__) + 1e-3 + ylim_inf
        pygame.init()
        pygame.display.init()
        font = pygame.font.SysFont('Helvetica', 32)
        if save: w = imageio.get_writer('Videos/carpy.mp4', format='FFMPEG', fps=60)

        # Defining auxiliar colors 
        # white = 0, 0, 0
        # line_colour = 0, 0, 0
        
        # Creating animation screen
        size = np.array([1920, 1080])
        screen = pygame.display.set_mode(size)
        
        # Preparing images
        background = pygame.image.load('../Animation/TrackT01.png')
        background = (pygame.transform.scale(background, (size)))
        
        # background = screen.fill([255, 255, 255])

        car = pygame.image.load('../Animation/Car.png')
        car = pygame.transform.scale(car, (128, 128))

        car_side = pygame.image.load('../Animation/Car side.png')
        car_side = pygame.transform.scale(car_side, (int(65 * car_side.get_width() / car_side.get_height()), 65))

        car_front = pygame.image.load('../Animation/Car front.png')
        car_front = pygame.transform.scale(car_front, (128, 128))

        # Initialize position vectors
        initial_position = np.array([0, size[1]-128])
        position = initial_position
        position_side = np.array([0, 128+52])
        position_front = np.array([0, 270])
        # last_position = position
        timeCount = 0
        scale = 3

        # Iteration in time
        while timeCount <= timeMax:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: sys.exit()

            # Get position, current time and vehicle velocity
            position = initial_position + self.mapi(self.x(timeCount), -self.y(timeCount), xlim, ylim, xlim_inf, ylim_inf, size-np.array([128,128]), realScale=True)     
            tempo = font.render("Tempo : {:.1f}s".format(timeCount), True, (255, 255, 255))
            velocidade = font.render("Velocidade : {:.1f}km/h".format(self.v(timeCount) * 3.6), True, (255, 255, 255))
            pitch = font.render("Arfagem : {:.1f}°".format(self.theta(timeCount)), True, (255, 255, 255))
            ang_roll = self.phi(timeCount) if abs(self.phi(timeCount)) > 0.1 else 0
            roll = font.render("Rolagem : {:.1f}°".format(ang_roll), True, (255, 255, 255))
            escala = font.render("Escala : {:.1f}x".format(scale), True, (255, 255, 255))

            # Blit the data into the screen
            screen.blit(background, (0, 0))
            # screen.fill([211, 211, 211])
            screen.blit(tempo, (5, 0))  
            screen.blit(velocidade, (5, 35))
            screen.blit(pitch, (5, 130)) 
            screen.blit(roll, (5, 255)) 
            screen.blit(escala, (5, 255 + 140))   
            # pygame.draw.aaline(screen, line_colour, last_position, position)
            self.blitRotateCenter(screen, car, position, self.psi(timeCount))
            self.blitRotateCenter(screen, car_side, position_side, -scale*self.theta(timeCount))
            self.blitRotateCenter(screen, car_front, position_front, scale*self.phi(timeCount))
            pygame.display.flip()

            timeCount += 3e-2
            if save:
                pygame.image.save(screen, 'Videos/temp.jpg')
                w.append_data(imageio.imread('Videos/temp.jpg'))
        pygame.display.quit()
        if save: w.close()