import numpy as np
from Function import Function
from scipy import integrate

# For the animation class
import sys, pygame

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
        self.atol = atol
        self.rtol = rtol
        self.maxStep = maxStep
        self.solver()
    
    def uDot(self, t, sol, post_process=False):
        '''Calculate the state space derivatives'''
        # Retrieve integration data
        x, y, vx, vy, psi, psi_dot, omega1, omega2, omega3, omega4 = sol

        # Get velocities in vehicle coordinate system
        vvx = vx * np.cos(psi) + vy * np.sin(psi)
        vvy = -vx * np.sin(psi) + vy * np.cos(psi)
        v = np.sqrt(vx**2 + vy**2)

        # Get control variables
        T1 = self.Driver.T1_control(v)
        T2 = self.Driver.T2_control(v)
        T3 = self.Driver.T3_control(v)
        T4 = self.Driver.T4_control(v)

        delta_1 = self.Driver.delta_1_control(t)
        delta_2 = self.Driver.delta_2_control(t)
        delta_3 = self.Driver.delta_3_control(t)
        delta_4 = self.Driver.delta_4_control(t)

        # Calculate tire forces
        ftx_1, fty_1, Mz_1 = self.Vehicle.Tire(vvx - psi_dot * self.Vehicle.wf/2, vvy + psi_dot * self.Vehicle.lf, delta_1, omega1, self.Vehicle.Fz / self.Vehicle.n_tires)
        ftx_2, fty_2, Mz_2 = self.Vehicle.Tire(vvx + psi_dot * self.Vehicle.wf/2, vvy + psi_dot * self.Vehicle.lf, delta_2, omega2, self.Vehicle.Fz / self.Vehicle.n_tires)
        ftx_3, fty_3, Mz_3 = self.Vehicle.Tire(vvx - psi_dot * self.Vehicle.wr/2, vvy - psi_dot * self.Vehicle.lr, delta_3, omega3, self.Vehicle.Fz / self.Vehicle.n_tires)
        ftx_4, fty_4, Mz_4 = self.Vehicle.Tire(vvx + psi_dot * self.Vehicle.wr/2, vvy - psi_dot * self.Vehicle.lr, delta_4, omega4, self.Vehicle.Fz / self.Vehicle.n_tires)

        # Convert tire forces to vehicle coordinate system
        fvx_1 = ftx_1 * np.cos(delta_1) - fty_1 * np.sin(delta_1) 
        fvy_1 = ftx_1 * np.sin(delta_1) + fty_1 * np.cos(delta_1)

        fvx_2 = ftx_2 * np.cos(delta_2) - fty_2 * np.sin(delta_2)
        fvy_2 = ftx_2 * np.sin(delta_2) + fty_2 * np.cos(delta_2)

        fvx_3 = ftx_3 * np.cos(delta_3) - fty_3 * np.sin(delta_3) 
        fvy_3 = ftx_3 * np.sin(delta_3) + fty_3 * np.cos(delta_3)

        fvx_4 = ftx_4 * np.cos(delta_4) - fty_4 * np.sin(delta_4)
        fvy_4 = ftx_4 * np.sin(delta_4) + fty_4 * np.cos(delta_4)

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
        psi_dot_dot = tau_z / self.Vehicle.Izz

        # Calculate accelerations in the wheels
        # Deveria aqui ser considerado o raio dinâmico?
        omega1_dot = (T1 - ftx_1 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz
        omega2_dot = (T2 - ftx_2 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz
        omega3_dot = (T3 - ftx_3 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz
        omega4_dot = (T4 - ftx_4 * self.Vehicle.Tire.tire_radius) / self.Vehicle.Tire.Jzz

        if not post_process:
            uDot = [
                vx,
                vy,
                ax,
                ay,
                psi_dot,
                psi_dot_dot,
                omega1_dot,
                omega2_dot,
                omega3_dot,
                omega4_dot
            ]
            return uDot   

        else:
            self.T1.append(T1)
            self.T2.append(T2)
            self.T3.append(T3)
            self.T4.append(T4)
            self.delta_1.append(delta_1)
            self.delta_2.append(delta_2)
            self.delta_3.append(delta_3)
            self.delta_4.append(delta_4)
            self.ftx_1.append(ftx_1)
            self.fty_1.append(fty_1)
            self.Mz_1.append(Mz_1)
            self.ftx_2.append(ftx_2)
            self.fty_2.append(fty_2)
            self.Mz_2.append(Mz_2)
            self.ftx_3.append(ftx_3)
            self.fty_3.append(fty_3)
            self.Mz_3.append(Mz_3)
            self.ftx_4.append(ftx_4)
            self.fty_4.append(fty_4)
            self.Mz_4.append(Mz_4)
            self.Frr.append(Frr)
            self.Fd.append(Fd)

            self.x.append(x)
            self.y.append(y)
            self.vx.append(vx)
            self.vy.append(vy)
            self.vvx.append(vvx)
            self.vvy.append(vvy)
            self.ax.append(ax)
            self.ay.append(ay)
            self.avx.append(avx)
            self.avy.append(avy)
            self.psi.append(psi)
            self.psi_dot.append(psi_dot)
            self.psi_dot_dot.append(psi_dot_dot)
            self.omega1.append(omega1)
            self.omega2.append(omega2)
            self.omega3.append(omega3)
            self.omega4.append(omega4)
            self.omega1_dot.append(omega1_dot)
            self.omega2_dot.append(omega2_dot)
            self.omega3_dot.append(omega3_dot)
            self.omega4_dot.append(omega4_dot)

    def solver(self):
        self.time = [0]
        self.solution = [[0, 0, 0, 0, 0, 0, 0, 0]]
        solver = integrate.LSODA(self.uDot, t0=0, y0=self.initialSolution, max_step=self.maxStep, t_bound=self.maxTime, rtol=self.rtol, atol=self.atol)
        while solver.status != 'finished':
            solver.step()
            self.time.append(solver.t)
            self.solution.append(solver.y)
        print("Solution Finished")

    def post_process(self):
        self.T1 = [0]
        self.T2 = [0]
        self.T3 = [0]
        self.T4 = [0]
        self.delta_1 = [0]
        self.delta_2 = [0]
        self.delta_3 = [0]
        self.delta_4 = [0]
        self.ftx_1 = [0]
        self.fty_1 = [0]
        self.Mz_1 = [0]
        self.ftx_2 = [0]
        self.fty_2 = [0]
        self.Mz_2 = [0]
        self.ftx_3 = [0]
        self.fty_3 = [0]
        self.Mz_3 = [0]
        self.ftx_4 = [0]
        self.fty_4 = [0]
        self.Mz_4 = [0]
        self.Frr = [0]
        self.Fd = [0]

        self.x = [0]
        self.y = [0]
        self.vx = [0]
        self.vy = [0]
        self.vvx = [0]
        self.vvy = [0]
        self.ax = [0]
        self.ay = [0]
        self.avx = [0]
        self.avy = [0]
        self.psi = [0]
        self.psi_dot = [0]
        self.psi_dot_dot = [0]
        self.omega1 = [0]
        self.omega2 = [0]
        self.omega3 = [0]
        self.omega4 = [0]
        self.omega1_dot = [0]
        self.omega2_dot = [0]
        self.omega3_dot = [0]
        self.omega4_dot = [0]

        for i in range(1, len(self.solution)):
            # print(self.solution)
            # print(self.time)
            self.uDot(self.time[i], self.solution[i], post_process=True)

        self.T1 = Function(self.time, self.T1, xlabel='Time (s)', ylabel='Front Left tire torque (Nm)', name='T1')
        self.T2 = Function(self.time, self.T2, xlabel='Time (s)', ylabel='Front Right tire torque (Nm)', name='T2')
        self.T3 = Function(self.time, self.T3, xlabel='Time (s)', ylabel='Rear Left tire torque (Nm)', name='T3')
        self.T4 = Function(self.time, self.T4, xlabel='Time (s)', ylabel='Rear Right tire torque (Nm)', name='T4')
        self.delta_1 = Function(self.time, 180 / np.pi * np.array(self.delta_1), xlabel='Time (s)', ylabel='Front Left tire steering angle (degrees)', name='δ1')
        self.delta_2 = Function(self.time, 180 / np.pi * np.array(self.delta_2), xlabel='Time (s)', ylabel='Front Right tire steering angle (degrees)', name='δ2')
        self.delta_3 = Function(self.time, 180 / np.pi * np.array(self.delta_3), xlabel='Time (s)', ylabel='Rear Left tire steering angle (degrees)', name='δ3')
        self.delta_4 = Function(self.time, 180 / np.pi * np.array(self.delta_4), xlabel='Time (s)', ylabel='Rear Right tire steering angle (degrees)', name='δ4')
        self.ftx_1 = Function(self.time, self.ftx_1, xlabel='Time (s)', ylabel='Front Left tire force in tire x axis (Nm)', name='ftx_1')
        self.fty_1 = Function(self.time, self.fty_1, xlabel='Time (s)', ylabel='Front Left tire force in tire y axis (Nm)', name='fty_1')
        self.Mz_1 = Function(self.time, self.Mz_1, xlabel='Time (s)', ylabel='Front Left tire moment in tire z axis (Nm)', name='Mz_1')
        self.ftx_2 = Function(self.time, self.ftx_2, xlabel='Time (s)', ylabel='Front Right tire force in tire x axis (Nm)', name='ftx_2')
        self.fty_2 = Function(self.time, self.fty_2, xlabel='Time (s)', ylabel='Front Right tire force in tire y axis (Nm)', name='fty_2')
        self.Mz_2 = Function(self.time, self.Mz_2, xlabel='Time (s)', ylabel='Front Right tire moment in tire z axis (Nm)', name='Mz_2')
        self.ftx_3 = Function(self.time, self.ftx_3, xlabel='Time (s)', ylabel='Rear Left tire force in tire x axis (Nm)', name='ftx_3')
        self.fty_3 = Function(self.time, self.fty_3, xlabel='Time (s)', ylabel='Rear Left tire force in tire y axis (Nm)', name='fty_3')
        self.Mz_3 = Function(self.time, self.Mz_3, xlabel='Time (s)', ylabel='Rear Left tire moment in tire z axis (Nm)', name='Mz_3')
        self.ftx_4 = Function(self.time, self.ftx_4, xlabel='Time (s)', ylabel='Rear Right tire force in tire x axis (Nm)', name='ftx_4')
        self.fty_4 = Function(self.time, self.fty_4, xlabel='Time (s)', ylabel='Rear Right tire force in tire y axis (Nm)', name='fty_4')
        self.Mz_4 = Function(self.time, self.Mz_4, xlabel='Time (s)', ylabel='Rear Right tire moment in tire z axis (Nm)', name='Mz_4')
        self.Frr = Function(self.time, self.Frr, xlabel='Time (s)', ylabel='Total roling resistance (N)', name='Frr')
        self.Fd = Function(self.time, self.Fd, xlabel='Time (s)', ylabel='Drag force (N)', name='Fd')

        self.r = Function(self.time, np.sqrt(np.array(self.x)**2 + np.array(self.y)**2), xlabel='Time (s)', ylabel='Distance from the origin (m)', name='r')
        self.v = Function(self.time, np.sqrt(np.array(self.vx)**2 + np.array(self.vy)**2), xlabel='Time (s)', ylabel='Velocity (m)', name='v')
        self.xy = Function(self.x, self.y, xlabel='Position in inertial frame x axis (m)', ylabel='Position in inertial frame y axis (m)', name='xy')
        self.x = Function(self.time, self.x, xlabel='Time (s)', ylabel='Position in inertial frame x axis (m)', name='x')
        self.y = Function(self.time, self.y, xlabel='Time (s)', ylabel='Position in inertial frame y axis (m)', name='y')
        self.vx = Function(self.time, self.vx, xlabel='Time (s)', ylabel='Velocity in inertial frame x axis (m/s)', name='vx')
        self.vy = Function(self.time, self.vy, xlabel='Time (s)', ylabel='Velocity in inertial frame y axis (m/s)', name='vy')
        self.vvx = Function(self.time, self.vvx, xlabel='Time (s)', ylabel='Velocity in vehicle frame x axis (m/s)', name='vvx')
        self.vvy = Function(self.time, self.vvy, xlabel='Time (s)', ylabel='Velocity in vehicle frame y axis (m/s)', name='vvy')
        self.ax = Function(self.time, self.ax, xlabel='Time (s)', ylabel='Acceleration in inertial frame x axis (m/s²)', name='ax')
        self.ay = Function(self.time, self.ay, xlabel='Time (s)', ylabel='Acceleration in inertial frame y axis (m/s²)', name='ay')
        self.avx = Function(self.time, self.avx, xlabel='Time (s)', ylabel='Acceleration in vehicle frame x axis (m/s²)', name='avx')
        self.avy = Function(self.time, self.avy, xlabel='Time (s)', ylabel='Acceleration in vehicle frame y axis (m/s²)', name='avy')
        self.psi = Function(self.time,  180 / np.pi * np.array(self.psi), xlabel='Time (s)', ylabel='Picht angle (degrees)', name='psi')
        self.psi_dot = Function(self.time, self.psi_dot, xlabel='Time (s)', ylabel='Picht angular velocity (rad/s)', name='psi_dot')
        self.psi_dot_dot = Function(self.time, self.psi_dot_dot, xlabel='Time (s)', ylabel='Picht angular acceleration (rad/s²)', name='psi_dot_dot')
        self.omega1= Function(self.time, self.omega1, xlabel='Time (s)', ylabel='Front Left tire angular velocity (rad/s)', name='ω1')
        self.omega2= Function(self.time, self.omega2, xlabel='Time (s)', ylabel='Front Right tire angular velocity (rad/s)', name='ω2')
        self.omega3 = Function(self.time, self.omega3, xlabel='Time (s)', ylabel='Rear Left tire angular velocity (rad/s)', name='ω3')
        self.omega4 = Function(self.time, self.omega4, xlabel='Time (s)', ylabel='Rear Right tire angular velocity (rad/s)', name='ω4')
        self.omega1_dot = Function(self.time, self.omega1_dot, xlabel='Time (s)', ylabel='Front Left tire angular acceleration (rad/s²)', name='ω1_dot')
        self.omega2_dot = Function(self.time, self.omega2_dot, xlabel='Time (s)', ylabel='Front Right tire angular acceleration (rad/s²)', name='ω2_dot')
        self.omega3_dot = Function(self.time, self.omega3_dot, xlabel='Time (s)', ylabel='Rear Left tire angular acceleration (rad/s²)', name='ω3_dot')
        self.omega4_dot = Function(self.time, self.omega4_dot, xlabel='Time (s)', ylabel='Rear Right tire angular acceleration (rad/s²)', name='ω4_dot')


    def mapi(self, x, y, xlim, ylim, size):
        x_factor = size[0] / xlim
        y_factor = size[1] / ylim
        return [x * x_factor, y * y_factor]

    def animate(self, timeMax):
        # Initialization
        xlim = max(self.x.__Y_source__)
        ylim = max(self.y.__Y_source__) + 1e-3
        pygame.init()
        # font = pygame.font.SysFont('Helvetica', 28)

        # Defining auxiliar colors 
        # white = 0, 0, 0
        
        # Creating animation screen
        size = 1920, 1080
        screen = pygame.display.set_mode(size)
        
        # Preparing images
        background = pygame.image.load('../Animation/TrackT01.png')
        background = pygame.transform.scale(background, (size))

        car = pygame.image.load('../Animation/Car.png')
        car = pygame.transform.scale(car, (128, 128))
        position = [0, size[1]-128*1.5]

        timeCount = 0
        while timeCount <= timeMax:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: sys.exit()

            screen.blit(background, [0, 0])

            position = [0, size[1]-128*1.5] + self.mapi(self.x(timeCount), -self.y(timeCount), xlim, ylim, size)          
            screen.blit(car, position)
            pygame.display.flip()
            timeCount += 1e-2