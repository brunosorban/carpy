from Function import Function
import numpy as np

class Driver:
    def __init__(self, accelerator='P', steering='steering'):
        self.set_accelerator(accelerator)
        self.set_steering(steering)
        self.plot_steering()
    
    def set_accelerator(self, steer):
        if steer == 'P': 
            self.T1_control = self.proportional_accelerator
            self.T2_control = self.proportional_accelerator
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V

        elif steer == 'PI': 
            self.T1_control = self.PI_throttle
            self.T2_control = self.PI_throttle
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V
            self.integrado = 0
            self.last_time = 0
            self.last_time_control = 0
            self.Vref = 20
            self.last_error = self.Vref

        else: 
            self.T1_control = lambda V, t: 0*V
            self.T2_control = lambda V, t: 0*V
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V

        self.torque_control = np.array([self.T1_control, self.T2_control, self.T3_control, self.T4_control])

    def set_steering(self, steering):
        if steering == 'steering': 
            print('steering defined.')
            self.delta_control = self.steering

        elif steering == 'steering2':
            print('steering2 defined.')
            self.delta_control = self.steering2

        elif steering == 'S': 
            print('S-Curve defined.')
            self.delta_control = self.Scurve

        elif steering == 'SS': 
            print('Smooth S-Curve defined.')
            self.delta_control = self.smoothScurve(name="δ1")

        elif type(steering)==int or type(steering)==float:
            print('Constant radius of {:.2f}m defined.'.format(steering))
            self.delta_control = self.constant_radius(steering)

        else: 
            print('Driver not found. Steering always 0 defined.')
            self.delta_control = lambda t: 0*t

        self.steering_control = self.delta_control

    def get_torque(self, V, t):
        return  self.torque_control[0](V, t), self.torque_control[1](V, t), self.torque_control[2](V, t), self.torque_control[3](V, t)
    
    def get_steering(self, sim_time):
        return self.steering_control(sim_time)
            
    def proportional_accelerator(self, V, t):
        K = 100
        Vref = 10
        dV = Vref - V
        return K * dV if  K * dV < 250 else 250

    def PI_throttle(self, V, t):
        dt = 0.01
        Kp = 1
        Ki = 0.1
        e = self.Vref - V
        if t - self.last_time_control > dt:
            self.integrado += (self.last_error + e) / 2 * dt
            self.last_time_control = t
            self.last_error = e
        return Kp * e + Ki * self.integrado

    def steering(self, sim_time):
        if sim_time <= 30: return 0
        elif 30 < sim_time <= 32: return (sim_time - 30) * 3 / 2 * np.pi / 180
        else: return 3 * np.pi / 180

    def steering2(self, sim_time):
        if sim_time <= 0.1: return 0
        else: return 3 * np.pi / 180
    
    def Scurve(self, sim_time, angle=3):
        if sim_time <= 30: return 0
        elif 30 < sim_time <= 50: return (angle *  np.pi / 180) * np.sin((sim_time - 30) * 2 * np.pi / 20)
        else: return 0

    def smoothScurve(self, angle=3, name=False):
        t = np.linspace(0, 100, 101)
        y = [self.Scurve(ti, angle=angle) for ti in t]
        curve = Function(t, y, xlabel='Time (s)', ylabel=name+' steering input (degrees)', method='cubicSpline', name=name)
        return curve

    def constant_radius(self, r):
        a = 1.07 + 1.605
        s = 1.517
        d1 = lambda t: t**0 * np.arctan(a / r)
        d2 = lambda t: t**0 * np.arctan(a / (r + s))
        return d1, d2

    def plot_steering(self):
        t = np.linspace(0, 100, 10001)
        y1 = [self.steering_control(ti) for ti in t]
        delta = Function(t, 180 / np.pi * np.array(y1), xlabel='Time (s)', ylabel='δ1 steering input (degrees)', method='linear', name='δsw')
        
    def find_radius(self, x, y):
        a=1