from Function import Function
import numpy as np
from simple_pid import PID

class Driver:
    def __init__(self, accelerator='P', steering='steering'):
        self.set_accelerator(accelerator)
        self.set_steering(steering)
        # self.plot_steering()
    
    def set_accelerator(self, accelerator):
        self.accelerator = accelerator
        
        if accelerator == 'P': 
            self.T1_control = self.proportional_accelerator
            self.T2_control = self.proportional_accelerator
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V

        elif accelerator == 'PI': 
            self.T1_control = self.PI_throttle
            self.T2_control = self.PI_throttle
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V
            self.integrado = 0
            self.last_time = 0
            self.last_time_control = 0
            self.Vref = 100 / 3.6
            self.last_error = self.Vref
            
        elif accelerator == 'PID':
            Kp = 150
            Ti = 7
            Td = 0.3
            V = 60 / 3.6
            self.dt = 0.01
            
            self.pid = PID(Kp, Kp * Ti, Kp * Td, setpoint=V, sample_time=None, output_limits = (-500, 500))
            self.last_time_control = 0
            self.output = 0
            self.T1_control = self.PID_throttle
            self.T2_control = self.PID_throttle
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V
            
        elif accelerator == 'PIDSSCD':
            Kp = 150
            Ti = 7
            Td = 0.3
            V = 0
            self.dt = 0.001
            
            self.pid = PID(Kp, Kp * Ti, Kp * Td, setpoint=V, sample_time=None, output_limits = (-500, 500))
            self.last_time_control = 0
            self.output = 0
            self.T1_control = self.PID_throttle_SSCD
            self.T2_control = self.PID_throttle_SSCD
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V  
            
        elif accelerator == 'PIDCSI':
            Kp = 150
            Ti = 7
            Td = 0.3
            V = 100 / 3.6
            self.dt = 0.001
            
            self.pid = PID(Kp, Kp * Ti, Kp * Td, setpoint=V, sample_time=None, output_limits = (-500, 500))
            self.last_time_control = 0
            self.output = 0
            self.T1_control = self.PID_throttle_CSI
            self.T2_control = self.PID_throttle_CSI
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V          

        else: 
            self.T1_control = lambda V, t: 0*V
            self.T2_control = lambda V, t: 0*V
            self.T3_control = lambda V, t: 0*V
            self.T4_control = lambda V, t: 0*V

        self.torque_control = np.array([self.T1_control, self.T2_control, self.T3_control, self.T4_control])

    def set_steering(self, steering):
        self.steering = steering
        
        if steering == 'steering1': 
            print('steering defined.')
            self.delta_control = self.steering1

        elif steering == 'steering2':
            print('steering2 defined.')
            self.delta_control = self.steering2

        elif steering == 'S': 
            print('S-Curve defined.')
            self.delta_control = self.Scurve

        elif steering == 'SS': 
            print('Smooth S-Curve defined.')
            self.delta_control = self.smoothScurve(name="δ1")

        elif steering == 'DLC':
            print('Double Lane-Change defined.')
            def brute_DLC(t, angle=6):
                if t <= 30: return 0
                elif 30 < t <= 31.93: return np.deg2rad(angle) * np.sin((t - 30) * 2 * np.pi /2)
                elif 31.93 < t < 33.4: return 0
                elif 33.4 < t <= 35.318: return -np.deg2rad(angle) * np.sin((t - 33.4) * 2 * np.pi / 2)
                else: return 0
            t = np.linspace(0, 100, 10001)
            y = [brute_DLC(ti) for ti in t]
            self.DLC_func = Function(t, y, method='cubicSpline', xlabel='Time (s)', ylabel='Steering angle (rad)')
            
            self.delta_control = self.DLC
            
            self.DLC_func.plot2D('Steering angle in time', lower=25, upper=40)
            
        elif steering == 'sinusoidal': 
            print('sinusoidal defined.')
            self.delta_control = self.sinusoidal

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
        Vref = 100 / 3.6
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

    def PID_throttle(self, V, t):
        dt = t - self.last_time_control
        if dt > self.dt: 
            self.output = self.pid(V, dt=dt)
            self.last_time_control = t
        return self.output
    
    def PID_throttle_SSCD(self, V, t):
        dt = t - self.last_time_control
        if dt > 0:
            self.pid.setpoint = 0.1 * t 
            self.output = self.pid(V, dt=dt)
            self.last_time_control = t
        return self.output
    
    def PID_throttle_CSI(self, V, t):
        dt = t - self.last_time_control
        if dt > 0 and t < 35:
            self.output = self.pid(V, dt=dt)
            self.last_time_control = t
        return self.output
    
    def steering1(self, sim_time, angle=np.deg2rad(2.9)):
        if sim_time <= 30: return 0
        elif 30 < sim_time <= 30.1: return ((sim_time - 30) / 0.1) * angle
        else: return angle

    def steering2(self, sim_time, angle=np.deg2rad(2.9)):
        if sim_time <= 10: return 0
        elif 10 < sim_time <= 10.1: return ((sim_time - 10) / 0.1) * angle
        else: return angle
    
    def Scurve(self, sim_time, angle=3):
        if sim_time <= 30: return 0
        elif 30 < sim_time <= 50: return (angle *  np.pi / 180) * np.sin((sim_time - 30) * 2 * np.pi / 20)
        else: return 0

    def DLC(self, sim_time):
        # O angulo deve ser definido antes
        return self.DLC_func(sim_time)
        
    def smoothScurve(self, angle=3, name=None):
        t = np.linspace(0, 100, 101)
        y = [self.Scurve(ti, angle=angle) for ti in t]
        curve = Function(t, y, xlabel='Time (s)', ylabel=name+' steering input (degrees)', method='cubicSpline', name=name)
        return curve
    
    def sinusoidal(self, sim_time, angle=np.deg2rad(2.9), freq=0.2):
        return 0 if sim_time < 40 else angle * np.sin(2 * np.pi * freq * sim_time)

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