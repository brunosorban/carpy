import numpy as np
from Function import Function

class Driver:
    def __init__(self, T1_control='P', T2_control='P', T3_control=None, T4_control=None, delta_1_control='steering', delta_2_control='steering', delta_3_control=None, delta_4_control=None):
        self.T1_control = T1_control
        self.T2_control = T2_control
        self.T3_control = T3_control
        self.T4_control = T4_control
        self.delta_1_control = delta_1_control
        self.delta_2_control = delta_2_control
        self.delta_3_control = delta_3_control
        self.delta_4_control = delta_4_control
        self.set_accelerator()
        self.set_steering()
    
    def set_accelerator(self):
        if self.T1_control == 'P':
            self.T1_control = self.proportional_accelerator
        else:
            self.T1_control = lambda V: 0

        if self.T2_control == 'P':
            self.T2_control = self.proportional_accelerator
        else:
            self.T2_control = lambda V: 0

        if self.T3_control == 'P':
            self.T3_control = self.proportional_accelerator
        else:
            self.T3_control = lambda V: 0

        if self.T4_control == 'P':
            self.T4_control = self.proportional_accelerator
        else:
            self.T4_control = lambda V: 0

    def set_steering(self):
        if self.delta_1_control == 'steering':
            self.delta_1_control = self.steering
        else:
            self.delta_1_control = lambda t: 0

        if self.delta_2_control == 'steering':
            self.delta_2_control = self.steering
        else:
            self.delta_2_control = lambda t: 0

        if self.delta_3_control == 'steering':
            self.delta_3_control = self.steering
        else:
            self.delta_3_control = lambda t: 0

        if self.delta_4_control == 'steering':
            self.delta_4_control = self.steering
        else:
            self.delta_4_control = lambda t: 0

    def proportional_accelerator(self, V):
        K = 100
        Vref = 10
        dV = Vref - V
        return K * dV

    def steering(self, time):
        if time <= 30: return 0
        elif 30 < time <= 32: return (time - 30) * 3 / 2 * np.pi / 180;
        else: return 3 * np.pi / 180