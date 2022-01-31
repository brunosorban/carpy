import numpy as np
from traitlets.traitlets import ValidateHandler
from Function import Function
import plotly.graph_objects as go

class Track:
    def __init__(self, z_func, vz_func):
        self.z_func = z_func
        self.vz_func = vz_func
        # self.noise_mean = 0
        # self.noise_variance = 0
    
    def __call__(self, x, y):
        return self.get_track_height(x, y)
    
    def get_track_height(self, x, y):
        return self.z_func(x, y) #+ self.get_noise()
    
    def get_track_speed(self, x, y, vx, vy):
        return self.vz_func(x, y, vx, vy)
    
    # def set_noise(self, mean, variance):
    #     self.noise_mean = mean
    #     self.noise_variance = variance
    #     return self.get_noise
    
    # def get_noise(self, x, y):
    #     return self.noise_mean + np.sqrt(self.noise_variance) * np.random.randn()
    
    def profile(self, x=(0, 100), y=(0, 100), resolution=0.1):
        '''Plots the track profile.
        Inputs: 
        x -> Range of the plot in x direction
        y -> Range of the plot in y direction
        resolution -> Distance, in meters, between the dots of the plots
        
        Outputs:
        None
        '''
        X = np.linspace(x[0], x[1], int((x[1] - x[0]) / resolution) + 1)
        Y = np.linspace(y[0], y[1], int((y[1] - y[0]) / resolution) + 1)
        
        [X, Y] = np.meshgrid(X, Y)
        
        # Read data from a csv
        Z = self.get_track_height(X, Y)

        fig = go.Figure(data=[go.Surface(x=X, y=Y, z=Z)])

        fig.update_traces(contours_z=dict(show=True, usecolormap=True,
                                  highlightcolor="limegreen", project_z=True))

        fig.update_layout(title='Track profile', autosize=True)
        # ,width=1500, height=700
        
        # fix the ratio in the plot to follow data size
        fig.update_layout(scene_aspectmode='data')


        fig.show()
        return None
        
# hx = 0.1
# hy = 0
# wx = 100
# wy = 100

# zf = lambda x, y: hx * np.sin(wx * x) + hy * np.sin(wy * y)
# vzf = lambda x, vx, y, vy: hx * wx * vx * np.cos(wx * x) + hy * wy * vy * np.cos(wy * y)

# T = Track(zf, vzf)
# T.profile(x=[0, 100], y=[0,100], resolution=1/10)