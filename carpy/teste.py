import numpy as np
phi, theta, psi = 0, 0, 60 * np.pi / 180
c11 = np.cos(theta) * np.cos(psi)
c12 = np.cos(theta) * np.sin(psi)
c13 = -np.sin(theta)
c21 = np.sin(phi) * np.sin(theta) * np.cos(psi) - np.cos(phi) * np.sin(psi)
c22 = np.sin(phi) * np.sin(theta) * np.sin(psi) + np.cos(phi) * np.cos(psi)
c23 = np.sin(phi) * np.cos(theta)
c31 = np.cos(phi) * np.sin(theta) * np.cos(psi) + np.sin(phi) * np.sin(psi)
c32 = np.cos(phi) * np.sin(theta) * np.sin(psi) - np.sin(phi) * np.cos(psi)
c33 = np.cos(phi) * np.cos(theta)

C = np.array([[c11, c12, c13], [c21, c22, c23], [c31, c32, c33]])

print(C @ np.array([0, 0, 1]))