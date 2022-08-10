import numpy as np


class Vehicle:
    def __init__(
        self,
        Tire,
        vehicle_mass,
        Ixx,
        Iyy,
        Izz,
        lf,
        lr,
        wf,
        wr,
        af,
        cd,
        CG_height,
        rho=1.225,
    ):
        self.Tire = Tire
        self.n_tires = 4
        self.vehicle_mass = vehicle_mass
        self.vehicle_equivalent_mass = (
            self.vehicle_mass + self.n_tires * Tire.Jzz / Tire.tire_radius ** 2
        )
        self.lf = lf
        self.lr = lr
        self.wf = wf
        self.wr = wr
        self.af = af
        self.cd = cd
        self.rho = rho
        self.CG_height = CG_height
        self.Ixx = Ixx
        self.Iyy = Iyy
        self.Izz = Izz
        self.K_arf = 0  # anti-roll bar parameters
        self.K_arr = 0  # anti-roll bar parameters
        self.gamma1 = 0
        self.gamma2 = 0
        self.gamma3 = 0
        self.gamma4 = 0
        self.gravity = 9.81

        # Pre-calculations
        self.h = CG_height - Tire.tire_radius
        self.Ff0 = (
            self.lr / (2 * (self.lf + self.lr)) * self.vehicle_mass * self.gravity
        )
        self.Fr0 = (
            self.lf / (2 * (self.lf + self.lr)) * self.vehicle_mass * self.gravity
        )

    def set_suspension(self, K_sf, K_sr, C_sf, C_sr):
        self.K_sf = K_sf
        self.K_sr = K_sr
        self.C_sf = C_sf
        self.C_sr = C_sr

    def set_anti_roll_bar(self, position, d=0, a=0, b=0, G=0, K_arz=False):
        K_arz = G * (np.pi * d ** 4 / 32) / (b * a ** 2) if not K_arz else K_arz

        if position == "f":
            self.K_arf = K_arz
            print("Anti-roll Bar (front) = {:.1f} Nm/৹".format(np.deg2rad(K_arz)))
        elif position == "r":
            self.K_arr = K_arz
            print("Anti-roll Bar (rear)  = {:.1f} Nm/৹".format(np.deg2rad(K_arz)))
        else:
            return print('Please insert a valid position. Position must be "f" or "r"')

        return K_arz

    def set_camber(self, gamma1, gamma2, gamma3, gamma4):
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.gamma3 = gamma3
        self.gamma4 = gamma4

    def get_vertical_load(
        self,
        z,
        vz,
        phi,
        theta,
        phi_dot,
        theta_dot,
        z1,
        z2,
        z3,
        z4,
        vz1,
        vz2,
        vz3,
        vz4,
        zc1,
        zc2,
        zc3,
        zc4,
        vzc1,
        vzc2,
        vzc3,
        vzc4,
    ):
        # Retrieve suspension properties
        Kt = self.Tire.cz
        Ct = 0  # tire without damping

        # Get static forces distribution
        F10 = F20 = self.Ff0
        F30 = F40 = self.Fr0

        # Calculate sines used
        sphi = np.sin(phi)
        stheta = np.sin(theta)
        cphi = np.cos(phi)
        ctheta = np.cos(theta)

        # Calculate tire forces
        Ft1 = F10 + Kt * (zc1 - z1) + Ct * (vzc1 - vz1)
        Ft2 = F20 + Kt * (zc2 - z2) + Ct * (vzc2 - vz2)
        Ft3 = F30 + Kt * (zc3 - z3) + Ct * (vzc3 - vz3)
        Ft4 = F40 + Kt * (zc4 - z4) + Ct * (vzc4 - vz4)

        # Calculate suspension forces
        Fs1 = (
            F10
            + self.K_sf * (z1 - z - self.wf / 2 * sphi + self.lf * stheta)
            + self.C_sf
            * (vz1 - vz - phi_dot * self.wf / 2 * cphi + theta_dot * self.lf * ctheta)
        )
        Fs2 = (
            F20
            + self.K_sf * (z2 - z + self.wf / 2 * sphi + self.lf * stheta)
            + self.C_sf
            * (vz2 - vz + phi_dot * self.wf / 2 * cphi + theta_dot * self.lf * ctheta)
        )
        Fs3 = (
            F30
            + self.K_sr * (z3 - z - self.wr / 2 * sphi - self.lr * stheta)
            + self.C_sr
            * (vz3 - vz - phi_dot * self.wr / 2 * cphi - theta_dot * self.lr * ctheta)
        )
        Fs4 = (
            F40
            + self.K_sr * (z4 - z + self.wr / 2 * sphi - self.lr * stheta)
            + self.C_sr
            * (vz4 - vz + phi_dot * self.wr / 2 * cphi - theta_dot * self.lr * ctheta)
        )

        # Calculate anti-roll bar forces
        Farf = self.K_arf * (z2 - z1)  # The force is positive in 1 and negative in 2
        Farr = self.K_arr * (z4 - z3)  # The force is positive in 3 and negative in 4

        return Ft1, Ft2, Ft3, Ft4, Fs1, Fs2, Fs3, Fs4, Farf, Farr

    def get_steer_wbw(self, delta_sw):
        # return the ackerman geometry wheel by wheel
        if delta_sw == 0:
            return 0, 0, 0, 0
        else:
            a = self.lf + self.lr
            s = self.wf
            r = a / np.tan(abs(delta_sw))
            d1 = np.arctan(a / r)
            d2 = np.arctan(a / (r + s))
            if delta_sw > 0:
                return d1, d2, 0, 0
            else:
                return -d2, -d1, 0, 0
