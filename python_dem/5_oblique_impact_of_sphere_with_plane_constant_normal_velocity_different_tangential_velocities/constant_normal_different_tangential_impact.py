import numpy as np
import matplotlib.pyplot as plt
from numpy import log, pi
import os


class Particle():
    def __init__(self, x, y, z, radius, rho):
        self.x = x
        self.y = y
        self.z = z
        self.u = 0.
        self.v = 0.
        self.w = 0.
        # angular velocity
        self.wx = 0.
        self.wy = 0.
        self.wz = 0.

        self.overlap = 0.
        self.fn = 0.
        self.fx = 0.
        self.fy = 0.
        self.fz = 0.
        self.tor_x = 0.
        self.tor_y = 0.
        self.tor_z = 0.

        self.rho = rho
        self.radius = radius
        self.m = 4. / 3. * np.pi * radius**3. * self.rho
        # moment of inertia
        self.moi = 2. / 5. * self.m * radius**2.
        self.nu = 0.
        self.E = 0.
        self.G = 0.


class Wall():
    def __init__(self, x, y, z, normal_x, normal_y, normal_z):
        self.x = x
        self.y = y
        self.z = z
        self.normal_x = normal_x
        self.normal_y = normal_y
        self.normal_z = normal_z


def update_velocity(particle, dt):
    particle.u += particle.fx/particle.m * dt
    particle.v += particle.fy/particle.m * dt
    particle.w += particle.fz/particle.m * dt
    particle.wx += particle.tor_x/particle.moi * dt
    particle.wy += particle.tor_y/particle.moi * dt
    particle.wz += particle.tor_z/particle.moi * dt


def update_position(particle, dt):
    particle.x += particle.u * dt
    particle.y += particle.v * dt
    particle.z += particle.w * dt

def make_forces_zero(particle):
    # compute the force
    particle.fx = 0.
    particle.fy = 0.
    particle.fz = 0.
    particle.tor_x = 0.
    particle.tor_y = 0.
    particle.tor_z = 0.

def body_force(particle, gx, gy, gz):
    # compute the force
    particle.fx += particle.m * gx
    particle.fy += particle.m * gy
    particle.fz += particle.m * gz


def distance(particle_1, particle_2):
    # compute the overlap
    xij = particle_1.x - particle_2.x
    yij = particle_1.y - particle_2.y
    zij = particle_1.z - particle_2.z
    return (xij**2. + yij**2. + zij**2.)**0.5


def distance(particle, wall):
    # compute the overlap
    xij = particle_1.x - particle_2.x
    yij = particle_1.y - particle_2.y
    zij = particle_1.z - particle_2.z
    return (xij**2. + yij**2. + zij**2.)**0.5


def compute_force_on_particle_due_to_wall(particle_1, wall, cor, fric_coeff, dt):
    # compute the overlap
    xij = particle_1.x - wall.x
    yij = particle_1.y - wall.y
    zij = particle_1.z - wall.z

    nij_x = wall.normal_x
    nij_y = wall.normal_y
    nij_z = wall.normal_z

    tmp = xij * nij_x  + yij * nij_y + zij * nij_z
    overlap = particle_1.radius - tmp

    # n \cross omega * R
    uij = particle_1.u + (nij_y * particle_1.wz - nij_z * particle_1.wy) * (particle_1.radius - overlap)
    vij = particle_1.v - (nij_x * particle_1.wz - nij_z * particle_1.wx) * (particle_1.radius - overlap)
    wij = particle_1.w + (nij_x * particle_1.wy - nij_y * particle_1.wx) * (particle_1.radius - overlap)

    vn = uij * nij_x + vij * nij_y + wij * nij_z
    vn_x = vn * nij_x
    vn_y = vn * nij_y
    vn_z = vn * nij_z

    if overlap > 0:
        ############################
        # normal force computation #
        ############################
        # Compute stiffness
        # effective Young's modulus
        tmp_1 = (1. - particle_1.nu**2.) / particle_1.E
        # tmp_2 = (1. - particle_2.nu**2.) / particle_2.E
        tmp_2 = 0.
        E_eff = 1. / (tmp_1 + tmp_2)
        tmp_1 = 1. / particle_1.radius
        # tmp_2 = 1. / particle_2.radius
        tmp_2 = 0.
        R_eff = 1. / (tmp_1 + tmp_2)
        # Eq 4 [1]
        kn = 4. / 3. * E_eff * R_eff**0.5

        # compute damping coefficient
        tmp_1 = log(cor)
        tmp_2 = log(cor)**2. + pi**2.
        alpha_1 = -tmp_1 * (5. / tmp_2)**0.5
        tmp_1 = 1. / particle_1.m
        # tmp_2 = 1. / particle_2.m
        tmp_2 = 0.
        m_eff = 1. / (tmp_1 + tmp_2)
        eta = alpha_1 * (m_eff * kn)**0.5 * overlap**0.25

        fn = kn * overlap**1.5
        fn_x = fn * nij_x - eta * vn_x
        fn_y = fn * nij_y - eta * vn_y
        fn_z = fn * nij_z - eta * vn_z

        particle_1.fn = fn
        particle_1.overlap = overlap
        particle_1.fx += fn_x
        particle_1.fy += fn_y
        particle_1.fz += fn_z

        # tangential unit vector
        vt_x = uij - vn * nij_x
        vt_y = vij - vn * nij_y
        vt_z = wij - vn * nij_z
        vij_t = (vt_x**2. + vt_y**2. + vt_z**2.)**0.5
        if vij_t > 1e-12:
            tij_x = vt_x / vij_t
            tij_y = vt_y / vij_t
            tij_z = vt_z / vij_t
        else:
            tij_x = 0.
            tij_y = 0.
            tij_z = 0.

        if vij_t < 1e-12:
            particle_1.ft_x = 0.
            particle_1.ft_y = 0.
            particle_1.ft_z = 0.

            # reset the spring length
            particle_1.tng_disp_x = 0.
            particle_1.tng_disp_y = 0.
            particle_1.tng_disp_z = 0.
        else:
            # compute the tangential force, by equation 27
            tng_disp_x = particle_1.tng_disp_x + uij * dt
            tng_disp_y = particle_1.tng_disp_y + vij * dt
            tng_disp_z = particle_1.tng_disp_z + wij * dt
            tng_disp_magn = tng_disp_x * tij_x + tng_disp_y * tij_y + tng_disp_z * tij_z
            particle_1.tng_disp_x = tng_disp_magn * tij_x
            particle_1.tng_disp_y = tng_disp_magn * tij_y
            particle_1.tng_disp_z = tng_disp_magn * tij_z

            # Compute the tangential stiffness
            tmp_1 = (2. - particle_1.nu) / particle_1.G
            # tmp_2 = (2. - s_nu[s_idx]) / s_G[s_idx]
            tmp_2 = 0.
            G_eff = 1. / (tmp_1 + tmp_2)
            kt = 8. * G_eff * (R_eff * overlap)**0.5
            tmp_1 = log(cor)
            tmp_2 = log(cor)**2. + pi**2.
            beta = -tmp_1 / tmp_2**0.5
            eta_t = 2. * beta * (5/6)**0.5 * (kt * m_eff)**0.5

            ft_x_star = - kt * tng_disp_x - eta_t * vt_x
            ft_y_star = - kt * tng_disp_y - eta_t * vt_y
            ft_z_star = - kt * tng_disp_z - eta_t * vt_z

            # apply coulomb law
            ft_magn = (ft_x_star**2. + ft_y_star**2. + ft_z_star**2.)**0.5
            fn_magn = fn

            ft_magn_star = min(fric_coeff * fn_magn, ft_magn)
            # compute the tangential force, by equation 27
            particle_1.ft_x = -ft_magn_star * tij_x
            particle_1.ft_y = -ft_magn_star * tij_y
            particle_1.ft_z = -ft_magn_star * tij_z
            particle_1.fx += particle_1.ft_x
            particle_1.fy += particle_1.ft_y
            particle_1.fz += particle_1.ft_z

            # reset the spring length
            particle_1.tng_disp_x = -particle_1.ft_x / kt
            particle_1.tng_disp_y = -particle_1.ft_y / kt
            particle_1.tng_disp_z = -particle_1.ft_z / kt

        # Compute torque
        particle_1.tor_x = -(nij_y * particle_1.ft_z - nij_z * particle_1.ft_y) * (particle_1.radius - overlap)
        particle_1.tor_y = (nij_x * particle_1.ft_z - nij_z * particle_1.ft_x) * (particle_1.radius - overlap)
        particle_1.tor_z = -(nij_x * particle_1.ft_y - nij_y * particle_1.ft_x) * (particle_1.radius - overlap)

    else:
        fn = 0.
        particle_1.fn = 0.
        particle_1.fx = 0.
        particle_1.fy = 0.
        particle_1.fz = 0.

        particle_1.tng_disp_x = 0.
        particle_1.tng_disp_y = 0.
        particle_1.tng_disp_z = 0.

        particle_1.ft_x = 0.
        particle_1.ft_y = 0.
        particle_1.ft_z = 0.
        particle_1.tor_x = 0.
        particle_1.tor_y = 0.
        particle_1.tor_z = 0.


def get_shear_modulus(E, nu):
    return E / (2. * (1. + nu))


def run():
    # vt = np.arange(0.1, 70., 10.)
    vt_1 = np.linspace(0.1, 12, 30)
    vt_2 = np.linspace(12.1, 30., 4)
    vt = np.concatenate((vt_1, vt_2))
    vn = 5
    ang_vel_dem = []
    # create the particle
    for _vt in vt:
        E = 2.08 * 1e11
        nu = 0.3
        G = get_shear_modulus(E, nu)
        fric_coeff = 0.3
        cor = 1.
        density = 7850
        radius = 1. * 1e-5
        velocity = 5
        particle_1 = Particle(0., radius + radius/10000, 0.0, radius, density)
        wall = Wall(0., 0., 0., 0., 1., 0.)
        particle_1.E = E
        particle_1.nu = nu
        particle_1.G = G
        # wall.E = E
        # wall.nu = nu
        # wall.G = G
        particle_1.u = _vt
        particle_1.v = -vn
        dt = 1e-9
        _t = 0
        tf = 10 * 1e-8
        gx = 0.
        gy = 0.
        gz = 0.
        t = []
        u = []
        y = []
        v = []
        fn = []
        overlap = []
        fy = []
        fx = []
        wz = []
        while _t < tf:
            y.append(particle_1.y)
            u.append(particle_1.u)
            v.append(particle_1.v)
            fn.append(particle_1.fn)
            overlap.append(particle_1.overlap)
            fy.append(particle_1.fy)
            fx.append(particle_1.fx)
            wz.append(particle_1.wz)
            t.append(_t)

            update_velocity(particle_1, dt)
            # update_velocity(particle_2, dt)
            make_forces_zero(particle_1)
            # make_forces_zero(particle_2)
            compute_force_on_particle_due_to_wall(particle_1, wall, cor, fric_coeff, dt)
            # compute_force_on_particles_due_to_particles(particle_2, particle_1, cor, fric_coeff, dt)
            body_force(particle_1, gx, gy, gz)
            # body_force(particle_2, gx, gy, gz)
            update_position(particle_1, dt)
            # update_position(particle_2, dt)
            # print(_t)
            _t += dt

        ang_vel_dem.append(particle_1.wz)
        # print("wz is", particle_1.wz)

    # convert list into numpy arrays
    ang_vel_dem = particle_1.radius * np.asarray(ang_vel_dem) / (fric_coeff * vn)
    thetas = vt / (fric_coeff * vn)

    plt.clf()
    path = os.path.abspath(__file__)
    directory = os.path.dirname(path)
    # print(path)
    data_wz_vs_angle_fem = np.loadtxt(
        os.path.join(
            directory,
            'fem_wu_2003.csv'),
        delimiter=','
    )
    angle_fem, wz_fem = data_wz_vs_angle_fem[:, 0], data_wz_vs_angle_fem[:, 1]
    plt.scatter(angle_fem, wz_fem, label='FEM')

    data_wz_vs_angle_dem_chung = np.loadtxt(
        os.path.join(
            directory,
            'dem_steel_chung.csv'),
        delimiter=','
    )
    angle_dem_chung, wz_dem_chung = data_wz_vs_angle_dem_chung[:, 0], data_wz_vs_angle_dem_chung[:, 1]
    plt.plot(angle_dem_chung, wz_dem_chung, label='DEM Chung')

    plt.plot(thetas, ang_vel_dem, "*-", label='DEM')
    # plt.scatter(cors, cor_dem, label='DEM')
    plt.xlabel('Vt/(fric x Vn)')
    plt.ylabel('r ang_vel / (fric x Vn)')
    plt.legend()
    fig = os.path.join(directory, "angle_vs_repose_omega.pdf")
    plt.savefig(fig, dpi=300)

run()
