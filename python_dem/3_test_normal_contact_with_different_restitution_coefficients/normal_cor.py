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

    uij = particle_1.u - (0)
    vij = particle_1.v - (0)
    wij = particle_1.w - (0)

    vn = uij * nij_x + vij * nij_y + wij * nij_z
    vn_x = vn * nij_x
    vn_y = vn * nij_y
    vn_z = vn * nij_z

    tmp = xij * nij_x  + yij * nij_y + zij * nij_z
    overlap = particle_1.radius - tmp

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

    else:
        fn = 0.
        particle_1.fn = 0.
        particle_1.fx = 0.
        particle_1.fy = 0.
        particle_1.fy = 0.


def get_shear_modulus(E, nu):
    return E / (2. * (1. + nu))


def run():
    cors = np.linspace(0.01, 1., 6)
    cor_dem = []
    # create the particle
    for _cor in cors:
        E = 3.8 * 1e11
        nu = 0.23
        G = get_shear_modulus(E, nu)
        fric_coeff = 0.0
        cor = _cor
        density = 4000
        radius = 0.0025
        velocity = 3.9
        particle_1 = Particle(0., radius + radius/10000, 0.0, radius, density)
        wall = Wall(0., 0., 0., 0., 1., 0.)
        particle_1.E = E
        particle_1.nu = nu
        particle_1.G = G
        # wall.E = E
        # wall.nu = nu
        # wall.G = G
        particle_1.v = -velocity
        dt = 1e-8
        _t = 0
        tf = 1000 * 1e-6
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

        final_vel = particle_1.v
        cor_computed = final_vel / velocity
        cor_dem.append(cor_computed)
        print("cor dem is", cor_computed)

    # convert list into numpy arrays
    cor_dem = np.asarray(cor_dem)

    plt.clf()
    path = os.path.abspath(__file__)
    directory = os.path.dirname(path)
    plt.plot(cors, cors, label='Analytical')
    plt.scatter(cors, cor_dem, label='DEM')
    plt.xlabel('Restitution coefficient')
    plt.ylabel('Computer simulation ratio')
    plt.legend()
    fig = os.path.join(directory, "cor_vs_cor_dem.pdf")
    plt.savefig(fig, dpi=300)

run()
