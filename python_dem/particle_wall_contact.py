import numpy as np
import matplotlib.pyplot as plt


class Particle():
    def __init__(self, x, y, radius, rho):
        self.x = x
        self.y = y
        self.u = 0.
        self.v = 0.
        # angular velocity
        self.wz = 0.
        self.fx = 0.
        self.fy = 0.
        self.tor = 0.
        self.rho = rho
        self.radius = radius
        self.m = np.pi * radius**2. * self.rho
        # moment of inertia
        self.moi = np.pi * radius**2. * self.rho


class Wall():
    def __init__(self, x, y, normal_x, normal_y):
        self.x = x
        self.y = y
        self.normal_x = normal_x
        self.normal_y = normal_y


def update_velocity(particle, dt):
    particle.u += particle.fx/particle.m * dt
    particle.v += particle.fy/particle.m * dt
    particle.wz += particle.tor/particle.moi * dt


def update_position(particle, dt):
    particle.x += particle.u * dt
    particle.y += particle.v * dt


def body_force(particle, gx, gy):
    # compute the force
    particle.fx += particle.m * gx
    particle.fy += particle.m * gy


def compute_force_on_particle_due_to_wall(particle, wall, kn, kt, fric_coeff, dt):
    # compute the overlap
    xij = particle.x - wall.x
    yij = particle.y - wall.y
    nij_x = wall.normal_x
    nij_y = wall.normal_y

    tmp = xij * nij_x  + yij * nij_y
    overlap = particle.radius - tmp

    if overlap > 0:
        # print("inside overlap")
        # compute the force
        # kn = kn
        fn = kn * overlap
        particle.fx = fn * nij_x
        particle.fy = fn * nij_y

        # tangential force
        # find the relative velocity
        vij_x = particle.u + nij_y * particle.wz * particle.radius - (0)
        vij_y = particle.v - nij_x * particle.wz * particle.radius - (0)
        vij_n = vij_x * nij_x + vij_y * nij_y

        # check if there is relative motion
        vij_magn = (vij_x**2. + vij_y**2.)**0.5

        # tangential unit vector
        _tij_x = vij_x - vij_n * nij_x
        _tij_y = vij_y - vij_n * nij_y
        _tij_magn = (_tij_x**2. + _tij_x**2.)**0.5
        if _tij_magn > 1e-12:
            tij_x = _tij_x / _tij_magn
            tij_y = _tij_y / _tij_magn
        else:
            tij_x = 0.
            tij_y = 0.

        if vij_magn < 1e-12:
            particle.ft_x = 0.
            particle.ft_y = 0.

            # reset the spring length
            particle.tng_disp_x = 0.
            particle.tng_disp_y = 0.

        else:
            # compute the tangential force, by equation 27
            tng_disp_x = particle.tng_disp_x + vij_x * dt
            tng_disp_y = particle.tng_disp_y + vij_y * dt
            tng_disp_magn = tng_disp_x * tij_x + tng_disp_y * tij_y
            particle.tng_disp_x = tng_disp_magn * tij_x
            particle.tng_disp_y = tng_disp_magn * tij_y

            ft_x_star = - kt * tng_disp_x
            ft_y_star = - kt * tng_disp_y

            # apply coulomb law
            ft_magn = (ft_x_star**2. + ft_y_star**2.)**0.5
            fn_magn = fn

            ft_magn_star = min(fric_coeff * fn_magn, ft_magn)
            # compute the tangential force, by equation 27
            particle.ft_x = -ft_magn_star * tij_x
            particle.ft_y = -ft_magn_star * tij_y
            particle.fx += particle.ft_x
            particle.fy += particle.ft_y

            # reset the spring length
            particle.tng_disp_x = -particle.ft_x / kt
            particle.tng_disp_y = -particle.ft_y / kt

        particle.tor = particle.ft_x * particle.radius

    else:
        fn = 0.
        particle.fx = 0.
        particle.fy = 0.

        particle.tng_disp_x = 0.
        particle.tng_disp_y = 0.

        particle.ft_x = 0.
        particle.ft_y = 0.
        particle.tor = 0.


def run():
    # create the particle
    particle = Particle(0.3, 0.2, 1e-1, 2000)
    particle.v = -10
    particle.u = 10
    wall = Wall(0., 0., 0., 1.)
    dt = 1e-4
    _t = 0
    tf = 0.04
    kn = 1e7
    kt = 1e5
    fric_coeff = 0.3
    gx = 0.
    gy = 0.
    t = []
    u = []
    y = []
    v = []
    fy = []
    fx = []
    wz = []
    while _t < tf:
        y.append(particle.y)
        u.append(particle.u)
        v.append(particle.v)
        fy.append(particle.fy)
        fx.append(particle.fx)
        wz.append(particle.wz)
        t.append(_t)

        update_velocity(particle, dt)
        compute_force_on_particle_due_to_wall(particle, wall, kn, kt, fric_coeff, dt)
        body_force(particle, gx, gy)
        update_position(particle, dt)
        # print(_t)
        _t += dt
        dt = 1e-4
    plt.plot(t, y)
    plt.savefig("y.jpg")
    plt.clf()
    plt.plot(t, v)
    plt.savefig("v.jpg")
    plt.clf()
    plt.plot(t, u)
    plt.savefig("u.jpg")
    plt.clf()
    plt.plot(t, fy)
    plt.savefig("fy.jpg")
    plt.clf()
    plt.plot(t, wz)
    plt.savefig("wz.jpg")
    plt.clf()
    plt.plot(t, fx)
    plt.savefig("fx.jpg")

run()
