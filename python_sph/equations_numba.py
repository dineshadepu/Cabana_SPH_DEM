from math import pi, sqrt, exp
import numba as nb
import numpy as np
from numba import njit, prange, jit

M_1_PI = 1.0 / pi
dim = 2


@jit
def compute_grad_w(xij, rij, h, dim, grad):
    if dim == 1:
        fac = 1.0 / 120.0
    elif dim == 2:
        fac = M_1_PI * 7.0 / 478.0
    elif dim == 3:
        fac = M_1_PI * 1.0 / 120.0

    h1 = 1. / h
    q = rij * h1

    # get the kernel normalizing factor
    if dim == 1:
        fac = fac * h1
    elif dim == 2:
        fac = fac * h1 * h1
    elif dim == 3:
        fac = fac * h1 * h1 * h1

    tmp3 = 3. - q
    tmp2 = 2. - q
    tmp1 = 1. - q

    # compute the gradient
    if (rij > 1e-12):
        if (q > 3.0):
            val = 0.0

        elif (q > 2.0):
            val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3

        elif (q > 1.0):
            val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
            val += 30.0 * tmp2 * tmp2 * tmp2 * tmp2
        else:
            val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
            val += 30.0 * tmp2 * tmp2 * tmp2 * tmp2
            val -= 75.0 * tmp1 * tmp1 * tmp1 * tmp1
    else:
        val = 0.0

    wdash = val * fac

    h1 = 1. / h

    # compute the gradient.
    if (rij > 1e-12):
        tmp = wdash * h1 / rij
    else:
        tmp = 0.0

    grad[0] = tmp * xij[0]
    grad[1] = tmp * xij[1]
    grad[2] = tmp * xij[2]


@jit
def compute_w_kernel(xij=[0., 0, 0], rij=1.0, h=1.0, dim=2):
    if dim == 1:
        fac = 1.0 / 120.0
    elif dim == 2:
        fac = M_1_PI * 7.0 / 478.0
    elif dim == 3:
        fac = M_1_PI * 1.0 / 120.0
    h1 = 1. / h
    q = rij * h1

    # get the kernel normalizing factor
    if dim == 1:
        fac = fac * h1
    elif dim == 2:
        fac = fac * h1 * h1
    elif dim == 3:
        fac = fac * h1 * h1 * h1

    tmp3 = 3. - q
    tmp2 = 2. - q
    tmp1 = 1. - q
    if (q > 3.0):
        val = 0.0

    elif (q > 2.0):
        val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3

    elif (q > 1.0):
        val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
        val -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2

    else:
        val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
        val -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2
        val += 15. * tmp1 * tmp1 * tmp1 * tmp1 * tmp1

    return val * fac


def fluid_stage_1_numba(pa, dt):
    dtb2 = 0.5 * dt
    d_u = pa.u
    d_v = pa.v
    d_w = pa.w

    d_au = pa.au
    d_av = pa.av
    d_aw = pa.aw

    d_u[:] += dtb2 * d_au[:]
    d_v[:] += dtb2 * d_av[:]
    d_w[:] += dtb2 * d_aw[:]


def fluid_stage_2_numba(pa, dt):
    dtb2 = 0.5 * dt
    d_x = pa.x
    d_y = pa.y
    d_z = pa.z

    d_u = pa.u
    d_v = pa.v
    d_w = pa.w

    d_rho = pa.rho
    d_arho = pa.arho

    d_p = pa.p
    d_ap = pa.ap

    d_rho[:] += dt * d_arho[:]
    d_p[:] += dt * d_ap[:]

    d_x[:] += dt * d_u[:]
    d_y[:] += dt * d_v[:]
    d_z[:] += dt * d_w[:]


def fluid_stage_3_numba(pa, dt):
    dtb2 = 0.5 * dt
    d_u = pa.u
    d_v = pa.v
    d_w = pa.w

    d_au = pa.au
    d_av = pa.av
    d_aw = pa.aw

    d_u[:] += dtb2 * d_au[:]
    d_v[:] += dtb2 * d_av[:]
    d_w[:] += dtb2 * d_aw[:]


def make_arho_zero_numba(fluid):
    d_arho = fluid.arho

    d_arho[:] = 0.


@njit(parallel=True)
def continuity_equation(d_x, d_y, d_z, d_u, d_v, d_w, d_arho, d_h,
                        s_x, s_y, s_z, s_u, s_v, s_w, s_h, s_m):
    for i in prange(len(d_x)):
        for j in range(len(s_x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = np.array([0., 0., 0.])
            XIJ[0] = d_x[i] - s_x[j]
            XIJ[1] = d_y[i] - s_y[j]
            XIJ[2] = d_z[i] - s_z[j]

            VIJ = np.array([0., 0., 0.])
            VIJ = [0., 0., 0.]
            VIJ[0] = d_u[i] - s_u[j]
            VIJ[1] = d_v[i] - s_v[j]
            VIJ[2] = d_w[i] - s_w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5

            hij = (d_h[i] + s_h[j]) / 2.

            DWIJ = np.array([0., 0., 0.])
            compute_grad_w(XIJ, RIJ, hij, dim, DWIJ)
            # =========================
            # Precompute variables ends
            # =========================

            # continuity_equation stuff starts here
            vijdotdwij = DWIJ[0]*VIJ[0] + DWIJ[1]*VIJ[1] + DWIJ[2]*VIJ[2]
            d_arho[i] += s_m[j]*vijdotdwij


def continuity_equation_numba(fluid, tank):
    d_x = fluid.x
    d_y = fluid.y
    d_z = fluid.z
    d_u = fluid.u
    d_v = fluid.v
    d_w = fluid.w
    d_arho = fluid.arho
    d_h = fluid.h

    s_x = tank.x
    s_y = tank.y
    s_z = tank.z
    s_u = tank.u
    s_v = tank.v
    s_w = tank.w
    s_h = tank.h
    s_m = tank.m

    continuity_equation(d_x, d_y, d_z, d_u, d_v, d_w, d_arho, d_h,
                        s_x, s_y, s_z, s_u, s_v, s_w, s_h, s_m)


def state_equation_numba(fluid, p0, rho0, b):
    for i in range(len(fluid.x)):
        fluid.p[i] = p0 * (fluid.rho[i] / rho0 - b)


@njit(parallel=True)
def set_wall_pressure_bc(d_x, d_y, d_z, d_au, d_av, d_aw, d_wij, d_p,
                         d_h, s_x, s_y, s_z, s_p, s_rho,
                         gx, gy, gz):
    for i in range(len(d_x)):
        d_p[i] = 0.
        d_wij[i] = 0.

    for i in range(len(d_x)):
        for j in range(len(s_x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = np.array([0., 0., 0.])
            XIJ[0] = d_x[i] - s_x[j]
            XIJ[1] = d_y[i] - s_y[j]
            XIJ[2] = d_z[i] - s_z[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5

            hij = d_h[i]

            WIJ = compute_w_kernel(XIJ, RIJ, hij, dim)
            # =========================
            # Precompute variables ends
            # =========================

            gdotxij = (gx - d_au[i])*XIJ[0] + \
                (gy - d_av[i])*XIJ[1] + \
                (gz - d_aw[i])*XIJ[2]

            d_p[i] += s_p[j]*WIJ + s_rho[j]*gdotxij*WIJ
            d_wij[i] += WIJ

    for i in range(len(d_x)):
        if d_wij[i] > 1e-14:
            d_p[i] /= d_wij[i]


def set_wall_pressure_bc_numba(tank, fluid, gx, gy, gz):
    d_x = tank.x
    d_y = tank.y
    d_z = tank.z
    d_au = tank.u
    d_av = tank.v
    d_aw = tank.w
    d_h = tank.h
    d_wij = tank.wij
    d_p = tank.p

    s_x = fluid.x
    s_y = fluid.y
    s_z = fluid.z
    s_p = fluid.p
    s_rho = fluid.rho

    set_wall_pressure_bc(d_x, d_y, d_z, d_au, d_av, d_aw, d_wij, d_p,
                         d_h, s_x, s_y, s_z, s_p, s_rho,
                         gx, gy, gz)


def body_force_equation_numba(fluid, gx, gy, gz):
    fluid.au[:] = gx
    fluid.av[:] = gy
    fluid.aw[:] = gz


def momentum_equation(
        d_x, d_y, d_z, d_u, d_v, d_w, d_arho, d_h, d_p, d_rho,
        d_au, d_av, d_aw,
        s_x, s_y, s_z, s_u, s_v, s_w, s_h, s_m, s_p, s_rho):
    for i in range(len(d_x)):
        for j in range(len(s_x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = np.array([0., 0., 0.])
            XIJ[0] = d_x[i] - s_x[j]
            XIJ[1] = d_y[i] - s_y[j]
            XIJ[2] = d_z[i] - s_z[j]

            VIJ = np.array([0., 0., 0.])
            VIJ = [0., 0., 0.]
            VIJ[0] = d_u[i] - s_u[j]
            VIJ[1] = d_v[i] - s_v[j]
            VIJ[2] = d_w[i] - s_w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5

            hij = (d_h[i] + s_h[j]) / 2.

            DWIJ = np.array([0., 0., 0.])
            compute_grad_w(XIJ, RIJ, hij, dim, DWIJ)
            # =========================
            # Precompute variables ends
            # =========================

            rhoi2 = d_rho[i] * d_rho[i]
            rhoj2 = s_rho[j] * s_rho[j]

            pij = d_p[i]/rhoi2 + s_p[j]/rhoj2

            tmp = -s_m[j] * pij

            d_au[i] += tmp * DWIJ[0]
            d_av[i] += tmp * DWIJ[1]
            d_aw[i] += tmp * DWIJ[2]


def momentum_equation_numba(fluid, tank):
    d_x = fluid.x
    d_y = fluid.y
    d_z = fluid.z
    d_u = fluid.u
    d_v = fluid.v
    d_w = fluid.w
    d_arho = fluid.arho
    d_h = fluid.h
    d_au = fluid.au
    d_av = fluid.av
    d_aw = fluid.aw
    d_p = fluid.p
    d_rho = fluid.rho

    s_x = tank.x
    s_y = tank.y
    s_z = tank.z
    s_u = tank.u
    s_v = tank.v
    s_w = tank.w
    s_h = tank.h
    s_m = tank.m
    s_p = tank.p
    s_rho = tank.rho

    momentum_equation(
        d_x, d_y, d_z, d_u, d_v, d_w, d_arho, d_h, d_p, d_rho,
        d_au, d_av, d_aw,
        s_x, s_y, s_z, s_u, s_v, s_w, s_h, s_m, s_p, s_rho)


@njit(parallel=True)
def momentum_equation_artificial_viscosity(
        d_x, d_y, d_z, d_u, d_v, d_w, d_arho, d_h, d_p, d_rho,
        d_au, d_av, d_aw,
        s_x, s_y, s_z, s_u, s_v, s_w, s_h, s_m, s_p, s_rho,
        alpha, c0):
    for i in prange(len(d_x)):
        for j in range(len(s_x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = np.array([0., 0., 0.])
            XIJ[0] = d_x[i] - s_x[j]
            XIJ[1] = d_y[i] - s_y[j]
            XIJ[2] = d_z[i] - s_z[j]

            VIJ = np.array([0., 0., 0.])
            VIJ = [0., 0., 0.]
            VIJ[0] = d_u[i] - s_u[j]
            VIJ[1] = d_v[i] - s_v[j]
            VIJ[2] = d_w[i] - s_w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5
            R2IJ = RIJ**2.

            hij = (d_h[i] + s_h[j]) / 2.

            DWIJ = np.array([0., 0., 0.])
            compute_grad_w(XIJ, RIJ, hij, dim, DWIJ)
            # =========================
            # Precompute variables ends
            # =========================

            # v_{ab} \cdot r_{ab}
            vijdotrij = VIJ[0]*XIJ[0] + VIJ[1]*XIJ[1] + VIJ[2]*XIJ[2]

            # scalar part of the accelerations Eq. (11)
            piij = 0.0
            if vijdotrij < 0:
                HIJ = (d_h[i] + s_h[j]) / 2.
                EPS = 0.01 * HIJ**2.
                RHOIJ1 = 1. / (d_rho[i] + s_rho[j])
                muij = (HIJ * vijdotrij)/(R2IJ + EPS)

                piij = -alpha*c0*muij
                piij = s_m[j] * piij*RHOIJ1

            d_au[i] += -piij * DWIJ[0]
            d_av[i] += -piij * DWIJ[1]
            d_aw[i] += -piij * DWIJ[2]


def momentum_equation_artificial_viscosity_numba(fluid, tank, alpha, c0):
    d_x = fluid.x
    d_y = fluid.y
    d_z = fluid.z
    d_u = fluid.u
    d_v = fluid.v
    d_w = fluid.w
    d_arho = fluid.arho
    d_h = fluid.h
    d_au = fluid.au
    d_av = fluid.av
    d_aw = fluid.aw
    d_p = fluid.p
    d_rho = fluid.rho

    s_x = tank.x
    s_y = tank.y
    s_z = tank.z
    s_u = tank.u
    s_v = tank.v
    s_w = tank.w
    s_h = tank.h
    s_m = tank.m
    s_p = tank.p
    s_rho = tank.rho

    momentum_equation_artificial_viscosity(
        d_x, d_y, d_z, d_u, d_v, d_w, d_arho, d_h, d_p, d_rho,
        d_au, d_av, d_aw,
        s_x, s_y, s_z, s_u, s_v, s_w, s_h, s_m, s_p, s_rho, alpha, c0)
