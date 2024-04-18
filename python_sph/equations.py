from math import pi, sqrt, exp
import numba as nb
from numba import njit, prange

M_1_PI = 1.0 / pi
dim = 2


def compute_grad_w(xij=[0., 0., 0.], rij=1.0, h=1.0, dim=2, grad=[0, 0, 0]):
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


# @njit(parallel=True)
# def fluid_stage_1_numba(u, v, w, au, av, aw, dtb2):
#     for i in prange(int(u.shape[0])):
#         u[i] += dtb2 * au[i]
#         v[i] += dtb2 * av[i]
#         w[i] += dtb2 * aw[i]


def fluid_stage_1(pa, dt):
    dtb2 = 0.5 * dt
    u = pa.u
    v = pa.v
    w = pa.w

    au = pa.au
    av = pa.av
    aw = pa.aw

    for i in range(len(pa.x)):
        u[i] += dtb2 * au[i]
        v[i] += dtb2 * av[i]
        w[i] += dtb2 * aw[i]
    # fluid_stage_1_numba(u, v, w, au, av, aw, dtb2)


def fluid_stage_2(pa, dt):
    for i in range(len(pa.x)):
        pa.rho[i] += dt * pa.arho[i]
        pa.p[i] += dt * pa.ap[i]

        pa.x[i] += dt * pa.u[i]
        pa.y[i] += dt * pa.v[i]
        pa.z[i] += dt * pa.w[i]


def fluid_stage_3(pa, dt):
    dtb2 = 0.5 * dt
    for i in range(len(pa.x)):
        pa.u[i] += dtb2 * pa.au[i]
        pa.v[i] += dtb2 * pa.av[i]
        pa.w[i] += dtb2 * pa.aw[i]


def make_arho_zero(fluid):
    fluid.arho[:] = 0.


def continuity_equation(fluid, tank):
    for i in range(len(fluid.x)):
        for j in range(len(tank.x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = [0., 0., 0.]
            XIJ[0] = fluid.x[i] - tank.x[j]
            XIJ[1] = fluid.y[i] - tank.y[j]
            XIJ[2] = fluid.z[i] - tank.z[j]

            VIJ = [0., 0., 0.]
            VIJ[0] = fluid.u[i] - tank.u[j]
            VIJ[1] = fluid.v[i] - tank.v[j]
            VIJ[2] = fluid.w[i] - tank.w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5

            hij = (fluid.h[i] + tank.h[j]) / 2.

            DWIJ = [0., 0., 0.]
            compute_grad_w(XIJ, RIJ, hij, dim, DWIJ)
            # =========================
            # Precompute variables ends
            # =========================

            # continuity_equation stuff starts here
            vijdotdwij = DWIJ[0]*VIJ[0] + DWIJ[1]*VIJ[1] + DWIJ[2]*VIJ[2]
            fluid.arho[i] += tank.m[j]*vijdotdwij


def state_equation(fluid, p0, rho0, b):
    for i in range(len(fluid.x)):
        fluid.p[i] = p0 * (fluid.rho[i] / rho0 - b)


def set_wall_pressure_bc(tank, fluid, gx, gy, gz):
    tank.p[:] = 0.
    tank.wij[:] = 0.
    for i in range(len(tank.x)):
        for j in range(len(fluid.x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = [0., 0., 0.]
            XIJ[0] = tank.x[i] - fluid.x[j]
            XIJ[1] = tank.y[i] - fluid.y[j]
            XIJ[2] = tank.z[i] - fluid.z[j]

            VIJ = [0., 0., 0.]
            VIJ[0] = tank.u[i] - fluid.u[j]
            VIJ[1] = tank.v[i] - fluid.v[j]
            VIJ[2] = tank.w[i] - fluid.w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5

            hij = tank.h[i]

            DWIJ = [0., 0., 0.]
            compute_grad_w(XIJ, RIJ, hij, dim, DWIJ)
            WIJ = compute_w_kernel(XIJ, RIJ, hij, dim)
            # =========================
            # Precompute variables ends
            # =========================

            gdotxij = (gx - tank.au[i])*XIJ[0] + \
                (gy - tank.av[i])*XIJ[1] + \
                (gz - tank.aw[i])*XIJ[2]

            tank.p[i] += fluid.p[j]*WIJ + fluid.rho[j]*gdotxij*WIJ
            tank.wij[:] += WIJ

    for i in range(len(tank.x)):
        if tank.wij[i] > 1e-14:
            tank.p[i] /= tank.wij[i]


def body_force_equation(fluid, gx, gy, gz):
    for i in range(len(fluid.x)):
        fluid.au[i] = gx
        fluid.av[i] = gy
        fluid.aw[i] = gz


def momentum_equation(fluid, tank):
    for i in range(len(fluid.x)):
        for j in range(len(tank.x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = [0., 0., 0.]
            XIJ[0] = fluid.x[i] - tank.x[j]
            XIJ[1] = fluid.y[i] - tank.y[j]
            XIJ[2] = fluid.z[i] - tank.z[j]

            VIJ = [0., 0., 0.]
            VIJ[0] = fluid.u[i] - tank.u[j]
            VIJ[1] = fluid.v[i] - tank.v[j]
            VIJ[2] = fluid.w[i] - tank.w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5

            hij = fluid.h[i]

            DWIJ = [0., 0., 0.]
            compute_grad_w(XIJ, RIJ, hij, dim, DWIJ)
            # =========================
            # Precompute variables ends
            # =========================

            rhoi2 = fluid.rho[i] * fluid.rho[i]
            rhoj2 = tank.rho[j] * tank.rho[j]

            pij = fluid.p[i]/rhoi2 + tank.p[j]/rhoj2

            tmp = -tank.m[j] * pij

            fluid.au[i] += tmp * DWIJ[0]
            fluid.av[i] += tmp * DWIJ[1]
            fluid.aw[i] += tmp * DWIJ[2]


def momentum_equation_artificial_viscosity(fluid, tank, alpha, c0):
    for i in range(len(fluid.x)):
        for j in range(len(tank.x)):
            # =========================
            # Precompute variables
            # =========================
            XIJ = [0., 0., 0.]
            XIJ[0] = fluid.x[i] - tank.x[j]
            XIJ[1] = fluid.y[i] - tank.y[j]
            XIJ[2] = fluid.z[i] - tank.z[j]

            VIJ = [0., 0., 0.]
            VIJ[0] = fluid.u[i] - tank.u[j]
            VIJ[1] = fluid.v[i] - tank.v[j]
            VIJ[2] = fluid.w[i] - tank.w[j]

            RIJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)**0.5
            R2IJ = (XIJ[0]**2. + XIJ[1]**2. + XIJ[2]**2.)

            HIJ = fluid.h[i]

            DWIJ = [0., 0., 0.]
            compute_grad_w(XIJ, RIJ, HIJ, dim, DWIJ)
            # =========================
            # Precompute variables ends
            # =========================

            # v_{ab} \cdot r_{ab}
            vijdotrij = VIJ[0]*XIJ[0] + VIJ[1]*XIJ[1] + VIJ[2]*XIJ[2]

            # scalar part of the accelerations Eq. (11)
            piij = 0.0
            if vijdotrij < 0:
                EPS = 0.01 * HIJ**2.
                RHOIJ1 = 1. / (fluid.rho[i] + tank.rho[j])
                muij = (HIJ * vijdotrij)/(R2IJ + EPS)

                piij = -alpha*c0*muij
                piij = tank.m[j] * piij*RHOIJ1

            fluid.au[i] += -piij * DWIJ[0]
            fluid.av[i] += -piij * DWIJ[1]
            fluid.aw[i] += -piij * DWIJ[2]
