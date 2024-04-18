import numpy as np
import matplotlib.pyplot as plt

# geometry imports
from sph_dem.geometry import (
    hydrostatic_tank_2d,
    translate_system_with_left_corner_as_origin,
    create_tank_2d_from_block_2d)

from sph_dem.fluids import (FluidsScheme,
                            get_particle_array_fluid,
                            get_particle_array_boundary)

from equations import (fluid_stage_1,
                       fluid_stage_2,
                       fluid_stage_3,
                       make_arho_zero,
                       continuity_equation,
                       state_equation,
                       set_wall_pressure_bc,
                       body_force_equation,
                       momentum_equation_artificial_viscosity,
                       momentum_equation)

from equations_numba import (fluid_stage_1_numba,
                             fluid_stage_2_numba,
                             fluid_stage_3_numba,
                             make_arho_zero_numba,
                             continuity_equation_numba,
                             state_equation_numba,
                             set_wall_pressure_bc_numba,
                             body_force_equation_numba,
                             momentum_equation_artificial_viscosity_numba,
                             momentum_equation_numba)

# GLOBAL VARIABLES
fluid_length = 1.
fluid_height = 0.5

tank_length = fluid_length
tank_height = 1.2 * fluid_height
tank_layers = 3

dx = 0.1
hdx = 1
h = hdx * dx
dim = 2
fluid_rho = 1000.
gx = 0.
gy = -9.81
gz = 0.
vref = np.sqrt(2. * abs(gy) * fluid_height)
c0 = 10 * vref
alpha = 0.1
p0 = fluid_rho * c0**2.
rho0 = fluid_rho
b = 1.

# numerical variable
dt = 1e-4
# tf = 100 * dt
tf = 1.
t = 0.
step = 0
print("total steps are", tf / dt)


def create_particles():
    xf, yf, xt, yt = hydrostatic_tank_2d(fluid_length, fluid_height,
                                         tank_height, tank_layers,
                                         dx, dx, False)
    if tank_length > fluid_length:
        xt, yt = create_tank_2d_from_block_2d(xf, yf,
                                              tank_length,
                                              tank_height,
                                              dx,
                                              tank_layers)

    zt = np.zeros_like(xt)
    zf = np.zeros_like(xf)

    # move fluid such that the left corner is at the origin of the
    # co-ordinate system
    translation = translate_system_with_left_corner_as_origin(xf, yf, zf)
    xt[:] = xt[:] - translation[0]
    yt[:] = yt[:] - translation[1]
    zt[:] = zt[:] - translation[2]

    # xf, yf, zf = np.array([0.02]), np.array([self.fluid_height]), np.array([0.])

    m = dx**dim * fluid_rho
    fluid = get_particle_array_fluid(name='fluid', x=xf, y=yf, z=zf, h=h, m=m, rho=fluid_rho)
    tank = get_particle_array_boundary(name='tank', x=xt, y=yt, z=zt, h=h, m=m, rho=fluid_rho)

    # set the pressure of the fluid
    fluid.p[:] = - fluid_rho * gy * (max(fluid.y) - fluid.y[:])
    # fluid.c0_ref[0] = self.c0
    # fluid.p0_ref[0] = self.p0
    return fluid, tank


# create particles
fluid, tank = create_particles()
# plt.scatter(fluid.x, fluid.y)
# plt.scatter(tank.x, tank.y)
# plt.show()

while t < tf:
    print(step)
    # # move velocities
    # fluid_stage_1(fluid, dt)

    # # compute acceleration 1
    # make_arho_zero(fluid)
    # continuity_equation(fluid, fluid)
    # continuity_equation(fluid, tank)

    # # move positions
    # fluid_stage_2(fluid, dt)

    # # compute acceleration 2
    # state_equation(fluid, p0, rho0, b)
    # set_wall_pressure_bc(tank, fluid, gx, gy, gz)
    # body_force_equation(fluid, gx, gy, gz)
    # momentum_equation(fluid, fluid)
    # momentum_equation(fluid, tank)
    # momentum_equation_artificial_viscosity(fluid, fluid, alpha, c0)
    # momentum_equation_artificial_viscosity(fluid, tank, alpha, c0)

    # # move velocities
    # fluid_stage_3(fluid, dt)

    # move velocities
    fluid_stage_1_numba(fluid, dt)

    # compute acceleration 1
    make_arho_zero_numba(fluid)
    continuity_equation_numba(fluid, fluid)
    continuity_equation_numba(fluid, tank)

    # move positions
    fluid_stage_2_numba(fluid, dt)

    # compute acceleration 2
    state_equation_numba(fluid, p0, rho0, b)
    set_wall_pressure_bc_numba(tank, fluid, gx, gy, gz)
    body_force_equation_numba(fluid, gx, gy, gz)
    momentum_equation_numba(fluid, fluid)
    momentum_equation_numba(fluid, tank)
    momentum_equation_artificial_viscosity_numba(fluid, fluid, alpha, c0)
    momentum_equation_artificial_viscosity_numba(fluid, tank, alpha, c0)

    # move velocities
    fluid_stage_3_numba(fluid, dt)

    t += dt
    if step % 10 == 0:
        plt.clf()
        plt.scatter(fluid.x, fluid.y)
        plt.scatter(tank.x, tank.y)
        # plt.show()
        name = "output/db_" + str(step) + ".png"
        plt.savefig(name)


    step += 1
