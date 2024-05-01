/****************************************************************************
 * Copyright (c) 2018-2022 by the Cabana authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the Cabana library. Cabana is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/
#include <Cabana_Core.hpp>
#include<math.h>

#include <iostream>


#define DIM 2
/*
  Start by declaring the types the particles will store. The first element
  will represent the coordinates, the second will be the particle's ID, the
  third velocity, and the fourth the radius of the particle.
*/

// rigid body data type
// cm, vcm, mass, force, indices limits
using RigidBodyDataType = Cabana::MemberTypes<
  int, // ids (0)
  double[3], // position(1)
  double[3], // velocity(2)
  double[3], // force(3)
  double, // mass (4)
  double // radius (5)
  >;


// auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
// auto rb_position = Cabana::slice<1>     ( rb,    "rb_position");
// auto rb_velocity = Cabana::slice<2>     ( rb,    "rb_velocity");
// auto rb_force = Cabana::slice<3>        ( rb,       "rb_force");
// auto rb_mass = Cabana::slice<4>        ( rb,        "rb_mass");
// auto rb_rad_s = Cabana::slice<5>        ( rb,        "rb_radius");

/*
  Next declare the data layout of the AoSoA. We use the host space here
  for the purposes of this example but all memory spaces, vector lengths,
  and member type configurations are compatible.
*/
const int VectorLength = 8;
using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecutionSpace::memory_space;

using DeviceType = Kokkos::Device<ExecutionSpace, MemorySpace>;
using RBAoSoA = Cabana::AoSoA<RigidBodyDataType, DeviceType, VectorLength>;

using ListAlgorithm = Cabana::FullNeighborTag;
using ListLayout = Cabana::VerletLayoutCSR;
using ListType = Cabana::VerletList<MemorySpace,ListAlgorithm,ListLayout>;


typedef Kokkos::View<double*>   ViewVectorType;
typedef Kokkos::View<double**>  ViewMatrixType;


void compute_force_on_dem_particles(RBAoSoA rb, double dt,
				    ListType * verlet_list, double fric_coeff){
  auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
  auto rb_position = Cabana::slice<1>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<2>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<3>        ( rb,       "rb_force");
  auto rb_mass = Cabana::slice<4>        ( rb,        "rb_mass");
  auto rb_rad_s = Cabana::slice<5>        ( rb,        "rb_radius");

  auto force_kernel = KOKKOS_LAMBDA( const int i, const int j )
    {
      const double mass_i = rb_mass( i );
      const double mass_j = rb_mass( j );
      const double radius_i = rb_rad_s( i );
      const double radius_j = rb_rad_s( j );

      const double x_i = rb_position( i, 0 );
      const double y_i = rb_position( i, 1 );
      const double z_i = rb_position( i, 2 );

      const double x_j = rb_position( j, 0 );
      const double y_j = rb_position( j, 1 );
      const double z_j = rb_position( j, 2 );

      const double xij = x_i - x_j;
      const double yij = y_i - y_j;
      const double zij = z_i - z_j;

      const double rsq = xij * xij + yij * yij + zij * zij;
      const double dist = sqrt(rsq);

      const double nij_x = xij / dist;
      const double nij_y = yij / dist;
      const double nij_z = zij / dist;

      const double overlap = radius_i + radius_j - dist;

      if (dist > 1e-12){

	// find the force if the particles are overlapping
	if (overlap > 0.) {
	  double kn = 1e6;
	  double fn_magn = kn * overlap;

	  double fn_x = fn_magn * nij_x;
	  double fn_y = fn_magn * nij_y;
	  double fn_z = fn_magn * nij_z;

	  // Delete these three once the tangential force is fixed
	  rb_force( i, 0 ) += fn_x;
	  rb_force( i, 1 ) += fn_y;
	  rb_force( i, 2 ) += fn_z;
	}

      }
    };

  Kokkos::RangePolicy<ExecutionSpace> policy(0, rb_mass.size());

  Cabana::neighbor_parallel_for( policy,
				 force_kernel,
				 *verlet_list,
				 Cabana::FirstNeighborsTag(),
				 Cabana::SerialOpTag(),
				 "CabanaDEM:Equations:ForceTorqueComputation" );
  Kokkos::fence();
}


void body_force_rigid_body(RBAoSoA & rb,
			   ViewVectorType & gravity,
			   double dt){
  auto rb_force = Cabana::slice<3>        ( rb,       "rb_force");
  auto rb_mass = Cabana::slice<4>        ( rb,        "rb_mass");

  auto half_dt = dt * 0.5;
  auto body_force_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );

      rb_force( i, 0 ) = gravity[0] * mass_i;
      rb_force( i, 1 ) = gravity[1] * mass_i;
      rb_force( i, 2 ) = gravity[2] * mass_i;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_mass.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBBodyForce", policy,
			body_force_lambda_func );
}

void euler_stage_1(RBAoSoA rb, double dt){
  auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
  auto rb_position = Cabana::slice<1>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<2>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<3>        ( rb,       "rb_force");
  auto rb_mass = Cabana::slice<4>        ( rb,        "rb_mass");
  auto rb_rad_s = Cabana::slice<5>        ( rb,        "rb_radius");

  auto half_dt = dt * 0.5;
  auto euler_stage_1_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );
      auto mass_i_1 = 1. / mass_i;

      rb_velocity( i, 0 ) += half_dt * rb_force( i, 0 ) * mass_i_1;
      rb_velocity( i, 1 ) += half_dt * rb_force( i, 1 ) * mass_i_1;
      rb_velocity( i, 2 ) += half_dt * rb_force( i, 2 ) * mass_i_1;

      rb_position( i, 0 ) += rb_velocity( i, 0 ) * dt;
      rb_position( i, 1 ) += rb_velocity( i, 1 ) * dt;
      rb_position( i, 2 ) += rb_velocity( i, 2 ) * dt;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:EulerStage1", policy,
			euler_stage_1_lambda_func );
}


void output_rb_data(RBAoSoA rb, int num_particles, int step, double time)
{
  // This is for setting HDF5 options
  auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
  auto rb_position = Cabana::slice<1>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<2>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<3>        ( rb,       "rb_force");
  auto rb_mass = Cabana::slice<4>        ( rb,        "rb_mass");
  auto rb_rad_s = Cabana::slice<5>        ( rb,        "rb_radius");

  Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
  Cabana::Experimental::HDF5ParticleOutput::
    writeTimeStep(
		  h5_config, "particles", MPI_COMM_WORLD,
		  step, time, num_particles,
		  rb_position,
		  rb_ids,
		  rb_velocity,
		  rb_force,
		  rb_mass);
}


//---------------------------------------------------------------------------//
// TODO: explain this function in short
//---------------------------------------------------------------------------//
void two_particles_colliding(const double t_final,
			     const int write_freq)
{
  int comm_rank = -1;
  MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

  if ( comm_rank == 0 )
    std::cout << "Cabana Rigid body solver example\n" << std::endl;


  /*
    ========================================================================
    Step 1: Define the dimentions and set the time step, time and pfreq
    ========================================================================
  */
  // Material properties of the fluid
  auto gx = 0.;
  auto gy = -9.81;
  auto gz = 0.;
  ViewVectorType gravity( "gravity", 3 );

  // Create host mirrors of device views.
  ViewVectorType::HostMirror h_gravity = Kokkos::create_mirror_view( gravity );

  h_gravity[0] = gx;
  h_gravity[1] = gy;
  h_gravity[2] = gz;

  // Deep copy host views to device views.
  Kokkos::deep_copy( gravity, h_gravity );
  // Numerical parameters of the SPH scheme for fluid

  // rigid body related properties
  auto rigid_body_rho = 2800.;
  auto rigid_body_radius = 0.01;
  auto mass_rb_i = 4/3. * M_PI * pow(rigid_body_radius, 3) * rigid_body_rho;
  auto _I = 2. / 5. * mass_rb_i * pow(rigid_body_radius, 2);
  auto I_inverse = 1. / _I;
  auto E_rb = 4.8 * 1e10;
  auto nu_rb = 0.2;
  auto u_rb = 10.;

  // integration related variables
  double h = 1. * rigid_body_radius;
  double dt = 1e-5;
  auto final_time = 0.002;
  auto time = 0.;
  int steps = final_time / dt;
  std::cout << "steps are :" << steps << std::endl;
  int print_freq = write_freq;
  /*
    ================================================
    End: Step 1
    ================================================
  */


  /*
    ================================================
    Step 2
    ================================================
    1. Create the particles
  */

  int no_rb_particles = 2;

  ViewVectorType x( "x", 2 );
  ViewVectorType y( "y", 2 );
  ViewVectorType z( "z", 2 );
  x[0] = 0.;
  x[1] = rigid_body_radius * 2 + rigid_body_radius / 10;

  y[0] = 0.;
  y[1] = 0.;

  z[0] = 0.;
  z[1] = 0.;

  /*
    ================================================
    End: Step 2
    ================================================
  */

  /*
    ================================================
    Step 3: Assign the particle properties to the aosoa
    ================================================
  */
  int total_no_particles = no_rb_particles;
  int num_particles = total_no_particles;
  std::cout << "Total no of particles: " << num_particles << std::endl;

  // std::vector<int> rigid_limits = {8, 16};
  RBAoSoA rb( "rb", total_no_particles );
  auto rb_position = Cabana::slice<1>( rb,    "rb_position");
  // create a mirror in the host space
  auto rb_host =
    Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), rb );

  auto rb_host_ids = Cabana::slice<0>          ( rb,         "rb_host_ids");
  auto rb_host_position = Cabana::slice<1>     ( rb,    "rb_host_position");
  auto rb_host_velocity = Cabana::slice<2>     ( rb,    "rb_host_velocity");
  auto rb_host_force = Cabana::slice<3>        ( rb,       "rb_host_force");
  auto rb_host_mass = Cabana::slice<4>        ( rb,        "rb_host_mass");
  auto rb_host_rad_s = Cabana::slice<5>        ( rb,        "rb_host_radius");

  for ( std::size_t i = 0; i < rb_host_position.size(); ++i )
    {

      std::cout << "Particle id i: " << i << std::endl;
      rb_host_ids ( i ) = i;

      rb_host_position ( i, 0 ) = x[i];
      rb_host_position ( i, 1 ) = y[i];
      rb_host_position ( i, 2 ) = z[i];

      for ( std::size_t j = 0; j < 3; ++j ) {
	rb_host_velocity ( i, j ) = 0.;

	rb_host_force ( i, j ) = 0.;
      }
      // This is only set explicitly for this example
      rb_host_velocity ( i, 0 ) = pow( -1, i) * u_rb;
      rb_host_mass ( i ) = mass_rb_i;
      rb_host_rad_s ( i ) = rigid_body_radius;
    }
  // copy it back to aosoa
  Cabana::deep_copy( rb, rb_host );

  output_rb_data(rb, num_particles, 0, time);
  // ================================================
  // ================================================
  // create the neighbor list
  // ================================================
  // ================================================

  double neighborhood_radius = 4. * rigid_body_radius;
  double grid_min[3] = { -10. * rigid_body_radius, -10. * rigid_body_radius, -neighborhood_radius - rigid_body_radius};
  double grid_max[3] = { 10. * rigid_body_radius, 10. * rigid_body_radius, neighborhood_radius + rigid_body_radius};
  // double grid_min[3] = { 0.0, -4.0, 0.0 };
  // double grid_max[3] = { 4.1, 4.0, 0.0 };
  double cell_ratio = 1.0;

  ListType verlet_list( rb_position, 0,
			rb_position.size(), neighborhood_radius,
			cell_ratio, grid_min, grid_max );

  // Main timestep loop
  // steps = 0;
  for ( int step = 0; step < steps; step++ )
    {
      verlet_list.build( rb_position, 0, rb_position.size(), neighborhood_radius,
			 cell_ratio, grid_min, grid_max );

      body_force_rigid_body(rb, gravity, dt);
      compute_force_on_dem_particles(rb, dt,
				     &verlet_list, 0.1);

      euler_stage_1(rb, dt);

      if ( step % print_freq == 0 )
	{

	  std::cout << "Time is:" << time << std::endl;
	  output_rb_data(rb, total_no_particles, step, time);
	}

      time += dt;

    }
}

int main( int argc, char* argv[] )
{

  MPI_Init( &argc, &argv );
  Kokkos::initialize( argc, argv );

  // check inputs and write usage
  if ( argc < 2 )
    {
      std::cerr << "Usage: ./DB2d body_spacing body_height "
	"body_length tank_height tank_length t_final write_freq\n";
      std::cerr << "      t_final           simulation end time\n";
      std::cerr
	<< "      write_freq      number of steps between output files\n";

      std::cerr << "\nfor example: ./NNPSBug 0.1 10\n";
      Kokkos::finalize();
      MPI_Finalize();
      return 0;
    }

  // end time.
  double t_final = std::atof( argv[1] );

  // write frequency
  int write_freq = std::atoi( argv[2] );

  // // device type
  // std::string device( argv[7] );

  // run the problem.
  two_particles_colliding(t_final, write_freq);

  Kokkos::finalize();

  MPI_Finalize();
  return 0;
}
