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

/*
  Start by declaring the types the particles will store. The first element
  will represent the coordinates, the second will be the particle's ID, the
  third velocity, and the fourth the radius of the particle.
*/

using DataTypes = Cabana::MemberTypes<double[3], // position(0) (x, y, z)
				      int, // ids(1)
				      double[3], // velocity(2) (u, v, w)
				      double[3], // acceleration(3) (au, av, aw)
				      double, // mass(4) (m)
				      double, // density(5) (rho)
				      double, // h (smoothing length) (6) (h)
				      double, // pressure (7) (p)
				      int, // is_fluid(8)
				      int, // is_boundary(9)
				      double, // rate of density change (10) (arho)
				      double, // rate of pressure change (11) (ap)
				      double, // sum_wij (12) (wij)
				      double[3], // velocity of wall (13)
				      double[3], // dummy velocity for computation (14)
				      int, // is_rb (15)
				      double, // rad_s (16)
				      int, // body_id (17)
				      double[3], // dx0, body frame local position (18)
				      double[3], // frc_dem, force on to rb particles (19)
				      double, // rb_m (20)
				      double // rb_rho (21)
				      >;

// rigid body data type
// cm, vcm, mass, force, indices limits
using RigidBodyDataType = Cabana::MemberTypes<int, // ids (0)
                                              int[2], // limits(1)
                                              double[3], // position(2)
                                              double[3], // velocity(3)
                                              double[3], // force(4)
                                              double[3], // torque(5)
                                              double[3], // linear acceleration (6)
                                              double[3], // angular acceleration (7)
                                              double[3], // angular momentum (8)
					      double[3], // angular velocity (9)
					      double[9], // Rotation matrix (10)
					      double, // mass (11)
					      double, // density (12)
					      double[9], // moment of inertia body frame (13)
					      double[9], // moment of inertia inverse body frame (14)
					      double[9], // moment of inertia global frame (15)
					      double[9], // moment of inertia inverse global frame (16)
					      double, // Rotation angle (17)
					      double // I_zz moi for 2d bodies (18)
					      >;

// auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
// auto aosoa_ids = Cabana::slice<1>          ( aosoa,    "ids");
// auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
// auto aosoa_acc = Cabana::slice<3>          ( aosoa,    "acc");
// auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
// auto aosoa_density = Cabana::slice<5>     ( aosoa,    "density");
// auto aosoa_h = Cabana::slice<6>           ( aosoa,    "h");
// auto aosoa_p = Cabana::slice<7>           ( aosoa,    "p");
// auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
// auto aosoa_is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");
// auto aosoa_density_acc = Cabana::slice<10>           ( aosoa,    "density_acc");
// auto aosoa_p_acc = Cabana::slice<11>           ( aosoa,    "p_acc");
// auto aosoa_sum_wij = Cabana::slice<12>           ( aosoa,    "sum_wij");
// auto aosoa_velocity_g = Cabana::slice<13>     ( aosoa,    "velocity_g");
// auto aosoa_velocity_f = Cabana::slice<14>     ( aosoa,    "velocity_f");
// auto aosoa_is_rb = Cabana::slice<15>       ( aosoa,    "is_rb");
// auto aosoa_rad_s = Cabana::slice<16>      ( aosoa,    "rad_s");
// auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
// auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");
// auto aosoa_frc_dem = Cabana::slice<19>     ( aosoa,    "frc_dem");
// auto aosoa_rb_m = Cabana::slice<20>     ( aosoa,    "aosoa_rb_m");
// auto aosoa_rb_rho = Cabana::slice<21>     ( aosoa,    "aosoa_rb_rho");


// auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
// auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
// auto rb_position = Cabana::slice<2>     ( rb,    "rb_position");
// auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
// auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
// auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
// auto rb_lin_acc = Cabana::slice<6>      ( rb,     "rb_lin_acc");
// auto rb_ang_acc = Cabana::slice<7>      ( rb,     "rb_ang_acc");
// auto rb_ang_mom = Cabana::slice<8>      ( rb,     "rb_ang_mom");
// auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
// auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");
// auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
// auto rb_density = Cabana::slice<12>     ( rb,     "rb_density");
// auto rb_body_moi = Cabana::slice<13>    ( rb,    "rb_body_moi");
// auto rb_inv_body_moi = Cabana::slice<14>( rb,    "rb_inv_body_moi");
// auto rb_global_moi = Cabana::slice<15>  ( rb,    "rb_global_moi");
// auto rb_inv_global_moi = Cabana::slice<16>( rb,    "rb_inv_global_moi");
// auto rb_rotation_angle = Cabana::slice<17>  ( rb,    "rb_rotation_angle");
// auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");

/*
  Next declare the data layout of the AoSoA. We use the host space here
  for the purposes of this example but all memory spaces, vector lengths,
  and member type configurations are compatible.
*/
const int VectorLength = 8;
using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = typename ExecutionSpace::memory_space;

using DeviceType = Kokkos::Device<ExecutionSpace, MemorySpace>;
using AoSoAType = Cabana::AoSoA<DataTypes, DeviceType, VectorLength>;
using RBAoSoA = Cabana::AoSoA<RigidBodyDataType, DeviceType, VectorLength>;


using ListAlgorithm = Cabana::FullNeighborTag;
using ListType =
  Cabana::VerletList<MemorySpace, ListAlgorithm, Cabana::VerletLayout2D>;
// using ListType =
//   Cabana::VerletList<MemorySpace, ListAlgorithm, Cabana::VerletLayoutCSR,
// 		     Cabana::TeamOpTag>;


typedef Kokkos::View<double*>   ViewVectorType;
typedef Kokkos::View<double**>  ViewMatrixType;


std::tuple<std::vector<double>,
	   std::vector<double>,
	   std::vector<double>> create_2d_block(double length, double height, double spacing){
  std::vector<double> x_vec;
  std::vector<double> y_vec;
  std::vector<double> z_vec;
  int x_no_points = length / spacing;
  int y_no_points = height / spacing;
  for(int i=0; i<=x_no_points; i++)
    {
      double x = i * spacing;
      for(int j=0; j<=y_no_points; j++)
	{
	  double y = j * spacing;
	  x_vec.push_back(x);
	  y_vec.push_back(y);
	  z_vec.push_back(0.);
	}
    }
  return std::make_tuple(x_vec, y_vec, z_vec);
}


// void set_body_limits_of_rb(AoSoATypeHost aosoa, RBAoSoAHost rb){
//   auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
//   auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");

//   /* =================================================================
//      =================================================================
//      ================================================================= */
//   /* start: compute center of mass and total mass of the rigid body */
//   /* =================================================================
//      =================================================================
//      ================================================================= */

//   for ( std::size_t i = 0; i < rb_limits.size(); ++i )
//     {
//       bool f_s_bool = false;
//       int found_start = 0;
//       int found_end = 0;

//       for ( std::size_t j = 0; j < aosoa_body_id.size(); ++j ) {
// 	if (aosoa_body_id(j) == i) {
// 	  if (f_s_bool == false) {
// 	    found_start = j;
// 	    f_s_bool = true;
// 	  }
// 	  else{
// 	    found_end = j;
// 	  }
// 	}
//       }
//       rb_limits( i, 0 ) = found_start;
//       rb_limits( i, 1 ) = found_end + 1;
//     }
// }



void setup_rigid_body_properties(AoSoAType aosoa, RBAoSoA rb, int * index_limits){
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
  auto aosoa_ids = Cabana::slice<1>          ( aosoa,    "ids");
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_acc = Cabana::slice<3>          ( aosoa,    "acc");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_density = Cabana::slice<5>     ( aosoa,    "density");
  auto aosoa_h = Cabana::slice<6>           ( aosoa,    "h");
  auto aosoa_p = Cabana::slice<7>           ( aosoa,    "p");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto aosoa_is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");
  auto aosoa_density_acc = Cabana::slice<10>           ( aosoa,    "density_acc");
  auto aosoa_p_acc = Cabana::slice<11>           ( aosoa,    "p_acc");
  auto aosoa_sum_wij = Cabana::slice<12>           ( aosoa,    "sum_wij");
  auto aosoa_velocity_g = Cabana::slice<13>     ( aosoa,    "velocity_g");
  auto aosoa_velocity_f = Cabana::slice<14>     ( aosoa,    "velocity_f");
  auto aosoa_is_rb = Cabana::slice<15>       ( aosoa,    "is_rb");
  auto aosoa_rad_s = Cabana::slice<16>      ( aosoa,    "rad_s");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");
  auto aosoa_frc_dem = Cabana::slice<19>     ( aosoa,    "frc_dem");
  auto aosoa_rb_m = Cabana::slice<20>     ( aosoa,    "aosoa_rb_m");
  auto aosoa_rb_rho = Cabana::slice<21>     ( aosoa,    "aosoa_rb_rho");


  auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_position = Cabana::slice<2>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_lin_acc = Cabana::slice<6>      ( rb,     "rb_lin_acc");
  auto rb_ang_acc = Cabana::slice<7>      ( rb,     "rb_ang_acc");
  auto rb_ang_mom = Cabana::slice<8>      ( rb,     "rb_ang_mom");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_density = Cabana::slice<12>     ( rb,     "rb_density");
  auto rb_body_moi = Cabana::slice<13>    ( rb,    "rb_body_moi");
  auto rb_inv_body_moi = Cabana::slice<14>( rb,    "rb_inv_body_moi");
  auto rb_global_moi = Cabana::slice<15>  ( rb,    "rb_global_moi");
  auto rb_inv_global_moi = Cabana::slice<16>( rb,    "rb_inv_global_moi");
  auto rb_rotation_angle = Cabana::slice<17>  ( rb,    "rb_rotation_angle");
  auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");

  /* =================================================================
     =================================================================
     ================================================================= */
  /* start: compute center of mass and total mass of the rigid body */
  /* =================================================================
     =================================================================
     ================================================================= */
  // auto com_total_mass_func = KOKKOS_LAMBDA( const int i, double& total_m, double& cm_x, double& cm_y, double& cm_z)
  //   {
  //     auto m_i = aosoa_mass( i );
  //     total_m += m_i;
  //     cm_x += m_i * aosoa_position(i, 0);
  //     cm_y += m_i * aosoa_position(i, 1);
  //     cm_z += m_i * aosoa_position(i, 2);
  //   };

  // for ( std::size_t i = 0; i < rb_position.size(); ++i )
  //   {
  //     double total_mass = 0.0;
  //     double cm_x = 0.0;
  //     double cm_y = 0.0;
  //     double cm_z = 0.0;

  //     std::cout << "rb limits are: " << rb_limits(i, 0) << " " << rb_limits(i, 1) << "\n";

  //     Kokkos::RangePolicy<ExecutionSpace> policy( rb_limits(i, 0), rb_limits(i, 1) );
  //     Kokkos::parallel_reduce(
  //                             "COM_and_total_mass_computation", policy,
  //                             com_total_mass_func, total_mass, cm_x, cm_y, cm_z);
  //     rb_mass(i) = total_mass;
  //     rb_position(i, 0) = cm_x / total_mass;
  //     rb_position(i, 1) = cm_y / total_mass;
  //     rb_position(i, 2) = cm_z / total_mass;
  //   }

  auto parallel_total_mass_func = KOKKOS_LAMBDA( const int i )
    {
      rb_mass(i) = 0.;
      rb_position(i, 0) = 0.;
      rb_position(i, 1) = 0.;
      rb_position(i, 2) = 0.;

      // The limits variable is used here. So the last index is not considerd
      for ( std::size_t j = rb_limits(i, 0); j < rb_limits(i, 1); ++j )
	{
	  auto m_j = aosoa_mass( j );
	  rb_mass(i) += m_j;
	  rb_position(i, 0) += m_j * aosoa_position(j, 0);
	  rb_position(i, 1) += m_j * aosoa_position(j, 1);
	  rb_position(i, 2) += m_j * aosoa_position(j, 2);
	}
      rb_position(i, 0) /= rb_mass(i);
      rb_position(i, 1) /= rb_mass(i);
      rb_position(i, 2) /= rb_mass(i);

    };
  Kokkos::RangePolicy<ExecutionSpace> policy_tm( 0, rb_mass.size());
  Kokkos::parallel_for( "CabanaRB:RBSetup:TotalMass", policy_tm,
                        parallel_total_mass_func );

  // for ( std::size_t i = 0; i < rb_position.size(); ++i )
  //   {
  //     std::cout << "\n";
  //     std::cout << "\n";
  //     std::cout << "rb mass: " << rb_mass(i) << "\n";
  //   }
  /* =================================================================
     end: compute center of mass and total mass of the rigid body
     ================================================================= */

  /* =================================================================
     start: compute moment of inertia
     ================================================================= */
  auto parallel_moi = KOKKOS_LAMBDA( const int i )
    {
      rb_I_zz(i) = 0.;
      // rb_position(i, 0) = 0.;
      // rb_position(i, 1) = 0.;
      // rb_position(i, 2) = 0.;

      // The limits variable is used here. So the last index is not considerd
      for ( std::size_t j = rb_limits(i, 0); j < rb_limits(i, 1); ++j )
	{
	  auto m_j = aosoa_mass( j );
	  double tmp_1 = pow((aosoa_position(j, 0) - rb_position(i, 0)), 2.);
	  double tmp_2 = pow((aosoa_position(j, 1) - rb_position(i, 1)), 2.);
	  rb_I_zz(i) += m_j * (tmp_1 + tmp_2);
	}
    };
  Kokkos::parallel_for( "CabanaRB:RBSetup:MOI", policy_tm,
                        parallel_moi );

  /* =================================================================
     end: compute center of mass and total mass of the rigid body
     ================================================================= */

  /* =================================================================
     start: save the initial positions of the particles in body frame axis
     ================================================================= */
  auto rb_save_initial_positions_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);
      aosoa_dx0(i, 0) = aosoa_position(i, 0) - rb_position( particle_body_id, 0 );
      aosoa_dx0(i, 1) = aosoa_position(i, 1) - rb_position( particle_body_id, 1 );
      aosoa_dx0(i, 2) = aosoa_position(i, 2) - rb_position( particle_body_id, 2 );
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( index_limits[0], index_limits[1] );
  Kokkos::parallel_for( "CabanaRB:RBSetup:SaveInitPos", policy,
                        rb_save_initial_positions_lambda_func );
  /* =================================================================
     end: save the initial positions of the particles in body frame axis
     ================================================================= */

  /* =================================================================
     start: save the initial positions of the particles in body frame axis
     ================================================================= */
  auto setup_all_properties_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      rb_ids ( i ) = i;
      rb_velocity ( i, 0 ) = 0.;
      rb_velocity ( i, 1 ) = 0.;
      rb_velocity ( i, 2 ) = 0.;

      rb_force ( i, 0 ) = 0.;
      rb_force ( i, 1 ) = 0.;
      rb_force ( i, 2 ) = 0.;

      rb_torque ( i, 0 ) = 0.;
      rb_torque ( i, 1 ) = 0.;
      rb_torque ( i, 2 ) = 0.;

      rb_lin_acc ( i, 0 ) = 0.;
      rb_lin_acc ( i, 1 ) = 0.;
      rb_lin_acc ( i, 2 ) = 0.;

      rb_ang_acc ( i, 0 ) = 0.;
      rb_ang_acc ( i, 1 ) = 0.;
      rb_ang_acc ( i, 2 ) = 0.;

      rb_ang_mom ( i, 0 ) = 0.;
      rb_ang_mom ( i, 1 ) = 0.;
      rb_ang_mom ( i, 2 ) = 0.;

      rb_ang_vel ( i, 0 ) = 0.;
      rb_ang_vel ( i, 1 ) = 0.;
      rb_ang_vel ( i, 2 ) = 0.;

      rb_rotation_matrix( i, 0 ) = 1.;
      rb_rotation_matrix( i, 1 ) = 0.;
      rb_rotation_matrix( i, 2 ) = 0.;

      rb_rotation_matrix( i, 3 ) = 0.;
      rb_rotation_matrix( i, 4 ) = 1.;
      rb_rotation_matrix( i, 5 ) = 0.;

      rb_rotation_matrix( i, 6 ) = 0.;
      rb_rotation_matrix( i, 7 ) = 0.;
      rb_rotation_matrix( i, 8 ) = 1.;

      rb_density ( i ) = 0.;

      rb_body_moi( i, 0 ) = 0.;
      rb_body_moi( i, 1 ) = 0.;
      rb_body_moi( i, 2 ) = 0.;

      rb_body_moi( i, 3 ) = 0.;
      rb_body_moi( i, 4 ) = 0.;
      rb_body_moi( i, 5 ) = 0.;

      rb_body_moi( i, 6 ) = 0.;
      rb_body_moi( i, 7 ) = 0.;
      rb_body_moi( i, 8 ) = 0.;

      rb_inv_body_moi( i, 0 ) = 0.;
      rb_inv_body_moi( i, 1 ) = 0.;
      rb_inv_body_moi( i, 2 ) = 0.;

      rb_inv_body_moi( i, 3 ) = 0.;
      rb_inv_body_moi( i, 4 ) = 0.;
      rb_inv_body_moi( i, 5 ) = 0.;

      rb_inv_body_moi( i, 6 ) = 0.;
      rb_inv_body_moi( i, 7 ) = 0.;
      rb_inv_body_moi( i, 8 ) = 0.;

      rb_global_moi( i, 0 ) = 0.;
      rb_global_moi( i, 1 ) = 0.;
      rb_global_moi( i, 2 ) = 0.;

      rb_global_moi( i, 3 ) = 0.;
      rb_global_moi( i, 4 ) = 0.;
      rb_global_moi( i, 5 ) = 0.;

      rb_global_moi( i, 6 ) = 0.;
      rb_global_moi( i, 7 ) = 0.;
      rb_global_moi( i, 8 ) = 0.;

      rb_inv_global_moi( i, 0 ) = 0.;
      rb_inv_global_moi( i, 1 ) = 0.;
      rb_inv_global_moi( i, 2 ) = 0.;

      rb_inv_global_moi( i, 3 ) = 0.;
      rb_inv_global_moi( i, 4 ) = 0.;
      rb_inv_global_moi( i, 5 ) = 0.;

      rb_inv_global_moi( i, 6 ) = 0.;
      rb_inv_global_moi( i, 7 ) = 0.;
      rb_inv_global_moi( i, 8 ) = 0.;

      rb_rotation_angle( i ) = 0.;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy_rb(0, rb_position.size() );
  Kokkos::parallel_for( "CabanaRB:RBSetup:SetupAllProperties", policy_rb,
                        setup_all_properties_lambda_func );
  /* =================================================================
     end: save the initial positions of the particles in body frame axis
     ================================================================= */
}

void set_linear_velocity_rigid_body(AoSoAType aosoa, RBAoSoA rb,
				    ViewVectorType & lin_vel,
				    int * index_limits){
  // get the properties of the particles and rigid body
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");

  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");

  // set the velocity of the particles of the rigid body
  auto rb_set_lin_vel_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      // set the linear velocity of the rigid body
      rb_velocity( i, 0) = lin_vel[3 * i + 0];
      rb_velocity( i, 1) = lin_vel[3 * i + 1];
      rb_velocity( i, 2) = lin_vel[3 * i + 2];
    };
  // range of particle indices for the current body no is
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_limits.size() );
  Kokkos::parallel_for( "CabanaRB:RB:SetAngVel", policy,
                        rb_set_lin_vel_func );

  // set the velocity of the particles of the rigid body
  auto aosoa_set_lin_vel_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      double dx = (rb_rotation_matrix( particle_body_id, 0 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 1 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 2 ) * aosoa_dx0( i, 3 ));
      double dy = (rb_rotation_matrix( particle_body_id, 3 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 4 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 5 ) * aosoa_dx0( i, 3 ));
      double dz = (rb_rotation_matrix( particle_body_id, 6 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 7 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 8 ) * aosoa_dx0( i, 3 ));

      double du = (rb_ang_vel( particle_body_id, 1 ) * dz -
		   rb_ang_vel( particle_body_id, 2 ) * dy);
      double dv = (rb_ang_vel( particle_body_id, 2 ) * dx -
		   rb_ang_vel( particle_body_id, 0 ) * dz);
      double dw = (rb_ang_vel( particle_body_id, 0 ) * dy -
		   rb_ang_vel( particle_body_id, 1 ) * dx);

      aosoa_velocity( i, 0 ) = rb_velocity( particle_body_id, 0 ) + du;
      aosoa_velocity( i, 1 ) = rb_velocity( particle_body_id, 1 ) + dv;
      aosoa_velocity( i, 2 ) = rb_velocity( particle_body_id, 2 ) + dw;

    };
  // range of particle indices for the current body no is
  // auto start = rb_limits(body_no, 0);
  // auto end = rb_limits(body_no, 1);
  Kokkos::RangePolicy<ExecutionSpace> policy_aosoa( index_limits[0], index_limits[1] );
  Kokkos::parallel_for( "CabanaRB:RB:AoSoASetLinVel", policy_aosoa,
                        aosoa_set_lin_vel_func );
}


void set_angular_velocity_rigid_body(AoSoAType aosoa, RBAoSoA rb,
				     ViewVectorType & ang_vel, int * index_limits){
  // get the properties of the particles and rigid body
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");

  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");

  // set the velocity of the particles of the rigid body
  auto rb_set_lin_vel_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      // set the linear velocity of the rigid body
      rb_ang_vel( i, 0) = ang_vel[3 * i + 0];
      rb_ang_vel( i, 1) = ang_vel[3 * i + 1];
      rb_ang_vel( i, 2) = ang_vel[3 * i + 2];
    };
  // range of particle indices for the current body no is
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_limits.size() );
  Kokkos::parallel_for( "CabanaRB:RB:SetAngVel", policy,
                        rb_set_lin_vel_func );

  // set the velocity of the particles of the rigid body
  auto aosoa_set_lin_vel_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      double dx = (rb_rotation_matrix( particle_body_id, 0 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 1 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 2 ) * aosoa_dx0( i, 3 ));
      double dy = (rb_rotation_matrix( particle_body_id, 3 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 4 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 5 ) * aosoa_dx0( i, 3 ));
      double dz = (rb_rotation_matrix( particle_body_id, 6 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 7 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 8 ) * aosoa_dx0( i, 3 ));

      double du = (rb_ang_vel( particle_body_id, 1 ) * dz -
		   rb_ang_vel( particle_body_id, 2 ) * dy);
      double dv = (rb_ang_vel( particle_body_id, 2 ) * dx -
		   rb_ang_vel( particle_body_id, 0 ) * dz);
      double dw = (rb_ang_vel( particle_body_id, 0 ) * dy -
		   rb_ang_vel( particle_body_id, 1 ) * dx);

      aosoa_velocity( i, 0 ) = rb_velocity( particle_body_id, 0 ) + du;
      aosoa_velocity( i, 1 ) = rb_velocity( particle_body_id, 1 ) + dv;
      aosoa_velocity( i, 2 ) = rb_velocity( particle_body_id, 2 ) + dw;

    };
  // range of particle indices for the current body no is
  // auto start = rb_limits(body_no, 0);
  // auto end = rb_limits(body_no, 1);
  Kokkos::RangePolicy<ExecutionSpace> policy_aosoa( index_limits[0], index_limits[1] );
  Kokkos::parallel_for( "CabanaRB:RB:AoSoASetLinVel", policy_aosoa,
                        aosoa_set_lin_vel_func );
}


void rigid_body_stage_1(RBAoSoA rb, double dt){
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_lin_acc = Cabana::slice<6>      ( rb,     "rb_lin_acc");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");

  auto half_dt = dt * 0.5;
  auto rb_stage_1_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );
      auto mass_i_1 = 1. / mass_i;
      rb_lin_acc( i, 0 ) = rb_force( i, 0 ) * mass_i_1;
      rb_lin_acc( i, 1 ) = rb_force( i, 1 ) * mass_i_1;
      rb_lin_acc( i, 2 ) = rb_force( i, 2 ) * mass_i_1;

      rb_velocity( i, 0 ) += rb_lin_acc( i, 0 ) * half_dt;
      rb_velocity( i, 1 ) += rb_lin_acc( i, 1 ) * half_dt;
      rb_velocity( i, 2 ) += rb_lin_acc( i, 2 ) * half_dt;

      rb_ang_vel( i, 2 ) += rb_torque( i, 2 ) / rb_I_zz( i ) * half_dt;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBStage1", policy,
                        rb_stage_1_lambda_func );
}


void rigid_body_particles_stage_1(AoSoAType aosoa, RBAoSoA rb, double dt, int * index_limits){
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");

  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");


  auto half_dt = dt * 0.5;
  auto rb_particles_stage_1_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      double dx = (rb_rotation_matrix( particle_body_id, 0 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 1 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 2 ) * aosoa_dx0( i, 2 ));
      double dy = (rb_rotation_matrix( particle_body_id, 3 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 4 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 5 ) * aosoa_dx0( i, 2 ));
      double dz = (rb_rotation_matrix( particle_body_id, 6 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 7 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 8 ) * aosoa_dx0( i, 2 ));

      double du = (rb_ang_vel( particle_body_id, 1 ) * dz -
		   rb_ang_vel( particle_body_id, 2 ) * dy);
      double dv = (rb_ang_vel( particle_body_id, 2 ) * dx -
		   rb_ang_vel( particle_body_id, 0 ) * dz);
      double dw = (rb_ang_vel( particle_body_id, 0 ) * dy -
		   rb_ang_vel( particle_body_id, 1 ) * dx);

      aosoa_velocity( i, 0 ) = rb_velocity( particle_body_id, 0 ) + du;
      aosoa_velocity( i, 1 ) = rb_velocity( particle_body_id, 1 ) + dv;
      aosoa_velocity( i, 2 ) = rb_velocity( particle_body_id, 2 ) + dw;
    };

  Kokkos::RangePolicy<ExecutionSpace> policy( index_limits[0], index_limits[1] );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBParticlesStage1", policy,
                        rb_particles_stage_1_lambda_func );
}


void rigid_body_stage_2(RBAoSoA rb, double dt){
  auto rb_position = Cabana::slice<2>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");
  auto rb_rotation_angle = Cabana::slice<17>  ( rb,    "rb_rotation_angle");

  auto rb_stage_2_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      rb_position( i, 0 ) += rb_velocity( i, 0 ) * dt;
      rb_position( i, 1 ) += rb_velocity( i, 1 ) * dt;
      rb_position( i, 2 ) += rb_velocity( i, 2 ) * dt;

      rb_rotation_angle( i ) += rb_ang_vel( i, 2 ) * dt;

      rb_rotation_matrix( i, 0 ) = cos( rb_rotation_angle( i ) );
      rb_rotation_matrix( i, 1 ) = - sin( rb_rotation_angle( i ) );
      rb_rotation_matrix( i, 2 ) = 0.;

      rb_rotation_matrix( i, 3 ) =  sin( rb_rotation_angle( i ) );
      rb_rotation_matrix( i, 4 ) =  cos( rb_rotation_angle( i ) );
      rb_rotation_matrix( i, 5 ) = 0.;

      rb_rotation_matrix( i, 6 ) = 0.;
      rb_rotation_matrix( i, 7 ) = 0.;
      rb_rotation_matrix( i, 8 ) = 1.;

    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBStage2", policy,
                        rb_stage_2_lambda_func );
}


void rigid_body_particles_stage_2(AoSoAType aosoa, RBAoSoA rb, double dt, int * index_limits){
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");

  auto rb_position = Cabana::slice<2>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");

  auto rb_particles_stage_2_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      double dx = (rb_rotation_matrix( particle_body_id, 0 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 1 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 2 ) * aosoa_dx0( i, 2 ));
      double dy = (rb_rotation_matrix( particle_body_id, 3 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 4 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 5 ) * aosoa_dx0( i, 2 ));
      double dz = (rb_rotation_matrix( particle_body_id, 6 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 7 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 8 ) * aosoa_dx0( i, 2 ));

      aosoa_position( i, 0 ) = rb_position( particle_body_id, 0 ) + dx;
      aosoa_position( i, 1 ) = rb_position( particle_body_id, 1 ) + dy;
      aosoa_position( i, 2 ) = rb_position( particle_body_id, 2 ) + dz;
    };

  Kokkos::RangePolicy<ExecutionSpace> policy( index_limits[0], index_limits[1] );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBParticlesStage2", policy,
                        rb_particles_stage_2_lambda_func );
}


void compute_force_on_rigid_body_particles(AoSoAType aosoa, double dt,
					   ListType * verlet_list_source,
					   int * index_limits){
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
  auto aosoa_ids = Cabana::slice<1>          ( aosoa,    "ids");
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_is_rb = Cabana::slice<15>       ( aosoa,    "is_rb");
  auto aosoa_rad_s = Cabana::slice<16>      ( aosoa,    "rad_s");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");
  auto aosoa_frc_dem = Cabana::slice<19>     ( aosoa,    "frc_dem");
  auto aosoa_rb_m = Cabana::slice<20>     ( aosoa,    "aosoa_rb_m");
  auto aosoa_rb_rho = Cabana::slice<21>     ( aosoa,    "aosoa_rb_rho");

  Cabana::deep_copy( aosoa_frc_dem, 0. );

  auto dem_contact_force_kernel = KOKKOS_LAMBDA( const int i, const int j )
    {
      if (aosoa_is_rb( j ) == 1 && aosoa_body_id( i ) != aosoa_body_id( j )){

	const double mass_i = aosoa_mass( i );
	const double mass_j = aosoa_mass( j );
	const double radius_i = aosoa_rad_s( i );
	const double radius_j = aosoa_rad_s( j );

	const double x_i = aosoa_position( i, 0 );
	const double y_i = aosoa_position( i, 1 );
	const double z_i = aosoa_position( i, 2 );

	const double x_j = aosoa_position( j, 0 );
	const double y_j = aosoa_position( j, 1 );
	const double z_j = aosoa_position( j, 2 );

	const double xij = x_i - x_j;
	const double yij = y_i - y_j;
	const double zij = z_i - z_j;

	const double rsq = xij * xij + yij * yij + zij * zij;
	const double dist = sqrt(rsq);

	const double nij_x = xij / dist;
	const double nij_y = yij / dist;
	const double nij_z = zij / dist;

	const double u_i = aosoa_velocity( i, 0 );
	const double v_i = aosoa_velocity( i, 1 );
	const double w_i = aosoa_velocity( i, 2 );

	const double u_j = aosoa_velocity( j, 0 );
	const double v_j = aosoa_velocity( j, 1 );
	const double w_j = aosoa_velocity( j, 2 );

	// const double wx_i = aosoa_ang_velocity( i, 0 );
	// const double wy_i = aosoa_ang_velocity( i, 1 );
	// const double wz_i = aosoa_ang_velocity( i, 2 );

	// const double wx_j = aosoa_ang_velocity( j, 0 );
	// const double wy_j = aosoa_ang_velocity( j, 1 );
	// const double wz_j = aosoa_ang_velocity( j, 2 );

	// const double u1 = u_i + (nij_y * wz_i - nij_z * wy_i) * radius_i;
	// const double u2 = -u_j + (nij_y * wz_j - nij_z * wy_j) * radius_j;
	// const double uij = u1 + u2;
	// const double v1 = v_i - (nij_x * wz_i - nij_z * wx_i) * radius_i;
	// const double v2 = -v_j - (nij_x * wz_j - nij_z * wx_j) * radius_j;
	// const double vij = v1 + v2;
	// const double w1 = w_i - (nij_x * wy_i - nij_y * wx_i) * radius_i;
	// const double w2 = -w_j - (nij_x * wy_j - nij_y * wx_j) * radius_j;
	// const double wij = w1 + w2;

	// const double vn = uij * nij_x + vij * nij_y + wij * nij_z;
	// const double vn_x = vn * nij_x;
	// const double vn_y = vn * nij_y;
	// const double vn_z = vn * nij_z;

	const double overlap = radius_i + radius_j - dist;

	if (overlap > 0){
	  /*
	    ############################
	    # normal force computation #
	    ############################
	  */
	  // // Compute stiffness
	  // // effective Young's modulus
	  // double tmp_1 = (1. - nu_i**2.) / E_i;
	  // double tmp_2 = (1. - nu_j**2.) / E_j;
	  // double E_eff = 1. / (tmp_1 + tmp_2);
	  // double tmp_1 = 1. / radius_i;
	  // double tmp_2 = 1. / radius_j;
	  // double R_eff = 1. / (tmp_1 + tmp_2);
	  // // Eq 4 [1];
	  // double kn = 4. / 3. * E_eff * R_eff**0.5;
	  double kn = 1e7;

	  // // compute damping coefficient
	  // double tmp_1 = log(cor);
	  // double tmp_2 = log(cor)**2. + pi**2.;
	  // double alpha_1 = -tmp_1 * (5. / tmp_2)**0.5;
	  // double tmp_1 = 1. / mass_i;
	  // double tmp_2 = 1. / mass_j;
	  // double m_eff = 1. / (tmp_1 + tmp_2);
	  // double eta = alpha_1 * (m_eff * kn)**0.5 * overlap**0.25;

	  double fn = kn * overlap * sqrt(overlap);
	  // double fn_x = fn * nij_x - eta * vn_x;
	  // double fn_y = fn * nij_y - eta * vn_y;
	  // double fn_z = fn * nij_z - eta * vn_z;
	  double fn_x = fn * nij_x;
	  double fn_y = fn * nij_y;
	  double fn_z = fn * nij_z;

	  // aosoa_frc_dem particle_1.fn = fn;
	  // particle_1.overlap = overlap;
	  aosoa_frc_dem (i, 0) += fn_x;
	  aosoa_frc_dem (i, 1) += fn_y;
	  aosoa_frc_dem (i, 2) += fn_z;
	}
      }
    };

  Kokkos::RangePolicy<ExecutionSpace> policy(index_limits[0], index_limits[1]);


  Cabana::neighbor_parallel_for( policy,
				 dem_contact_force_kernel,
				 *verlet_list_source,
				 Cabana::FirstNeighborsTag(),
				 Cabana::SerialOpTag(),
				 "dem_contact_force_loop" );
  Kokkos::fence();
}


void compute_effective_force_and_torque_on_rigid_body(RBAoSoA rb,  AoSoAType aosoa){
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_frc_dem = Cabana::slice<19>     ( aosoa,    "frc_dem");

  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_position = Cabana::slice<2>     ( rb,    "rb_position");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");

  auto total_force_torque_func = KOKKOS_LAMBDA( const int i )
    {
      // rb_force(i, 0) = 0.;
      // rb_force(i, 1) = 0.;
      // rb_force(i, 2) = 0.;
      auto mass_i = rb_mass( i );
      // rb_force( i, 0 ) = gravity[0] * mass_i;
      // rb_force( i, 1 ) = gravity[1] * mass_i;
      // rb_force( i, 2 ) = gravity[2] * mass_i;

      rb_torque(i, 0) = 0.;
      rb_torque(i, 1) = 0.;
      rb_torque(i, 2) = 0.;

      for ( std::size_t j = rb_limits(i, 0); j < rb_limits(i, 1); ++j )
	{
	  double fx_j = aosoa_frc_dem(j, 0);
	  double fy_j = aosoa_frc_dem(j, 1);
	  double fz_j = aosoa_frc_dem(j, 2);

	  rb_force(i, 0) += fx_j;
	  rb_force(i, 1) += fy_j;
	  rb_force(i, 2) += fz_j;

	  double dx = aosoa_position( j, 0 ) - rb_position( i, 0 );
	  double dy = aosoa_position( j, 1 ) - rb_position( i, 1 );
	  double dz = aosoa_position( j, 2 ) - rb_position( i, 2 );

	  rb_torque(i, 0) += dy * fz_j - dz * fy_j;
	  rb_torque(i, 1) += dz * fx_j - dx * fz_j;
	  rb_torque(i, 2) += dx * fy_j - dy * fx_j;
	}
    };

  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_mass.size());
  Kokkos::parallel_for( "CabanaRB:RB:TotalForceTorque", policy,
                        total_force_torque_func );
}


void body_force_rigid_body(RBAoSoA & rb,
			   ViewVectorType & gravity,
			   double dt){
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_lin_acc = Cabana::slice<6>      ( rb,     "rb_lin_acc");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");

  auto half_dt = dt * 0.5;
  auto body_force_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );

      rb_force( i, 0 ) = gravity[0] * mass_i;
      rb_force( i, 1 ) = gravity[1] * mass_i;
      rb_force( i, 2 ) = gravity[2] * mass_i;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBBodyForce", policy,
                        body_force_lambda_func );
}

void rigid_body_stage_3(RBAoSoA rb, double dt){
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_lin_acc = Cabana::slice<6>      ( rb,     "rb_lin_acc");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");

  auto half_dt = dt * 0.5;
  auto rb_stage_3_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );
      auto mass_i_1 = 1. / mass_i;
      rb_lin_acc( i, 0 ) = rb_force( i, 0 ) * mass_i_1;
      rb_lin_acc( i, 1 ) = rb_force( i, 1 ) * mass_i_1;
      rb_lin_acc( i, 2 ) = rb_force( i, 2 ) * mass_i_1;

      rb_velocity( i, 0 ) += rb_lin_acc( i, 0 ) * half_dt;
      rb_velocity( i, 1 ) += rb_lin_acc( i, 1 ) * half_dt;
      rb_velocity( i, 2 ) += rb_lin_acc( i, 2 ) * half_dt;

      rb_ang_vel( i, 2 ) += rb_torque( i, 2 ) / rb_I_zz( i ) * half_dt;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBStage3", policy,
                        rb_stage_3_lambda_func );
}


void rigid_body_particles_stage_3(AoSoAType aosoa, RBAoSoA rb, double dt, int * index_limits){
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto aosoa_dx0 = Cabana::slice<18>         ( aosoa,    "dx0");

  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");

  auto half_dt = dt * 0.5;
  auto rb_particles_stage_3_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto particle_body_id = aosoa_body_id(i);

      double dx = (rb_rotation_matrix( particle_body_id, 0 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 1 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 2 ) * aosoa_dx0( i, 2 ));
      double dy = (rb_rotation_matrix( particle_body_id, 3 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 4 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 5 ) * aosoa_dx0( i, 2 ));
      double dz = (rb_rotation_matrix( particle_body_id, 6 ) * aosoa_dx0( i, 0 ) +
		   rb_rotation_matrix( particle_body_id, 7 ) * aosoa_dx0( i, 1 ) +
		   rb_rotation_matrix( particle_body_id, 8 ) * aosoa_dx0( i, 2 ));

      double du = (rb_ang_vel( particle_body_id, 1 ) * dz -
		   rb_ang_vel( particle_body_id, 2 ) * dy);
      double dv = (rb_ang_vel( particle_body_id, 2 ) * dx -
		   rb_ang_vel( particle_body_id, 0 ) * dz);
      double dw = (rb_ang_vel( particle_body_id, 0 ) * dy -
		   rb_ang_vel( particle_body_id, 1 ) * dx);

      aosoa_velocity( i, 0 ) = rb_velocity( particle_body_id, 0 ) + du;
      aosoa_velocity( i, 1 ) = rb_velocity( particle_body_id, 1 ) + dv;
      aosoa_velocity( i, 2 ) = rb_velocity( particle_body_id, 2 ) + dw;
    };

  Kokkos::RangePolicy<ExecutionSpace> policy( index_limits[0], index_limits[1] );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBParticlesStage3", policy,
                        rb_particles_stage_3_lambda_func );
}


void output_data(AoSoAType aosoa, int num_particles, int step, double time)
{
  // This is for setting HDF5 options
  auto position = Cabana::slice<0>     ( aosoa,    "position");
  auto ids = Cabana::slice<1>          ( aosoa,    "ids");
  auto velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto acc = Cabana::slice<3>          ( aosoa,    "acc");
  auto mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto density = Cabana::slice<5>     ( aosoa,    "density");
  auto h = Cabana::slice<6>           ( aosoa,    "h");
  auto p = Cabana::slice<7>           ( aosoa,    "p");
  auto is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");
  auto density_acc = Cabana::slice<10>           ( aosoa,    "density_acc");
  auto p_acc = Cabana::slice<11>           ( aosoa,    "p_acc");
  auto sum_wij = Cabana::slice<12>           ( aosoa,    "sum_wij");
  auto velocity_g = Cabana::slice<13>     ( aosoa,    "velocity_g");
  auto velocity_f = Cabana::slice<14>     ( aosoa,    "velocity_f");
  auto is_rb = Cabana::slice<15>       ( aosoa,    "is_rb");
  auto rad_s = Cabana::slice<16>      ( aosoa,    "rad_s");
  auto body_id = Cabana::slice<17>     ( aosoa,    "body_id");
  auto dx0 = Cabana::slice<18>         ( aosoa,    "dx0");
  auto frc_dem = Cabana::slice<19>     ( aosoa,    "frc_dem");
  auto rb_m = Cabana::slice<20>     ( aosoa,    "aosoa_rb_m");
  auto rb_rho = Cabana::slice<21>     ( aosoa,    "aosoa_rb_rho");

  Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
  Cabana::Experimental::HDF5ParticleOutput::
    writeTimeStep(
		  h5_config, "particles", MPI_COMM_WORLD,
		  step, time, num_particles, position,
		  ids, velocity, density, rad_s, body_id, mass, frc_dem);
}


void output_rb_data(RBAoSoA rb, int num_particles, int step, double time)
{
  // This is for setting HDF5 options
  auto rb_ids = Cabana::slice<0>          ( rb,         "rb_ids");
  auto rb_position = Cabana::slice<2>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_rotation_angle = Cabana::slice<17>  ( rb,    "rb_rotation_angle");
  auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");

  Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
  Cabana::Experimental::HDF5ParticleOutput::
    writeTimeStep(
		  h5_config, "rigid_bodies", MPI_COMM_WORLD,
		  step, time, num_particles,
		  rb_position,
		  rb_ids,
		  rb_velocity,
		  rb_force,
		  rb_torque,
		  rb_ang_vel,
		  rb_mass,
		  rb_rotation_angle,
		  rb_I_zz);
}


//---------------------------------------------------------------------------//
// TODO: explain this function in short
//---------------------------------------------------------------------------//
void rb_freely_rotating_square_bodies(const double body_spacing,
				      const double body_height,
				      const double body_length,
				      const int no_of_bodies,
				      const double t_final,
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
  auto gy = 0.;
  auto gz = 0.;
  auto fluid_rho = 1000.;
  auto vref = sqrt(2. * abs(gy) * body_height);
  std::cout << "vref is :" << vref << std::endl;
  auto c0 = 10. * vref;
  // auto fluid_viscosity = 1e-4;
  auto alpha = 0.05;
  auto p0 = c0*c0 * fluid_rho;
  auto rho0 = fluid_rho;
  ViewVectorType b_rho0_p0( "b_rho0_p0", 3 );
  ViewVectorType gravity( "gravity", 3 );
  ViewVectorType alpha_c0( "alpha_c0", 2 );

  // Create host mirrors of device views.
  ViewVectorType::HostMirror h_b_rho0_p0 = Kokkos::create_mirror_view( b_rho0_p0 );
  ViewVectorType::HostMirror h_gravity = Kokkos::create_mirror_view( gravity );
  ViewVectorType::HostMirror h_alpha_c0 = Kokkos::create_mirror_view( alpha_c0 );

  h_b_rho0_p0[0] = 1.;
  h_b_rho0_p0[1] = fluid_rho;
  h_b_rho0_p0[2] = p0;

  h_gravity[0] = gx;
  h_gravity[1] = gy;
  h_gravity[2] = gz;

  h_alpha_c0[0] = alpha;
  h_alpha_c0[1] = c0;

  // Deep copy host views to device views.
  Kokkos::deep_copy( b_rho0_p0, h_b_rho0_p0 );
  Kokkos::deep_copy( gravity, h_gravity );
  Kokkos::deep_copy( alpha_c0, h_alpha_c0 );
  // Numerical parameters of the SPH scheme for fluid

  // rigid body related properties
  auto rigid_body_rho = 2000.;

  // integration related variables
  double h = 1. * body_spacing;
  double dt = 0.25 * h / (c0 + vref);
  dt = 1e-4;
  std::cout << "dt is :" << dt << std::endl;
  auto final_time = t_final;
  auto time = 0.;
  int steps = final_time / dt;
// int steps = 10000;
   // int steps = 1;
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
    1. Create the particles of the rigid bodies
  */
  // double fluid_spacing = fluid_spacing;
  // double fluid_height = fluid_height;
  // double fluid_length = fluid_length;
  // double tank_height = tank_height;
  // double tank_length = tank_length;
  // int tank_layers = 3;

  int no_rb_particles = 0;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<int> is_fluid;
  std::vector<int> is_bdry;
  std::vector<int> is_rb;
  std::vector<int> body_id;

  for ( std::size_t i = 0; i < no_of_bodies; ++i ) {

    std::vector<double> x_rb, y_rb, z_rb;
    std::tie(x_rb, y_rb, z_rb) = create_2d_block(body_length, body_height,
						 body_spacing);


    for ( std::size_t j = 0; j < x_rb.size(); ++j )
      {
	int fac = i % 10;
	int fac_1 = int (i / 4);
	x.push_back(x_rb[j] + fac * body_length + fac * 3. * body_spacing);

	// y.push_back(y_rb[j] + fac_1 * body_length + fac_1 * 8. * body_spacing);
	y.push_back(y_rb[j] + (i * 4. * body_spacing));

	z.push_back(z_rb[j]);

	no_rb_particles += 1;

	is_fluid.push_back(0);
	is_bdry.push_back(0);
	is_rb.push_back(1);
	body_id.push_back(i);
      }
  }

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
  int no_fluid_particles = 0;
  int no_bdry_particles = 0;
  int total_no_particles = no_fluid_particles + no_bdry_particles + no_rb_particles;
  int num_particles = total_no_particles;

  AoSoAType aosoa( "particles", total_no_particles );
  // setup_aosoa(aosoa, &x, &y, &z, &u);
  auto aosoa_position = Cabana::slice<0>( aosoa,    "position");

  // create a mirror in the host space
  auto aosoa_host =
    Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), aosoa );

  auto aosoa_host_position = Cabana::slice<0>     ( aosoa_host,    "position");
  auto aosoa_host_ids = Cabana::slice<1>          ( aosoa_host,    "ids");
  auto aosoa_host_velocity = Cabana::slice<2>     ( aosoa_host,    "velocity");
  auto aosoa_host_acc = Cabana::slice<3>          ( aosoa_host,    "acc");
  auto aosoa_host_mass = Cabana::slice<4>        ( aosoa_host,    "mass");
  auto aosoa_host_density = Cabana::slice<5>     ( aosoa_host,    "density");
  auto aosoa_host_h = Cabana::slice<6>           ( aosoa_host,    "h");
  auto aosoa_host_p = Cabana::slice<7>           ( aosoa_host,    "p");
  auto aosoa_host_is_fluid = Cabana::slice<8>    ( aosoa_host,    "is_fluid");
  auto aosoa_host_is_boundary = Cabana::slice<9>    ( aosoa_host,    "is_boundary");
  auto aosoa_host_density_acc = Cabana::slice<10>           ( aosoa_host,    "density_acc");
  auto aosoa_host_p_acc = Cabana::slice<11>           ( aosoa_host,    "p_acc");
  auto aosoa_host_sum_wij = Cabana::slice<12>           ( aosoa_host,    "sum_wij");
  auto aosoa_host_velocity_g = Cabana::slice<13>     ( aosoa_host,    "velocity_g");
  auto aosoa_host_velocity_f = Cabana::slice<14>     ( aosoa_host,    "velocity_f");
  auto aosoa_host_is_rb = Cabana::slice<15>       ( aosoa_host,    "is_rb");
  auto aosoa_host_rad_s = Cabana::slice<16>      ( aosoa_host,    "rad_s");
  auto aosoa_host_body_id = Cabana::slice<17>     ( aosoa_host,    "body_id");
  auto aosoa_host_dx0 = Cabana::slice<18>         ( aosoa_host,    "dx0");
  auto aosoa_host_frc_dem = Cabana::slice<19>     ( aosoa_host,    "frc_dem");
  auto aosoa_host_rb_m = Cabana::slice<20>     ( aosoa_host,    "aosoa_host_rb_m");
  auto aosoa_host_rb_rho = Cabana::slice<21>     ( aosoa_host,    "aosoa_host_rb_rho");

  for ( std::size_t i = 0; i < aosoa_host_position.size(); ++i )
    {
      aosoa_host_position ( i, 0 ) = x[i];
      aosoa_host_position ( i, 1 ) = y[i];
      aosoa_host_position ( i, 2 ) = z[i];

      aosoa_host_ids ( i ) = i;

      aosoa_host_velocity ( i, 0 ) = 0.;
      aosoa_host_velocity ( i, 1 ) = 0.;
      aosoa_host_velocity ( i, 2 ) = 0.;

      aosoa_host_acc ( i, 0 ) = 0.;
      aosoa_host_acc ( i, 1 ) = 0.;
      aosoa_host_acc ( i, 2 ) = 0.;

      aosoa_host_mass ( i ) = fluid_rho * pow(body_spacing, DIM);
      aosoa_host_density ( i ) = fluid_rho;
      aosoa_host_h ( i ) = 1. * body_spacing;

      aosoa_host_p ( i ) = 0.;
      aosoa_host_is_fluid ( i ) = is_fluid[i];
      aosoa_host_is_boundary ( i ) = is_bdry[i];
      aosoa_host_density_acc ( i ) = 0.;
      aosoa_host_p_acc ( i ) = 0.;
      aosoa_host_sum_wij ( i ) = 0.;

      aosoa_host_velocity_g ( i, 0 ) = 0.;
      aosoa_host_velocity_g ( i, 1 ) = 0.;
      aosoa_host_velocity_g ( i, 2 ) = 0.;

      aosoa_host_velocity_f ( i, 0 ) = 0.;
      aosoa_host_velocity_f ( i, 1 ) = 0.;
      aosoa_host_velocity_f ( i, 2 ) = 0.;

      // Only rigid body corresponding variables
      aosoa_host_is_rb ( i ) = is_rb[i];

      aosoa_host_rad_s ( i ) = body_spacing * 0.5;

      aosoa_host_body_id ( i ) = body_id[i];

      aosoa_host_dx0 ( i, 0 ) = 0.;
      aosoa_host_dx0 ( i, 1 ) = 0.;
      aosoa_host_dx0 ( i, 2 ) = 0.;

      aosoa_host_frc_dem ( i, 0 ) = 0.;
      aosoa_host_frc_dem ( i, 1 ) = 0.;
      aosoa_host_frc_dem ( i, 2 ) = 0.;

      aosoa_host_rb_m ( i ) = rigid_body_rho * pow(body_spacing, DIM);
      aosoa_host_rb_rho ( i ) = rigid_body_rho;
    }
  // copy it back to aosoa
  Cabana::deep_copy( aosoa, aosoa_host );

  // int fluid_limits[2] = {0, no_fluid_particles};

  // std::cout << "Fluid limits are :" << fluid_limits[0] << ", " << fluid_limits[1] << std::endl;
  // int bdry_limits[2] = {no_fluid_particles, no_fluid_particles + no_bdry_particles};
  // std::cout << "boundary limits are :" << bdry_limits[0] << ", " << bdry_limits[1] << std::endl;

  int rigid_limits[2] = {0, no_rb_particles};
  // std::cout << "rigid limits are :" << rigid_limits[0] << ", " << rigid_limits[1] << std::endl;

  output_data(aosoa, num_particles, 0, time);
  // output_data(aosoa, num_particles, 100, time);

  /*
    ================================
    rigid body starts
    ================================
  */
  /*
    Create the rigid body data type.
  */
  // std::vector<int> rigid_limits = {8, 16};
  RBAoSoA rb( "rb", no_of_bodies );
  {
    /*
      ================================
      set body limits of rb
      ================================
    */
    auto aosoa_body_id = Cabana::slice<17>     ( aosoa,    "body_id");
    auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
    auto set_limits_func = KOKKOS_LAMBDA( const int i )
      {
	bool f_s_bool = false;
	int found_start = 0;
	int found_end = 0;

	for ( std::size_t j = 0; j < aosoa_body_id.size(); ++j ) {
	  if (aosoa_body_id(j) == i) {
	    if (f_s_bool == false) {
	      found_start = j;
	      f_s_bool = true;
	    }
	    else{
	      found_end = j;
	    }
	  }
	}
	rb_limits( i, 0 ) = found_start;
	rb_limits( i, 1 ) = found_end + 1;
      };
    Kokkos::RangePolicy<ExecutionSpace> policy_set_limits( 0, rb_limits.size());
    Kokkos::parallel_for( "CabanaRB:RBSetup:SetLimits", policy_set_limits,
			  set_limits_func );

  }

  // get the rigid body limits
  // set_body_limits_of_rb(aosoa_host, rb_host);
  // auto rb_limits = Cabana::slice<1>       ( rb_host,      "rb_limits");

  // for ( std::size_t j = 0; j < rb.size(); ++j )
  //   {
  //     std::cout << "rb_limits of body : " << j << " are: " << rb_limits(j, 0) << ", " << rb_limits(j, 1) << std::endl;

  //   }

  // setup the master rigid body properties
  setup_rigid_body_properties(aosoa, rb, rigid_limits);

  // create a mirror in the host space
  // auto rb_host =
  //   Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), rb );

  // Set the velocity of the bodies
  ViewVectorType lin_vel( "lin_vel", 3*no_of_bodies );
  ViewVectorType ang_vel( "ang_vel", 3*no_of_bodies );

  // Create host mirrors of device views.
  ViewVectorType::HostMirror h_lin_vel = Kokkos::create_mirror_view( lin_vel );
  ViewVectorType::HostMirror h_ang_vel = Kokkos::create_mirror_view( ang_vel );

  for ( std::size_t i = 0; i < rb.size(); ++i )
    {
      h_lin_vel[3 * i + 0] = pow(-1., i) * 2.;
      h_lin_vel[3 * i + 1] = 0.;
      h_lin_vel[3 * i + 2] = 0.;
      h_ang_vel[3 * i + 0] = 0.;
      h_ang_vel[3 * i + 1] = 0.;
      h_ang_vel[3 * i + 2] = 0.;
    }

  // Deep copy host views to device views.
  Kokkos::deep_copy( lin_vel, h_lin_vel );
  Kokkos::deep_copy( ang_vel, h_ang_vel );

  set_linear_velocity_rigid_body(aosoa, rb, lin_vel, rigid_limits);
  set_angular_velocity_rigid_body(aosoa, rb, ang_vel, rigid_limits);

  // copy it back to aosoa
  // Cabana::deep_copy( rb, rb_host );
  // Cabana::deep_copy( aosoa, aosoa_host );
  /*
    ================================================
    End: Step 3
    ================================================
    */


    // output the data
    // output_data(aosoa, num_particles, 0, 0.);
    // output_data(aosoa, num_particles, 0, 0.);

    // ================================================
    // ================================================
    // create the neighbor list
    // ================================================
    // ================================================
  double neighborhood_radius = 3. * body_spacing;
  double grid_min[3] = { -500.0, -500.0, -neighborhood_radius };
  double grid_max[3] = { 500.0, 500.0, 1. * neighborhood_radius };
  // double grid_min[3] = { 0.0, -4.0, 0.0 };
    // double grid_max[3] = { 4.1, 4.0, 0.0 };
    double cell_ratio = 1.0;

    ListType verlet_list( aosoa_position, 0,
			  aosoa_position.size(), neighborhood_radius,
			  cell_ratio, grid_min, grid_max );

    // Main timestep loop
    for ( int step = 0; step < steps; step++ )
    {
      rigid_body_stage_1(rb, dt);
      rigid_body_particles_stage_1(aosoa, rb, dt, rigid_limits);

      rigid_body_stage_2(rb, dt);
      rigid_body_particles_stage_2(aosoa, rb, dt, rigid_limits);

      // Compute the neighbours
      ListType verlet_list( aosoa_position, 0,
			    aosoa_position.size(), neighborhood_radius,
			    cell_ratio, grid_min, grid_max );

      compute_force_on_rigid_body_particles(aosoa, dt,
					    &verlet_list,
					    rigid_limits);

      body_force_rigid_body(rb, gravity, dt);
      compute_effective_force_and_torque_on_rigid_body(rb,
						       aosoa);

      rigid_body_stage_3(rb, dt);
      rigid_body_particles_stage_3(aosoa, rb, dt, rigid_limits);

      if ( step % print_freq == 0 )
	{
	  std::cout << "Time is:" << time << std::endl;
	  output_data(aosoa, num_particles, step, time);
	  output_rb_data(rb, no_of_bodies, step, time);
	}

      time += dt;

    }
}

int main( int argc, char* argv[] )
{

  MPI_Init( &argc, &argv );
  Kokkos::initialize( argc, argv );

  // check inputs and write usage
  if ( argc < 7 )
    {
      std::cerr << "Usage: ./DB2d body_spacing body_height "
	"body_length tank_height tank_length t_final write_freq\n";
      std::cerr << "\nwhere body_spacing       spacing between the particles"
	"\n";
      std::cerr
	<< "      body_height  length of the body block\n";
      std::cerr
	<< "      body_length  length of the body block\n";
      std::cerr
	<< "      no_of_bodies  no of rigid bodies \n";

      std::cerr << "      t_final           simulation end time\n";
      std::cerr
	<< "      write_freq      number of steps between output files\n";

      std::cerr << "\nfor example: ./01RBFreelyRotatingBodies 0.1 1.0 1.0 6 3.0 100\n";
      Kokkos::finalize();
      MPI_Finalize();
      return 0;
    }

  // cell size
  double body_spacing = std::atof( argv[1] );

  // body height
  double body_height = std::atof( argv[2] );

  // body length
  double body_length = std::atof( argv[3] );

  // body height
  double no_of_bodies = std::atoi( argv[4] );

  // end time.
  double t_final = std::atof( argv[5] );

  // write frequency
  int write_freq = std::atoi( argv[6] );

  // // device type
  // std::string device( argv[7] );

  // run the problem.
  rb_freely_rotating_square_bodies( body_spacing, body_height, body_length,
				    no_of_bodies,
				    t_final, write_freq);

  Kokkos::finalize();

  MPI_Finalize();
  return 0;
}
