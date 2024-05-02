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
using RigidBodyDataType = Cabana::MemberTypes<double[3], // position(0)
					      int[2], // limits(1)
					      int, // ids (2)
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
					      double, // moi for spherical particles (18)
					      double, // radius of the particle (19)
					      double, // inverse of moi for spherical particles (20)
					      double, // Youngs mod (21)
					      double, // Poisson's ratio (22)
					      int, // Total no of current contacts a dem particle is in (23)
					      int[10], // Indices of the contacts particle has (24)
					      double[10][3], // tangential force of dem particles (25)
					      double[10][3], // tangential disp of dem particles (26)
					      double[10], // normal force of the contact (27)
					      double[10], // normal overlap of the contact (28)
					      double // Shear modulus (29)
					      >;


// auto rb_position = Cabana::slice<0>     ( rb,    "rb_position");
// auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
// auto rb_ids = Cabana::slice<2>          ( rb,         "rb_ids");
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
// auto rb_rad_s = Cabana::slice<19>( rb,    "rb_rad_s");
// auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");
// auto rb_E = Cabana::slice<21>( rb,    "rb_E");
// auto rb_nu = Cabana::slice<22>( rb,    "rb_nu");
// auto rb_contacts_count = Cabana::slice<23>( rb,    "rb_contact_count");
// auto rb_contact_idx = Cabana::slice<24>( rb,    "rb_contact_idx");
// auto rb_contact_tng_frc = Cabana::slice<25>( rb,    "rb_contact_tng_frc");
// auto rb_contact_tng_disp = Cabana::slice<26>( rb,    "rb_contact_tng_disp");
// auto rb_contact_fn_magn = Cabana::slice<27>( rb,    "rb_contact_fn_magn");
// auto rb_contact_normal_overlap = Cabana::slice<28>( rb,    "rb_contact_normal_overlap");
// auto rb_G = Cabana::slice<29>( rb,    "rb_G");

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


void rigid_body_stage_1(RBAoSoA rb, double dt){
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");

  auto half_dt = dt * 0.5;
  auto rb_stage_1_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );
      auto mass_i_1 = 1. / mass_i;

      rb_velocity( i, 0 ) += half_dt * rb_force( i, 0 ) * mass_i_1;
      rb_velocity( i, 1 ) += half_dt * rb_force( i, 1 ) * mass_i_1;
      rb_velocity( i, 2 ) += half_dt * rb_force( i, 2 ) * mass_i_1;

      rb_ang_vel( i, 0 ) += half_dt * rb_torque( i, 0 ) * rb_moi_1( i );
      rb_ang_vel( i, 1 ) += half_dt * rb_torque( i, 1 ) * rb_moi_1( i );
      rb_ang_vel( i, 2 ) += half_dt * rb_torque( i, 2 ) * rb_moi_1( i );
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBStage1", policy,
			rb_stage_1_lambda_func );
}


void rigid_body_stage_2(RBAoSoA rb, double dt){
  auto rb_position = Cabana::slice<0>     ( rb,    "rb_position");
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");

  auto rb_stage_2_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      rb_position( i, 0 ) += rb_velocity( i, 0 ) * dt;
      rb_position( i, 1 ) += rb_velocity( i, 1 ) * dt;
      rb_position( i, 2 ) += rb_velocity( i, 2 ) * dt;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBStage2", policy,
			rb_stage_2_lambda_func );
}


void update_tangential_contacts(RBAoSoA rb, double dt, int * limits){
  auto rb_position = Cabana::slice<0>     ( rb,    "rb_position");
  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_ids = Cabana::slice<2>          ( rb,         "rb_ids");
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
  auto rb_rad_s = Cabana::slice<19>( rb,    "rb_rad_s");
  auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");
  auto rb_E = Cabana::slice<21>( rb,    "rb_E");
  auto rb_nu = Cabana::slice<22>( rb,    "rb_nu");
  auto rb_contacts_count = Cabana::slice<23>( rb,    "rb_contact_count");
  auto rb_contact_idx = Cabana::slice<24>( rb,    "rb_contact_idx");
  auto rb_contact_tng_frc = Cabana::slice<25>( rb,    "rb_contact_tng_frc");
  auto rb_contact_tng_disp = Cabana::slice<26>( rb,    "rb_contact_tng_disp");
  auto rb_contact_fn_magn = Cabana::slice<27>( rb,    "rb_contact_fn_magn");
  auto rb_contact_normal_overlap = Cabana::slice<28>( rb,    "rb_contact_normal_overlap");
  auto rb_G = Cabana::slice<29>( rb,    "rb_G");


  auto half_dt = dt * 0.5;
  auto update_tangential_contacts_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      int count = 0;
      int k = 0;
      int idx_total_ctcs = rb_contacts_count( i );
      int last_idx_tmp = rb_contacts_count( i ) - 1;
      int sidx = -1;
      // loop over all the contacts of particle d_idx
      while (count < idx_total_ctcs){
	// The index of the particle with which
	// d_idx in contact is
	sidx = rb_contact_idx( i, k );
	if (sidx == -1){
	  break;
	}
	else {
	  double pos_i[3] = {rb_position( i, 0 ),
	    rb_position( i, 1 ),
	    rb_position( i, 2 )};

	  double pos_j[3] = {rb_position( sidx, 0 ),
	    rb_position( sidx, 1 ),
	    rb_position( sidx, 2 )};

	  double pos_ij[3] = {rb_position( i, 0 ) - rb_position( sidx, 0 ),
	    rb_position( i, 1 ) - rb_position( sidx, 1 ),
	    rb_position( i, 2 ) - rb_position( sidx, 2 )};

	  // squared distance
	  double r2ij = pos_ij[0] * pos_ij[0] + pos_ij[1] * pos_ij[1] + pos_ij[2] * pos_ij[2];
	  // distance between i and j
	  double rij = sqrt(r2ij);
	  // Find the overlap amount
	  double overlap =  rb_rad_s( i ) + rb_rad_s( sidx ) - rij;

	  if (overlap <= 0.) {
	    // if the swap index is the current index then
	    // simply make it to null contact.
	    if (k == last_idx_tmp){
	      rb_contact_idx( i, k ) = -1;
	      rb_contact_tng_frc( i, k, 0 ) = 0.;
	      rb_contact_tng_frc( i, k, 1 ) = 0.;
	      rb_contact_tng_frc( i, k, 2 ) = 0.;
	      rb_contact_tng_disp( i, k, 0 ) = 0.;
	      rb_contact_tng_disp( i, k, 1 ) = 0.;
	      rb_contact_tng_disp( i, k, 2 ) = 0.;
	    }
	    else {
	      // swap the current tracking index with the final
	      // contact index
	      rb_contact_idx( i, k ) = rb_contact_idx( i, last_idx_tmp );
	      rb_contact_idx( i, last_idx_tmp ) = -1;

	      // swap tangential x frc
	      rb_contact_tng_frc( i, k, 0 ) = rb_contact_tng_frc( i, last_idx_tmp, 0 );
	      rb_contact_tng_frc( i, last_idx_tmp, 0 ) = 0.;

	      // swap tangential y frc
	      rb_contact_tng_frc( i, k, 1 ) = rb_contact_tng_frc( i, last_idx_tmp, 1 );
	      rb_contact_tng_frc( i, last_idx_tmp, 1 ) = 0.;

	      // swap tangential z frc
	      rb_contact_tng_frc( i, k, 2 ) = rb_contact_tng_frc( i, last_idx_tmp, 2 );
	      rb_contact_tng_frc( i, last_idx_tmp, 2 ) = 0.;

	      // swap tangential x displacement
	      rb_contact_tng_disp( i, k, 0 ) = rb_contact_tng_disp( i, last_idx_tmp, 0 );
	      rb_contact_tng_disp( i, last_idx_tmp, 0 ) = 0.;

	      // swap tangential y displacement
	      rb_contact_tng_disp( i, k, 1 ) = rb_contact_tng_disp( i, last_idx_tmp, 1 );
	      rb_contact_tng_disp( i, last_idx_tmp, 1 ) = 0.;

	      // swap tangential z displacement
	      rb_contact_tng_disp( i, k, 2 ) = rb_contact_tng_disp( i, last_idx_tmp, 2 );
	      rb_contact_tng_disp( i, last_idx_tmp, 2 ) = 0.;

	      // swap normal force and magnitude
	      rb_contact_fn_magn( i, k ) = rb_contact_fn_magn( i, last_idx_tmp );
	      rb_contact_fn_magn( i, last_idx_tmp ) = 0.;

	      rb_contact_normal_overlap( i, k ) = rb_contact_normal_overlap( i, last_idx_tmp );
	      rb_contact_normal_overlap( i, last_idx_tmp ) = 0.;

	      // decrease the last_idx_tmp, since we swapped it to
	      // -1
	      last_idx_tmp -= 1;
	    }

	    // decrement the total contacts of the particle
	    rb_contacts_count( i ) -= 1;
	  }
	  else
	    {
	      k = k + 1;
	    }
	}
	count += 1;
      }
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_mass.size());
  Kokkos::parallel_for( "CabanaSPH:Integrator:UpdateTngCnts", policy,
			update_tangential_contacts_lambda_func );
}


void compute_force_on_dem_particles(RBAoSoA rb, double dt,
				    ListType * verlet_list, double fric_coeff){
  auto rb_position = Cabana::slice<0>     ( rb,    "rb_position");
  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_ids = Cabana::slice<2>          ( rb,         "rb_ids");
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
  auto rb_rad_s = Cabana::slice<19>( rb,    "rb_rad_s");
  auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");
  auto rb_E = Cabana::slice<21>( rb,    "rb_E");
  auto rb_nu = Cabana::slice<22>( rb,    "rb_nu");
  auto rb_contacts_count = Cabana::slice<23>( rb,    "rb_contact_count");
  auto rb_contact_idx = Cabana::slice<24>( rb,    "rb_contact_idx");
  auto rb_contact_tng_frc = Cabana::slice<25>( rb,    "rb_contact_tng_frc");
  auto rb_contact_tng_disp = Cabana::slice<26>( rb,    "rb_contact_tng_disp");
  auto rb_contact_fn_magn = Cabana::slice<27>( rb,    "rb_contact_fn_magn");
  auto rb_contact_normal_overlap = Cabana::slice<28>( rb,    "rb_contact_normal_overlap");
  auto rb_G = Cabana::slice<29>( rb,    "rb_G");

  auto dem_contact_force_kernel = KOKKOS_LAMBDA( const int i, const int j )
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

      const double u_i = rb_velocity( i, 0 );
      const double v_i = rb_velocity( i, 1 );
      const double w_i = rb_velocity( i, 2 );

      const double u_j = rb_velocity( j, 0 );
      const double v_j = rb_velocity( j, 1 );
      const double w_j = rb_velocity( j, 2 );

      const double wx_i = rb_ang_vel( i, 0 );
      const double wy_i = rb_ang_vel( i, 1 );
      const double wz_i = rb_ang_vel( i, 2 );

      const double wx_j = rb_ang_vel( j, 0 );
      const double wy_j = rb_ang_vel( j, 1 );
      const double wz_j = rb_ang_vel( j, 2 );

      const double u1 = u_i + (nij_y * wz_i - nij_z * wy_i) * radius_i;
      const double u2 = -u_j + (nij_y * wz_j - nij_z * wy_j) * radius_j;
      const double uij = u1 + u2;
      const double v1 = v_i - (nij_x * wz_i - nij_z * wx_i) * radius_i;
      const double v2 = -v_j - (nij_x * wz_j - nij_z * wx_j) * radius_j;
      const double vij = v1 + v2;
      const double w1 = w_i - (nij_x * wy_i - nij_y * wx_i) * radius_i;
      const double w2 = -w_j - (nij_x * wy_j - nij_y * wx_j) * radius_j;
      const double wij = w1 + w2;

      const double vn = uij * nij_x + vij * nij_y + wij * nij_z;
      const double vn_x = vn * nij_x;
      const double vn_y = vn * nij_y;
      const double vn_z = vn * nij_z;

      const double overlap = radius_i + radius_j - dist;

      if (dist > 1e-12){

	// find the force if the particles are overlapping
	if (overlap > 0.) {
	  // Check if the particle is being tracked

	  // if the particle is being tracked then update the
	  // tracked tangential index
	  auto found = 0;
	  auto found_at = -1;
	  for ( std::size_t k = 0; k < rb_contacts_count( i ); ++k )
	    {
	      if (j == rb_contact_idx( i, k )) {
		found = 1;
		found_at = k;
	      }
	    }

	  if (found == 1)
	    {
	      // don't do anything
	    }
	  else{
	    // implies this is a new contact
	    // so add it to the contacting indices and increase the total
	    // count of the contacting indices

	    // create a new index at the end of the tracking
	    rb_contacts_count ( i ) += 1;
	    found_at = rb_contacts_count ( i ) - 1;
	    rb_contact_idx ( i, found_at ) = j;

	    rb_contact_tng_disp ( i, found_at, 0 ) = 0.;
	    rb_contact_tng_disp ( i, found_at, 1 ) = 0.;
	    rb_contact_tng_disp ( i, found_at, 2 ) = 0.;

	    rb_contact_tng_frc ( i, found_at, 0 ) = 0.;
	    rb_contact_tng_frc ( i, found_at, 1 ) = 0.;
	    rb_contact_tng_frc ( i, found_at, 2 ) = 0.;

	    rb_contact_fn_magn ( i, found_at ) = 0.;
	    rb_contact_normal_overlap ( i, found_at ) = 0.;
	  }

	  // =====================================================
	  // Start: Compute the stiffness constants and damping constants
	  // =====================================================
	  double tmp_1 = (1. - pow(rb_nu ( i ), 2.)) / rb_E ( i );
	  double tmp_2 = (1. - pow(rb_nu ( j ), 2.)) / rb_E ( j );
	  double E_eff = 1. / (tmp_1 + tmp_2);

	  tmp_1 = 1. / radius_i;
	  tmp_2 = 1. / radius_j;
	  double R_eff = 1. / (tmp_1 + tmp_2);
	  double kn = 4. / 3. * E_eff * sqrt(R_eff);
	  tmp_1 = (2. - rb_nu ( i )) / rb_G ( i );
	  tmp_2 = (2. - rb_nu ( j )) / rb_G ( j );
	  double G_eff = 1. / (tmp_1 + tmp_2);
	  double kt = 8. * G_eff * sqrt(R_eff * overlap);
	  // =====================================================
	  // End: Compute the stiffness constants and damping constants
	  // =====================================================


	  double fn_magn = kn * overlap;

	  rb_contact_fn_magn ( i, found_at ) = fn_magn;
	  rb_contact_normal_overlap ( i, found_at ) = overlap;

	  double fn_x = fn_magn * nij_x;
	  double fn_y = fn_magn * nij_y;
	  double fn_z = fn_magn * nij_z;

	  // Incremenet the tangential force and add it to the total force
	  // check if there is relative motion
	  double vij_magn = sqrt(uij*uij + vij*vij + wij*wij);
	  if (vij_magn < 1e-12) {
	    // make the tangential displacement
	    rb_contact_tng_disp ( i, found_at, 0 ) = 0.;
	    rb_contact_tng_disp ( i, found_at, 1 ) = 0.;
	    rb_contact_tng_disp ( i, found_at, 2 ) = 0.;

	    // make the tangential tangential force
	    rb_contact_tng_frc ( i, found_at, 0 ) = 0.;
	    rb_contact_tng_frc ( i, found_at, 1 ) = 0.;
	    rb_contact_tng_frc ( i, found_at, 2 ) = 0.;

	    // // # repulsive force
	    // d_fn_x[t2] = fn_x;
	    // d_fn_y[t2] = fn_y;
	    // d_fn_z[t2] = fn_z;
	  }
	  else{
	    // the tangential vector is
	    double tx_tmp = uij - vn_x;
	    double ty_tmp = vij - vn_y;
	    double tz_tmp = wij - vn_z;

	    double ti_magn = pow(tx_tmp*tx_tmp + ty_tmp*ty_tmp + tz_tmp*tz_tmp, 0.5);

	    double ti_x = 0.;
	    double ti_y = 0.;
	    double ti_z = 0.;

	    if (ti_magn > 1e-12){
	      ti_x = tx_tmp / ti_magn;
	      ti_y = ty_tmp / ti_magn;
	      ti_z = tz_tmp / ti_magn;
	    }

	    // // save the normals to output and view in viewer
	    // d_ti_x[t2] = ti_x;
	    // d_ti_y[t2] = ti_y;
	    // d_ti_z[t2] = ti_z;

	    // this is correct
	    double delta_lt_x_star = rb_contact_tng_disp ( i, found_at, 0 ) + uij * dt;
	    double delta_lt_y_star = rb_contact_tng_disp ( i, found_at, 1 ) + vij * dt;
	    double delta_lt_z_star = rb_contact_tng_disp ( i, found_at, 2 ) + wij * dt;

	    double delta_lt_dot_ti = (delta_lt_x_star * ti_x +
				      delta_lt_y_star * ti_y +
				      delta_lt_z_star * ti_z);

	    rb_contact_tng_disp ( i, found_at, 0 ) = delta_lt_dot_ti * ti_x;
	    rb_contact_tng_disp ( i, found_at, 1 ) = delta_lt_dot_ti * ti_y;
	    rb_contact_tng_disp ( i, found_at, 2 ) = delta_lt_dot_ti * ti_z;

	    double ft_x_star = -kt * rb_contact_tng_disp ( i, found_at, 0 );
	    double ft_y_star = -kt * rb_contact_tng_disp ( i, found_at, 1 );
	    double ft_z_star = -kt * rb_contact_tng_disp ( i, found_at, 2 );

	    double ft_magn = pow(ft_x_star*ft_x_star + ft_y_star*ft_y_star + ft_z_star*ft_z_star, 0.5);
	    double fn_magn = pow(fn_x*fn_x + fn_y*fn_y + fn_z*fn_z, 0.5);

	    double ft_magn_star = ft_magn;
	    if (ft_magn_star > fric_coeff * fn_magn){
	      ft_magn_star = fric_coeff * fn_magn;
	    }
	    // compute the tangential force, by equation 27
	    rb_contact_tng_frc ( i, found_at, 0 ) = -ft_magn_star * ti_x;
	    rb_contact_tng_frc ( i, found_at, 1 ) = -ft_magn_star * ti_y;
	    rb_contact_tng_frc ( i, found_at, 2 ) = -ft_magn_star * ti_z;

	    // reset the spring length
	    double modified_delta_lt_x = -rb_contact_tng_frc ( i, found_at, 0 ) / kt;
	    double modified_delta_lt_y = -rb_contact_tng_frc ( i, found_at, 1 ) / kt;
	    double modified_delta_lt_z = -rb_contact_tng_frc ( i, found_at, 2 ) / kt;

	    // repulsive force
	    rb_contact_fn_magn ( i, found_at ) = fn_magn;
	  }

	  double ft_x = rb_contact_tng_frc ( i, found_at, 0 );
	  double ft_y = rb_contact_tng_frc ( i, found_at, 1 );
	  double ft_z = rb_contact_tng_frc ( i, found_at, 2 );

	  rb_force( i, 0 ) += fn_x + ft_x;
	  rb_force( i, 1 ) += fn_y + ft_y;
	  rb_force( i, 2 ) += fn_z + ft_z;

	  // torque = n cross F;
	  rb_torque( i, 0 ) += (nij_y * ft_z - nij_z * ft_y) * radius_i;
	  rb_torque( i, 1 ) += (nij_z * ft_x - nij_x * ft_z) * radius_i;
	  rb_torque( i, 2 ) += (nij_x * ft_y - nij_y * ft_x) * radius_i;
	}

      }
    };

  Kokkos::RangePolicy<ExecutionSpace> policy(0, rb_mass.size());

  Cabana::neighbor_parallel_for( policy,
				 dem_contact_force_kernel,
				 *verlet_list,
				 Cabana::FirstNeighborsTag(),
				 Cabana::SerialOpTag(),
				 "CabanaDEM:Equations:ForceTorqueComputation" );
  Kokkos::fence();
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

      rb_torque( i, 0 ) = 0.;
      rb_torque( i, 1 ) = 0.;
      rb_torque( i, 2 ) = 0.;
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBBodyForce", policy,
			body_force_lambda_func );
}

void rigid_body_stage_3(RBAoSoA rb, double dt){
  auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
  auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
  auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
  auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
  auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
  auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");

  auto half_dt = dt * 0.5;
  auto rb_stage_3_lambda_func = KOKKOS_LAMBDA( const int i )
    {
      auto mass_i = rb_mass( i );
      auto mass_i_1 = 1. / mass_i;

      rb_velocity( i, 0 ) += half_dt * rb_force( i, 0 ) * mass_i_1;
      rb_velocity( i, 1 ) += half_dt * rb_force( i, 1 ) * mass_i_1;
      rb_velocity( i, 2 ) += half_dt * rb_force( i, 2 ) * mass_i_1;

      rb_ang_vel( i, 0 ) += half_dt * rb_torque( i, 0 ) * rb_moi_1( i );
      rb_ang_vel( i, 1 ) += half_dt * rb_torque( i, 1 ) * rb_moi_1( i );
      rb_ang_vel( i, 2 ) += half_dt * rb_torque( i, 2 ) * rb_moi_1( i );
    };
  Kokkos::RangePolicy<ExecutionSpace> policy( 0, rb_velocity.size() );
  Kokkos::parallel_for( "CabanaRB:Integrator:RBStage3", policy,
			rb_stage_3_lambda_func );
}


void output_rb_data(RBAoSoA rb, int num_particles, int step, double time)
{
  // This is for setting HDF5 options
  auto rb_position = Cabana::slice<0>     ( rb,    "rb_position");
  auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
  auto rb_ids = Cabana::slice<2>          ( rb,         "rb_ids");
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
  auto rb_rad_s = Cabana::slice<19>( rb,    "rb_rad_s");
  auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");
  auto rb_E = Cabana::slice<21>( rb,    "rb_E");
  auto rb_nu = Cabana::slice<22>( rb,    "rb_nu");
  auto rb_contacts_count = Cabana::slice<23>( rb,    "rb_contact_count");
  auto rb_contact_idx = Cabana::slice<24>( rb,    "rb_contact_idx");
  auto rb_contact_tng_frc = Cabana::slice<25>( rb,    "rb_contact_tng_frc");
  auto rb_contact_tng_disp = Cabana::slice<26>( rb,    "rb_contact_tng_disp");
  auto rb_contact_fn_magn = Cabana::slice<27>( rb,    "rb_contact_fn_magn");
  auto rb_contact_normal_overlap = Cabana::slice<28>( rb,    "rb_contact_normal_overlap");
  auto rb_G = Cabana::slice<29>( rb,    "rb_G");


  Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
  Cabana::Experimental::HDF5ParticleOutput::
    writeTimeStep(
		  h5_config, "particles", MPI_COMM_WORLD,
		  step, time, num_particles,
		  rb_position,
		  rb_ids,
		  rb_limits,
		  rb_velocity,
		  rb_force,
		  rb_torque,
		  rb_lin_acc,
		  rb_ang_acc,
		  rb_ang_mom,
		  rb_ang_vel,
		  rb_rotation_matrix,
		  rb_mass,
		  rb_density,
		  rb_body_moi,
		  rb_inv_body_moi,
		  rb_global_moi,
		  rb_inv_global_moi,
		  rb_rotation_angle,
		  rb_I_zz,
		  rb_rad_s,
		  rb_moi_1,
		  rb_E,
		  rb_nu,
		  rb_contacts_count,
		  rb_contact_idx,
		  rb_contact_tng_frc,
		  rb_contact_tng_disp,
		  rb_contact_fn_magn,
		  rb_contact_normal_overlap,
		  rb_G);
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
  auto rigid_body_rho = 2800.;
  auto rigid_body_radius = 0.04;
  auto rigid_body_diameter = 2. * rigid_body_radius;
  auto mass_rb_i = 4/3. * M_PI * pow(rigid_body_radius, 3) * rigid_body_rho;
  auto _I = 2. / 5. * mass_rb_i * pow(rigid_body_radius, 2);
  auto I_inverse = 1. / _I;
  auto E_rb = 4.8 * 1e5;
  auto nu_rb = 0.2;
  auto u_rb = 1.;

  // integration related variables
  double h = 1. * body_spacing;
  double dt = 0.25 * h / (c0 + vref);
  dt = 1e-4;
  std::cout << "dt is :" << dt << std::endl;
  auto final_time = t_final;
  auto time = 0.;
  int steps = final_time / dt;
  std::cout << "steps are :" << steps << std::endl;
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

  int no_rb_particles = 6;
  std::vector<int> is_fluid;
  std::vector<int> is_bdry;
  std::vector<int> is_rb;
  std::vector<int> body_id;

  std::vector<double> x = {
    0.,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5
  };

  std::vector<double> y = {
    0.,
    0.,
    0.,
    0.,
    0.,
    0.,
  };
  std::vector<double> z = {0., 0., 0., 0., 0., 0.};
  std::vector<double> u = {1., 0., -1., 1., 0., -1.};

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
  std::cout << "Total no of particles: " << num_particles << std::endl;

  // std::vector<int> rigid_limits = {8, 16};
  RBAoSoA rb( "rb", total_no_particles );
  auto rb_position = Cabana::slice<0>( rb,    "position");
  // create a mirror in the host space
  auto rb_host =
    Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), rb );

  auto rb_host_position = Cabana::slice<0>     ( rb_host,    "rb_host_position");
  auto rb_host_limits = Cabana::slice<1>       ( rb_host,      "rb_host_limits");
  auto rb_host_ids = Cabana::slice<2>          ( rb_host,         "rb_host_ids");
  auto rb_host_velocity = Cabana::slice<3>     ( rb_host,    "rb_host_velocity");
  auto rb_host_force = Cabana::slice<4>        ( rb_host,       "rb_host_force");
  auto rb_host_torque = Cabana::slice<5>       ( rb_host,      "rb_host_torque");
  auto rb_host_lin_acc = Cabana::slice<6>      ( rb_host,     "rb_host_lin_acc");
  auto rb_host_ang_acc = Cabana::slice<7>      ( rb_host,     "rb_host_ang_acc");
  auto rb_host_ang_mom = Cabana::slice<8>      ( rb_host,     "rb_host_ang_mom");
  auto rb_host_ang_vel = Cabana::slice<9>      ( rb_host,     "rb_host_ang_vel");
  auto rb_host_rotation_matrix = Cabana::slice<10>     ( rb_host,     "rb_host_rotation_matrix");
  auto rb_host_mass = Cabana::slice<11>        ( rb_host,        "rb_host_mass");
  auto rb_host_density = Cabana::slice<12>     ( rb_host,     "rb_host_density");
  auto rb_host_body_moi = Cabana::slice<13>    ( rb_host,    "rb_host_body_moi");
  auto rb_host_inv_body_moi = Cabana::slice<14>( rb_host,    "rb_host_inv_body_moi");
  auto rb_host_global_moi = Cabana::slice<15>  ( rb_host,    "rb_host_global_moi");
  auto rb_host_inv_global_moi = Cabana::slice<16>( rb_host,    "rb_host_inv_global_moi");
  auto rb_host_rotation_angle = Cabana::slice<17>  ( rb_host,    "rb_host_rotation_angle");
  auto rb_host_I_zz = Cabana::slice<18>( rb_host,    "rb_host_I_zz");
  auto rb_host_rad_s = Cabana::slice<19>( rb_host,    "rb_host_rad_s");
  auto rb_host_moi_1 = Cabana::slice<20>( rb_host,    "rb_host_moi_1");
  auto rb_host_E = Cabana::slice<21>( rb_host,    "rb_host_E");
  auto rb_host_nu = Cabana::slice<22>( rb_host,    "rb_host_nu");
  auto rb_host_contacts_count = Cabana::slice<23>( rb_host,    "rb_host_contact_count");
  auto rb_host_contact_idx = Cabana::slice<24>( rb_host,    "rb_host_contact_idx");
  auto rb_host_contact_tng_frc = Cabana::slice<25>( rb_host,    "rb_host_contact_tng_frc");
  auto rb_host_contact_tng_disp = Cabana::slice<26>( rb_host,    "rb_host_contact_tng_disp");
  auto rb_host_contact_fn_magn = Cabana::slice<27>( rb_host,    "rb_host_contact_fn_magn");
  auto rb_host_contact_normal_overlap = Cabana::slice<28>( rb_host,    "rb_host_contact_normal_overlap");
  auto rb_host_G = Cabana::slice<29>( rb_host,    "rb_host_G");

  for ( std::size_t i = 0; i < rb_host_position.size(); ++i )
    {

      std::cout << "Particle id i: " << i << std::endl;
      rb_host_ids ( i ) = i;

      rb_host_limits ( i, 0 ) = 0;
      rb_host_limits ( i, 1 ) = 0;

      rb_host_position ( i, 0 ) = x[i];
      rb_host_position ( i, 1 ) = y[i];
      rb_host_position ( i, 2 ) = z[i];

      for ( std::size_t j = 0; j < 3; ++j ) {
	rb_host_velocity ( i, j ) = 0.;

	rb_host_force ( i, j ) = 0.;

	rb_host_torque ( i, j ) = 0.;

	rb_host_lin_acc ( i, j ) = 0.;

	rb_host_ang_acc ( i, j ) = 0.;

	rb_host_ang_mom ( i, j ) = 0.;

	rb_host_ang_vel ( i, j ) = 0.;
      }
      // This is only set explicitly for this example
      rb_host_velocity ( i, 0 ) = u[i];
      // rb_host_velocity ( i, 1 ) = 0.;
      // rb_host_velocity ( i, 2 ) = pow( -1, i+1) * u_rb / 10.;

      for ( std::size_t j = 0; j < 9; ++j ) {
	rb_host_rotation_matrix( i, j ) = 0.;
	rb_host_body_moi( i, j ) = 0.;
	rb_host_inv_body_moi( i, j ) = 0.;
	rb_host_global_moi( i, j ) = 0.;
	rb_host_inv_global_moi( i, j ) = 0.;
      }

      rb_host_mass ( i ) = mass_rb_i;
      rb_host_density ( i ) = rigid_body_rho;

      rb_host_rotation_angle( i ) = 0.;

      rb_host_I_zz ( i ) = _I;
      rb_host_rad_s ( i ) = rigid_body_radius;
      rb_host_moi_1 ( i ) = I_inverse;
      rb_host_E ( i ) = E_rb;
      rb_host_nu ( i ) = nu_rb;

      rb_host_contacts_count ( i ) = 0;

      for ( std::size_t j = 0; j < 10; ++j ) {
	rb_host_contact_idx ( i, j ) = -1;
	rb_host_contact_tng_frc ( i, j, 0 ) = 0.;
	rb_host_contact_tng_frc ( i, j, 1 ) = 0.;
	rb_host_contact_tng_frc ( i, j, 2 ) = 0.;
	rb_host_contact_tng_disp ( i, j, 0 ) = 0.;
	rb_host_contact_tng_disp ( i, j, 1 ) = 0.;
	rb_host_contact_tng_disp ( i, j, 2 ) = 0.;

	rb_host_contact_fn_magn ( i, j ) = 0.;
	rb_host_contact_normal_overlap ( i, j ) = 0.;
      }

      rb_host_G ( i ) = E_rb / (2. * (1. + nu_rb));
    }
  // copy it back to aosoa
  Cabana::deep_copy( rb, rb_host );

  int rigid_limits[2] = {0, 2};

  output_rb_data(rb, num_particles, 0, time);
  // ================================================
  // ================================================
  // create the neighbor list
  // ================================================
  // ================================================
  // double neighborhood_radius_rb = 1. * rigid_body_radius;
  // double grid_min_rb[3] = { -20. * rigid_body_radius, -20. * rigid_body_radius, -neighborhood_radius_rb };
  // double grid_max_rb[3] = {  20. * rigid_body_radius, 20. * rigid_body_radius, neighborhood_radius_rb};
  // double cell_ratio = 1.0;
  // ListType verlet_list_rb( rb_position, 0,
  // 			   rb_position.size(), neighborhood_radius_rb,
  // 			   cell_ratio, grid_min_rb, grid_max_rb );

  std::cout << "Making a liked list now:===================" << std::endl;
  double neighborhood_radius = 4. * rigid_body_radius;
  double grid_min[3] = { -10. * rigid_body_radius, -10. * rigid_body_radius, -neighborhood_radius - rigid_body_radius};
  double grid_max[3] = { 10. * rigid_body_radius, 10. * rigid_body_radius, neighborhood_radius + rigid_body_radius};
  // double grid_min[3] = { 0.0, -4.0, 0.0 };
  // double grid_max[3] = { 4.1, 4.0, 0.0 };
  double cell_ratio = 1.0;

  ListType verlet_list( rb_position, 0,
			rb_position.size(), neighborhood_radius,
			cell_ratio, grid_min, grid_max );
  std::cout << "A working linked list is produced :" << std::endl;

  // Main timestep loop
  // steps = 0;
  for ( int step = 0; step < steps; step++ )
    {
      rigid_body_stage_1(rb, dt);

      rigid_body_stage_2(rb, dt);

      verlet_list.build( rb_position, 0, rb_position.size(), neighborhood_radius,
			 cell_ratio, grid_min, grid_max );

      body_force_rigid_body(rb, gravity, dt);

      update_tangential_contacts(rb, dt, rigid_limits);
      compute_force_on_dem_particles(rb, dt,
				     &verlet_list, 0.1);

      rigid_body_stage_3(rb, dt);

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

      std::cerr << "\nfor example: ./TestDEM01MultipleContactsAtOnce  0.1 1.0 1.0 6 0.1 10\n";
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


// import h5py
// import numpy as np
// import matplotlib.pyplot as plt
// import os


// # start: Get the files
// files = [filename for filename in os.listdir('.') if filename.startswith("particles") and filename.endswith("xmf") ]
// files.sort()
// files_num = []
// for f in files:
//     f_last = f[10:]
//     files_num.append(int(f_last[:-4]))
// files_num.sort()

// sorted_files = []
// for num in files_num:
//     sorted_files.append("particles_" + str(num) + ".h5")
// print(sorted_files)
// # ends: Get the files
// files = sorted_files

// contact_count = []
// for f in files:
//     f = h5py.File(f, "r")
//     # print(np.array(f["radius"]))
//     x = np.array(f["rb_position"][:, 0])
//     y = np.array(f["rb_position"][:, 1])
//     idx_check = 1
//     contact_count.append(f["rb_contact_count"][idx_check])

//     cnt = f["rb_contact_count"][idx_check]
//     if cnt > 0:
//         # print(f["rb_contact_idx"][1])
//         print("----------------------------")
//         # print(f["rb_contact_tng_frc"][0])
//         # print(f["rb_contact_tng_frc"][1])
//         print(f["rb_contact_idx"][idx_check])
//         # print(f["rb_contact_tng_frc"][idx_check])
//         print(f["rb_contact_fn_magn"][idx_check])
//         print("----------------------------")

// # print(x)
// # print(y)
// print(contact_count)

// # plt.scatter(x, y, label="SPH appr")
// # plt.legend()
// # plt.savefig("colliding_fluid_blocks.png")
// # plt.show()
// # plt.plot()
