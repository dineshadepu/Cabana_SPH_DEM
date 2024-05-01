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

/****************************************************************************
 * We follow this article's code for implementation
 * https://github.com/pmocz/sph-python/blob/master/sph.py *
 ****************************************************************************/

#include <Cabana_Core.hpp>
#include <math.h>

#include <iostream>

#define DIM 2


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
				      double[3] // dummy velocity for computation (14)
				      >;

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

using ListAlgorithm = Cabana::FullNeighborTag;
using ListType =
  Cabana::VerletList<MemorySpace, ListAlgorithm, Cabana::VerletLayout2D>;

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

std::tuple<std::vector<double>,
	   std::vector<double>,
	   std::vector<double>> create_2d_tank(double length, double height,
					       double spacing,
					       int no_layers){
  std::vector<double> x_bottom;
  std::vector<double> y_bottom;
  std::vector<double> z_bottom;
  int x_no_points = length / spacing;
  int y_no_points = (no_layers * spacing) / spacing;
  for(int i=0; i<=x_no_points; i++)
    {
      double x = i * spacing;
      for(int j=0; j<=y_no_points; j++)
	{
	  double y = j * spacing;
	  x_bottom.push_back(x);
	  y_bottom.push_back(y);
	  z_bottom.push_back(0.);
	}
    }

  std::vector<double> x_left;
  std::vector<double> y_left;
  std::vector<double> z_left;
  x_no_points = (no_layers * spacing) / spacing;
  y_no_points = height / spacing + (no_layers * spacing) / spacing;
  for(int i=0; i<=x_no_points; i++)
    {
      double x = i * spacing;
      for(int j=0; j<=y_no_points; j++)
	{
	  double y = j * spacing;
	  x_left.push_back(x);
	  y_left.push_back(y);
	  z_left.push_back(0.);
	}
    }

  std::vector<double> x_right;
  std::vector<double> y_right;
  std::vector<double> z_right;
  x_no_points = (no_layers * spacing) / spacing;
  y_no_points = height / spacing + (no_layers * spacing) / spacing;
  for(int i=0; i<=x_no_points; i++)
    {
      double x = i * spacing;
      for(int j=0; j<=y_no_points; j++)
	{
	  double y = j * spacing;
	  x_right.push_back(x);
	  y_right.push_back(y);
	  z_right.push_back(0.);
	}
    }

  /*
  // move the blocks to their appropriate positions
  */
  // first move left block to the left
  auto x_bottom_min = *std::min_element(x_bottom.begin(), x_bottom.end());
  auto x_left_max = *std::max_element(x_left.begin(), x_left.end());
  for(double& d : x_left)
    d -= (x_left_max - x_bottom_min) + 1. * spacing;
  // now move it to the correct y location
  auto y_bottom_min = *std::min_element(y_bottom.begin(), y_bottom.end());
  auto y_left_min = *std::min_element(y_left.begin(), y_left.end());
  for(double& d : y_left)
    d += (y_bottom_min - y_left_min);

  // next move the right block to the right
  auto x_bottom_max = *std::max_element(x_bottom.begin(), x_bottom.end());
  auto x_right_min = *std::min_element(x_right.begin(), x_right.end());
  for(double& d : x_right)
    d += (x_bottom_max - x_right_min) + spacing * 1;
  // now move it to the correct y location
  auto y_right_min = *std::min_element(y_right.begin(), y_right.end());
  for(double& d : y_right)
    d += (y_bottom_min - y_right_min);

  // The tank with all three blocks combined
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  for(int i=0; i<x_bottom.size(); i++)
    {
      x.push_back(x_bottom[i]);
      y.push_back(y_bottom[i]);
      z.push_back(z_bottom[i]);
    }
  for(int i=0; i<x_left.size(); i++)
    {
      x.push_back(x_left[i]);
      y.push_back(y_left[i]);
      z.push_back(z_left[i]);
    }
  for(int i=0; i<x_right.size(); i++)
    {
      x.push_back(x_right[i]);
      y.push_back(y_right[i]);
      z.push_back(z_right[i]);
    }

  return std::make_tuple(x, y, z);
}


KOKKOS_INLINE_FUNCTION
void compute_quintic_wij(double rij, double h, double *result){
  double h1 =  1. / h;
  double q =  rij * h1;
  double fac = M_1_PI *  7. / 478. * h1 * h1;

  double tmp3 = 3. - q;
  double tmp2 = 2. - q;
  double tmp1 = 1. - q;

  double val = 0.;
  if (q > 3.) {
    val = 0.;
  } else if ( q > 2.) {
    val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3;
  } else if ( q > 1.) {
    val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3;
    val -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2;
  } else {
    val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3;
    val -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2;
    val += 15. * tmp1 * tmp1 * tmp1 * tmp1 * tmp1;
  }

  *result = val * fac;
}


KOKKOS_INLINE_FUNCTION
void compute_quintic_gradient_wij(double *xij, double rij, double h, double *result){
  double h1 =  1. / h;
  double q =  rij * h1;

  double fac = M_1_PI *  7. / 478. * h1 * h1;

  double tmp3 = 3. - q;
  double tmp2 = 2. - q;
  double tmp1 = 1. - q;

  double val = 0.;
  if (rij > 1e-12){
    if (q > 3.) {
      val = 0.;
    } else if ( q > 2.) {
      val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3;
      val *= h1 / rij;
    } else if ( q > 1.) {
      val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3;
      val += 30.0 * tmp2 * tmp2 * tmp2 * tmp2;
      val *= h1 / rij;
    } else {
      val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3;
      val += 30.0 * tmp2 * tmp2 * tmp2 * tmp2;
      val -= 75.0 * tmp1 * tmp1 * tmp1 * tmp1;
      val *= h1 / rij;
    }
  } else {
    val = 0.;
  }

  double tmp = val * fac;
  result[0] = tmp * xij[0];
  result[1] = tmp * xij[1];
  result[2] = tmp * xij[2];
}


void fluid_stage_1(AoSoAType & aosoa, double dt, int * limits){
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_acc = Cabana::slice<3>          ( aosoa,    "acc");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto aosoa_is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");

  auto half_dt = dt * 0.5;

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      aosoa_velocity(i, 0) += aosoa_acc(i, 0) * half_dt;
      aosoa_velocity(i, 1) += aosoa_acc(i, 1) * half_dt;
      aosoa_velocity(i, 2) += aosoa_acc(i, 2) * half_dt;
    }
  }
}


void fluid_stage_2(AoSoAType & aosoa, double dt, int * limits){
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_density = Cabana::slice<5>     ( aosoa,    "density");
  auto aosoa_p = Cabana::slice<7>           ( aosoa,    "p");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto aosoa_density_acc = Cabana::slice<10>           ( aosoa,    "density_acc");
  auto aosoa_p_acc = Cabana::slice<11>           ( aosoa,    "p_acc");

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      aosoa_density( i ) += aosoa_density_acc( i ) * dt;
      aosoa_p( i ) += aosoa_p_acc( i ) * dt;

      aosoa_position(i, 0) += aosoa_velocity( i, 0 ) * dt;
      aosoa_position(i, 1) += aosoa_velocity( i, 1 ) * dt;
      aosoa_position(i, 2) += aosoa_velocity( i, 2 ) * dt;
    }
  }
}


void fluid_stage_3(AoSoAType & aosoa, double dt, int * limits){
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_acc = Cabana::slice<3>          ( aosoa,    "acc");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto aosoa_is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");

  auto half_dt = dt * 0.5;

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      aosoa_velocity(i, 0) += aosoa_acc(i, 0) * half_dt;
      aosoa_velocity(i, 1) += aosoa_acc(i, 1) * half_dt;
      aosoa_velocity(i, 2) += aosoa_acc(i, 2) * half_dt;
    }
  }
}


void continuity_equation(AoSoAType &aosoa, double dt,
			 ListType * verlet_list,
			 int * limits){
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
  auto aosoa_ids = Cabana::slice<1>          ( aosoa,    "ids");
  auto aosoa_velocity = Cabana::slice<2>     ( aosoa,    "velocity");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_density = Cabana::slice<5>     ( aosoa,    "density");
  auto aosoa_h = Cabana::slice<6>           ( aosoa,    "h");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto aosoa_is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");
  auto aosoa_density_acc = Cabana::slice<10>           ( aosoa,    "density_acc");

  // Make the acceleration zero first
  Cabana::deep_copy( aosoa_density_acc, 0. );


  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      for (int j = 0; j < aosoa_mass.size() ; j++) {
	/*
	  Common to all equations in SPH.

	  We compute:
	  1.the vector passing from j to i
	  2. Distance between the points i and j
	  3. Distance square between the points i and j
	  4. Velocity vector difference between i and j
	  5. Kernel value
	  6. Derivative of kernel value
	*/
	double pos_i[3] = {aosoa_position( i, 0 ),
	  aosoa_position( i, 1 ),
	  aosoa_position( i, 2 )};

	double pos_j[3] = {aosoa_position( j, 0 ),
	  aosoa_position( j, 1 ),
	  aosoa_position( j, 2 )};

	double pos_ij[3] = {aosoa_position( i, 0 ) - aosoa_position( j, 0 ),
	  aosoa_position( i, 1 ) - aosoa_position( j, 1 ),
	  aosoa_position( i, 2 ) - aosoa_position( j, 2 )};

	double vel_ij[3] = {aosoa_velocity( i, 0 ) - aosoa_velocity( j, 0 ),
	  aosoa_velocity( i, 1 ) - aosoa_velocity( j, 1 ),
	  aosoa_velocity( i, 2 ) - aosoa_velocity( j, 2 )};

	// wij and dwij
	// double wij = 0.;
	double dwij[3] = {0., 0., 0.};

	// h value of particle i
	double h_i = aosoa_h( i );

	// squared distance
	double r2ij = pos_ij[0] * pos_ij[0] + pos_ij[1] * pos_ij[1] + pos_ij[2] * pos_ij[2];
	// distance between i and j
	double rij = sqrt(r2ij);
	// compute the kernel wij
	// compute_quintic_wij(rij, h_i, &wij);
	// compute the gradient of kernel dwij
	compute_quintic_gradient_wij(pos_ij, rij, h_i, dwij);

	/*
	  ====================================
	  End: common to all equations in SPH.
	  ====================================
	*/

	// const double mass_i = aosoa_mass( i );
	const double mass_j = aosoa_mass( j );

	double vijdotdwij = dwij[0]*vel_ij[0] + dwij[1]*vel_ij[1] + dwij[2]*vel_ij[2];
	aosoa_density_acc (i) += mass_j * vijdotdwij;
      }
    }
  }
}


void state_equation(AoSoAType & aosoa, ViewVectorType & b_rho0_p0, double dt, int * limits){
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_density = Cabana::slice<5>     ( aosoa,    "aosoa_density");
  auto aosoa_p = Cabana::slice<7>           ( aosoa,    "aosoa_p");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      double tmp = b_rho0_p0[2] * (aosoa_density( i ) / b_rho0_p0[1] - b_rho0_p0[0]);
      aosoa_p( i ) = tmp;
    }
  }
}



void solid_wall_pressure_bc(AoSoAType & aosoa, ViewVectorType & gravity, double dt,
			    ListType * verlet_list,
			    int * limits){
  auto aosoa_ids = Cabana::slice<1>          ( aosoa,    "ids");
  auto aosoa_position = Cabana::slice<0>     ( aosoa,    "position");
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

  Cabana::deep_copy( aosoa_sum_wij, 0. );

  // Divide by wij as the end
  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_boundary(i) == 1 ){
      aosoa_p( i ) = 0.;
    }
  }


  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_boundary(i) == 1 ){
      for (int j = 0; j < aosoa_mass.size() ; j++) {
	if (aosoa_is_fluid(j) == 1 ){
	  /*
	    Common to all equations in SPH.

	    We compute:
	    1.the vector passing from j to i
	    2. Distance between the points i and j
	    3. Distance square between the points i and j
	    4. Velocity vector difference between i and j
	    5. Kernel value
	    6. Derivative of kernel value
	  */
	  double pos_i[3] = {aosoa_position( i, 0 ),
	    aosoa_position( i, 1 ),
	    aosoa_position( i, 2 )};

	  double pos_j[3] = {aosoa_position( j, 0 ),
	    aosoa_position( j, 1 ),
	    aosoa_position( j, 2 )};

	  double pos_ij[3] = {aosoa_position( i, 0 ) - aosoa_position( j, 0 ),
	    aosoa_position( i, 1 ) - aosoa_position( j, 1 ),
	    aosoa_position( i, 2 ) - aosoa_position( j, 2 )};

	  double vel_ij[3] = {aosoa_velocity( i, 0 ) - aosoa_velocity( j, 0 ),
	    aosoa_velocity( i, 1 ) - aosoa_velocity( j, 1 ),
	    aosoa_velocity( i, 2 ) - aosoa_velocity( j, 2 )};

	  // wij and dwij
	  double wij = 0.;
	  // double dwij[3] = {0., 0., 0.};

	  // h value of particle i
	  double h_i = aosoa_h( i );

	  // squared distance
	  double r2ij = pos_ij[0] * pos_ij[0] + pos_ij[1] * pos_ij[1] + pos_ij[2] * pos_ij[2];
	  // distance between i and j
	  double rij = sqrt(r2ij);
	  // compute the kernel wij
	  compute_quintic_wij(rij, h_i, &wij);
	  // // compute the gradient of kernel dwij
	  // compute_quintic_gradient_wij(pos_ij, rij, h_i, dwij);

	  /*
	    ====================================
	    End: common to all equations in SPH.
	    ====================================
	  */
	  double tmp1 = (gravity[0] - aosoa_acc( i, 0 )) * pos_ij [ 0 ];
	  double tmp2 = (gravity[1] - aosoa_acc( i, 1 )) * pos_ij [ 1 ];
	  double tmp3 = (gravity[2] - aosoa_acc( i, 2 )) * pos_ij [ 2 ];
	  double gdotxij = tmp1 +  tmp2 + tmp3;

	  // pressure of the boundary particle
	  aosoa_p ( i ) += aosoa_p ( j ) * wij + aosoa_density ( j ) * gdotxij * wij;

	  // sum the wij (source number density)
	  aosoa_sum_wij ( i ) += wij;
	}
      }
    }
  }

  // Divide by wij as the end
  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_boundary(i) == 1 ){
      if (aosoa_sum_wij( i ) > 1e-12){
	aosoa_p( i ) /= aosoa_sum_wij( i );
      }
    }
  }
}


void clamp_wall_pressure(AoSoAType & aosoa, double dt, int * limits){
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_p = Cabana::slice<7>           ( aosoa,    "p");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");
  auto aosoa_is_boundary = Cabana::slice<9>    ( aosoa,    "is_boundary");

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_boundary(i) == 1 ){
      if (aosoa_p( i ) < 0. ){
	aosoa_p( i ) = 0.;
      }
    }
  }
}


void momentum_equation(AoSoAType & aosoa, ViewVectorType & gravity,
		       ViewVectorType & alpha_c0,
		       double dt,
		       ListType * verlet_list,
		       int * limits){
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

  // Cabana::deep_copy( aosoa_acc, 0. );

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      for (int j = 0; j < aosoa_mass.size() ; j++) {
	/*
	  Common to all equations in SPH.

	  We compute:
	  1.the vector passing from j to i
	  2. Distance between the points i and j
	  3. Distance square between the points i and j
	  4. Velocity vector difference between i and j
	  5. Kernel value
	  6. Derivative of kernel value
	*/
	double pos_i[3] = {aosoa_position( i, 0 ),
	  aosoa_position( i, 1 ),
	  aosoa_position( i, 2 )};

	double pos_j[3] = {aosoa_position( j, 0 ),
	  aosoa_position( j, 1 ),
	  aosoa_position( j, 2 )};

	double pos_ij[3] = {aosoa_position( i, 0 ) - aosoa_position( j, 0 ),
	  aosoa_position( i, 1 ) - aosoa_position( j, 1 ),
	  aosoa_position( i, 2 ) - aosoa_position( j, 2 )};

	double vel_ij[3] = {aosoa_velocity( i, 0 ) - aosoa_velocity( j, 0 ),
	  aosoa_velocity( i, 1 ) - aosoa_velocity( j, 1 ),
	  aosoa_velocity( i, 2 ) - aosoa_velocity( j, 2 )};

	// wij and dwij
	// double wij = 0.;
	double dwij[3] = {0., 0., 0.};

	// h value of particle i
	double h_i = aosoa_h( i );

	// squared distance
	double r2ij = pos_ij[0] * pos_ij[0] + pos_ij[1] * pos_ij[1] + pos_ij[2] * pos_ij[2];
	// distance between i and j
	double rij = sqrt(r2ij);
	// compute the kernel wij
	// compute_quintic_wij(rij, h_i, &wij);
	// compute the gradient of kernel dwij
	compute_quintic_gradient_wij(pos_ij, rij, h_i, dwij);

	/*
	  ====================================
	  End: common to all equations in SPH.
	  ====================================
	*/
	double rhoi2 = 1. / (aosoa_density( i ) * aosoa_density( i ));
	double rhoj2 = 1. / (aosoa_density( j ) * aosoa_density( j ));
	double pij = aosoa_p ( i ) * rhoi2 + aosoa_p ( j ) * rhoj2;
	double tmp = - aosoa_mass ( j ) * pij;

	// pressure acceleration
	aosoa_acc ( i, 0 ) += tmp * dwij[0];
	aosoa_acc ( i, 1 ) += tmp * dwij[1];
	aosoa_acc ( i, 2 ) += tmp * dwij[2];

	// artificial viscosity
	const double mass_i = aosoa_mass( i );
	const double mass_j = aosoa_mass( j );

	double vijdotrij = vel_ij[0]*pos_ij[0] + vel_ij[1]*pos_ij[1] + vel_ij[2]*pos_ij[2];

	double piij = 0.;

	if (vijdotrij < 0.) {
	  double muij = (h_i * vijdotrij)/(r2ij + 0.01 * h_i * h_i);
	  double rhoij = 0.5*(aosoa_density( i ) + aosoa_density( j ));
	  double rhoij1 = 1. / rhoij;

	  piij = -alpha_c0[0] * alpha_c0[1] * muij;
	  piij = aosoa_mass( j ) * piij * rhoij1;
	}

	aosoa_acc ( i, 0 ) += -piij * dwij[0];
	aosoa_acc ( i, 1 ) += -piij * dwij[1];
	aosoa_acc ( i, 2 ) += -piij * dwij[2];

	// // Real viscosity
	// double etai = nu[0] * aosoa_density( i );
	// double etaj = nu[0] * aosoa_density( j );

	// double etaij = 4 * (etai * etaj)/(etai + etaj);

	// double xdotdij = dwij[0]*pos_ij[0] + dwij[1]*pos_ij[1] + dwij[2]*pos_ij[2];

	// double tmp = mass_j / (aosoa_density( i ) * aosoa_density( j ));
	// double fac = tmp * etaij * xdotdij/(r2ij + 0.01 * h_i * h_i);

	// aosoa_acc ( i, 0 ) += fac * vel_ij[0];
	// aosoa_acc ( i, 1 ) += fac * vel_ij[1];
	// aosoa_acc ( i, 2 ) += fac * vel_ij[2];
      }
    }
  }
};


void body_force(AoSoAType & aosoa, ViewVectorType & gravity, double dt,
		ListType * verlet_list, int * limits){
  auto aosoa_acc = Cabana::slice<3>          ( aosoa,    "acc");
  auto aosoa_mass = Cabana::slice<4>        ( aosoa,    "mass");
  auto aosoa_is_fluid = Cabana::slice<8>    ( aosoa,    "is_fluid");

  for (int i = 0; i < aosoa_mass.size() ; i++) {
    if (aosoa_is_fluid(i) == 1 ){
      aosoa_acc( i, 0 ) = gravity[0];
      aosoa_acc( i, 1 ) = gravity[1];
      aosoa_acc( i, 2 ) = gravity[2];
    }
  }
}


void setup_aosoa(AoSoAType aosoa, double *x, double *y, double *z, double *u){
  auto aosoa_position = Cabana::slice<0>( aosoa,    "position");
  auto aosoa_ids = Cabana::slice<1>( aosoa,    "ids");
  auto aosoa_velocity = Cabana::slice<2>( aosoa,    "velocity");
  auto aosoa_acc = Cabana::slice<3>( aosoa,    "acc");
  auto aosoa_mass = Cabana::slice<4>( aosoa,    "mass");
  auto aosoa_density = Cabana::slice<5>( aosoa,    "density");
  auto aosoa_h = Cabana::slice<6>( aosoa,    "h");
  auto aosoa_p = Cabana::slice<7>( aosoa,    "p");
  auto aosoa_is_fluid = Cabana::slice<8>( aosoa,    "is_fluid");
  auto aosoa_is_boundary = Cabana::slice<9>( aosoa,    "is_boundary");
  auto aosoa_density_acc = Cabana::slice<10>( aosoa,    "density_acc");
  auto aosoa_p_acc = Cabana::slice<11>( aosoa,    "p_acc");
  auto aosoa_sum_wij = Cabana::slice<12>( aosoa,    "wij");

  // auto aosoa_init_func = KOKKOS_LAMBDA( const int i )
  //   {
  //     aosoa_position ( i, 0 ) = x[i];
  //     aosoa_position ( i, 1 ) = y[i];
  //     aosoa_position ( i, 2 ) = z[i];

  //     aosoa_ids ( i ) = i;

  //     aosoa_velocity ( i, 0 ) = u[i];
  //     aosoa_velocity ( i, 1 ) = 0.;
  //     aosoa_velocity ( i, 2 ) = 0.;

  //     aosoa_acc ( i, 0 ) = 0.;
  //     aosoa_acc ( i, 1 ) = 0.;
  //     aosoa_acc ( i, 2 ) = 0.;

  //     aosoa_mass ( i ) = fluid_rho * pow(fluid_spacing, DIM);
  //     aosoa_density ( i ) = fluid_rho;
  //     aosoa_h ( i ) = 1. * fluid_spacing;

  //     aosoa_p ( i ) = 0.;
  //     // aosoa_is_fluid ( i ) = Cabana::slice<8>( aosoa,    "is_fluid");
  //     // aosoa_is_boundary ( i ) = Cabana::slice<9>( aosoa,    "is_boundary");
  //     aosoa_density_acc ( i ) = 0.;
  //     aosoa_p_acc ( i ) = 0.;
  //     aosoa_sum_wij ( i ) = 0.;
  //   };

  // Kokkos::RangePolicy<ExecutionSpace> policy( 0, aosoa_position.size() );
  // Kokkos::parallel_for( "init_aosoa",  policy,
  // 			aosoa_init_func );
}


void output_data(AoSoAType & aosoa, int num_particles, int step, double time)
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


  Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
  Cabana::Experimental::HDF5ParticleOutput::
    writeTimeStep(
		  h5_config, "particles", MPI_COMM_WORLD,
		  step, time, num_particles,
		  position,
		  ids ,
		  velocity,
		  acc,
		  mass,
		  density,
		  h,
		  p,
		  is_fluid,
		  is_boundary,
		  density_acc,
		  p_acc,
		  sum_wij
		  );
}


//---------------------------------------------------------------------------//
// TODO: explain this function in short
//---------------------------------------------------------------------------//
void damBreak(const double fluid_spacing,
	      const double fluid_height,
	      const double fluid_length,
	      const double tank_height,
	      const double tank_length,
	      const double t_final,
	      const int write_freq)
{
  /* This simulation is setup in three parts

     1. Create the particles (Both fluid and boundary)
     2. Assign the particle properties to the aosoa
     2a. This includes setting the limits of the particles
     3. Set the time step, time and pfreq
     4. Execute the run
     4a. Set up the neighbour list
    */

    int comm_rank = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

    if ( comm_rank == 0 )
      std::cout << "Cabana Rigid body solver example\n" << std::endl;

    /*
      ================================================
      Step 1
      ================================================
      1. Create the particles (Both fluid and boundary)
    */
    // double fluid_length = 1.;
    // double fluid_height = 0.6;
    // double fluid_spacing = 0.1;
    // double tank_length = 1.;
    // double tank_height = 1.7 * fluid_height;
    int tank_layers = 3;


    std::vector<double> x_fluid, y_fluid, z_fluid;
    std::tie(x_fluid, y_fluid, z_fluid) = create_2d_block(fluid_length, fluid_height,
							  fluid_spacing);
    std::vector<double> u_fluid(x_fluid.size(), 0);

    std::vector<double> x_tank, y_tank, z_tank;
    std::tie(x_tank, y_tank, z_tank) = create_2d_tank(tank_length, tank_height,
						      fluid_spacing, tank_layers);
    for ( std::size_t i = 0; i < x_fluid.size(); ++i )
      {
	x_fluid[i] += (tank_layers - 2) * fluid_spacing;
	y_fluid[i] += (tank_layers + 1) * fluid_spacing;
	x_fluid[i] -= (1) * fluid_spacing;
      // y_fluid[i] += (2) * fluid_spacing;
    }
  // for ( std::size_t i = 0; i < x_fluid.size(); ++i )
  //   {
  //     x_fluid[i] += fluid_spacing;
  //     y_fluid[i] += fluid_spacing;
  //   }


  int no_fluid_particles = x_fluid.size();
  int no_bdry_particles = x_tank.size();

  int total_no_particles = no_fluid_particles + no_bdry_particles;
  int num_particles = total_no_particles;

  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> u;
  std::vector<int> is_fluid;
  std::vector<int> is_bdry;
  for ( std::size_t i = 0; i < x_fluid.size(); ++i )
    {
      x.push_back(x_fluid[i]);
      y.push_back(y_fluid[i]);
      z.push_back(z_fluid[i]);
      u.push_back(0.);
      is_fluid.push_back(1);
      is_bdry.push_back(0);
    }

  for ( std::size_t i = 0; i < x_tank.size(); ++i )
    {
      x.push_back(x_tank[i]);
      y.push_back(y_tank[i]);
      z.push_back(z_tank[i]);
      u.push_back(0.);
      is_fluid.push_back(0);
      is_bdry.push_back(1);
    }
  /*
    ================================================
    End: Step 1
    ================================================
  */

  /*
    ================================================
    Step 2: Set the time step, time and pfreq
    ================================================
  */
  // Material properties of the fluid
  auto gx = 0.;
  auto gy = -9.81;
  auto gz = 0.;
  auto fluid_rho = 1000.;
  auto vref = sqrt(2. * abs(gy) * fluid_height);
  auto c0 = 10. * vref;
  // auto fluid_viscosity = 1e-4;
  auto alpha = 0.01;
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

  // integration related variables
  double h = 1. * fluid_spacing;
  // double dt = 0.25 * h / (c0 + vref);
  double dt = 1e-5;
  std::cout << "dt is :" << dt << std::endl;
  auto final_time = t_final;
  auto time = 0.;
  int steps = final_time / dt;
  // int steps = 10000;
  // int steps = 1;
  int print_freq = write_freq;
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
  AoSoAType aosoa( "particles", total_no_particles );
  // setup_aosoa(aosoa, &x, &y, &z, &u);
  auto aosoa_position = Cabana::slice<0>( aosoa,    "position");

  // create a mirror in the host space
  auto aosoa_host =
    Cabana::create_mirror_view_and_copy( Kokkos::HostSpace(), aosoa );
  auto aosoa_host_position = Cabana::slice<0>( aosoa_host,    "position");
  auto aosoa_host_ids = Cabana::slice<1>( aosoa_host,    "ids");
  auto aosoa_host_velocity = Cabana::slice<2>( aosoa_host,    "velocity");
  auto aosoa_host_acc = Cabana::slice<3>( aosoa_host,    "acc");
  auto aosoa_host_mass = Cabana::slice<4>( aosoa_host,    "mass");
  auto aosoa_host_density = Cabana::slice<5>( aosoa_host,    "density");
  auto aosoa_host_h = Cabana::slice<6>( aosoa_host,    "h");
  auto aosoa_host_p = Cabana::slice<7>( aosoa_host,    "p");
  auto aosoa_host_is_fluid = Cabana::slice<8>( aosoa_host,    "is_fluid");
  auto aosoa_host_is_boundary = Cabana::slice<9>( aosoa_host,    "is_boundary");
  auto aosoa_host_density_acc = Cabana::slice<10>( aosoa_host,    "density_acc");
  auto aosoa_host_p_acc = Cabana::slice<11>( aosoa_host,    "p_acc");
  auto aosoa_host_sum_wij = Cabana::slice<12>( aosoa_host,    "wij");

  for ( std::size_t i = 0; i < aosoa_host_position.size(); ++i )
    {
      aosoa_host_position ( i, 0 ) = x[i];
      aosoa_host_position ( i, 1 ) = y[i];
      aosoa_host_position ( i, 2 ) = z[i];

      aosoa_host_ids ( i ) = i;

      aosoa_host_velocity ( i, 0 ) = u[i];
      aosoa_host_velocity ( i, 1 ) = 0.;
      aosoa_host_velocity ( i, 2 ) = 0.;

      aosoa_host_acc ( i, 0 ) = 0.;
      aosoa_host_acc ( i, 1 ) = 0.;
      aosoa_host_acc ( i, 2 ) = 0.;

      aosoa_host_mass ( i ) = fluid_rho * pow(fluid_spacing, DIM);
      aosoa_host_density ( i ) = fluid_rho;
      aosoa_host_h ( i ) = 1. * fluid_spacing;

      aosoa_host_p ( i ) = 0.;
      aosoa_host_is_fluid ( i ) = is_fluid[i];
      aosoa_host_is_boundary ( i ) = is_bdry[i];
      aosoa_host_density_acc ( i ) = 0.;
      aosoa_host_p_acc ( i ) = 0.;
      aosoa_host_sum_wij ( i ) = 0.;
    }
  // copy it back to aosoa
  Cabana::deep_copy( aosoa, aosoa_host );

  int fluid_limits[2] = {0, no_fluid_particles};

  std::cout << "Fluid limits is :" << fluid_limits[0] << ", " << fluid_limits[1] << std::endl;
  int bdry_limits[2] = {no_fluid_particles, no_fluid_particles + no_bdry_particles};
  output_data(aosoa, num_particles, 0, time);
  // output_data(aosoa, num_particles, 100, time);
  /*
    ================================================
    End: Step 3
    ================================================
  */

  // ================================================
  // ================================================
  // create the neighbor list
  // ================================================
  // ================================================
  double neighborhood_radius = 4. * fluid_spacing;
  double grid_min[3] = { -0.4, 0.0, -neighborhood_radius };
  double grid_max[3] = { 3.4, 1.6, 0. * neighborhood_radius };
  // double grid_min[3] = { 0.0, -4.0, 0.0 };
  // double grid_max[3] = { 4.1, 4.0, 0.0 };
  double cell_ratio = 1.0;

  ListType verlet_list( aosoa_position, 0,
			aosoa_position.size(), neighborhood_radius,
			cell_ratio, grid_min, grid_max );

  // Main timestep loop
  for ( int step = 0; step < steps; step++ )
    {
      fluid_stage_1(aosoa, dt, fluid_limits);

      continuity_equation(aosoa, dt,
			  &verlet_list,
			  fluid_limits);

      fluid_stage_2(aosoa, dt, fluid_limits);

      state_equation(aosoa, b_rho0_p0, dt, fluid_limits);
      // // Compute the neighbours
      // ListType verlet_list( aosoa_position, 0,
      // 			    aosoa_position.size(), neighborhood_radius,
      // 			    cell_ratio, grid_min, grid_max );

      solid_wall_pressure_bc(aosoa, gravity, dt,
			     &verlet_list,
			     bdry_limits);
      // clamp_wall_pressure(aosoa, dt, bdry_limits);
      // set_wall_velocity(aosoa, dt,
      // 			&verlet_list,
      // 			bdry_limits);
      body_force(aosoa, gravity, dt,
		 &verlet_list,
		 fluid_limits);

      momentum_equation(aosoa, gravity,
			alpha_c0,
			dt,
			&verlet_list,
			fluid_limits);

      fluid_stage_3(aosoa, dt, fluid_limits);

      // output
      if ( step % print_freq == 0 )
	{

	  std::cout << "Time is:" << time << std::endl;
	  output_data(aosoa, num_particles, step, time);
	}

      time += dt;

    }
}

int main( int argc, char* argv[] )
{

  MPI_Init( &argc, &argv );
  Kokkos::initialize( argc, argv );

  // check inputs and write usage
  if ( argc < 8 )
    {
      std::cerr << "Usage: ./DB2d fluid_spacing fluid_height "
	"fluid_length tank_height tank_length t_final write_freq\n";
      std::cerr << "\nwhere fluid_spacing       spacing between the particles"
	"\n";
      std::cerr
	<< "      fluid_height  length of the fluid block\n";
      std::cerr
	<< "      fluid_length  length of the fluid block\n";
      std::cerr
	<< "      tank_height  length of the tank \n";
      std::cerr
	<< "      tank_length  length of the tank \n";

      std::cerr << "      t_final           simulation end time\n";
      std::cerr
	<< "      write_freq      number of steps between output files\n";

      std::cerr << "\nfor example: ./DB2d 0.1 1.0 1.0 1.4 2.0 1 100\n";
      Kokkos::finalize();
      MPI_Finalize();
      return 0;
    }

  // cell size
  double fluid_spacing = std::atof( argv[1] );

  // fluid height
  double fluid_height = std::atof( argv[2] );

  // fluid length
  double fluid_length = std::atof( argv[3] );

  // fluid height
  double tank_height = std::atof( argv[4] );

  // fluid length
  double tank_length = std::atof( argv[5] );

  // end time.
  double t_final = std::atof( argv[6] );

  // write frequency
  int write_freq = std::atoi( argv[7] );

  // // device type
  // std::string device( argv[7] );

  // run the problem.
  damBreak( fluid_spacing, fluid_height, fluid_length,
	    tank_height, tank_length,
	    t_final, write_freq);

  Kokkos::finalize();

  MPI_Finalize();
  return 0;
}
