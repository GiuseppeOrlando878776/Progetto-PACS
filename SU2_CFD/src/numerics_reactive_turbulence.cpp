#include "../include/numerics_reactive_turbulence.hpp"
#include "../include/variable_reactive_turbulence.hpp"

#include "../../Common/include/su2_assert.hpp"
#include "../include/variable_reactive.hpp"

#include <algorithm>

//
//
/*!
 * \brief Convective term class constructor
 */
//
//
CUpwReactiveTurb::CUpwReactiveTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig* config):
                  CNumerics(val_nDim, val_nVar, config) {
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  nSpecies = CReactiveEulerVariable::GetnSpecies();

  /*--- Allocate auxiliary vectors ---*/
  Velocity_i.resize(nDim);
  Velocity_j.resize(nDim);
}

//
//
/*!
 * \brief Compute residual convective term
 */
//
//
void CUpwReactiveTurb::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) {
  SU2_Assert(val_residual != NULL,"The array of residual for convective flux for turbulence has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive variables at node i for turbulence has not been allocated");
  SU2_Assert(V_j != NULL,"The array of primitive variables at node j for turbulence has not been allocated");

  AD::StartPreacc();
  AD::SetPreaccIn(TurbVar_i,2);
  AD::SetPreaccIn(TurbVar_j,2);
  AD::SetPreaccIn(Normal, nDim);
  if(grid_movement)
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);

  AD::SetPreaccIn(V_i, nSpecies + nDim + 5);
  AD::SetPreaccIn(V_j, nSpecies + nDim + 5);

  Density_i = V_i[CReactiveEulerVariable::RHO_INDEX_PRIM];
  Density_j = V_j[CReactiveEulerVariable::RHO_INDEX_PRIM];

  su2double q_ij = 0.0;
  if(grid_movement) {
    for(unsigned short iDim = 0; iDim < nDim; ++iDim) {
      Velocity_i[iDim] = V_i[CReactiveEulerVariable::VX_INDEX_PRIM + iDim] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[CReactiveEulerVariable::VX_INDEX_PRIM + iDim] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim] + Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    std::copy(V_i + CReactiveEulerVariable::VX_INDEX_PRIM, V_i + CReactiveEulerVariable::VX_INDEX_PRIM + nDim, Velocity_i.begin());
    std::copy(V_j + CReactiveEulerVariable::VX_INDEX_PRIM, V_j + CReactiveEulerVariable::VX_INDEX_PRIM + nDim, Velocity_j.begin());
    for(unsigned short iDim = 0; iDim < nDim; ++iDim)
      q_ij += 0.5*(Velocity_i[iDim] + Velocity_j[iDim])*Normal[iDim];
  }

  su2double a0 = 0.5*(q_ij + std::abs(q_ij));
  su2double a1 = 0.5*(q_ij - std::abs(q_ij));

  val_residual[CReactiveTurbVariable::TURB_KINE_INDEX] = a0*Density_i*TurbVar_i[CReactiveTurbVariable::TURB_KINE_INDEX] +
                                                         a1*Density_j*TurbVar_j[CReactiveTurbVariable::TURB_KINE_INDEX];
  val_residual[CReactiveTurbVariable::TURB_OMEGA_INDEX] = a0*Density_i*TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX] +
                                                          a1*Density_j*TurbVar_j[CReactiveTurbVariable::TURB_OMEGA_INDEX];

  if(implicit) {
    SU2_Assert(val_Jacobian_i != NULL,"The matrix of Jacobian at node i for turbulence has not been allocated");
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      SU2_Assert(val_Jacobian_i[iVar] != NULL,"The matrix of Jacobian at node i for turbulence has not been allocated");
    SU2_Assert(val_Jacobian_j != NULL,"The matrix of Jacobian at node i at node j for turbulence has not been allocated");
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      SU2_Assert(val_Jacobian_j[iVar] != NULL,"The matrix of Jacobian at node j for turbulence has not been allocated");

    val_Jacobian_i[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = a0;
    val_Jacobian_i[CReactiveTurbVariable::TURB_KINE_INDEX][1] = 0.0;
    val_Jacobian_i[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = 0.0;
    val_Jacobian_i[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = a0;

    val_Jacobian_j[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = a1;
    val_Jacobian_j[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = 0.0;
    val_Jacobian_j[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = 0.0;
    val_Jacobian_j[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = a1;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

//
//
/*!
 * \brief Viscous term class constructor
 */
//
//
CAvgGrad_ReactiveTurb::CAvgGrad_ReactiveTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig* config):
                       CNumerics(val_nDim, val_nVar, config) {
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  nSpecies = CReactiveEulerVariable::GetnSpecies();

  /*--- Allocate auxiliary vectors and matrices ---*/
  Edge_Vector.resize(nDim);
  Mean_GradTurbVar = new su2double* [nVar];
  for(unsigned short iVar = 0; iVar < nVar; ++iVar)
    Mean_GradTurbVar[iVar] = new su2double[nDim];
  Proj_Mean_GradTurbVar_Normal.resize(nVar);
}

//
//
/*!
 * \brief Viscous term class destructor
 */
//
//
CAvgGrad_ReactiveTurb::~CAvgGrad_ReactiveTurb() {
  if(Mean_GradTurbVar != NULL) {
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      delete[] Mean_GradTurbVar[iVar];
    delete[] Mean_GradTurbVar;
  }
}


//
//
/*!
 * \brief Compute residual viscous term for turbulence
 */
//
//
void CAvgGrad_ReactiveTurb::ComputeResidual(su2double* val_residual, su2double** Jacobian_i, su2double** Jacobian_j, CConfig* config) {
  SU2_Assert(val_residual != NULL,"The array of residual for convective flux for turbulence has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive variables at node i for turbulence has not been allocated");
  SU2_Assert(V_j != NULL,"The array of primitive variables at node j for turbulence has not been allocated");

  unsigned short iDim, iVar;
  su2double dist_ij_2, proj_vector_ij;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega, diff_kine, diff_omega;

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(V_i, nSpecies + nDim + 5); AD::SetPreaccIn(V_j, nSpecies + nDim + 5);

  Density_i = V_i[CReactiveEulerVariable::RHO_INDEX_PRIM];
  Density_j = V_j[CReactiveEulerVariable::RHO_INDEX_PRIM];

  /*--- Compute mean effective viscosity ---*/
  const su2double omega = 0.5;
  const su2double omega_star = 0.6;
  diff_i_kine  = Laminar_Viscosity_i + omega*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + omega*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + omega_star*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + omega_star*Eddy_Viscosity_j;

  diff_kine  = 2.0*(1.0/diff_i_kine + 1.0/diff_j_kine);
  diff_omega = 2.0*(1.0/diff_i_omega + 1.0/diff_j_omega);

  /*--- Compute vector going from iPoint to jPoint ---*/
  std::transform(Coord_j, Coord_j + nDim, Coord_i, Edge_Vector.begin(), std::less<su2double>());
  dist_ij_2 = std::inner_product(Edge_Vector.cbegin(), Edge_Vector.cend(), Edge_Vector.cbegin(), 0.0);
  proj_vector_ij = std::inner_product(Edge_Vector.cbegin(), Edge_Vector.cend(), Normal, 0.0);
  if(dist_ij_2 < EPS)
    proj_vector_ij = 0.0;
  else
    proj_vector_ij /= dist_ij_2;

  /*--- Mean gradient approximation and projection of the mean gradient in the direction of the edge ---*/
  for(iVar = 0; iVar < nVar; ++iVar) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    for(iDim = 0; iDim < nDim; ++iDim) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }

  val_residual[CReactiveTurbVariable::TURB_KINE_INDEX] = diff_kine*Proj_Mean_GradTurbVar_Normal[CReactiveTurbVariable::TURB_KINE_INDEX];
  val_residual[CReactiveTurbVariable::TURB_OMEGA_INDEX] = diff_omega*Proj_Mean_GradTurbVar_Normal[CReactiveTurbVariable::TURB_OMEGA_INDEX];

  /*--- For Jacobians -> Use of TSL(Thin Shear Layer) approximation to compute derivatives of the gradients ---*/
  if(implicit) {
    SU2_Assert(Jacobian_i != NULL,"The matrix of Jacobian at node i for turbulence has not been allocated");
    for(iVar = 0; iVar < nVar; ++iVar)
      SU2_Assert(Jacobian_i[iVar] != NULL,"The matrix of Jacobian at node i for turbulence has not been allocated");
    SU2_Assert(Jacobian_j != NULL,"The matrix of Jacobian at node i at node j for turbulence has not been allocated");
    for(iVar = 0; iVar < nVar; ++iVar)
      SU2_Assert(Jacobian_j[iVar] != NULL,"The matrix of Jacobian at node j for turbulence has not been allocated");

    Jacobian_i[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = -diff_kine*proj_vector_ij/Density_i;
    Jacobian_i[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = 0.0;
    Jacobian_i[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = 0.0;
    Jacobian_i[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = -diff_omega*proj_vector_ij/Density_i;

    Jacobian_j[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = diff_kine*proj_vector_ij/Density_j;
    Jacobian_j[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = 0.0;
    Jacobian_j[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = 0.0;
    Jacobian_j[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = diff_omega*proj_vector_ij/Density_j;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

//
//
/*!
 * \brief Source term class constructor
 */
//
//
CSourcePieceWise_ReactiveTurb::CSourcePieceWise_ReactiveTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig* config):
                               CNumerics(val_nDim, val_nVar, config) {
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  nSpecies = CReactiveEulerVariable::GetnSpecies();
}

//
//
/*!
 * \brief Compute residual source term
 */
//
//
void CSourcePieceWise_ReactiveTurb::ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i,
                                                    su2double** val_Jacobian_j, CConfig* config) {
  SU2_Assert(val_residual != NULL,"The array of residual for convective flux for turbulence has not been allocated");
  SU2_Assert(V_i != NULL,"The array of primitive variables at node i for turbulence has not been allocated");

  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume);
  AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim + 1, nDim);
  AD::SetPreaccIn(V_i, nSpecies + nDim + 5);

  unsigned short iDim;
  su2double diverg, pk, pw;
  const su2double alfa = 13.0/25.0;
  const su2double beta  = 0.5;
  const su2double beta_star  = 0.09;

  Density_i = V_i[CReactiveEulerVariable::RHO_INDEX_PRIM];

  /*--- Computation of blended constants for the source terms---*/
  if(dist_i > 1.0e-10) {
    /*--- Production ---*/
    diverg = 0.0;
    for(iDim = 0; iDim < nDim; ++iDim)
      diverg += PrimVar_Grad_i[CReactiveEulerVariable::VX_INDEX_GRAD + iDim][iDim];

    pk = Eddy_Viscosity_i*StrainMag_i*StrainMag_i - TWO3*Density_i*TurbVar_i[CReactiveTurbVariable::TURB_KINE_INDEX]*diverg;
    pw = StrainMag_i*StrainMag_i - TWO3*TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX]*diverg;

    val_residual[CReactiveTurbVariable::TURB_KINE_INDEX] += pk*Volume;
    val_residual[CReactiveTurbVariable::TURB_OMEGA_INDEX] += alfa*Density_i*pw*Volume;

    /*--- Dissipation ---*/
    val_residual[CReactiveTurbVariable::TURB_KINE_INDEX] -= beta_star*Density_i*TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX]*
                                                            TurbVar_i[CReactiveTurbVariable::TURB_KINE_INDEX]*Volume;
    val_residual[CReactiveTurbVariable::TURB_OMEGA_INDEX] -= beta*Density_i*TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX]*
                                                             TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX]*Volume;

    /*--- Implicit part ---*/
    if(implicit) {
      SU2_Assert(val_Jacobian_i != NULL,"The matrix of Jacobian at node i for turbulence has not been allocated");
      for(unsigned short iVar = 0; iVar < nVar; ++iVar)
        SU2_Assert(val_Jacobian_i[iVar] != NULL,"The matrix of Jacobian at node i for turbulence has not been allocated");

      val_Jacobian_i[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] =
      -beta_star*TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX]*Volume;

      val_Jacobian_i[CReactiveTurbVariable::TURB_KINE_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] = 0.0;

      val_Jacobian_i[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_KINE_INDEX] = 0.0;
      val_Jacobian_i[CReactiveTurbVariable::TURB_OMEGA_INDEX][CReactiveTurbVariable::TURB_OMEGA_INDEX] =
      -2.0*beta*TurbVar_i[CReactiveTurbVariable::TURB_OMEGA_INDEX]*Volume;
    }
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}
