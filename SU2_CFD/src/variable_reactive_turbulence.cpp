#include "../include/variable_reactive_turbulence.hpp"

//
//
/*!
 * \brief Class constructor
 */
//
//
CReactiveTurbVariable::CReactiveTurbVariable(su2double val_kine, su2double val_omega, su2double val_muT, unsigned short val_nDim,
                                             unsigned short val_nvar, CConfig* config): CTurbVariable(val_nDim, val_nvar, config) {
  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  /*--- Initialization of k-\omega variables ---*/
  Solution[TURB_KINE_INDEX] = val_kine;       Solution[TURB_OMEGA_INDEX] = val_omega;
  Solution_Old[TURB_KINE_INDEX] = val_kine;   Solution_Old[TURB_OMEGA_INDEX] = val_omega;

  /*--- Initialization of eddy viscosity ---*/
  muT = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[TURB_KINE_INDEX]  = val_kine;   Solution_time_n[TURB_OMEGA_INDEX]  = val_omega;
    Solution_time_n1[TURB_KINE_INDEX] = val_kine;   Solution_time_n1[TURB_OMEGA_INDEX] = val_omega;
  }
}
