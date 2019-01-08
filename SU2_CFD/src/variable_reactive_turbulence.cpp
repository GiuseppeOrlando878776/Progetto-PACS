#include "../include/variable_reactive_turbulence.hpp"

CReactiveTurbVariable::CReactiveTurbVariable(su2double val_kine, su2double val_omega, su2double val_muT, unsigned short val_nDim,
                                             unsigned short val_nvar, CConfig *config): CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = config->GetUnsteady_Simulation() == DT_STEPPING_1ST || config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  /*--- Initialization of k -w variables ---*/
  Solution[0] = val_kine;     Solution_Old[0] = val_kine;
  Solution[1] = val_omega;    Solution_Old[1] = val_omega;

  /*--- Initialization of eddy viscosity ---*/
  muT = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_kine;   Solution_time_n[1]  = val_omega;
    Solution_time_n1[0] = val_kine;   Solution_time_n1[1] = val_omega;
  }
}
