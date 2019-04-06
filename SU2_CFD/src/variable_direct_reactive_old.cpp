#include "../include/variable_reactive.hpp"
#include "../../Common/include/reacting_model_library.hpp"

unsigned short CReactiveEulerVariable::nSpecies = 0;
CReactiveEulerVariable::LibraryPtr CReactiveEulerVariable::library = NULL;

//
//
/*!
 *\brief Class default constructor
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable():CVariable() {
  nPrimVar = 0;
  nPrimVarGrad = 0;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nPrimVarLim = 0;
}

//
//
/*!
 *\brief Class constructor
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                               unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                               CConfig* config): CVariable(val_nDim,val_nvar,config), Cp() {
  nSpecies = val_nSpecies;
  nPrimVar = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;
  nSecondaryVar = 0;
  nSecondaryVarGrad = 0;
  nPrimVarLim = val_nprimvarlim;

  library = LibraryPtr(new Framework::ReactingModelLibrary(config->GetConfigLibFile()));

  Primitive.resize(nPrimVar);
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for(unsigned short iVar = 0; iVar < nPrimVarGrad; ++iVar)
    Gradient_Primitive[iVar] = new su2double[nDim];
  Limiter_Primitive.resize(nPrimVarLim);

  Limiter = new su2double [nVar];
  Solution_Max = new su2double[nPrimVarLim];
  Solution_Min = new su2double[nPrimVarLim];

  Ys.resize(nSpecies);
}

//
//
/*!
 *\brief Class overloaded constructor (initialization values)
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                                               const su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, CConfig* config):
                                               CReactiveEulerVariable(val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                      val_nprimvarlim,config) {
  /*--- Rename and initialize for convenience ---*/
  su2double T = val_temperature;   // Translational-rotational temperature [K]
  su2double P = val_pressure;

  su2double rho,rhoE;

  /*--- Compute mixture density ---*/
  //rho = library->ComputeDensity(T,P,val_massfrac);
  rho *= config->GetGas_Constant_Ref();

  su2double dim_temp = T*config->GetTemperature_Ref();
  bool US_System = (config->GetSystemMeasurements() == US);
  if(US_System)
    dim_temp *= 5.0/9.0;

  /*--- Compute energy (RHOE) from supplied primitive quanitites ---*/
  /*
  su2double sqvel = std::inner_product(val_velocity.cbegin(), val_velocity.cend(), val_velocity.cbegin());
  su2double e_tot = library->ComputeEnergy(dim_temp,val_massfrac)/config->GetEnergy_Ref();
  rhoE = rho*(0.5*sqvel + e_tot);
  */

  /*--- Initialize Solution and Solution_Old vectors ---*/
  /*--- Initialize mixture density and partial density ---*/
  Solution[RHO_INDEX_SOL] = Solution_Old[RHO_INDEX_SOL] = rho;
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Solution[RHOS_INDEX_SOL + iSpecies] = Solution_Old[RHOS_INDEX_SOL + iSpecies] = rho*val_massfrac[iSpecies];

  /*--- Initialize momentum ---*/
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    Solution[RHOVX_INDEX_SOL + iDim] = Solution_Old[RHOVX_INDEX_SOL + iDim] = rho*val_velocity[iDim];

  /*--- Initialize energy contribution ---*/
  Solution[RHOE_INDEX_SOL] = Solution_Old[RHOE_INDEX_SOL] = rhoE;

  /*--- Assign primitive variables ---*/
  Primitive.at(T_INDEX_PRIM) = T;
  Primitive.at(P_INDEX_PRIM) = P;

  /*--- Set specific heat at constant pressure ---*/
  //Cp = library->ComputeCP(dim_temp,val_massfrac);
  Cp /= config->GetEnergy_Ref()*config->GetTemperature_Ref();
  if(US_System)
    Cp *= 3.28084*3.28084*9.0/5.0;
}

//
//
/*!
 *\brief Class overloaded constructor (initialization vector)
 */
//
//
CReactiveEulerVariable::CReactiveEulerVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                               unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                               unsigned short val_nprimvarlim, CConfig* config):
                                               CReactiveEulerVariable(val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                      val_nprimvarlim,config) {
  /*--- Initialize Solution and Solution_Old vectors ---*/
  SU2_Assert(Solution != NULL,"The array Solution has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array Solution_Old has not been allocated");

  std::copy(val_solution.cbegin(),val_solution.cend(),Solution);
  std::copy(val_solution.cbegin(),val_solution.cend(),Solution_Old);

  /*--- Initialize T and P to the free stream for Secant method ---*/
  Primitive.at(T_INDEX_PRIM) = config->GetTemperature_FreeStream();
  Primitive.at(P_INDEX_PRIM) = config->GetPressure_FreeStream();

  /*--- Set specific heat at constant pressure ---*/
  su2double dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Ys[iSpecies] = Solution[RHOS_INDEX_SOL]/Solution[RHO_INDEX_SOL];
  bool US_System = (config->GetSystemMeasurements() == US);
  if(US_System)
    dim_temp *= 5.0/9.0;
  //Cp = library->ComputeCP(dim_temp,Ys);
  Cp /= config->GetEnergy_Ref()*config->GetTemperature_Ref();
  if(US_System)
    Cp *= 3.28084*3.28084*9.0/5.0;
}

//
//
/*!
 *\brief Class destructor
 */
//
//
CReactiveEulerVariable::~CReactiveEulerVariable() {
  if(Gradient_Primitive != NULL) {
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      delete[] Gradient_Primitive[iVar];
    delete[] Gradient_Primitive;
  }
}

//
//
/*!
 *\brief Set primitive variables
 */
//
//
bool CReactiveEulerVariable::SetPrimVar(CConfig* config) {
  /*--- Convert conserved to primitive variables ---*/
  bool nonPhys = Cons2PrimVar(config, Solution, Primitive.data());
  if(nonPhys) {
    std::copy(Solution_Old,Solution_Old + nVar,Solution);
    bool check_old = Cons2PrimVar(config, Solution, Primitive.data());
    SU2_Assert(check_old == true, "Neither the old solution is feasible to set primitive variables: problem unsolvable");
  }

  /*--- Set specific heat at constant pressure ---*/
  su2double dim_temp = Primitive.at(T_INDEX_PRIM)*config->GetTemperature_Ref();
  bool US_System = (config->GetSystemMeasurements() == US);
  if(US_System)
    dim_temp *= 5.0/9.0;
  //Cp = library->ComputeCP(dim_temp,GetMassFractions());
  Cp /= config->GetEnergy_Ref()*config->GetTemperature_Ref();
  if(US_System)
    Cp *= 3.28084*3.28084*9.0/5.0;

  return nonPhys;
}

//
//
/*!
 *\brief Pass from conserved to primitive variables
 */
//
//
bool CReactiveEulerVariable::Cons2PrimVar(CConfig* config, su2double* U, su2double* V) {
  SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
  SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");

  bool NRconvg, Bconvg, nonPhys;
	unsigned short iDim, iSpecies, iIter, maxBIter, maxNIter;
  su2double rho, rhoE;
  su2double sqvel;
  su2double f, df, NRtol, Btol;
  su2double T, Told, Tnew, hs, hs_old, Tmin, Tmax;

  /*--- Conserved and primitive vector layout ---*/
  // U:  [rho, rhou, rhov, rhow, rhoE, rho1, ..., rhoNs]^T
  // V: [T, u, v, w, P, rho, h, a, rho1, ..., rhoNs,]^T

  /*--- Set booleans ---*/
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0;
  Tmax   = 6.0e4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0e-6;    // Tolerance for the Secant method
  Btol     = 1.0e-4;    // Tolerance for the Bisection method
  maxNIter = 7;        // Maximum Secant method iterations
  maxBIter = 32;        // Maximum Bisection method iterations

  /*--- Rename variables forconvenience ---*/
  rhoE  = U[RHOE_INDEX_SOL];          // Density * total energy [J/m3]

  /*--- Assign species mass fraction and mixture density ---*/
  // NOTE: If any species densities are < 0, these values are re-assigned
  //       in the conserved vector to ensure positive density
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(U[RHOS_INDEX_SOL + iSpecies] < EPS) {
      U[RHOS_INDEX_SOL + iSpecies] = EPS;
      nonPhys = true;
    }
  }

  if(U[RHO_INDEX_SOL] < EPS) {
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL] = EPS;
    nonPhys = true;
  }
  else
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL];

  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    V[RHOS_INDEX_PRIM + iSpecies] = U[RHOS_INDEX_SOL]/U[RHO_INDEX_SOL];

  /*--- Checking sum of mass fraction ---*/
  std::copy(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  nonPhys = nonPhys || (std::abs(std::accumulate(Ys.cbegin(),Ys.cend(),0.0) - 1.0) > EPS);

  /*--- Rename for convenience ---*/
  rho = U[RHO_INDEX_SOL];

  /*--- Assign mixture velocity ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
  sqvel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Translational-Rotational Temperature ---*/
  const su2double Rgas = library->ComputeRgas(Ys)/config->GetGas_Constant_Ref();
  const su2double C1 = (-rhoE + rho*sqvel)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  NRconvg = false;
  bool US_System = (config->GetSystemMeasurements() == US);
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    /*--- Execute a secant root-finding method to find T ---*/
    su2double dim_temp, dim_temp_old;
    dim_temp = T*config->GetTemperature_Ref();
    dim_temp_old = Told*config->GetTemperature_Ref();
    if(US_System) {
      dim_temp *= 5.0/9.0;
      dim_temp_old *= 5.0/9.0;
    }
    //hs_old = library->ComputeEnthalpy(dim_temp_old,Ys)/config->GetEnergy_Ref();
    //hs = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
    if(US_System) {
      hs_old *= 3.28084*3.28084;
      hs *= 3.28084*3.28084;
    }

    f = T - C1 - C2*hs;
    df = T - Told + C2*(hs_old-hs);
    Tnew = T - f*(T-Told)/df;

    /*--- Check for convergence ---*/
    if(std::abs(Tnew - T) < NRtol) {
      NRconvg = true;
      break;
    }
    else {
      Told = T;
      T = Tnew;
    }
  }

  /*--- If the secant method has converged, assign the value of T.
        Otherwise execute a bisection root-finding method ---*/
  if(NRconvg)
    V[T_INDEX_PRIM] = T;
  else {
    Bconvg = false;
    su2double Ta = Tmin;
    su2double Tb = Tmax;
    for(iIter = 0; iIter < maxBIter; ++iIter) {
      T = (Ta + Tb)/2.0;
      su2double dim_temp = T*config->GetTemperature_Ref();;
      if(US_System)
        dim_temp *= 5.0/9.0;
      //hs = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
      if(US_System)
        hs *= 3.28084*3.28084;
      f = T - C1 - C2*hs;

      if(std::abs(f) < Btol) {
        V[T_INDEX_PRIM] = T;
        Bconvg = true;
        break;
      }
      else {
        if(f > 0)
          Ta = T;
        else
          Tb = T;
      }

    /*--- If absolutely no convergence, then something is going really wrong ---*/
    if(!Bconvg)
      throw std::runtime_error("Convergence not achieved for bisection method");
    }
  }

  if(V[T_INDEX_PRIM] < Tmin) {
    V[T_INDEX_PRIM] = Tmin;
    nonPhys = true;
  }
  else if(V[T_INDEX_PRIM] > Tmax) {
    V[T_INDEX_PRIM] = Tmax;
    nonPhys = true;
  }

  /*--- Pressure ---*/
  //V[P_INDEX_PRIM] = library->ComputePressure(V[RHO_INDEX_PRIM],V[T_INDEX_PRIM],Ys)/config->GetGas_Constant_Ref();
  if(V[P_INDEX_PRIM] < EPS) {
    V[P_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Sound speed ---*/
  su2double dim_temp = V[T_INDEX_PRIM]*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  //V[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp,Ys,V[P_INDEX_PRIM],V[RHO_INDEX_PRIM]);
  if(V[A_INDEX_PRIM] < EPS) {
    V[A_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Enthalpy ---*/
  V[H_INDEX_PRIM] = (U[RHOE_INDEX_SOL] + V[P_INDEX_PRIM])/rho;

  return nonPhys;
}

//
//
/*!
 *\brief Pass from primitive to conserved variables
 */
//
//
void CReactiveEulerVariable::Prim2ConsVar(CConfig* config, su2double* V, su2double* U) {
  SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
  SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");

  /*--- Set mixture density and species density ---*/
  U[RHO_INDEX_SOL] = V[RHO_INDEX_PRIM];
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    U[RHOS_INDEX_SOL + iSpecies] = V[RHO_INDEX_PRIM]*V[RHOS_INDEX_PRIM];

  /*--- Set momentum ---*/
  for(unsigned short iDim = 0; iDim < nDim; ++iDim)
    U[RHOVX_INDEX_SOL + iDim] = V[RHO_INDEX_PRIM]*V[VX_INDEX_PRIM + iDim];

  /*--- Set energies ---*/
  U[RHOE_INDEX_SOL] = V[RHO_INDEX_PRIM]*V[H_INDEX_PRIM] - V[P_INDEX_PRIM];
}

//
//
/*!
 *\brief Set density
 */
//
//
inline bool CReactiveEulerVariable::SetDensity(void) {
  SU2_Assert(Solution != NULL,"The array of solution variables has not been allocated");
  Primitive.at(RHO_INDEX_PRIM) = Solution[RHO_INDEX_SOL];

  if(Primitive.at(RHO_INDEX_PRIM) < EPS)
    return true;

  return false;
}

//
//
/*!
 *\brief Set pressure (NOTE: it requires SetDensity() call)
 */
//
//
bool CReactiveEulerVariable::SetPressure(CConfig* config) {
  /*--- Compute this gas mixture property from library ---*/
  su2double Pressure;
  //su2double Pressure = library->ComputePressure(Primitive.at(T_INDEX_PRIM),Primtive.at(RHO_INDEX_PRIM), GetMassFractions());
  //Pressure /= config->GetGas_Constant_Ref();

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(P_INDEX_PRIM) = Pressure;
  if(Pressure < EPS)
    return true;

  return false;
}

//
//
/*!
 *\brief Set sound speed (requires SetDensity() call)
 */
//
//
bool CReactiveEulerVariable::SetSoundSpeed(CConfig* config) {
  su2double dim_temp = Primitive.at(T_INDEX_PRIM)*config->GetTemperature_Ref();
  if(config->GetSystemMeasurements() == US)
    dim_temp *= 5.0/9.0;
  su2double Sound_Speed;
  //su2double Sound_Speed = library->ComputeFrozenSoundSpeed(dim_temp, GetMassFractions(),
  //                                                         Primitive.at(P_INDEX_PRIM), Primitive.at(RHO_INDEX_PRIM));

  /*--- Store computed value and check for a physical solution ---*/
  Primitive.at(A_INDEX_PRIM) = Sound_Speed;
  if(Sound_Speed < EPS)
    return true;

  return false;
}

//
//
/*!
 *\brief Set the velocity vector from old solution
 */
//
//
inline void CReactiveEulerVariable::SetVelocity_Old(su2double* val_velocity) {
  SU2_Assert(val_velocity != NULL,"The array of velocity val_velocity has not been allocated");
  SU2_Assert(Solution_Old != NULL,"The array of solution variables at previous step has not been allocated");

  for(unsigned short iDim=0; iDim < nDim; ++iDim)
    Solution_Old[RHOVX_INDEX_SOL+iDim] = val_velocity[iDim]*Primitive.at(RHO_INDEX_PRIM);
}

//
//
/*!
 *\brief Class constructor
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(unsigned short val_nDim, unsigned short val_nvar, unsigned short val_nSpecies,
                                         unsigned short val_nprimvar, unsigned short val_nprimvargrad, unsigned short val_nprimvarlim,
                                         CConfig* config):
                                         CReactiveEulerVariable(val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                val_nprimvarlim,config) {
  Laminar_Viscosity = 0.0;
  Thermal_Conductivity = 0.0;
  Diffusion_Coeffs.resize(nSpecies,nSpecies);

  /*--- Turbulence part ---*/
  Eddy_Viscosity = 0.0;
  Vorticity = {};
  StrainMag = 0.0;
}

//
//
/*!
 *\brief Class overloaded constructor (initialization values)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(const su2double val_pressure, const RealVec& val_massfrac, const RealVec& val_velocity,
                                         const su2double val_temperature, unsigned short val_nDim, unsigned short val_nvar,
                                         unsigned short val_nSpecies,unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                         unsigned short val_nprimvarlim, CConfig* config):
                                         CReactiveEulerVariable(val_pressure,val_massfrac,val_velocity,val_temperature,val_nDim,val_nvar,
                                                                val_nSpecies,val_nprimvar,val_nprimvargrad,val_nprimvarlim,config) {
  su2double dim_temp, dim_press;
  bool US_System = (config->GetSystemMeasurements() == US);
  dim_temp = val_temperature*config->GetTemperature_Ref();
  dim_press = val_pressure*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }

  /*--- Compute transport properties --- */
  //Laminar_Viscosity = library->GetLambda(dim_temp,val_massfrac)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,val_massfrac)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);

  /*--- Turbulence part ---*/
  Eddy_Viscosity = 0.0;
  Vorticity = {};
  StrainMag = 0.0;
}

//
//
/*!
 *\brief Class overloaded constructor (initialization vector)
 */
//
//
CReactiveNSVariable::CReactiveNSVariable(const RealVec& val_solution, unsigned short val_nDim, unsigned short val_nvar,
                                         unsigned short val_nSpecies, unsigned short val_nprimvar, unsigned short val_nprimvargrad,
                                         unsigned short val_nprimvarlim, CConfig* config):
                                         CReactiveEulerVariable(val_solution,val_nDim,val_nvar,val_nSpecies,val_nprimvar,val_nprimvargrad,
                                                                val_nprimvarlim,config) {
  su2double dim_temp, dim_press;
  bool US_System = (config->GetSystemMeasurements() == US);
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }

  /*--- Compute transport properties --- */
  Ys = GetMassFractions();
  //Laminar_Viscosity = library->GetLambda(dim_temp,Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);

  /*--- Turbulence part ---*/
  Eddy_Viscosity = 0.0;
  Vorticity = {};
  StrainMag = 0.0;
}

//
//
/*!
 *\brief Set primitive variables
 */
//
//
bool CReactiveNSVariable::SetPrimVar(CConfig* config) {
  /*--- Convert conserved to primitive variables using Euler version since primitives are the same ---*/
  bool nonPhys = CReactiveEulerVariable::SetPrimVar(config);

  su2double dim_temp, dim_press;
  bool US_System = (config->GetSystemMeasurements() == US);
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }

  /*--- Compute transport properties --- */
  Ys = GetMassFractions();
  //Laminar_Viscosity = library->GetLambda(dim_temp,Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);

  return nonPhys;
}

//
//
/*!
 *\brief Pass from conserved to primitive variables including turbulence effects
 */
//
//
bool CReactiveNSVariable::Cons2PrimVar(CConfig* config, su2double* U, su2double* V, su2double turb_ke) {
  SU2_Assert(U != NULL,"The array of conserved varaibles has not been allocated");
  SU2_Assert(V != NULL,"The array of primitive varaibles has not been allocated");

  bool NRconvg, Bconvg, nonPhys;
	unsigned short iDim, iSpecies, iIter, maxBIter, maxNIter;
  su2double rho, rhoE;
  su2double sqvel;
  su2double f, df, NRtol, Btol;
  su2double T, Told, Tnew, hs, hs_old, Tmin, Tmax;

  /*--- Conserved and primitive vector layout ---*/
  // U:  [rho, rhou, rhov, rhow, rhoE, rho1, ..., rhoNs]^T
  // V: [T, u, v, w, P, rho, h, a, rho1, ..., rhoNs,]^T

  /*--- Set booleans ---*/
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0;
  Tmax   = 6.0e4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0e-6;    // Tolerance for the Secant method
  Btol     = 1.0e-4;    // Tolerance for the Bisection method
  maxNIter = 7;        // Maximum Secant method iterations
  maxBIter = 32;        // Maximum Bisection method iterations

  /*--- Rename variables forconvenience ---*/
  rhoE  = U[RHOE_INDEX_SOL];          // Density * total energy [J/m3]

  /*--- Assign species mass fraction and mixture density ---*/
  // NOTE: If any species densities are < 0, these values are re-assigned
  //       in the conserved vector to ensure positive density
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(U[RHOS_INDEX_SOL + iSpecies] < EPS) {
      U[RHOS_INDEX_SOL + iSpecies] = EPS;
      nonPhys = true;
    }
  }

  if(U[RHO_INDEX_SOL] < EPS) {
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL] = EPS;
    nonPhys = true;
  }
  else
    V[RHO_INDEX_PRIM] = U[RHO_INDEX_SOL];

  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    V[RHOS_INDEX_PRIM + iSpecies] = U[RHOS_INDEX_SOL]/U[RHO_INDEX_SOL];

  /*--- Checking sum of mass fraction ---*/
  std::copy(V + RHOS_INDEX_PRIM, V + (RHOS_INDEX_PRIM + nSpecies), Ys.begin());
  nonPhys = nonPhys || (std::abs(std::accumulate(Ys.cbegin(),Ys.cend(),0.0) - 1.0) > EPS);

  /*--- Rename for convenience ---*/
  rho = U[RHO_INDEX_SOL];

  /*--- Assign mixture velocity ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    V[VX_INDEX_PRIM + iDim] = U[RHOVX_INDEX_SOL + iDim]/rho;
  sqvel = std::inner_product(V + VX_INDEX_PRIM, V + (VX_INDEX_PRIM + nDim), V + VX_INDEX_PRIM, 0.0);

  /*--- Translational-Rotational Temperature ---*/
  const su2double Rgas = library->ComputeRgas(Ys)/config->GetGas_Constant_Ref();
  const su2double C1 = (-rhoE + rho*sqvel + rho*turb_ke)/(rho*Rgas);
  const su2double C2 = 1.0/Rgas;

  T = V[T_INDEX_PRIM];
  Told = T + 1.0;
  NRconvg = false;
  bool US_System = (config->GetSystemMeasurements() == US);
  for(iIter = 0; iIter < maxNIter; ++iIter) {
    /*--- Execute a secant root-finding method to find T ---*/
    su2double dim_temp, dim_temp_old;
    dim_temp = T*config->GetTemperature_Ref();
    dim_temp_old = Told*config->GetTemperature_Ref();
    if(US_System) {
      dim_temp *= 5.0/9.0;
      dim_temp_old *= 5.0/9.0;
    }
    //hs_old = library->ComputeEnthalpy(dim_temp_old,Ys)/config->GetEnergy_Ref();
    //hs = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
    if(US_System) {
      hs_old *= 3.28084*3.28084;
      hs *= 3.28084*3.28084;
    }

    f = T - C1 - C2*hs;
    df = T - Told + C2*(hs_old-hs);
    Tnew = T - f*(T-Told)/df;

    /*--- Check for convergence ---*/
    if(std::abs(Tnew - T) < NRtol) {
      NRconvg = true;
      break;
    }
    else {
      Told = T;
      T = Tnew;
    }
  }

  /*--- If the secant method has converged, assign the value of T.
        Otherwise execute a bisection root-finding method ---*/
  if(NRconvg)
    V[T_INDEX_PRIM] = T;
  else {
    Bconvg = false;
    su2double Ta = Tmin;
    su2double Tb = Tmax;
    for(iIter = 0; iIter < maxBIter; ++iIter) {
      T = (Ta + Tb)/2.0;
      su2double dim_temp = T*config->GetTemperature_Ref();;
      if(US_System)
        dim_temp *= 5.0/9.0;
      //hs = library->ComputeEnthalpy(dim_temp,Ys)/config->GetEnergy_Ref();
      if(US_System)
        hs *= 3.28084*3.28084;
      f = T - C1 - C2*hs;

      if(std::abs(f) < Btol) {
        V[T_INDEX_PRIM] = T;
        Bconvg = true;
        break;
      }
      else {
        if(f > 0)
          Ta = T;
        else
          Tb = T;
      }

    /*--- If absolutely no convergence, then something is going really wrong ---*/
    if(!Bconvg)
      throw std::runtime_error("Convergence not achieved for bisection method");
    }
  }

  if(V[T_INDEX_PRIM] < Tmin) {
    V[T_INDEX_PRIM] = Tmin;
    nonPhys = true;
  }
  else if(V[T_INDEX_PRIM] > Tmax) {
    V[T_INDEX_PRIM] = Tmax;
    nonPhys = true;
  }

  /*--- Pressure ---*/
  //V[P_INDEX_PRIM] = library->ComputePressure(V[RHO_INDEX_PRIM],V[T_INDEX_PRIM],Ys)/config->GetGas_Constant_Ref();
  if(V[P_INDEX_PRIM] < EPS) {
    V[P_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Sound speed ---*/
  su2double dim_temp = V[T_INDEX_PRIM]*config->GetTemperature_Ref();
  if(US_System)
    dim_temp *= 5.0/9.0;
  //V[A_INDEX_PRIM] = library->ComputeFrozenSoundSpeed(dim_temp,Ys,V[P_INDEX_PRIM],V[RHO_INDEX_PRIM]);
  if(V[A_INDEX_PRIM] < EPS) {
    V[A_INDEX_PRIM] = EPS;
    nonPhys = true;
  }

  /*--- Enthalpy ---*/
  V[H_INDEX_PRIM] = (U[RHOE_INDEX_SOL] + V[P_INDEX_PRIM])/rho;

  return nonPhys;
}

//
//
/*!
 *\brief Set primitive variables including turbu variables
 */
//
//
bool CReactiveNSVariable::SetPrimVar(su2double eddy_visc, su2double turb_ke, CConfig* config) {
  /*--- Convert conserved to primitive variables including turbulence effect---*/
  bool nonPhys = Cons2PrimVar(config, Solution, Primitive.data(), turb_ke);
  if(nonPhys) {
    std::copy(Solution_Old,Solution_Old + nVar,Solution);
    bool check_old = Cons2PrimVar(config, Solution, Primitive.data(), turb_ke);
    SU2_Assert(check_old == true, "Neither the old solution is feasible to set primitive variables: problem unsolvable");
  }

  su2double dim_temp, dim_press;
  bool US_System = (config->GetSystemMeasurements() == US);
  dim_temp = Primitive[T_INDEX_PRIM]*config->GetTemperature_Ref();
  dim_press = Primitive[P_INDEX_PRIM]*config->GetPressure_Ref()/101325.0;
  if(US_System) {
    dim_temp *= 5.0/9.0;
    dim_press *= 47.8803;
  }

  Ys = GetMassFractions();

  /*--- Compute specific heat at constant pressure --- */
  //Cp = library->ComputeCP(dim_temp,Ys);
  Cp /= config->GetEnergy_Ref()*config->GetTemperature_Ref();
  if(US_System)
    Cp *= 3.28084*3.28084*9.0/5.0;

  /*--- Compute transport properties --- */
  //Laminar_Viscosity = library->GetLambda(dim_temp,Ys)/config->GetViscosity_Ref();
  if(US_System)
    Laminar_Viscosity *= 0.02088553108;
  //Thermal_Conductivity = library->GetEta(dim_temp,Ys)/config->GetConductivity_Ref();
  if(US_System)
    Thermal_Conductivity *= 0.12489444444;
  //Diffusion_Coeffs = library->GetDij_SM(dim_temp,dim_press)/(config->GetVelocity_Ref()*config->GetLength_Ref()*1e-4);

  /*--- Set eddy viscoisty ---*/
  Eddy_Viscosity = eddy_visc;

  return nonPhys;
}


//
//
/*!
 *\brief Set the vorticity
 */
//
//
bool CReactiveNSVariable::SetVorticity(bool val_limiter) {
  Vorticity[0] = 0.0;
  Vorticity[1] = 0.0;
  Vorticity[2] = Gradient_Primitive[VX_INDEX_GRAD + 1][0] - Gradient_Primitive[VX_INDEX_GRAD][1];

  if (nDim == 3) {
    Vorticity[0] = Gradient_Primitive[VX_INDEX_GRAD + 2][1] - Gradient_Primitive[VX_INDEX_GRAD + 1][2];
    Vorticity[1] = Gradient_Primitive[VX_INDEX_GRAD][2] - Gradient_Primitive[VX_INDEX_GRAD + 2][0];
  }
  return false;
}

//
//
/*!
 *\brief Set strain rate tensor magnitude
 */
//
//
bool CReactiveNSVariable::SetStrainMag(bool val_limiter) {
  su2double div_vel;
  unsigned short iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive, nDim + 1, nDim);

  div_vel = 0.0;
  for(iDim = 0; iDim < nDim; ++iDim)
    div_vel += Gradient_Primitive[VX_INDEX_GRAD + iDim][iDim];

  StrainMag = 0.0;

  /*--- Add diagonal part ---*/
  for(iDim = 0; iDim < nDim; ++iDim)
    StrainMag += (Gradient_Primitive[VX_INDEX_GRAD + iDim][iDim] - 1.0/3.0*div_vel)*
                 (Gradient_Primitive[VX_INDEX_GRAD + iDim][iDim] - 1.0/3.0*div_vel);

  /*--- Add off diagonals ---*/
  StrainMag += 0.5*(Gradient_Primitive[VX_INDEX_GRAD][1] + Gradient_Primitive[VX_INDEX_GRAD + 1][0])*
                   (Gradient_Primitive[VX_INDEX_GRAD][1] + Gradient_Primitive[VX_INDEX_GRAD + 1][0]);

  if (nDim == 3) {
    StrainMag += 0.5*(Gradient_Primitive[VX_INDEX_GRAD][2] + Gradient_Primitive[VX_INDEX_GRAD + 2][0])*
                     (Gradient_Primitive[VX_INDEX_GRAD][2] + Gradient_Primitive[VX_INDEX_GRAD + 2][0]);
    StrainMag += 0.5*(Gradient_Primitive[VX_INDEX_GRAD + 1][2] + Gradient_Primitive[VX_INDEX_GRAD + 2][1])*
                     (Gradient_Primitive[VX_INDEX_GRAD + 1][2] + Gradient_Primitive[VX_INDEX_GRAD + 2][1]);
  }

  StrainMag = std::sqrt(2.0*StrainMag);

  AD::SetPreaccOut(StrainMag);
  AD::EndPreacc();

  return false;
}
