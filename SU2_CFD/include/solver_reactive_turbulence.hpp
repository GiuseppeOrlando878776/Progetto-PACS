#ifndef SU2_SOLVER_REACTIVE_TURBULENCE
#define SU2_SOLVER_REACTIVE_TURBULENCE

#include "solver_structure.hpp"

#include "variable_reactive_turbulence.hpp"

/*!
 * \class CReactiveTurbSolver
 * \brief Main class for defining the k-\omega turbulence model solver.
 * \author G. Orlando
 */
class CReactiveTurbSolver: public CTurbSolver {
protected:
  unsigned short nSpecies;    /*!< \brief Number of species. */

  su2double kine_Inf,         /*!< \brief Free-stream turbulent kinetic energy. */
            omega_Inf;        /*!< \brief Free-stream specific dissipation. */

public:
  /*!
   * \brief Constructor of the class.
   */
  CReactiveTurbSolver(): CTurbSolver(), nSpecies(), kine_Inf(), omega_Inf() {}

  /*!
   * \overload Constructor of the class
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  CReactiveTurbSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CReactiveTurbSolver() {}

  /*!
   * \brief Load file for restart
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] rank - Rank of current process.
   * \param[in] muT_Inf - Freestream Eddy viscosity.
   */
  void LoadRestart(CGeometry* geometry, CConfig* config, int rank, su2double muT_Inf);

  /*!
   * \brief Restart residual and compute gradients.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh, unsigned short iRKStep,
                     unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Computes the eddy viscosity.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Postprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Index of the current Runge-Kutta step.
   */
  void Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config,
                        unsigned short iMesh, unsigned short iRKStep) override;

  /*!
   * \brief Source term computation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] second_numerics - Description of the second numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                       CNumerics* second_numerics, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the Navier-Stokes wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                          CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the Far Field boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                    CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the inlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Inlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose the outlet boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                 CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Set the solution using the freestream values.
   * \param[in] config - Definition of the particular problem.
   */
  inline void SetFreeStream_Solution(CConfig* config) override {
    for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
      node[iPoint]->SetSolution(CReactiveTurbVariable::TURB_KINE_INDEX, kine_Inf);
      node[iPoint]->SetSolution(CReactiveTurbVariable::TURB_OMEGA_INDEX, omega_Inf);
    }
  }

};

#endif
