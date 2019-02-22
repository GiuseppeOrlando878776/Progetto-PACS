#ifndef SU2_SOLVER_REACTIVE
#define SU2_SOLVER_REACTIVE

#include "solver_structure.hpp"
#include "variable_reactive.hpp"

/*! \class CReactiveEulerSolver
 *  \brief Main class for defining a solver for chemically reacting inviscid flows.
 *  \author G. Orlando.
 */
class CReactiveEulerSolver: public CSolver {
public:
  using RealVec = CReactiveEulerVariable::RealVec;
  using RealMatrix = CReactiveNSVariable::RealMatrix;
  using LibraryPtr = CReactiveEulerVariable::LibraryPtr;

protected:
  static LibraryPtr library; /*!< \brief Smart pointer to the library that computes physical-chemical properties. */

  unsigned short nSpecies; /*!< \brief Total number of species. */
  unsigned short nPrimVarLim; /*!< \brief Number of primitive variables to limit. */

  bool  space_centered,  /*!< \brief True if space centered scheme used. */
        implicit,      /*!< \brief True if euler implicit scheme used. */
        grid_movement, /*!< \brief True if grid movement is used. */
        least_squares,  /*!< \brief True if computing gradients by least squares. */
        second_order, /*!< \brief True if second order recosntruction is applied. */
        limiter; /*!< \brief True if limiting strategy is applied. */

  bool  turbulent,   /*!< \brief True if a turbulence model is used. */
        tkeNeeded;   /*!< \brief True if a turbulent kinetic energy model is used. */

  RealVec   Lower_Limit,   /*!< \brief Lower limit conserved variables. */
            Upper_Limit;   /*!< \brief Upper limit conserved variables. */

  su2double Density_Inf,       /*!< \brief Free stream density. */
            Pressure_Inf,		  /*!< \brief Free stream pressure. */
	          Temperature_Inf;   /*!< \brief Trans.-rot. free stream temperature. */

  RealVec   Velocity_Inf,  /*!< \brief Free stream flow velocity. */
            MassFrac_Inf;  /*!< \brief Free stream species mass fraction. */

  RealVec   PrimVar_i,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at point i. */
            PrimVar_j,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at point j. */
            PrimVar_Vertex,  /*!< \brief Auxiliary nPrimVarGrad vector for storing primitive at boundary node. */
            PrimVar_Average,  /*!< \brief Auxiliary nPrimVarGrad vector for storing average of primitive.. */
            Partial_Res;  /*!< \brief Auxiliary nPrimVarGrad vector. */

  RealVec   Prim_i, /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at point i. */
            Prim_j, /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at point j. */
            Primitive; /*!< \brief Auxiliary nPrimVarLim vector for storing primitive at boundary node. */

  RealMatrix C_Mat,S_Mat; /*!< \brief Auxiliary matrices for least squares computation. */

  RealVec Ys_i,Ys_j;  /*!< \brief Auxiliary vectors to store mass fractions at node i and j. */
  RealVec Ys;         /*!< \brief Auxiliary vector to store mass fractions. */

public:

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveEulerSolver();

	/*!
	 * \overloaded constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
	 */
	CReactiveEulerSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveEulerSolver() {}

  /*!
   * \brief Set the simulation to explicit
   */
  inline void SetExplicit(void) {
    implicit = false;
  }

  /*!
 	 * \brief Looking for non physical points in the initial solution
 	 * \param[in] config - Definition of the particular problem.
 	 */
 	void Check_FreeStream_Solution(CConfig* config);

  /*!
	 * \brief Reading files in case of restart
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] val_filename - Name of the file for the restart
	 */
  void Load_Restart(CGeometry* geometry, CConfig* config, std::string val_filename);

  /*!
   * \brief Set primitive variables in each point reporting non physical data
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) override;

  /*!
   * \brief Set gradient primitive variables static const unsigned Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set gradient primitive variables static const unsigned weighted least squares.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Compute the limiter of the primitive variables.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetPrimitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set the fluid solver nondimensionalization.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Call MPI to set solution in case of parallel simulation
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Solution(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Call MPI to set limiter of primitive variables in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Primitive_Limiter(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Call MPI to set gradient of primitive variables in case of parallel simulation.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void Set_MPI_Primitive_Gradient(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
  void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh,
                     unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Compute the time step for solving the Euler equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
  void SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh, unsigned long Iteration) override;

  /*!
   * \brief Compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
   void SetResidual_DualTime(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                             unsigned short iRKStep, unsigned short iMesh, unsigned short RunTime_EqSystem) override;

  /*!
   * \brief Compute the spatial integration static const unsigned a centered scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Centered_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) override;

  /*!
   * \brief Compute the spatial integration static const unsigned a upwind scheme.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   */
  void Upwind_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config, unsigned short iMesh) override;

   /*!
    * \brief Source term integration.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] second_numerics - Description of the second numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iMesh - Index of the mesh in multigrid computations.
    */
   void Source_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CNumerics* second_numerics,
                        CConfig* config, unsigned short iMesh) override;

   /*!
    * \brief Set the free-stream solution all over the domain.
    * \param[in] config - Definition of the particular problem.
    */
    void SetFreeStream_Solution(CConfig* config) override;

   /*!
    * \brief Impose via the residual the Euler wall boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Euler_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics, CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose the far-field boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Far_Field(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                     CConfig* config, unsigned short val_marker) override {}

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
    * \brief Impose a supersonic inlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Supersonic_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                            CConfig *config, unsigned short val_marker) override;

   /*!
    * \brief Impose a supersonic outlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method.
    * \param[in] visc_numerics - Description of the numerical method.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Supersonic_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                             CConfig *config, unsigned short val_marker) override {}

   /*!
    * \brief Impose the outlet boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Outlet(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics, CNumerics* visc_numerics,
                  CConfig* config, unsigned short val_marker) override;

   /*!
    * \brief Impose the symmetry plane boundary condition.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] conv_numerics - Description of the numerical method for convective term.
    * \param[in] visc_numerics - Description of the numerical method for viscous term.
    * \param[in] config - Definition of the particular problem.
    * \param[in] val_marker - Surface marker where the boundary condition is applied.
    */
   void BC_Sym_Plane(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                     CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override {}

    /*!
     * \brief Set the initial conditions.
     * \param[in] geometry - Geometrical definition of the problem.
     * \param[in] solver_container - Container with all the solutions.
     * \param[in] config - Definition of the particular problem.
     * \param[in] ExtIter - External iteration.
    */
   void SetInitialCondition(CGeometry** geometry, CSolver*** solver_container, CConfig *config, unsigned long ExtIter) override;


   /*!
    * \brief Update the solution static const unsigned a Runge-Kutta scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
    */
   void ExplicitRK_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iRKStep) override;

   /*!
    * \brief Update the solution static const unsigned the explicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
   void ExplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

   /*!
    * \brief Update the solution static const unsigned an implicit Euler scheme.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] solver_container - Container vector with all the solutions.
    * \param[in] config - Definition of the particular problem.
    */
   void ImplicitEuler_Iteration(CGeometry* geometry, CSolver** solver_container, CConfig* config) override;

   /*!
    * \brief Get turbulent Prandtl number.
    */
   virtual su2double GetPrandtl_Turb() {}
};

/*! \class CReactiveNSSolver
 *  \brief Main class for defining a solver for chemically reacting viscous flows.
 *  \author G. Orlando.
 */
class CReactiveNSSolver:public CReactiveEulerSolver {
protected:

  su2double Viscosity_Inf;	/*!< \brief Viscosity at the infinity. */

  su2double Tke_Inf;	/*!< \brief Turbulent kinetic energy at infinity. */

  su2double Prandtl_Lam,    /*!< \brief Laminar Prandtl number. */
            Prandtl_Turb;   /*!< \brief Turbulent Prandtl number. */

public:

  /*!
	 * \brief Default constructor of the class.
	 */
  CReactiveNSSolver(): CReactiveEulerSolver(), Viscosity_Inf() {}

	/*!
	 * \overloaded Constructor of the class
	 * \param[in] geometry - Geometrical definition of the problem.
	 * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
	 */
	CReactiveNSSolver(CGeometry* geometry, CConfig* config, unsigned short iMesh);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CReactiveNSSolver() {}

  /*!
   * \brief Set gradient primitive variables static const unsigned Green Gauss.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_GG(CGeometry* geometry, CConfig* config) override;

  /*!
   * \brief Set gradient primitive variables static const unsigned weighted least squares.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
   void SetPrimitive_Gradient_LS(CGeometry* geometry, CConfig* config) override;

   /*!
    * \brief Set the fluid solver nondimensionalization.
    * \param[in] geometry - Geometrical definition of the problem.
    * \param[in] config - Definition of the particular problem.
    */
   void SetNondimensionalization(CGeometry* geometry, CConfig* config, unsigned short iMesh) override;

  /*!
   * \brief Preprocessing.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   * \param[in] RunTime_EqSystem - System of equations which is going to be solved.
   * \param[in] Output - boolean to determine whether to print output.
   */
   void Preprocessing(CGeometry* geometry, CSolver** solver_container, CConfig* config, unsigned short iMesh, unsigned short iRKStep,
                      unsigned short RunTime_EqSystem, bool Output) override;

  /*!
   * \brief Set primitive variables in each point reporting non physical data
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] Output - boolean to determine whether to print output.
   * \return - The number of non-physical points.
   */
  unsigned long SetPrimitive_Variables(CSolver** solver_container, CConfig* config, bool Output) override;

  /*!
   * \brief Compute the time step for solving the Navier-Stokes equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] Iteration - Index of the current iteration.
   */
   void SetTime_Step(CGeometry* geometry, CSolver** solver_container, CConfig* config,
                     unsigned short iMesh, unsigned long Iteration) override;

  /*!
   * \brief Compute the viscous residuals.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
   void Viscous_Residual(CGeometry* geometry, CSolver** solver_container, CNumerics* numerics,
                         CConfig* config, unsigned short iMesh, unsigned short iRKStep) override;

  /*!
   * \brief Impose the Navier-Stokes boundary condition (strong).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method for convective term.
   * \param[in] visc_numerics - Description of the numerical method for viscous term.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
  */
  void BC_HeatFlux_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                        CNumerics* visc_numerics, CConfig* config, unsigned short val_marker) override;

  /*!
   * \brief Impose an isothermal wall boundary condition (no-slip).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   */
  void BC_Isothermal_Wall(CGeometry* geometry, CSolver** solver_container, CNumerics* conv_numerics,
                          CNumerics* visc_numerics, CConfig* config,unsigned short val_marker) override;

  /*!
   * \brief Get turbulent kinetic energy at infinity.
   */
  inline su2double GetTke_Inf(void) override {
    return Tke_Inf;
  }

  /*!
   * \brief Get turbulent Prandtl number.
   */
  inline su2double GetPrandtl_Turb(void) override {
    return Prandtl_Turb;
  }

};

#endif
