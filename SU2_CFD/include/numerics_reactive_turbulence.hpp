#ifndef SU2_NUMERICS_REACTIVE_TURBULENCE
#define SU2_NUMERICS_REACTIVE_TURBULENCE

#include "numerics_structure.hpp"

/*!
 * \class CUpwReactiveTurb
 * \brief Class for doing a scalar upwind solver for the k-\omega turbulence model equations.
 * \author G.Orlando
 */
class CUpwReactiveTurb: public CNumerics {
public:
  typedef std::vector<su2double> RealVec;

protected:
  bool implicit, /*!< \brief True if euler implicit scheme used. */
       grid_movement; /*!< \brief True if grid movement is used. */

  unsigned short nSpecies;  /*!< \brief Number of species. */

  RealVec Velocity_i, /*!< \brief Velocity at node i. */
          Velocity_j; /*!< \brief Velcoity at node j. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwReactiveTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwReactiveTurb() {}

  /*!
   * \brief Set number of species.
   * \param[in] val_ns - Number of species of the problem.
   */
  inline void SetnSpecies(unsigned short val_ns) {
    nSpecies = val_ns;
  }

  /*!
   * \brief Compute the scalar upwind flux between two nodes i and j.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) override;
};


/*!
 * \class CAvgGrad_ReactiveTurb
 * \brief Class for computing viscous term using average of gradient (k-\omega turbulence model).
 * \author G.Orlando
 */
class CAvgGrad_ReactiveTurb: public CNumerics {
public:
  typedef std::vector<su2double> RealVec;
  typedef su2double** SU2Matrix;

protected:
  bool implicit;                        /*!< \brief True if euler implicit scheme used. */
  unsigned short nSpecies;              /*!< \brief Number of species. */

  RealVec Edge_Vector;                  /*!< \brief Vector from node i to node j. */

  SU2Matrix Mean_GradTurbVar;              /*!< \brief Average of gradients at cell face */
  RealVec Proj_Mean_GradTurbVar_Normal;   /*!< \brief Mean_gradTurbVar DOT normal */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_ReactiveTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_ReactiveTurb();

  /*!
   * \brief Set number of species.
   * \param[in] val_ns - Number of species of the problem.
   */
  inline void SetnSpecies(unsigned short val_ns) {
    nSpecies = val_ns;
  }

  /*!
   * \brief Compute the viscous turbulent residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** Jacobian_i, su2double** Jacobian_j, CConfig* config) override;
};


/*!
 * \class CSourcePieceWise_ReactiveTurb
 * \brief Class for integrating the source terms of the k-\omega turbulence model equations.
 * \author G.Orlando
 */
class CSourcePieceWise_ReactiveTurb: public CNumerics {
protected:
  bool implicit;                        /*!< \brief True if euler implicit scheme used. */
  unsigned short nSpecies;              /*!< \brief Number of species. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_ReactiveTurb(unsigned short val_nDim, unsigned short val_nVar, CConfig* config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_ReactiveTurb() {}

  /*!
   * \brief Set number of species.
   * \param[in] val_ns - Number of species of the problem.
   */
  inline void SetnSpecies(unsigned short val_ns) {
    nSpecies = val_ns;
  }

  /*!
   * \brief Compute the residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double* val_residual, su2double** val_Jacobian_i, su2double** val_Jacobian_j, CConfig* config) override;
};

#endif
