#ifndef SU2_REACTIVE_TURBULENT
#define SU2_REACTIVE_TURBULENT

#include "variable_structure.hpp"
//#include "../../Common/include/su2_assert.hpp"
//#include "../../Common/include/datatypes/vectorT.hpp"

/*!
 * \class CReactiveTurbVariable
 * \brief Main class for defining the variables of the k - \omega turbulence model.
 * \author G. Orlando.
 */

class CReactiveTurbVariable : public CTurbVariable {
public:
  /*!
   * \brief Default constructor of the class.
   */
  CReactiveTurbVariable(): CTurbVariable() {}

  /*!
   * \overload Class constructor
   * \param[in] val_rho_kine - Turbulent kinetic energy value (initialization value).
   * \param[in] val_rho_omega - Specific dissipation rate (initialization value).
   * \param[in] val_muT - Eddy viscosity (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CReactiveTurbVariable(su2double val_rho_kine, su2double val_rho_omega, su2double val_muT, unsigned short val_nDim,
                        unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CReactiveTurbVariable() = default;

};

#endif
