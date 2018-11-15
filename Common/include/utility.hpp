#ifndef SU2_UTILITY_HPP
#define SU2_UTILITY_HPP

#include <set>
#include <map>
#include "datatype/vectorT.hpp"
#include "su2_assert.hpp"

namespace Utility {

  typedef std::map<std::string,unsigned short> MyMap;
  typedef std::set<std::string> MySet;
  using RealVec = Common::RealVec;

  /*!
    * \brief Read a single chemical reaction
    * \param[in] line - line of the reaction
    * \param[in] map_names - map with all species name
    * \param[out] stoich - vector of stechiometric coefficients
    * \param[out] stoich_exp - vector of stechiometric coefficients at the exponent of concetration
    * \param[out] species_names - detected symbols of the species
  */
  void Parse_Terms(std::string& line, MyMap& map_names, RealVec& stoich, RealVec& stoich_exp, MySet& species_names);

  /*!
    * \brief Preprocess to generate data tables for transport properties
    * \param[in] f_name - name of the file with the coefficients
  */
  void PreprocessTransp(const std::string& f_name);

  /*!
    * \brief Preprocess to generate data tables for thermodynamic properties
    * \param[in] f_name - name of the file with the coefficients
  */
  void PreprocessThermo(const std::string& f_name);

  /*!
    * \brief Preprocess to generate data tables for transport/thermodynamic properties for all species
    * \param[in] f_name - name of the file with the name of sub-file to genrate tables
  */
  void Preprocess(const std::string& f_name);


  //RealVec As1Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec Bs1Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec Cs1Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec Ds1Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec As2Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec Bs2Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec Cs2Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec Ds2Mu; /*!< \brief Coefficient for viscosity computation. */

  //RealVec As1Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec Bs1Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec Cs1Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec Ds1Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec As2Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec Bs2Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec Cs2Kappa; /*!< \brief Coefficient for conductivity computation. */

  //RealVec Ds2Kappa; /*!< \brief Coefficient for conductivity computation. */

}

#endif
