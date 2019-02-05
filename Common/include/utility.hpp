#ifndef SU2_UTILITY_HPP
#define SU2_UTILITY_HPP

#include <set>
#include <map>
#include "su2_assert.hpp"
#include "../../externals/Eigen/Dense"

namespace Utility {

  typedef std::map<std::string,unsigned short> MyMap;
  typedef std::set<std::string> MySet;
  using RealMatrix = Eigen::MatrixXd;

  /*!
    * \brief Read a single chemical reaction
    * \param[in] line - line of the reaction
    * \param[in] n_reac - current number of reaction
    * \param[in] is_elem - flag for elementary reaction
    * \param[in] map_names - map with all species name
    * \param[out] stoich - vector of stechiometric coefficients
    * \param[out] stoich_exp - vector of stechiometric coefficients at the exponent of concetration
    * \param[out] species_names - detected symbols of the species
  */
  void Parse_Terms(std::string& line, unsigned n_reac, bool is_elem, const MyMap& map_names, RealMatrix& stoich,
                   RealMatrix& stoich_exp, MySet& species_names);

  /*!
    * \brief Read a single chemical reaction
    * \param[in] line - line of the reaction
    * \param[in] n_reac - current number of reaction
    * \param[in] is_elem - flag for elementary reaction
    * \param[in] map_names - map with all species name
    * \param[out] stoich - vector of stechiometric coefficients
    * \param[out] stoich_exp - vector of stechiometric coefficients at the exponent of concetration
    * \param[out] species_names - detected symbols of the species
  */
  void CompleteBackwardRate(std::string& line, unsigned n_reac, bool is_elem, const MyMap& map_names, const RealMatrix& stoich,
                            RealMatrix& stoich_exp);

}

#endif
