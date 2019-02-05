#include "../include/utility.hpp"

#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace Utility {
  //
  //
  /*--- This function reads a reaction from line ---*/
  void Parse_Terms(std::string& line, unsigned n_reac, bool is_elem, const MyMap& map_names, RealMatrix& stoich_coeff,
                   RealMatrix& stoich_coeff_exp, MySet& species_names) {
     auto size = line.size();
     std::size_t idx = 0;
     std::string coeff,symbol,exp_coeff,sub;

     /*--- Find first stoichiometric coefficient or first letter of species name ---*/
     while(!std::isdigit(line[idx]) and !std::isalpha(line[idx]))
       idx++;

     /*--- Extract leading coefficient ---*/
     while(std::isdigit(line[idx]) || std::ispunct(line[idx])) {
       coeff += line[idx];
       idx++;
     }

     do {
       /*--- Extract chemical symbol ---*/
       while(std::isalpha(line[idx]) || std::isdigit(line[idx])) {
         symbol += line[idx];
         idx++;
       }
     } while(!std::ispunct(line[idx]) and !std::isspace(line[idx]) and line[idx] != '+' and idx < size);

     auto it = map_names.find(symbol);
     SU2_Assert(it!=map_names.end(),std::string("The symbol " + symbol + " is not in the mixture list"));
     species_names.insert(symbol);

     /*--- Retrieve the index of the species detected ---*/
     unsigned short idx_species = it->second;

     /*--- Saving stoichiometric coefficient ---*/
     double coefficient;
     if(coeff.empty())
      coefficient = 1.0;
     else {
       std::istringstream str_to_coeff(coeff);
       str_to_coeff>>coefficient;
     }
     stoich_coeff(idx_species,n_reac - 1) = coefficient;

     /*--- Saving possibly coefficient at the exponent of coencentration in reaction rate ---*/
     if(std::ispunct(line[idx])) {
       idx++;
       while(std::isdigit(line[idx]) || std::ispunct(line[idx])) {
         exp_coeff += line[idx];
         idx++;
       }
     }
     double exp_coefficient = 0.0;
     if(!exp_coeff.empty()) {
       std::istringstream str_to_expcoeff(exp_coeff);
       str_to_expcoeff>>exp_coefficient;
       stoich_coeff_exp(n_reac - 1,idx_species) = exp_coefficient;
     }
     /*--- If the reaction is elementary we need absolutely the exponent so if it is not present we put it ---*/
     if(is_elem && exp_coefficient == 0.0)
       stoich_coeff_exp(n_reac - 1,idx_species) = stoich_coeff(idx_species,n_reac - 1);

     /*--- Any more terms to extract? ---*/
     if(idx == size)
      return;

     sub = line.substr((idx + 1), (size - idx));
     if(sub.empty())
      return;
     else
      Parse_Terms(sub,map_names,stoich_coeff,stoich_coeff_exp,species_names);
  }

  //
  //
  /*--- This function sets the rate exponent for reactants side ---*/
  void CompleteBackwardRate(std::string& line, unsigned n_reac, bool is_elem, const MyMap& map_names, const RealMatrix& stoich_coeff,
                            RealMatrix& stoich_coeff_exp) {
    if(!is_elem)
      return;

    auto size = line.size();
    std::size_t idx = 0;
    std::string coeff,symbol,exp_coeff,sub;

    /*--- Find first stoichiometric coefficient or first letter of species name ---*/
    while(!std::isdigit(line[idx]) and !std::isalpha(line[idx]))
      idx++;

    /*--- Extract leading coefficient ---*/
    while(std::isdigit(line[idx]) || std::ispunct(line[idx])) {
      coeff += line[idx];
      idx++;
    }

    do {
      /*--- Extract chemical symbol ---*/
      while(std::isalpha(line[idx]) || std::isdigit(line[idx])) {
        symbol += line[idx];
        idx++;
      }
    } while(!std::ispunct(line[idx]) and !std::isspace(line[idx]) and line[idx] != '+' and idx < size);

    auto it = map_names.find(symbol);

    /*--- Retrieve the index of the species detected ---*/
    unsigned short idx_species = it->second;
    if(stoich_coeff_exp(n_reac - 1,idx_species) != 0.0)
      return;

    stoich_coeff_exp(n_reac - 1,idx_species) = 1.0 - stoich_coeff(idx_species,n_reac - 1);
  }

} /*--- End of namespace Utility ---*/
