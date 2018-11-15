#include "../include/utility.hpp"

#include <fstream>
#include <sstream>
#include <cmath>
#include <functional>
#include <exception>
#include <algorithm>
#include <iomanip>

namespace Utility {
  //
  //
  /*--- This function reads a reaction from line ---*/
  void Parse_Terms(std::string& line,MyMap& map_names,RealVec& stoich_coeff,RealVec& stoich_coeff_exp,MySet& species_names) {
     auto size = line.size();
     std::size_t idx = 0;
     std::string coeff,symbol,exp_coeff,sub;

     while(!std::isdigit(line[idx]) and !std::isalpha(line[idx]))
        idx++;

     //Extract leading coefficient
     while(std::isdigit(line[idx]) || std::ispunct(line[idx])) {
          coeff += line[idx];
          idx++;
     }

     do {
          //Extract chemical symbol
          while(std::isalpha(line[idx]) || std::isdigit(line[idx])) {
               symbol += line[idx];
               idx++;
          }
     } while(!std::ispunct(line[idx]) and !std::isspace(line[idx]) and line[idx] != '+' and idx < size);

     SU2_Assert(map_names.find(symbol)!=map_names.end(),std::string("The symbol " + symbol + " is not in the mixture list"));
     species_names.insert(symbol);

     //Saving stoichiometric coefficient
     double coefficient;
     if(coeff.empty())
      coefficient = 1.0;
     else {
       std::istringstream str_to_coeff(coeff);
       str_to_coeff>>coefficient;
     }
     stoich_coeff.push_back(coefficient);

     //Saving possibly coefficient at the exponent of coencentration in reaction rate
     if(std::ispunct(line[idx])) {
       idx++;
       while(std::isdigit(line[idx]) || std::ispunct(line[idx])) {
            exp_coeff += line[idx];
            idx++;
       }
     }
     double exp_coefficient;
     if(!exp_coeff.empty()) {
       std::istringstream str_to_expcoeff(exp_coeff);
       str_to_expcoeff>>exp_coefficient;
       stoich_coeff_exp.push_back(exp_coefficient);
     }

     //Any more terms to extract?
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
  /*--- Preprocess transport properties ---*/
  void PreprocessTransp(const std::string& f_name) {
    std::string line;
    unsigned n_line = 0;

    std::vector<std::string> read_strings;

    std::ifstream transpfile(f_name);
    if(transpfile.is_open()) {
      while(transpfile.good() and !transpfile.eof()) {
        std::getline(transpfile,line);
        // We clearly avoid reading and empty lines
        if(!line.empty() and !std::ispunct(line[0])) {
          if(n_line == 0)
            SU2_Assert(std::isalpha(line[0]),"Species not specified in reading data");
          else
            SU2_Assert(std::isdigit(line[0]),"Unknown parameter in the file");
          read_strings.push_back(line);
        }
      }
    transpfile.close();
    }
    else {
      std::cerr<<"Unable to open the file "<<f_name<<std::endl;
      std::exit(1);
    }

    //Now we read coefficients
    std::string species = read_strings[0];
    if(species == "C4H6")
      SU2_Assert(read_strings.size() == 5,std::string("Wrong number of lines in the file " + f_name));
    else
      SU2_Assert(read_strings.size() == 8,std::string("Wrong number of lines in the file " + f_name));

    double step,temp_min,temp_max;
    double A1Mu,B1Mu,C1Mu,D1Mu,A2Mu,B2Mu,C2Mu,D2Mu;
    double A1Kappa,B1Kappa,C1Kappa,D1Kappa,A2Kappa,B2Kappa,C2Kappa,D2Kappa;
    unsigned short Butidene_ID = 0;
    double Butidene_MM = 0.0;

    if(species != "C4H6") {
      for(auto idx = 1; idx < 9; ++idx) {
        std::istringstream curr_line(read_strings[idx]);
        if(idx == 1)
          curr_line>>step;
        if(idx == 2)
          curr_line>>temp_min;
        if(idx == 3)
          curr_line>>temp_max;
        if(idx == 4) {
          curr_line>>A1Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing A1Mu coefficient for " + species));
          curr_line>>B1Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing B1Mu coefficient for " + species));
          curr_line>>C1Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing C1Mu coefficient for " + species));
          curr_line>>D1Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing D1Mu coefficient for " + species));
        }
        if(idx == 5) {
          curr_line>>A2Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing A2Mu coefficient for " + species));
          curr_line>>B2Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing B2Mu coefficient for " + species));
          curr_line>>C2Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing C2Mu coefficient for " + species));
          curr_line>>D2Mu;
          SU2_Assert(!curr_line.fail(),std::string("Missing D2Mu coefficient for " + species));
        }
        if(idx == 6) {
          curr_line>>A1Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing A1Mu coefficient for " + species));
          curr_line>>B1Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing B1Mu coefficient for " + species));
          curr_line>>C1Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing C1Mu coefficient for " + species));
          curr_line>>D1Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing D1Mu coefficient for " + species));
        }
        if(idx == 7) {
          curr_line>>A2Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing A2Mu coefficient for " + species));
          curr_line>>B2Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing B2Mu coefficient for " + species));
          curr_line>>C2Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing C2Mu coefficient for " + species));
          curr_line>>D2Kappa;
          SU2_Assert(!curr_line.fail(),std::string("Missing D2Mu coefficient for " + species));
        }
      }
    }
    else {
      for(auto idx = 1; idx < 9; ++idx) {
        std::istringstream curr_line(read_strings[idx]);
        if(idx == 1) {
          curr_line>>Butidene_ID;
          SU2_Assert(std::isdigit(Butidene_ID),"Wrong format of index for butidene");
        }
        if(idx == 2)
          curr_line>>Butidene_MM;
        if(idx == 3)
          curr_line>>step;
        if(idx == 4)
          curr_line>>temp_min;
        if(idx == 5)
          curr_line>>temp_max;
      }
    }
    SU2_Assert(temp_max > temp_min,"The range of temperature is swapped");
    if(temp_min < 200 or temp_max > 6000)
      throw std::out_of_range("Temperature is out of available range");
    //Now we generate data
    std::string out_name = species + "_transp.txt";
    std::vector<double> temp_data;
    for(double curr_temp = temp_min; curr_temp <= temp_max; curr_temp = std::min(curr_temp + step,temp_max))
      temp_data.push_back(curr_temp);
    //First we focus on viscosity
    std::function<double(double)> compute_mu,compute_kappa;
    std::vector<double> mu_data;
    std::vector<double> kappa_data;
    if(species == "C4H6") {
      const double a = 1.16145;
      const double b = 0.14874;
      const double c = 0.52487;
      const double d = 0.77320;
      const double e = 2.16178;
      const double f = 2.43787;

      const double AF = 0.192;
      const double Vc = 220;
      const double Fc = 1.0 - 0.2756*AF;
      auto compute_mu_tmp = [=](const double temp) {
        const double Tstar = 0.003*temp;
        const double Omega_v = a*std::pow(Tstar, -b) + c*std::exp(-d*Tstar) + e*std::exp(-f*Tstar);
        return 40.785*Fc*std::sqrt(Butidene_MM*temp)/(std::pow(Vc,2.0/3)*Omega_v);
      };
      compute_mu = std::move(compute_mu_tmp);
    }
    else {
      auto compute_mu_tmp = [=](const double temp) {
        if(temp >= 200 and temp <= 1000)
          return std::exp(A1Mu*std::log(temp) + B1Mu/temp + C1Mu/(temp*temp) + D1Mu);
        else
          return std::exp(A2Mu*std::log(temp) + B2Mu/temp + C2Mu/(temp*temp) + D2Mu);
      };
      compute_mu = std::move(compute_mu_tmp);
    }
    std::transform(temp_data.cbegin(),temp_data.cend(),std::back_inserter(mu_data),compute_mu);
    //Now we focus on conductivity
    std::function<double(double)> compute_mu,compute_kappa;
    std::vector<double> mu_data;
    std::vector<double> kappa_data;
    if(species == "C4H6") {
      const double a = 1.16145;
      const double b = 0.14874;
      const double c = 0.52487;
      const double d = 0.77320;
      const double e = 2.16178;
      const double f = 2.43787;

      const double AF = 0.192;
      const double Vc = 220;
      const double Fc = 1.0 - 0.2756*AF;
      auto compute_mu_tmp = [=](const double temp) {
        const double Tstar = 0.003*temp;
        const double Omega_v = a*std::pow(Tstar, -b) + c*std::exp(-d*Tstar) + e*std::exp(-f*Tstar);
        return 40.785*Fc*std::sqrt(Butidene_MM*temp)/(std::pow(Vc,2.0/3)*Omega_v);
      };
      compute_mu = std::move(compute_mu_tmp);
    }
    else {
      auto compute_mu_tmp = [=](const double temp) {
        if(temp >= 200 and temp <= 1000)
          return std::exp(A1Mu*std::log(temp) + B1Mu/temp + C1Mu/(temp*temp) + D1Mu);
        else
          return std::exp(A2Mu*std::log(temp) + B2Mu/temp + C2Mu/(temp*temp) + D2Mu);
      };
      compute_mu = std::move(compute_mu_tmp);
    }
    std::transform(temp_data.cbegin(),temp_data.cend(),std::back_inserter(mu_data),compute_mu);
    //Save data
    std::ofstream output(out_name);
    for(auto idx = 0; idx < temp_data.size(); ++idx)
      output<<temp_data[idx]<<std::setw(10)<<mu_data[idx]<<std::setw(10)<<kappa_data[idx]<<std::endl;
  }

  //
  //
  /*--- Preprocess thermodynamic properties ---*/
  void PreprocessThermo(const std::string& f_name) {

  }

  //
  //
  /*--- Preprocess all properties ---*/
  void Preprocess(const std::string& f_name) {

  }



//
//
/*--- Computing molecular viscosity of each species---*/
/*
void ComputeViscosities(const su2double temp) {
  unsigned short iSpecies;
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(temp > 200 and temp <= 1000)
      Viscosities[iSpecies] = std::exp(As1Mu[iSpecies]*std::log(temp) + Bs1Mu[iSpecies]/temp + Cs1Mu[iSpecies]/(temp*temp) + Ds1Mu[iSpecies]);
    else if(temp > 1000 and temp <= 6000)
      Viscosities[iSpecies] = std::exp(As2Mu[iSpecies]*std::log(temp) + Bs2Mu[iSpecies]/temp + Cs2Mu[iSpecies]/(temp*temp) + Ds2Mu[iSpecies]);
  }

  auto is_butadiene = Species_Names.find("C4H6");
  if(is_butadiene != Species_Names.end()) {
    const su2double a = 1.16145;
    const su2double b = 0.14874;
    const su2double c = 0.52487;
    const su2double d = 0.77320;
    const su2double e = 2.16178;
    const su2double f = 2.43787;

    const su2double AF = 0.192;
    const su2double Vc = 220;
    const su2double Tstar = 0.003*temp;
    const su2double Omega_v = a*std::pow(Tstar, -b) + c*std::exp(-d*Tstar) + e*std::exp(-f*Tstar);
    const su2double Fc = 1.0 - 0.2756*AF;

    const unsigned idx_butadiene = is_butadiene->second - 1;
    Viscosities[idx_butadiene] = 40.785*Fc*std::sqrt(mMasses[idx_butadiene]*temp)/(std::pow(Vc,2.0/3)*Omega_v);
  }
}
*/
//
//
/*--- Computing thermal conductivity of each species ---*/
/*
void ReactingModelLibrary::ComputeConductivities(const su2double temp) {
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(temp > 200 and temp <= 1000)
      Thermal_Conductivities[iSpecies] = std::exp(As1Kappa[iSpecies]*std::log(temp) + Bs1Kappa[iSpecies]/temp +
                                                  Cs1Kappa[iSpecies]/(temp*temp) + Ds1Kappa[iSpecies]);
    else if(temp > 1000 and temp <= 6000)
      Thermal_Conductivities[iSpecies] = std::exp(As2Kappa[iSpecies]*std::log(temp) + Bs2Kappa[iSpecies]/temp +
                                                  Cs2Kappa[iSpecies]/(temp*temp) + Ds2Kappa[iSpecies]);
  }

  auto is_butadiene = Species_Names.find("C4H6");
  if(is_butadiene != Species_Names.end()) {
    const su2double a = 1.16145;
    const su2double b = 0.14874;
    const su2double c = 0.52487;
    const su2double d = 0.77320;
    const su2double e = 2.16178;
    const su2double f = 2.43787;

    const su2double AF = 0.192;
    const su2double Vc = 220;
    const su2double Tc = 425.17;
    const su2double Tstar = 0.003*temp;
    const su2double Omega_v = a*std::pow(Tstar, -b) + c*std::exp(-d*Tstar) + e*std::exp(-f*Tstar);
    const su2double Fc = 1.0 - 0.2756*AF;

    const unsigned idx_butadiene = is_butadiene->second - 1;
    const su2double mu_butadiene = 40.785*Fc*std::sqrt(mMasses[idx_butadiene]*temp)/(std::pow(Vc,2.0/3)*Omega_v);
    const su2double CP = ComputeCP(temp);
    const su2double CV = CP - Rgas;
    const su2double CV_ov_MM = mu_butadiene*CV/mMasses[idx_butadiene];
    switch(Butidene_ID) {
      case 0:
        Thermal_Conductivities[idx_butadiene] = 1.32*CV_ov_MM + 1.77*Rgas/CV;
        break;
      case 1:
        Thermal_Conductivities[idx_butadiene] = 1.15*CV_ov_MM + 2.03*Rgas/CV;
        break;
      case 2: {
        const su2double alpha = CV/Rgas - 1.5;
        const su2double beta = 0.7682 - 0.7109*AF + 1.3168*AF*AF;
        const su2double Tbar = temp/Tc;
        const su2double zeta = 2.0 + 10.5*Tbar*Tbar;
        const su2double psi = 1.0 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*zeta) / (0.6366 + beta*zeta + 1.061*alpha*beta));
        Thermal_Conductivities[idx_butadiene] = CV_ov_MM + 3.77*Rgas*psi/CV;
      }
        break;
      default:
        throw std::out_of_range("Unknown index for computing Butidene conductivity");

    }

  }
}
*/

}
