#include "../include/reacting_model_library.hpp"
#include "../include/not_setup_exception.hpp"
#include "../include/spline.hpp"
#include "../include/utility.hpp"

#include <experimental/filesystem>
#include <fstream>
#include <sstream>
#include <algorithm>

namespace Framework {

  //
  //
  /*--- Setting gas constant for each species ---*/
  void ReactingModelLibrary::SetRiGas(void) {
    Ri.clear();
    Ri.reserve(nSpecies);
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Ri.push_back(R_ungas/mMasses[iSpecies]);
  }

  //
  //
  /*--- Setting gas constant for mixture ---*/
  void ReactingModelLibrary::SetRgas(const Vec& ys) {
    SU2_Assert(ys.size() == nSpecies,"The dimension of vector rs doesn't match nSpecies");
    Rgas = std::inner_product(ys.cbegin(),ys.cend(),Ri.cbegin(),0.0);
  }

  //
  //
  /*--- Computing gas constant for mixture ---*/
  double ReactingModelLibrary::ComputeRgas(const Vec& ys) {
    SetRgas(ys);
    return Rgas;
  }

  //
  //
  /*--- Setting molar fraction for each species ---*/
  void ReactingModelLibrary::SetMolarFractions(const Vec& xs) {
    SU2_Assert(xs.size() == nSpecies,"The dimension of vector xs doesn't match nSpecies");
    Xs = (RealVec) xs;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      if(Xs[iSpecies] < 0.0)
        Xs[iSpecies] = 1E-20;

      SU2_Assert(Xs[iSpecies]<=1.0, std::string("The molar fraction of species number " + std::to_string(iSpecies) + "is greater than 1"));
    }
  }

  //
  //
  /* This function sets the molar fractions from mass fractions */
  void ReactingModelLibrary::SetMolarFromMass(const Vec& ys) {
    SU2_Assert(ys.size() == nSpecies,"The dimension of vector ys doesn't match nSpecies");

    unsigned short iSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      SU2_Assert(ys[iSpecies] <= 1.0,std::string("The mass fraction of species number " + std::to_string(iSpecies) + "is greater than one"));
      Xs[iSpecies] = ys[iSpecies]/mMasses[iSpecies];
    }

    double massTot = std::accumulate(ys.cbegin(),ys.cend(),0.0)/std::accumulate(Xs.cbegin(),Xs.cend(),0.0);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Xs[iSpecies] *= massTot;
  }

  //
  //
  /* This function returns the molar fractions from mass fractions */
  ReactingModelLibrary::Vec ReactingModelLibrary::GetMolarFromMass(const Vec& ys) {
    SetMolarFromMass(ys);
    return Vec(Xs);
  }

  //
  //
  /*--- Setting mass fraction for each species ---*/
  void ReactingModelLibrary::SetMassFractions(const Vec& ys) {
    SU2_Assert(ys.size() == nSpecies,"The dimension of vector ys doesn't match nSpecies");
    Ys = (RealVec) ys;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      if(Ys[iSpecies] < 0.0)
        Ys[iSpecies] = 1E-20;

      SU2_Assert(Ys[iSpecies]<=1.0,std::string("The mass fraction of species number " + std::to_string(iSpecies) + "is greater than 1"));
    }
  }

  //
  //
  /* This function sets the mass fractions from molar fractions */
  void ReactingModelLibrary::SetMassFromMolar(const Vec& xs) {
    SU2_Assert(xs.size() == nSpecies,"The dimension of vector xs doesn't match nSpecies");

    unsigned short iSpecies;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      SU2_Assert(xs[iSpecies] <= 1.0,std::string("The molar fraction of species number " + std::to_string(iSpecies) + "is greater than 1"));
      Ys[iSpecies] = xs[iSpecies]*mMasses[iSpecies];
    }

    double massTot = std::accumulate(Ys.cbegin(),Ys.cend(),0.0)/std::accumulate(xs.cbegin(),xs.cend(),0.0);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Ys[iSpecies] /= massTot;
  }

  //
  //
  /* This function returns the mass fractions from molar fractions */
  ReactingModelLibrary::Vec ReactingModelLibrary::GetMassFromMolar(const Vec& xs) {
    SetMassFromMolar(xs);
    return Vec(Ys);
  }

  //
  //
  /* This function computes gamma and the frozen sound speed */
  void ReactingModelLibrary::Gamma_FrozenSoundSpeed(const double temp, const Vec& ys, double& gamma, double& sound_speed) {
    double Cp = ComputeCP(temp,ys);
    SetRgas(ys);
    double Cv = Cp - Rgas;
    gamma = Cp/Cv;
    sound_speed = std::sqrt(gamma*Rgas*temp);
  }

  //
  //
  /* This function computes the frozen gamma */
  double ReactingModelLibrary::ComputeFrozenGamma(const double temp, const Vec& ys) {
    double Cp = ComputeCP(temp,ys);
    SetRgas(ys);
    double Cv = Cp - Rgas;
    return Cp/Cv;
  }

  //
  //
  /* This function computes the frozen sound speed */
  inline double ReactingModelLibrary::ComputeFrozenSoundSpeed(const double temp, const Vec& ys) {
    double gamma = ComputeFrozenGamma(temp,ys);
    return std::sqrt(gamma*Rgas*temp);
  }

  //
  //
  /* This function computes the frozen sound speed */
  inline double ReactingModelLibrary::ComputeFrozenSoundSpeed(const double temp, const Vec& ys, const double press, const double rho) {
    double gamma = ComputeFrozenGamma(temp,ys);
    return std::sqrt(gamma*press/rho);
  }

  //
  //
  /* This function computes pressure at given temperature and density */
  inline double ReactingModelLibrary::ComputePressure(const double temp, const double rho, const Vec& ys) {
    SetRgas(ys);
    return rho*temp*Rgas;
  }

  //
  //
  /* This function computes density at given temperature and pressure */
  inline double ReactingModelLibrary::ComputeDensity(const double temp, const double pressure, const Vec& ys) {
    SetRgas(ys);
    return pressure/(temp*Rgas);
  }

  //
  //
  /* This function computes internal energy at given temperature and pressure */
  inline double ReactingModelLibrary::ComputeEnergy(const double temp, const Vec& ys) {
    double enthalpy = ComputeEnthalpy(temp,ys);
    SetRgas(ys);
    return enthalpy - temp*Rgas;
  }

  //
  //
  /* This function computes density,internal energy and static enthalpy at given temperature and pressure */
  inline void ReactingModelLibrary::Density_Enthalpy_Energy(const double temp, const double pressure, const Vec& ys, RealVec& dhe) {
    dhe.resize(3);
    dhe[0] = ComputeDensity(temp,pressure,ys);
    dhe[1] = ComputeEnthalpy(temp,ys);
    dhe[2] = dhe[1] - pressure/dhe[0];
  }

  //
  //
  /* This function computes the static enthalpy for each species */
  ReactingModelLibrary::Vec ReactingModelLibrary::ComputePartialEnthalpy(const double temp) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Enthalpies[iSpecies] = MathTools::GetSpline(*std::get<0>(Enth_Spline[iSpecies]),std::get<1>(Enth_Spline[iSpecies]),
                                                   std::get<2>(Enth_Spline[iSpecies]),temp)/mMasses[iSpecies];

    return Vec(Enthalpies);
  }

  //
  //
  /* This function computes the static enthalpy of the mixture  */
  double ReactingModelLibrary::ComputeEnthalpy(const double temp, const Vec& ys) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Enthalpies[iSpecies] = MathTools::GetSpline(*std::get<0>(Enth_Spline[iSpecies]),std::get<1>(Enth_Spline[iSpecies]),
                                                   std::get<2>(Enth_Spline[iSpecies]),temp)/mMasses[iSpecies];

    SetMassFractions(ys);
    return std::inner_product(Ys.cbegin(),Ys.cend(),Enthalpies.cbegin(),0.0);
  }


  //
  //
  /* This function computes the specific heat at constant pressure */
  double ReactingModelLibrary::ComputeCP(const double temp, const Vec& ys) {
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      CPs[iSpecies] = MathTools::GetSpline(*std::get<0>(Cp_Spline[iSpecies]),std::get<1>(Cp_Spline[iSpecies]),
                                            std::get<2>(Cp_Spline[iSpecies]),temp)/mMasses[iSpecies];

    SetMassFractions(ys);
    return std::inner_product(Ys.cbegin(),Ys.cend(),CPs.cbegin(),0.0);
  }

  //
  //
  /*--- Computing molecular viscosity of each species---*/
  void ReactingModelLibrary::ComputeViscosities(const double temp) {
    SU2_Assert(Mu_Spline.size() == nSpecies,"The dimension of Mu_Spline doesn't match nSpecies");
    SU2_Assert(Viscosities.size() == nSpecies,"The dimension of Viscosities doesn't match nSpecies");

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Viscosities[iSpecies] = MathTools::GetSpline(*std::get<0>(Mu_Spline[iSpecies]),std::get<1>(Mu_Spline[iSpecies]),
                                                    std::get<2>(Mu_Spline[iSpecies]),temp);
  }

  //
  //
  /*--- Computing viscosity of the mixture---*/
  double ReactingModelLibrary::ComputeEta(const double temp, const Vec& ys) {
    ComputeViscosities(temp);
    SU2_Assert(ys.size() == nSpecies,"The dimension of ys doesn't match nSpecies");

    unsigned short iSpecies, jSpecies;
    double phi;
    RealVec ys_over_mm;
    ys_over_mm.reserve(nSpecies);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      ys_over_mm.push_back(ys[iSpecies]/mMasses[iSpecies]);

    /*--- Mixture viscosity calculation as sum weighted over PHI ---*/
    double Viscosity_Mixture = 0.0;
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)  {
      phi = 0.0;
      for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
        phi += ys_over_mm[jSpecies]/std::sqrt(8.0*(1.0 + mMasses[iSpecies]/mMasses[jSpecies]))*
               (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25))*
               (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25));

      Viscosity_Mixture += Viscosities[iSpecies]*ys_over_mm[iSpecies]/phi;
    }
    return Viscosity_Mixture;
  }

  //
  //
  /*--- Computing thermal conductivity of each species ---*/
  void ReactingModelLibrary::ComputeConductivities(const double temp) {
    SU2_Assert(Kappa_Spline.size() == nSpecies,"The dimension of Kappa_Spline doesn't match nSpecies");
    SU2_Assert(Thermal_Conductivities.size() == nSpecies,"The dimension of Thermal_Conductivities doesn't match nSpecies");

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      Thermal_Conductivities[iSpecies] = MathTools::GetSpline(*std::get<0>(Kappa_Spline[iSpecies]),std::get<1>(Kappa_Spline[iSpecies]),
                                                               std::get<2>(Kappa_Spline[iSpecies]),temp);
  }

  //
  //
  /*--- Computing thermal conductivity of the mixtures---*/
  double ReactingModelLibrary::ComputeLambda(const double temp, const Vec& ys) {
    ComputeConductivities(temp);
    ComputeViscosities(temp);
    SU2_Assert(ys.size() == nSpecies,"The dimension of ys doesn't match nSpecies");

    double Thermal_Conductivity_Mixture = 0.0;
    unsigned short iSpecies,jSpecies;
    double phi;
    RealVec ys_over_mm;
    ys_over_mm.reserve(nSpecies);
    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
      ys_over_mm.push_back(ys[iSpecies]/mMasses[iSpecies]);


    for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      phi = 0.0;
	    for(jSpecies = 0; jSpecies < nSpecies && jSpecies != iSpecies; ++jSpecies)
        if(jSpecies != iSpecies)
          phi += 1.065*ys_over_mm[jSpecies]/std::sqrt(8.0*(1.0 + mMasses[iSpecies]/mMasses[jSpecies]))*
                      (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25))*
                      (1.0 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25));

      phi += ys_over_mm[iSpecies];
      Thermal_Conductivity_Mixture += Thermal_Conductivities[iSpecies]*ys_over_mm[iSpecies]/phi;
    }
	  return Thermal_Conductivity_Mixture;
  }

  //
  //
  /* This function computes the species diffusion in case of a constant Lewis number */
  ReactingModelLibrary::Vec ReactingModelLibrary::GetRhoUdiff(const double temp, const double rho, const Vec& ys) {
    double kappa = ComputeLambda(temp,ys);
    double Cp = ComputeCP(temp,ys);
    Vec rhoUdiff(nSpecies);
    std::fill(rhoUdiff.begin(),rhoUdiff.end(),kappa/(rho*Cp*Le));
    return rhoUdiff;
  }

  //
  //
  /* This function computes the binary diffusion coefficients for Stefan-Maxwell diffusion with an empirical formula */
  RealMatrix ReactingModelLibrary::GetDij_SM(const double temp, const double pressure) {
    RealMatrix result(nSpecies,nSpecies);
    double Mij,diff_vol_i,molar_mass_i,diff_vol_j;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      molar_mass_i = mMasses[iSpecies];
      diff_vol_i = std::cbrt(Diff_Volumes[iSpecies]);
      for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
        Mij = std::sqrt(2.0*(molar_mass_i*mMasses[jSpecies])/(molar_mass_i + mMasses[jSpecies]));
        diff_vol_j = std::cbrt(Diff_Volumes[jSpecies]);
        result(iSpecies,jSpecies) = 0.0143*std::pow(temp,1.75)/(pressure*Mij*(diff_vol_i + diff_vol_j)*(diff_vol_i + diff_vol_j));
      }
    }
    return result;
  }

  //
  //
  /* This function computes the effective diffusion coefficients for each species */
  ReactingModelLibrary::Vec ReactingModelLibrary::GetDiffCoeffs(const double temp, const double pressure, const Vec& ys) {
    Vec result;
    result.reserve(nSpecies);
    double tmp;

    RealMatrix Dij = GetDij_SM(temp,pressure);
    SetMolarFromMass(ys);
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      tmp = 0.0;
      for(unsigned short jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
        if(iSpecies != jSpecies)
          tmp += Xs[jSpecies]/Dij(iSpecies,jSpecies);
      }
      result.push_back((1.0 - Xs[iSpecies])/tmp);
    }

    return result;
  }

  //
  //
  /* This function computes the species concetration */
  inline void ReactingModelLibrary::SetConcentration(const double rho, const Vec& ys) {
    SU2_Assert(ys.size() == nSpecies, "The dimension of vector with mass fractinos doesn't match nSpecies");
    SU2_Assert(Cs.size() == nSpecies, "The dimension of vector with partial concetrations doesn't match nSpecies");

    std::transform(ys.cbegin(),ys.cend(),mMasses.cbegin(),Cs.begin(),[rho](double ys_elem,double mass_elem){return rho*ys_elem/mass_elem;});
  }

  //
  //
  /* This function computes the equilibrium constants for a specific reaction */
  std::pair<double,double> ReactingModelLibrary::GetKeq(const double temp,const unsigned short iReac) {
    SU2_Assert(iReac < nReactions, "The index of reaction exceeds the number of reactions detected");

    double kf,kb;
    unsigned short iSpecies;

    kf = As[iReac]*std::pow(temp,Ns[iReac])*std::exp(-Temps_Activation[iReac]/(temp));

    if(Elementary_Reactions[iReac])
      kb = 0.0;
    else {
      Vec coeff_reac = Stoich_Coeffs_Reactants.col(iReac);
      Vec coeff_prod = Stoich_Coeffs_Products.col(iReac);
      double dG = 0.0;
      int dnu = 0;
      for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        int dcoeff = coeff_reac[iSpecies] - coeff_prod[iSpecies];
        if(dcoeff != 0)
          dG += dcoeff*
                (MathTools::GetSpline(*std::get<0>(Enth_Spline[iSpecies]),std::get<1>(Enth_Spline[iSpecies]), std::get<2>(Enth_Spline[iSpecies]),temp) -
                 temp*
                 MathTools::GetSpline(*std::get<0>(Entr_Spline[iSpecies]),std::get<1>(Entr_Spline[iSpecies]), std::get<2>(Entr_Spline[iSpecies]),temp));
        dnu += dcoeff;
      }
      const double RT = R_ungas*temp;
      const double Keq = std::exp(-dG/RT)*std::pow(RT,-dnu);
      kb = kf/Keq;
    }

    return std::make_pair(kf,kb);
  }

  //
  //
  /* This function computes the omega term */
  ReactingModelLibrary::Vec ReactingModelLibrary::GetMassProductionTerm(const double temp, const double rho, const Vec& ys) {
    SetConcentration(rho,ys);
    Vec omega;
    omega.reserve(nSpecies);

    double kf_con,kb_con;
    double omega_temp;
    RealVec stoich_reac,stoich_prod,cs_exp,cs_prod;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      omega_temp = 0.0;
      for(unsigned iReac = 0; iReac < nReactions; ++iReac) {
        cs_exp = Cs;
        cs_prod = Cs;
        stoich_reac = Stoich_Coeffs_Reactants.GetRow<RealVec>(iSpecies);
        std::transform(cs_exp.begin(), cs_exp.end(), stoich_reac.cbegin(), cs_exp.begin(),[](double base,double expon){return std::pow(base,expon);});
        kf_con = std::accumulate(cs_exp.cbegin(),cs_exp.cend(),1,std::multiplies<double>());
        stoich_prod = Stoich_Coeffs_Products.GetRow<RealVec>(iSpecies);
        std::transform(cs_prod.begin(), cs_prod.end(), stoich_prod.cbegin(), cs_prod.begin(),[](double base,double expon){return std::pow(base,expon);});
        kf_con = std::accumulate(cs_prod.cbegin(),cs_prod.cend(),1,std::multiplies<double>());
        auto Keq = GetKeq(temp,iReac);
        kf_con *= Keq.first;
        kb_con *= Keq.second;
        omega_temp += (stoich_reac[iReac] - stoich_reac[iReac])*(kf_con - kb_con);
      }
      omega.push_back(omega_temp*mMasses[iSpecies]);
    }
    return omega;
  }

  //
  //
  /* This functions reads the species and sets their order */
  void ReactingModelLibrary::ReadDataMixture(const std::string& f_name) {
    std::string line;
    std::string curr_species;
    unsigned short n_line = 0;

    std::ifstream mixfile(f_name);
    if(mixfile.is_open()) {
      /*--- Clear vectors for safety ---*/
      Species_Names.clear();
      mMasses.clear();
      Formation_Enthalpies.clear();
      Diff_Volumes.clear();
      /*---- Read lines ---*/
      while(mixfile.good() && !mixfile.eof()) {
        std::getline(mixfile,line);
        /*--- Check if we encounter the termination character ---*/
        if(line == "STOP")
          break;
        /*--- We avoid clearly reading comments and empty lines in the file ---*/
        if(!line.empty() && !std::ispunct(line[0])) {
          if(n_line == 0) {
            std::istringstream curr_line(line);
            curr_line>>nSpecies;
            SU2_Assert(!curr_line.fail(),"You have to specify the number of species before proceding");
            mMasses.reserve(nSpecies);
            Formation_Enthalpies.reserve(nSpecies);
            Diff_Volumes.reserve(nSpecies);
            n_line++;
          }
          else {
            SU2_Assert(std::isalpha(line[0]),std::string("Empty species field at line " + std::to_string(n_line)));
            std::string curr_species;
            double curr_mass,curr_enth,curr_vol;
            std::istringstream curr_line(line);
            curr_line>>curr_species;
            SU2_Assert(!curr_line.fail(),std::string("Empty species field at line " + std::to_string(n_line)));
            /*--- Read molar mass ---*/
            curr_line>>curr_mass;
            SU2_Assert(!curr_line.fail(),std::string("The molar mass of species " + curr_species + " is missing"));
            mMasses.push_back(curr_mass);
            /*--- Read formation enthalpy ---*/
            curr_line>>curr_enth;
            SU2_Assert(!curr_line.fail(),std::string("The formation enthalpy of species " + curr_species + " is missing"));
            Formation_Enthalpies.push_back(curr_enth);
            /*--- Read formation enthalpy ---*/
            curr_line>>curr_vol;
            SU2_Assert(!curr_line.fail(),std::string("The diffusion of species " + curr_species + " is missing"));
            Diff_Volumes.push_back(curr_vol);
            /*--- Insert in the map ---*/
            auto res = Species_Names.insert(std::make_pair(curr_species,n_line - 1));
            SU2_Assert(res.second==true,std::string("The species " + curr_species + " has already been declared"));
            n_line++;
          }
        }
        mixfile.close();
        SU2_Assert(Species_Names.size() == nSpecies,"The number of species detected doesn't match nSpecies");

        /*--- Resize vector that wil be often used ---*/
        Ys.resize(nSpecies);
        Xs.resize(nSpecies);
        Cs.resize(nSpecies);
        Viscosities.resize(nSpecies);
        Enthalpies.resize(nSpecies);
        CPs.resize(nSpecies);
        Thermal_Conductivities.resize(nSpecies);
      }
      mixfile.close();
    }
    else {
      std::cerr<<"Unable to open the mixture file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- Reading data about chemistry ---*/
  void ReactingModelLibrary::ReadDataChem(const std::string& f_name) {
    std::string line;
    unsigned n_line = 0;
    unsigned n_reac_read = 0;

    std::ifstream chemfile(f_name);
    if(chemfile.is_open()) {
      std::set<std::string> Species_Reactions;
      /*--- Clear for safety ---*/
      Stoich_Coeffs_Reactants.clear();
      Stoich_Coeffs_Reactants.clear();
      Stoich_Coeffs_Products_Exp.clear();
      Stoich_Coeffs_Reactants_Exp.clear();
      Elementary_Reactions.clear();
      As.clear();
      Ns.clear();
      Temps_Activation.clear();
      while(chemfile.good() && !chemfile.eof()) {
        std::getline(chemfile,line);
        /*--- Check if we encounter the termination character ---*/
        if(line == "STOP")
          break;

        /*--- We avoid clearly reading empty lines and comments in the file ---*/
        if(!line.empty() && !std::ispunct(line[0])) {
          if(n_line == 0) {
            std::istringstream curr_line(line);
            curr_line>>nReactions;
            SU2_Assert(!curr_line.fail(),"You have to specify the number of reactions before proceding");
            Stoich_Coeffs_Reactants.resize(nSpecies,nReactions);
            Stoich_Coeffs_Products.resize(nSpecies,nReactions);
            Stoich_Coeffs_Reactants_Exp.resize(nReactions,nSpecies);
            Stoich_Coeffs_Products_Exp.resize(nReactions,nSpecies);
            Elementary_Reactions.reserve(nReactions);
            As.reserve(nReactions);
            Ns.reserve(nReactions);
            Temps_Activation.reserve(nReactions);
            n_line++;
          }
          else {
            bool is_elementary;
            if(n_line % 2 == 1) {
              is_elementary = line.find('<') != std::string::npos;
              Elementary_Reactions.push_back(is_elementary);
              n_reac_read++;
              ReadReactSpecies(line,is_elementary,n_reac_read,Species_Reactions);
            }
            else {
              Elementary_Reactions.push_back(is_elementary);
              ReadChemCoefs(line);
            }
          }
          n_line++;
        }
      }
      SU2_Assert(n_reac_read == nReactions, "The number of reactions detected doesn't match nReactions");
      SU2_Assert(Species_Reactions.size() == Species_Names.size(),"Some species declared in the mixture are not in the reactions");
      chemfile.close();
    }
    else {
      std::cerr<<"Unable to open the chemical file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- This function reads chemical reactions and store stoichiometric coefficients ---*/
  void ReactingModelLibrary::ReadReactSpecies(const std::string& line, bool is_elem, unsigned n_reac, std::set<std::string>& Species_Reactions) {
    auto minor_pos = line.find('<');
    auto major_pos = line.find('>');
    SU2_Assert(major_pos != std::string::npos,"No reaction in this line");
    SU2_Assert(line.find('>',major_pos+1),"Already detected > symbol for reactions");

    std::string reactants_side,products_side;

    if(is_elem) {
      SU2_Assert(minor_pos+2==major_pos,"Incorrect symbol to detect reactions");
      SU2_Assert(line.find('<',minor_pos+1),"Already detected < symbol for reactions");
      reactants_side = line.substr(0,minor_pos);
    }
    else {
      auto bar_pos = line.find('-');
      SU2_Assert(bar_pos != std::string::npos,"No reaction in this line");
      SU2_Assert(bar_pos+1 == major_pos,"Incorrect symbol to detect reactions");
      reactants_side = line.substr(0,bar_pos);
    }
    products_side = line.substr(major_pos + 1);

    Utility::Parse_Terms(reactants_side,n_reac,Species_Names,Stoich_Coeffs_Reactants,Stoich_Coeffs_Reactants_Exp,Species_Reactions);
    Utility::Parse_Terms(products_side,n_reac,Species_Names,Stoich_Coeffs_Products,Stoich_Coeffs_Products_Exp,Species_Reactions);
  }

  //
  //
  /*--- This function reads coefficients for equilibrium constants ---*/
  void ReactingModelLibrary::ReadChemCoefs(const std::string& line) {
    double A,Ta;
    int n;
    std::istringstream curr_line(line);
    curr_line>>A;
    SU2_Assert(!curr_line.fail(),"No exponential prefactor after a reaction");
    As.push_back(A);
    curr_line>>n;
    SU2_Assert(!curr_line.fail(),"No temperature exponent after a reaction");
    Ns.push_back(n);
    curr_line>>Ta;
    SU2_Assert(!curr_line.fail(),"No activation temperature after a reaction");
    Temps_Activation.push_back(Ta);
  }

  //
  //
  /*--- Reading data about transport properties ---*/
  void ReactingModelLibrary::ReadDataTransp(const std::string& f_name,const unsigned short iSpecies) {
    SU2_Assert(iSpecies < nSpecies, std::string("The index " + std::to_string(iSpecies) + " is out of range"));

    std::string line;
    unsigned n_line = 0;
    double curr_temp,curr_visc,curr_cond;
    std::shared_ptr<RealVec> temp_data(new RealVec());
    RealVec mu_data,kappa_data;

    std::ifstream transpfile(f_name);

    if(transpfile.is_open()) {
      while(transpfile.good() && !transpfile.eof()) {
        std::getline(transpfile,line);
        /*--- We avoid clearly reading empty lines and comments ---*/
        if(!line.empty() && !std::ispunct(line[0])) {
          std::istringstream curr_line(line);
          /*--- Reading temperature ---*/
          curr_line>>curr_temp;
          SU2_Assert(!curr_line.fail(),std::string("Empty Temperature field at line " + std::to_string(n_line + 1) +
                                                   " for species " + std::to_string(iSpecies)));
          temp_data->push_back(curr_temp);
          /*--- Reading viscosity ---*/
          curr_line>>curr_visc;
          SU2_Assert(!curr_line.fail(),std::string("Empty Mu field at line " + std::to_string(n_line + 1) +
                                                   " for species " + std::to_string(iSpecies)));
          mu_data.push_back(curr_visc);
          /*--- Reading conductivity ---*/
          curr_line>>curr_cond;
          SU2_Assert(!curr_line.fail(),std::string("Empty Kappa field at line " + std::to_string(n_line + 1) +
                                                   " for species " + std::to_string(iSpecies)));
          kappa_data.push_back(curr_cond);

          n_line++;
        }
      }

      RealVec y2_mu,y2_kappa;

      MathTools::SetSpline(*temp_data,mu_data,0.0,0.0,y2_mu);
      Mu_Spline[iSpecies] = std::make_tuple(temp_data,std::move_if_noexcept(mu_data),std::move_if_noexcept(y2_mu));

      MathTools::SetSpline(*temp_data,kappa_data,0.0,0.0,y2_kappa);
      Kappa_Spline[iSpecies] = std::make_tuple(temp_data,std::move_if_noexcept(kappa_data),std::move_if_noexcept(y2_kappa));

      transpfile.close();
    }
    else {
      std::cerr<<"Unable to open the species file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- Reading data around thermodynamical properties ---*/
  void ReactingModelLibrary::ReadDataThermo(const std::string& f_name,const unsigned short iSpecies) {
    SU2_Assert(iSpecies < nSpecies, std::string("The index " + std::to_string(iSpecies) + " is out of range"));

    std::string line;
    unsigned n_line = 0;

    double curr_temp,curr_enth,curr_Cp,curr_entr;

    std::shared_ptr<RealVec> temp_data(new RealVec());
    RealVec cp_data,enth_data,entr_data;

    std::ifstream thermofile(f_name);

    if(thermofile.is_open()) {
      while(thermofile.good() && !thermofile.eof()) {
      std::getline(thermofile,line);
      /*--- We clearly avoid reading and empty lines ---*/
      if(!line.empty() && !std::ispunct(line[0])) {
        std::istringstream curr_line(line);
        /*--- Reading temperature ---*/
        curr_line>>curr_temp;
        SU2_Assert(!curr_line.fail(),std::string("Empty Temperature field at line " + std::to_string(n_line + 1) +
                                                 " for species " + std::to_string(iSpecies)));
        temp_data->push_back(curr_temp);
        /*--- Reading Cp ---*/
        curr_line>>curr_Cp;
        SU2_Assert(!curr_line.fail(),std::string("Empty Cp field at line " + std::to_string(n_line + 1) +
                                                  " for species " + std::to_string(iSpecies)));
        cp_data.push_back(curr_Cp);
        /*--- Reading Enthalpy ---*/
        curr_line>>curr_enth;
        SU2_Assert(!curr_line.fail(),std::string("Empty Enthalpy field at line " + std::to_string(n_line + 1) +
                                                 " for species " + std::to_string(iSpecies)));
        enth_data.push_back(curr_enth);

        curr_line>>curr_entr;
        SU2_Assert(!curr_line.fail(),std::string("Empty Entropy field at line " + std::to_string(n_line + 1) +
                                                 " for species " + std::to_string(iSpecies)));
        entr_data.push_back(curr_entr);

        n_line++;
       }
      }

      RealVec y2_cp,y2_enth,y2_entr;

      MathTools::SetSpline(*temp_data,cp_data,0.0,0.0,y2_cp);
      Cp_Spline[iSpecies] = std::make_tuple(temp_data,std::move_if_noexcept(cp_data),std::move_if_noexcept(y2_cp));

      MathTools::SetSpline(*temp_data,enth_data,0.0,0.0,y2_enth);
      Enth_Spline[iSpecies] = std::make_tuple(temp_data,std::move_if_noexcept(enth_data),std::move_if_noexcept(y2_enth));

      MathTools::SetSpline(*temp_data,entr_data,0.0,0.0,y2_entr);
      Entr_Spline[iSpecies] = std::make_tuple(temp_data,std::move_if_noexcept(entr_data),std::move_if_noexcept(y2_entr));

      thermofile.close();
    }
    else {
      std::cerr<<"Unable to open the thermo file: "<<f_name<<std::endl;
      std::exit(1);
    }
  }

  //
  //
  /*--- Setup library ---*/
  void ReactingModelLibrary::Setup(void) {
    if(!Lib_Setup) {
      Le = 1.0;
      /*--- If nobody has configured the library path, we try to do it here with a default value ---*/
      if(Lib_Path=="") {
        std::cout<<"Library path set to default"<<std::endl;
        auto base_dir = std::experimental::filesystem::current_path().string();
        Lib_Path = base_dir + "/../../Common/include";
      }

      std::vector<std::string> list_file;
      std::ifstream config_file(Config_File);
      if(config_file.is_open()) {
        while(config_file.good() && !config_file.eof()) {
          std::string curr_line;
          std::getline(config_file,curr_line);
          if(!curr_line.empty() && !std::ispunct(curr_line[0]))
            list_file.push_back(curr_line);
        }
      }
      else {
        std::cerr<<"Unable to open the specified file with all the file names for setting library."<<std::endl;
        std::exit(1);
      }

      using size_type = std::vector<std::string>::size_type;
      size_type n_file = 2*nSpecies + 2;
      SU2_Assert(list_file.size() == n_file,"The number of files present in the configuration file is wrong");

      std::string file_mix = list_file[0];
      ReadDataMixture(file_mix);
      std::cout<<"Mixture Data read"<<std::endl;
      SetRiGas();

      std::string file_chem = list_file[1];
      ReadDataChem(file_chem);
      std::cout<<"Chemical Reactions read"<<std::endl;
      /*--- We assume that the order of declaration of files is the same as the order of species
            (I can't check the content so it seems reasonable) ---*/
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        std::string file_species = list_file[iSpecies*2 + 2];
        ReadDataTransp(file_species,iSpecies);
        std::cout<<"Transport Data Species "<<iSpecies<<" read"<<std::endl;
        std::string file_thermo = list_file[iSpecies*2 + 3];
        ReadDataThermo(file_thermo,iSpecies);
        std::cout<<"Thermo Data Species "<<iSpecies<<" read"<<std::endl;
      }

      Lib_Setup = true;
      std::cout<<"Library set."<<std::endl;
   }

   else
      throw Common::NotSetup("Trying to setup again without calling unsetup first.");
  }

  //
  //
  /*--- Unsetup the library ---*/
  void ReactingModelLibrary::Unsetup(void) {
    if(Lib_Setup) {
      Species_Names.clear();

      Ri.clear();
      mMasses.clear();
      Diff_Volumes.clear();
      Ys.clear();
      Xs.clear();
      Cs.clear();
      Viscosities.clear();
      Enthalpies.clear();
      CPs.clear();
      Thermal_Conductivities.clear();
      Formation_Enthalpies.clear();

      Stoich_Coeffs_Products.clear();
      Stoich_Coeffs_Products_Exp.clear();
      Stoich_Coeffs_Reactants.clear();
      Stoich_Coeffs_Reactants_Exp.clear();
      As.clear();
      Ns.clear();
      Temps_Activation.clear();
      Elementary_Reactions.clear();

      Enth_Spline.clear();
      Cp_Spline.clear();
      Mu_Spline.clear();
      Kappa_Spline.clear();
      Entr_Spline.clear();

      Lib_Setup = false;
    }

    else
      throw Common::NotSetup("Trying to unsetup without calling setup first.");

  }

} /*--- End of namespace Common ---*/
