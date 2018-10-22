#include "../include/reacting_model_library.hpp"
#include "../include/not_setup_exception.hpp"
#include "../include/spline.hpp"
#include "../include/utility.hpp"

#include <experimental/filesystem>
//#include <locale>
#include <fstream>
#include <sstream>
#include <exception>
#include <algorithm>
#include <limits>

namespace Framework {

//
//
/*--- Setting gas constant for each species ---*/
void ReactingModelLibrary::SetRiGas(void) {
  SU2_Assert(mMasses.size()==nSpecies,"The dimension of vector mMasses doesn't match nSpecies");
  Ri.resize(nSpecies);
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Ri[iSpecies] = R_ungas/mMasses[iSpecies];
}

//
//
/*--- Setting gas constant for mixture ---*/
void ReactingModelLibrary::SetRgas(const RealVec& ys) {
  SU2_Assert(ys.size()==nSpecies,"The dimension of vector rs doesn't match nSpecies");
  Rgas = std::inner_product(ys.cbegin(),ys.cend(),Ri.cbegin(),0.0);
}

//
//
/*--- Computing gas constant for mixture ---*/
su2double ReactingModelLibrary::ComputeRgas(const RealVec& ys) {
  SetRgas(ys);
  return GetRgas();
}

//
//
/*--- Setting molar fraction for each species ---*/
void ReactingModelLibrary::SetMolarFractions(const RealVec& xs) {
  SU2_Assert(xs.size() == nSpecies,"The dimension of vector xs doesn't match nSpecies");
  Xs = xs;

  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(Xs[iSpecies] < 0.0)
      Xs[iSpecies] = 1E-20;

    SU2_Assert(Xs[iSpecies]<=1.0, std::string("The molar fraction of species number " + std::to_string(iSpecies) + "is greater than 1"));
  }
}

//
//
/*--- Setting moleculesIds for each species ---*/
/*
void ReactingModelLibrary::SetMoleculesIDs(std::vector<unsigned>& Ids) {
  Ids.reserve(nSpecies);

  for(auto i=0;i<nSpecies;++i) {
    if(Atoms[i]>1)
      Ids.push_back(i);
  }
}
*/

//
//
/*--- Setting mass fraction for each species ---*/
void ReactingModelLibrary::SetMassFractions(const RealVec& ys) {
  SU2_Assert(ys.size() == nSpecies,"The dimension of vector ys doesn't match nSpecies");
  Ys = ys;

  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    if(Ys[iSpecies] < 0.0)
      Ys[iSpecies] = 1E-20;

    SU2_Assert(Ys[iSpecies]<=1.0,std::string("The mass fraction of species number " + std::to_string(iSpecies) + "is greater than 1"));
  }
}

//
//
/* This file sets the molar fractions from mass fractions */
void ReactingModelLibrary::GetMolarFractions(const RealVec& ys, RealVec& xs) {
  SU2_Assert(ys.size() == nSpecies,"The dimension of vector ys doesn't match nSpecies");
  xs.clear();
  xs.reserve(nSpecies);
  su2double massTot=0.0;

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    SU2_Assert (ys[iSpecies] <= 1.0,std::string("The mass fraction of species number " + std::to_string(iSpecies) + "is greater than one"));

    const su2double mm = ys[iSpecies]/mMasses[iSpecies];
    massTot += mm;
    xs.push_back(mm);
  }
  std::for_each(xs.begin(),xs.end(),[massTot](double& elem){elem /= massTot;});
}

//
//
/* This file sets the mass fractions from molar fractions */
void ReactingModelLibrary::GetMassFractions(const RealVec& xs, RealVec& ys) {
  SU2_Assert(xs.size() == nSpecies,"The dimension of vector xs doesn't match nSpecies");
  ys.clear();
  ys.reserve(nSpecies);
  su2double massTot=0.0;

  for (unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    SU2_Assert(xs[iSpecies] <= 1.0,std::string("The molar fraction of species number " + std::to_string(iSpecies) + "is greater than 1"));
    const su2double mm = xs[iSpecies] * mMasses[iSpecies];
    massTot += mm;
    ys.push_back(mm);
  }
  std::for_each(ys.begin(),ys.end(),[massTot](double& elem){elem /= massTot;});

}

//
//
/* This function computes the specific heat at constant pressure */
su2double ReactingModelLibrary::ComputeCP(const su2double temp) {
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    CPs[iSpecies] = MathTools::GetSpline(*std::get<0>(Cp_Spline[iSpecies]),std::get<1>(Cp_Spline[iSpecies]),
                                       std::get<2>(Cp_Spline[iSpecies]),temp);

  SU2_Assert(Ys.size() == nSpecies,"The dimension of vector Ys doesn't match nSpecies");
  return std::inner_product(Ys.begin(),Ys.end(),CPs.begin(),0.0);

}

//
//
/* This function computes gamma and the frozen sound speed */
void ReactingModelLibrary::Gamma_FrozenSoundSpeed(const su2double temp,su2double& gamma,su2double& sound_speed) {
  RealVec CPs_Over_MM;
  SU2_Assert(CPs.size() == nSpecies,"The dimension of vector CPs doesn't match nSpecies");
  SU2_Assert(mMasses.size() == nSpecies,"The dimension of vector mMasses doesn't match nSpecies");
  std::transform(CPs.cbegin(),CPs.cend(),mMasses.cbegin(),std::back_inserter(CPs_Over_MM),std::divides<su2double>());
  SU2_Assert(Ys.size() == nSpecies,"The dimension of vector Ys doesn't match nSpecies");
  su2double Cp = std::inner_product(Ys.cbegin(),Ys.cend(),CPs_Over_MM.cbegin(),0.0);
  su2double Cv = Cp - Rgas;
  gamma = Cp/Cv;
  sound_speed = std::sqrt(gamma*Rgas*temp);
}

//
//
/*--- Computing molecular viscosity of each species---*/
void ReactingModelLibrary::ComputeViscosities(const su2double temp) {
  SU2_Assert(Mu_Spline.size() == nSpecies,"The dimension of Mu_Spline doesn't match nSpecies");
  SU2_Assert(Viscosities.size() == nSpecies,"The dimension of Viscosities doesn't match nSpecies");
  for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Viscosities[iSpecies] = MathTools::GetSpline(*std::get<0>(Mu_Spline[iSpecies]),std::get<1>(Mu_Spline[iSpecies]),
                                                  std::get<2>(Mu_Spline[iSpecies]),temp);
}

//
//
/*--- Computing viscosity of the mixture---*/
su2double ReactingModelLibrary::Eta(const su2double temp) {
  ComputeViscosities(temp);
  unsigned short iSpecies, jSpecies;
  su2double phi;
  // Mixture viscosity calculation as sum weighted over PHI.
  su2double Viscosity_Mixture = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies)  {
      phi=0.0;
      for(jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
        phi +=  Xs[jSpecies]*std::pow(1+std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25),2) /
                             std::sqrt(8*(1 + mMasses[iSpecies]/mMasses[jSpecies]));

        Viscosity_Mixture += Xs[iSpecies]*Viscosities[iSpecies]/phi;
  }

  return Viscosity_Mixture;

}


//
//
/*--- Computing thermal conductivity of each species ---*/
void ReactingModelLibrary::ComputeConductivities(const su2double temp) {
  SU2_Assert(Mu_Spline.size() == nSpecies,"The dimension of Mu_Spline doesn't match nSpecies");
  SU2_Assert(Thermal_Conductivities.size() == nSpecies,"The dimension of Thermal_Conductivities doesn't match nSpecies");
    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
    Thermal_Conductivities[iSpecies] = MathTools::GetSpline(*std::get<0>(Kappa_Spline[iSpecies]),std::get<1>(Kappa_Spline[iSpecies]),
                                                            std::get<2>(Kappa_Spline[iSpecies]),temp);
}


//
//
/*--- Computing thermal conductivity of the mixtures---*/
su2double ReactingModelLibrary::Lambda(const su2double temp) {
  ComputeConductivities(temp);
  ComputeViscosities(temp);
  su2double Thermal_Conductivity_Mixture = 0.0;
  unsigned short iSpecies,jSpecies;
  su2double phi;
  for(iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
    phi = 0.0;
	  for(jSpecies = 0; jSpecies < nSpecies and jSpecies != iSpecies; ++jSpecies)
      phi += Xs[jSpecies]*std::pow(1 + std::sqrt(Viscosities[iSpecies]/Viscosities[jSpecies])*std::pow(mMasses[jSpecies]/mMasses[iSpecies],0.25),2) /
                          std::sqrt(8.0*(1 + mMasses[iSpecies]/mMasses[jSpecies]));

    Thermal_Conductivity_Mixture += Thermal_Conductivities[iSpecies]/(1.0 + phi/Xs[iSpecies]);
  }

	return Thermal_Conductivity_Mixture;

}

su2double ReactingModelLibrary::GetRhoUdiff(const su2double temp,const su2double rho) {
  su2double kappa = Lambda(temp);
  su2double Cp = ComputeCP(temp);
  return kappa/(rho*Cp*Le);
}

//
//
/* This file reads the species and sets their order */
void ReactingModelLibrary::ReadDataMixture(const std::string& f_name) {
  Species_Names.clear();
  mMasses.clear();
  Formation_Enthalpies.clear();

  std::string line;
  std::string curr_species;
  double curr_enth,curr_mass;
  unsigned short n_line = 0;

  std::ifstream mixfile(f_name);
  if(mixfile.is_open()) {
    while(mixfile.good() and !mixfile.eof()) {
      // We assume to have a titled table ( to be fixed)
      mixfile.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
      std::getline(mixfile,line);
      // We avoid clearly reading comments and empty lines in the file
      if(!line.empty() and !std::ispunct(line[0])) {
        SU2_Assert(std::isalpha(line[0]),std::string("Empty species field at line " + std::to_string(n_line + 1)));
        std::istringstream curr_line(line);
        curr_line>>curr_species;
        SU2_Assert(!curr_line.fail(),std::string("Empty species field at line " + std::to_string(n_line + 1)));
        auto res = Species_Names.insert(std::make_pair(curr_species,n_line));
        SU2_Assert(res.second==true,std::string("The species " + curr_species + " has already been declared"));
        curr_line>>curr_mass;
        SU2_Assert(!curr_line.fail(),std::string("The molar mass of species " + curr_species + " is missing"));
        mMasses.push_back(curr_mass);
        curr_line>>curr_enth;
        Formation_Enthalpies.push_back(curr_enth);
        n_line++;
      }
      nSpecies = n_line;

      // Resize vector that wil be often used
      Ri.resize(nSpecies);
      Ys.resize(nSpecies);
      Xs.resize(nSpecies);
      Viscosities.resize(nSpecies);
      //Internal_Energies.resize(nSpecies);
      Enthalpies.resize(nSpecies);
      //Heat_Capacities.resize(nSpecies);
      CPs.resize(nSpecies);
      //CVs.resize(nSpecies);
      Thermal_Conductivities.resize(nSpecies);
      //Sound_Speeds.resize(nSpecies);
      //Atoms.resize(nSpecies);

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

  std::ifstream chemfile(f_name);
  if(chemfile.is_open()) {
    std::set<std::string> Species_Reactions;
    while(chemfile.good() and !chemfile.eof()) {
      std::getline(chemfile,line);
      // We avoid clearly reading empty lines and comments in the file
      if(!line.empty() and !std::ispunct(line[0])) {
          std::istringstream curr_line(line);
          if(n_line % 2 == 0)
            ReadReactSpecies(line,Species_Reactions);
          else
            ReadChemCoefs(line);

          n_line++;
      }
    }
    SU2_Assert(Species_Reactions.size() == Species_Names.size(),"Some species declared in the mixture are not in the reactions");
    SetNReactions(n_line);
    chemfile.close();
  }
  else {
    std::cerr<<"Unable to open the chemical file: "<<f_name<<std::endl;
    std::exit(1);
  }
}

void ReactingModelLibrary::ReadReactSpecies(const std::string& line,std::set<std::string>& Species_Reactions) {
  auto minor_pos = line.find('<');
  bool is_reversible = minor_pos!= std::string::npos;

  auto major_pos = line.find('>');
  SU2_Assert(major_pos != std::string::npos,"No reaction in this line");
  SU2_Assert(line.find('>',major_pos+1),"Already detected > symbol for reactions");

  std::string reactants_side,products_side;

  if(is_reversible) {
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

  Stoich_Coeffs_Reactants.clear();
  Stoich_Coeffs_Reactants.clear();
  Stoich_Coeffs_Products_Exp.clear();
  Stoich_Coeffs_Reactants_Exp.clear();
  Utility::Parse_Terms(reactants_side,Species_Names,Stoich_Coeffs_Reactants,Stoich_Coeffs_Reactants_Exp,Species_Reactions);
  Utility::Parse_Terms(products_side,Species_Names,Stoich_Coeffs_Products,Stoich_Coeffs_Products_Exp,Species_Reactions);

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
    while(transpfile.good() and !transpfile.eof()) {
      std::getline(transpfile,line);
      // We avoid clearly reading empty lines and comments
      if(!line.empty() and std::isdigit(line[0])) {
        std::istringstream curr_line(line);
        curr_line>>curr_temp;
        SU2_Assert(!curr_line.fail(),std::string("Empty Temperature field at line " + std::to_string(n_line + 1) +
                                                 " for species " + std::to_string(iSpecies)));
        temp_data->push_back(curr_temp);
        curr_line>>curr_visc;
        SU2_Assert(!curr_line.fail(),std::string("Empty Mu field at line " + std::to_string(n_line + 1) +
                                                  " for species " + std::to_string(iSpecies)));
        mu_data.push_back(curr_visc);
        curr_line>>curr_cond;
        SU2_Assert(!curr_line.fail(),std::string("Empty Kappa field at line " + std::to_string(n_line + 1) +
                                                 " for species " + std::to_string(iSpecies)));
        kappa_data.push_back(curr_cond);

        n_line++;
      }
    }
    RealVec y2_mu,y2_kappa;

    MathTools::SetSpline(*temp_data,mu_data,1.0,1.0,y2_mu);
    Mu_Spline[iSpecies] = std::make_tuple(temp_data,std::move( mu_data),std::move(y2_mu));

    MathTools::SetSpline(*temp_data,kappa_data,1.0,1.0,y2_kappa);
    Kappa_Spline[iSpecies] = std::make_tuple(temp_data,std::move(kappa_data),std::move(y2_kappa));

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
    while(thermofile.good() and !thermofile.eof()) {
      std::getline(thermofile,line);
      // We clearly avoid reading and empty lines
      if(!line.empty() and std::isdigit(line[0])) {
        std::istringstream curr_line(line);
        curr_line>>curr_temp;
        SU2_Assert(!curr_line.fail(),std::string("Empty Temperature field at line " + std::to_string(n_line + 1) +
                                                 " for species " + std::to_string(iSpecies)));
        temp_data->push_back(curr_temp);
        curr_line>>curr_Cp;
        SU2_Assert(!curr_line.fail(),std::string("Empty Cp field at line " + std::to_string(n_line + 1) +
                                                  " for species " + std::to_string(iSpecies)));
        cp_data.push_back(curr_Cp);
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

    MathTools::SetSpline(*temp_data,cp_data,1.0,1.0,y2_cp);
    Cp_Spline[iSpecies] = std::make_tuple(temp_data,std::move(cp_data),std::move(y2_cp));

    MathTools::SetSpline(*temp_data,enth_data,1.0,1.0,y2_cp);
    Enth_Spline[iSpecies] = std::make_tuple(temp_data,std::move(enth_data),std::move(y2_enth));

    MathTools::SetSpline(*temp_data,entr_data,1.0,1.0,y2_cp);
    Entr_Spline[iSpecies] = std::make_tuple(temp_data,std::move(entr_data),std::move(y2_entr));

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

    //AF = 0.192;
    //T_ref = 298.16;
    Le = 1.0;

    // if nobody has configured the library path, we try to do it here with a default value
    if(Lib_Path=="") {
      std::cout<<"Library path set to default"<<std::endl;
      auto base_dir = std::experimental::filesystem::current_path().string();
      Lib_Path = base_dir + "../../Common/include";
    }

    std::vector<std::string> list_file;
    std::ifstream config_file(Config_File);
    if(config_file.is_open()) {
      while(config_file.good() and !config_file.eof()) {
        std::string curr_line;
        std::getline(config_file,curr_line);
        if(curr_line.size()!=0 and !std::ispunct(curr_line[0]))
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

    std::string file_chem = list_file[1];
    ReadDataChem(file_chem);
    std::cout<<"Chemical Reactions read"<<std::endl;

    for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      std::string file_species = list_file[iSpecies*2 + 2];
      ReadDataTransp(file_species,iSpecies);
      std::cout<<"Transport Data Species "<<iSpecies<<" read"<<std::endl;
      std::string file_thermo = list_file[iSpecies*2 + 3];
      ReadDataThermo(file_thermo,iSpecies);
      std::cout<<"Thermo Data Species "<<iSpecies<<" read"<<std::endl;
    }

    SetRiGas();
    //SetMassFractions(Ys);
    //GetMolarFractions(Ys,Xs);
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
    //Atoms.clear();
    Ri.clear();
    mMasses.clear();
    Ys.clear();
    Xs.clear();
    Viscosities.clear();
    //Internal_Energies.clear();
    Enthalpies.clear();
    //Heat_Capacities.clear();
    CPs.clear();
    //CVs.clear();
    Thermal_Conductivities.clear();
    //Sound_Speeds.clear();
    Formation_Enthalpies.clear();

    Stoich_Coeffs_Products.clear();
    Stoich_Coeffs_Products_Exp.clear();
    Stoich_Coeffs_Reactants.clear();
    Stoich_Coeffs_Reactants_Exp.clear();
    Chem_Coeffs.clear();

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


}
