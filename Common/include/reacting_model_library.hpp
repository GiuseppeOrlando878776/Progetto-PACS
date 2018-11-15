#ifndef SU2_REACTING_MODEL_LIBRARY
#define SU2_REACTING_MODEL_LIBRARY

#include "physical_chemical_library.hpp"

#include <memory>
#include <vector>
#include <map>
#include <set>
#include <tuple>

namespace Framework {

  /*!
   * \brief Provides a particular library definition to compute the physical and chemical properties.
   */

  class ReactingModelLibrary: public  Framework::PhysicalChemicalLibrary {

  public:

    typedef std::tuple<std::shared_ptr<RealVec>,RealVec,RealVec> MyTuple;

  public: // functions

    /*!
     * \brief Constructor with the name of the library.
     */
    explicit ReactingModelLibrary(const std::string& name):PhysicalChemicalLibrary(),Lib_Name(name),Rgas(),Le() {}


    /*!
     * \brief Default destructor.
     */
    ~ReactingModelLibrary() = default;

    /*!
     * \brief Setups the library name.
     */
    inline void SetLibName(const std::string& library_name) {
      Lib_Name = library_name;
    }

    /*!
     * \brief Check if library is setup.
     */
    inline bool IsSetup(void) const {
      return Lib_Setup;
    }

    /*!
     * \brief Setups the data of the library.
     */
    void Setup(void);

    /*!
     * \brief Unsetups the data of the library.
     */
    void Unsetup(void);

    /*!
     * \brief Set the constant of gases for each species [J/(Kg*K)]
     */
    void SetRiGas(void) override;

    /*!
     * \brief Get the constant of perfect gases [J/(Kg*K)]
     */
    inline su2double GetRgas(void) const override {
      return Rgas;
    }

    /*!
     * \brief Set the constant of perfect gases for the mixture[J/(Kg*K)]
     * \param[in] ys - The vector of the mass fractions of species
     */
    void SetRgas(const RealVec& ys) override;

    /*!
     * Compute the constant of perfect gases for the mixture [J/(Kg*K)]
     * \param[in] ys - The vector of the mass fractions of species
     */
    virtual su2double ComputeRgas(const RealVec& ys) = 0;

    /*!
     * \brief Get the molar masses of the species
     */
    inline void GetMolarMasses(RealVec& mm) const override {
      mm=mMasses;
    }

    /*!
     * \brief Gets the molar fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
     * \param[out] xs - The vector of the molar fractions of species (output)
    */
    void GetMolarFractions(const RealVec& ys, RealVec& xs) override;

    /*!
     * Sets the molar fractions of elements Xn.This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] xn - The vector of the mass fractions of elements
    */
    void SetMolarFractions(const RealVec& xs) override;

    /*!
     * \brief Gets the mass fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
     * \param[out] ys - The vector of the mass fractions of species (output)
     */
    void GetMassFractions(const RealVec& xs, RealVec& ys) override;

    /*!
     * \brief Sets the mass fractions. This function should be called before getting
     * \brief thermodynamic quantities or transport properties.
     * \param[in] ys The vector of the mass fractions of species
    */
    void SetMassFractions(const RealVec& ys) override;

    /*!
     * \brief Computes the frozen specific heat ratio and the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[out] gamma - specific heat ratio (output)
     * \param[out] sound_speed - speed of sound (output)
    */
    void Gamma_FrozenSoundSpeed(const su2double temp, const RealVec& ys,su2double& gamma, su2double& sound_speed) override;

    /*!
     * \brief Computes the density, the enthalpy and the internal energy
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
     * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
     */
    void Density_Enthalpy_Energy(const su2double temp, const su2double pressure, const RealVec& ys, RealVec& dhe) override;

    /*!
     * \brief Computes the density at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    su2double ComputePressure(const su2double temp, const su2double rho, const RealVec& ys) override;

    /*!
     * \brief Computes the density at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
    */
    su2double ComputeDensity(const su2double temp, const su2double pressure, const RealVec& ys) override;

    /*!
     * \brief Computes the internal energy per unit of mass at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
     */
    su2double ComputeEnergy(const su2double temp, const su2double pressure, const RealVec& ys) override;

    /*!
     * \brief Returns the formation enthalpies per unit mass of species
     * \param[out] hs - species formation enthalpies (output)
    */
    inline void GetFormationEnthalpies(RealVec& hs) override {
      hs = Formation_Enthalpies;
    }

    /*!
     * \brief Returns the static enthalpy per unit of mass
     * \param[in] temp - the mixture temperature
     * \param[in] ys - mass fractions
     * \return hsTot - species total enthalpy
    */
    su2double ComputeEnthalpy(const su2double temp,const RealVec& ys) override;

    /*!
      * \brief Computes the specific heat at constant volume
      * \param[in] temp - temperature
      * \param[in] ys - mass fractions
      * \return Specific heat at constant volume
    */
    su2double ComputeCV(const su2double temp,const RealVec& ys) override {
      return ComputeCP(temp,ys) - ComputeRgas(ys);
    }

    /*!
      * \brief Computes the specific heat at constant pressure
      * \param[in] temp - temperature
      * \param[in] ys - mass fractions
      * \return Specific heat at constant pressure
    */
    su2double ComputeCP(const su2double temp,const RealVec& ys) override;

    /*!
     * \brief Computes the thermal conductivity for each species at given temperature
     * \param[in] temp - temperature
     * \return Thermal conductivity for each species
     */
    void ComputeConductivities(const su2double temp);

    /*!
      * \brief Computes the thermal conductivity at given temperature
      * \param[in] temp - temperature
      * \param[in] ys - mass fractions
      * \return Mixture thermal conductivity
      */
    su2double ComputeLambda(const su2double temp,const RealVec& ys) override;

    /*!
      * \brief Computes the dynamic viscosity for each species at given temperature
      * \param[in] temp - temperature
      * \return Molecular viscosity for each species
      */
    void ComputeViscosities(const su2double temp);

    /*!
      * \brief Computes the dynamic viscosity at given temperature
      * \param[in] temp - temperature
      * \param[in] ys - mass fractions
      * \return Molecular viscosity for the mixture in thermal non equilibrium
      */
    su2double ComputeEta(const su2double temp,const RealVec& ys) override;

    /*!
     * \brief Returns the mass production/destruction terms [kg m^-3 s^-1] in CHEMICAL
     * \brief NONEQUILIBRIUM based on Arrhenius's formula.
     * \param[in] pressure - the mixture pressure
     * \param[in] temp - the mixture temperature
     * \param[in] ys - the species mass fractions
     * \param[in] omega - the mass production terms
     * \param[in] jacobian - the Jacobian matrix of the mass production terms
    */
    void GetMassProductionTerm(const su2double temp, const su2double pressure, const su2double rho,
                               RealVec& ys, RealVec& omega, RealMatrix& jacobian) override;

    /*!
     * Returns the source terms species continuity and energy equations
     * \param[in] temp - the mixture temperature
     * \param[in] pressure - the mixture pressure
     * \param[in] ys - the species mass fractions
     * \param[in] rho - the mixture density
     * \param[in] omega - the mass producrtion term
     * \param[in] omegav - the source term
    */
    //void GetSourceTerm(const su2double temp, const su2double pressure, const su2double rho,
    //                   RealVec& ys,RealVec& omega,RealVec& omegav) override;

   /*!
    * Returns the diffusion velocities of species multiplied by the species
    * densities for nonequilibrium computations
    * \param[in] temp - the mixture temperature
    * \param[in] rho  - the mixture density
    * \param[in] ys - the species mass fractions
    * \return Diffusion coefficient with constant Lewis number for each species
    */
    RealVec GetRhoUdiff(const su2double temp, const su2double rho, const RealVec& ys) override;

   /*!
    * Returns the diffusion velocities of species multiplied by the species
    * densities for nonequilibrium computations
    * \param[in] temp - the mixture temperature
    * \param[in] pressure - the mixture pressure
    */
    void GetDij_fick(RealVec& dx, const su2double pressure, const su2double temperature, RealMatrix& Dij) override;

  private:

    /*!
    * \brief Returns the reaction contribution to source term of species continuity equations.
    * \param[in] temp - the mixture temperature
    * \param[in] rho the mixture density
    * \param[in] ys - the species mass fractions
    * \param[in] mMasses -  molar masses
    * \param[out] omega the mass production term
    */
    void OmegaContribution(const su2double temp,const su2double rho,
                           const RealVec& ys,const RealVec& mMasses,RealVec& omega);

   /*!
    * \brief Returns the forward reaction rate coefficient.
    * \param[in] temp - the corrected temperature Tq
    */
    su2double Kf(const su2double temp);

    /*!
     * \brief Returns the backward reaction rate coefficient.
     * \param[in] temp the mixture temperature
    */
    su2double Kb(const su2double temp);

    /*!
      * \brief Preprocess to generate data tables for thermodynamic/transport properties
      * \param[in] f_name - name of the file with the properties
      * \param[in] iSpecies - index of the desired species
    */
    void Preprocess(const std::string& f_name,const unsigned short iSpecies);

    /*!
      * \brief Read mixture data
      */
    void ReadDataMixture(const std::string& f_name);

    /*!
    * \brief Read chemistry data
    * \param[in] f_name - name of the file with the nvolved reactions
    */
    void ReadDataChem(const std::string& f_name);

    /*!
      * \brief Read transport data
      * \param[in] f_name - name of the file with the properties
      * \param[in] iSpecies - index of the desired species
    */
    void ReadDataTransp(const std::string& f_name, const unsigned short iSpecies);

    /*!
      * \brief Read thermodynamical data
      * \param[in] f_name - name of the file with the properties
      * \param[in] iSpecies - index of the desired species
    */
    void ReadDataThermo(const std::string& f_name, const unsigned short iSpecies);

    /*!
      * \brief Read species involved in a chemical reaction from line
      * \param[in] line - line to be read
      * \param[out] Species_Reactions - names of species form reactions
    */
    void ReadReactSpecies(const std::string& line,std::set<std::string>& Species_Reactions);

    /*!
      * Read coefficients to compute reaction rates from line
      * \param[in] idx - index of position in "line" separating reaction definition and coefficients
      * \param[in] ChemCoefs - Coefficients to compute reaction rates
      * \param[in] line - line to be read
    */
    void ReadChemCoefs(const std::string& line);

  protected:

    //bool Lib_Setup; /*!\brief Flag to check setting library */

    std::string Lib_Name;  /*!< \brief Name of the library. */

    //unsigned short Butidene_ID; /*!< \brief Identification chemical model. */

    std::string Config_File; /*!< \brief Name of the file for reading file names thermodynamical and transport properties. */
                              // where to add it?

    std::map<std::string,unsigned short> Species_Names;  /*!< \brief Names of species in the mixture. */

    su2double Rgas; /*!< \brief Gas constant of the mixture. */

    su2double Le; /*!< \brief Lewis number. */

    //su2double Lam_Pr; /*!< \brief Laminar Prandtl number. */

    //su2double T_ref; /*!< \brief Reference temperature. */

    //su2double Viscosity_Mixture;  /*!< \brief Viscosity of the mixture. */

    //su2double Thermal_Conductivity_Mixture; /*!< \brief Laminar Prandtl number. */

    //su2double Cp_Mixture; /*!< \brief Specific heat constant pressure for the mxiture. */

    //su2double Gamma; /*!< \brief Specific heat ratio. */

    //su2double Sound_Speed; /*!< \brief Frozen Sound speed. */

    //su2double AF; /*!< \brief Acentric factor. */

    //su2double Tc; /*!< \brief Critical temperature. */

    //su2double Vc; /*!< \brief Critical volume. */

    //su2double EpsKb; /*!< \brief Critical volume. */

    RealVec mMasses; /*!< \brief Molar mass for each species. */

    //RealVec Tref;   /*!< \brief Reference tempeature for each species. */

    RealVec Ri; /*!< \brief Specific gas constant for each species. */

    RealVec Ys;    /*!<  \brief Mass fraction for each species. */

    RealVec Xs;    /*!<  \brief Molar fraction for each species. */

    RealVec Viscosities; /*!< \brief Viscosity for each species. */

    RealVec Thermal_Conductivities; /*!< \brief Thermal conductivity for each species. */

    std::map<unsigned short,MyTuple> Mu_Spline; /*!< \brief Spline interpolation coefficient for viscosity computation. */

    std::map<unsigned short,MyTuple> Kappa_Spline; /*!< \brief Spline interpolation coefficient for conductivity computation. */

    std::map<unsigned short,MyTuple> Entr_Spline; /*!< \brief Spline interpolation coefficient for entropy computation. */

    std::map<unsigned short,MyTuple> Cp_Spline; /*!< \brief Spline interpolation coefficient for specific heat at constant pressure computation. */

    std::map<unsigned short,MyTuple> Enth_Spline; /*!< \brief Spline interpolation coefficient for enthalpy computation. */

    //RealVec Internal_Energies; /*!< \brief Internal energy for each species. */

    RealVec Enthalpies; /*!< \brief Enthalpy for each species. */

    //RealVec Heat_Capacities; /*!< \brief Heat Capacity for each species. */

    RealVec CPs; /*!< \brief Specific heat at constant pressure for each species (Cp). */

    //RealVec CVs; /*!< \brief Specific heat at constant volume for each species (Cv). */

    //RealVec Sound_Speeds; /*!< \brief Sound speed for each species. */

    RealVec Formation_Enthalpies; /*!< \brief Formation enthalpy for each species. */

    RealVec Stoich_Coeffs_Reactants; /*!< \brief Stochiometric coefficents vector. */

    RealVec Stoich_Coeffs_Products; /*!< \brief Stochiometric coefficents vector. */

    RealVec Stoich_Coeffs_Reactants_Exp; /*!< \brief Stochiometric coefficents vector. */

    RealVec Stoich_Coeffs_Products_Exp; /*!< \brief Stochiometric coefficents vector. */

    RealVec Chem_Coeffs;     /*!< \brief Vector with coefficients to estimate reaction rates. */

  }; /*-- End of class ReactingModelLibrary ---*/

} /*-- End of Namespace Framework ---*/

#endif
