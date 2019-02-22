#ifndef SU2_REACTING_MODEL_LIBRARY
#define SU2_REACTING_MODEL_LIBRARY

#include "physical_chemical_library.hpp"

#include "../../externals/Eigen/Dense"
#include "su2_assert.hpp"
#include <memory>
#include <map>
#include <set>
#include <tuple>

namespace Framework {

  /*!
   * \brief Provides a particular library definition to compute the physical and chemical properties.
   */
  class ReactingModelLibrary: public Framework::PhysicalChemicalLibrary<RealVec,Eigen::MatrixXd> {

  public:

    typedef std::tuple<std::shared_ptr<RealVec>,RealVec,RealVec> MyTuple;
    using RealMatrix = Eigen::MatrixXd;
    using Vec = Eigen::VectorXd;

  public:

    /*!
     * \brief Constructor with the name of the comfiguration file of the library.
     */
    ReactingModelLibrary(const std::string& name): PhysicalChemicalLibrary(name), Rgas(), Le() {}

    /*!
     * \brief Default destructor.
     */
    ~ReactingModelLibrary() = default;

    /*!
     * \brief Check if library is setup.
     */
    inline bool IsSetup(void) const {
      return PhysicalChemicalLibrary<RealVec,RealMatrix>::Lib_Setup;
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
     * \brief Set the constant of gases for each species [J/(Kmol*K)]
     */
    void SetRiGas(void) override;

    /*!
     * \brief Get the constant of perfect gases [J/(Kmol*K)]
     */
    inline double GetRgas(void) const override {
      return Rgas;
    }

    /*!
     * \brief Set the constant of perfect gases for the mixture[J/(Kmol*K)]
     * \param[in] ys - The vector of the mass fractions of species
     */
    void SetRgas(const RealVec& ys) override;

    /*!
     * Compute the constant of perfect gases for the mixture [J/(Kmol*K)]
     * \param[in] ys - The vector of the mass fractions of species
     */
    double ComputeRgas(const RealVec& ys) override;

    /*!
     * \brief Get the molar masses of the species
     */
    inline RealVec GetMolarMasses(void) const override {
      SU2_Assert(mMasses.size() == nSpecies,"The number of elements in the vector of molar masses doesn't match nSpecies");
      return mMasses;
    }

    /*!
     * \brief Get the molar fractions from the mass fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    RealVec GetMolarFromMass(const RealVec& ys) override;

    /*!
     * Set the molar fractions of elements Xs. This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] xs - The vector of the molar fractions of elements
    */
    void SetMolarFractions(const RealVec& xs) override;

    /*!
     * \brief Set the molar fractions from mass fractions.
     * \param[in] ys - The vector of the mass fractions of species (input)
    */
    void SetMolarFromMass(const RealVec& ys) override;

    /*!
     * \brief Get the mass fractions from the molar fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
     */
    RealVec GetMassFromMolar(const RealVec& xs) override;

    /*!
     * \brief Set the mass fractions. This function should be called before getting
     * \brief thermodynamic quantities or transport properties.
     * \param[in] ys The vector of the mass fractions of species
    */
    void SetMassFractions(const RealVec& ys) override;

    /*!
     * \brief Get the mass fractions from molar fractions.
     * \param[in] xs - The vector of the molar fractions of species (input)
    */
    void SetMassFromMolar(const RealVec& xs) override;

    /*!
     * \brief Computes the frozen specific heat ratio and the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[out] gamma - specific heat ratio (output)
     * \param[out] sound_speed - speed of sound (output)
    */
    void Gamma_FrozenSoundSpeed(const double temp, const RealVec& ys, double& gamma, double& sound_speed) override;

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     */
    double ComputeFrozenGamma(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
    */
    double ComputeFrozenSoundSpeed(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The vector of the mass fractions of species
     * \param[in] press - pressure
     * \param[in] rho - density
    */
    double ComputeFrozenSoundSpeed(const double temp, const RealVec& ys, const double press, const double rho) override;

    /*!
     * \brief Computes the density, the enthalpy and the internal energy
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
     * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
     */
    void Density_Enthalpy_Energy(const double temp, const double pressure, const RealVec& ys, RealVec& dhe) override;

    /*!
     * \brief Computes the density at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    double ComputePressure(const double temp, const double rho, const RealVec& ys) override;

    /*!
     * \brief Computes the density at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
    */
    double ComputeDensity(const double temp, const double pressure, const RealVec& ys) override;

    /*!
     * \brief Computes the internal energy per unit of mass at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     */
    double ComputeEnergy(const double temp, const RealVec& ys) override;

    /*!
     * \brief Returns the formation enthalpies per unit mass of species
    */
    inline RealVec GetFormationEnthalpies(void) const override {
      return Formation_Enthalpies;
    }

    /*!
     * \brief Returns the static enthalpy per unit of mass
     * \param[in] temp - the mixture temperature
     * \param[in] ys - mass fractions
     * \return Mixture static enthalpy
    */
    double ComputeEnthalpy(const double temp, const RealVec& ys) override;

    /*!
     * \brief Returns the static enthalpy per unit of mass of each species
     * \param[in] temp - the mixture temperature
    */
    RealVec ComputePartialEnthalpy(const double temp) override;

    /*!
     * \brief Computes the mixture total concentration
     * \param[in] rho - density
     */
    inline double ComputeConcentration(const double rho, const RealVec& ys) override {
      double totMass = 0.0;
      for(unsigned short iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        totMass += ys[iSpecies]/mMasses[iSpecies];
      return rho*totMass;
    }

    /*!
     * \brief Computes the specific heat at constant pressure
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Specific heat at constant pressure
    */
    double ComputeCP(const double temp,const RealVec& ys) override;

    /*!
     * \brief Computes the specific heat at constant volume
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Specific heat at constant volume
    */
    double ComputeCV(const double temp, const RealVec& ys) override {
      return ComputeCP(temp,ys) - ComputeRgas(ys);
    }

    /*!
     * \brief Computes the thermal conductivity for each species at given temperature
     * \param[in] temp - temperature
     * \return Thermal conductivity for each species
    */
    void ComputeConductivities(const double temp);

    /*!
     * \brief Computes the thermal conductivity at given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Mixture thermal conductivity
    */
    double ComputeLambda(const double temp, const RealVec& ys) override;

    /*!
     * \brief Computes the dynamic viscosity for each species at given temperature
     * \param[in] temp - temperature
     * \return Molecular viscosity for each species
    */
    void ComputeViscosities(const double temp);

    /*!
     * \brief Computes the dynamic viscosity at given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Molecular viscosity for the mixture in thermal non equilibrium
    */
    double ComputeEta(const double temp, const RealVec& ys) override;

    /*!
     * Returns the mass production/destruction terms [kg m^-3 s^-1] in chemical
     * non-equilibrium based on Arrhenius's formula.
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \param[in] ys - the species mass fractions
     * \return Mass production terms
    */
    RealVec GetMassProductionTerm(const double temp, const double rho, const RealVec& ys) override;

   /*!
    * Returns the diffusion velocities of species multiplied by the species
    * densities for nonequilibrium computations
    * \param[in] temp - the mixture temperature
    * \param[in] rho  - the mixture density
    * \param[in] ys - the species mass fractions
    * \return Diffusion coefficient with constant Lewis number for each species
    */
    RealVec GetRhoUdiff(const double temp, const double rho, const RealVec& ys) override;

   /*!
    * Returns the binary diffusion coefficients
    * \param[in] pressure - the mixture pressure
    * \param[in] temp - the mixture temperature
    */
    RealMatrix GetDij_SM(const double pressure, const double temp) override;

    /*!
     * Returns thematrix of Stefan-Maxwell equations
     * \param[in] rho - the mixture density
     * \param[in] xs - current molar fractions
     * \param[in] ys - current mass fractions
     * \param[in] val_Dij - current binary diffusion coefficients
     */
    RealMatrix GetGamma(const double rho, const RealVec& xs, const RealVec& ys, const RealMatrix& val_Dij) override;

    /*!
     * \brief Returns the effective diffusion coefficients to solve Stefan-Maxwell equation
     * \param[in] temp - the mixture temperature
     * \param[in] pressure - the mixture pressure
     * \param[in] ys - mass fractions in the mixture
     */
    RealVec GetDiffCoeffs(const double temp, const double pressure, const RealVec& ys) override;

  private:

    /*!
     * \brief Set the actual concetration for each species
     * \param[in] rho - the mixture density
     * \param[in] ys - the actual mass fractions
    */
    void SetConcentration(const double rho, const RealVec& ys);

    /*!
     * \brief Returns the forward and backward reaction rate coefficients.
     * \param[in] temp - the mixture temperature
     * \param[in] iReac - index of the desired reaction
    */
    std::pair<double,double> GetKeq(const double temp, const unsigned short iReac);

    /*!
      * \brief Read transport data
      * \param[in] f_name - name of the file with the properties
    */
    void ReadDataTransp(const std::string& f_name);

    /*!
      * \brief Read thermodynamical data
      * \param[in] f_name - name of the file with the properties
    */
    void ReadDataThermo(const std::string& f_name);

    /*!
      * \brief Read mixture data
      * \param[in] f_name - name of the file with the involved species
      */
    void ReadDataMixture(const std::string& f_name);

    /*!
    * \brief Read chemistry data
    * \param[in] f_name - name of the file with the involved reactions
    */
    void ReadDataChem(const std::string& f_name);

    /*!
      * \brief Read species involved in a chemical reaction from line
      * \param[in] line - line to be read
      * \param[in] is_elem - flag whether the reaction is elementary or not
      * \parm[in]  n_reac - the cuurent number of reactions detected
      * \param[out] Species_Reactions - names of species form reactions
    */
    void ReadReactSpecies(const std::string& line, bool is_elem, unsigned n_reac);

    /*!
      * Read coefficients to compute reaction rates from line
      * \param[in] line - line to be read
    */
    void ReadChemCoefs(const std::string& line);

  protected:

    std::map<std::string,unsigned short> Species_Names;  /*!< \brief Names of species in the mixture. */

    double Rgas; /*!< \brief Gas constant of the mixture. */

    double Le; /*!< \brief Lewis number. */

    RealVec mMasses; /*!< \brief Molar mass for each species. */

    RealVec Diff_Volumes; /*!< \brief Molecular diffusion volume for each species. */

    RealVec Ri; /*!< \brief Specific gas constant for each species. */

    RealVec Ys;    /*!<  \brief Mass fraction for each species. */

    RealVec Xs;    /*!<  \brief Molar fraction for each species. */

    RealVec Cs;    /*!<  \brief Actual concentration for each species. */

    RealVec Viscosities; /*!< \brief Viscosity for each species. */

    RealVec Thermal_Conductivities; /*!< \brief Thermal conductivity for each species. */

    std::vector<MyTuple> Mu_Spline; /*!< \brief Spline interpolation coefficient for viscosity computation. */

    std::vector<MyTuple> Kappa_Spline; /*!< \brief Spline interpolation coefficient for conductivity computation. */

    std::vector<MyTuple> Entr_Spline; /*!< \brief Spline interpolation coefficient for entropy computation. */

    std::vector<MyTuple> Cp_Spline; /*!< \brief Spline interpolation coefficient for specific heat at constant pressure computation. */

    std::vector<MyTuple> Enth_Spline; /*!< \brief Spline interpolation coefficient for enthalpy computation. */

    RealVec Enthalpies; /*!< \brief Enthalpy for each species. */

    RealVec CPs; /*!< \brief Specific heat at constant pressure for each species (Cp). */

    RealVec Formation_Enthalpies; /*!< \brief Formation enthalpy for each species. */

    RealMatrix Stoich_Coeffs_Reactants; /*!< \brief Stochiometric coefficents vector. */

    RealMatrix Stoich_Coeffs_Products; /*!< \brief Stochiometric coefficents vector. */

    RealMatrix Stoich_Coeffs_Reactants_Exp; /*!< \brief Stochiometric coefficents vector. */

    RealMatrix Stoich_Coeffs_Products_Exp; /*!< \brief Stochiometric coefficents vector. */

    std::vector<bool> Elementary_Reactions; /*!< \brief Vector to check if a reaction is elementary. */

    RealVec As;  /*!< \brief Vector with exponential pre factor. */

    std::vector<int> Ns;  /*!< \brief Vector with temperature exponent. */

    RealVec Temps_Activation;  /*!< \brief Vector with activation temperatures to estimate reaction rates for each reaction. */

    RealVec ys_over_mm; /*!< \brief Auxiliary vector to compute viscosity and thermal conductivity. */

    RealVec rhoUdiff; /*!< \brief Auxiliary vector for mass production term in case of constant Lewis number. */

    RealVec Dm_coeffs; /*!< \brief Auxiliary vector for effective diffusion coefficients. */

    RealVec omega; /*!< \brief Auxiliary vector for mass production term. */

    RealMatrix Dij; /*!< \brief Auxiliary matrix for diffusion binary coefficients. */

    RealMatrix Gamma; /*!< \brief Auxiliary matrix for Stefan-Maxwell equations. */

  private:

    enum {T_DATA_SPLINE = 0, X_DATA_SPLINE = 1, Y_DATA_SPLINE = 2}; /*!< \brief Enumerator for spline indexes. */

  }; /*-- End of class ReactingModelLibrary ---*/

} /*-- End of Namespace Framework ---*/

#endif
