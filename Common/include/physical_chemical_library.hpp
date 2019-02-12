#ifndef SU2_PHYSICAL_CHEMICAL_LIBRARY
#define SU2_PHYSICAL_CHEMICAL_LIBRARY

#include "concrete_provider.hpp"

namespace Framework  {

  typedef std::vector<double> RealVec;

  /*!
   * /brief Provides an abstract interface for libraries that compute the physical properties.
  */
  template<class Vector,class Matrix>
  class PhysicalChemicalLibrary: public Common::NotCopyable<PhysicalChemicalLibrary<Vector,Matrix>> {

  public:
    typedef Common::ConcreteProvider<PhysicalChemicalLibrary> Provider;
    typedef const std::string& Arg1;

  public:

    /*!
     * \brief Class constructor.
     */
    PhysicalChemicalLibrary(const std::string& name): Config_File(name), Lib_Setup(false), nSpecies(), nReactions() {};

    /*!
     * \brief Default destructor.
     */
    virtual ~PhysicalChemicalLibrary() = default;

    /*!
     * \brief Get the name of the base class of library.
    */
    inline static std::string GetBaseName(void) {
      return "PhysicalChemicalLibrary";
    }

    /*!
     *\brief Setup the path library name.
     *\param[in] lib_path_name - Path of the library
    */
    inline void SetLibPathName(const std::string& lib_path_name) {
      Lib_Path = lib_path_name;
    }

    /*!
     * Get the number of species in the mixture
     */
    inline const unsigned short GetNSpecies(void) const {
      return nSpecies;
    }

    /*!
     * Set the number of species in the mixture
    */
    inline void SetNSpecies(const unsigned short ns) {
      nSpecies = ns;
    }

    /*!
     * Get the number of reactions in the mixture
    */
    inline const unsigned GetNReactions(void) const {
      return nReactions;
    }

    /*!
     * Set the number of reactions in the mixture
    */
    inline void SetNReactions(const unsigned short nr) {
      nReactions = nr;
    }

    /*!
     * Set the constant of gases for each species [J/(Kmol*K)]
    */
    virtual void SetRiGas(void) = 0;

    /*!
     * \brief Get the constant of perfect gases for the mixture [J/(Kmol*K)]
    */
    virtual double GetRgas(void) const = 0;

    /*!
     * \brief Set the constant of perfect gases for the mixture [J/(Kmol*K)]
     * \param[in] ys - The Vectortor of the mass fractions of species (input)
    */
    virtual void SetRgas(const Vector& ys) = 0;

    /*!
     * \brief Compute the constant of perfect gases for the mixture [J/(Kmol*K)]
     * \param[in] ys - The Vectortor of the mass fractions of species (input)
    */
    virtual double ComputeRgas(const Vector& ys) = 0;

    /*!
     * \brief Get the molar masses of the species
    */
    virtual Vector GetMolarMasses(void) const = 0;

    /*!
     * Set the molar fractions of elements Xn.This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] xs - The Vectortor of the mass fractions of elements
    */
    virtual void SetMolarFractions(const Vector& xs) = 0;

    /*!
     * Set the mass fractions. This function should be called before getting
     * thermodynamic quantities or transport properties.
     * \param[in] ys - The Vectortor of the mass fractions of species
    */
    virtual void SetMassFractions(const Vector& ys) = 0;

    /*!
     * \brief Set the molar fractions from mass fractions.
     * \param[in] ys - The Vectortor of the mass fractions of species (input)
    */
    virtual void SetMolarFromMass(const Vector& ys) = 0;

    /*!
     * \brief Get the molar fractions from mass fractions.
     * \param[in] ys - The Vectortor of the mass fractions of species (input)
    */
    virtual Vector GetMolarFromMass(const Vector& ys) = 0;

    /*!
     * \brief Get the mass fractions from molar fractions.
     * \param[in] xs - The Vectortor of the molar fractions of species (input)
    */
    virtual void SetMassFromMolar(const Vector& xs) = 0;

    /*!
     * \brief Get the mass fractions from molar fractions.
     * \param[in] xs - The Vectortor of the molar fractions of species (input)
    */
    virtual Vector GetMassFromMolar(const Vector& xs) = 0;

    /*!
     * \brief Compute the specific heat ratio and the speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The Vectortor of the mass fractions of species
     * \param[out] gamma - specific heat ratio
     * \param[out] sound_speed - speed of sound
    */
    virtual void Gamma_FrozenSoundSpeed(const double temp, const Vector& ys, double& gamma, double& sound_speed) = 0;

    /*!
     * \brief Computes the frozen specific heat ratio.
     * \param[in] temp - temperature
     * \param[in] ys - The Vectortor of the mass fractions of species
     */
    virtual double ComputeFrozenGamma(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The Vectortor of the mass fractions of species
    */
    virtual double ComputeFrozenSoundSpeed(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Computes the frozen speed of sound.
     * \param[in] temp - temperature
     * \param[in] ys - The Vectortor of the mass fractions of species
     * \param[in] press - pressure
     * \param[in] rho - density
    */
    virtual double ComputeFrozenSoundSpeed(const double temp, const Vector& ys, const double press, const double rho) = 0;

    /*!
     * \brief Compute the density, the enthalpy and the internal energy
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
     * \param[out] dhe - Vectortor with density, enthalpy, energy (output) for thermal equilibrium
    */
    virtual void Density_Enthalpy_Energy(const double temp, const double pressure, const Vector& ys, RealVec& dhe) = 0;

    /*!
     * \brief Compute the density given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    virtual double ComputePressure(const double temp, const double rho, const Vector& ys) = 0;

    /*!
     * \brief Compute the density given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] pressure - pressure
     * \param[in] ys - mass fractions
    */
    virtual double ComputeDensity(const double temp, const double pressure, const Vector& ys) = 0;

    /*!
     * \brief Compute the internal energy per unit of mass at given temperature and pressure.
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     */
    virtual double ComputeEnergy(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Returns the formation enthalpies per unit mass of species
    */
    virtual Vector GetFormationEnthalpies(void) const = 0;

    /*!
     * \brief Returns the static enthalpy per unit of mass
     * \param[in] temp - the mixture temperature
     * \param[in] ys - mass fractions
     * \return Mixture static enthalpy (output)
    */
    virtual double ComputeEnthalpy(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Returns the static enthalpy per unit of mass of each species
     * \param[in] temp - the mixture temperature
    */
    virtual Vector ComputePartialEnthalpy(const double temp) = 0;

    /*!
     * \brief Computes the mixture total concentration
     * \param[in] rho - density
     * \param[in] ys - mass fractions
    */
    virtual double ComputeConcentration(const double rho, const Vector& ys) = 0;

    /*!
     * \brief Compute the specific heat at constant pressure
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Cp - specific heat at constant pressure
     */
    virtual double ComputeCP(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Compute the specific heat at constant volume
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Cv - specific heat at constant volume
    */
    virtual double ComputeCV(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Compute the thermal conductivity given temperature
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Total thermal conductivity for thermal non equilibrium
     */
    virtual double ComputeLambda(const double temp, const Vector& ys) = 0;

    /*!
     * \brief Compute the dynamic viscosity given temperature and
     * \param[in] temp - temperature
     * \param[in] ys - mass fractions
     * \return Total molecular viscosity for thermal non equilibrium
    */
    virtual double ComputeEta(const double temp, const Vector& ys) = 0;

    /*!
     * Returns the mass production/destruction terms [kg m^-3 s^-1] in chemical
     * non-equilibrium based on Arrhenius's formula.
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \param[in] ys - the species mass fractions
     * \return Mass production terms
    */
    virtual Vector GetMassProductionTerm(const double temp, const double rho, const Vector& ys) = 0;

    /*!
     * Returns the diffusion velocities of species multiplied by the species
     * densities for nonequilibrium computations
     * \param[in] temp - the mixture temperature
     * \param[in] rho - the mixture density
     * \param[in] ys - the species mass fractions
     * \return Diffusion coefficient with constant Lewis number for each species
    */
    virtual Vector GetRhoUdiff(const double temp, const double rho, const Vector& ys) = 0;

    /*!
     * Returns the diffusion velocities of species multiplied by the species
     * densities for nonequilibrium computations
     * \param[in] pressure - the mixture pressure
     * \param[in] temp - the mixture temperature
    */
    virtual Matrix GetDij_SM(const double pressure, const double temp) = 0;

    /*!
     * Returns thematrix of Stefan-Maxwell equations
     * \param[in] rho - the mixture density
     * \param[in] xs - current molar fractions
     * \param[in] ys - current mass fractions
     * \param[in] val_Dij - current binary diffusion coefficients
     */
    virtual Matrix GetGamma(const double rho, const Vector& xs, const Vector& ys, const Matrix& val_Dij) = 0;

    /*!
     * \brief Returns the effective diffusion coefficients to solve Stefan-Maxwell equation using Sutton algorithm
     * \param[in] temp - the mixture temperature
     * \param[in] pressure - the mixture pressure
     * \param[in] ys - mass fractions in the mixture
     */
    virtual Vector GetDiffCoeffs(const double temp, const double pressure, const Vector& ys) = 0;

  protected:

    std::string Config_File; /*!\brief File with all the configuration options */

    std::string Lib_Path; /*!\brief Path to some input library data */

    bool Lib_Setup; /*!\brief Bool to check if library is setup */

    unsigned short nSpecies; /*!< \brief Number of species. */

    unsigned short nReactions; /*!< \brief Number of reactions. */

  public:

    static constexpr double NA = 6.02214129*1e23; /*!< \brief Avogadro number. */

    static constexpr double KB = 1.3806488*1e-23; /*!< \brief Boltzmann constant. */

    static constexpr double R_ungas = NA*KB*1e3; /*!< \brief Universal constant of perfect gas. (J/kmol*K) */

    static constexpr double R_ungas_scal = 1.9858775; /*!< \brief Universal gas constant in cal/(mol K). */

  }; /*-- End of class PhysicalChemicalLibrary ---*/

} /*-- End of namespace Framework ---*/

namespace Common {
  template<class Vector,class Matrix>
  class Common::Factory<Framework::PhysicalChemicalLibrary<Vector,Matrix>>;
}

#endif
