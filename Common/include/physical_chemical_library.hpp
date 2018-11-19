#ifndef SU2_PHYSICAL_CHEMICAL_LIBRARY
#define SU2_PHYSICAL_CHEMICAL_LIBRARY

#include "physical_property_library.hpp"
#include "datatypes/matrixT.hpp"

namespace Framework  {

    typedef std::vector<double> RealVec;
    using RealMatrix = Common::RealMatrix;

    /*!
      * /brief Provides an abstract interface for libraries that compute the physical properties.
    */
    class PhysicalChemicalLibrary: public Framework::PhysicalPropertyLibrary {

    public:
      typedef Common::ConcreteProvider<PhysicalChemicalLibrary> Provider;
      typedef const std::string& Arg1;

    public:

      /*!
       * \brief Default constructor.
       */
       PhysicalChemicalLibrary():PhysicalPropertyLibrary(),nSpecies(),nReactions() {}

      /*!
       * \brief Class constructor.
       */
       explicit PhysicalChemicalLibrary(const std::string& name);

      /*!
       * \brief Default destructor.
       */
       virtual ~PhysicalChemicalLibrary() = default;

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
         * Set the constant of gases for each species [J/(Kg*K)]
         */
       virtual void SetRiGas(void) = 0;

       /*!
         * \brief Get the constant of perfect gases [J/(Kg*K)]
         */
       virtual double GetRgas(void) const = 0;

       /*!
         * \brief Set the constant of perfect gases for the mixture [J/(Kg*K)]
         * \param[in] ys The vector of the mass fractions of species (input)
         */
       virtual void SetRgas(const RealVec& ys) = 0;

       /*!
         * \brief Compute the constant of perfect gases for the mixture [J/(Kg*K)]
         * \param[in] ys The vector of the mass fractions of species (input)
         */
       virtual double ComputeRgas(const RealVec& ys) = 0;

       /*!
         * \brief Get the molar masses of the species
         */
       virtual RealVec GetMolarMasses(void) const = 0;

      /*!
        * \brief Gets the molar fractions from mass fractions.
        * \param[in] ys The vector of the mass fractions of species (input)
        */
       virtual void GetMolarFractions(const RealVec& ys) = 0;

       /*!
        * Sets the molar fractions of elements Xn.This function should be called before getting
        * thermodynamic quantities or transport properties.
        * \param[in] xn - The vector of the mass fractions of elements
        */
       virtual void SetMolarFractions(const RealVec& xs) = 0;

      /*!
        * \brief Gets the mass fractions from molar fractions.
        * \param[in] xs The vector of the molar fractions of species (input)
        */
       virtual void GetMassFractions(const RealVec& xs) = 0;

       /*!
        * Sets the mass fractions. This function should be called before getting
        * thermodynamic quantities or transport properties.
        * \param[in] ys The vector of the mass fractions of species
       */
       virtual void SetMassFractions(const RealVec& ys) = 0;

       /*!
        * \brief Computes the specific heat ratio and the speed of sound.
        * \param[in] temp - temperature
        * \param[in] ys - The vector of the mass fractions of species
        * \param[out] gamma - specific heat ratio
        * \param[out] sound_speed - speed of sound
        */
       virtual void Gamma_FrozenSoundSpeed(const double temp, const RealVec& ys, double& gamma, double& sound_speed) = 0;

       /*!
        * \brief Computes the density, the enthalpy and the internal energy
        * \param[in] temp - temperature
        * \param[in] pressure - pressure
        * \param[in] ys - mass fractions
        * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
        */
       virtual void Density_Enthalpy_Energy(const double temp, const double pressure, const RealVec& ys, RealVec& dhe) = 0;

       /*!
        * \brief Computes the density given temperature and pressure.
        * \param[in] temp - temperature
        * \param[in] rho - density
        * \param[in] ys - mass fractions
        */
       virtual double ComputePressure(const double temp, const double rho, const RealVec& ys) = 0;

       /*!
        * \brief Computes the density given temperature and pressure.
        * \param[in] temp - temperature
        * \param[in] pressure - pressure
        * \param[in] ys - mass fractions
        */
       virtual double ComputeDensity(const double temp, const double pressure, const RealVec& ys) = 0;

       /*!
        * \brief Computes the internal energy per unit of mass at given temperature and pressure.
        * \param[in] temp - temperature
        * \param[in] pressure - pressure
        * \param[in] ys - mass fractions
        */
       virtual double ComputeEnergy(const double temp, const double pressure, const RealVec& ys) = 0;

       /*!
        * \brief Returns the formation enthalpies per unit mass of species
        */
       virtual RealVec GetFormationEnthalpies(void) const = 0;

       /*!
        * \brief Returns the static enthalpy per unit of mass
        * \param[in] temp - the mixture temperature
        * \param[in] ys - mass fractions
        * \return Mixture static enthalpy (output)
        */
       virtual double ComputeEnthalpy(const double temp, const RealVec& ys) = 0;

       /*!
        * \brief Returns the static enthalpy per unit of mass of each species
        * \param[in] temp - the mixture temperature
        */
       virtual RealVec ComputePartialEnthalpy(const double temp) = 0;

       /*!
        * \brief Computes the specific heat at constant pressure
        * \param[in] temp - temperature
        * \param[in] ys - mass fractions
        * \return Cp - specific heat at constant pressure
        */
       virtual double ComputeCP(const double temp, const RealVec& ys) = 0;

       /*!
        * \brief Computes the specific heat at constant volume
        * \param[in] temp - temperature
        * \param[in] ys - mass fractions
        * \return Cv - specific heat at constant volume
        */
       virtual double ComputeCV(const double temp, const RealVec& ys) = 0;

       /*!
        * \brief Computes the mixture total concentration
        * \param[in] rho - density
        */
       virtual double ComputeConcentration(const double rho) = 0;

       /*!
         * \brief Computes the thermal conductivity given temperature
         * \param[in] temp - temperature
         * \param[in] ys - mass fractions
         * \return Total thermal conductivity for thermal non equilibrium
         */
       virtual double ComputeLambda(const double temp, const RealVec& ys) = 0;

       /*!
        * \brief Computes the dynamic viscosity given temperature and
        * \param[in] temp - temperature
        * \param[in] ys - mass fractions
        * \return Total molecular viscosity for thermal non equilibrium
        */
       virtual double ComputeEta(const double temp, const RealVec& ys) = 0;

       /*!
        * Returns the mass production/destruction terms [kg m^-3 s^-1] in chemical
        * non-equilibrium based on Arrhenius's formula.
        * \param[in] temp - the mixture temperature
        * \param[in] rho - the mixture density
        * \param[in] ys - the species mass fractions
        * \return Mass production terms
        */
       virtual RealVec GetMassProductionTerm(const double temp, const double rho, const RealVec& ys) = 0;

      /*!
       * Returns the diffusion velocities of species multiplied by the species
       * densities for nonequilibrium computations
       * \param[in] temp - the mixture temperature
       * \param[in] rho - the mixture density
       * \param[in] ys - the species mass fractions
       * \return Diffusion coefficient with constant Lewis number for each species
       */
      virtual RealVec GetRhoUdiff(const double temp, const double rho, const RealVec& ys) = 0;

      /*!
       * \brief Returns the binary diffusion coefficients for Stefan-Maxwell equation
       * \param[in] temp the mixture temperature
       * \param[in] pressure the mixture pressure
       */
      virtual RealMatrix GetDij_SM(const double temp, const double pressure) = 0;

      /*!
       * \brief Returns the effective diffusion coefficients to solve Stefan-Maxwell equation
       * \param[in] temp - the mixture temperature
       * \param[in] pressure - the mixture pressure
       * \param[in] ys - mass fractions in the mixture
       */
       virtual RealVec GetDiffCoeffs(const double temp, const double pressure, const RealVec& ys) = 0;

    protected:

       unsigned short nSpecies; /*!< \brief Number of species. */

       unsigned short nReactions; /*!< \brief Number of reactions. */

     public:

       static constexpr double NA = 6.02214129*1e23; /*!< \brief Avogadro number. */

       static constexpr double KB = 1.3806488*1e-23; /*!< \brief Boltzmann constant. */

       static constexpr double R_ungas = 8.31446215; /*!< \brief Universal constant of perfect gas. */

       static constexpr double R_ungas_scal = 1.9858775; /*!< \brief Universal gas constant in cal/(mol K). */

  }; /*-- End of class PhysicalChemicalLibrary ---*/

} /*-- End of namespace Framework ---*/

#endif
