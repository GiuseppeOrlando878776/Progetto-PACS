#ifndef SU2_PHYSICAL_CHEMICAL_LIBRARY
#define SU2_PHYSICAL_CHEMICAL_LIBRARY

#include "physical_property_library.hpp"
#include "datatypes/vectorT.hpp"
#include "datatypes/matrixT.hpp"

namespace Framework  {

  /*!
    * /brief Provides an abstract interface for libraries that compute the physical properties.
    */
    class PhysicalChemicalLibrary: public Framework::PhysicalPropertyLibrary {

    public:
      typedef Common::ConcreteProvider<PhysicalChemicalLibrary> Provider;
      typedef const std::string& Arg1;

      using RealVec = Common::RealVec;
      using RealMatrix = Common::RealMatrix;

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
       inline void SetNReactions(const unsigned nr) {
         nReactions = nr;
       }

       /*!
         * Set the constant of gases for each species [J/(Kg*K)]
         */
       virtual void SetRiGas(void) = 0;

       /*!
         * \brief Get the constant of perfect gases [J/(Kg*K)]
         */
       virtual su2double GetRgas(void) const = 0;

       /*!
         * \brief Set the constant of perfect gases for the mixture [J/(Kg*K)]
         * \param[in] ys The vector of the mass fractions of species (input)
         */
       virtual void SetRgas(const RealVec& ys) = 0;

       /*!
         * \brief Compute the constant of perfect gases for the mixture [J/(Kg*K)]
         * \param[in] ys The vector of the mass fractions of species (input)
         */
       virtual su2double ComputeRgas(const RealVec& ys) = 0;

       /*!
         * \brief Get the molar masses of the species
         */
       virtual void GetMolarMasses(RealVec& mm) const = 0;

      /*!
        * \brief Gets the molar fractions.
        * \param[in] ys The vector of the mass fractions of species (input)
        * \param[out] xs The vector of the molar fractions of species (output)
        */
       virtual void GetMolarFractions(const RealVec& ys, RealVec& xs) = 0;

       /*!
        * \brief Sets the molar fractions of elements Xn.This function should be called before getting
        * \brief thermodynamic quantities or transport properties.
        * \param[in] xn - The vector of the mass fractions of elements
        */
       virtual void SetMolarFractions(const RealVec& xs) = 0;

      /*!
        * \brief Gets the mass fractions.
        * \param[in] xs The vector of the molar fractions of species (input)
        * \param[out] ys The vector of the mass fractions of species (output)
        */
       virtual void GetMassFractions(const RealVec& xs, RealVec& ys) = 0;

       /*!
        * \brief Sets the mass fractions. This function should be called before getting
        * \brief thermodynamic quantities or transport properties.
        * \param[in] ys The vector of the mass fractions of species
       */
       virtual void SetMassFractions(const RealVec& ys) = 0;

      /*!
        * \brief Sets the mole fractions of elements Xn starting from the given species mass fractions Yn
        * \param[in] ys - The vector of the mass fractions of species
       */
      // virtual void SetMolarFromMass(const RealVec& ys) = 0;

       /*!
        * \brief Computes the specific heat ratio and the speed of sound.
        * \param[in] temp - temperature
        * \param[out] gamma - specific heat ratio
        * \param[out] sound_speed - speed of sound
        */
      virtual void Gamma_FrozenSoundSpeed(const su2double temp, su2double& gamma, su2double& sound_speed) = 0;

       /*!
        * \brief Computes the density, the enthalpy and the internal energy
        * \param[in] temp - temperature
        * \param[in] pressure - pressure
        * \param[in] ys - mass fractions
        * \param[out] dhe - Vector with density, enthalpy, energy (output) for thermal equilibrium
        */
      virtual void Density_Enthalpy_Energy(const su2double temp, const su2double pressure, const RealVec& ys, RealVec& dhe) = 0;

       /*!
        * \brief Computes the density given temperature and pressure.
        * \param[in] temp - temperature
        * \param[in] pressure - pressure
        * \param[in] ys - mass fractions
        */
      virtual su2double ComputeDensity(const su2double temp, const su2double pressure, const RealVec& ys) = 0;

       /*!
        * \brief Computes the internal energy per unit of mass at given temperature and pressure.
        * \param[in] temp - temperature
        * \param[in] pressure - pressure
        * \param[in] ys - mass fractions
        */
      virtual su2double ComputeEnergy(const su2double temp, const su2double pressure, const RealVec& ys) = 0;

       /*!
        * \brief Returns the formation enthalpies per unit mass of species
        * \param[out] hs - species formation enthalpies (output)
        */
      virtual void GetFormationEnthalpies(RealVec& hs) = 0;

       /*!
        * \brief Returns the static enthalpy per unit of mass
        * \param[in] temp - the mixture temperature
        * \param[in] ys - mass fractions
        * \return Mixture static enthalpy (output)
        */
      virtual su2double ComputeEnthalpy(const su2double temp, const RealVec& ys) = 0;

       /*!
        * \brief Computes the specific heat at constant pressure
        * \param[in] temp - temperature
        * \param[in] ys - mass fractions
        * \return Cp - specific heat at constant pressure
        */
      virtual su2double ComputeCP(const su2double temp, const RealVec& ys) = 0;

      /*!
        * \brief Computes the thermal conductivity given temperature
        * \param[in] temp - temperature
        * \param[in] ys - mass fractions
        * \return Total thermal conductivity for thermal non equilibrium
        */
      virtual su2double ComputeLambda(const su2double temp, const RealVec& ys) = 0;

       /*!
        * \brief Computes the dynamic viscosity given temperature and
        * \param[in] temp - temperature
        * \param[in] ys - mass fractions
        * \return Total molecular viscosity for thermal non equilibrium
        */
      virtual su2double ComputeEta(const su2double temp, const RealVec& ys) = 0;

       /*!
        * Returns the mass production/destruction terms [kg m^-3 s^-1] in chemical
        * brief non-equilibrium based on Arrhenius's formula.
        * \param[in] pressure - the mixture pressure
        * \param[in] temp - the mixture temperature
        * \param[in] ys - the species mass fractions
        * \param[in] omega - the mass production terms
        * \param[in] jacobian - the Jacobian matrix of the mass production terms
       */
      virtual void GetMassProductionTerm(const su2double temp, const su2double pressure, const su2double rho,
                                         RealVec& ys, RealVec& omega, RealMatrix& jacobian) = 0;

       /*!
        * \brief Returns the source terms species continuity and energy equations
        * \param[in] temp - the mixture temperature
        * \param[in] pressure - the mixture pressure
        * \param[in] ys - the species mass fractions
        * \param[in] rho - the mixture density
        * \param[in] omega - the mass production term
        * \param[in] omegav - the source term
       */
       //virtual void GetSourceTerm(const su2double temp, const su2double pressure, const su2double rho,
      //                            RealVec& ys, RealVec& omega, RealVec& omegav) = 0;

      /*!
       * Returns the diffusion velocities of species multiplied by the species
       * densities for nonequilibrium computations
       * \param[in] temp - the mixture temperature
       * \param[in] rho - the mixture density
       * \return rhoUdiff  - Diffusion coefficient with constant Lewis number for each species
       */
     virtual RealVec GetRhoUdiff(const su2double temp, const su2double rho) = 0;

     /*!
      * Returns the diffusion velocities of species multiplied by the species
      * densities for nonequilibrium computations
      * \param[in] temp the mixture temperature
      * \param[in] pressure the mixture pressure
      */
      virtual void GetDij_fick(RealVec& dx, const su2double pressure, const su2double temperature, RealMatrix& Dij) = 0;

    protected:

       unsigned short nSpecies; /*!< \brief Number of species. */

       unsigned nReactions; /*!< \brief Number of reactions. */

     public:

       static constexpr su2double NA = 6.02214129*1e23; /*!< \brief Avogadro number. */

       static constexpr su2double KB = 1.3806488*1e-23; /*!< \brief Boltzmann constant. */

       static constexpr su2double R_ungas = 8.31446215; /*!< \brief Universal constant of perfect gas. */

  }; /*-- End of class PhysicalChemicalLibrary ---*/

} /*-- End of namespace Framework ---*/

#endif
