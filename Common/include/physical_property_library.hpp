#ifndef SU2_PHYSICAL_PROPERTY_LIBRARY
#define SU2_PHYSICAL_PROPERTY_LIBRARY

#include "builder_provider.hpp"
#include "datatype_structure.hpp"

namespace Framework  {

  /*!
    * /brief Provides an abstract interface for libraries that compute the physical properties.
    */
    class PhysicalPropertyLibrary {

    public:

      /*!
       * \brief Default constructor.
       */
       PhysicalPropertyLibrary():Lib_Setup(false) {}

       /*!
        * \brief Class constructor.
        */
       explicit PhysicalPropertyLibrary(const std::string& name);


      /*!
       * \brief Default destructor.
       */
       virtual ~PhysicalPropertyLibrary() = default;

       /*!
        *\brief Setups the path library name.
       */
       inline void SetLibPathName(const std::string& lib_path_name)  {
         Lib_Path = lib_path_name;
       }

       /*!
       * \brief Get the name of the library.
       */
       inline static std::string GetBaseName(void) {
          return "PhysicalPropertyLibrary";
       }

    protected:

       std::string Lib_Path; /*!\brief Path to some input library data */

       bool Lib_Setup; /*!\brief Setup library */

  }; /*-- End of class PhysicalPropertyLibrary ---*/

} /*-- End of namespace Framework ---*/

#endif
