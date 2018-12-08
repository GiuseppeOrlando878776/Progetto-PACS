#ifndef SU2_PHYSICAL_PROPERTY_LIBRARY
#define SU2_PHYSICAL_PROPERTY_LIBRARY

#include "singleton.hpp"
#include "concrete_provider.hpp"

#include <string>

namespace Framework  {

  /*!
   * /brief Provides an abstract interface for libraries that compute the physical properties.
  */
  class PhysicalPropertyLibrary: public Common::Singleton<PhysicalPropertyLibrary>  {

  public:
    typedef Common::ConcreteProvider<PhysicalPropertyLibrary> Provider;
    typedef const std::string& Arg1;

  public:

    /*!
     * \brief Default constructor.
    */
    PhysicalPropertyLibrary() {}


    /*!
     * \brief Default destructor.
    */
    virtual ~PhysicalPropertyLibrary() = default;

    /*!
     *\brief Setup the path library name.
     *\param[in] lib_path_name - Path of the library
    */
    virtual void SetLibPathName(const std::string& lib_path_name) = 0;

    /*!
     * \brief Get the name of the library.
    */
    inline static std::string GetBaseName(void) {
        return "PhysicalPropertyLibrary";
    }

  }; /*-- End of class PhysicalPropertyLibrary ---*/

} /*-- End of namespace Framework ---*/

template class Common::Factory<Framework::PhysicalPropertyLibrary>;

#endif
