#ifndef SU2_ABSTRACT_FACTORY
#define SU2_ABSTRACT_FACTORY

#include "singleton.hpp"

#include <string>
#include <vector>

/*!
  * This namespace provides a factory class to load run-time the library that
  * will compute the physical and chemical properties of the considered mixture
*/

namespace Common {

  /*!
   * \class AbstractFactory
   * \brief Interface for factory to load libraries at run-time.
   * \author G. Orlando
  */
  class AbstractFactory: public Common::Singleton<AbstractFactory> {
  public:

    /*!
      * \brief Class constructor
    */
    AbstractFactory() {}

    /*!
      * \brief Virtual destructor
    */
    virtual ~AbstractFactory() {}

    /*
     * \brief Get the name of the Base class
    */
    virtual const std::string GetBaseName() const = 0;

    /*
     * \brief Get all the providers in this Factory in a std::vector
    */
    virtual std::vector<std::string> GetAllProviders() const = 0;

  };

}

#endif
