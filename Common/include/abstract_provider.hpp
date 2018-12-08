#ifndef SU2_ABSTRACT_PROVIDER
#define SU2_ABSTRACT_PROVIDER

#include "singleton.hpp"

#include <string>

namespace Common {

/*!
  * \brief Abstract class for provider types
  * \author G. Orlando
  */

  class AbstractProvider: public Common::Singleton<AbstractProvider> {

  public:

    /*!
      * \brief Class constructor
    */
    AbstractProvider() {}

    /*!
      * \brief Virtual destructor
    */
    virtual ~AbstractProvider() {}

    /*!
      * \brief Get the name of this provider
    */
    virtual std::string GetProviderName(void) const = 0;

  }; /*-- End of class AbstractProvider ---*/

} /*-- End of namespace Common ---*/

#endif
