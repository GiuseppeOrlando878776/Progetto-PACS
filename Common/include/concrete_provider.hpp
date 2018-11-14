#ifndef SU2_CONCRETE_PROVIDER
#define SU2_CONCRETE_PROVIDER

#include "abstract_provider.hpp"

#include <string>
#include <cassert>

namespace Common {

/*!
  * \brief Concrete class for provider types
  * \author G. Orlando
  */

template<class Base>
class ConcreteProvider: public Common::AbstractProvider<Base> {

public:

  /*!
    * \brief Constructor of the class
    */
  explicit ConcreteProvider(const std::string& name): Common::AbstractProvider<Base>(name) {}

  /*!
    * \brief Virtual destructor
    */
  ~ConcreteProvider() {}

  /*!
    *\brief Free an instance created by the factory
    *\@param ptr pointer to be freed
    */
  virtual std::unique_ptr<Base>  Create(typename Base::Arg1 arg) = 0;


}; /*-- End of class ConcreteProvider ---*/


} /*-- End of namespace Common ---*/

#endif
