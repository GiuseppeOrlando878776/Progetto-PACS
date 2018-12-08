#ifndef SU2_CONCRETE_PROVIDER
#define SU2_CONCRETE_PROVIDER

#include "provider.hpp"

namespace Common {

  /*!
    * \brief Concrete class for provider types
    * \author G. Orlando
  */

  template<class Base>
  class ConcreteProvider: public Common::Provider<Base> {

  public:

    /*!
      * \brief Constructor of the class
    */
    explicit ConcreteProvider(const std::string& name): Common::Provider<Base>(name) {}

    /*!
      * \brief Virtual destructor
    */
    ~ConcreteProvider() {}

    /*!
      *\brief Create an instance of provider: it must take exactly one argument
      *\param[in] arg - argument to construct the desired provider
    */
    virtual std::unique_ptr<Base> Create(typename Base::Arg1 arg) = 0;

  }; /*-- End of class ConcreteProvider ---*/


} /*-- End of namespace Common ---*/

#endif
