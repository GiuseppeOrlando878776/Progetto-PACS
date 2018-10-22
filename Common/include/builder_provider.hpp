#ifndef SU2_BUILDER_PROVIDER
#define SU2_BUILDER_PROVIDER

#include <type_traits>
#include <memory>

#include "concrete_provider.hpp"

namespace Common {

/*!
  * \brief Abstract builder class. It forms the base for builders to be used with factories.
  * \author G. Orlando
  */
template<class Abstract,class Concrete>
class ProviderBuilder {
 public:
  static_assert(std::is_base_of<Abstract, Concrete>::value,
  		         "Builder requires that Abstract be a base of Concrete");

  /*!
   * \brief Default constructor
   */
  ProviderBuilder() = default;

  /*!
   * \brief Class destructor;
   */
  ~ProviderBuilder() {}

  /*!
   * \brief Create objects of type Concrete with zero arguments;
   */
  std::unique_ptr<Abstract> Create(typename Abstract::Arg1 arg) override {
    return std::unique_ptr<Abstract>(new Concrete(arg));
  }

  /*!
    *\brief Free an instance created by the factory
    *\@param ptr pointer to be freed
    */
  void FreeInstance(void* ptr) = 0;


}; /*-- End of class ProviderBuilder ---*/

template<class Abstract,class Concrete>
void ProviderBuilder<Abstract,Concrete>::FreeInstance(void* ptr)  {

    assert(ptr != NULL);
    auto obj = reinterpret_cast<ConcreteProvider<Abstract>*>(ptr);

    assert(obj != NULL);
    delete obj;
}



} /*-- End of namespace Common ---*/
#endif
