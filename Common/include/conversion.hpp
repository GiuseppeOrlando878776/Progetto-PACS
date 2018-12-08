#ifndef SU2_CONVERSION_POINTER
#define SU2_CONVERSION_POINTER

#include <memory>
#include "su2_assert.hpp"

namespace Common {
  template <class Library,class Interface>
  Library* GetLibrary(Interface* interface_ptr) {
    SU2_Assert(interface_ptr != NULL, "The interface for the library has not been instancieted");
    return dynamic_cast<Library*>(interface_ptr);
  }
}

#endif
