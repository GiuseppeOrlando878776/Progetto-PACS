#ifndef SU2_DEFAULT_INITIALIZATION
#define SU2_DEFAULT_INITIALIZATION

#include <cstdlib>
#include <utility>
#include <tuple>

namespace Common {

  template<typename T,std::size_t... Is>
  std::tuple<T> repeat_impl(const T& x) {
    return std::make_tuple((Is,x),...);
  }

  template<std::size_t N,typename T>
  auto repeat(const T& x) {
    return repeat_impl(x,std::make_index_sequence<N>());
  }

}

#endif
