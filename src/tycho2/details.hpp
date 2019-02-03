namespace details{

#include "details/make_array.hpp" 

} // namespace details

template < class D = void, class... Types>
constexpr details::return_type<D, Types...> make_array(Types&&... t) {
  return {std::forward<Types>(t)... };
}



