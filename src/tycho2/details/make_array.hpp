template<class> struct is_ref_wrapper : std::false_type {};
template<class T> struct is_ref_wrapper<std::reference_wrapper<T>> : std::true_type {};
 
template<class T>
using not_ref_wrapper = std::negation<is_ref_wrapper<std::decay_t<T>>>;
 
template <class D, class...> struct return_type_helper { using type = D; };
template <class... Types>
struct return_type_helper<void, Types...> : std::common_type<Types...> {
  static_assert(std::conjunction_v<not_ref_wrapper<Types>...>,
		"Types cannot contain reference_wrappers when D is void");
};
 
template <class D, class... Types>
using return_type = std::array<typename return_type_helper<D, Types...>::type,
			       sizeof...(Types)>;
}
