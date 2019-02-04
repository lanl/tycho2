template<class...>
struct conjunction : std::true_type { };

template<class B1>
struct conjunction<B1> : B1 { };

template<class B1, class... Bn>
struct conjunction<B1, Bn...> 
    : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {};

template <bool B>
using bool_constant = std::integral_constant<bool, B>;

template<class B>
struct negation : bool_constant<!bool(B::value)> { };

template<class>
struct is_ref_wrapper : std::false_type {};

template<class T>
struct is_ref_wrapper<std::reference_wrapper<T>> : std::true_type {};
 
template<class T>
using not_ref_wrapper = negation<is_ref_wrapper<std::decay_t<T>>>;
 
template <class D, class...> struct return_type_helper { using type = D; };
template <class... Types>
struct return_type_helper<void, Types...> : std::common_type<Types...> {
  static_assert(conjunction<not_ref_wrapper<Types>...>::value,
		"Types cannot contain reference_wrappers when D is void");
};
 
template <class D, class... Types>
using return_type = std::array<typename return_type_helper<D, Types...>::type,
			       sizeof...(Types)>;
