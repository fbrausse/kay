/*
 * kay/bits.hh
 *
 * Copyright 2017-2019 Franz Brauße <brausse@informatik.uni-trier.de>
 *
 * This file is part of kay.
 * See the LICENSE file for terms of distribution.
 */

#ifndef KAY_BITS_HH
#define KAY_BITS_HH

#include <cstdint>
#include <climits>	/* CHAR_BIT */
#include <cstddef>	/* size_t */
#include <type_traits>	/* std::integral_constant */
#include <variant>	/* std::monostate */
#include <vector>	/* std::vector */

namespace kay {

/* Fowler/Noll/Vo hash */
template <typename T> struct fnv1_hash_params;

/* from <http://www.isthe.com/chongo/tech/comp/fnv/index.html#FNV-param> */
template <> struct fnv1_hash_params<uint32_t> {
	constexpr static const uint32_t offset_basis = 2166136261U;
	constexpr static const uint32_t fnv_prime = 16777619U;
};

template <> struct fnv1_hash_params<uint64_t> {
	constexpr static const uint64_t offset_basis = 14695981039346656037U;
	constexpr static const uint64_t fnv_prime = 1099511628211U;
};

template <typename T> struct fnv1_hash : fnv1_hash_params<T> {

	constexpr static T combine(T hsh, const T &v) noexcept
	{
		hsh *= fnv1_hash_params<T>::fnv_prime;
		hsh ^= v;
		return hsh;
	}
};

template <typename T> struct fnv1a_hash : fnv1_hash_params<T> {

	constexpr static T combine(T hsh, const T &v) noexcept
	{
		hsh ^= v;
		hsh *= fnv1_hash_params<T>::fnv_prime;
		return hsh;
	}
};

template <typename T, typename H>
struct hash
{
	static_assert(H::template is_std<T>, "specialization of kay::hash required");

	constexpr size_t operator()(const T &x) const noexcept
	{
		return H(x).v;
	}
};

template <template <typename> typename Std>
struct hash_base : private fnv1a_hash<size_t> {

	template <typename T> static const constexpr bool is_std = Std<T>::value;

	size_t v;
	template <typename T>
	hash_base(const T &t)
	{
		if constexpr (Std<T>::value) {
			v = std::hash<T>{}(t);
		} else
			v = hash<T,hash_base<Std>>{}(t);
	}
	hash_base() : v(offset_basis) {}
	friend hash_base operator^(const hash_base &a, const hash_base &b)
	{
		return fnv1a_hash<size_t>::combine(a.v, b.v);
	}
	friend hash_base & operator^=(hash_base &a, const hash_base &b)
	{
		return a = a ^ b;
	}
};

template <typename... Ts, typename H> struct hash<std::tuple<Ts...>,H> {

	template <size_t... Is>
	constexpr size_t do_hash(const std::tuple<Ts...> &t,
	                         std::index_sequence<Is...>) const noexcept
	{
		return (H() ^ ... ^ H(std::get<Is>(t))).v;
	}

	constexpr size_t operator()(const std::tuple<Ts...> &t) const noexcept
	{
		return do_hash(t, std::make_index_sequence<sizeof...(Ts)>{});
	}
};

template <typename T, typename H> struct hash<std::vector<T>,H> {

	size_t operator()(const std::vector<T> &v) const noexcept
	{
		H h(size(v));
		for (const T &x : v)
			h ^= H(x);
		return h.v;
	}
};

template <typename T> struct hash_std_default
: std::disjunction<std::is_integral<std::remove_cv_t<T>>, std::is_pointer<T>>
{};

template <typename T, template <typename> typename Std = hash_std_default>
static inline size_t do_hash(const T &v) { return hash<T,hash_base<Std>>{}(v); }



template <size_t n> struct ceil_log2;
template <size_t n> constexpr size_t ceil_log2_v = ceil_log2<n>::value;
template <> struct ceil_log2<1> : std::integral_constant<size_t,0> {};
template <size_t n> struct ceil_log2 : std::integral_constant<size_t,1+ceil_log2_v<(n+1)/2>> {};

/* total bits used */
template <typename T> struct type_bits : std::integral_constant<size_t,CHAR_BIT*sizeof(T)> {};
template <typename T> constexpr size_t type_bits_v = type_bits<T>::value;

template <typename T> struct cardinality;
template <typename T> constexpr size_t cardinality_v = cardinality<T>::value;

/* max #lower bits used */
template <typename T> struct max_bits : ceil_log2<cardinality_v<T>> {};
template <typename T> constexpr size_t max_bits_v = max_bits<T>::value;

using unit = std::monostate;

template <typename T> struct cardinality
: std::enable_if_t<(type_bits_v<size_t> >= max_bits_v<T>)
                  ,std::integral_constant<size_t,((size_t)1 << max_bits_v<T>)>> {};

/* cardinalities for special types that either don't use all their bits or have
 * no bits at all */
template <> struct cardinality<void> : std::integral_constant<size_t,0> {};
template <> struct cardinality<unit> : std::integral_constant<size_t,1> {};
template <> struct cardinality<bool> : std::integral_constant<size_t,2> {};

static_assert(max_bits_v<unit> == 0);
static_assert(max_bits_v<bool> == 1);


template <size_t n> struct integral_at_least_log_bytes;
template <> struct integral_at_least_log_bytes<0> { using type = uint8_t; };
template <> struct integral_at_least_log_bytes<1> { using type = uint16_t; };
template <> struct integral_at_least_log_bytes<2> { using type = uint32_t; };
template <> struct integral_at_least_log_bytes<3> { using type = uint64_t; };
template <> struct integral_at_least_log_bytes<4> { using type = unsigned __int128_t; };

template <size_t n>
struct integral_at_least_bits
: integral_at_least_log_bytes<ceil_log2_v<(n+7)/8>> {};

template <size_t n>
using integral_at_least_bits_t = typename integral_at_least_bits<n>::type;

template <> struct integral_at_least_bits<0> { using type = unit; };
template <> struct integral_at_least_bits<1> { using type = bool; };

template <size_t tag_sz, typename I, typename T,
          typename = std::enable_if<(type_bits_v<I> > tag_sz)>>
struct tagged_idx_base {

	static constexpr size_t idx_size = type_bits_v<I>;
	static constexpr size_t tag_size = tag_sz;

	union {
		I v;
		struct {
			I idx : type_bits_v<I> - tag_sz;
			T cat : tag_sz;
		};
	};

	tagged_idx_base() {};
	tagged_idx_base(I idx, T cat) : idx(idx), cat(cat) {}

	constexpr friend bool operator==(const tagged_idx_base &a,
	                                 const tagged_idx_base &b)
	{ return a.v == b.v; }

	constexpr friend bool operator!=(const tagged_idx_base &a,
	                                 const tagged_idx_base &b)
	{ return a.v != b.v; }

	constexpr friend bool operator< (const tagged_idx_base &a,
	                                 const tagged_idx_base &b)
	{ return a.v <  b.v; }
};

template <typename I, typename T>
struct tagged_idx_base<0,I,T> {

	union {
		I v;
		struct { I idx; };
	};

	constexpr friend bool operator==(const tagged_idx_base &a,
	                                 const tagged_idx_base &b)
	{ return a.v == b.v; }

	constexpr friend bool operator!=(const tagged_idx_base &a,
	                                 const tagged_idx_base &b)
	{ return a.v != b.v; }

	constexpr friend bool operator< (const tagged_idx_base &a,
	                                 const tagged_idx_base &b)
	{ return a.v <  b.v; }
};

template <typename T, typename I = uint32_t>
using tagged_idx = std::enable_if_t<sizeof(tagged_idx_base<max_bits_v<T>, I, T>) == sizeof(I)
                                   ,tagged_idx_base<max_bits_v<T>, I, T>>;

}

#endif
