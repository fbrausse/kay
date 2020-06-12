
#ifndef KAY_NUM_BITS_HH
#define KAY_NUM_BITS_HH

#include <kay/numbers.hh>

namespace kay {

namespace detail {
template <typename T, bool flt, bool sgned> struct fund_type_bit_cnt;
}

template <typename T> struct type_bit_cnt
: std::enable_if_t<std::is_fundamental_v<T>
                  ,detail::fund_type_bit_cnt<T,std::is_floating_point_v<T>
                                              ,std::is_signed_v<T>
                                            >
                  >
{};

template <typename T>
constexpr inline size_t flt_prec(const T &v)
{
	return type_bit_cnt<T>{}(v);
}

namespace detail {

template <typename T> struct fund_type_bit_cnt<T,true,true> {

	static_assert(FLT_RADIX == 2);

	static_assert(std::is_same_v<uint32_t,integral_at_least_bits_t<FLT_MANT_DIG>>);
	static_assert(std::is_same_v<uint64_t,integral_at_least_bits_t<DBL_MANT_DIG>>);
	static_assert(std::is_same_v<uint64_t,integral_at_least_bits_t<LDBL_MANT_DIG>>);

	size_t operator()(T v) const
	{
		int exp;
		switch (std::fpclassify(v)) {
		case FP_NAN:
		case FP_INFINITE:
		case FP_ZERO:
			return 0;
		case FP_SUBNORMAL:
		case FP_NORMAL:
			static_assert(std::numeric_limits<T>::radix == 2);
			constexpr size_t N = std::numeric_limits<T>::digits;
			return flt_prec(static_cast<integral_at_least_bits_t<N>>(
			        std::ldexp(std::frexp(std::fabs(v), &exp), N)));
		}
	}
};

template <typename T> struct type_u_ul_ull_bits0;

template <> struct type_u_ul_ull_bits0<unsigned> {

	static constexpr size_t clz(unsigned v) { return __builtin_clz(v); }
	static constexpr size_t ctz(unsigned v) { return __builtin_ctz(v); }
};

template <> struct type_u_ul_ull_bits0<unsigned long> {

	static constexpr size_t clz(unsigned long v) { return __builtin_clzl(v); }
	static constexpr size_t ctz(unsigned long v) { return __builtin_ctzl(v); }
};

template <> struct type_u_ul_ull_bits0<unsigned long long> {

	static constexpr size_t clz(unsigned long long v) { return __builtin_clzll(v); }
	static constexpr size_t ctz(unsigned long long v) { return __builtin_ctzll(v); }
};

template <typename T> struct fund_type_bit_cnt<T,false,false> {

	constexpr size_t operator()(T v) const
	{
		return v ? _bits(v) - type_u_ul_ull_bits0<T>::ctz(v) : 0;
	}

private:
	static constexpr size_t _bits(T v)
	{
		return type_bits_v<T> - type_u_ul_ull_bits0<T>::clz(v);
	}
};

template <typename T> struct fund_type_bit_cnt<T,false,true> {

	constexpr size_t operator()(T v) const
	{
		if constexpr (is_twos_complement)
			if (v == std::numeric_limits<T>::min())
				return type_bits_v<T>;
		return flt_prec(static_cast<std::make_unsigned_t<T>>(std::abs(v)));
	}

private:
	static constexpr bool is_twos_complement =
		std::numeric_limits<T>::min() < -std::numeric_limits<T>::max();

	static constexpr size_t _bits(T v)
	{
		if constexpr (is_twos_complement)
			if (v == std::numeric_limits<T>::min())
				return type_bits_v<T>;
		return bits(static_cast<std::make_unsigned_t<T>>(std::abs(v)));
	}
};
}

template <> struct type_bit_cnt<Z> {

	size_t operator()(const Z &v) const
	{
		return v ? bits(v) - ctz(v) : 0;
	}
};

}

#endif
