/*
 * numbers.hh
 *
 * Copyright 2018-2019 Franz Brauße <brausse@informatik.uni-trier.de>
 *
 * This file is part of kay.
 * See the LICENSE file for terms of distribution.
 */

#ifndef KAY_NUMBERS_HH
#define KAY_NUMBERS_HH

#include <kay/bits.hh>
#include <kay/gmpxx.hh>
#include <kay/flintxx.hh>
#include <kay/compiletime.hh>

#if !defined(KAY_USE_FLINT) && !defined(KAY_USE_GMPXX)
# if KAY_HAVE_FLINT
#  define KAY_USE_FLINT 1
# elif KAY_HAVE_GMPXX
#  define KAY_USE_GMPXX 1
# endif
#endif

#if (KAY_USE_FLINT-0) && (KAY_USE_GMPXX-0)
# error "cannot use both, flint and gmpxx"
#elif (KAY_USE_FLINT-0)

# include <kay/flintxx.hh>

namespace kay {
using Z = flintxx::Z;
using Q = flintxx::Q;

using flintxx::ui_pow_ui;

inline mpz_class to_mpz_class(const Z &z) { return static_cast<mpz_class>(z); }
inline mpq_class to_mpq_class(const Q &q) { return static_cast<mpq_class>(q); }

}

#elif (KAY_USE_GMPXX-0)

# include <kay/gmpxx.hh>

namespace kay {
using Z = ::mpz_class;
using Q = ::mpq_class;

inline Z ui_pow_ui(unsigned long base, unsigned long exp)
{
	mpz_class r;
	mpz_ui_pow_ui(r.get_mpz_t(), base, exp);
	return r;
}

inline Q pow(Q x, unsigned long n)
{
	mpz_pow_ui(x.get_num_mpz_t(), x.get_num_mpz_t(), n);
	mpz_pow_ui(x.get_den_mpz_t(), x.get_den_mpz_t(), n);
	/* x was canonicalized before, so it is now */
	return x;
}

inline mp_bitcnt_t sizeinbase(const Z &a, int base)
{
	return mpz_sizeinbase(a.get_mpz_t(), base);
}

inline void canonicalize(Q &v)
{
	v.canonicalize();
}

inline const mpz_class & to_mpz_class(const Z &z) { return z; }
inline       mpz_class & to_mpz_class(      Z &z) { return z; }

inline const mpq_class & to_mpq_class(const Q &q) { return q; }
inline       mpq_class & to_mpq_class(      Q &q) { return q; }

inline const mpq_class   inv(mpq_class q) { return 1/q; }

inline       mp_bitcnt_t ctz(const mpz_class &v) { return mpz_scan1(v.get_mpz_t(), 0); }

}

#else
# error "neither KAY_USE_FLINT nor KAY_USE_GMPXX; cannot define numbers"
#endif

namespace kay {

/* parses stuff like "0.85" */
static Q Q_from_str(char *rep, unsigned base = 10)
{
	char *dot = strchr(rep, '.');
	Q divisor = 1;
	if (dot) {
		unsigned frac_len = strlen(dot+1);
		divisor = ui_pow_ui(base, frac_len);
		memmove(dot, dot+1, frac_len+1);
	}
	Q r(rep, base);
	r /= divisor;
	return r;
}


namespace integral_literal_support {

template <size_t v> using cnst = std::integral_constant<size_t,v>;

template <char min, char max, char c> struct alph_help
: std::enable_if_t<(min <= c && c <= max),cnst<c-min>> {};

template <size_t b,char c> struct alph;
template <char c> struct alph<2 ,c> : alph_help<'0','1',c> {};
template <char c> struct alph<8 ,c> : alph_help<'0','7',c> {};
template <char c> struct alph<10,c> : alph_help<'0','9',c> {};
template <char c> struct alph<16,c>
: std::conditional_t<(c >= 'A' && c <= 'F'),cnst<(size_t)(c-'A'+10)>,
  std::conditional_t<(c >= 'a' && c <= 'f'),cnst<(size_t)(c-'a'+10)>,
  alph<10,c>>> {};

template <size_t b,char... cs> struct base : cnst<1> {};
template <size_t b,char c, char... cs>
struct base<b,c,cs...> : cnst<b*base<b,cs...>::value> {};

template <size_t b,char... cs> struct parse2;
template <char... cs> struct parse2<0,'0','x',cs...> : parse2<16,cs...> {};
template <char... cs> struct parse2<0,'0','X',cs...> : parse2<16,cs...> {};
template <char... cs> struct parse2<0,'0','b',cs...> : parse2<2,cs...> {};
template <char... cs> struct parse2<0,'0','B',cs...> : parse2<2,cs...> {};
template <char... cs> struct parse2<0,'0',cs...> : parse2<8,cs...> {};
template <char c, char... cs> struct parse2<0,c,cs...> : parse2<10,c,cs...> {};
template <size_t b,char c, char... cs>
struct parse2<b,c,cs...>
: cnst<parse2<b,cs...>::value+alph<b,c>::value*base<b,cs...>::value> {};
template <size_t b> struct parse2<b> : cnst<0> {};

template <char... cs> struct parse : parse2<0,cs...> {};
template <char... cs> constexpr size_t parse_v = parse<cs...>::value;

template <size_t v, size_t b>
struct digits : cnst<(v<b) ? 1 : 1+digits<v/b,b>::value> {};

template <size_t b> struct digits<0,b> : cnst<1> {};

static_assert(digits<255,16>::value == 2);
static_assert(digits<256,16>::value == 3);

static_assert(digits<9,10>::value == 1);
static_assert(digits<10,10>::value == 2);
static_assert(digits<~(uint64_t)0,10>::value == 20);
static_assert(digits<(~(uint64_t)0 >> 2),2>::value == 62);

template <char... cs> struct num_nonzero : cnst<sizeof...(cs)> {};
template <char... cs> struct num_nonzero<'0',cs...> : num_nonzero<cs...> {};


template <size_t, int, char...> struct fits3;
template <size_t max> struct fits3<max,0,'0'> : std::true_type {};

template <size_t max, char... cs>
struct fits3<max,0,'0','x',cs...> : fits3<max,16,cs...> {};

template <size_t max, char... cs>
struct fits3<max,0,'0','X',cs...> : fits3<max,16,cs...> {};

template <size_t max, char... cs>
struct fits3<max,0,'0','b',cs...> : fits3<max,2,cs...> {};

template <size_t max, char... cs>
struct fits3<max,0,'0','B',cs...> : fits3<max,2,cs...> {};

template <size_t max, char... cs>
struct fits3<max,0,'0',cs...> : fits3<max,8,cs...> {};

template <size_t max, char c, char... cs>
struct fits3<max,0,c,cs...> : fits3<max,10,c,cs...> {};

template <size_t max, int b, char c, char... cs>
struct fits3<max,b,c,cs...>
: std::bool_constant<(num_nonzero<c,cs...>::value <= digits<max,b>::value) &&
                     ((max - parse2<b,cs...>::value) / base<b,cs...>::value
                      >= alph<b,c>::value)> {};


template <typename UL, char... cs>
constexpr bool fits_v = fits3<std::numeric_limits<UL>::max(),0,cs...>::value;

static_assert(fits_v<unsigned,'0'>);
static_assert(fits_v<unsigned,'0','0'>);
static_assert((UINT32_MAX-parse_v<'2','9','4','9','6','7','2','9','6'>) == 4e9-1);
static_assert(fits_v<uint32_t,'4','2','9','4','9','6','7','2','9','5'>);
static_assert(fits_v<uint32_t,'0','3','7','7','7','7','7','7','7','7','7','7'>);
static_assert(!fits_v<uint32_t,'0','4','0','0','0','0','0','0','0','0','0','0'>);
static_assert(fits_v<uint32_t,'0','x','f','f','f','f','f','f','f','f'>);
static_assert(!fits_v<uint32_t,'0','x','1','0','0','0','0','0','0','0','0'>);
static_assert(!fits_v<uint32_t,'4','2','9','4','9','6','7','2','9','6'>);

static_assert(parse_v<'1','8','4','4','6','7','4','4','0','7','3','7','0','9','5','5','1','6','1','5'> == 18446744073709551615U);
static_assert(parse_v<'1','8','4','4','6','7','4','4','0','7','3','7','0','9','5','5','1','6','1','6'> == 0);


template <char... cs>
constexpr std::enable_if_t<fits_v<unsigned long,cs...>,Z> parse_Z()
{ return parse_v<cs...>; }

template <char... cs>
inline std::enable_if_t<!fits_v<unsigned long,cs...>,Z> parse_Z()
{ char str[] = {cs...,'\0'}; return Z(str); }

};

namespace literals {

template <char... cs>
inline Z operator""_Z() { return integral_literal_support::parse_Z<cs...>(); }

}

}

#endif
