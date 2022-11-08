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

#if defined(__GNUC__)
/* officially supported by GCC since 4.5, all GNUC-compliant compilers with
 * C++11-support know about this built-in */
# define kay_unreachable()	__builtin_unreachable()
#else
/* https://stackoverflow.com/questions/6031819/emulating-gccs-builtin-unreachable */
[[noreturn]] static inline void kay_unreachable() { /* intentional */ }
#endif

#include <charconv>	/* std::from_chars_result */
#include <cassert>

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

inline Z pow(Z x, unsigned long n)
{
	mpz_pow_ui(x.get_mpz_t(), x.get_mpz_t(), n);
	return x;
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
static inline Q Q_from_str(char *rep, unsigned base = 10)
{
	const char *elit = base == 10 ? "eE" : base == 16 ? "p" : nullptr;
	char *e = rep + (elit ? strcspn(rep, elit) : 0);
	if (*e)
		*e++ = '\0';
	else
		e = NULL;
	char *dot = strchr(rep, '.');
	Q divisor = 1;
	if (dot) {
		unsigned frac_len = strlen(dot+1);
		divisor = ui_pow_ui(base, frac_len);
		memmove(dot, dot+1, frac_len+1);
	}
	Q r(rep, base);
	r /= divisor;
	if (e) {
		long g = strtol(e, NULL, base);
		Z f = pow(Z(base), labs(g));
		if (g < 0)
			r /= f;
		else
			r *= f;
	}
	return r;
}

inline std::from_chars_result
from_chars(const char *rep, const char *end, Z &v, int base = 0,
           bool incl_sign = true, bool incl_prefix = true)
{
	assert(base > 1);
	assert(base < 36);
	if (rep == end)
		return { rep, std::errc::invalid_argument };
	bool is_neg = false;
	const char *st = rep;
	if (incl_sign) {
		is_neg = *st == '-';
		if (strchr("+-", *st))
			st++;
	}
	if (st == end || !isdigit(*st))
		return { rep, std::errc::invalid_argument };
	if (incl_prefix) {
		int new_base;
		if (end - rep >= 2 && st[0] == '0' && tolower(st[1]) == 'x') {
			new_base = 16;
			st += 2;
		} else if (st[0] == '0') {
			new_base = 8;
			st += 1;
		} else
			new_base = 10;
		if (base && base != new_base)
			return { rep, std::errc::invalid_argument };
		base = new_base;
	}
	if (st == end || !isdigit(*st))
		return { rep, std::errc::invalid_argument };
	v = 0;
	const char *beg = st;
	for (; st < end; st++) {
		char c = tolower(*st);
		int digit;
		if ('0' <= c && c <= '9')
			digit = c - '0';
		else if ('a' <= c && c <= 'z')
			digit = 10 + (c - 'a');
		else
			break;
		if (digit >= base)
			break;
		v *= base;
		v += digit;
	}
	if (beg == st)
		return { rep, std::errc::invalid_argument };
	if (is_neg)
		neg(v);
	return { st, std::errc {} };
}

namespace detail {

inline std::from_chars_result
from_chars_Q_component(const char *rep, const char *end, Q &v, int base)
{
	const char *beg = rep;
	bool is_neg = *rep == '-';
	if (strchr("+-", *rep))
		rep++;
	v = 0;
	auto [ze,zr] = from_chars(rep, end, v.get_num(), base, true, false);
	if (zr != std::errc {})
		return { beg, zr };
	if (ze == end) {
		if (is_neg)
			neg(v);
		return { ze, zr };
	}
	assert(ze < end);
	const char *st = ze;
	if (*st == '.' && isdigit(st[1])) {
		Z f;
		auto [fe,fr] = from_chars(st + 1, end, f, base, false, false);
		if (fr != std::errc {}) {
			end = st;
		} else {
			size_t frac_len = fe - (ze + 1);
			v += Q(v < 0 ? -f : f, ui_pow_ui(base, frac_len));
			st = fe;
		}
	}
	const char *elit = base == 10 ? "eE" : base == 16 ? "p" : nullptr;
	if (st < end && elit && strchr(elit, *st)) {
		long e;
		auto [ee,er] = std::from_chars(st + 1, end, e, base);
		if (er != std::errc {}) {
			end = st;
		} else {
			Z f = pow(Z(base), labs(e));
			if (e < 0)
				v /= f;
			else
				v *= f;
			st = ee;
		}
	}
	if (is_neg)
		neg(v);
	return { st, std::errc {} };
}
} /* namespace detail */

inline std::from_chars_result
from_chars(const char *rep, const char *end, Q &v, int base = 10)
{
	auto r = detail::from_chars_Q_component(rep, end, v, base);
	if (r.ec != std::errc {})
		return { rep, r.ec };
	const char *s = r.ptr;
	if (*s == '/') {
		Q d;
		r = detail::from_chars_Q_component(s+1, end, d, base);
		if (r.ec == std::errc {}) {
			s = r.ptr;
			v /= d;
		}
	}
	return { s, {} };
}

inline Q scale(Q v, ssize_t n)
{
	if (n > 0)
		v <<= n;
	else if (n)
		v >>= -n;
	return v;
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
static_assert(digits<~(size_t)0,10>::value == (sizeof(size_t) == 4 ? 10 : 20));
static_assert(digits<(~(size_t)0 >> 2),2>::value == (sizeof(size_t) * CHAR_BIT - 2));

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


static_assert(sizeof(size_t) != 8 || parse_v<'1','8','4','4','6','7','4','4','0','7','3','7','0','9','5','5','1','6','1','5'> == 18446744073709551615U);
static_assert(sizeof(size_t) != 4 || parse_v<'4','2','9','4','9','6','7','2','9','5'> == 4294967295U);
static_assert(sizeof(size_t) > 8 || parse_v<'1','8','4','4','6','7','4','4','0','7','3','7','0','9','5','5','1','6','1','6'> == 0);


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
