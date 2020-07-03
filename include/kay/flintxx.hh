/*
 * kay-flintxx.hh
 *
 * Copyright 2018-2019 Franz Brauße <brausse@informatik.uni-trier.de>
 *
 * This file is part of kay.
 * See the LICENSE file for terms of distribution.
 */

#ifndef KAY_FLINTXX_HH
#define KAY_FLINTXX_HH

#if __has_include(<flint/fmpz.h>) && __has_include(<flint/fmpq.h>)
# define KAY_HAVE_FLINT 1
#endif

#include <kay/gmpxx.hh>

#if KAY_HAVE_FLINT && KAY_HAVE_GMPXX

#include <utility>	/* std::swap */
#include <ostream>

#include <flint/fmpz.h>
#include <flint/fmpq.h>

namespace kay::flintxx {

class Z {

	fmpz z;

public:
	Z()           noexcept   { fmpz_init(get_fmpz_t()); }
	Z(const Z &v) noexcept   { fmpz_init_set(get_fmpz_t(), v.get_fmpz_t()); }
	Z(Z &&v)      noexcept   : z(v.z) { v.z = 0; /* flint internals: 0 is not alloc'ed */ }

#if COEFF_MIN <= INT_MIN && INT_MAX <= COEFF_MAX
	Z(signed v)              : z(v) {}
#else
	Z(signed v)              : Z(static_cast<signed long>(v)) {}
#endif
#if COEFF_MIN <= UINT_MIN && UINT_MAX <= COEFF_MAX
	Z(unsigned v)            : z(v) {}
#else
	Z(unsigned v)            : Z(static_cast<unsigned long>(v)) {}
#endif
	Z(signed long v)         : Z() { fmpz_set_si(get_fmpz_t(), v); }
	Z(unsigned long v)             { fmpz_init_set_ui(get_fmpz_t(), v); }
	explicit Z(mpz_srcptr v) : Z() { fmpz_set_mpz(get_fmpz_t(), v); }
	Z(const mpz_class &v)    : Z(v.get_mpz_t()) {}

	explicit Z(const char *s, int base=10)
	: Z() { fmpz_set_str(get_fmpz_t(), s, base); }

	~Z() { fmpz_clear(get_fmpz_t()); }

	friend void swap(Z &a, Z &b)
	{
		fmpz_swap(a.get_fmpz_t(), b.get_fmpz_t());
	}

	Z & operator=(const Z &v) noexcept { fmpz_set(get_fmpz_t(), v.get_fmpz_t()); return *this; }
	Z & operator=(Z &&v)      noexcept { swap(*this, v); return *this; }

	      fmpz * get_fmpz_t()       { return &z; }
	const fmpz * get_fmpz_t() const { return &z; }

	explicit operator bool() const { return !fmpz_is_zero(get_fmpz_t()); }

	explicit operator mpz_class() const
	{
		mpz_class r;
		fmpz_get_mpz(r.get_mpz_t(), get_fmpz_t());
		return r;
	}

	friend void neg(Z &a)
	{
		/* flint impl fmpz_neg(tgt, src) is not efficient when
		 * tgt == src */
		if (COEFF_IS_MPZ(a.z)) {
			mpz_ptr p = COEFF_TO_PTR(a.z);
			mpz_neg(p, p);
			a.z = PTR_TO_COEFF(p);
		} else {
			a.z = -a.z;
		}
	}

	friend Z operator+(Z a) { return a; }
	friend Z operator~(Z a) { fmpz_complement(a.get_fmpz_t(), a.get_fmpz_t()); return a; }
	friend Z operator-(Z a) { neg(a); return a; }

	Z & operator++()    { *this += 1; return *this; }
	Z   operator++(int) { Z old = *this; ++*this; return old; }

	Z & operator--()    { *this -= 1; return *this; }
	Z   operator--(int) { Z old = *this; --*this; return old; }

	friend Z & operator+=(Z &a, const Z &b)
	{ fmpz_add(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator+ (Z  a, const Z &b) { a += b; return a; }

	friend Z & operator-=(Z &a, const Z &b)
	{ fmpz_sub(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator- (Z  a, const Z &b) { a -= b; return a; }

	friend Z & operator*=(Z &a, const Z &b)
	{ fmpz_mul(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator* (Z  a, const Z &b) { a *= b; return a; }

	friend Z & operator/=(Z &a, const Z &b)
	{ fmpz_tdiv_q(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator/ (Z  a, const Z &b) { a /= b; return a; }

	friend Z & operator%=(Z &a, const Z &b)
	{ fmpz_mod(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator% (Z  a, const Z &b) { a %= b; return a; }

	template <typename T
	         ,typename = std::enable_if_t<std::is_unsigned_v<std::remove_cv_t<T>> &&
	                                      (type_bits_v<std::remove_cv_t<T>> <=
	                                       type_bits_v<::ulong>) &&
	                                      std::is_convertible_v<T,::ulong> &&
	                                      std::is_convertible_v<::ulong,T>>>
	friend T   operator% (Z  a, const T &b)
	{
		/* flint naming sometimes is strange, this is the modulo
		 * operation, though. */
		return fmpz_fdiv_ui(a.get_fmpz_t(), b);
	}

	friend Z & operator<<=(Z &a, mp_bitcnt_t e)
	{ fmpz_mul_2exp(a.get_fmpz_t(), a.get_fmpz_t(), e); return a; }
	friend Z   operator<< (Z  a, mp_bitcnt_t e) { a <<= e; return a; }

	friend Z & operator>>=(Z &a, mp_bitcnt_t e)
	{ fmpz_tdiv_q_2exp(a.get_fmpz_t(), a.get_fmpz_t(), e); return a; }
	friend Z   operator>> (Z  a, mp_bitcnt_t e) { a >>= e; return a; }

	friend Z & operator&=(Z &a, const Z &b)
	{ fmpz_and(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator& (Z  a, const Z &b) { a &= b; return a; }

	friend Z & operator|=(Z &a, const Z &b)
	{ fmpz_or(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator| (Z  a, const Z &b) { a |= b; return a; }

	friend Z & operator^=(Z &a, const Z &b)
	{ fmpz_xor(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }
	friend Z   operator^ (Z  a, const Z &b) { a ^= b; return a; }

	friend int cmp(const Z &a, const Z &b) { return fmpz_cmp(a.get_fmpz_t(), b.get_fmpz_t()); }

	friend int sgn(const Z &a)
	{
		/* fmpz_sgn is not inlined */
		if (COEFF_IS_MPZ(a.z)) {
			return mpz_sgn(COEFF_TO_PTR(a.z));
		} else {
			return a.z < 0 ? -1 : a.z > 0 ? +1 : 0;
		}
	}

	friend Z   pow(Z a, unsigned long x) { fmpz_pow_ui(a.get_fmpz_t(), a.get_fmpz_t(), x); return a; }
	friend Z   abs(Z a)                  { fmpz_abs(a.get_fmpz_t(), a.get_fmpz_t()); return a; }
	friend Z   gcd(Z a, const Z &b)      { fmpz_gcd(a.get_fmpz_t(), a.get_fmpz_t(), b.get_fmpz_t()); return a; }

	friend size_t sizeinbase(const Z &a, int base)
	{
		return fmpz_sizeinbase(a.get_fmpz_t(), base);
	}

	friend mp_bitcnt_t bits(const Z &a) { return fmpz_bits(a.get_fmpz_t()); }

	friend mp_bitcnt_t ctz(const Z &a) { return fmpz_val2(a.get_fmpz_t()); }

	friend bool operator==(const Z &a, const Z &b) { return cmp(a, b) == 0; }
	friend bool operator!=(const Z &a, const Z &b) { return cmp(a, b) != 0; }
	friend bool operator<=(const Z &a, const Z &b) { return cmp(a, b) <= 0; }
	friend bool operator< (const Z &a, const Z &b) { return cmp(a, b) <  0; }
	friend bool operator>=(const Z &a, const Z &b) { return cmp(a, b) >= 0; }
	friend bool operator> (const Z &a, const Z &b) { return cmp(a, b) >  0; }

	std::string get_str(int base=10) const { return fmpz_get_str(NULL, base, get_fmpz_t()); }

	friend std::ostream & operator<<(std::ostream &os, const Z &v)
	{
		/* TODO: inefficient */
		std::ios_base::fmtflags flags = os.flags();
		int base = 10;
		if (flags & os.oct) base = 8;
		if (flags & os.hex) base = 16;
		if (flags & os.showpos && sgn(v) >= 0)
			os << "+";
		if (flags & os.showbase)
			switch (base) {
			case 8: os << "0";
			case 16: os << (flags & os.uppercase ? "0X" : "0x");
			}
		char *s = fmpz_get_str(NULL, base, v.get_fmpz_t());
		os << s;
		free(s);
		return os;
	}
};

/* 0^0 yields 1 */
inline Z ui_pow_ui(unsigned long base, unsigned long exp)
{
	return pow(Z(base), exp);
}

struct Q {

	Z num;
	Z den;

	Q()           noexcept : num(), den(1U) {}
	Q(const Q &v) noexcept = default;
	Q(Q &&v)      noexcept = default;
	Q(Z num)               : num(std::move(num))
	                       , den(1) {}

	Q(const Z &num, const Z &den)
	: Q()
	{ fmpq_set_fmpz_frac(get_fmpq_t(), num.get_fmpz_t(), den.get_fmpz_t()); }

	Q(const Z &num, Z &&den)
	: num(num)
	, den(std::move(den))
	{ canonicalize(*this); }

	Q(Z &&num, Z den)
	: num(std::move(num))
	, den(std::move(den))
	{ canonicalize(*this); }

	Q(signed int num)  : num(num), den(1) {}
	Q(signed long num) : num(num), den(1) {}

	Q(signed long num, unsigned long den)
	: Q()
	{ fmpq_set_si(get_fmpq_t(), num, den); }

	Q(double d) : Q(mpq_class(d)) {}

	explicit Q(const char *s, int base=10)
	: Q()
	{
		using std::string;
		if (const char *delim = strchr(s, '/'))
			*this = Q(Z(string(s, delim-s).c_str(), base),
			          Z(delim+1, base));
		else
			*this = Z(s, base);
	}

	explicit Q(mpq_srcptr v) : Q() { fmpq_set_mpq(get_fmpq_t(), v); }
	Q(const mpq_class &v)    : Q(v.get_mpq_t()) {}

	~Q() {}

	friend void swap(Q &a, Q &b) { fmpq_swap(a.get_fmpq_t(), b.get_fmpq_t()); }

	Q & operator=(const Q &v) { fmpq_set(get_fmpq_t(), v.get_fmpq_t()); return *this; }
	Q & operator=(Q &&v) = default;

	constexpr       Z & get_num()       { return num; }
	constexpr const Z & get_num() const { return num; }
	constexpr       Z & get_den()       { return den; }
	constexpr const Z & get_den() const { return den; }

	      fmpq * get_fmpq_t()       { return reinterpret_cast<fmpq *>(this); }
	const fmpq * get_fmpq_t() const { return reinterpret_cast<const fmpq *>(this); }

	friend void canonicalize(Q &a) { fmpq_canonicalise(a.get_fmpq_t()); }

	explicit operator bool() const { return !fmpq_is_zero(get_fmpq_t()); }

	explicit operator mpq_class() const
	{
		mpq_class r;
		fmpq_get_mpq(r.get_mpq_t(), get_fmpq_t());
		return r;
	}

	friend void neg(Q &a) { neg(a.num); }

	friend Q operator+(Q a) { return a; }
	friend Q operator-(Q a) { neg(a); return a; }

	Q & operator++()    { num += den; return *this; }
	Q   operator++(int) { Q old = *this; ++*this; return old; }

	Q & operator--()    { num -= den; return *this; }
	Q   operator--(int) { Q old = *this; --*this; return old; }

	friend Q & operator+=(Q &a, const Q &b)
	{ fmpq_add(a.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); return a; }
	friend Q   operator+ (Q  a, const Q &b) { a += b; return a; }

	friend Q & operator-=(Q &a, const Q &b)
	{ fmpq_sub(a.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); return a; }
	friend Q   operator- (Q  a, const Q &b) { a -= b; return a; }

	friend Q & operator*=(Q &a, const Q &b)
	{ fmpq_mul(a.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); return a; }
	friend Q   operator* (Q  a, const Q &b) { a *= b; return a; }

	friend Q & operator/=(Q &a, const Q &b)
	{ fmpq_div(a.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); return a; }
	friend Q   operator/ (Q  a, const Q &b) { a /= b; return a; }

	friend Q & operator<<=(Q &a, mp_bitcnt_t e)
	{ fmpq_mul_2exp(a.get_fmpq_t(), a.get_fmpq_t(), e); return a; }
	friend Q   operator<< (Q  a, mp_bitcnt_t e) { a <<= e; return a; }

	friend Q & operator>>=(Q &a, mp_bitcnt_t e)
	{ fmpq_div_2exp(a.get_fmpq_t(), a.get_fmpq_t(), e); return a; }
	friend Q   operator>> (Q  a, mp_bitcnt_t e) { a >>= e; return a; }

	friend void fma(Q &r, const Q &a, const Q &b)
	{ fmpq_addmul(r.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); }
	friend void fms(Q &r, const Q &a, const Q &b)
	{ fmpq_submul(r.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); }

	friend int  sgn(const Q &a)
	{
		/* fmpq_sgn just delegates to fmpz_sgn */
		return sgn(a.num);
	}

	friend int  cmp(const Q &a, const Q &b) { return fmpq_cmp(a.get_fmpq_t(), b.get_fmpq_t()); }
	friend Q    inv(Q a)                    { fmpq_inv(a.get_fmpq_t(), a.get_fmpq_t()); return a; }
	friend Q    abs(Q a)                    { fmpq_abs(a.get_fmpq_t(), a.get_fmpq_t()); return a; }
	friend Q    gcd(Q a, const Q &b)        { fmpq_gcd(a.get_fmpq_t(), a.get_fmpq_t(), b.get_fmpq_t()); return a; }
	friend Q    pow(Q a, signed long e)     { fmpq_pow_si(a.get_fmpq_t(), a.get_fmpq_t(), e); return a; }

	friend bool operator==(const Q &a, const Q &b) { return cmp(a, b) == 0; }
	friend bool operator!=(const Q &a, const Q &b) { return cmp(a, b) != 0; }
	friend bool operator<=(const Q &a, const Q &b) { return cmp(a, b) <= 0; }
	friend bool operator< (const Q &a, const Q &b) { return cmp(a, b) <  0; }
	friend bool operator>=(const Q &a, const Q &b) { return cmp(a, b) >= 0; }
	friend bool operator> (const Q &a, const Q &b) { return cmp(a, b) >  0; }

	/* TODO: rational reconstruction; flint-2.5.2 manual 25.10 */
	/* TODO: continued fractions; flint-2.5.2 manual 25.12 */

	std::string get_str(int base=10) const { return fmpq_get_str(NULL, base, get_fmpq_t()); }

	/* truncates, i.e. rounds towards zero */
	double get_d() const
	{
		/* TODO: inefficient */
		return static_cast<mpq_class>(*this).get_d();
	}

	friend int mpfr_set_q(mpfr_t dest, const Q &src, mpfr_rnd_t rnd)
	{
		return fmpq_get_mpfr(dest, src.get_fmpq_t(), rnd);
	}

	friend int mpfr_sub_q(mpfr_t r, mpfr_t a, const Q &b, mpfr_rnd_t rnd)
	{
		/* TODO: inefficient */
		return mpfr_sub_q(r, a, static_cast<mpq_class>(b).get_mpq_t(), rnd);
	}

	friend Z floor(const Q &q)
	{
		Z r;
		fmpz_fdiv_q(r.get_fmpz_t(), q.get_num().get_fmpz_t(),
		                            q.get_den().get_fmpz_t());
		return r;
	}

	friend Z ceil(const Q &q)
	{
		Z r;
		fmpz_cdiv_q(r.get_fmpz_t(), q.get_num().get_fmpz_t(),
		                            q.get_den().get_fmpz_t());
		return r;
	}

	friend Z round(const Q &q)
	{
		return floor(Q(1,2)+q);
	}

	friend std::ostream & operator<<(std::ostream &os, const Q &v)
	{
		/* TODO: inefficient */
		/* TODO: obey os.flags() */
		char *s = fmpq_get_str(NULL, 10, v.get_fmpq_t());
		os << s;
		free(s);
		return os;
	}
};

static_assert(sizeof(Q) == sizeof(fmpq));

static_assert(std::is_standard_layout_v<Z>);
static_assert(std::is_standard_layout_v<Q>);

/* clang-7 thinks these expressions are not constant
static_assert((void *)&((Q *)nullptr)->get_num() == (void *)&((fmpq *)nullptr)->num);
static_assert((void *)&((Q *)nullptr)->get_den() == (void *)&((fmpq *)nullptr)->den);
*/

/* TODO: algebraic numbers using fmpz_poly */

}

namespace std {

template <>
struct hash<kay::flintxx::Z> : protected hash<mpz_class> {

	size_t operator()(const fmpz *v) const noexcept
	{
		if (!COEFF_IS_MPZ(*v))
			return *v;
		return hash<mpz_class>::operator()(COEFF_TO_PTR(*v));
	}

	size_t operator()(const kay::flintxx::Z &v) const noexcept
	{
		return operator()(v.get_fmpz_t());
	}

protected:
	size_t combine(size_t r, const fmpz *v) const noexcept
	{
		if (!COEFF_IS_MPZ(*v))
			return kay::fnv1a_hash<size_t>::combine(r, *v);
		return hash<mpz_class>::combine(r, COEFF_TO_PTR(*v));
	}
};

template <>
struct hash<kay::flintxx::Q> : protected hash<kay::flintxx::Z> {

	size_t operator()(const kay::flintxx::Q &v) const noexcept
	{
		return combine(offset_basis, v.get_fmpq_t());
	}

protected:
	size_t combine(size_t r, const fmpq *v) const noexcept
	{
		r = hash<kay::flintxx::Z>::combine(r, &v->num);
		r = hash<kay::flintxx::Z>::combine(r, &v->den);
		return r;
	}
};

}

#endif /* KAY_HAVE_FLINT && KAY_HAVE_GMPXX */

#endif
