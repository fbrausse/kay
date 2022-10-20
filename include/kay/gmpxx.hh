/*
 * kay/gmpxx.hh
 *
 * Copyright 2018-2019 Franz Brauße <brausse@informatik.uni-trier.de>
 *
 * This file is part of kay.
 * See the LICENSE file for terms of distribution.
 */

#ifndef KAY_GMPXX_HH
#define KAY_GMPXX_HH

#if __has_include(<gmpxx.h>)
# define KAY_HAVE_GMPXX 1
#endif

#if __has_include(<mpfr.h>)
# define KAY_HAVE_MPFR 1
#endif

#if KAY_HAVE_GMPXX
# include <gmpxx.h>
#endif
#if KAY_HAVE_MPFR
# include <mpfr.h>	/* mpfr_t */
#endif

#if KAY_HAVE_GMPXX

#include <kay/bits.hh>

namespace kay {

inline void neg(mpz_class &v) { v = -v; }

inline void neg(mpq_class &v) { v = -v; }

inline void fma(mpq_class &r, const mpq_class &a, const mpq_class &b)
{
	r += a * b;
}

/* unconditionally define these for GMP as long as cont-frac is not converted */

inline void floor(mpz_class &i, const mpq_class &q)
{
	mpz_fdiv_q(i.get_mpz_t(), q.get_num_mpz_t(), q.get_den_mpz_t());
}

inline mpz_class floor(const mpq_class &q)
{
	mpz_class z;
	floor(z, q);
	return z;
}

inline void ceil(mpz_class &i, const mpq_class &q)
{
	mpz_cdiv_q(i.get_mpz_t(), q.get_num_mpz_t(), q.get_den_mpz_t());
}

inline mpz_class ceil(const mpq_class &q)
{
	mpz_class z;
	ceil(z, q);
	return z;
}

inline mpz_class round(const mpq_class &q)
{
	return kay::floor(q+mpq_class(1,2));
}

inline size_t bits(const mpz_class &v)
{
	return mpz_sizeinbase(v.get_mpz_t(), 2);
}

#if KAY_HAVE_MPFR

inline int mpfr_set_q(mpfr_t dest, const mpq_class &src, mpfr_rnd_t rnd)
{
	return mpfr_set_q(dest, src.get_mpq_t(), rnd);
}

inline int mpfr_sub_q(mpfr_t r, mpfr_t a, const mpq_class &b, mpfr_rnd_t rnd)
{
	return mpfr_sub_q(r, a, b.get_mpq_t(), rnd);
}

#endif /* KAY_HAVE_MPFR */

}

namespace std {

template <>
struct hash<mpz_class> : enable_if_t<sizeof(mp_limb_t) == sizeof(size_t),
                                     kay::fnv1a_hash<size_t>> {

	size_t operator()(const mpz_class &v) const noexcept
	{
		return operator()(v.get_mpz_t());
	}

	size_t operator()(mpz_srcptr v) const noexcept
	{
		return combine(offset_basis, v);
	}

protected:
	size_t combine(size_t r, mpz_srcptr z) const noexcept
	{
		int s = z->_mp_size;
		unsigned n;
		if (s < 0) {
			r = kay::fnv1a_hash<size_t>::combine(r, 1);
			n = -s;
		} else
			n = s;
		for (unsigned i=0; i<n; i++)
			r = kay::fnv1a_hash<size_t>::combine(r, z->_mp_d[i]);
		return r;
	}
};

template <>
struct hash<mpq_class> : hash<mpz_class> {

	size_t operator()(const mpq_class &q) const noexcept
	{
		size_t r = offset_basis;
		r = combine(r, q.get_num_mpz_t());
		r = combine(r, q.get_den_mpz_t());
		return r;
	}
};

}

#endif /* KAY_HAVE_GMPXX */

#endif
