/*
 * ival.hh
 *
 * Copyright 2018-2022 Franz Brauße <fb@paxle.org>
 *
 * See the LICENSE file for terms of distribution.
 */

#ifndef KAY_DBL_IVAL_HH
#define KAY_DBL_IVAL_HH

#include <cfenv>	/* fe[gs]etround() */
#include <cmath>	/* INFINITY */
#include <sstream>
#include <cassert>
#include <kay/numbers.hh>
#include <kay/numbits.hh>

namespace kay::dbl {

using kay::Z;
using kay::Q;

/* double intervals */

/* no known compiler claims to understands these pragmas required by the
 * C standard for our use-case
#pragma STDC FENV_ACCESS ON
#pragma STDC FP_CONTRACT OFF
// C2x pragma to set the "static rounding direction" of a part of code
#pragma STDC FENV_ROUND FE_DOWNWARD
 */

class rounding_mode {

	const int old;
	const bool changed;

public:
	rounding_mode(int mode)
	: old(fegetround())
	, changed(old != mode)
	{
		if (!changed)
			return;
		if (int r = fesetround(mode)) {
			std::stringstream ss;
			ss << "fesetround(" << mode << ") failed with code "
			   << r;
			throw std::runtime_error(ss.str());
		}
	}

	rounding_mode(rounding_mode &&) = delete;

	~rounding_mode() { if (changed) fesetround(old); }

	rounding_mode & operator=(rounding_mode) = delete;
};

enum ieee1788_cmp {
	BOTH_EMPTY    = 1 <<  0, // - BOTH_EMPTY
	FIRST_EMPTY   = 1 <<  1, // - SECOND_EMPTY
	SECOND_EMPTY  = 1 <<  2, // - FIRST_EMPTY
	BEFORE        = 1 <<  3, // - AFTER
	MEETS         = 1 <<  4, // - MET_BY
	OVERLAPS      = 1 <<  5, // - OVERLAPPED_BY
	STARTS        = 1 <<  6, // - STARTED_BY
	CONTAINED_BY  = 1 <<  7, // - CONTAINS
	FINISHES      = 1 <<  8, // - FINISHED_BY
	EQUALS        = 1 <<  9, // - EQUALS
	FINISHED_BY   = 1 << 10,
	CONTAINS      = 1 << 11,
	STARTED_BY    = 1 << 12,
	OVERLAPPED_BY = 1 << 13,
	MET_BY        = 1 << 14,
	AFTER         = 1 << 15,
};

enum ival_pos {
	IVAL_LT  = BEFORE,  /* hl < 0 */
	IVAL_LE  = MEETS,  /* hl == 0 && ll < 0 && hh < 0 */
	IVAL_LO  = OVERLAPS, /* ll < 0 && hl > 0 && hh < 0 */
	IVAL_SUB = STARTS | CONTAINED_BY | FINISHES, /* ll == 0 && hh < 0 ||
	                                              * ll > 0 && hh < 0 ||
	                                              * ll > 0 && hh == 0 */
	/* ll >= 0 && hh <= 0 && ll != hh */
	IVAL_EQ  = EQUALS,  /* ll == 0 && hh == 0 */
	IVAL_SUP = FINISHED_BY | CONTAINS | STARTED_BY, /* ll < 0 && hh == 0 ||
	                                                 * ll < 0 && hh > 0 ||
	                                                 * ll == 0 && hh > 0 */
	/* ll <= 0 && hh >= 0 && ll != hh */
	IVAL_GO  = OVERLAPPED_BY, /* ll > 0 && lh < 0 && hh > 0 */
	IVAL_GE  = MET_BY, /* lh == 0 && hh > 0 && ll > 0 */
	IVAL_GT  = AFTER,  /* lh > 0 */
};

/* Assumes neither of a and b are NaN; infinities of the same sign compare equal */
constexpr inline int cmp(double a, double b)
{
	return a < b ? -1 : a > b ? +1 : 0;
}

/* -sgn([a,b]) == sgn(-[a,b]) */
enum ival_sgn : int32_t {
	NEG     = -1,
	ZERO    =  0,
	POS     =  1,
	OV_ZERO = INT32_MIN,
};

using kay::flt_prec;

template <typename C, typename R>
struct cnt_rad { C c; R r; };

template <typename C, typename R> cnt_rad(C,R) -> cnt_rad<C,R>;

struct endpts { double l, u; };

template <typename X, typename... Ts>
constexpr bool is_any = (std::is_same_v<X,Ts> || ...);

template <typename X, typename... Ts>
using op_compat_t = std::enable_if_t<is_any<std::remove_cv_t<X>,Ts...>,X>;

/* All operations on ival objects except creation from anything but
 * cnt_rad<double,double> require rounding_mode(FE_DOWNWARD).
 *
 * Represents intervals with double-precision endpoints:
 *  - non-empty (by construction)
 *  - point-intervals supported
 *  - endpoints must be finite or infinite, not NaN
 */
class ival {

	double lo_pos, hi_neg;

	constexpr ival(double lo_pos, double hi_neg)
	: lo_pos(lo_pos)
	, hi_neg(hi_neg)
	{
		assert(lo_pos <= -hi_neg);
		assert(!isempty(*this));
	}

	explicit ival(double d, int sgn)
	: ival {   sgn >= 0 ? d : nextafter(d, -INFINITY),
	         -(sgn <= 0 ? d : nextafter(d,  INFINITY)) }
	{
		assert(sgn || ispoint(*this));
	}

public:
	explicit ival(int32_t v=0) : ival { (double)v, -(double)v } {}
	explicit ival(int64_t v) : ival(v, flt_prec(v) <= DBL_MANT_DIG ? 0 : v < 0 ? -1 : v > 0 ? +1 : 0) {}
	explicit ival(const Z &v) : ival(Q(v).get_d(), flt_prec(v) <= DBL_MANT_DIG ? 0 : sgn(v)) {}
	explicit ival(const Q &v) : ival(v.get_d(), sgn(v)) {}
	explicit ival(double d) : ival { d, -d } {}
	constexpr ival(endpts e) : ival { e.l, -e.u } {}
	constexpr ival(cnt_rad<double,double> v) : ival { v.c - v.r, -v.c - v.r } {}

	friend constexpr double lo(const ival &v) { return  v.lo_pos; }
	friend constexpr double hi(const ival &v) { return -v.hi_neg; }

	/* always true by construction */
	friend constexpr bool   isempty(const ival &v) { return lo(v) > hi(v); }
	friend bool   isNaI(const ival &v) { return std::isnan(v.lo_pos) || std::isnan(v.hi_neg); }
	friend bool   ispoint(const ival &v) { return std::isfinite(lo(v)) && lo(v) == hi(v); }
	friend bool   isentire(const ival &v) { return !std::isfinite(lo(v)) && v.lo_pos == v.hi_neg; }
	friend bool   isbounded(const ival &v) { return std::isfinite(lo(v)) && std::isfinite(hi(v)); }

	       bool   contains(double d) const { return lo(*this) <= d && d <= hi(*this); }

	/* requires v non-empty */
	friend double inf(const ival &v) { return lo(v); }
	/* requires v non-empty */
	friend double sup(const ival &v) { return hi(v); }
	/* requires v non-empty */
	friend double mag(const ival &v) { return std::max(std::abs(lo(v)), std::abs(hi(v))); }
	/* requires v non-empty */
	friend double mig(const ival &v)
	{
		return lo(v) >= 0 ?  lo(v) : hi(v) <= 0 ? -hi(v) : 0;
	}
	/* requires v non-empty and bounded */
	friend ival   mid_enc(const ival &v) { return {  ( v.lo_pos-v.hi_neg)/2, (v.hi_neg-v.lo_pos)/2 }; }
	/* requires v non-empty */
	friend ival   wid_enc(const ival &v) { return { -(-v.lo_pos-v.hi_neg)  , (v.hi_neg+v.lo_pos)   }; }
	/* requires v non-empty and bounded */
	friend ival   rad_enc(const ival &v) { return { -(-v.lo_pos-v.hi_neg)/2, (v.hi_neg+v.lo_pos)/2 }; }

	friend double mid(const ival &v)
	{
		if (isempty(v))
			return NAN;
		if (isentire(v))
			return 0;
		if (std::isinf(lo(v)))
			return DBL_MIN;
		if (std::isinf(hi(v)))
			return DBL_MAX;
		return lo(mid_enc(v));
	}

	friend double rad(const ival &v)
	{
		if (isempty(v))
			return NAN;
		if (!isbounded(v))
			return INFINITY;
		return hi(rad_enc(v));
	}

	friend double wid(const ival &v) { return isempty(v) ? NAN : hi(wid_enc(v)); }

	friend ival   intersect(const ival &a, const ival &b)
	{
		using std::max;
		return { max(a.lo_pos, b.lo_pos), max(a.hi_neg, b.hi_neg) };
	}

	friend ival   convex_hull(const ival &a, const ival &b)
	{
		using std::min;
		return { min(a.lo_pos, b.lo_pos), min(a.hi_neg, b.hi_neg) };
	}

	friend ival   operator- (const ival &a) { return { a.hi_neg, a.lo_pos }; }

	friend void neg(ival &a) { using std::swap; swap(a.lo_pos, a.hi_neg); }

	friend ival & operator+=(ival &a, const ival &b)
	{ a.lo_pos += b.lo_pos; a.hi_neg += b.hi_neg; return a; }
	friend ival   operator+ (ival  a, const ival &b) { a += b; return a; }

	friend ival   operator-=(ival &a, const ival &b) { a += -b; return a; }
	friend ival   operator- (ival  a, const ival &b) { a -= b; return a; }

	template <typename L, typename = op_compat_t<L,float,double>>
	friend ival   operator+ (const ival &a, const L &b)
	{
		return { a.lo_pos + b, a.hi_neg - b };
	}

	template <typename L, typename = op_compat_t<L,float,double>>
	friend ival   operator* (const L &a, const ival &b)
	{
		if (a >= 0)
			return { a * b.lo_pos, a * b.hi_neg };
		else
			return { -a * b.hi_neg, -a * b.lo_pos };
	}

#if 0
	/* x*y + z */
	template <typename L, typename = op_compat_t<L,float,double>>
	friend ival   fma(const L &x, const ival &y, const ival &z)
	{
		if (x >= 0)
			return {
				fma(x, y.lo_pos, z.lo_pos),
				fma(x, y.hi_neg, z.hi_neg),
			};
		else
			return {
				fma(-x, y.hi_neg, z.lo_pos),
				fma(-x, y.lo_pos, z.hi_neg),
			};
	}

	template <typename L, typename = op_compat_t<L,float,double>>
	friend ival   fma(const L &x, const ival &y, const L &z)
	{
		if (x >= 0)
			return {
				fma(x, y.lo_pos, z),
				fma(x, y.hi_neg, z),
			};
		else
			return {
				fma(-x, y.hi_neg, z),
				fma(-x, y.lo_pos, z),
			};
	}
#endif

	friend ival   operator* (const ival &a, const ival &b)
	{
		if (lo(a) >= 0 && lo(b) >= 0) {
			/* both non-negative */
			return { a.lo_pos * b.lo_pos, -a.hi_neg * b.hi_neg };
		} else if (hi(a) <= 0 && hi(b) <= 0) {
			/* both non-positive */
			return { a.hi_neg * b.hi_neg, -a.lo_pos * b.lo_pos };
		} else if (hi(a) <= 0 && lo(b) >= 0) {
			/* a non-positive, b non-negative */
			return { a.lo_pos * -b.hi_neg, a.hi_neg * b.lo_pos };
		} else if (lo(a) >= 0 && hi(b) <= 0) {
			/* a non-negative, b non-positive */
			return { -a.hi_neg * b.lo_pos, a.lo_pos * b.hi_neg };
		} else {
			/* at least one contains zero */
			using std::min;
			return {
				min(-a.hi_neg * b.lo_pos, a.lo_pos * -b.hi_neg),
				min(-a.hi_neg * b.hi_neg, a.lo_pos * -b.lo_pos),
			};
		}
	}
	friend ival & operator*=(ival &a, const ival &b) { a = a * b; return a; }

	friend ival   operator/ (const ival &a, const ival &b)
	{
		if (b.lo_pos > 0) {
			if (a.lo_pos > 0)
				return {  a.lo_pos / -b.hi_neg,  a.hi_neg /  b.lo_pos };
			else if (a.hi_neg > 0)
				return {  a.lo_pos /  b.lo_pos,  a.hi_neg / -b.hi_neg };
			else
				return {  a.lo_pos /  b.lo_pos,  a.hi_neg /  b.lo_pos };
		} else if (b.hi_neg > 0) {
			if (a.lo_pos > 0)
				return {  a.hi_neg /  b.hi_neg, -a.lo_pos /  b.lo_pos };
			else if (a.hi_neg > 0)
				return { -a.hi_neg /  b.lo_pos,  a.lo_pos /  b.hi_neg };
			else
				return {  a.hi_neg /  b.hi_neg,  a.lo_pos /  b.hi_neg };
		} else {
			// contains zero
			return { -INFINITY, -INFINITY };
		}
	}
	friend ival & operator/=(ival &a, const ival &b) { a = a / b; return a; }

	friend ival square(const ival &i)
	{
		double lp = i.lo_pos, hn = i.hi_neg;
		using std::min;
		switch (sgn(i)) {
		case POS: return { lp * lp, -hn * hn };
		case NEG: return { hn * hn, -lp * lp };
		case ZERO: return i;
		case OV_ZERO: return { 0, min(-lp * lp, -hn * hn) };
		}
		kay_unreachable();
	}

	// TODO: use fma(3)
	friend void fma(ival &r, const ival &a, const ival &b) { r += a * b; }

	friend ival_pos cmp_detailed(const ival &a, const ival &b)
	{
		int ll = cmp(lo(a), lo(b));
		int hl = cmp(hi(a), lo(b));
		int lh = cmp(lo(a), hi(b));
		int hh = cmp(hi(a), hi(b));

		if (hl < 0) return IVAL_LT;
		if (ll < 0 && hh < 0)
			return !hl ? IVAL_LE : IVAL_LO;
		if (lh > 0) return IVAL_GT;
		if (ll > 0 && hh > 0)
			return !lh ? IVAL_GE : IVAL_GO;
		if (ll == hh)
			return IVAL_EQ;
		return ll > hh ? IVAL_SUB : IVAL_SUP;
	}

	/* -1 if all points in a are smaller than points in b
	 *  0 if a and b share at least one point
	 * +1 if all points in a are larger than points in b */
	friend int cmp(const ival &a, const ival &b)
	{
		if (hi(a) < lo(b))
			return -1;
		if (lo(a) > hi(b))
			return +1;
		return 0;
	}

	friend ival_sgn sgn(const ival &a)
	{
		if (lo(a) > 0)
			return POS;
		if (hi(a) < 0)
			return NEG;
		if (ispoint(a))
			return ZERO;
		return OV_ZERO;
	}

	friend ival max(const ival &a, double b)
	{
		using std::max;
		return endpts { max(lo(a), b), max(hi(a), b) };
	}

	friend ival min(const ival &a, double b)
	{
		using std::min;
		return endpts { min(lo(a), b), min(hi(a), b) };
	}

	friend ival tanh(const ival &a)
	{
		/* tanh(x) is monotonic */
		return endpts { tanh(lo(a)), tanh(hi(a)) };
	}

	friend bool issubset(const ival &a, const ival &b)
	{
		return lo(a) >= lo(b) && hi(a) <= hi(b);
	}

	friend std::ostream & operator<<(std::ostream &os, const ival &a)
	{
		if (isempty(a))
			os << "[]";
		else if (ispoint(a))
			os << "[" << lo(a) << "]";
		else {
			if (std::isinf(lo(a)))
				os << "(-infty";
			else
				os << "[" << lo(a);
			os << ",";
			if (std::isinf(hi(a)))
				os << "infty)";
			else
				os << hi(a) << "]";
		}
		return os;
	}
};

}

#endif
