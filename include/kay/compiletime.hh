/*
 * compiletime.hh
 *
 * Copyright 2019 Franz Brauße <brausse@informatik.uni-trier.de>
 *
 * This file is part of kay.
 * See the LICENSE file for terms of distribution.
 */

#ifndef KAY_COMPILETIME_HH
#define KAY_COMPILETIME_HH

namespace kay::compiletime {

/* contains strings termed 'L' (using namespace limbs) and arbitrary sized
 * natural numbers (using namespace N) */

/* --------------------------------------------------------------------------
 * strings over alphabet uint64_t
 * -------------------------------------------------------------------------- */
namespace limbs {

/* L is a string of limbs */
template <uint64_t... Limbs> using L = std::integer_sequence<uint64_t,Limbs...>;

namespace infix {

/* equality of strings */
template <uint64_t... As> constexpr bool operator==(L<As...>,L<As...>) { return true; }
template <uint64_t... As> constexpr bool operator!=(L<As...>,L<As...>) { return false; }

}

/* concat1: prepend one symbol to a string */
template <uint64_t A, typename B> struct concat1;
template <uint64_t A, typename B> using concat1_t = typename concat1<A,B>::type;
template <uint64_t A, uint64_t... Bs> struct concat1<A,L<Bs...>> {
	using type = L<A,Bs...>;
};

/* concat: concatenate two strings */
template <typename A, typename B> struct concat;
template <typename A, typename B> using concat_t = typename concat<A,B>::type;
template <typename B, uint64_t a, uint64_t... As> struct concat<L<a,As...>,B> {
	using type = concat1_t<a,concat_t<L<As...>,B>>;
};
template <typename B> struct concat<L<>,B> {
	using type = B;
};

static_assert(std::is_same_v<concat_t<L<0,1>,L<>>,L<0,1>>);
static_assert(std::is_same_v<concat_t<L<0,1>,L<2,3>>,L<0,1,2,3>>);

/* length of a string */
template <typename A> constexpr size_t length_v = A::size();

/* ileave2: interleave two strings A0 A1... and B0 B1... limb by limb:
 *          A0 B0 A1 B1 ...; extending the shorter one by '0' limbs */
template <typename A, typename B> struct ileave2;
template <typename A, typename B> using ileave2_t = typename ileave2<A,B>::type;

/* ileave1: helper to ileave2 */
template <typename A, uint64_t B, typename BB> struct ileave1;
template <typename A, uint64_t B, typename BB> using ileave1_t = typename ileave1<A,B,BB>::type;
template <uint64_t B, typename BB, uint64_t A, uint64_t... As> struct ileave1<L<A,As...>,B,BB> {
	using type = concat1_t<A,concat1_t<B,ileave2_t<L<As...>,BB>>>;
};
template <uint64_t B, typename BB> struct ileave1<L<>,B,BB> {
	using type = concat1_t<0,concat1_t<B,ileave2_t<L<>,BB>>>;
};

template <typename A, uint64_t B, uint64_t... Bs> struct ileave2<A,L<B,Bs...>> {
	using type = ileave1_t<A,B,L<Bs...>>;
};
template <typename A> struct ileave2<A,L<>> {
	using type = ileave1_t<A,0,L<>>;
};
template <> struct ileave2<L<>,L<>> {
	using type = L<>;
};

}

/* --------------------------------------------------------------------------
 * notation of natural numbers using limbs (strings over alphabet size_t)
 * -------------------------------------------------------------------------- */
namespace N {

/* exports instances of L<uint64_t...>:
 *
 * - N<uint64_t...>
 * - add_t<L<uint64_t...>,L<uint64_t...>>
 * - mul_t<L<uint64_t...>,L<uint64_t...>>
 *
 * and cmp, cmp::{LT,EQ,GT}, cmp_v<L<uint64_t...>,L<uint64_t...>>
 * and 'namespace infix' providing constexpr versions of
 *
 * - operator==(L<uint64_t...>,L<uint64_t...>)
 * - operator!=(L<uint64_t...>,L<uint64_t...>)
 * - operator< (L<uint64_t...>,L<uint64_t...>)
 * - operator<=(L<uint64_t...>,L<uint64_t...>)
 * - operator> (L<uint64_t...>,L<uint64_t...>)
 * - operator>=(L<uint64_t...>,L<uint64_t...>)
 * - operator+(L<uint64_t...>,L<uint64_t...>)
 * - operator*(L<uint64_t...>,L<uint64_t...>)
 */

namespace detail {

using namespace limbs;

template <uint64_t n, typename A> struct strip_leading_0_rev;
template <uint64_t n, typename A> using strip_leading_0_rev_t = typename strip_leading_0_rev<n,A>::type;
template <> struct strip_leading_0_rev<0,L<>> {
	using type = L<>;
};
template <uint64_t n, uint64_t... As> struct strip_leading_0_rev<n,L<As...>> {
	using type = concat1_t<n,L<As...>>;
};

template <typename A> struct strip_leading_0;
template <typename A> using strip_leading_0_t = typename strip_leading_0<A>::type;
template <uint64_t n, uint64_t... As> struct strip_leading_0<L<n,As...>> {
	using type = strip_leading_0_rev_t<n,strip_leading_0_t<L<As...>>>;
};
template <> struct strip_leading_0<L<>> {
	using type = L<>;
};

}

/* N: a natural number is a string of limbs with leading zeroes stripped */
template <uint64_t... As> using N = detail::strip_leading_0_t<limbs::L<As...>>;

namespace detail {

static_assert(std::is_same_v<N<0,0>,L<>>);
static_assert(std::is_same_v<N<1,0,0>,L<1>>);

/* iadd: addition of two interleaved natural numbers */
template <bool carry, typename A> struct iadd;
template <bool carry, typename A> using iadd_t = typename iadd<carry,A>::type;
template <bool carry, uint64_t A, uint64_t B, uint64_t... Cs> struct iadd<carry,L<A,B,Cs...>> {
	using type = concat1_t<A+B+carry,iadd_t<(A+B<A || A+B+carry<A+B),L<Cs...>>>;
};
template <uint64_t C> struct iadd<false,L<C>> {
	using type = L<C>;
};
template <> struct iadd<false,L<>> {
	using type = L<>;
};
template <> struct iadd<true,L<>> {
	using type = L<1>;
};

}

/* add: addition of two natural numbers represented by strings */
template <typename A, typename B> using add_t = detail::iadd_t<false,limbs::ileave2_t<A,B>>;

static_assert(std::is_same_v<add_t<N<>                     ,N<>              >,N<>          >);
static_assert(std::is_same_v<add_t<N<0,1>                  ,N<3>             >,N<3,1>       >);
static_assert(std::is_same_v<add_t<N<UINT64_MAX>           ,N<0,0>           >,N<UINT64_MAX>>);
static_assert(std::is_same_v<add_t<N<UINT64_MAX>           ,N<1,UINT64_MAX>  >,N<0,0,1>     >);
static_assert(std::is_same_v<add_t<N<UINT64_MAX,UINT64_MAX>,N<1>             >,N<0,0,1>     >);
static_assert(std::is_same_v<add_t<N<1,UINT64_MAX>         ,N<UINT64_MAX>    >,N<0,0,1>     >);
static_assert(std::is_same_v<add_t<N<1,0,1>                ,N<UINT64_MAX,0,1>>,N<0,1,2>     >);
static_assert(std::is_same_v<add_t<N<1>                    ,N<2>             >,N<3>         >);

namespace infix {

using limbs::infix::operator==;

/* infix addition */
template <typename A, typename B> constexpr add_t<A,B> operator+(A,B) { return {}; }

static_assert(N<1>{} + N<2>{} == N<3,0,0,0,0>{});

}


namespace detail {

/* mul1: multiply limb by natural number */
template <size_t a, typename B> struct imul1;
template <size_t a, typename B> using imul1_t = typename imul1<a,B>::type;
template <uint64_t a> struct imul1<a,L<>> {
	using type = N<>;
};
template <uint64_t a, uint64_t b, size_t... Bs> struct imul1<a,L<b,Bs...>> {
private:
	static constexpr uint64_t lo(uint64_t v) { return v & 0xffff'ffff; }
	static constexpr uint64_t hi(uint64_t v) { return v >> 32; }

	/* largest value (B-1)^2 = B^2-2B+1 for B=2^32 */
	static constexpr uint64_t ll = lo(a) * lo(b);
	static constexpr uint64_t lh = lo(a) * hi(b);
	static constexpr uint64_t hl = hi(a) * lo(b);
	static constexpr uint64_t hh = hi(a) * hi(b);

	/* no overflows in intermediate sums */
	static constexpr uint64_t x = hl + hi(ll); /* B^2-2B+1 + (B-1) = B^2-B */
	static constexpr uint64_t y = lo(x) + lh;  /* (B-1) + B^2-2B+1 = B^2-B */
	static constexpr uint64_t c = hi(x) + hi(y) + hh; /* 2(B-1) + B^2-2B+1 = B^2-1 */
	static constexpr uint64_t r = lo(ll) | lo(y) << 32;

public:
	using type = concat1_t<r,add_t<N<c>,imul1_t<a,L<Bs...>>>>;
};

static_assert(std::is_same_v<imul1_t<2,N<1>>,N<2>>);
static_assert(std::is_same_v<imul1_t<UINT64_MAX,N<2>>,N<(UINT64_MAX << 1),1>>);
static_assert(std::is_same_v<imul1_t<2,N<UINT64_MAX>>,N<(UINT64_MAX << 1),1>>);
static_assert(std::is_same_v<imul1_t<4,N<UINT64_MAX>>,N<(UINT64_MAX << 2),3>>);
static_assert(std::is_same_v<imul1_t<(1ULL<<32),N<UINT64_MAX>>,N<(UINT64_MAX << 32),(UINT64_MAX >> 32)>>);
static_assert(std::is_same_v<imul1_t<(1ULL<<63),N<UINT64_MAX>>,N<(UINT64_MAX << 63),(UINT64_MAX >> 1)>>);
static_assert(std::is_same_v<imul1_t<UINT64_MAX,N<UINT64_MAX>>,N<1,UINT64_MAX-1>>);
static_assert(std::is_same_v<imul1_t<UINT64_MAX,N<UINT64_MAX,UINT64_MAX>>,N<1,UINT64_MAX,UINT64_MAX-1>>);

/* mul: multiply natural numbers (schoolbook) */
template <typename A, typename B> struct mul;
template <typename A, typename B> using mul_t = typename mul<A,B>::type;
template <typename B> struct mul<L<>,B> {
	using type = N<>;
};
template <size_t a, size_t... As, typename B> struct mul<L<a,As...>,B> {
	using type = add_t<imul1_t<a,B>,concat1_t<0,mul_t<L<As...>,B>>>;
};

}

using detail::mul_t;

static_assert(std::is_same_v<mul_t<N<UINT64_MAX>,N<UINT64_MAX>>,N<1,UINT64_MAX-1>>);
static_assert(std::is_same_v<mul_t<N<UINT64_MAX,UINT64_MAX>,N<UINT64_MAX>>,N<1,UINT64_MAX,UINT64_MAX-1>>);
static_assert(std::is_same_v<mul_t<N<UINT64_MAX,UINT64_MAX>,N<UINT64_MAX,UINT64_MAX>>,N<1,0,UINT64_MAX-1,UINT64_MAX>>);

namespace infix {

/* infix multiplication */
template <typename A, typename B> constexpr mul_t<A,B> operator*(A,B) { return {}; }

static_assert(N<2>{} * N<2>{} * N<2>{} == N<8,0,0,0,0>{});

}

/* cmp cmp_v: compare natural numbers */
enum class cmp { LT = -1, EQ = 0, GT = +1 };

namespace detail {

template <uint64_t a, uint64_t b, cmp c> struct icmp_merge
: std::integral_constant<cmp,c> {};
template <uint64_t a, uint64_t b> struct icmp_merge<a,b,cmp::EQ> {
	static constexpr cmp value = a < b ? cmp::LT : a > b ? cmp::GT : cmp::EQ;
};

template <typename A> struct icmp;
template <typename A> static constexpr cmp icmp_v = icmp<A>::value;
template <> struct icmp<L<>> : std::integral_constant<cmp,cmp::EQ> {};
template <uint64_t a, uint64_t b, uint64_t... Cs> struct icmp<L<a,b,Cs...>>
: icmp_merge<a,b,icmp_v<L<Cs...>>> {};

}

template <typename A, typename B> static constexpr cmp cmp_v = detail::icmp_v<limbs::ileave2_t<A,B>>;

static_assert(cmp_v<N<1,2>,N<2,1>> == cmp::GT);
static_assert(cmp_v<N<>,N<>> == cmp::EQ);
static_assert(cmp_v<N<>,N<2>> == cmp::LT);
static_assert(cmp_v<N<1>,N<1,1>> == cmp::LT);

namespace infix {

template <typename A, typename B> constexpr bool operator< (A,B) { return cmp_v<A,B> == cmp::LT; }
template <typename A, typename B> constexpr bool operator<=(A,B) { return cmp_v<A,B> != cmp::GT; }
template <typename A, typename B> constexpr bool operator> (A,B) { return cmp_v<A,B> == cmp::GT; }
template <typename A, typename B> constexpr bool operator>=(A,B) { return cmp_v<A,B> != cmp::LT; }

}

} /* namespace N */

} /* namespace kay::compiletime */

#endif
