#pragma once

#include <cassert>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <sstream>

typedef std::pair<NTL::ZZ_p, NTL::ZZ_p> fp2_elem;

//// Printing
std::string StringFp2(fp2_elem const &a) {
    std::ostringstream oss;
    oss << a.first << " + i*" << a.second;
    return oss.str();
}

fp2_elem Fp2_zero() {
    return {NTL::ZZ_p(0), NTL::ZZ_p(0)};
}

fp2_elem Fp2_one() {
    return {NTL::ZZ_p(1), NTL::ZZ_p(0)};
}

fp2_elem Fp2_i() {
    return {NTL::ZZ_p(0), NTL::ZZ_p(1)};
}

fp2_elem Fp2_add(fp2_elem const &a, fp2_elem const &b) {
    return {a.first + b.first, a.second + b.second};
}

fp2_elem Fp2_sub(fp2_elem const &a, fp2_elem const &b) {
    return {a.first - b.first, a.second - b.second};
}

fp2_elem Fp2_mul(fp2_elem const &a, fp2_elem const &b) {
    NTL::ZZ_p t0 = a.first * b.first;
    NTL::ZZ_p t1 = a.second * b.second;
    NTL::ZZ_p t2 = a.first + a.second;
    NTL::ZZ_p t3 = b.first + b.second;

    NTL::ZZ_p res0 = t0 - t1;
    t0 = t0 + t1;
    t1 = t2 * t3;

    return {res0, t1 - t0};
}

fp2_elem Fp2_sqr(fp2_elem const &a) {
    NTL::ZZ_p t0 = a.first + a.second;
    NTL::ZZ_p t1 = a.first - a.second;
    NTL::ZZ_p t2 = a.first + a.first;

    return {t0*t1, t2*a.second};
}

fp2_elem Fp2_inv(fp2_elem const &a) {
    NTL::ZZ_p t0 = a.first * a.first;
    NTL::ZZ_p t1 = a.second * a.second;
    t0 = 1/(t0 + t1);

    return {a.first * t0, -(a.second * t0)};
}

fp2_elem Fp2_pow(fp2_elem a, int e) {
    fp2_elem res = Fp2_one();
    while (e > 0) {
        if (e & 1) {
            res = Fp2_mul(res, a);
        }
        a = Fp2_sqr(a);
        e>>= 1;
    }
    return res;
}

void Fp2_neg(fp2_elem &a) {
    //negate in place
    a.first = -a.first;
    a.second = -a.second;
}

fp2_elem Fp2_negative(fp2_elem const &a) {
    //negate with return
    return {-a.first,-a.second};
}

fp2_elem Fp2_div(fp2_elem const &a, fp2_elem const &b) {
    return Fp2_mul(a, Fp2_inv(b));
}

bool Fp2_issquare(fp2_elem const &a) {
    NTL::ZZ_p t0 = a.first * a.first;
	NTL::ZZ_p t1 = a.second * a.second;
	t0 = t0 + t1;

    NTL::ZZ p = NTL::ZZ_p::modulus();
    NTL::ZZ e = (p - 1)/2;

    return NTL::PowerMod(rep(t0), e, p) == 1;
}

fp2_elem Fp2_sqroot(fp2_elem const &a) {

    assert (Fp2_issquare(a));

    NTL::ZZ p = NTL::ZZ_p::modulus();
	if (a.second == 0) {
		if (NTL::Jacobi(rep(a.first), p) == 1) {
            NTL::ZZ_p b;
            NTL::conv(b, NTL::SqrRootMod(rep(a.first), p));
            return {b, NTL::ZZ_p(0)};
        } else {
            NTL::ZZ_p b;
            NTL::conv(b, NTL::SqrRootMod(rep(-a.first), p));
            return {NTL::ZZ_p(0), b};
        }
    }

    NTL::ZZ_p sdelta = a.first*a.first;
	NTL::ZZ_p t1 = a.second*a.second;
	sdelta = sdelta + t1;
	NTL::conv(sdelta, NTL::SqrRootMod(rep(sdelta), p));

    NTL::ZZ_p inv2(1);
    inv2 /= NTL::ZZ_p(2);

	NTL::ZZ_p re = a.first + sdelta;
	re *= inv2;

	if (NTL::Jacobi(rep(re), p) != 1) {
        re = a.first - sdelta;
        re *= inv2;
    }

    NTL::conv(re, NTL::SqrRootMod(rep(re), p));
    NTL::ZZ_p im(1);
    im /= re;
    im *= inv2;
	im *= a.second;

	return {re, im};
}

std::pair<fp2_elem, fp2_elem> QuadraticRoot(fp2_elem const &b, fp2_elem const &c) {
    // Find a root of the polynomial x^2 + bx + c
    // (-b \pm \sqrt(b^2 - 4c))/2
    fp2_elem disc = Fp2_sqr(b);
    fp2_elem t0 = Fp2_add(c, c);
    t0 = Fp2_add(t0, t0);
    disc = Fp2_sub(disc, t0);

    disc = Fp2_sqroot(disc);

    fp2_elem r1 = Fp2_sub(disc, b);
    Fp2_neg(disc);
    fp2_elem r2 = Fp2_sub(disc, b);

    fp2_elem inv2 = {NTL::ZZ_p(2), NTL::ZZ_p(0)};
    inv2 = Fp2_inv(inv2);

    r1 = Fp2_mul(r1, inv2);
    r2 = Fp2_mul(r2, inv2);

    return {r1, r2};
}

fp2_elem Fp2_randomElement() {
    return {NTL::random_ZZ_p(), NTL::random_ZZ_p()};
}

bool IsZero(fp2_elem const &a) {
    return (a.first == 0) && (a.second == 0);
}

bool IsOne(fp2_elem const &a) {
    return (a.first == 1) && (a.second == 0);
}

bool Fp2_equal(fp2_elem const &a, fp2_elem const &b) {
    return IsZero(Fp2_sub(a, b));
}