#pragma once

#include <cassert>
#include <array>
#include <string>
#include <sstream>
#include <functional>
#include <NTL/ZZ.h>
#include <NTL/HNF.h>

class quatlat;

struct quatalg
{
    NTL::ZZ p, q;

    // returns matrices of i and j modulo ell
    std::pair<NTL::mat_ZZ_p,NTL::mat_ZZ_p> splitting(NTL::ZZ const &ell) const
    {
        if (ell == 2)
            throw;  // not implemented
        if (q % ell == 0)
            throw;  // not implemented
        assert(NTL::ZZ_p::modulus() == ell);

        NTL::ZZ_p pmod, qmod, qinv, b, a;
        NTL::conv(pmod, p);
        NTL::conv(qmod, q);
        NTL::inv(qinv, qmod);
        NTL::set(b);
        for (size_t i = 0; i < 999; ++i, ++b) {
            auto a2 = -pmod - b*b*qinv;
            if (NTL::Jacobi(NTL::rep(a2), ell) == 1) {
                NTL::conv(a, NTL::SqrRootMod(NTL::rep(a2), ell));
                assert(a*a == a2);
                break;
            }
        }
        assert(!NTL::IsZero(a));

        std::pair<NTL::mat_ZZ_p,NTL::mat_ZZ_p> ret;
        ret.first.SetDims(2, 2);
        ret.second.SetDims(2, 2);
        ret.first[0][1] = -qmod;
        NTL::set(ret.first[1][0]);
        ret.second[0][0] = a;
        ret.second[0][1] = b;
        ret.second[1][0] = (-pmod-a*a)/b;
        ret.second[1][1] = -a;
        assert(NTL::IsIdent(ret.first*ret.first*NTL::inv(-qmod), 2));
        assert(NTL::IsIdent(ret.second*ret.second*NTL::inv(-pmod), 2));
        assert(ret.first*ret.second == -ret.second*ret.first);
        return ret;
    }    
    
    quatlat maximal_order() const; // return a standard choice of maximal order
};

class quat
{
private:
    std::array<NTL::ZZ,5> coeffs; //the four coefficients, and the denominator

public:
    quatalg const &alg;

    quat(quatalg const &alg_) : alg{alg_} { NTL::set(coeffs[4]); };
    quat(std::array<NTL::ZZ,5> const &coeffs_, quatalg const &alg_)
        : coeffs{coeffs_}, alg{alg_}
    {
        normalize();
    }

    void normalize()
    {
        assert(!NTL::IsZero(coeffs[4]));
        if (coeffs[4] < 0)
            for (auto &c: coeffs)
                c = -c;
        NTL::ZZ g = coeffs[4];
        for (unsigned i = 0; i < 4 && !NTL::IsOne(g); ++i)
            g = NTL::GCD(g, coeffs[i]);
        assert(g > 0);
        if (!NTL::IsOne(g))
            for (auto &c: coeffs)
                c /= g;
    }

    quat const conjugate() const
    {
        return {{coeffs[0], -coeffs[1], -coeffs[2], -coeffs[3], coeffs[4]}, alg};
    }

    // in-place
    void invert()
    {
        auto conj = conjugate();
        auto nrm = norm();
        coeffs[0] = conj[0]*nrm.second;
        coeffs[1] = conj[1]*nrm.second;
        coeffs[2] = conj[2]*nrm.second;
        coeffs[3] = conj[3]*nrm.second;
        coeffs[4] = conj[4]*nrm.first;
        normalize();
    }

    NTL::ZZ &operator[](size_t index)
    {
        if (index > 4)
            throw std::out_of_range("Index out of bounds");
        return coeffs[index];
    }

    NTL::ZZ const &operator[](size_t index) const
    {
        if (index > 4)
            throw std::out_of_range("Index out of bounds");
        return coeffs[index];
    }

    quat operator+(quat const &other) const
    {
        if (&other.alg != &alg) // object identity
            throw;
        auto const &[a,b,c,d,e] = coeffs;
        auto const &[A,B,C,D,E] = other.coeffs;
        return {{a*E+A*e, b*E+B*e, c*E+C*e, d*E+D*e, e*E}, alg};
    };

    quat &operator+=(quat const &other)
    {
        auto const sum = *this + other;
        std::copy(sum.coeffs.begin(), sum.coeffs.end(), this->coeffs.begin());
        return *this;
    };

    quat operator*(quat const &other) const
    {
        if (&other.alg != &alg) // object identity
            throw;
        auto const &[a,b,c,d,e] = this->coeffs;
        auto const &[A,B,C,D,E] = other.coeffs;
        auto const r = a*A - alg.q*b*B - alg.p*c*C - alg.p*alg.q*d*D;
        auto const s = a*B + A*b + alg.p*c*D - alg.p*C*d;
        auto const t = a*C + A*c - alg.q*b*D + alg.q*B*d;
        auto const u = a*D + A*d + b*C - B*c;
        return {{r,s,t,u, e*E}, alg};
    };

    quat operator*(NTL::ZZ const &n) const { return {{n*coeffs[0], n*coeffs[1], n*coeffs[2], n*coeffs[3], coeffs[4]}, alg}; }
    quat operator*(long n) const { return {{n*coeffs[0], n*coeffs[1], n*coeffs[2], n*coeffs[3], coeffs[4]}, alg}; }

    std::pair<NTL::ZZ,NTL::ZZ> norm() const  // returns numerator and denominator
    {
        auto n = (*this) * this->conjugate();
        return {n[0], n[4]};
    }

    std::pair<NTL::ZZ,NTL::ZZ> trace() const  // returns numerator and denominator
    {
        auto num = 2*coeffs[0], denom = coeffs[4];
        auto g = NTL::GCD(num, denom);
        if (!NTL::IsOne(g)) {
            num /= g;
            denom /= g;
        }
        return {num, denom};
    }

    std::pair<NTL::ZZ,NTL::ZZ> trace_pairing(quat const &other) const
    {
        if (&other.alg != &alg) // object identity
            throw;
        return (this->conjugate() * other).trace();
    }

    friend std::ostream& operator<<(std::ostream& o, quat const &alpha)
    {
            return o
            << alpha[0] << "/" << alpha[4] << " + "
            << alpha[1] << "/" << alpha[4] << "*i + "
            << alpha[2] << "/" << alpha[4] << "*j + "
            << alpha[3] << "/" << alpha[4] << "*k";
    }

    bool is_one() const
    {
        auto const &[a,b,c,d,e] = coeffs;
        assert(e != 0);
        return a == e && b == 0 && c == 0 && d == 0;
    }
};


#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

class quatlat
{
public:
    quatalg const &alg;

    NTL::mat_ZZ basis;
    NTL::ZZ denom;

private:
    mutable std::pair<NTL::ZZ,NTL::ZZ> the_norm {};

public:
    quatlat(NTL::mat_ZZ const &gens, NTL::ZZ const &denom_, quatalg const &alg_) : alg{alg_}, basis{gens}, denom{denom_}
    {
        assert(basis.NumCols() == 4);

        normalize();
    }

    quatlat(quatlat const &L) : alg{L.alg} { *this = L; }

    void normalize()
    {
        NTL::ZZ det;

        size_t numrows = basis.NumRows();
        size_t rank = NTL::LLL(det, basis);
        if (rank != numrows)
            for (size_t i = 0; i < rank; ++i)
                std::swap(basis[i], basis[numrows-rank+i]);
        basis.SetDims(rank, 4);

        assert(basis.NumRows() == 4);
        assert(basis.NumCols() == 4);

        assert(denom != 0);
        if (denom < 0) {
            basis = -basis;
            denom = -denom;
        }

        NTL::ZZ g = denom;
        for (size_t i = 0; i < (size_t) basis.NumRows(); ++i) {
            if (NTL::IsOne(g)) break;
            for (size_t j = 0; j < (size_t) basis.NumCols(); ++j) {
                if (NTL::IsOne(g)) break;
                g = NTL::GCD(g, basis[i][j]);
            }
        }

        assert(g > 0);
        if (!NTL::IsOne(g)) {
            denom /= g;
            for (size_t i = 0; i < (size_t) basis.NumRows(); ++i)
                for (size_t j = 0; j < (size_t) basis.NumCols(); ++j)
                    basis[i][j] /= g;
        }

        the_norm.first = the_norm.second = 0;
    }

    quatlat conjugate() const
    {
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(4, 4);
        for (unsigned i = 0; i < 4; ++i) {
            newbasis[i][0] =  basis[i][0];
            newbasis[i][1] = -basis[i][1];
            newbasis[i][2] = -basis[i][2];
            newbasis[i][3] = -basis[i][3];
        }
        return {newbasis, denom, alg};
    }

    quatlat inverse() const
    {
        auto conj = conjugate();
        auto nrm = norm();
        return {conj.basis*nrm.second, denom*nrm.first, alg};
    }

    quatlat operator+(quatlat const &other) const
    {
        size_t n = this->basis.NumRows();
        size_t m = other.basis.NumRows();
        NTL::mat_ZZ newbasis = this->basis * other.denom;
        newbasis.SetDims(n+m, 4);
        for (size_t i = 0; i < m; ++i)
            newbasis[n+i] = other.basis[i] * this->denom;
        auto newdenom = this->denom * other.denom;
        return {newbasis, newdenom, alg};
    }

    quatlat operator*(quatlat const &other) const
    {
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(16, 4);
        quat a{alg}, b{alg};
        a[4] = this->denom;
        b[4] = other.denom;
        auto newdenom = this->denom * other.denom;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t k = 0; k < 4; ++k)
                a[k] = this->basis[i][k];
            for (size_t j = 0; j < 4; ++j) {
                for (size_t k = 0; k < 4; ++k)
                    b[k] = other.basis[j][k];

                auto c = a * b;
//std::cerr << "(" << a << ") * (" << b << ") = " << c << std::endl;
                for (unsigned k = 0; k < 4; ++k)
                    newbasis[4*i+j][k] = c[k] * newdenom / c[4];
            }
        }
        return {newbasis, newdenom, alg};
    }

    quatlat operator*(quat const &other) const
    {
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(4, 4);
        quat a{alg};
        a[4] = this->denom;
        auto newdenom = this->denom * other[4];
        for (size_t i = 0; i < 4; ++i) {
            for (size_t k = 0; k < 4; ++k) {
                a[k] = this->basis[i][k];

                auto c = a * other;
                for (unsigned k = 0; k < 4; ++k)
                    newbasis[i][k] = c[k] * newdenom / c[4];
            }
        }
        return {newbasis, newdenom, alg};
    }

    quatlat operator*(NTL::ZZ const &other) const
    {
        NTL::mat_ZZ newbasis = basis;
        newbasis *= other;
        return {newbasis, denom, alg};
    }

    friend quatlat operator*(quat const &left, quatlat const &right)
    { return (right.conjugate()*left.conjugate()).conjugate(); }  //TODO write faster version

    friend quatlat operator*(NTL::ZZ const &left, quatlat const &right)
    { return right * left; }

    quatlat& operator=(const quatlat& other) {
        if (this != &other) {  // Avoid self-assignment
            // Copy the values
            if (&this->alg != &other.alg) {
                throw std::runtime_error("LOUD ERROR WHEN ASSIGNING QUATLAT");
            }
            basis = other.basis;
            denom = other.denom;
            the_norm = other.the_norm;
        }
        return *this;
    }

    std::pair<NTL::ZZ,NTL::ZZ> norm() const  // returns numerator and denominator
    {
        if (NTL::IsZero(the_norm.second)) {
            auto L = (*this) * this->conjugate();
            NTL::vec_ZZ v, x;
            v.SetLength(4);
            NTL::set(v[0]);
            NTL::ZZ num;
            NTL::solve1(num, x, L.basis, v);
            auto g = NTL::GCD(num, L.denom);
            the_norm = std::make_pair(num/g, L.denom/g);
        }
        return the_norm;
    }

    // in-place
    void _intersect(quatlat const &other)
    {
        auto bas1 = NTL::transpose(basis * other.denom);
        auto bas2 = NTL::transpose(other.basis * denom);
        NTL::ZZ det1, det2;
        NTL::mat_ZZ ker1, ker2;
        NTL::inv(det1, ker1, bas1);
        NTL::inv(det2, ker2, bas2);
        assert(ker1.NumRows() == 4 && ker2.NumRows() == 4);

        NTL::mat_ZZ ker;
        ker.SetDims(8, 4);
        for (unsigned j = 0, i = 0; j < 4; ++j) {
            for (i = 0; i < ker1.NumRows(); ++i)
                ker[i][j] = ker1[i][j] * det2;
            for (i = 0; i < ker2.NumRows(); ++i)
                ker[4+i][j] = ker2[i][j] * det1;
        }

        NTL::ZZ det;
        size_t rank = NTL::image(det, ker);
        assert(rank == 4);
        for (size_t i = 0; i < rank; ++i)
            std::swap(ker[i], ker[4+i]);
        ker.SetDims(4, 4);

        NTL::transpose(ker, ker);

        NTL::inv(det, basis, ker);
        denom *= other.denom * det;
        assert(denom % (det1 * det2) == 0);
        denom /= det1 * det2;

        normalize();
    }

    // set of all x such that xL <= L or Lx <= L
    quatlat _compute_order(bool right_order=false) const
    {
        auto nrm = norm();
        quatlat L(basis*nrm.second, denom*nrm.first, alg);

        for (unsigned k = 0; k < 4; ++k) {
            quat a{alg};
            for (unsigned j = 0; j < 4; ++j)
                a[j] = basis[k][j];
            a[4] = denom;
            a.invert();

            NTL::mat_ZZ mat;
            mat.SetDims(4, 4);
            NTL::ZZ d(1);
            for (unsigned i = 0; i < 4; ++i) {
                quat b{alg};
                for (unsigned j = 0; j < 4; ++j)
                    b[j] = basis[i][j];
                b[4] = denom;
                auto c = right_order ? a * b : b * a;
                mat *= c[4];
                for (unsigned j = 0; j < 4; ++j)
                    mat[i][j] = c[j] * d;
                d *= c[4];
            }

            quatlat N(mat, d, alg);
//std::cerr << "N:\n" << N << "\n";
            L._intersect(N);
//std::cerr << "L:\n" << L << "\n";
        }

        assert(L.is_order());
        return L;
    }
    quatlat left_order() const { return _compute_order(false); }
    quatlat right_order() const { return _compute_order(true); }

    bool is_order() const { auto [n,d] = norm(); assert(!NTL::IsZero(d)); return n == d; }

    // https://math.dartmouth.edu/~jvoight/articles/73446.pdf Lemma 7.2
    void right_ideals_of_norm(NTL::ZZ const &ell, std::function<void(quatlat const &)> const &fun)
    {
        if (!is_order())
            throw std::logic_error("not an order");

        NTL::ZZ_pPush push(ell);

        auto mat1 = NTL::ident_mat_ZZ_p(2);
        auto [mati, matj] = alg.splitting(ell);
        auto matk = mati * matj;
        NTL::mat_ZZ_p const *mats[4] = {&mat1, &mati, &matj, &matk};

        // can ignore denominator for now because it's all projective anyway

        NTL::mat_ZZ mat;
        mat.SetDims(4, 4);
        for (unsigned k = 0; k < 4; ++k) {
            NTL::mat_ZZ cur;
            cur.SetDims(2,2);
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j)
                    cur[j/2][j%2] += basis[k][i] * NTL::rep((*mats[i])[j/2][j%2]);
//            std::cerr << cur << std::endl;
            for (unsigned j = 0; j < 4; ++j)
                mat[k][j] = cur[j/2][j%2];
//            std::cerr << mat[k] << std::endl;
        }

        // we do need to be careful when the denominator is divisible by ell
        //XXX currently it fails in this case; TODO fix
//std::cerr << NTL::determinant(mat) << std::endl;
//std::cerr << mat << std::endl;

        NTL::mat_ZZ_p matmod, matinv;
        matmod.SetDims(4,4);
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                NTL::conv(matmod[i][j], mat[i][j]);
        NTL::inv(matinv, matmod);

        NTL::vec_ZZ_p rhs;
        rhs.SetLength(4);
        NTL::ZZ_p &x = rhs[0], &y = rhs[2];

        NTL::set(x);
        for (NTL::ZZ i; i <= ell; ++i) {

            auto sol = rhs * matinv;
            assert(sol * matmod == rhs);

            quat elt{alg};
            NTL::set(elt[4]);
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j)
                    elt[j] += NTL::rep(sol[i]) * basis[i][j];
            assert(NTL::IsOne(elt.norm().second));
            assert(NTL::IsZero(elt.norm().first % ell));

//std::cerr << rhs << std::endl;

            auto I = ell*(*this) + elt*(*this);
            assert(I.norm().first == ell && NTL::IsOne(I.norm().second));

            fun(I);

            if (NTL::IsZero(i)) {
                NTL::clear(x);
                NTL::set(y);
            }
            else
                ++x;
        }
    }

    void left_ideals_of_norm(NTL::ZZ const &ell, std::function<void(quatlat const &)> const &fun)
    {
        auto const fun2 = [&](quatlat const &lat) {
            fun(lat.conjugate());
        };
        right_ideals_of_norm(ell, fun2);
    }

    friend std::ostream& operator<<(std::ostream& o, quatlat const &L)
    {
        return o << L.basis << "/" << L.denom << "  (norm " << L.norm().first << "/" << L.norm().second << ")";
    }

    std::string sage() const
    {
        std::ostringstream s;
        s << "i.parent().ideal([";
        for (unsigned i = 0; i < 4; ++i) {
            s << "[";
            for (unsigned j = 0; j < 4; ++j)
                s << basis[i][j] << "/" << denom << ",";
            s << "],";
        }
        s << "])";
        return s.str();
    }

    NTL::ZZ get_denom() const {
        return denom;
    }

    NTL::mat_ZZ get_basis() const {
        return basis;
    }

    NTL::mat_ZZ HNF_basis() const 
    {
        NTL::mat_ZZ M = basis;
        NTL::mat_ZZ H;
        NTL::HNF(H, NTL::transpose(M), NTL::determinant(M));
        return H;
    }

    friend class quatlatenum;
};

