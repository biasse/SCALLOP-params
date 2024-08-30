#include "quaternions.hpp"

quatlat quatalg::maximal_order() const 
{
    NTL::mat_ZZ mat;
    NTL::ZZ denom;
    mat.SetDims(4, 4);
    NTL::ZZ q(this->q);
    NTL::ZZ two_q = q*2;

    if (q == 1) {
        denom = 2;
        mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 1; mat[0][3] = 0;
        mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = 1;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 2; mat[2][3] = 0;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 2;
    }

    else if (q == 2) {
        denom = 4;
        mat[0][0] = 2; mat[0][1] = 0; mat[0][2] = 2; mat[0][3] = 2;
        mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 2; mat[1][3] = 1;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 4; mat[2][3] = 0;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 4;
    }

    else {
        NTL::ZZ p(this->p);
        assert (q % 4 == 3);
        NTL::ZZ a(0);
        while (!((a*a*p + 1) % q == 0)) {
            a += 1;
        }
        denom = q*2;
        mat[0][0] = q; mat[0][1] = q; mat[0][2] = 0; mat[0][3] = 0;
        mat[1][0] = 0; mat[1][1] = 2; mat[1][2] = 0; mat[1][3] = a*2;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = q; mat[2][3] = q;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = denom;
    }

    quatlat O(mat, denom, *this); 
    assert (O.is_order());
    return O;
};