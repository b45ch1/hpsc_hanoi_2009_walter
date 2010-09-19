/*        Generated by TAPENADE     (INRIA, Tropics team)
    Tapenade 3.4 (r3375) - 10 Feb 2010 15:08
*/
#include "GlobalDeclarations_b.c"

/*
  Differentiation of qr in reverse (adjoint) mode:
   gradient     of useful results: *r *qt
   with respect to varying inputs: *qt *a
*/
void qr_b(double *a, double *ab, double *qt, double *qtb, double *r, double *
        rb, int na) {
    double at, bt, rt, ct, st, rnk, qtnk;
    double atb, btb, rtb, ctb, stb, rnkb, qtnkb;
    int n, m, k, tmp;
    int res;
    int res0;
    int res1;
    int res2;
    int res3;
    int res4;
    int res5;
    int res6;
    double tmp0;
    int res7;
    int res8;
    double tmp1;
    int res9;
    int res10;
    int res11;
    double tmp2;
    int res12;
    int res13;
    double tmp3;
    int adFrom;
    double tmp0b;
    double tmp3b;
    double tempb;
    double tmp2b;
    double tmp1b;
    /* prepare Q and r */
    for (n = 0; n <= na-1; ++n)
        for (m = 0; m <= na-1; ++m) {
            if (n == m)
                tmp = 1;
            else
                tmp = 0;
            res = myindex(n, m, na);
            qt[res] = tmp;
        }
    for (n = 0; n <= na-1; ++n)
        for (m = 0; m <= na-1; ++m) {
            res0 = myindex(n, m, na);
            res1 = myindex(n, m, na);
            r[res0] = a[res1];
        }
    /* MAIN ALGORITHM */
    for (n = 0; n <= na-1; ++n) {
        adFrom = n + 1;
        for (m = adFrom; m <= na-1; ++m) {
            /* defining coefficients of the Givens rotation */
            res2 = myindex(n, n, na);
            pushreal8_(at);
            at = r[res2];
            res3 = myindex(m, n, na);
            pushreal8_(bt);
            bt = r[res3];
            pushreal8_(rt);
            rt = sqrt(at*at + bt*bt);
            pushreal8_(ct);
            ct = at/rt;
            pushreal8_(st);
            st = bt/rt;
            for (k = 0; k <= na-1; ++k) {
                /* update r */
                res4 = myindex(n, k, na);
                pushreal8_(rnk);
                rnk = r[res4];
                res5 = myindex(n, k, na);
                res6 = myindex(m, k, na);
                tmp0 = ct*rnk + st*r[res6];
                pushreal8_(r[res5]);
                r[res5] = tmp0;
                res7 = myindex(m, k, na);
                res8 = myindex(m, k, na);
                tmp1 = -st*rnk + ct*r[res8];
                pushreal8_(r[res7]);
                r[res7] = tmp1;
                /* update Q */
                res9 = myindex(n, k, na);
                pushreal8_(qtnk);
                qtnk = qt[res9];
                res10 = myindex(n, k, na);
                res11 = myindex(m, k, na);
                tmp2 = ct*qtnk + st*qt[res11];
                pushreal8_(qt[res10]);
                qt[res10] = tmp2;
                res12 = myindex(m, k, na);
                res13 = myindex(m, k, na);
                tmp3 = -st*qtnk + ct*qt[res13];
                pushreal8_(qt[res12]);
                qt[res12] = tmp3;
            }
        }
        pushinteger4_(adFrom);
    }
    for (n = na-1; n >= 0; --n) {
        popinteger4_(&adFrom);
        for (m = na-1; m >= adFrom; --m) {
            stb = 0.0;
            ctb = 0.0;
            for (k = na-1; k >= 0; --k) {
                popreal8_(&qt[res12]);
                tmp3b = qtb[res12];
                qtb[res12] = 0.0;
                qtb[res13] = qtb[res13] + ct*tmp3b;
                tmp2b = qtb[res10];
                ctb = ctb + qtnk*tmp2b + qt[res13]*tmp3b;
                stb = stb - qtnk*tmp3b;
                qtnkb = ct*tmp2b - st*tmp3b;
                popreal8_(&qt[res10]);
                qtb[res10] = 0.0;
                popreal8_(&qtnk);
                popreal8_(&r[res7]);
                tmp1b = rb[res7];
                rb[res7] = 0.0;
                rb[res8] = rb[res8] + ct*tmp1b;
                tmp0b = rb[res5];
                ctb = ctb + rnk*tmp0b + r[res8]*tmp1b;
                popreal8_(&r[res5]);
                stb = stb + r[res6]*tmp0b - rnk*tmp1b + qt[res11]*tmp2b;
                qtb[res11] = qtb[res11] + st*tmp2b;
                qtb[res9] = qtb[res9] + qtnkb;
                rnkb = ct*tmp0b - st*tmp1b;
                rb[res5] = 0.0;
                rb[res6] = rb[res6] + st*tmp0b;
                popreal8_(&rnk);
                rb[res4] = rb[res4] + rnkb;
            }
            popreal8_(&st);
            rtb = -(at*ctb/(rt*rt)) - bt*stb/(rt*rt);
            if (at*at + bt*bt == 0.0)
                tempb = 0.0;
            else
                tempb = rtb/(2.0*sqrt(at*at+bt*bt));
            btb = 2*bt*tempb + stb/rt;
            popreal8_(&ct);
            atb = 2*at*tempb + ctb/rt;
            popreal8_(&rt);
            popreal8_(&bt);
            rb[res3] = rb[res3] + btb;
            popreal8_(&at);
            rb[res2] = rb[res2] + atb;
        }
    }
    *ab = 0.0;
    for (n = na-1; n >= 0; --n)
        for (m = na-1; m >= 0; --m) {
            ab[res1] = ab[res1] + rb[res0];
            rb[res0] = 0.0;
        }
    for (n = na-1; n >= 0; --n)
        for (m = na-1; m >= 0; --m)
            qtb[res] = 0.0;
}
