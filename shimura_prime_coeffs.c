#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

void eis_series(fmpz_poly_t *eis, ulong k, ulong prec);
void eta_series(fmpz_poly_t *eta, ulong prec);

int main(int argc, char *argv[])
{
    ulong r, s, p, prec, s_prec;
    ulong i, m, lambda, char1, char2;
    fmpz_poly_t eis, eta, eta_r;
    fmpz_t p1, c1, c2, ss1, ss, pp, ppl, pplc;

    // input r, s, and the precision desired
    if (argc == 4) {
        r = atol(argv[1]);
        s = atol(argv[2]);
        prec = atol(argv[3]);
        prec = n_nextprime(prec, 1) - 1;
    }

    if (argc != 4) {
        printf("usage: shimura_coeffs r s n\n");
        return EXIT_FAILURE;
    }

    // Compute the weight of the form f_r_s
    lambda = r / 2 + s;

    //printf("%d\n", lambda);

    // This is the necessary precision of the input
    s_prec = (prec + r) * prec * r / 24;
    //printf("Necessary precision: %d terms\n", s_prec);

    // First compute eta
    //printf("Computing eta\n");
    fmpz_poly_init(eta);
    eta_series(&eta, s_prec);

    // Next compute eta^r
    //printf("Computing eta^%d\n", r);
    fmpz_poly_init(eta_r);
    fmpz_poly_pow_trunc(eta_r, eta, r, s_prec);
    fmpz_poly_clear(eta);

    // Next compute E_s
    //printf("Computing E_%d\n", s);
    
    fmpz_poly_init(eis);
    if(s != 0) {
        eis_series(&eis, s, s_prec);
    }

    //printf("Computing coefficients\n");

    // Instead of computing the product, we will compute each coefficient as
    // needed. This is potentially slower, but saves memory.
    //
    fmpz_init(p1);
    fmpz_init(c1);
    fmpz_init(c2);
    fmpz_init(ss);
    fmpz_init(ss1);

    // compute alpha(p^2*r)
    p = 2;
    m = r*p*p;
    while(p <= prec) {
        fmpz_set_ui(ss1, 0);

        if(s != 0) {
            for(i = r; i <= m; i += 24) {
                if((m - i) % 24 == 0) {
                    fmpz_poly_get_coeff_fmpz(c1, eta_r, (i - r) / 24);
                    fmpz_poly_get_coeff_fmpz(c2, eis, (m - i) / 24);
                    fmpz_addmul(ss1, c1, c2);
                }
            }
        }
        else {
            if((m - r) % 24 == 0)
                fmpz_poly_get_coeff_fmpz(ss1, eta_r, (m - r) / 24);
        }

        /*
        if(p == 199) {
            fmpz_print(ss1);
            printf("\n");
        }
        */

        // compute character X12(p)
        char1 = n_jacobi(12, p);
        
        // compute character ((-1)^(lambda) * n | p)
        if(lambda % 2 == 0) {
            char2 = n_jacobi(r, p);
        }
        else {
            char2 = n_jacobi(-r, p);
        }

        // compute p^(lambda - 1)
        fmpz_set_ui(pp, p);
        fmpz_pow_ui(ppl, pp, lambda - 1);
        fmpz_mul_si(pplc, ppl, char1*char2);

        fmpz_add(ss, ss1, pplc);

        printf("%d:\t", p);
        fmpz_print(ss);
        printf("\n");

        p = n_nextprime(p, 1);
        m = r*p*p;
    }

    fmpz_clear(p1);
    fmpz_clear(c1);
    fmpz_clear(c2);
    fmpz_clear(ss);
    fmpz_clear(ss1);

    return EXIT_SUCCESS;
}


void eis_series(fmpz_poly_t *eis, ulong k, ulong prec)
{
    ulong i, p, ee, p_pow, ind;
    fmpz_t last, last_m1, mult, fp, t_coeff, term, term_m1;

    long bern_fac[15] = {0, 0, 0, 0, 240, 0, -504, 0, 480, 0, -264, 0, 0, 0, -24};

    // Now, we compute the eisenstein series
    // First, compute primes by sieving.

    // allocate memory for the array.
    /*
    char *primes;
    primes = calloc(prec + 1, sizeof(char));

    for(i = 0; i <= prec; i++) {
        primes[i] = 1;
    }

    primes[0] = 0;
    primes[1] = 0;

    p = 2;
    while(p*p <= prec) {
        j = p*p;
        while(j <= prec) {
            primes[j] = 0;
            j += p;
        }

        p += 1;
        while(primes[p] != 1)
            p += 1;
    }
    */

    // Now, create the eisenstein series.
    // We build up each coefficient using the multiplicative properties of the
    // divisor sum function.

    fmpz_poly_set_coeff_si(*eis, 0, 1);

    for(i = 1; i <= prec; i++)
        fmpz_poly_set_coeff_si(*eis, i, bern_fac[k]);

    fmpz_init(last);
    fmpz_init(last_m1);
    fmpz_init(mult);
    fmpz_init(fp);
    fmpz_init(t_coeff);
    fmpz_init(term);
    fmpz_init(term_m1);

    ee = k - 1;
    p = 2;

    while(p <= prec) {
        p_pow = p;
        fmpz_set_ui(fp, p);
        fmpz_pow_ui(mult, fp, ee);
        fmpz_pow_ui(term, mult, 2);
        fmpz_set(last, mult);

        while(p_pow <= prec) {
            ind = p_pow;
            fmpz_sub_ui(term_m1, term, 1);
            fmpz_sub_ui(last_m1, last, 1);

            while(ind <= prec) {
                fmpz_poly_get_coeff_fmpz(t_coeff, *eis, ind);
                fmpz_mul(t_coeff, t_coeff, term_m1);
                fmpz_divexact(t_coeff, t_coeff, last_m1);
                fmpz_poly_set_coeff_fmpz(*eis, ind, t_coeff);
                ind += p_pow;
            }

            p_pow *= p;
            fmpz_set(last, term);
            fmpz_mul(term, term, mult);
        }
        p = n_nextprime(p, 1);
    }

    fmpz_clear(last);
    fmpz_clear(last_m1);
    fmpz_clear(mult);
    fmpz_clear(fp);
    fmpz_clear(t_coeff);
    fmpz_clear(term);
    fmpz_clear(term_m1);
}


void eta_series(fmpz_poly_t *eta, ulong prec)
{
    ulong i, e;
    int no[4] = {1, -1, 1, -1};

    fmpz_poly_set_coeff_si(*eta, 0, 1);

    i = 1;
    e = i*(3*i + 1) / 2;
    while (e <= prec) {
        fmpz_poly_set_coeff_si(*eta, e, no[i % 4]);
        e += (3*i + 2);
        i += 1;
    }

    i = 1;
    e = i*(3*i - 1) / 2;
    while (e <= prec) {
        fmpz_poly_set_coeff_si(*eta, e, no[i % 4]);
        e += (3*i + 1);
        i += 1;
    }
}

