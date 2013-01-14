import numpy

def wt2_eisenstein_coeffs(prec = 100):
    #coeff_list = [ -24 ] * prec
    #coeff_list = -24 * numpy.ones(prec, dtype = object)
    coeff_list = -24 * numpy.ones(prec, dtype = numpy.int64)

    for p in primes(1, prec + 1):
        ppow = p
        term = p**2
        last = p

        while ppow < prec:
            ind = ppow
            term_m1 = term - 1
            last_m1 = last - 1
            coeff_list[ind::ppow] *= term_m1
            coeff_list[ind::ppow] /= last_m1

            ppow *= p
            last = term
            term *= p

    coeff_list[0] = 1

    return coeff_list

################################################################################

@cached_function
def E2_newform(r, s, prec):
    """
    E2_newform(r, s, prec):
    Given r, s, this computes the newform associated to f_r_s to precision prec
    by writing it in terms of E2 forms.

    Examples:

    sage: f = E2_newform(7, 6, 100)
    sage: f
    """
    # this is larger than the sturm bound in the weight 50 case, so identifying
    # all forms to this precision identifies them uniquely
    sturm_bd = 52
    wt = r // 2 + s

    file_name = "f_{0}_{1}.sobj".format(r, s)
    nf = load(nf_dir + file_name)
    nf_trunc = nf.q_expansion(sturm_bd)

    gens = wtk_level6_generators(2*wt, sturm_bd)
    coeffs = in_span(nf_trunc, gens)

    gen_set = wtk_level6_generators(2*wt, prec)

    linear_comb = 0
    for (i, gen) in enumerate(gen_set):
        linear_comb += coeffs[i] * gen

    return linear_comb

################################################################################

def E2_newform_modp_twist_lowmem(r, s, p, prec):
    """
    E2_newform(r, s, prec):
    Given r, s, this computes the newform associated to f_r_s to precision prec
    by writing it in terms of E2 forms.

    Examples:

    sage: f = E2_newform(7, 6, 100)
    sage: f
    """
    # this is larger than the sturm bound in the weight 50 case, so identifying
    # all forms to this precision identifies them uniquely
    sturm_bd = 52
    wt = r // 2 + s

    file_name = "f_{0}_{1}.sobj".format(r, s)
    nf = load(nf_dir + file_name)
    nf_trunc = nf.q_expansion(sturm_bd)

    gen_set = wtk_level6_generators(2*wt, sturm_bd)
    print "Finding linear relation"
    coeffs = in_span(nf_trunc, gen_set)

    # it may be the case that some of these coefficients have denominators
    # divisible by p. If so, compute the coefficients exactly, then reduce.

    # construct the eisenstein series
    Fp = Integers(p)
    R = PowerSeriesRing(Fp, 'q')
    #E21 = eisenstein_series_qexp_normalized(2, q, prec)
    #E22 = vd_operator_trunc(E21, 2)
    #E23 = vd_operator_trunc(E21, 3)
    #E26 = vd_operator_trunc(E21, 6)

    print "Generating Eisenstein Series"
    E21 = wt2_eisenstein_coeffs(prec)

    print "Constructing shifts"
    E22 = numpy.zeros(prec, dtype = numpy.int64)
    shift = 0 if prec % 2 == 0 else 1
    E22[0::2] = E21[0:prec // 2 + shift]

    E23 = numpy.zeros(prec, dtype = numpy.int64)
    shift = 0 if prec % 3 == 0 else 1
    E23[0::3] = E21[0:prec // 3 + shift]

    E26 = numpy.zeros(prec, dtype = numpy.int64)
    shift = 0 if prec % 6 == 0 else 1
    E26[0::6] = E21[0:prec // 6 + shift]

    # do this so that
    # f1 = -1 - 24*q - 48*q^2 - 48*q^3 - 96*q^4 - 144*q^5 - 72*q^6 + ...
    # f2 = -2 - 24*q - 24*q^2 - 120*q^3 - 24*q^4 - 144*q^5 - 120*q^6 + ...
    # f3 = -24*q - 24*q^2 - 72*q^3 - 24*q^4 - 144*q^5 - 72*q^6 + ...
    #f3, f2, f1 = [ c1*E21 + c2*E22 + c3*E23 + c6*E26 for (c1, c2, c3, c6) in all_coeffs ]
    all_coeffs = [(1, -2, -1, 2), (1, -2, 1, -2), (1, -1, -2, 1)]

    print "Constructing generators"
    f1 = E21 - E22 - 2*E23 + E26
    f2 = E21 - 2*E22 + E23 - 2*E26
    f3 = E21 - 2*E22 - E23 + 2*E26

    ff1 = R(list(f1)).add_bigoh(prec)
    ff2 = R(list(f2)).add_bigoh(prec)
    ff3 = R(list(f3)).add_bigoh(prec)

    exps = lex_exponents(wt)
    modp_coeffs = [ Fp(c) for c in coeffs ]

    print "Constructing linear combination"
    linear_comb = 0

    i = 0
    for (a, b, c) in exps:
        if c <= 1:
            coeff = modp_coeffs[i]
            if coeff != 0:
                print (a, b, c)
                linear_comb += coeff * ff1**a * ff2**b * ff3**c
            i += 1

    print "Twisting"
    return twist_operator(linear_comb)

################################################################################

def E2_newform_modp_twist(r, s, p, prec):
    """
    E2_newform(r, s, prec):
    Given r, s, this computes the newform associated to f_r_s to precision prec
    by writing it in terms of E2 forms.

    Examples:

    sage: f = E2_newform(7, 6, 100)
    sage: f
    """
    # this is larger than the sturm bound in the weight 50 case, so identifying
    # all forms to this precision identifies them uniquely
    sturm_bd = 52
    wt = r // 2 + s

    file_name = "f_{0}_{1}.sobj".format(r, s)
    nf = load(nf_dir + file_name)
    nf_trunc = nf.q_expansion(sturm_bd)

    gen_set = wtk_level6_generators(2*wt, sturm_bd)
    print "Finding linear relation"
    coeffs = in_span(nf_trunc, gen_set)

    # it may be the case that some of these coefficients have denominators
    # divisible by p. If so, compute the coefficients exactly, then reduce.

    try:
        modp_coeffs = [ Integers(p)(c) for c in coeffs ]
        print "Creating generators (mod p)"
        modp_gens = wtk_level6_generators_mod(2*wt, p, prec)

        print "Constructing linear combination"
        linear_comb = 0
        for (i, gen) in enumerate(modp_gens):
            linear_comb += modp_coeffs[i] * gen
    except:
        print "Creating generators"
        gens = wtk_level6_generators(2*wt, prec)

        print "Constructing linear combination"
        linear_comb = 0
        for (i, gen) in enumerate(gens):
            linear_comb += coeffs[i] * gen

        print "Reducing modulo p"
        linear_comb = linear_comb % p

    return twist_operator(linear_comb)

################################################################################

def lex_exponents(n):
    """
    lex_exponents(n):
    This returns a list of the exponents of the homogeneous polynomial in three
    variables of degree n in lexicographic order.

    Examples:
    sage: S = SFAHomogeneous(QQ)
    sage: S([2]).expand(3)
    x0^2 + x0*x1 + x1^2 + x0*x2 + x1*x2 + x2^2
    sage: lex_exponents(2)
    [(2, 0, 0), (1, 1, 0), (0, 2, 0), (1, 0, 1), (0, 1, 1), (0, 0, 2)]
    """
    parts = range(n + 1)
    triples = []

    def backtrack(total, index, used):
        for part in parts:
            # don't add parts greater than what's left to fill
            if part > total:
                break
            else:
                t_part = (part,) + used
                if part == total:
                    # if the part we're adding is all that's left to fill, then
                    # append this triple with the appropriate number of zeros
                    # appended to the beginning
                    triples.append((0,) * (2 - index) + t_part)
                elif index < 2:
                    # this is the case where part < total and the index < 2
                    # recurse!
                    backtrack(total - part, index + 1, t_part)

    # backtrack from the beginning situation where no parts are used
    backtrack(n, 0, tuple())

    return triples

################################################################################

def wtk_level6_generators(k, prec):
    """
    wtk_level6_generators(k):
    This returns a generator over level 6 forms of weight k.

    Examples:

    sage: M = ModularForms(6, 4)
    sage: M
    Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(6) of
    weight 4 over Rational Field
    sage: N = M.newforms()
    sage: nf = N[0].q_expansion(10)
    sage: nf
    q - 2*q^2 - 3*q^3 + 4*q^4 + 6*q^5 + O(q^6)
    sage: L = wtk_level6_generators(4, 6)
    sage: for e in L:
              print e

    1 + 48*q + 672*q^2 + 2400*q^3 + 4800*q^4 + 9504*q^5 + O(q^6)
    2 + 72*q + 696*q^2 + 1944*q^3 + 5400*q^4 + 10224*q^5 + O(q^6)
    4 + 96*q + 672*q^2 + 1632*q^3 + 6432*q^4 + 7488*q^5 + O(q^6)
    24*q + 600*q^2 + 1800*q^3 + 4056*q^4 + 7632*q^5 + O(q^6)
    48*q + 624*q^2 + 1296*q^3 + 5232*q^4 + 6048*q^5 + O(q^6)

    sage: C = in_span(nf, L)
    sage: sum([ C[i] * L[i] for i in range(len(L)) ])
    q - 2*q^2 - 3*q^3 + 4*q^4 + 6*q^5 + O(q^6)
    sage: L = wtk_level6_generators(4, 20)
    sage: lin_comb = sum([ C[i] * L[i] for i in range(len(L)) ])
    sage: lin_comb
    q - 2*q^2 - 3*q^3 + 4*q^4 + 6*q^5 + 6*q^6 - 16*q^7 - 8*q^8 + 9*q^9
    - 12*q^10 + 12*q^11 - 12*q^12 + 38*q^13 + 32*q^14 - 18*q^15 + 16*q^16
    - 126*q^17 - 18*q^18 + 20*q^19 + O(q^20)
    sage: nf = N[0].q_expansion(10)
    sage: nf == lin_comb
    True
    """
    all_coeffs = [(1, -2, -1, 2), (1, -2, 1, -2), (1, -1, -2, 1)]

    # construct the eisenstein series
    R.<q> = QQ[[]]
    E21 = eisenstein_series_qexp_normalized(2, q, prec)
    E22 = vd_operator_trunc(E21, 2)
    E23 = vd_operator_trunc(E21, 3)
    E26 = vd_operator_trunc(E21, 6)

    # do this so that
    # f1 = -1 - 24*q - 48*q^2 - 48*q^3 - 96*q^4 - 144*q^5 - 72*q^6 + ...
    # f2 = -2 - 24*q - 24*q^2 - 120*q^3 - 24*q^4 - 144*q^5 - 120*q^6 + ...
    # f3 = -24*q - 24*q^2 - 72*q^3 - 24*q^4 - 144*q^5 - 72*q^6 + ...
    f3, f2, f1 = [ c1*E21 + c2*E22 + c3*E23 + c6*E26 for (c1, c2, c3, c6) in all_coeffs ]

    # the weight k forms are the products of f1, f2, f3 of total degree k / 2
    f1_pows = {}
    f1_pows[0] = 1
    f1_pows[1] = f1

    f2_pows = {}
    f2_pows[0] = 1
    f2_pows[1] = f2

    f3_pows = {}
    f3_pows[0] = 1
    f3_pows[1] = f3

    # we need f1^a and f2^b for 0 <= a, b <= k / 2
    # we only need f3^c for 0 <= c <= 1
    for k in range(2, k + 1):
        f1_pows[k] = f1 * f1_pows[k - 1]
        f2_pows[k] = f2 * f2_pows[k - 1]

    exps = lex_exponents(k // 2)

    for (a, b, c) in exps:
        if c <= 1:
            yield f1_pows[a] * f2_pows[b] * f3_pows[c]

################################################################################

def wtk_level6_generators_mod(k, p, prec):
    """
    wtk_level6_generators(k, p, prec):
    This returns a generator over level 6 newforms of weight k.

    Examples:

    sage: M = ModularForms(6, 4)
    sage: M
    Modular Forms space of dimension 5 for Congruence Subgroup Gamma0(6) of
    weight 4 over Rational Field
    sage: N = M.newforms()
    sage: nf = N[0].q_expansion(10)
    sage: nf
    q - 2*q^2 - 3*q^3 + 4*q^4 + 6*q^5 + O(q^6)
    sage: L = wtk_level6_generators(4, 6)
    sage: for e in L:
              print e

    1 + 48*q + 672*q^2 + 2400*q^3 + 4800*q^4 + 9504*q^5 + O(q^6)
    2 + 72*q + 696*q^2 + 1944*q^3 + 5400*q^4 + 10224*q^5 + O(q^6)
    4 + 96*q + 672*q^2 + 1632*q^3 + 6432*q^4 + 7488*q^5 + O(q^6)
    24*q + 600*q^2 + 1800*q^3 + 4056*q^4 + 7632*q^5 + O(q^6)
    48*q + 624*q^2 + 1296*q^3 + 5232*q^4 + 6048*q^5 + O(q^6)

    sage: C = in_span(nf, L)
    sage: sum([ C[i] * L[i] for i in range(len(L)) ])
    q - 2*q^2 - 3*q^3 + 4*q^4 + 6*q^5 + O(q^6)
    sage: L = wtk_level6_generators(4, 20)
    sage: lin_comb = sum([ C[i] * L[i] for i in range(len(L)) ])
    sage: lin_comb
    q - 2*q^2 - 3*q^3 + 4*q^4 + 6*q^5 + 6*q^6 - 16*q^7 - 8*q^8 + 9*q^9
    - 12*q^10 + 12*q^11 - 12*q^12 + 38*q^13 + 32*q^14 - 18*q^15 + 16*q^16
    - 126*q^17 - 18*q^18 + 20*q^19 + O(q^20)
    sage: nf = N[0].q_expansion(10)
    sage: nf == lin_comb
    True
    """
    all_coeffs = [(1, -2, -1, 2), (1, -2, 1, -2), (1, -1, -2, 1)]

    # construct the eisenstein series
    R.<q> = Integers(p)[[]]
    E21 = eisenstein_series_qexp_normalized(2, q, prec)
    E22 = vd_operator_trunc(E21, 2)
    E23 = vd_operator_trunc(E21, 3)
    E26 = vd_operator_trunc(E21, 6)

    # do this so that
    # f1 = -1 - 24*q - 48*q^2 - 48*q^3 - 96*q^4 - 144*q^5 - 72*q^6 + ...
    # f2 = -2 - 24*q - 24*q^2 - 120*q^3 - 24*q^4 - 144*q^5 - 120*q^6 + ...
    # f3 = -24*q - 24*q^2 - 72*q^3 - 24*q^4 - 144*q^5 - 72*q^6 + ...
    f3, f2, f1 = [ c1*E21 + c2*E22 + c3*E23 + c6*E26 for (c1, c2, c3, c6) in all_coeffs ]

    # the weight k forms are the products of f1, f2, f3 of total degree k / 2
    f1_pows = {}
    f1_pows[0] = 1
    f1_pows[1] = f1

    f2_pows = {}
    f2_pows[0] = 1
    f2_pows[1] = f2

    f3_pows = {}
    f3_pows[0] = 1
    f3_pows[1] = f3

    # we need f1^a and f2^b for 0 <= a, b <= k / 2
    # we only need f3^c for 0 <= c <= 1
    for l in range(2, k + 1):
        f1_pows[l] = f1 * f1_pows[l - 1]
        f2_pows[l] = f2 * f2_pows[l - 1]

    exps = lex_exponents(k // 2)

    for (a, b, c) in exps:
        if c <= 1:
            yield f1_pows[a] * f2_pows[b] * f3_pows[c]

################################################################################

