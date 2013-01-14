# Miscellaneous tools for computations with modular forms.                     

import collections
import glob
import itertools
import subprocess
from functools import wraps

################################################################################
# Relevant directories
################################################################################

base_dir = "/home/noname/grad/modular_forms/code/"
ap_dir = base_dir + "ap_data/"
image_dir = base_dir + "images/"
nf_dir = base_dir + "newforms/"
nf_mod_dir = base_dir + "newforms_modp/"
ext_dir = base_dir + "extended_newforms/"
shimura_binary = base_dir + "shimura_prime_coeffs"

r_list = [1, 5, 7, 11, 13, 17, 19, 23]
s_list = [0, 4, 6, 8, 10, 14]

################################################################################
# function to attach all files
################################################################################

def attach_all():
    to_attach = glob.glob(base_dir + '*.sage')
    for f in to_attach:
        print "attaching", f
        attach(f)

################################################################################
# memoization wrapper
################################################################################

def cached_function(func):
    cache = {}
    @wraps(func)
    def wrap(*args):
        if args not in cache:
            cache[args] = func(*args)
        return cache[args]
    return wrap

################################################################################
# Utilities to generate forms of specific types
################################################################################

@cached_function
def delta_series_qexp(q, prec = 1000):
    """
    delta_series_qexp(q, prec = 1000):
    Returns the q-expansion of the Ramanujan Delta function to precision p.
    """
    return q * eta_series_qexp(q, prec - 1)**24

################################################################################

def find_eta_d_r_level(d, r):
    m = 1
    mul = r
    while mul % 24 != 0:
        mul += r
        m += 1
    return m * d

#@cached_function
def eisenstein_series_qexp_normalized(k, v, prec = 1000):
    """
    eisenstein_series_qexp_normalized(k, v, prec = 1000):
    Returns the q-expansion of the normalized weight-k Eisenstein series
    on SL_2(ZZ) to precision prec in the variable v. (The normalization chosen
    here is the one that forces the constant term to be 1.)

    Examples:

    >>> eisenstein_series_qexp_normalized(6, q, 5)
    1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 + O(q^5)
    """
    R = v.parent()
    cc = -2*k / bernoulli(k)
    ee = k - 1
    coeff_list = [ cc ] * prec

    for p in primes(1, prec + 1):
        ppow = p
        mult = p**ee
        term = mult**2
        last = mult

        while ppow < prec:
            ind = ppow
            term_m1 = term - 1
            last_m1 = last - 1
            while ind < prec:
                coeff_list[ind] *= term_m1
                coeff_list[ind] /= last_m1
                ind += ppow

            ppow *= p
            last = term
            term *= mult

    coeff_list[0] = 1

    return R(coeff_list).add_bigoh(prec)

################################################################################

@cached_function
def eisenstein24_series_qexp(k, v, prec = 1000):
    """
    eisenstein24_series_qexp(k, v, prec = 1000):
    Returns the q-expansion of E_k(24z) to precision p, where E_k is the weight
    k Eisenstein series.

    Examples:

    >>> R.<q> = QQ[[]]
    >>> eisenstein24_series_qexp(6, q, 100)
    1 - 504*q^24 - 16632*q^48 - 122976*q^72 - 532728*q^96 + O(q^100)
    """
    R = v.parent()
    cc = -2*k / bernoulli(k)
    ee = k - 1
    coeff_list = [ cc * int(i % 24 == 0) for i in xrange(prec) ]

    for p in primes(prec // 24 + 1):
        ppow = p
        mult = p**ee
        term = mult**2
        last = mult

        while 24*ppow < prec:
            ind = ppow
            term_m1 = term - 1
            last_m1 = last - 1
            while 24*ind < prec:
                coeff_list[24*ind] *= term_m1
                coeff_list[24*ind] /= last_m1
                ind += ppow

            ppow *= p
            last = term
            term *= mult

    coeff_list[0] = 1

    return R(coeff_list).add_bigoh(prec)

################################################################################

@cached_function
def eta_series_qexp(v, prec = 1000):
    """
    eta_series_qexp(z, prec = 1000):
    Returns the q-expansion for eta(z) to precision p, without the q^(1/24).

    Examples:

    >>> R.<x> = QQ[]
    >>> eta_series_qexp(x, 40)
    1 - x - x^2 + x^5 + x^7 - x^12 - x^15 + x^22 + x^26 + O(x^40)
    """
    def exps():
        n = 1
        while True:
            yield n * (3*n - 1) // 2
            yield n * (3*n + 1) // 2
            n += 1

    R = v.parent()
    coeff_list = [0] * prec
    coeff_list[0] = 1

    for term, sign in itertools.izip(exps(), itertools.cycle((-1, -1, 1, 1))):
        if term >= prec:
            break
        coeff_list[term] = sign

    return R(coeff_list).add_bigoh(prec)

################################################################################

@cached_function
def eta24_series_qexp(v, prec = 1000):
    """
    eta24_series_qexp(z, prec = 1000):
    Returns the q-expansion of eta(24z) to precision prec.

    Examples:

    >>> R.<q> = PowerSeriesRing(QQ)
    >>> eta24_series_qexp(q, 300)
    q - q^25 - q^49 + q^121 + q^169 - q^289 + O(q^300)
    """
    def exps():
        k = 1
        while True:
            yield (1 - 6*k)**2
            yield (6*k + 1)**2
            k += 1

    R = v.parent()
    coeff_list = [0] * prec
    coeff_list[1] = 1

    for term, sign in itertools.izip(exps(), itertools.cycle((-1, -1, 1, 1))):
        if term >= prec:
            break
        coeff_list[term] = sign

    return R(coeff_list).add_bigoh(prec)

################################################################################

@cached_function
def eta8_cubed_series_qexp(v, prec = 1000):
    """
    eta8_cubed_series_qexp(z, prec = 1000):
    Returns the q-expansion of eta(8z)^3 to precision prec.

    Examples:

    >>> R.<q> = PowerSeriesRing(QQ)
    >>> eta24_series_qexp(q, 300)
    q - q^25 - q^49 + q^121 + q^169 - q^289 + O(q^300)
    """
    R = v.parent()
    coeff_list = [0] * prec
    coeff_list[1] = 1
    
    coeff = 3
    coeff2 = 9
    while coeff2 < prec:
        if coeff % 4 == 1:
            coeff_list[coeff2] = coeff
        else:
            coeff_list[coeff2] = -coeff
        coeff += 2
        coeff2 = coeff**2

    return R(coeff_list).add_bigoh(prec)

################################################################################

def eta_d_r_series_qexp(v, d, r, prec=1000):
    """
    Returns the q-expansion of eta(d*z)^r to precision prec.
    """
    if d*r % 24 != 0:
        raise ValueError("d*r must be divisible by 24")

    ee = vd_operator_trunc(eta_series_qexp(v, prec), d)
    pp = v**(d*r // 24) * ee**r
    return pp.O(prec)

################################################################################

@cached_function
def pr_series_qexp(v, r, prec = 1000):
    """
    pr_series_qexp(v, r, prec = 1000):
    Returns the q-expansion of the rth power of the partition function.
    Pr(n) = \sum_{n = 0}^{\infty} p_r(n) q^n
          = \prod_{n = 1}^{\infty} (1 - q^n)^(-1)
          = 1 + rq + ...

    Examples:

    >>> R.<q> = QQ[[]]
    >>> pr_series_qexp(5, q, 5)
    1 + 5*q + 20*q^2 + 65*q^3 + 190*q^4 + O(q^5)
    """
    pp = eta_series_qexp(v, prec)**(-r)
    return pp

################################################################################
# Scripts related to the paper on rth power partitions by Boylan
################################################################################

def partition_number(n, P = {0: 1, 1: 1}):
    """
    partition_number(n):
    Given a non-negative integer n, this returns the number of integer partitions
    of n; i.e. the number of ways to write n as a non-increasing sum of positive
    integers

    Examples:
    
    >>> partition_number(4)
    5
    >>> integer_partitions(5)
    [[1, 1, 1, 1, 1], [2, 1, 1, 1], [2, 2, 1], [3, 1, 1], [3, 2], [4, 1], [5]]
    >>> len(_)
    5
    """
    # if already computed, return the cached value
    if n in P:
        return P[n]
    if n < 0:
        return 0

    # since P(m) = 0 for m < 0, we don't need to check all pentagonal numbers
    max_k = int((1 + math.sqrt(1 + 24*n)) / 6) + 1

    # this is the sign in the alternating sum
    s = (-1)**(max_k + 1)

    # this comes from the pentagonal number theorem
    ss = 0
    for k in xrange(max_k, 0, -1):
        t1 = partition_number(n - k * (3*k - 1) // 2)
        t2 = partition_number(n - k * (3*k + 1) // 2)
        ss += s * (t1 + t2)
        s *= (-1)

    P[n] = ss
    return ss

################################################################################

def rth_power_partition_number(r, n):
    r"""
    rth_power_partition_number(r, n):
    This returns p_r(n); the nth coefficient of the rth power of the generating
    function for p(n).

    Examples:

    >>> e = eta_series_qexp(q, 8)^(-5)
    >>> e
    1 + 5*q + 20*q^2 + 65*q^3 + 190*q^4 + 506*q^5 + 1265*q^6 + 2990*q^7 + O(q^8)
    >>> rth_power_partition_number(5, 3)
    65
    >>> rth_power_partition_number(5, 6)
    1265
    """
    if n == 0:
        return 1
    if r == 1:
        return partition_number(n)

    P = [ partition_number(k) for k in xrange(0, n + 1) ]
    Pr = [0] * (n + 1)
    Pr[0] = 1

    # this bit uses the j.c.p. pure power recurrence to compute the rth
    # power coefficients of the partition generating function
    for k in xrange(1, n + 1):
        ss = 0
        for i in xrange(1, k + 1):
            ss += P[i] * ((r + 1) * i - k) * Pr[k - i]
        Pr[k] = ss // k

    return Pr[n]

################################################################################
# Utilities related to operators on integral-weight modular forms
################################################################################

def drl_operator(F, r, l):
    """
    drl_operator(f, r, l):
    This returns f | D_r(l) := (f * phi_l(z)^r) | U(l)
    where phi_l(z) = eta(l^2 z) / eta(z)
    """
    R = F.parent()
    q = R.gen()

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
    else:
        prec = len(list(F))

    ee = eta_series_qexp(q, prec)
    eel2 = ud_operator(ee, l**2)
    quo = (ee / eel2)**r
    pp = F * quo

    return ud_operator(pp, l)

################################################################################

def hecke_tpk_operator(F, p, k, chi=None):
    """
    hecke_tpk_operator(F, p, k):
    This returns the action of the Hecke operator T_{p, k, chi_12} on the
    weight k form F.
    """
    if not is_prime(p):
        raise ValueError('p must be prime')

    if chi is None:
        chi = lambda x: 1

    R = F.parent()
    coeffs = list(F)
    prec = len(coeffs)

    new_coeffs = []
    n = 0
    while n*p < prec:
        cc = coeffs[n*p]
        if n % p == 0:
            cc += chi(p) * p**(k - 1) * coeffs[n // p]
        new_coeffs.append(cc)
        n += 1

    return R(new_coeffs).add_bigoh(n)

################################################################################

def theta_operator(F, n = 1):
    """
    theta_op(F, n = 1):
    Given a q-expansion F, this applies the theta operator n times.
    """
    R = F.parent()
    coeffs = list(F)
    lc = len(coeffs)

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
    else:
        prec = lc

    for i in xrange(lc):
        coeffs[i] *= i**n

    return R(coeffs).add_bigoh(prec)

################################################################################

def twist_operator(F, chi=None):
    """
    twist_op(F, chi = None):
    Given a q-series F and a character chi, this returns the twist by chi_12
    """
    R = F.parent()
    coeffs = list(F)
    lc = len(coeffs)

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
    else:
        prec = lc

    if chi is None:
        chi = kronecker_character(3)

    for i in xrange(lc):
        coeffs[i] *= chi(i)

    return R(coeffs).add_bigoh(prec)

################################################################################

def ud_operator(F, d):
    """
    ud_operator(F, d, n = 1):
    This applies the hecke U(d) operator to F (n times).
    Given a q-expansion F = \sum c(n) q^n and a positive integer d, the U_d
    operator is defined by F | U_d = \sum c(d n) q^{n}.

    If d | N, then
        U_d : M_k(Gamma0(N), chi) -> M_k(Gamma0(N), chi)

    and U_d maps cusp forms to cusp forms.
    """
    R = F.parent()
    coeffs = list(F)

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
    else:
        prec = len(coeffs)

    new_prec = prec // d
    new_coeffs = coeffs[0::d]

    return R(new_coeffs).add_bigoh(new_prec)

################################################################################

def vd_operator(F, d):
    """
    vd_operator(F, d):
    Given a q-expansion F = \sum c(n) q^n and a positive integer d, the V_d
    operator is defined by F | V_d = \sum c(n) q^{d n}.

    We have
        V_d : M_k(Gamma0(N), chi) -> M_k(Gamma0(N*d), chi)
    
    and V_d maps cusp forms to cusp forms.
    """
    R = F.parent()
    coeffs = list(F)

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
        coeffs.extend([0] * (prec - len(coeffs)))
    else:
        prec = len(coeffs)

    new_prec = prec * d
    new_coeffs = [0] * new_prec
    new_coeffs[0::d] = coeffs

    return R(new_coeffs).add_bigoh(new_prec)


################################################################################

def vd_operator_trunc(F, d):
    """
    vd_operator(F, d):
    Given a q-expansion F = \sum c(n) q^n and a positive integer d, the V_d
    operator is defined by F | V_d = \sum c(n) q^{d n}.

    We have
        V_d : M_k(Gamma0(N), chi) -> M_k(Gamma0(N*d), chi)
    
    and V_d maps cusp forms to cusp forms.

    Examples:

    >>> E2 = weight2_eisenstein_series_qexp(q, 6)
    >>> E2
    1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
    >>> vd_operator_trunc(E2, 2)
    1 - 24*q^2 - 72*q^4 + O(q^6)
    >>> vd_operator_trunc(E2, 3)
    1 - 24*q^3 + O(q^6)
    """
    R = F.parent()
    coeffs = list(F)

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
        coeffs.extend([0] * (prec - len(coeffs)))
    else:
        prec = len(coeffs)

    new_coeffs = [0] * prec
    
    i = 0
    while d*i < prec:
        new_coeffs[d*i] = coeffs[i]
        i += 1

    return R(new_coeffs).add_bigoh(prec)

################################################################################
# Utilities related to half-integral weight modular forms
################################################################################

def hecke_tp2_operator(F, p, l, flag = False):
    """
    hecke_tp2_operator(F, p, l, flag = False):
        - F is the q-expansion of a modular form of weight l + 1/2
        - p is a prime

    This returns F | T(p^2, l, chi)
    If flag is True, this returns the eigenvalue, else returns the form.

    Examples:

    >>> f_17_4 = eta24_series_qexp(q)^17 * eisenstein24_series_qexp(4, q)
    >>> F_17_4 = shimura_map(f_17_4, 17, 12)
    >>> F_17_4
    q + 9019770*q^5 - 515282432*q^7 - 855114401460*q^11 - 8296664277034*q^13 + ...
    >>> hecke_tp2(f_17_4, 12, 5, True)
    9019770
    >>> hecke_tp2(F_17_4, 12, 7, True)
    -515282432
    >>> hecke_tp2(F_17_4, 12, 7)
    -515282432*q^17 - 114907982336*q^41 + 928023660032*q^65 + 952241934336*q^89  + ...
    """
    R = F.parent()
    coeffs = list(F)
    prec = len(coeffs)
    new_coeffs = []
    n = 0
    while n * p**2 < prec:
        cc = coeffs[p**2 * n]
        cc += kronecker_character(n * (-1)**l * 12)(p) * p**(l - 1) * coeffs[n]
        if n % p**2 == 0:
            cc += kronecker_character((-1)**l * 12)(p) * p**(2*l - 1) * coeffs[n // p**2]
        if flag and n == l:
            return cc
        new_coeffs.append(cc)
        n += 1

    return R(new_coeffs).add_bigoh(n)

################################################################################

def shimura_map(F, t, l, chi_list=None):
    """
    shimura_coefficient(F, t, l, chi_list)
        - F is the q-expansion of a modular form of weight l + 1/2
        - t, l = t, lambda
        - chi_list is list of characters attached to F

    Here, F is a half-integral weight cusp form of weight l + 1/2.

    Examples:

    >>> f_7 = eta24_series_qexp(q, 10000)^7
    >>> shimura_map(f_7, 7, 3)
    q + 66*q^5 - 176*q^7 - 60*q^11 - 658*q^13 + 414*q^17 - 956*q^19 \
            + 600*q^23 + 1231*q^25 - 5574*q^29 + 3592*q^31 - 11616*q^35
    >>> shimura_map(f_7, 55, 3) / 14
    q + 66*q^5 - 176*q^7 - 60*q^11
    >>> f_7_u5 = ud_operator(f_7, 5)
    # Recall that the U5 operator changes the character
    >>> shimura_map(f_7_u5, 11, 3, [X5, X12]) / 14
    q + 66*q^5 - 176*q^7 - 60*q^11 - 658*q^13 + 414*q^17 - 956*q^19 \
            + 600*q^23 + 1231*q^25 - 5574*q^29 + 3592*q^31 - 11616*q^35 \
            - 8458*q^37 - 19194*q^41
    """
    R = F.parent()
    bd = integer_floor(sqrt(F.degree() / t))

    if chi_list is None:
        chi_list = [ kronecker_character(12),
                     kronecker_character((-1)**l),
                     kronecker_character(t)
                   ]

    if not isinstance(chi_list, list):
        chi_list = [ chi_list ]

    def psi_t(n):
        return prod([ X(n) for X in chi_list ])

    # we compute the coefficients of the shimura image using a dirichlet
    # convolution.
    L = [0]
    for k in xrange(1, bd):
        ss = sum([ psi_t(d) * d**(l - 1) * F[t*(k // d)**2] for d in divisors(k) ])
        L.append(ss)

    return R(L)

################################################################################
# Utilities related to mod p modular forms and level-lowering
################################################################################

def in_span(F, B):
    """
    in_span(F, B, coeffs = True):
    This tries to find a linear combination of elements of B that give F.
    These need to be power series, all with the same error term.
    """
    coeffs = list(F)
    prec = len(coeffs)
    L = []
    R = F.base_ring()

    # construct L to be a list of lists of the coefficients of the elements of B
    for e in B:
        L.append([ e[i] for i in xrange(prec) ])

    M = Matrix(R, L)
    X = Matrix(R, coeffs)

    return M.solve_left(X).list()

################################################################################

def cuspform_congruent_to_lower_level_modp(F, level, weight, p):
    """
    congruent_to_lower_level(f, m):
    Swinnerton-Dyer proved a result which says if p >= 5 is prime and h is a modular
    form of weight k, then theta(h) is congruent (mod p) to a modular form of weight
    k + p + 1.
    Given a q-expansion f, a space M, and a prime p, this attempts to find a
    lower-level space in which f is congruent to a modular form in this space.

    Examples:

    >>> E
    Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 + 37*x + 281 over Rational
    Field
    >>> H = E.modular_form()
    >>> H.parent()
    Modular Forms space of dimension 42 for Congruence Subgroup Gamma0(150) of
    weight 2 over Rational Field

    Is H congruent to a form of lower level mod 5?

    >>> Hq5 = H.q_expansion(100) % 5
    >>> B
    """
    R = F.parent()
    coeffs = list(F)

    if hasattr(F, 'prec') and F.prec() < Infinity:
        prec = F.prec()
    else:
        prec = len(coeffs)
    prec = f.prec()

    divs = [ d for d in divisors(level) if d > 1 ]

    for t_level in divs:
        print t_level
        space = CuspForms(t_level, weight)
        qexp_basis = space.q_expansion_basis(f_prec)
        basis_mod_p = [ b % p for b in qexp_basis ]
        try:
            lin_comb = in_span(f % p, new_B)
            print "Level:", t_level
            return lin_comb
        except:
            pass

    return False

################################################################################
