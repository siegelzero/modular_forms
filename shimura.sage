def shimura_image(r, s):
    """
    shimura_image(r, s):
    Given r, s, this computes the shimura image of f_r,s to enough precision to
    identify it in S_{2 lambda}(Gamma0(288)).

    Our hope is to identify it in S_{2 lambda)(Gamma0(6)) twist X_12.
    This returns the desired element of S_{2 lambda)(Gamma0(6)).

    Examples:

    """

    # compute lambda
    l = r // 2 + s

    # compute sturm bd; add two
    bd = 96 * l + 2

    # this is the command to compute the coefficients at primes
    args = shimura_binary + " {0} {1} {2}".format(r, s, bd)

    # this calls the external program to compute the shimura prime coefficients
    print "Computing image"
    pp = subprocess.Popen(args, shell = True, stdout = subprocess.PIPE)
    ss = pp.communicate()[0]
    lines = ss.split("\n")

    # Now, reconstruct the coefficients
    print "Reconstructing coefficients"
    #F = reconstruct_coeffs(output, l)
    D = {}

    # now read off the remaining lines
    for line in lines:
        try:
            a, b = line.split(':')
            p, c = ZZ(a), ZZ(b)
            D[p] = c
        except:
            pass

    # we want to reconstruct as far as possible, so we go out to right before
    # the next prime.
    max_p = max(D.keys())
    coeffs = [0] * bd

    # result wil be a cusp form, so go ahead and make the q term have coeff 1
    coeffs[1] = 1

    # only interested on nonzero coefficients
    p_list = [ p for p in D if D[p] != 0 ]

    # copy over the values at primes
    for p in p_list:
        coeffs[p] = D[p]

    # build up the coefficients using the fact that
    # A(p**r) = A(p) * A(p**(r - 1)) - p**(2*l - 1) * A(p**(r - 2))
    for p in p_list:
        k = 2
        while p**k < bd:
            coeffs[p**k] = coeffs[p] * coeffs[p**(k - 1)] - p**(2*l - 1) * coeffs[p**(k - 2)]
            k += 1

    # now build up the remaining coefficients using multiplicativity
    for n in range(4, bd):
        if coeffs[n] == 0:
            n_fac = factor(n)
            coeffs[n] = prod([ coeffs[p**e] for (p, e) in n_fac ])

    R = PowerSeriesRing(QQ, 'q')

    return R(coeffs).add_bigoh(bd)

################################################################################

def sage_shimura_image(r, s):
    """
    Give r, s, this computes the shimura image of f_r_s
    """
    # compute lambda
    l = r // 2 + s

    # compute sturm bd, and throw in a few extra for good measure
    bd = 96 * l + 12

    print "Computing image"

    R = PowerSeriesRing(QQ, 'q')
    q = R.gen()
    s_prec = bd**2 * r // 24
    eta_pow = eta_series_qexp(q, s_prec + 1)**r
    eis_fac = eisenstein_series_qexp_normalized(s, q, s_prec + 1)

    D = {}
    p = 2
    m = r*p*p
    
    while p < bd:
        ss = 0
        first = r
        while first <= bd and (m - first) % 24 != 0:
            first += 24

        for i in xrange(first, m + 1, 24):
            ss += eta_pow[(i - r) // 24] * eis_fac[(m - i) // 24]

        D[p] = ss + kronecker_symbol(12, p) * kronecker_symbol((-1)**l * r, p) * p**(l - 1)
        p = next_prime(p)
        m = r*p*p

    # we want to reconstruct as far as possible, so we go out to right before
    # the next prime.
    max_p = max(D.keys())
    coeffs = [0] * bd

    # result wil be a cusp form, so go ahead and make the q term have coeff 1
    coeffs[1] = 1

    # only interested on nonzero coefficients
    p_list = [ p for p in D if D[p] != 0 ]

    # copy over the values at primes
    for p in p_list:
        coeffs[p] = D[p]

    # build up the coefficients using the fact that
    # A(p**r) = A(p) * A(p**(r - 1)) - p**(2*l - 1) * A(p**(r - 2))
    for p in p_list:
        k = 2
        while p**k < bd:
            coeffs[p**k] = coeffs[p] * coeffs[p**(k - 1)] - p**(2*l - 1) * coeffs[p**(k - 2)]
            k += 1

    # now build up the remaining coefficients using multiplicativity
    for n in range(4, bd):
        if coeffs[n] == 0:
            n_fac = factor(n)
            coeffs[n] = prod([ coeffs[p**e] for (p, e) in n_fac ])

    R = PowerSeriesRing(QQ, 'q')

    return R(coeffs).add_bigoh(bd)

################################################################################

def reconstruct_coeffs(file_name, l):
    """
    Given hecke tp2 eigenvalues for primes p, this reconstructs all 
    coefficients using multiplicativity.

    Here, ff is a file named f_r_s, where the lines are formatted as 
    p:  coeff

    l is the weight lambda
    """
    ff = open(file_name, 'r')
    D = {}

    # now read off the remaining lines
    for line in ff:
        a, b = line.split(':')
        p, c = ZZ(a), ZZ(b)
        D[p] = c

    # we want to reconstruct as far as possible, so we go out to right before
    # the next prime.
    max_p = max(D.keys())
    bd = next_prime(max_p)
    coeffs = [0] * bd

    # result wil be a cusp form, so go ahead and make the q term have coeff 1
    coeffs[1] = 1

    # only interested on nonzero coefficients
    p_list = [ p for p in D if D[p] != 0 ]

    # copy over the values at primes
    for p in p_list:
        coeffs[p] = D[p]

    # build up the coefficients using the fact that
    # A(p**r) = A(p) * A(p**(r - 1)) - p**(2*l - 1) * A(p**(r - 2))
    for p in p_list:
        k = 2
        while p**k < bd:
            coeffs[p**k] = coeffs[p] * coeffs[p**(k - 1)] - p**(2*l - 1) * coeffs[p**(k - 2)]
            k += 1

    # now build up the remaining coefficients using multiplicativity
    for n in range(4, bd):
        if coeffs[n] == 0:
            n_fac = factor(n)
            coeffs[n] = prod([ coeffs[p**e] for (p, e) in n_fac ])

    R = PowerSeriesRing(QQ, 'q')

    return R(coeffs).add_bigoh(bd)

################################################################################

def generate_images():
    """
    compute_images():
    This computes all shimura images and saves their power series into the
    current directory. The fourier series are computed to enough precision to
    identify them in the space S_{2 lambda}(Gamma_0(288)).
    """
    for r in r_list:
        for s in s_list:
            t = r // 2 + s
            if t > 0:
                print (r, s)
                f = shimura_image(r, s)
                file_name = "f_{0}_{1}.sobj".format(r, s)
                save(f, image_dir + file_name)

################################################################################

def generate_newforms():
    """
    generate_newforms():
    This generates the corresponding newforms and names them according to the
    form they correspond to.
    """
    D = collections.defaultdict(list)

    for r in r_list:
        for s in s_list:
            l = r // 2 + s
            if l > 0:
                D[l] += [ (r, s) ]

    for l in D:
        bd = 96 * l + 2
        M = ModularForms(6, 2*l)
        M_new = M.newforms('a')
        M_new_zz = [ f for f in M_new if f.base_ring() is QQ ]

        for (r, s) in D[l]:
            f = load(image_dir + "f_{0}_{1}.sobj".format(r, s))
            f_trunc = f.add_bigoh(30)

            for nf in M_new_zz:
                nf_trunc = nf.q_expansion(30)
                if twist_operator(nf_trunc) == f_trunc:
                    print (r, s)
                    print nf
                    #nf.q_expansion(bd)
                    print
                    save(nf, nf_dir + "f_{0}_{1}.sobj".format(r, s))

################################################################################

def shimura_families(p):
    """
    shimura_families(p):
    This generates all equivalence classes of shimura images modulo p; i.e.
    this returns a list of classes of forms whose shimura images are
    congruent modulo p
    """
    D = collections.defaultdict(list)

    for r in r_list:
        for s in s_list:
            if r // 2 + s > 0:
                f = load(image_dir + "f_{0}_{1}.sobj".format(r, s))
                f_trunc = f.add_bigoh(50)
                f_mod = f_trunc % p
                D[f_mod].append((r, s))

    # now remove redundancies from the list

    families = D.values()

    return families

    for family in families:
        for (r, s) in reversed(family):
            if (r, s - p + 1) in family:
                family.remove((r, s))

    return families

################################################################################

def family_conj(p):
    """
    family_conj(p):
    Suppose S_r1(f_r1_0) = S_r2(f_r2_0) (mod p).
    Then S_r1(f_r1_k) = S_r2(f_r2_k) (mod p) for k = 4, 6, 8, 10, 14.
    """
    for r1 in r_list[1:]:
        f_r1 = load(image_dir + 'f_{0}_0.sobj'.format(r1))
        f_r1_mod = f_r1.O(100) % p
        for r2 in r_list[1:]:
            # we want r1 > r2
            if r2 <= r1:
                continue

            f_r2 = load(image_dir + 'f_{0}_0.sobj'.format(r2))
            f_r2_mod = f_r2.O(100) % p

            if f_r1_mod == f_r2_mod:
                print "f_{0}_0 = f_{1}_0 (mod {2})".format(r1, r2, p)
                for s in [4, 6, 8, 10, 14]:
                    new_f_r1 = load(image_dir + 'f_{0}_{1}.sobj'.format(r1, s))
                    new_f_r2 = load(image_dir + 'f_{0}_{1}.sobj'.format(r2, s))
                    new_f_r1_mod = new_f_r1.O(100) % p
                    new_f_r2_mod = new_f_r2.O(100) % p

                    if new_f_r1_mod != new_f_r2_mod:
                        print "Failure at f_{0}_{1}, f_{2}_{1}".format(r1, s, r2)
                    else:
                        print "Success at f_{0}_{1}, f_{2}_{1}".format(r1, s, r2)

                print

################################################################################
