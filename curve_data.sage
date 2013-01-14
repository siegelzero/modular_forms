def find_curve(M, p):
    """
    Given a modular form M and a prime p, this attempts to find an elliptic curve
    with conductor divisible by p with ap values congruent to the coefficients
    of the power series M, after twisting and taking the appropriate number of
    derivatives.

    Examples:

    sage:
    """
    C = CremonaDatabase()
    my_primes = [ pp for pp in M.exponents() if is_prime(pp) ]
    agree = len(my_primes)

    conductor_bound = C.largest_conductor()
    good_curves = []

    # the conjecture is that we only need to check curves of conductor
    # 6*p^l, where p is the prime we're interested in
    # even more, we only look while l <= 2
    cond = 6*p
    k = 1
    while k <= 2 and cond <= conductor_bound:
        # we look at each optimal curve of the given conductor
        for curve in C.iter_optimal([cond]):
            # for each curve, check twists and thetas
            for a in range(p - 1):
                count = 0
                for pp in my_primes:
                    #t = compute_ap(curve, pp) * pp**a * kronecker_symbol(12, pp) % p
                    t = curve.an(pp) * pp**a * kronecker_symbol(12, pp) % p

                    if t == M[pp] % p:
                        count += 1
                    else:
                        break

                if count >= agree:
                    good_curves.append((a, curve))
        cond *= p
        k += 1

    return good_curves

################################################################################

def compute_all_data(p):
    """
    compute_all_data(p):
    Corollary 3 of the Yang paper says that eta(24z)^r * f(24z) is a Hecke
    eigenform whenever 0 < r < 24 is an odd integer, and f(24z) is one of the
    functions 1, E_4, E_6, E_4^2 = E_8, E_4 * E_6 = E_10, and E_4^2 * E_6 = E_14.

    We will generate these forms, and compute their integer-weight images under
    the Shimura map S_r. Then, we will look for elliptic curves with associated
    weight 2 cusp forms congruent to these images.

    Given a prime p, we will compute these products, their images, and look for
    associated curves.
    """
    file_name = 'prime{0}.txt'.format(p)
    data = open(file_name, 'w')

    # this is just a separator for data blocks
    lb = '#' * 80

    data.write(lb)
    data.write('\n')
    data.write('# p = {0}\n'.format(p))
    data.write(lb)
    data.write('\n\n')

    # curve_data is a dict holding all curves for each form
    # curve_data[r, s] is a list of all curves congruent to the form f_r,s
    curve_data = {}
    all_forms = glob.glob(ext_dir + 'f_*.sobj')

    for r in r_list:
        for s in s_list:
            print (r, s)
            l = r // 2 + s
            t = r
            file_name = ext_dir + "f_{0}_{1}.sobj".format(r, s)

            if file_name not in all_forms:
                print file_name, "not in directory"
                continue

            f_r_s = load(file_name)
            F_r_s = twist_operator(f_r_s)

            # generate a list of all curves with forms congruent mod p
            curve_list = find_curve(F_r_s, p)

            # for each curve, add this half integer weight form to the list of forms
            # we will print out all forms associated with each curve to look for
            # shimura pairs
            for (d, curve) in curve_list:
                if curve not in curve_data:
                    curve_data[curve] = [(r, s, d)]
                else:
                    curve_data[curve] += [(r, s, d)]

    for curve in curve_data:
        # if a curve is a twist of a minimal curve, make sure to print this
        twist = curve.minimal_quadratic_twist()
        if twist[0] in curve_data and twist[1] != 1:
            #print curve
            #print "Twist of", twist[0], "by", twist[1]

            line = "{0}\n".format(curve)
            data.write(line)

            line = "Twist of {0} by {1}\n".format(twist[0], twist[1])
            data.write(line)

        elif twist[0] in curve_data and twist[1] == 1:
            #print curve
            line = "{0}\n".format(curve)
            data.write(line)

        # this is the useful data we want (for now)
        # it's easy to add more info here if we need
        #print
        #print "Conductor:", curve.conductor()
        #print "Rank:", curve.rank()
        #print "Torsion:", curve.torsion_order()
        #print "Complex Multiplication:", curve.has_cm()
        #print
        line = "\nConductor: {0}\nRank: {1}\nTorsion: {2}\nComplex Multiplication: {3}\n\n"
        data.write(line.format(curve.conductor(), curve.rank(), curve.torsion_order(), curve.has_cm()))

        # this is sort of ugly.
        # create a dict of forms, where r_dict[r] is a list of all eisenstein factors
        # that appear with this eta^r
        r_dict = {}
        for (r, s, d) in curve_data[curve]:
            if r not in r_dict:
                r_dict[r] = [ (s, d) ]
            else:
                r_dict[r] += [ (s, d) ]

        # now, iterate through a sorted list of the eta powers
        for r in sorted(r_dict.keys()):
            f_list = r_dict[r]
            # print f_r_s
            if len(f_list) > 0:
                #print "eta^{0}*E_k".format(r),
                #print "for k in", [ e[0] for e in f_list ]
                line = "f_{0}_s for s in {1}".format(r, [ e[0] for e in f_list ])
                data.write(line)
            data.write("\n")
                
        #print
        #print lb
        #print
        data.write("\n{0}\n\n".format(lb))

    data.close()

################################################################################

def compute_all_data_with_thetas(p):
    """
    compute_all_data(p):
    Corollary 3 of the Yang paper says that eta(24z)^r * f(24z) is a Hecke
    eigenform whenever 0 < r < 24 is an odd integer, and f(24z) is one of the
    functions 1, E_4, E_6, E_4^2 = E_8, E_4 * E_6 = E_10, and E_4^2 * E_6 = E_14.

    We will generate these forms, and compute their integer-weight images under
    the Shimura map S_r. Then, we will look for elliptic curves with associated
    weight 2 cusp forms congruent to these images.

    Given a prime p, we will compute these products, their images, and look for
    associated curves.
    """
    file_name = 'prime{0}.txt'.format(p)
    data = open(file_name, 'w')

    # this is just a separator for data blocks
    lb = '#' * 80

    data.write(lb)
    data.write('\n')
    data.write('# p = {0}\n'.format(p))
    data.write(lb)
    data.write('\n\n')

    # curve_data is a dict holding all curves for each form
    # curve_data[r, s] is a list of all curves congruent to the form f_r,s
    curve_data = {}
    all_forms = glob.glob(image_dir + 'f_*.sobj')

    for r in r_list:
        for s in s_list:
            print (r, s)
            l = r // 2 + s
            t = r
            file_name = image_dir + "f_{0}_{1}.sobj".format(r, s)

            if file_name not in all_forms:
                print file_name, "not in directory"
                continue

            F_r_s = load(file_name)

            # generate a list of all curves with forms congruent mod p
            curve_list = find_curve(F_r_s, p)

            # for each curve, add this half integer weight form to the list of forms
            # we will print out all forms associated with each curve to look for
            # shimura pairs
            for (d, curve) in curve_list:
                if curve not in curve_data:
                    curve_data[curve] = [(r, s, d)]
                else:
                    curve_data[curve] += [(r, s, d)]

    for curve in curve_data:
        # if a curve is a twist of a minimal curve, make sure to print this
        twist = curve.minimal_quadratic_twist()
        if twist[0] in curve_data and twist[1] != 1:
            #print curve
            #print "Twist of", twist[0], "by", twist[1]

            line = "{0}\n".format(curve)
            data.write(line)

            line = "Twist of {0} by {1}\n".format(twist[0], twist[1])
            data.write(line)

        elif twist[0] in curve_data and twist[1] == 1:
            #print curve
            line = "{0}\n".format(curve)
            data.write(line)

        # this is the useful data we want (for now)
        # it's easy to add more info here if we need
        #print
        #print "Conductor:", curve.conductor()
        #print "Rank:", curve.rank()
        #print "Torsion:", curve.torsion_order()
        #print "Complex Multiplication:", curve.has_cm()
        #print
        line = "\nConductor: {0}\nRank: {1}\nTorsion: {2}\nComplex Multiplication: {3}\n"
        data.write(line.format(curve.conductor(), curve.rank(), curve.torsion_order(), curve.has_cm()))

        # this is sort of ugly.
        # create a dict of forms, where r_dict[r] is a list of all eisenstein factors
        # that appear with this eta^r
        r_dict = {}
        for (r, s, d) in curve_data[curve]:
            if r not in r_dict:
                r_dict[r] = [ (s, d) ]
            else:
                r_dict[r] += [ (s, d) ]

        # now, iterate through a sorted list of the eta powers
        for r in sorted(r_dict.keys()):
            f_list = r_dict[r]
            # print f_r_s
            if len(f_list) > 0:
                #print "eta^{0}*E_k".format(r),
                #print "for k in", [ e[0] for e in f_list ]
                line = "\nf_{0}_s for s in {1}\n".format(r, [ e[0] for e in f_list ])
                data.write(line)
                for e in f_list:
                    line = "\ts = {0}: theta^{1}\n".format(e[0], e[1])
                    data.write(line)
                
        #print
        #print lb
        #print
        data.write("\n{0}\n\n".format(lb))

    data.close()

################################################################################

def prove_congruence(r, s, p):
    """
    This proves the congruence for f_r_s
    """
    file_name = image_dir + "f_{0}_{1}.sobj".format(r, s)

    if file_name not in glob.glob(image_dir + '*.sobj'):
        return "Nothing to prove"

    f_r_s = load(file_name)
    curves = find_curve(f_r_s, p)
    
    Fp = Integers(p)

    print
    print "#"*80
    print

    for (a, curve) in curves:
        print curve
        label = curve.label()
        print "Label:", label
        print
        cond = curve.conductor()
        old_weight = r // 2 + s
        new_weight = r - 1 + 2*s
        #print "f_{0}_{1} has weight {2}".format(r, s, old_weight)
        print "F_{0}_{1} has weight {2}".format(r, s, new_weight)
        print "M lives in M_2(Gamma0({0}))".format(cond)
        theta_weight = 2 + a*(p + 1) 
        print "Theta^{0}(M) is congruent (mod {1}) to a form of weight {2}".format(a, p, theta_weight)
        new_level = lcm(cond, 12^2)
        print "X_12 * Theta^{0}(M) is congruent (mod {1}) to a form of weight {2} and level {3}".format(a, p, theta_weight, new_level)
        print
        
        if theta_weight > new_weight:
            diff = theta_weight - new_weight
            if diff % (p - 1) != 0:
                print "Error"
            pp = diff // (p - 1)
            print "Trying to prove:"
            print "F_{0}_{1} (E_{2})^{3} = X_12 * Theta^{4}(M) (mod {5})".format(r, s, p - 1, pp, a, p)
            max_weight = theta_weight

        elif theta_weight < new_weight:
            diff = new_weight - theta_weight
            if diff % (p - 1) != 0:
                print "Error"
            pp = diff // (p - 1)
            print "Trying to prove:"
            print "F_{0}_{1} = X_12 * Theta^{4}(M) (E_{2})^{3} (mod {5})".format(r, s, p - 1, pp, a, p)
            max_weight = new_weight
        else:
            # diff == 0
            print "Trying to prove:"
            print "F_{0}_{1} = X_12 * Theta^{2}(M) (mod {3})".format(r, s, a, p)
            max_weight = theta_weight

        bound = sturm_bound(new_level, max_weight)

        print "Level:", new_level
        print "Weight:", max_weight
        print "Sturm bound:", bound
        print

        # list of names of newforms in the newform_modp dir.
        # if the file is already computed, read it in to save repeating
        # the computation
        list1 = glob.glob(nf_mod_dir + "*.sobj")
        
        # the file name is f_r_s_p, where r and s are as usual, and p is the
        # prime f_r_s is reduced modulo
        base_name ="f_{0}_{1}_{2}_{3}.sobj".format(r, s, p, label) 
        file_name1 = nf_mod_dir + base_name

        if file_name1 in list1:
            print "Using computed newform"
            # ensure that all files in modp folder are computed to sufficient
            # precision
            f_r_s_p = load(file_name1)
            save1 = False
        else:
            print "Generating twisted newform"
            # else compute the twisted newform mod p
            #f_r_s_p = E2_newform_modp_twist(r, s, p, bound + 20)
            f_r_s_p = E2_newform_modp_twist_lowmem(r, s, p, bound + 20)
            save1 = True
            #f_r_s_p = twist_operator(E2_newform(r, s, bound + 12)) % p

        print

        lf = list(f_r_s_p)
        len_lf = len(lf)

        list2 = glob.glob(ap_dir + "*.sobj")
        file_name2 = ap_dir + base_name

        if file_name2 in list2:
            print "Using computed a_n values"
            an_s_p = load(file_name2)
            save2 = False
        else:
            print "Generating a_n values"
            #an_s_p = [ kronecker_symbol(12, i) * i**a * curve.an(i) % p for i in xrange(len_lf) ]
            an_s_p = [ Fp(kronecker_symbol(12, i) * i**a * curve.an(i)) for i in xrange(len_lf) ]
            save2 = True

        print

        len_al = len(an_s_p)

        if len_al > bound and lf == an_s_p:
            print "Success."
            print "Verified to precision", len_al

            if save1:
                print "Saving twisted newform"
                save(f_r_s_p, file_name1)

            if save2:
                print "Saving a_n values"
                save(an_s_p, file_name2)
        else:
            raise ValueError("Curve does not match. Increase precision of image.")

        print
        print "#"*80
        print

