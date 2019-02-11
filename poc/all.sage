load("primes.sage")

Hex = lambda x: map(hex,map(int,x))

def my_abs(a):
    lim = (p-1)//2
    if ZZ(a) > ZZ(lim):
        a = -a
    return a

def my_sqrt(a):
    sqrt_of_minusone = F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
    exp = ZZ( (p+3)//8 )
    b = a**exp
    if b**2 != a:
        # print("case neg\n")
        b = b*sqrt_of_minusone
    # print("root: ",hex(int(b)))
    return my_abs(b)


def my_elligator2(r,u):
    v = -A/(1+u*r**2)
    e = legendre_symbol(v**3+A*v**2+B*v,p)
    # print("case e=%+d\n"%e)
    x = e*v-(1-e)*A/F(2)
    y = -e*my_sqrt(x**3+A*x**2+B*x)
    return E([x,y])

def ell2curve255(r):
    sqrt_of_minusone = F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
    u = F(2)
    v = -A/(1+u*r**2)
    l = v**3+A*v**2+B*v
    power = l**ZZ(((p+3)//8))
    if power**4 == l**2 :
        # print("case e=+1\n")
        x = v
        y = power
        if power**2 != l:
            # print("case e=+1 neg\n")
            y = y*sqrt_of_minusone
        y = -my_abs(y)
    else:
        # print("case e=-1\n")
        x = v*u*r**2
        l = l*u*r**2
        y = power*r*u**ZZ((p+3)//8)
        if y**2 != l:
            # print("case e=-1 neg\n")
            y = y*sqrt_of_minusone
        y = my_abs(y)
    return E([x,y])


def my_swu(t,u):
    x1 = u
    g1 = x1**3+A*x1+B
    x2 = -B/A*(1+(1/(t**4*g1**2+t**2*g1)))
    g2 = x2**3+A*x2+B
    x3 = t**2*g1*x2
    g3 = x3**3+A*x3+B

    if g1.is_square():
        print("case1\n")
        x = x1
        y = g1**((p+1)//4)
    elif g2.is_square():
        print("case2\n")
        x = x2
        y = g2**((p+1)//4)
    else:
        print("case3\n")
        x = x3
        y = g3**((p+1)//4)
    return E([x,y])

def my_swu_opt(t,u):
    x1 = u
    g1 = x1**3+A*x1+B
    y1 = g1**((p+1)//4)
    if y1**2 == g1:
        print("case1\n")
        x = x1
        y = y1
    else:
        x2 = -B/A*(1+(1/(t**4*g1**2+t**2*g1)))
        g2 = x2**3+A*x2+B
        y2 = g2**((p+1)//4)
        if y2**2 == g2:
            print("case2\n")
            x = x2
            y = y2
        else:
            print("case3\n")
            x = t**2*g1*x2
            y = legendre_symbol(t,p)*(t*y1)**3*y2
    return E([x,y])

def my_sswu(u):
    x2 = -B/A*(1+(1/(u**4-u**2)))
    x3 = -u**2*x2
    g2 = x2**3+A*x2+B
    g3 = x3**3+A*x3+B

    if g2.is_square():
        x = x2
        y = g2**((p+1)//4)
    else:
        x = x3
        y = g3**((p+1)//4)
    return E([x,y])

def my_sswu_opt(u):
    """
    Optimizations
    1) If the character of the input is known
       one expo can be saved.
    2) Squaring the input, thus we know it is
       a square.
    2) The character only determines the sign
       of the y-coord, so we can disregard it.
    """
    t0 = u**2
    t1 = t0-1
    x  = (-B/A)*(t1**2+t0)/(t0*t1)
    g  = (x**2+A)*x+B

    y = g**((p+1)//4)
    if y**2 != g:
        x = -u**2*x
        y = (F(-1)**((p+1)//4))*legendre_symbol(u,p)*u**3*y
    return E([x,y])



def my_icart(u):
    v = (3*A-u**4)/(6*u)
    x = (v**2-B-u**6/27)**((2*p-1)//3)+u**2/3
    y = u*x+v
    return E([x,y])

def my_icart_opt(u, Jacobian=True):
    t = u**2
    v = t*(t**3+18*A*t+108*B)-27*A**2

    z = 6*u
    x = t-v*(4*t*v**2)**((p-2)//3)
    y = -t**2+2*x*t+3*A
    if Jacobian is True:
        x = 12*t*x
        y = 36*t*y
        return E([x/z**2,y/z**3])
    else: # homogeneous == True:
        x = 2*u*x
        return E([x,y,z])

def compare(f0,f1,n):
    l = []
    for i in range(2,n):
        r = F.random_element()
        if f0(r) != f1(r):
            l += [r]
    return l

def doub_mont(P):
    X,Y,Z = P
    # print("X:",hex(int(X)))
    # print("Y:",hex(int(Y)))
    # print("Z:",hex(int(Z)))
    l0 = 3*X**2+2*A*X*Z+Z**2
    l1 = 2*B*Y*Z
    t0 = B*l0**2*Z-(2*X+A*Z)*l1**2
    X2 = t0*l1
    Y2 = -l0*t0+l1**2*(X*l0-Y*l1)
    Z2 = l1**3*Z
    print("X:",hex(int(X2)))
    print("Y:",hex(int(Y2)))
    print("Z:",hex(int(Z2)))
    #return [X2,Y2,Z2]
    return E([X2,Y2,Z2])

def add_mont(P,Q):
    X1,Y1,Z1 = P
    X2,Y2,Z2 = Q
    # print("X:",hex(int(X)))
    # print("Y:",hex(int(Y)))
    # print("Z:",hex(int(Z)))
    l0 = Y2*Z1-Y1*Z2
    l1 = X2*Z1-X1*Z2
    t0 = B*l0**2*Z1*Z2-(X1*Z2+X2*Z1+A*Z1*Z2)*l1**2
    print("l0:",hex(int(l0)))
    print("l1:",hex(int(l1)))
    print("t0:",hex(int(t0)))
    X3 = t0*l1
    Y3 = -l0*t0+l1**2*Z2*(X1*l0-Y1*l1)
    Z3 = l1**3*Z1*Z2
    print("X:",hex(int(X3)))
    print("Y:",hex(int(Y3)))
    print("Z:",hex(int(Z3)))
    #return [X2,Y2,Z2]
    #return E.lift_x(X3/Z3)
    return E([X3,Y3,Z3])



#EOF
