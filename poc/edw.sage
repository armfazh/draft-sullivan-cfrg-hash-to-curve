A = F(486662)
B = F(1)
d = F(37095705934669439343138083508754565189542113879843219016388785533085940283555)
a = F(-1)
r = 2^252 + 0x14def9dea2f79cd65812631a5cf5d3ed
sA = sqrt(F(-(A+2)))
ZeroEdw = [F(0),F(1),F(0),F(1)]

assert is_square(a) == true , "a must be square"
assert is_square(d) == false , "d must be non-square"
assert d == F(-(A-2)/(A+2)) , "d wrong"
assert a == F(-1) , "a wrong"


def toMont(x,y):
    return ((1+y)/(1-y), sqrt(-(A+2))*u/x)

def toEdw(P):
    if P == E.point(0):
        return ZeroEdw
    u,v,z = P
    x = sA*u/v
    y = (u-z)/(u+z)
    return [x,y,x*y,1]

def isEdwardsPoint(P):
    x,y,t,z = P
    return (a*x**2+y**2 == z**2+d*t**2) and (x*y==z*t)

def areEqual(P,Q):
    x0,y0,t0,z0 = P
    x1,y1,t1,z1 = Q
    assert isEdwardsPoint(P), "P is not a point in Edw"
    assert isEdwardsPoint(Q), "Q is not a point in Edw"
    return (x0*z1 == x1*z0)  \
       and (y0*z1 == y1*z0)  \
       and (t0*z1 == t1*z0)

def addEdw(P,Q):
    x1,y1,t1,z1 = P
    x2,y2,t2,z2 = Q
    A = x1*x2
    B = y1*y2
    C = d*t1*t2
    D = z1*z2
    E = (x1+y1)*(x2+y2)-A-B
    F = D-C
    G = D+C
    H = B+A
    x3 = E*F
    y3 = G*H
    t3 = E*H
    z3 = F*G
    return x3,y3,t3,z3

def dedaddEdw(P,Q):
    x1,y1,t1,z1 = P
    x2,y2,t2,z2 = Q
    A = (y1-x1)*(y2+x2)
    B = (y1+x1)*(y2-x2)
    C = 2*z1*t2
    D = 2*t1*z2
    E = D+C
    F = B-A
    G = B+A
    H = D-C
    x3 = E*F
    y3 = G*H
    t3 = E*H
    z3 = F*G
    return x3,y3,t3,z3

def smulEdw(k,P):
    K = ZZ(k).digits(2)
    Q = ZeroEdw
    for ki in reversed(K):
        Q = addEdw(Q,Q)
        if ki == 1:
            Q = addEdw(Q,P)
    return Q

def normalEdw(P):
    assert isEdwardsPoint(P)
    x1,y1,t1,z1 = P
    return x1/z1, y1/z1, x1*y1/z1**2, 1

def isZeroEdw(P):
    return areEqual(P,ZeroEdw)

def orderEdw(P):
    for ki in [1,2,4,8,r,2*r,4*r,8*r]:
        if isZeroEdw(smulEdw(ki,P)):
            return ki

def my_abs(z):
    lim = (p-1)//2
    if ZZ(z) > ZZ(lim):
        z = -z
    return z

def my_sqrt(n):
    sqrt_of_minusone = F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
    exp = ZZ( (p+3)//8 )
    b = n**exp
    if b**2 != n:
        # print("case neg\n")
        b = b*sqrt_of_minusone
    # print("root: ",hex(int(b)))
    return my_abs(b)

def ell2Mon(r):
    u = F(2)
    v = -A/(1+u*r**2)
    e = legendre_symbol(v**3+A*v**2+B*v,p)
    # print("case e=%+d\n"%e)
    x = e*v-(1-e)*A/F(2)
    y = -e*my_sqrt(x**3+A*x**2+B*x)
    return E([x,y])

def ell2MonOpt(r):
    u = F(2)
    sqrt_of_minusone = F(0x2b8324804fc1df0b2b4d00993dfbd7a72f431806ad2fe478c4ee1b274a0ea0b0)
    sqrt_of_u = u**ZZ((p+3)//8)
    v = -A/(1+u*r**2)
    l = v**3+A*v**2+B*v
    power = l**ZZ(((p+3)//8))
    if power**4 == l**2 :
        # print("case e=+1\n")
        t = power
        x = v
        y = t
        if t**2 != l:
            # print("case e=+1 neg\n")
            y = y*sqrt_of_minusone
        y = -my_abs(y)
    else:
        # print("case e=-1\n")
        t = power*sqrt_of_u
        l = l*u
        x = v*u*r**2
        y = t*r
        if t**2 != l:
            # print("case e=-1 neg\n")
            y = y*sqrt_of_minusone
        y = my_abs(y)
    return E([x,y])

inputs = [1,2,3,4,5, 7, 13, 1<<7, 1<<8, 1<<64, 1<<64-1, p-1, p+1]
tts = [(alpha, elligator2(alpha), ell2Mon(alpha), ell2MonOpt(alpha)) for alpha in inputs]
for pair in tts:
    assert pair[1] == pair[2] == pair[3], "error in ell2mont t:{0}".format(pair[0])


def ell2Edw(t):
    return normalEdw(toEdw(ell2MonOpt(t)))

def ell2EdwOpt(t):
    U,V,W = ell2MonOpt(t)
    e = sA*U
    f = U+W
    g = V
    h = U-W
    return [e*f,g*h,e*h,f*g]

inputs = [1, 7, 13, 1<<7, 1<<8, 1<<64, 1<<64-1, p-1, p+1]
tts = [(alpha, ell2Edw(alpha), ell2EdwOpt(alpha)) for alpha in inputs]
for pair in tts:
    assert areEqual(pair[1],pair[2]),  "error in ell2edw t:{0}".format(pair[0])

for i in range(2**10):
    t = F.random_element()
    P0 = ell2Edw(t)
    P1 = ell2EdwOpt(t)
    assert areEqual(P0,P1), "error test {0}".format(t)
