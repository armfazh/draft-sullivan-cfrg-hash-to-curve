A = F(486662)
B = F(1)
d = F(37095705934669439343138083508754565189542113879843219016388785533085940283555)
a = F(-1)
r = 2^252 + 0x14def9dea2f79cd65812631a5cf5d3ed

assert is_square(a) == true , "a must be square"
assert is_square(d) == false , "d must be non-square"
assert d == F(-(A-2)/(A+2)) , "d wrong"
assert a == F(-1) , "a wrong"


def toMont(x,y):
    return ((1+y)/(1-y), sqrt(-(A+2))*u/x)

def toEdw(P):
    u,v,z = P
    x = sqrt(F(-(A+2)))*u/v
    y = (u-1)/(u+1)
    return [x,y,x*y,1]

def isEdwardsPoint(P):
    x,y,t,z = P
    return (a*x**2+y**2 == z**2+d*t**2) and (x*y==z*t)

def areEqual(P,Q):
    x0,y0,t0,z0 = P
    x1,y1,t1,z1 = Q
    return x0*y0*t1*z1 == x1*y1*t0*z0

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
    Q = [0,1,0,1]
    for ki in reversed(K):
        Q = addEdw(Q,Q)
        if ki ==1:
            Q = addEdw(Q,P)
    return Q

def normalEdw(P):
    assert isEdwardsPoint(P)
    x1,y1,t1,z1 = P
    return x1/z1, y1/z1, x1*y1/z1**2, 1
    
    
    
    
    

