def CMOV(x, y, b):
    """
    Returns x if b=False; otherwise returns y
    """
    return int(not(bool(b)))*x + int(bool(b))*y

def CTEQ(x, y):
    """
    Returns True if x=y; otherwise returns False
    """
    assert len(x) == len(y)
    r = True
    for xi,yi in zip(x,y):
        r = r & (xi == yi)
    return r

def is_QR(x, p):
    """
    Returns True if x is a quadratic residue; otherwise returns False
    """
    z = x**((p-1)//2)
    if z == (-1%p):
        return False
    else:
        return True

def mult_inv(x, p):
    """
    Returns the inverse of x modulo p. If x is zero, returns 0.
    """
    return x**(p-2)

def absolute(x, p):
    """
    Returns |x|=x if x =< (p-1)/2, ohterwise returns -x modulo p.
    """
    lim = (p-1)//2
    if ZZ(x) > ZZ(lim):
        x = -x%p
    return x

def sq_root(x, p):
    """
    Returns the principal square root defined through fixed formulas.
    """
    if p%4 == 3:
        return x**((p+1)//4)
    elif p%8 == 5:
        z = x**((p+3)//8)
        if z**2 == -x:
            sqrt_of_minusone = sqrt(F(-1))
            z = z*sqrt_of_minusone
        return absolute(z, p)
    else:
        raise("cannot handle this square root")
