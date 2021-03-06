from operator import contains, eq, ge, gt, le, lt, ne
try:
    from numpy import isclose
except ImportError:
    def isclose(a, b, rtol=1e-05, atol=1e-08):
        return abs(a - b) <= (atol + rtol * abs(b))


def isin(a, b):
    return contains(b, a)


def endswith(a, b):
    return a.endswith(b)


def startswith(a, b):
    return a.startswith(b)


def ncontains(a, b):
    return not contains(a, b)


def isnotclose(a, b):
    return not isclose(a, b)


def isnotin(a, b):
    return not contains(b, a)


def nendswith(a, b):
    return not a.endswith(b)


def nstartswith(a, b):
    return not a.startswith(b)
