
# Online Python - IDE, Editor, Compiler, Interpreter
import math
from decimal import Decimal

gamma = 0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495
A = 1.28242712910062263687534256886979172776768892732500119206374002174040630885882646112973649195820237439420646120
# z2 = 0.26106396943075527347386628720298336605780785784492507780922552042229678210188446789076863129427510035210829679669627
# z3 = 0.19143408063318455797017055106920434168973878478494635223844132270658416908847711002251672208423541888096
# z4 = 0.15081566749163775274492006017567606638806650405324149790965818916962925627277165142927830810539831713073
# z5 = 0.124416230941457303658905426513360274809113396225794667397111652701079007147537732979071749824727453921012

# Bernoulli numbers until 29
Bernoulli = [1 ,-1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66, 0, -691/2730, 0, 7/6, 0, -3617/510, 0, 43867/798, 0, -174611/330, 0, 854513/138, 0, -236364091/2730, 0, 8553103/6, 0, -23749461029/870, 0]

# ascending power
def asc_pow(s, j):
    result = 1
    for i in range(0, j):
        result *= s+i
    return result

# 14th 100 degree accurancy zeta
def zeta(s):
    m = 14
    n = 1000
    Sn = 0
    for k in range(1, n+1):
        Sn += 1 / k**s
    B_sum = 0
    for j in range(1, m+1):
        B_sum += Bernoulli[2*j] * asc_pow(s, 2*j-1) / math.factorial(2*j) * n**(-s-2*j+1)
    E = -n**(1-s)/(1-s) - n**(-s)/2 + B_sum

    return Sn + E

# Multiple Aitken Algorithm
def MMA(list):
    l = len(list)
    if l == 1:
        return list[0]
    elif l == 2:
        return (list[0] + list[1]) / 2

    result = []
    for n in range(0, l-2):
        if list[n] == list[n+1]:
            return list[n]
        elif list[n+1] == list[n+2]:
            return list[n+1]
        r = (list[n+1] - list[n+2]) / (list[n] - list[n+1])
        if r == 1:
            return list[n]
        result.append((list[n+1] - r*list[n]) / (1-r))
    return MMA(result)

# last of list of float
def last(list):
    if len(list) == 0:
        return 0
    else:
        return list[-1]

# 200-ple accurancy z_r
def z(r):
    list = []
    for m in range(2, 201):
        pm = 1 if m%2 == 0 else -1
        list.append(last(list) + pm * (1 - 2**(-m)) * zeta(m) / (m + r - 1))
    return MMA(list)

# for r in range(2, 6):
#     print(z(r))

# A3 = 2**(-1/42) * (math.pi)**(-1/14) * A**(4/7) * math.exp(gamma/42 - z(3)/7)
# A4 = 2**(29/1800) * (math.pi)**(1/30) * A**(-2/5) * A3**(4/5) * math.exp(-gamma/120 + z(4)/15)
# A5 = 2**(-3/310) * (math.pi)**(-1/62) * A**(8/31) * A3**(-24/31) * A4**(32/31) * math.exp(gamma/310 - z(5)/31)

# print(Decimal(A3))
# print(Decimal(A4))
# print(Decimal(A5))

A3 = 1.0309167521973920944589053760864771902561187744140625
A4 = 0.9795555269428444233881236868910491466522216796875
A5 = 0.99204797452504001054052196195698343217372894287109375


def hyper_factorial(r, n):
    result = 0
    for k in range(1, n):
        result += (k ** (r-1)) * math.log10(k)
    return result

# r=4
def relative_error(z):
    log10_A = math.log10( A5 )
    log10_U = (-13*z/360 + z**3/12 - z**5/25) * math.log10(math.e)
    log10_B = (-z/30 + z**3/3 - z**4/2 + z**5/5)*math.log10(z)

    log10_measured_value = log10_A + log10_U + log10_B
    log10_true_value = hyper_factorial(5, z)
    return 10**(log10_measured_value - log10_true_value) - 1

for i in range(2, 10):
    print(relative_error(i))
