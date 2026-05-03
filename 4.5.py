import math

from scipy import stats
def read(ind: int) -> list[float]:
    a = []
    with open("resources/RGR2_A-1_X1-X4.csv", "r") as file:
        for i in file.readlines():
            a.append(float(i.split(';')[ind]))
    return a

def avg(a: list[float]):
    n = len(a)
    sum = 0.0
    for i in range(n):
        sum += a[i]
    return sum / n

def deviation(a: list[float]):
    average = avg(a)
    sum = 0.0
    n = len(a)
    for i in range(n):
        sum += (a[i] - average)**2
    return sum / float(n)

x = read(3)
n = len(x)
x.sort()
l = 1/avg(x)

a = 0.05
k = 5
p_k = 0.2
cum_probs = [i * p_k for i in range(1, k)]

bounds = []
for cp in cum_probs:
    bound = -math.log(1 - cp) / l
    bounds.append(bound)

observed_freq = [0] * k
for value in x:
    if value <= bounds[0]:
        observed_freq[0] += 1
    elif value <= bounds[1]:
        observed_freq[1] += 1
    elif value <= bounds[2]:
        observed_freq[2] += 1
    elif value <= bounds[3]:
        observed_freq[3] += 1
    else:
        observed_freq[4] += 1

expected = n * p_k
chi = sum((obs - expected)**2 / expected for obs in observed_freq)
print(chi)
r = k - 1 - 1
chi_crit = stats.chi2.ppf(1 - a, r)
if chi < chi_crit:
    print("Принимаем H0")
else:
    print("Отвергаем HO")