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

a = 0.05
mu_0 = 62.55
x = read(2)
n = 108
s = deviation(x)
mu = avg(x)
t = (mu - mu_0)/((s/n)**0.5)
t_crit = stats.t.ppf(1 - a/2, n-1)
print(f"t = {t}")
print(f"s = {s}")
print(f"avg(x) = {mu}")
print(f"t_crit = {t_crit}")
if( abs(t) > t_crit):
    print("Отвергаем H0")
else:
    print("Принимаем H0")