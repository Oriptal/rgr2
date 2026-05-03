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

x = read(0)
y = read(1)
n = 108
m = 108
a = 0.05
k = n + m - 2

s_x = deviation(x)
s_y = deviation(y)
s = (s_x*(n - 1) + s_y*(m - 1))/(k)

mu_x = avg(x)
mu_y = avg(y)

T = (mu_x - mu_y)/(s*(1/n + 1/m))**0.5
p_value = 2 * (1 - stats.t.cdf(abs(T), k))
t_crit = stats.t.ppf(1 - a/2, k)

print(f"Критерий Стьюдента")
print(f"T = {T}")
print(f"p_value = {p_value}")
print(f"t = {t_crit}")
print(f"O = ({-t_crit}, {t_crit})")
print(f"W = (-∞, {-t_crit}) ∪ ({t_crit}, ∞)")
if abs(T) >= t_crit:
    print("T принадлежит W")
    print("Отвергаем H0: мат. ожидания не равны")
else:
    print("T принадлежит O")
    print("Принимаем H0: мат. ожидания равны")