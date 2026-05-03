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
new_x = [[cur, 0] for cur in x]
new_y = [[cur, 1] for cur in y]
all = new_x + new_y
all.sort()
sums = [0, 0]
count = {}
for i in range(len(all)):
    if all[i][0] in count:
        count[all[i][0]].append(i + 1)
    else:
        count[all[i][0]] = [i + 1]
for i in range(len(all)):
    sums[all[i][1]] += sum(count[all[i][0]])/len(count[all[i][0]])

u_x = sums[0] - n*(n + 1)/2
u_y = sums[1] - m*(m + 1)/2
u = min(u_x, u_y)

res = stats.mannwhitneyu(x, y)

u_statistic = res.statistic
p_value = res.pvalue

print(f"U статистика: {u}")
print(f"p-value: {p_value}")

alpha = 0.05
if p_value < alpha:
    print("Отвергаем H0: распределения выборок различаются.")
else:
    print("Нет оснований отвергнуть H0: распределения выборок статистически значимо не различаются.")