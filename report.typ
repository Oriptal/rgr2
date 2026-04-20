#set text(size: 12pt, lang: "ru")
#set page(numbering: "1")
#import "@preview/lilaq:0.5.0" as lq
#set par(justify: true)

#show: lq.set-grid(
  stroke-sub: .5pt + luma(90),
)

#let data = json("artifacts/results.json")
#let ds = csv("artifacts/descriptive_statistics.csv")
#let pearson = csv("artifacts/x4_pearson_intervals.csv")
#let tests = data.at("tests")
#let t12 = tests.at("Параметрический тест для X1 и X2")
#let x3 = tests.at("Одновыборочный тест для X3")
#let mw = tests.at("Критерий Манна–Уитни для X1 и X2")
#let x4 = tests.at("Критерий согласия Пирсона для X4")

#let fmt(x, digits: 4) = {
  let factor = calc.pow(10, digits)
  if type(x) == float or type(x) == int {
    str(calc.round(x * factor) / factor)
  } else {
    str(x)
  }
}

#align(center)[
  #text(weight: "bold")[Университет ИТМО]
  #linebreak()
  #text(weight: "bold")[Мегафакультет компьютерных технологий и управления]
  #linebreak()
  #text(weight: "bold")[Факультет программной инженерии и компьютерной техники]

  #v(2em)
  #v(10em)

  #text(weight: "bold", size: 15pt)[ОТЧЕТ ПО РГР №2]
  #linebreak()
  #text(weight: "bold", size: 15pt)[Проверка статистических гипотез]

  #v(10em)
]

#align(right + bottom)[
  #table(
    stroke: none,
    columns: 2,
    align: left,
    [Преподаватель:], [*Селина Елена Георгиевна*],
    table.cell(rowspan: 3)[Подготовили:], [*Деревянко Владимир Владимирович*],
    [], [*Лоскутов Прохор Александрович*],
    [], [*Стивен Франклин Фейвор*],
  )
]

#pagebreak()

= Введение и исходные данные

В работе используются данные варианта #data.at("meta").at("variant") из файла `#data.at("meta").at("data_path")`.
После чтения файла получена выборка объёма #data.at("meta").at("sample_size"), уровень значимости принят равным #data.at("meta").at("alpha"), единицы измерения --- #data.at("meta").at("units").
Проверялись гипотезы о равенстве средних для X1 и X2, о соответствии среднего X3 значению 62.55 и о согласии X4 с показательным распределением с параметром lambda = 0.085.

== Описательная статистика

#table(
  columns: 9,
  inset: 4pt,
  stroke: .5pt + black,
  table.header(
    [Показатель], [n], [mean], [std], [min], [q25], [median], [q75], [max],
  ),
  for row in ds.slice(1) {
    [#row.at(0)]
    [#row.at(1)]
    [#row.at(2)]
    [#row.at(3)]
    [#row.at(4)]
    [#row.at(5)]
    [#row.at(6)]
    [#row.at(7)]
    [#row.at(8)]
  },
)

== Графики исходных данных

#figure(
  image("figures/x1_hist.png", width: 80%),
  caption: [Гистограмма выборки X1],
)

#figure(
  image("figures/x2_hist.png", width: 80%),
  caption: [Гистограмма выборки X2],
)

#figure(
  image("figures/x1_x2_comparison_hist.png", width: 90%),
  caption: [Сравнительные гистограммы X1 и X2],
)

#figure(
  image("figures/x1_x2_boxplot.png", width: 70%),
  caption: [Boxplot для X1 и X2],
)

= Проверка гипотезы о равенстве математических ожиданий (X1, X2)

Проверяются гипотезы: #t12.at("null_hypothesis"); #t12.at("alternative_hypothesis").
Используется #t12.at("criterion"). Основание выбора: #t12.at("justification")

Наблюдаемое значение статистики равно #fmt(t12.at("statistic"), digits: 6), число степеней свободы равно #fmt(t12.at("degrees_of_freedom"), digits: 0), p-value = #fmt(t12.at("p_value"), digits: 6).
Контрольный расчёт через SciPy даёт то же значение статистики #fmt(t12.at("statistic_scipy"), digits: 6) и p-value #fmt(t12.at("p_value_scipy"), digits: 6).
Итоговый вывод: #t12.at("conclusion")
Ошибка I рода в данном пункте означает: #t12.at("type_i_error")

= Проверка гипотезы по X3

Проверяются гипотезы: #x3.at("null_hypothesis"); #x3.at("alternative_hypothesis").
Используется #x3.at("criterion"). Основание выбора: #x3.at("justification")

Наблюдаемое значение статистики равно #fmt(x3.at("statistic"), digits: 6), число степеней свободы равно #fmt(x3.at("degrees_of_freedom"), digits: 0), p-value = #fmt(x3.at("p_value"), digits: 6).
Контрольный расчёт через SciPy даёт то же значение статистики #fmt(x3.at("statistic_scipy"), digits: 6) и p-value #fmt(x3.at("p_value_scipy"), digits: 6).
Итоговый вывод: #x3.at("conclusion")
Ошибка I рода в данном пункте означает: #x3.at("type_i_error")

#figure(
  image("figures/x3_hist.png", width: 80%),
  caption: [Гистограмма выборки X3],
)

= Непараметрический критерий для X1 и X2

Проверяются гипотезы: #mw.at("null_hypothesis"); #mw.at("alternative_hypothesis").
Используется #mw.at("criterion"). Основание выбора: #mw.at("justification")

Наблюдаемое значение статистики равно #fmt(mw.at("statistic"), digits: 6), p-value = #fmt(mw.at("p_value"), digits: 6).
Итоговый вывод: #mw.at("conclusion")
Сравнение с параметрическим тестом: #mw.at("details").at("comparison_with_t_test")
Ошибка I рода в данном пункте означает: #mw.at("type_i_error")

= Критерий согласия Пирсона для X4

Проверяются гипотезы: #x4.at("null_hypothesis"); #x4.at("alternative_hypothesis").
Используется #x4.at("criterion"). Основание выбора: #x4.at("justification")

Начальное разбиение строилось по равным по ширине интервалам от нуля до верхней границы, покрывающей наблюдаемый диапазон и правый хвост распределения.
После этого соседние интервалы объединялись до выполнения условия применимости критерия Пирсона, то есть до достижения ожидаемой частоты не менее 5 в каждом интервале.

Наблюдаемое значение статистики равно #fmt(x4.at("statistic"), digits: 6), число степеней свободы равно #fmt(x4.at("degrees_of_freedom"), digits: 0), p-value = #fmt(x4.at("p_value"), digits: 6).
Итоговый вывод: #x4.at("conclusion")
Ошибка I рода в данном пункте означает: #x4.at("type_i_error")

== Таблица интервалов

#table(
  columns: 5,
  inset: 4pt,
  stroke: .5pt + black,
  table.header(
    [Интервал], [Наблюдаемые], [Ожидаемые], [Вероятность], [Вклад в chi-square],
  ),
  for row in pearson.slice(1) {
    [#row.at(1)]
    [#row.at(4)]
    [#row.at(5)]
    [#row.at(6)]
    [#row.at(7)]
  },
)

#figure(
  image("figures/x4_hist.png", width: 80%),
  caption: [Гистограмма выборки X4],
)

#figure(
  image("figures/x4_pearson_intervals.png", width: 95%),
  caption: [Сравнение наблюдаемых и ожидаемых частот по интервалам для X4],
)

= Итоговый вывод

#data.at("summary")
