from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


BASE_DIR = Path(__file__).resolve().parent.parent
DATA_PATH = BASE_DIR / "resources" / "RGR2_A-1_X1-X4.csv"
FIGURES_DIR = BASE_DIR / "figures"
ARTIFACTS_DIR = BASE_DIR / "artifacts"
ALPHA = 0.05
MU0_X3 = 62.55
LAMBDA_X4 = 0.085
EXPECTED_COLUMNS = ["X1", "X2", "X3", "X4"]
UNITS = "мс"
VARIANT = "A-1"
GROUP = "A"
PEARSON_MIN_EXPECTED = 5.0


@dataclass
class TestResult:
    name: str
    null_hypothesis: str
    alternative_hypothesis: str
    criterion: str
    justification: str
    statistic_name: str
    statistic: float
    statistic_scipy: float | None
    p_value: float
    p_value_scipy: float | None
    degrees_of_freedom: float | None
    critical_value: float | None
    alpha: float
    reject_null: bool
    conclusion: str
    type_i_error: str
    details: dict


def ensure_output_dirs() -> None:
    FIGURES_DIR.mkdir(exist_ok=True)
    ARTIFACTS_DIR.mkdir(exist_ok=True)


def load_dataset() -> pd.DataFrame:
    df = pd.read_csv(DATA_PATH, sep=";")
    missing = [column for column in EXPECTED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"В CSV отсутствуют обязательные столбцы: {missing}")
    numeric_df = df[EXPECTED_COLUMNS].apply(pd.to_numeric, errors="raise")
    if numeric_df.isnull().any().any():
        raise ValueError("Обнаружены пропуски в обязательных столбцах.")
    return numeric_df


def build_descriptive_statistics(df: pd.DataFrame) -> pd.DataFrame:
    stats_df = pd.DataFrame(
        {
            "count": df.count(),
            "mean": df.mean(),
            "std": df.std(ddof=1),
            "min": df.min(),
            "q25": df.quantile(0.25),
            "median": df.median(),
            "q75": df.quantile(0.75),
            "max": df.max(),
        }
    )
    stats_df.index.name = "variable"
    return stats_df.round(4)


def save_descriptive_statistics(stats_df: pd.DataFrame) -> None:
    stats_df.to_csv(ARTIFACTS_DIR / "descriptive_statistics.csv", encoding="utf-8-sig")


def perform_two_sample_t_test(x1: pd.Series, x2: pd.Series) -> TestResult:
    n1 = len(x1)
    n2 = len(x2)
    mean1 = float(x1.mean())
    mean2 = float(x2.mean())
    var1 = float(x1.var(ddof=1))
    var2 = float(x2.var(ddof=1))
    pooled_variance = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
    standard_error = np.sqrt(pooled_variance * (1 / n1 + 1 / n2))
    statistic = (mean1 - mean2) / standard_error
    degrees_of_freedom = n1 + n2 - 2
    p_value = 2 * stats.t.sf(abs(statistic), df=degrees_of_freedom)
    critical_value = stats.t.ppf(1 - ALPHA / 2, df=degrees_of_freedom)
    scipy_result = stats.ttest_ind(x1, x2, equal_var=True, alternative="two-sided")
    reject_null = bool(p_value < ALPHA)
    conclusion = (
        "На уровне значимости 0.05 гипотеза о равенстве математических ожиданий отвергается."
        if reject_null
        else "На уровне значимости 0.05 нет оснований отвергать гипотезу о равенстве математических ожиданий."
    )
    return TestResult(
        name="Параметрический тест для X1 и X2",
        null_hypothesis="H0: mu_X1 = mu_X2",
        alternative_hypothesis="H1: mu_X1 != mu_X2",
        criterion="Двухвыборочный t-критерий Стьюдента для независимых выборок",
        justification=(
            "По условию рассматриваются две независимые нормальные совокупности. "
            "Использован классический t-критерий с объединённой оценкой дисперсии."
        ),
        statistic_name="t",
        statistic=float(statistic),
        statistic_scipy=float(scipy_result.statistic),
        p_value=float(p_value),
        p_value_scipy=float(scipy_result.pvalue),
        degrees_of_freedom=float(degrees_of_freedom),
        critical_value=float(critical_value),
        alpha=ALPHA,
        reject_null=reject_null,
        conclusion=conclusion,
        type_i_error=(
            "Ошибка I рода: отклонить H0 и заявить о различии математических ожиданий, "
            "хотя в генеральных совокупностях средние на самом деле равны."
        ),
        details={
            "n1": n1,
            "n2": n2,
            "mean1": mean1,
            "mean2": mean2,
            "variance1": var1,
            "variance2": var2,
            "pooled_variance": float(pooled_variance),
            "standard_error": float(standard_error),
            "mean_difference": float(mean1 - mean2),
        },
    )


def perform_one_sample_t_test(x3: pd.Series) -> TestResult:
    n = len(x3)
    mean_value = float(x3.mean())
    std_value = float(x3.std(ddof=1))
    standard_error = std_value / np.sqrt(n)
    statistic = (mean_value - MU0_X3) / standard_error
    degrees_of_freedom = n - 1
    p_value = 2 * stats.t.sf(abs(statistic), df=degrees_of_freedom)
    critical_value = stats.t.ppf(1 - ALPHA / 2, df=degrees_of_freedom)
    scipy_result = stats.ttest_1samp(x3, popmean=MU0_X3, alternative="two-sided")
    reject_null = bool(p_value < ALPHA)
    conclusion = (
        "На уровне значимости 0.05 гипотеза H0: mu = 62.55 отвергается."
        if reject_null
        else "На уровне значимости 0.05 нет оснований отвергать гипотезу H0: mu = 62.55."
    )
    return TestResult(
        name="Одновыборочный тест для X3",
        null_hypothesis="H0: mu_X3 = 62.55",
        alternative_hypothesis="H1: mu_X3 != 62.55",
        criterion="Одновыборочный t-критерий Стьюдента",
        justification=(
            "По условию используется критерий для среднего нормального распределения "
            "при неизвестной дисперсии."
        ),
        statistic_name="t",
        statistic=float(statistic),
        statistic_scipy=float(scipy_result.statistic),
        p_value=float(p_value),
        p_value_scipy=float(scipy_result.pvalue),
        degrees_of_freedom=float(degrees_of_freedom),
        critical_value=float(critical_value),
        alpha=ALPHA,
        reject_null=reject_null,
        conclusion=conclusion,
        type_i_error=(
            "Ошибка I рода: отклонить H0 и сделать вывод, что математическое ожидание X3 "
            "отличается от 62.55, хотя на самом деле mu = 62.55."
        ),
        details={
            "n": n,
            "sample_mean": mean_value,
            "sample_std": std_value,
            "standard_error": float(standard_error),
            "mu0": MU0_X3,
        },
    )


def perform_mann_whitney_test(x1: pd.Series, x2: pd.Series, t_test_result: TestResult) -> TestResult:
    scipy_result = stats.mannwhitneyu(x1, x2, alternative="two-sided", method="asymptotic")
    reject_null = bool(scipy_result.pvalue < ALPHA)
    if reject_null == t_test_result.reject_null:
        comparison = "Вывод совпадает с параметрическим t-критерием."
    else:
        comparison = (
            "Вывод отличается от параметрического t-критерия, что может быть связано "
            "с разной чувствительностью критериев к форме распределения и выбросам."
        )
    return TestResult(
        name="Критерий Манна–Уитни для X1 и X2",
        null_hypothesis="H0: распределения X1 и X2 имеют одинаковое положение",
        alternative_hypothesis="H1: распределения X1 и X2 различаются по положению",
        criterion="Непараметрический критерий Манна–Уитни",
        justification=(
            "Критерий использован как непараметрическая альтернатива для сравнения "
            "двух независимых выборок без жёсткой опоры на нормальность."
        ),
        statistic_name="U",
        statistic=float(scipy_result.statistic),
        statistic_scipy=float(scipy_result.statistic),
        p_value=float(scipy_result.pvalue),
        p_value_scipy=float(scipy_result.pvalue),
        degrees_of_freedom=None,
        critical_value=None,
        alpha=ALPHA,
        reject_null=reject_null,
        conclusion=(
            "На уровне значимости 0.05 нулевая гипотеза отвергается."
            if reject_null
            else "На уровне значимости 0.05 нет оснований отвергать нулевую гипотезу."
        ),
        type_i_error=(
            "Ошибка I рода: заявить о различии положений распределений X1 и X2, "
            "хотя в действительности различий нет."
        ),
        details={
            "comparison_with_t_test": comparison,
            "t_test_reject_null": t_test_result.reject_null,
        },
    )


def _initial_exponential_bins(x4: pd.Series) -> np.ndarray:
    bin_count = int(np.ceil(np.log2(len(x4)) + 1))
    finite_upper = float(max(x4.max(), stats.expon.ppf(0.995, scale=1 / LAMBDA_X4)))
    edges = np.linspace(0.0, finite_upper, bin_count)
    return np.concatenate([edges, [np.inf]])


def _merge_small_expected_bins(observed: np.ndarray, expected: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    observed = observed.astype(float)
    expected = expected.astype(float)
    edges = edges.astype(float)
    while len(expected) > 1 and np.any(expected < PEARSON_MIN_EXPECTED):
        idx = int(np.where(expected < PEARSON_MIN_EXPECTED)[0][0])
        if idx == 0:
            merge_with = 1
        elif idx == len(expected) - 1:
            merge_with = idx - 1
        else:
            merge_with = idx - 1 if expected[idx - 1] <= expected[idx + 1] else idx + 1
        left = min(idx, merge_with)
        right = max(idx, merge_with)
        observed[left] = observed[left] + observed[right]
        expected[left] = expected[left] + expected[right]
        observed = np.delete(observed, right)
        expected = np.delete(expected, right)
        edges = np.delete(edges, right)
    return observed, expected, edges


def perform_pearson_chi_square_test(x4: pd.Series) -> tuple[TestResult, pd.DataFrame]:
    edges = _initial_exponential_bins(x4)
    observed, _ = np.histogram(x4, bins=edges)
    cdf_values = stats.expon.cdf(edges, scale=1 / LAMBDA_X4)
    probabilities = np.diff(cdf_values)
    expected = len(x4) * probabilities
    observed_merged, expected_merged, edges_merged = _merge_small_expected_bins(observed, expected, edges)
    chi_components = (observed_merged - expected_merged) ** 2 / expected_merged
    statistic = float(np.sum(chi_components))
    degrees_of_freedom = len(observed_merged) - 1
    p_value = float(stats.chi2.sf(statistic, df=degrees_of_freedom))
    critical_value = float(stats.chi2.ppf(1 - ALPHA, df=degrees_of_freedom))
    reject_null = bool(p_value < ALPHA)

    rows = []
    for idx, (left, right, obs, exp) in enumerate(
        zip(edges_merged[:-1], edges_merged[1:], observed_merged, expected_merged),
        start=1,
    ):
        probability = float(exp / len(x4))
        interval = (
            f"[{left:.3f}; +inf)"
            if np.isinf(right)
            else f"[{left:.3f}; {right:.3f})"
        )
        rows.append(
            {
                "interval_index": idx,
                "interval": interval,
                "left": round(float(left), 6),
                "right": None if np.isinf(right) else round(float(right), 6),
                "observed": int(obs),
                "expected": round(float(exp), 6),
                "probability": round(probability, 6),
                "chi_component": round(float((obs - exp) ** 2 / exp), 6),
            }
        )
    table_df = pd.DataFrame(rows)
    table_df.to_csv(ARTIFACTS_DIR / "x4_pearson_intervals.csv", index=False, encoding="utf-8-sig")

    result = TestResult(
        name="Критерий согласия Пирсона для X4",
        null_hypothesis="H0: X4 имеет показательное распределение с lambda = 0.085",
        alternative_hypothesis="H1: распределение X4 отличается от Exp(0.085)",
        criterion="Критерий согласия Пирсона chi-square",
        justification=(
            "Параметр lambda задан заранее и не оценивается по данным, поэтому число "
            "степеней свободы определяется как k - 1 после объединения интервалов."
        ),
        statistic_name="chi-square",
        statistic=statistic,
        statistic_scipy=None,
        p_value=p_value,
        p_value_scipy=None,
        degrees_of_freedom=float(degrees_of_freedom),
        critical_value=critical_value,
        alpha=ALPHA,
        reject_null=reject_null,
        conclusion=(
            "На уровне значимости 0.05 гипотеза о показательном распределении отвергается."
            if reject_null
            else "На уровне значимости 0.05 нет оснований отвергать гипотезу о показательном распределении."
        ),
        type_i_error=(
            "Ошибка I рода: отвергнуть гипотезу о показательном распределении Exp(0.085), "
            "хотя на самом деле выборка X4 получена из такого распределения."
        ),
        details={
            "lambda": LAMBDA_X4,
            "bin_count_after_merge": int(len(observed_merged)),
            "initial_bin_count": int(len(observed)),
            "merged_bins": int(len(observed) - len(observed_merged)),
            "chi_square_components": [float(value) for value in chi_components],
        },
    )
    return result, table_df


def save_plot(fig: plt.Figure, filename: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)


def create_plots(df: pd.DataFrame, pearson_table: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["X1"], bins="sturges", color="#4C72B0", edgecolor="black", alpha=0.85)
    ax.set_title("Гистограмма X1")
    ax.set_xlabel(f"X1, {UNITS}")
    ax.set_ylabel("Частота")
    save_plot(fig, "x1_hist.png")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["X2"], bins="sturges", color="#55A868", edgecolor="black", alpha=0.85)
    ax.set_title("Гистограмма X2")
    ax.set_xlabel(f"X2, {UNITS}")
    ax.set_ylabel("Частота")
    save_plot(fig, "x2_hist.png")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["X1"], bins="sturges", alpha=0.55, label="X1", color="#4C72B0", edgecolor="black")
    ax.hist(df["X2"], bins="sturges", alpha=0.55, label="X2", color="#55A868", edgecolor="black")
    ax.set_title("Сравнительные гистограммы X1 и X2")
    ax.set_xlabel(f"Значение, {UNITS}")
    ax.set_ylabel("Частота")
    ax.legend()
    save_plot(fig, "x1_x2_comparison_hist.png")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.boxplot([df["X1"], df["X2"]], labels=["X1", "X2"], patch_artist=True)
    ax.set_title("Boxplot для X1 и X2")
    ax.set_ylabel(f"Значение, {UNITS}")
    save_plot(fig, "x1_x2_boxplot.png")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["X3"], bins="sturges", color="#C44E52", edgecolor="black", alpha=0.85)
    ax.axvline(MU0_X3, color="black", linestyle="--", label=f"mu0 = {MU0_X3}")
    ax.set_title("Гистограмма X3")
    ax.set_xlabel(f"X3, {UNITS}")
    ax.set_ylabel("Частота")
    ax.legend()
    save_plot(fig, "x3_hist.png")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["X4"], bins="sturges", color="#8172B2", edgecolor="black", alpha=0.85)
    ax.set_title("Гистограмма X4")
    ax.set_xlabel(f"X4, {UNITS}")
    ax.set_ylabel("Частота")
    save_plot(fig, "x4_hist.png")

    fig, ax = plt.subplots(figsize=(10, 5))
    x_positions = np.arange(len(pearson_table))
    ax.bar(x_positions - 0.2, pearson_table["observed"], width=0.4, label="Наблюдаемые", color="#4C72B0")
    ax.bar(x_positions + 0.2, pearson_table["expected"], width=0.4, label="Ожидаемые", color="#DD8452")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(pearson_table["interval"], rotation=35, ha="right")
    ax.set_title("Интервалы критерия Пирсона для X4")
    ax.set_xlabel("Интервалы")
    ax.set_ylabel("Частоты")
    ax.legend()
    save_plot(fig, "x4_pearson_intervals.png")


def build_summary(
    t_test_result: TestResult,
    x3_result: TestResult,
    mann_whitney_result: TestResult,
    pearson_result: TestResult,
) -> str:
    verdict = lambda result: "H0 отвергается" if result.reject_null else "нет оснований отвергать H0"
    lines = [
        "В работе были проверены четыре статистические гипотезы по выборкам X1, X2, X3 и X4.",
        "Для X1 и X2 использованы параметрический двухвыборочный t-критерий Стьюдента и непараметрический критерий Манна–Уитни.",
        f"По t-критерию для X1 и X2: {verdict(t_test_result)} (p-value = {t_test_result.p_value:.4f}).",
        f"По критерию Манна–Уитни для X1 и X2: {verdict(mann_whitney_result)} (p-value = {mann_whitney_result.p_value:.4f}).",
        f"Для X3 применён одновыборочный t-критерий Стьюдента: {verdict(x3_result)} (p-value = {x3_result.p_value:.4f}).",
        f"Для X4 применён критерий согласия Пирсона: {verdict(pearson_result)} (p-value = {pearson_result.p_value:.4f}).",
        "Содержательно результаты показывают, насколько выборочные данные согласуются с заданными предположениями о средних значениях и виде распределения.",
        "Все выводы сделаны при уровне значимости alpha = 0.05.",
    ]
    return "\n".join(lines)


def print_console_report(
    descriptive_stats: pd.DataFrame,
    results: list[TestResult],
    pearson_table: pd.DataFrame,
    sample_size: int,
) -> None:
    print(f"Вариант: {VARIANT}")
    print(f"Группа: {GROUP}")
    print(f"Размер выборки по файлу: {sample_size}")
    print(f"Уровень значимости alpha = {ALPHA}")
    print()
    print("Описательная статистика:")
    print(descriptive_stats.to_string())
    print()
    for result in results:
        print(result.name)
        print(f"  {result.null_hypothesis}")
        print(f"  {result.alternative_hypothesis}")
        print(f"  Критерий: {result.criterion}")
        print(f"  Обоснование: {result.justification}")
        print(
            f"  {result.statistic_name} = {result.statistic:.6f}"
            + (
                f" (scipy: {result.statistic_scipy:.6f})"
                if result.statistic_scipy is not None
                else ""
            )
        )
        if result.degrees_of_freedom is not None:
            print(f"  Степени свободы: {result.degrees_of_freedom:.4f}")
        if result.critical_value is not None:
            print(f"  Критическое значение: {result.critical_value:.6f}")
        print(
            f"  p-value = {result.p_value:.6f}"
            + (
                f" (scipy: {result.p_value_scipy:.6f})"
                if result.p_value_scipy is not None
                else ""
            )
        )
        print(f"  Вывод: {result.conclusion}")
        print(f"  Ошибка I рода: {result.type_i_error}")
        if result.details:
            print("  Детали:")
            for key, value in result.details.items():
                print(f"    {key}: {value}")
        print()
    print("Таблица интервалов критерия Пирсона:")
    print(pearson_table.to_string(index=False))


def save_results_json(
    sample_size: int,
    descriptive_stats: pd.DataFrame,
    results: list[TestResult],
    summary_text: str,
) -> None:
    payload = {
        "meta": {
            "variant": VARIANT,
            "group": GROUP,
            "sample_size": sample_size,
            "alpha": ALPHA,
            "mu0_x3": MU0_X3,
            "lambda_x4": LAMBDA_X4,
            "units": UNITS,
            "data_path": str(DATA_PATH.relative_to(BASE_DIR)),
        },
        "descriptive_statistics": descriptive_stats.reset_index().to_dict(orient="records"),
        "tests": {result.name: asdict(result) for result in results},
        "summary": summary_text,
    }
    (ARTIFACTS_DIR / "results.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def main() -> None:
    ensure_output_dirs()
    df = load_dataset()
    sample_size = len(df)
    descriptive_stats = build_descriptive_statistics(df)
    save_descriptive_statistics(descriptive_stats)

    t_test_result = perform_two_sample_t_test(df["X1"], df["X2"])
    x3_result = perform_one_sample_t_test(df["X3"])
    mann_whitney_result = perform_mann_whitney_test(df["X1"], df["X2"], t_test_result)
    pearson_result, pearson_table = perform_pearson_chi_square_test(df["X4"])

    create_plots(df, pearson_table)

    results = [t_test_result, x3_result, mann_whitney_result, pearson_result]
    summary_text = build_summary(t_test_result, x3_result, mann_whitney_result, pearson_result)
    save_results_json(sample_size, descriptive_stats, results, summary_text)
    print_console_report(descriptive_stats, results, pearson_table, sample_size)


if __name__ == "__main__":
    main()
