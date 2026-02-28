import sqlite3
import csv
from pathlib import Path
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent
DB_PATH = BASE_DIR / "immune_trial.db"
OUTPUTS_DIR = BASE_DIR / "outputs"
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")

# Allowed deviation from 100.0 when summing per-sample percentages (float rounding)
PCT_SUM_TOLERANCE = 0.01


def get_connection() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def run_frequency_metrics(conn: sqlite3.Connection) -> None:
    """Generate outputs/sample_population_frequencies.csv with one row per (sample, population)."""
    query = """
        SELECT
            cc.sample_id AS sample,
            SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS total_count,
            cc.population,
            cc.count
        FROM cell_counts cc
        ORDER BY cc.sample_id, cc.population
    """
    rows = conn.execute(query).fetchall()

    # --- Build output records and validate before writing ---
    records = []
    pct_sums: dict[str, float] = {}
    pop_counts: dict[str, int] = {}

    for row in rows:
        sample = row["sample"]
        total_count = row["total_count"]
        count = row["count"]

        if not total_count:
            raise RuntimeError(f"Zero total_count for sample {sample}")

        pct = count / total_count * 100.0
        records.append((sample, total_count, row["population"], count, round(pct, 6)))

        pct_sums[sample] = pct_sums.get(sample, 0.0) + pct
        pop_counts[sample] = pop_counts.get(sample, 0) + 1

    # Row count: must equal 5 * number_of_samples
    sample_count = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    expected_rows = sample_count * len(POPULATIONS)
    if len(records) != expected_rows:
        raise RuntimeError(
            f"Frequency metrics row count mismatch: got {len(records)}, expected {expected_rows} "
            f"(samples={sample_count}, populations={len(POPULATIONS)})"
        )

    # Each sample must have exactly 5 population rows
    bad_counts = {s: n for s, n in pop_counts.items() if n != len(POPULATIONS)}
    if bad_counts:
        raise RuntimeError(
            f"Samples with wrong population count (expected {len(POPULATIONS)}): {bad_counts}"
        )

    # Per-sample percentages must sum to ~100
    bad_sums = {s: v for s, v in pct_sums.items() if abs(v - 100.0) > PCT_SUM_TOLERANCE}
    if bad_sums:
        worst = max(bad_sums.items(), key=lambda x: abs(x[1] - 100.0))
        raise RuntimeError(
            f"Percentage sums deviate from 100 beyond tolerance={PCT_SUM_TOLERANCE} "
            f"in {len(bad_sums)} sample(s); worst: {worst[0]}={worst[1]:.6f}"
        )

    # --- All checks passed; write output ---
    out_path = OUTPUTS_DIR / "sample_population_frequencies.csv"
    with out_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "total_count", "population", "count", "percentage"])
        writer.writerows(records)

    print(f"[Frequency metrics] {out_path} written ({len(records)} rows, {sample_count} samples)")


def _bh_correction(p_values: list[float]) -> list[float]:
    """Benjamini–Hochberg FDR correction. Returns q-values in the same order as input."""
    n = len(p_values)
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    q_values = [0.0] * n
    running_min = 1.0
    for rev_rank, (orig_idx, p) in enumerate(reversed(indexed)):
        rank = n - rev_rank  # 1-based rank (descending order)
        q = p * n / rank
        running_min = min(running_min, q)
        q_values[orig_idx] = running_min
    return q_values


def run_responder_analysis(conn: sqlite3.Connection) -> None:
    """Responder vs non-responder differential analysis for melanoma+miraclib+PBMC."""
    # Query percentages for the filtered cohort directly from the DB
    query = """
        SELECT
            sub.response,
            cc.population,
            CAST(cc.count AS REAL) * 100.0 / SUM(cc.count) OVER (PARTITION BY cc.sample_id) AS percentage
        FROM samples sa
        JOIN subjects sub ON sa.subject_id = sub.subject_id
        JOIN cell_counts cc ON sa.sample_id = cc.sample_id
        WHERE sub.condition = 'melanoma'
          AND sub.treatment = 'miraclib'
          AND sa.sample_type = 'PBMC'
          AND sub.response IN ('yes', 'no')
        ORDER BY cc.population, sub.response
    """
    rows = conn.execute(query).fetchall()

    if not rows:
        raise RuntimeError("No rows returned for melanoma+miraclib+PBMC cohort. Check data filters.")

    # Group percentages by (population, response)
    from collections import defaultdict
    groups: dict[str, dict[str, list[float]]] = {
        pop: {"yes": [], "no": []} for pop in POPULATIONS
    }
    for row in rows:
        groups[row["population"]][row["response"]].append(row["percentage"])

    # Run Mann–Whitney U tests
    results = []
    for pop in POPULATIONS:
        yes_vals = groups[pop]["yes"]
        no_vals = groups[pop]["no"]
        n_yes = len(yes_vals)
        n_no = len(no_vals)
        if n_yes == 0 or n_no == 0:
            raise RuntimeError(
                f"Population '{pop}' has empty group: n_yes={n_yes}, n_no={n_no}. "
                "Cannot run Mann–Whitney U test."
            )
        stat, p_val = mannwhitneyu(yes_vals, no_vals, alternative="two-sided")
        results.append({"population": pop, "n_yes": n_yes, "n_no": n_no, "p_value": p_val})

    # Benjamini–Hochberg FDR correction
    p_values = [r["p_value"] for r in results]
    q_values = _bh_correction(p_values)
    for r, q in zip(results, q_values):
        r["q_value"] = q

    # Write CSV
    csv_path = OUTPUTS_DIR / "responder_differential_analysis.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["population", "n_yes", "n_no", "p_value", "q_value"])
        writer.writeheader()
        writer.writerows(results)
    print(f"[Responder analysis] {csv_path} written ({len(results)} populations)")

    # --- Boxplot ---
    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(4 * len(POPULATIONS), 5), sharey=False)
    for ax, pop in zip(axes, POPULATIONS):
        yes_vals = groups[pop]["yes"]
        no_vals = groups[pop]["no"]
        ax.boxplot([yes_vals, no_vals], labels=["Responder\n(yes)", "Non-responder\n(no)"],
                   patch_artist=True,
                   boxprops=dict(facecolor="#a8c8f0"),
                   medianprops=dict(color="black", linewidth=2))
        ax.set_title(pop.replace("_", " ").title(), fontsize=10)
        ax.set_ylabel("Percentage (%)")
        ax.tick_params(axis="x", labelsize=8)

    fig.suptitle("Responder vs Non-responder — Melanoma / Miraclib / PBMC", fontsize=12, y=1.02)
    plt.tight_layout()
    plot_path = OUTPUTS_DIR / "responder_stratification_boxplot.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[Responder analysis] {plot_path} saved")


def main() -> None:
    OUTPUTS_DIR.mkdir(exist_ok=True)

    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found: {DB_PATH}. Run python load_data.py first."
        )

    with get_connection() as conn:
        run_frequency_metrics(conn)
        run_responder_analysis(conn)


if __name__ == "__main__":
    main()