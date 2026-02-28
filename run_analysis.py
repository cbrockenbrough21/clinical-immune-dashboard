import sqlite3
import csv
from pathlib import Path

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
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sample", "total_count", "population", "count", "percentage"])
        writer.writerows(records)

    print(f"[Frequency metrics] {out_path} written ({len(records)} rows, {sample_count} samples)")


def main() -> None:
    OUTPUTS_DIR.mkdir(exist_ok=True)

    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found: {DB_PATH}. Run python load_data.py first."
        )

    with get_connection() as conn:
        run_frequency_metrics(conn)


if __name__ == "__main__":
    main()