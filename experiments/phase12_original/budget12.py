"""Phase 12 runtime budget from the pilot (results/pilot/*.jsonl): mean CPU s per task x confirmatory task count / 4."""
import json
from pathlib import Path

from experiments.phase12_original import run12 as R

HERE = Path(__file__).resolve().parent
CONF_COUNTS = {"b1": 20, "b2": 10 * 20 * 3, "b4a": 400, "b4v": 4 * 210, "b4s": 900, "b4t": 439, "b4w": 988,
               "pm97": 33, "pm97nl": 33}
# pilot -> confirmatory scale per task (pilot PM97 blocks are 2,000 beats; confirmatory 20,000)
SCALE = {"pm97": 10.0, "pm97nl": 10.0}


def main():
    rows, tot = [], 0.0
    for part in R.TIER1B + R.TIER1C:
        R_ = [json.loads(line) for line in open(HERE / "results" / "pilot" / f"{part}.jsonl")]
        ok = [r for r in R_ if not r.get("error") and "excluded" not in r]
        cpu = sum(r["t_cpu"] for r in ok) / max(1, len(ok)) * SCALE.get(part, 1.0)
        wall = cpu * CONF_COUNTS[part] / 4 / 3600
        tot += wall
        rows.append(f"| {part} | {len(ok)} | {cpu:.2f} | {CONF_COUNTS[part]} | {wall:.2f} h |")
    md = ["| part | pilot tasks | CPU s / task (scaled) | confirmatory tasks | projected wall (4 workers) |",
          "|---|---|---|---|---|", *rows, f"| **total** | | | | **{tot:.2f} h** |"]
    (HERE / "results" / "pilot" / "budget.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
