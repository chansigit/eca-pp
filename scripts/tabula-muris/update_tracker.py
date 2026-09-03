"""Update the "scRNA-seq Dataset Tracker" Google Sheet with eca-pp run results.

Rows are keyed on the h5ad_path column: existing rows are updated in place (only the
eca_pp_* / claude_* / status / last_updated / notes columns are touched, other columns
such as eca_rsi_* are preserved), new files are appended.

IO/auth: /scratch/users/chensj16/gsheets/gsheets_client.py (gspread; token in
~/.config/gsheets/token.json, re-auth with gsheets_auth.py on invalid_grant). One read and one
batch write per run; the header is normalised to the canonical column list (duplicates
collapsed, missing columns added, extra columns kept at the end).

    /scratch/users/chensj16/venvs/dl2025/.venv/bin/python scripts/tabula-muris/update_tracker.py \
        --datasets DS [DS ...] [--dry-run] \
        [--skipped DS=reason ...] [--note DS=text ...] [--note DS/SAMPLE=text ...]
Env: ECA_TRACKER_SESSION_DIR / ECA_TRACKER_SESSION_ID identify the Claude session.
"""
import glob, json, os, re, sys, datetime
from collections import Counter

sys.path.insert(0, "/scratch/users/chensj16/gsheets")
from gsheets_client import open_sheet  # noqa: E402

OAK_BASE = "/home/users/chensj16/oak/data/sc"
RUNS = os.path.join(os.environ["SCRATCH"], "eca-pp-runs")
# the Claude session doing the run — override per session via env
SESSION_DIR = os.environ.get("ECA_TRACKER_SESSION_DIR",
    "/home/users/chensj16/.claude/projects/-scratch-users-chensj16-projects-eca-pp")
SESSION_ID = os.environ.get("ECA_TRACKER_SESSION_ID", "4993b76c-6277-42f2-bd70-806aecd1a47f")
HEADER = ["dataset_name", "h5ad_path", "eca_pp_done", "eca_pp_output_dir", "eca_rsi_done",
          "eca_rsi_run_dir", "eca_rsi_release_path", "claude_session_dir", "claude_session_id",
          "status", "last_updated", "notes"]
NOW = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

# per-file caveats worth surfacing in the sheet
CAVEATS = {
    ("tabula-muris-drop", "Mammary_Gland"): "identify-columns cell_type=null (eca-pp bug: agent returned explicit null, cli.py:319 does not fall back); cell_ontology_class is the obvious column",
}
FACS_NOTE = ("caution: 1-63% cells have empty-string tissue/mouse.id/plate.barcode/mouse.sex "
             "(unannotated in TM); profiler treats them as a group, batch verdict may be biased. "
             "Filled with 'missing' in standardized.h5ad")


def status_of(d):
    st = os.path.join(d, "status.txt")
    txt = open(st).read() if os.path.isfile(st) else ""
    g = lambda k: (re.findall(rf"^{k} exit=(\d+)$", txt, re.M) or [None])[-1]
    return g("standardize"), g("identify-columns"), g("fill-missing")


def build_rows(datasets, skipped, extra_notes=()):
    """extra_notes: [(dataset, sample_or_"*", text)] appended to matching rows."""
    rows = {}
    for d in sorted(glob.glob(f"{RUNS}/*/*/")):
        ds = os.path.basename(os.path.dirname(d.rstrip("/")))
        if ds not in datasets:
            continue
        organ = os.path.basename(d.rstrip("/"))
        s1, s2, s3 = status_of(d)
        if s1 is None:
            continue
        r1 = r2 = None
        try: r1 = json.load(open(os.path.join(d, "standardize/result.json")))
        except Exception: pass
        try: r2 = json.load(open(os.path.join(d, "identify_columns/result.json")))
        except Exception: pass
        src = (r1 or {}).get("src") or ""
        if not src:
            continue
        done = s1 == "0" and s2 == "0"
        notes = []
        std_status = (r1 or {}).get("status")
        if r1:
            m = r1["metrics"]
            h = m.get("harmonization") or {}
            notes.append(f"cells={m.get('n_cells')} species={(r1.get('species') or {}).get('code')} "
                         f"counts={(r1.get('counts') or {}).get('source') or m.get('counts_source')} "
                         f"genes_kept={h.get('genes_kept')} dropped_frac={h.get('genes_dropped_frac')}")
            if r1.get("reasons"):
                notes.append(f"standardize {std_status}: " + " | ".join(r1["reasons"]))
        if r2:
            b = (r2.get("columns") or {}).get("batch")
            ct = (r2.get("columns") or {}).get("cell_type")
            llm = r2["metrics"].get("llm") or {}
            notes.append(f"batch={(b or {}).get('value')} correction={(b or {}).get('correction')} "
                         f"cell_type={(ct or {}).get('value')} probes={len(r2.get('trials') or [])} "
                         f"agent=${llm.get('cost_usd', 0):.2f} ({','.join(llm.get('models') or [])})")
        if s3 is not None:
            notes.append("empty batch strings -> 'missing' (in job)")
        if ds == "tabula-muris-facs":
            notes.append(FACS_NOTE)
        if (ds, organ) in CAVEATS:
            notes.append(CAVEATS[(ds, organ)])
        if ds in ("tabula-muris-drop", "tabula-muris-facs"):
            notes.append("empty batch strings -> 'missing' (post-hoc fill_missing.py)")
        for nds, nsample, text in extra_notes:
            if nds == ds and nsample in ("*", organ):
                notes.append(text)
        if done:
            status = "eca_pp ok" if std_status == "ok" else f"eca_pp ok (standardize {std_status})"
        else:
            status = f"eca_pp failed (standardize exit {s1}, identify-columns exit {s2})"
            if r1 and r1.get("reasons"):
                notes.append("failure: " + " | ".join(r1["reasons"]))
        rows[src] = {
            "dataset_name": ds, "h5ad_path": src, "eca_pp_done": "TRUE" if done else "FALSE",
            "eca_pp_output_dir": f"{OAK_BASE}/{ds}/eca-pp/{organ}",
            "claude_session_dir": SESSION_DIR, "claude_session_id": SESSION_ID,
            "status": status, "last_updated": NOW, "notes": "; ".join(notes),
        }
    for ds, why in skipped:
        for f in sorted(glob.glob(f"{OAK_BASE}/{ds}/*.h5ad")):
            rows[f] = {
                "dataset_name": ds, "h5ad_path": f, "eca_pp_done": "FALSE", "eca_pp_output_dir": "",
                "claude_session_dir": SESSION_DIR, "claude_session_id": SESSION_ID,
                "status": "eca_pp skipped", "last_updated": NOW, "notes": why,
            }
    return rows


def main():
    dry = "--dry-run" in sys.argv
    def after(flag):
        if flag not in sys.argv: return []
        out = []
        for a in sys.argv[sys.argv.index(flag) + 1:]:
            if a.startswith("--"): break
            out.append(a)
        return out
    datasets = after("--datasets")
    skipped = [tuple(spec.partition("=")[::2]) for spec in after("--skipped")]
    # --note DS=text  or  --note DS/SAMPLE=text
    extra_notes = []
    for spec in after("--note"):
        target, _, text = spec.partition("=")
        nds, _, nsample = target.partition("/")
        extra_notes.append((nds, nsample or "*", text))
    if not datasets and not skipped:
        sys.exit("usage: --datasets DS [DS ...] [--skipped DS=reason ...] [--dry-run]")
    ws = open_sheet()
    raw_rows = ws.get_all_values()
    raw_header = raw_rows[0] if raw_rows else list(HEADER)
    # Normalise the header: canonical columns first (in canonical order), then any
    # genuinely extra columns; duplicated names (e.g. after an accidental drag-fill in
    # the UI) are collapsed and reported.
    dups = [h for i, h in enumerate(raw_header) if h in raw_header[:i]]
    extra = [h for h in dict.fromkeys(raw_header) if h not in HEADER]
    header = list(HEADER) + extra
    if dups or [h for h in HEADER if h not in raw_header]:
        print(f"WARNING sheet header was {raw_header}; duplicates={dups}; normalising to {header}",
              file=sys.stderr)
    existing = []
    for vals in raw_rows[1:]:
        if not any(v.strip() for v in vals):
            continue
        row = {}
        for h, v in zip(raw_header, vals):
            if h in row and h in dups:
                continue  # keep the first occurrence of a duplicated column
            row[h] = v
        existing.append(row)
    new = build_rows(set(datasets), skipped, extra_notes)
    out, seen, n_upd = [], set(), 0
    for row in existing:
        key = (row.get("h5ad_path") or "").strip()
        if key in new:
            row.update(new[key]); seen.add(key); n_upd += 1
        out.append(row)
    n_add = 0
    for key, row in new.items():
        if key not in seen:
            out.append(row); n_add += 1
    table = [header] + [[str(row.get(h, "")) for h in header] for row in out]
    print("rows per dataset:", dict(Counter(r["dataset_name"] for r in new.values())))
    print(f"existing rows={len(existing)} updated={n_upd} appended={n_add} total={len(out)}")
    if dry:
        for r in table[:3]: print(r)
        return
    ws.update(range_name="A1", values=table, value_input_option="RAW")
    if len(raw_rows) > len(table):   # rows removed: clear what is left below the new table
        ws.batch_clear([f"A{len(table)+1}:{chr(64+len(header))}{len(raw_rows)}"])
    back = ws.get_all_values()
    print(f"verified: sheet now has {len(back)-1} data rows, header ok={back[0]==header}")


if __name__ == "__main__":
    main()
