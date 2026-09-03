"""Upsert eca-pp run results for a dataset into the scRNA-seq Dataset Tracker sheet,
one row per h5ad, matched on the h5ad_path column (gsheets_client, authorized).

    ECA_TRACKER_SESSION_DIR=... ECA_TRACKER_SESSION_ID=... \
    /scratch/users/chensj16/venvs/dl2025/.venv/bin/python scripts/tome-mouse/update_tracker.py \
        --datasets tome-mouse [--note DS=text] [--note DS/SAMPLE=text] [--done-only] [--dry-run]

Row content (notes, status, eca_pp_done, ...) is built by scripts/tabula-muris/update_tracker.py's
build_rows(); this script only swaps the write path for gsheets_client.upsert_row, which also
stamps last_updated.
"""
import os, sys
sys.path.insert(0, "/scratch/users/chensj16/gsheets")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tabula-muris"))
from gsheets_client import open_sheet
import update_tracker as base   # tabula-muris version: build_rows + arg parsing helpers

WRITE_COLS = ["dataset_name", "eca_pp_done", "eca_pp_output_dir", "claude_session_dir",
              "claude_session_id", "status", "notes"]


def main():
    argv = sys.argv
    dry = "--dry-run" in argv
    def after(flag):
        if flag not in argv: return []
        out = []
        for a in argv[argv.index(flag) + 1:]:
            if a.startswith("--"): break
            out.append(a)
        return out
    datasets = after("--datasets")
    skipped = [tuple(s.partition("=")[::2]) for s in after("--skipped")]
    notes = []
    for spec in after("--note"):
        target, _, text = spec.partition("=")
        nds, _, nsample = target.partition("/")
        notes.append((nds, nsample or "*", text))
    if not datasets and not skipped:
        sys.exit("usage: --datasets DS [DS ...] [--skipped DS=reason ...] [--note DS[/SAMPLE]=text ...] [--dry-run]")

    rows = base.build_rows(set(datasets), skipped, notes)
    if "--done-only" in argv:  # jobs still running would otherwise be written as failed
        rows = {k: r for k, r in rows.items() if r["eca_pp_done"] == "TRUE"}

    # One read, then at most two writes (Sheets API: 60 reads/min/user — per-row
    # upsert_row blew through that at ~3 reads per row).
    import datetime
    import gspread
    ws = open_sheet()
    vals = ws.get_all_values()
    h = vals[0]
    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    col = {c: j for j, c in enumerate(h)}
    by_path = {}
    for i, r in enumerate(vals[1:], start=2):
        r = r + [""] * (len(h) - len(r))
        by_path.setdefault(r[col["h5ad_path"]].strip(), i)
    last = max((i for i, r in enumerate(vals, 1) if any(c.strip() for c in r)), default=1)
    cells, block = [], []
    for key, row in rows.items():
        upd = {c: row.get(c, "") for c in WRITE_COLS}
        upd["last_updated"] = stamp
        tag = "update" if key in by_path else "append"
        print(f"{tag:6s} {row['status']:12s} {key.split('/')[-1]:28s} {upd['notes'][:100]}")
        if key in by_path:
            i = by_path[key]
            cells += [gspread.Cell(i, col[c] + 1, str(v)) for c, v in upd.items() if c in col]
        else:
            full = {**upd, "h5ad_path": key}
            block.append([str(full.get(c, "")) for c in h])
    print(f"sheet data rows={last - 1}; updating {len({c.row for c in cells})} rows, "
          f"appending {len(block)} rows at row {last + 1}")
    if dry:
        print("(dry run, nothing written)"); return
    if cells:
        ws.update_cells(cells, value_input_option="USER_ENTERED")
    if block:
        if last + len(block) > ws.row_count:
            ws.add_rows(last + len(block) - ws.row_count)
        ws.update(range_name=f"A{last + 1}", values=block, value_input_option="USER_ENTERED")


if __name__ == "__main__":
    main()
