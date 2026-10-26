# subset 步骤 · 需求规格

> eca-pp 0.5.6 · 2026-10-26
> 从一个已标准化的产物里挑出一部分细胞，写成同样合格的标准形。
> 起因：tabula-sapiens-clean（28 个器官）是 2026-09-07 用临时脚本从 tabula-sapiens 挑出来的，
> X 里存的是原始计数，result.json 的数字也是从源产物照抄的；identify-columns 的探针做 HVG
> 时对 X 取 expm1，结果溢出成无穷大而崩溃。从一个产物派生另一个产物，以后只走这一步。

## 1. 输入

- `src`：standardize 或 subset 的产物，可以给 OUTDIR，也可以给其中的 `standardized.h5ad`；
  同目录下的 `result.json` 的 status 必须是 `ok` 或 `needs_review`。
- 源产物必须有整数的 `layers["counts"]`（I1），并且 `uns["eca_pp_standardize"]["species"]` 存在。
  任一条件不满足，就以 rejected 结束（exit 2，`rejected_at: input`）。

## 2. 筛选（全部在 `layers["counts"]` 上算，按严格大于判断）

| 参数 | 保留条件 | 被丢弃时记下的理由 |
|---|---|---|
| `--obs COL=VALUE`（可多次） | 每一个给出的列都等于对应的值 | `obs COL!=VALUE` |
| `--min-genes N` | 检出基因数 > N | `n_genes<=N` |
| `--min-counts N` | 总计数 > N | `n_counts<=N` |

一个细胞同时不满足几条时，几条理由用 `+` 连起来。
`OUTDIR/dropped_cells.csv.gz` 有 `cell_id` 和 `reason` 两列，每个被丢弃的细胞占一行。

## 3. 输出契约

与 [standardize](standardize-spec.md) 的标准形相同：

| # | 不变量 | 本步的做法 |
|---|---|---|
| I1 | `layers["counts"]` = 整数原始 counts | 源产物 counts 中被保留的那些行，原样保留 |
| I2 | `X` = log1p(normalize_total(counts, 1e4))，float32 | 在子集上重新计算（`build.finalize`） |
| I3 | 基因空间 | 不变：不丢基因，`var` 原样保留 |
| I4 | 四个 QC 列 | 先删掉源产物的 QC 列，再在子集上重新计算（所以不会产生 `__original` 备份列） |
| I7 | 通过硬 QC 门 | 沿用源产物 result.json 的 `params`（min_cells、min_genes、no_gate）重新检查；不通过则 rejected（`final_gate`），不写 h5ad |
| I8 | 溯源 | `uns["eca_pp_standardize"]` 沿用源产物；`uns["eca_pp_subset"]` 记录本步版本、源和筛选条件；result.json 见 §4 |

- `obsp` 一律去掉：它们是全量产物上建的近邻图，不描述子集。
- `--counts-only`：只留 counts 层，同时去掉 `obsm` 和 `varm`。

## 4. result.json

格式与 standardize 相同（schema_version 2），区别如下：

- `step` 为 `subset`；
- `species` 照抄源产物的值；
- `metrics` 在子集上重新计算：`n_cells`、`n_vars`、`n_genes_detected`、`counts_integer`、`x_normalization`（`stancounts.detect_normalization`）、`qc`；
- `status` 沿用源产物的 `ok` / `needs_review`，源产物的 `reasons` 一并带上；
- 新增 `derived_from` 一项，含：`src`、`src_sha256`、`src_step`、`src_step_version`、`filters`、`n_source_cells`、`n_dropped`。

下游的 ECA-RSI 0.4.18 起认 `step: subset` 的产物，与 standardize 的产物同等对待。

## 5. 输入检查（integration-probe 与 identify-columns）

两个步骤读入数据后，先用 `stancounts.detect_normalization` 检查 X：如果 X 不是 log1p
（比如存的是原始计数），就以 rejected 结束（exit 2），并写明原因，让人先跑 standardize 或 subset。
identify-columns 在调用模型之前就做这项检查。

## 6. CLI

```
eca-pp-subset SRC -o OUTDIR [--obs COL=VALUE ...] [--min-genes N] [--min-counts N] [--counts-only]
bash run.sh subset SRC -o OUTDIR ...
```

退出码与 standardize 相同：0 表示 ok / needs_review，1 表示出错，2 表示 rejected。
