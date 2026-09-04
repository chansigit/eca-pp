# eca-pp-standardize 教程(Tabula Muris 实战)

把一个来源不明的单样本 `.h5ad` 变成下游可信赖的标准形。本教程所有输出都来自
真实运行:`data/tabula_muris_droplet/`(小鼠 12 器官,软链接,已 gitignore)。

## 1. 跑一个真样本

在**计算节点**上:

```bash
cd /scratch/users/chensj16/projects/eca-pp
bash run.sh standardize data/tabula_muris_droplet/Heart_and_Aorta.h5ad \
    -o data/out/Heart_and_Aorta
echo $?        # → 0
```

约 20 秒,产出两个文件:

- `data/out/Heart_and_Aorta/standardized.h5ad` — 标准形;
- `data/out/Heart_and_Aorta/result.json` — 判定全记录(失败也会写,永远先看它)。

## 2. 读 result.json(真实内容)

```
status : ok                    ← reasons 为空;存疑点才会变 needs_review
species: mouse | inferred | 0.9 | rule=symbol_overlap
         ← 该数据没有 Ensembl ID 列,是靠"与内置参考的基因清单交集率"判出小鼠的
counts : source=X              ← X 本身就是整数 counts,直接采用
cells  : 624
harmonization:                       ← 键名自带量纲:数的全是基因,不是细胞
  genes_kept=23143  genes_dropped=290 (genes_dropped_frac=1.2%)
    unmapped=183  ambiguous=15  non_gene_feature=92   ← 那 92 个是 ERCC spike-in
  注:standardize 只裁剪基因维度;细胞要么整样本过门、要么整样本被拒(exit 2),
      永远不会部分删除细胞
qc:
  n_mt_genes=0  n_hb_genes=10  median_pct_counts_hb=0.0 ...
  ← n_mt_genes=0 已实查:输入 23433 个基因名里 0 个 mt- 基因——是数据发布方
    在上游就把线粒体基因剔掉了(stangene 小鼠参考含全部 37 个 mt- 基因,
    F4 也没丢过任何 mt 基因)。这是正常情况,**不触发 needs_review**;
    正确解读是"pct_counts_mt 在此数据集上不可用",判断依据就在 metrics.qc
```

**每步耗时**也在 result.json 里(`metrics.timings`,Marrow 样本真实值):

| 步骤 | 耗时 | 占比 |
|---|---|---|
| f1_validate(校验+预门控) | 0.00s | — |
| load(读 h5ad) | 0.58s | 3% |
| f2_counts(counts 定位) | 0.06s | — |
| f3_gates(QC 硬门) | 0.06s | — |
| f4a_species(物种推断) | 0.30s | 1% |
| **f4_harmonize(23k 基因名统一)** | **16.38s** | **79%** |
| f5_qc(QC 列计算) | 0.04s | — |
| f7_build_write(lognorm + 写盘) | 3.27s | 16% |
| **total** | **20.70s** | |

结论:耗时大头是基因名统一(和基因数成正比),其余全部亚秒级。
`--species mm` 只省 0.3 秒——指定物种的意义是**可复现**,不是提速。

## 3. 退出码 → 下一步动作

| 码 | 含义 | 动作 |
|---|---|---|
| 0 | 成功(status=`needs_review` 时附带存疑点) | 用输出;看 `reasons` |
| 2 | 数据永久性问题(细胞/基因太少、counts 不可恢复) | 放弃该样本,重试无用 |
| 3 | **等你拍板** | 看证据 → 补参数重跑 |
| 1 | 意外错误 | 原样重跑一次 |

**真实的 exit 3 长这样**(故意指错 counts 层):

```bash
bash run.sh standardize data/tabula_muris_droplet/Liver.h5ad \
    -o data/out/_demo_blocked --counts-layer nope
echo $?        # → 3
# result.json → reasons:
#   "--counts-layer 'nope' not found; layers present: []"
```

处置:按 reasons 提示补对参数重跑(这里 Liver 的 counts 就在 X,去掉
`--counts-layer` 即可)。物种判不出的 exit 3 同理:看
`species.evidence`(前缀计票 / 命名惯例 / 各物种交集率)→ `--species CODE` 重跑。

## 4. 常用参数

```bash
--species mm          # 指定物种(hs/mm/rn/dr/dm/ce/cyno/rhesus/marmoset/lemur)
--counts-layer NAME   # 指定 counts 所在 layer
--keep-unmapped       # 保留映射不上的基因(默认丢弃;上面 290 个会被留下、维持原名)
--llm                 # 物种推断失败时允许一次 LLM 兜底(默认关)
--min-cells 100 --min-genes 5000   # QC 硬门:低于阈值整样本拒绝(exit 2,不产出 h5ad)
--no-gate                          # 跳过上述两道门(如珍稀小样本);指标照算照记

```

## 5. 批量 12 个器官(真实运行记录)

```bash
for f in data/tabula_muris_droplet/*.h5ad; do
  bash run.sh standardize "$f" -o "data/out/$(basename "$f" .h5ad)"
  case $? in 3) echo "$f" >> pending.txt ;; 2) echo "$f" >> rejected.txt ;; esac
done
```

| organ | status | cells | genes_kept | genes_dropped | 耗时 |
|---|---|---|---|---|---|
| Bladder | ok | 2351 | 23143 | 290 | 39.5s |
| Heart_and_Aorta | ok | 624 | 23143 | 290 | 17.5s |
| Kidney | ok | 2781 | 23143 | 290 | 21.6s |
| Limb_Muscle | ok | 4536 | 23143 | 290 | 18.7s |
| Liver | ok | 1845 | 23143 | 290 | 17.4s |
| Lung | ok | 4956 | 23143 | 290 | 59.1s |
| Mammary_Gland | ok | 4481 | 23143 | 290 | 20.4s |
| Marrow | ok | 3652 | 23143 | 290 | 20.7s |
| Spleen | ok | 9552 | 23143 | 290 | 21.0s |
| Thymus | ok | 1429 | 23143 | 290 | 18.1s |
| Tongue | ok | 7538 | 23143 | 290 | 24.7s |
| Trachea | ok | 11269 | 23143 | 290 | 22.9s |

物种 12/12 全判对(mouse,symbol_overlap),12/12 status=ok。
**丢的是基因不是细胞**(细胞 0 丢失);12 个器官共用同一套 2018 版基因注释
(var 完全相同),所以被丢的是同一批 290 个 feature:92 个 ERCC spike-in
(non_gene_feature)+ 183 个 unmapped(LOC 临时名、克隆名 4932411L15 之类、
假基因 *-ps、miRNA 旧名)+ 15 个 ambiguous(旧符号对应多个现行基因)。
耗时与细胞数基本无关(大头在基因数);个别样本偏慢(Lung 59s)是共享文件系统
抖动,timings 块能看出慢在哪一步。

## 6. 用输出

```python
import anndata as ad
A = ad.read_h5ad("data/out/Marrow/standardized.h5ad")

A.layers["counts"]                  # 整数 counts(下游 doublets/HVG 用这个)
A.X                                 # log1p(normalize_total(counts, 1e4)), float32
A.obs[["total_counts", "n_genes_by_counts",
       "pct_counts_mt", "pct_counts_hb"]]      # 本步权威计算的 QC 列
A.var[["original_feature_name", "mapping_status"]]   # 改名溯源
A.uns["eca_pp_standardize"]       # {step_version, species, counts_source, ...}
```

## 7. 别的机器 / 容器

代码零集群依赖:任意机器

```bash
pip install git+https://github.com/chansigit/stancounts \
            git+https://github.com/chansigit/stangene
pip install ".[probe,agent]"        # 本 repo 的 checkout
```

Sherlock 容器环境用 `bash scripts/test-in-container.sh` 一键验证
(python:3.12-slim + Apptainer)。

## 8. 识别批次列 / 细胞类型列(identify-columns)

### 8.1 配置 agent harness(不配也能确定性降级)

默认是 DSH + Doubao。Sherlock 上使用源码构建的 dsh CLI:

```bash
export ARK_API_KEY=...
export DSH_BIN=$SCRATCH/tools/deepseek-harness-src/apps/cli/lib/bin.js
# DSH_BIN 不设时会自动尝试上面这个路径
```

要切回 Claude Agent SDK:

```bash
export HARNESS=claude
export ANTHROPIC_API_KEY=sk-ant-...  # 或先用 Claude Code CLI 登录
export ECA_PP_CLAUDE_CLI=/path/to/claude  # 可选
pip install ".[claude]"
```

两种后端共用 submit-tool 协议:模型必须调用提交工具,Python handler 校验
结构与候选合法性;无效提交会在同一 session 内返回错误并允许修正。

**选择模型**:`--model MODEL_ID`(或环境变量 `ECA_PP_AGENT_MODEL`)。
不指定时随后端使用 `doubao-seed-2-1-turbo-260628` 或
`claude-sonnet-5`。每轮实际模型记录在 `decisions[].usage.model`。

**不配置会怎样**:不报错——自动使用确定性 policy 跑完;无法安全判断时
输出 null 和结构化 warning(exit 0)。

### 8.2 跑一个真样本(Marrow,3652 细胞)

```bash
bash run.sh identify-columns data/out/Marrow/standardized.h5ad \
    -o data/out/Marrow_idc
```

真实结果(约 90 秒,agent 一轮试验即判定):

```
columns.batch     = channel   correction=unnecessary
  ← agent 从列间关系图看出 channel ≡ mouse.id,探一列即覆盖两个候选;
    整合前 iLISI 0.92(批次本已混合)、校正增益仅 0.013 → 无需校正
columns.cell_type = cell_ontology_class
```

查证据:`result.json` 的 `decisions`(每轮决策 + agent 回复全文)、
`trials`(指标)、`trial_1_umap.png`(图)。下游消费:
`--batch-col channel`;`correction=unnecessary` 时 integration 应跳过校正。

**费用可见**:每轮的 token 用量与费用记录在 `decisions[].usage`,汇总在
`metrics.llm`。DSH/Ark 与 Claude 的账户级消费分别在各自控制台查询;
`result.json` 的 `billing_url` 会随当前后端给出入口。
