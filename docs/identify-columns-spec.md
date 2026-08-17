# identify-columns 环节 · 需求规格

> eca-pp · 2026-08-17 · v0.3 · 状态:待终审(代码未动工)
> 前置:standardize 已交付;本环节从其范围中移出(原 F6),从零实现,不依赖 stanmetacols。

## 1. 目标与定位

从 standardized.h5ad 的 obs 中识别**批次列**(供 doublets 分组与 integration
校正)与**细胞类型列**,并以试验证据支撑结论:批次列选错往往到整合阶段才暴露
(批次过多/病态时 Harmony 不收敛或产出劣质嵌入),因此识别内置**小规模整合
试验**做经验验证,而不只看列名和值分布。

本环节服务于**图谱构建**:整合的目的是把不同实验中的同类细胞身份对应起来。
由此推出与差异分析场景不同的教义(§4)。

**这是本项目第一个 agent 型环节。** 对调用方它仍是普通 CLI(h5ad 进 →
result.json 出 → 退出码汇报);内部由 Agent SDK 驱动"提出候选 → 试验 →
评估 → 换候选"的循环。红线:

> **确定性环节内不放 agent loop。agent 型环节须显式声明(本环节是第一个),
> 其所有工具调用必须是确定性 CLI——每轮试验可复现、可比较、可审计;
> agent 只承担候选选择与终止判断,每次选择连同理由记入 trials。**

## 2. 判定空间

批次结论四分,均为一等公民:

| 判定 | 含义 | 状态 |
|---|---|---|
| ① 批次列 + **建议校正** | 找到列,试验证实校正有益 | ok |
| ② 批次列 + **无需校正** | 找到列,但整合前 iLISI 已高 / 批次效应可忽略 | ok(`correction: "unnecessary"`) |
| ③ **无批次结构** | 单批次数据:候选列均不存在/常量/结构性不合格,证据充分 | ok(`batch: null`) |
| ④ 定不了 | 证据不足或候选试尽无合格者 | exit 3,附全部证据 |

细胞类型结论:找到列 / 未找到(null)/ 存疑(needs_review)。

## 3. obs 画像(确定性,ecasteps 内实现)

三层证据,整体写入 result.json,agent 与人工复核看同一份:

1. **逐列统计 + 抽样值**:列名、dtype、基数、缺失率、是否常量、是否每细胞
   唯一、去重值样本(≤20 个,**附各值细胞数**)。值是主信号,列名会骗人。
2. **分组体检**(候选分组列):组数、每组细胞数 min/median/max、归一化熵
   (均衡度)、微组(<25 细胞)占比。
3. **列间关系图**(候选分组列两两):**嵌套**(A 的每组恰落入 B 的一组 →
   A 细分 B)与**等价**(同划分不同标签)。groupby-nunique 即得,构成层级
   搜索的格。

另做派生候选枚举:barcode 前/后缀拆分(obs_names 分隔符结构)、双列复合
(限基数乘积合理者)。

## 4. 教义(进 agent prompt,全文见附录 A)

1. **先分类再选择**。分组列归类:技术分组(lane/channel/library/run/pool/
   hash)、供体分组(donor/mouse/patient)、实验条件(disease/treatment/
   timepoint/genotype)、注释列、QC 数值列、标识符。
2. **候选池**:批次候选 = 任何"希望跨其对齐细胞身份"的分组——技术分组、
   供体、**实验条件都算**(图谱场景:跨条件对齐同类细胞正是目的;条件间
   生物学差异不丢——counts 原样在产物里,差异分析回 counts 层)。
   采纳条件性分组时 evidence 须标注:"校正将合并条件驱动的表达偏移
   (图谱对齐之目的);差异分析请回 counts 层"。
3. **结构性不合格**(永不做批次):注释列、QC 数值列、每细胞唯一标识符、
   常量列。
4. **嵌套时自底向上**:批次 = 最细的可行技术层。典型设计里 library 嵌套于
   condition,**在最细技术层校正已顺带实现跨条件对齐**;细层病态(§5 预检)
   → 上移一层。条件列亲自做批次,主要发生在"数据未记录任何技术分组"时。
5. **正交分组**:v0.3 只选一列——选试验表现好者,另一者的存在如实记录。
   多协变量校正后置。
6. **宁 exit 3 不硬猜**。

## 5. integration-probe · CLI 契约(确定性工具)

```bash
ecasteps-integration-probe SRC.h5ad -o OUTDIR --batch-col SPEC \
    [--cell-type-col SPEC] [--n-cells 5000] [--n-hvg 2000] [--seed 0]
```

- `SPEC` = obs 列名,或派生列值文件路径(TSV:`cell_id<TAB>value`)。
- **抽样**:固定 seed 均匀抽样 ≤ n-cells。同输入同 seed → 同细胞子集,
  **一次识别运行内所有试验共用同一子集**——候选间指标差异只来自批次列本身。
  不按被试列分层(微批次在均匀抽样下指标难看是真实信号,不掩盖)。
- 流程:抽样 → HVG → PCA → Harmony → 指标 → UMAP 图。
- **指标**(经 harmonypy `compute_lisi`,零新增依赖):
  - **iLISI**(批次混合度):整合前 / 后,原始值 + 归一值
    `(iLISI−1)/(n_batches−1)`。整合前 iLISI 已高 → 判定 ②(无需校正);
    整合后提升量是校正收益的核心信号。
  - **cLISI**(细胞类型纯度):整合前 / 后,归一 `(n_types−cLISI)/(n_types−1)`。
    守"同类聚齐、异类不混";无细胞类型候选列时,在**整合前的 kNN 图上跑
    Leiden**(固定 seed 与 resolution)取伪标签——图社区才贴合转录组的
    簇结构,`clisi_labels: "annotated" | "pseudo"` 如实标注。
  - `harmony_converged`:**不收敛/崩溃是合法观测结果**(status=ok、如实
    记录)——病态批次的最强信号。
  - 预检随附:n_batches、每组细胞数分布、微组占比、批次×细胞类型混杂度、
    批次对前若干主成分的回归方差(辅助 ② 判定)。
- **UMAP 图**:每轮固定 seed 出双面板 PNG(按批次着色 / 按细胞类型候选
  着色),供 agent 视觉复核与人工复盘。纪律:**指标定案,图只做旁证与
  否决**(以图否决须写明理由入 trials)。
- 退出码:0 = 试验完成(无论整合好坏);2 = 输入不合法;1 = 意外错误。
- 依赖走 extras `[probe]`:scanpy、harmonypy、scikit-learn、umap-learn、
  leidenalg(伪标签 Leiden 用)。

**为何固定 Harmony 而非 scVI**:probe 是诊断仪器,不是生产整合器——要灵敏
不要鲁棒。Harmony 秒级、CPU、seed 稳定,且对病态批次"崩得响亮"(此脆弱性
即诊断信号);scVI 分钟级、盼 GPU、训练随机,其鲁棒性反而掩盖病态。正式
integration 环节自行选方法(Harmony / scVI / MrVI);批次列的身份判定与
层级选择是方法无关的。已知偏差如实声明:对数百微批次但结构规整的数据
(如 plate-based),probe 会保守地推上一层,likelihood 类整合器或可承受
更细层——遇此情形 result.json 留注,供正式整合改判。

## 6. 病态批次:三重防线

1. **预检阈值**:微组(<25 细胞)占比 <5% 细胞 → 列可采纳,附注;
   **多数组为微组** → 列判结构病态,不进试验(多半是过细复合列/准标识符);
2. **试验实证**:漏网者在 probe 里现形(不收敛、iLISI 崩),如实记录;
3. **层级逃生**:自底向上搜索自然上移一层。

## 7. identify-columns · 环节契约

```bash
ecasteps-identify-columns SRC.h5ad -o OUTDIR \
    [--max-probes 6] [--n-cells 5000] [--no-probe] [--seed 0]
```

流程:

```
① obs 画像(确定性,§3)
② 候选生成:启发式枚举 ∪ agent 从画像补提(教义 §4 约束)
③ 确定性预检:结构性不合格与病态者出局(§4.3、§6.1)
④ 试验循环(≤ max-probes):agent 按自底向上顺序选候选 → probe → 读指标
   与 UMAP → 采纳 / 判"无需校正" / 换下一个
⑤ 细胞类型列:值词汇判断(agent)+ cLISI 佐证(与 ④ 共用试验)
⑥ 定案写 result.json;选中派生列时另写值文件(§8)
```

- **终止条件**:出现"收敛 + iLISI 提升显著 + cLISI 不劣化"的候选 → ①定案;
  整合前 iLISI 已高 → ②定案;候选池空且证据充分 → ③定案;
  试尽或达 max-probes 无合格者 → ④ exit 3。
- `--no-probe`:仅画像 + 排名(降级模式,快而弱)。
- **无 API key / Agent SDK 不可用**:确定性降级——画像 + 启发式排名照常
  产出,status=needs_review,exit 3。
- 退出码沿用全项目契约:0 / 3 / 2 / 1。

## 8. 产物

| 文件 | 内容 |
|---|---|
| `OUTDIR/result.json` | 画像、候选、每轮试验(指标 + 取舍理由 + UMAP 路径)、判定(§2) |
| `OUTDIR/trial_<n>_umap.png` | 每轮试验的双面板 UMAP |
| `OUTDIR/batch.tsv` | 仅当选中**派生**批次列:`cell_id<TAB>value`;真实 obs 列直接给列名 |

下游工具(doublets、integration、probe 自身)的 `--batch-col` 统一接受
"obs 列名或 TSV 路径";输入 h5ad 永不被修改。

```jsonc
{
  "schema_version": 1,
  "step": "identify_columns",
  "status": "ok | needs_review | rejected | error",
  "exit_code": 0,
  "src": "…/standardized.h5ad",
  "params": { "max_probes": 6, "n_cells": 5000, "no_probe": false, "seed": 0 },
  "profile": {
    "columns": [ { "column": "", "dtype": "", "n_unique": 0, "entropy": 0.0,
                   "missing_frac": 0.0, "examples": {"值": 123},
                   "group_sizes": { "min": 0, "median": 0, "max": 0,
                                    "tiny_frac": 0.0 } } ],
    "relations": [ { "finer": "", "coarser": "", "kind": "nested | equivalent" } ],
    "derived": [ { "label": "", "kind": "barcode | composite", "n_groups": 0 } ]
  },
  "candidates": { "batch": [], "cell_type": [] },
  "trials": [ { "batch_col": "",
                "metrics": { "ilisi_pre": 0.0, "ilisi_post": 0.0,
                             "ilisi_norm_pre": 0.0, "ilisi_norm_post": 0.0,
                             "clisi_norm_pre": 0.0, "clisi_norm_post": 0.0,
                             "clisi_labels": "annotated | pseudo",
                             "harmony_converged": true,
                             "n_batches": 0, "pc_regression_r2": 0.0 },
                "umap": "trial_1_umap.png",
                "verdict": "adopted | rejected | correction_unnecessary",
                "reason": "" } ],
  "columns": {
    "batch": { "value": "obs 列名 或 batch.tsv 路径 或 null",
               "kind": "existing | derived | null",
               "correction": "recommended | unnecessary | null",
               "confidence": 0.0, "evidence": "" },
    "cell_type": { "value": "", "kind": "existing", "confidence": 0.0,
                   "evidence": "" }
  },
  "metrics": { "timings": {} }
}
```

## 9. 非功能要求

- **每轮留痕**:agent 的每次选择与理由入 trials;同输入同 seed 下 probe 指标
  逐位可复现(agent 的选择顺序可异,每条证据可独立复算)。
- **费用有界**:agent 调用次数受 max-probes 与流程结构约束,无开放循环。
- 打包:extras `[probe]`(scanpy/harmonypy/sklearn/umap-learn)、
  `[agent]`(claude-agent-sdk);核心包不引入重依赖。
- 双环境验收(原生 + 容器)照旧;agent 相关测试以 mock 驱动,不真调 API。

## 10. 验收测试(草案)

| 用例 | 期望 |
|---|---|
| 画像:常量 / 每细胞唯一 / 均衡分组列 | 字段正确,前两者预检出局 |
| 嵌套列对(donor ⊃ lane 合成) | relations 正确;搜索先试 lane |
| barcode 前缀批次(合成) | 派生候选被枚举,选中时物化为 TSV |
| probe:两批次合成数据,真批次列 vs 打乱列 | 真列 iLISI 提升显著高于打乱列 |
| probe:批次效应注入 vs 无效应数据 | 后者整合前 iLISI 已高 → correction_unnecessary |
| probe:数百微批次(病态) | 不崩,harmony_converged=false 如实记录 |
| 单批次数据(无分组列) | 判定 ③,batch=null,exit 0 |
| identify-columns 降级(无 key / --no-probe) | exit 3,画像 + 排名完整 |
| agent 流程(mock SDK) | trials 留痕完整,终止条件逻辑正确 |
| Tabula Muris 实跑 | 识别出 channel/mouse.id 类批次列,cell_type 列命中 |

## 11. 不做的事(non-goals)

- organ / tissue 等其余角色——后置;
- 写规范 obs 列、修改任何输入 h5ad——只产出结论与证据;
- 多协变量校正(harmony 多 key)——后置;
- 整合方法选型(scVI/MrVI probe)——probe 固定 Harmony(§5 已论证),
  正式整合是另一环节的职责;
- 细胞类型重注释——只识别"哪列是作者的标注",不生产新标注。

## 附录 A · agent prompt(草案,随实现定稿)

```
You identify two roles among the obs columns of a standardized scRNA-seq
dataset: the BATCH column (for integration) and the CELL TYPE column. You
work for an atlas-building pipeline: the purpose of integration is to align
cell identities across experiments.

Evidence provided: a three-layer profile (per-column stats with sampled
values and per-value cell counts; group-size health for grouping columns;
a nesting/equivalence graph among grouping columns), plus the metrics and
UMAP image of every probe trial run so far.

Doctrine:
1. Classify grouping columns first: technical (lane/channel/library/run/
   pool/hash), donor, experimental condition (disease/treatment/timepoint/
   genotype), annotation, QC numeric, identifier.
2. Batch candidates: ANY grouping across which cell identities should be
   aligned - technical, donor, AND experimental condition (atlas setting:
   merging condition-driven expression shifts is the goal; count data stay
   untouched for downstream differential analysis). When you adopt a
   condition-like column, state this consequence in your evidence.
3. Never batch: annotation columns, QC numeric columns, per-cell-unique
   identifiers, constant columns.
4. Nested groupings: prefer the finest viable technical level (correcting
   at library level already aligns across conditions when libraries nest
   within conditions). If a level is pathological (mostly tiny groups, or
   a probe shows non-convergence / poor iLISI), move one level up.
5. Orthogonal groupings: choose ONE column - the one that probes better;
   record the other's existence.
6. Decide from probe metrics first (iLISI gain, cLISI preservation,
   convergence); use the UMAP image only as supporting evidence or a veto,
   and state the reason when vetoing. If pre-integration iLISI is already
   high, conclude "correction unnecessary". If nothing qualifies, stop and
   report undecidable rather than guessing.
Every action you take must name the candidate, the reason, and the
expected signal, in structured form.
```
