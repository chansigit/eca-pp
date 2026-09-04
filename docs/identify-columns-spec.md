# identify-columns 环节 · 需求规格

> eca-pp · 2026-09-03 · v0.4 · 状态:✅ clean-container 回归验证通过
> 前置:standardize 已交付;本环节从其范围中移出(原 F6),独立实现,不依赖 stanmetacols。

## 1. 目标与定位

从 standardized.h5ad 的 obs 中识别**批次列**(供 doublets 分组与 integration
校正)与**细胞类型列**,并以试验证据支撑结论:批次列选择不当往往到整合阶段
才暴露(批次过多或病态时 Harmony 不收敛或产出劣质嵌入),因此识别过程内置
**小规模整合试验**进行经验验证,而非仅依据列名与取值分布。

本环节服务于**图谱构建**:整合的目的是对齐不同实验中同类细胞的身份。由此
推出与差异分析场景不同的判定准则(§4)。

**这是本项目第一个 agent 型环节。** 对调用方它是普通 CLI(h5ad 输入 →
result.json 输出 → 退出码汇报);内部由可切换 harness 驱动"提出候选 → 试验 →
评估 → 更换候选"的循环。架构约束:

> **确定性环节内不引入 agent loop。agent 型环节须显式声明(本环节是第一个),
> 其所有工具调用必须是确定性 CLI——每轮试验可复现、可比较、可审计;
> agent 仅承担候选选择与终止判断,每次选择连同理由记入审计记录。**

## 2. 判定空间

批次判定分为四类,均为正式判定结果:

| 判定 | 含义 | 状态 |
|---|---|---|
| ① 批次列 + **建议校正** | 已识别批次列,试验证实校正有益 | ok |
| ② 批次列 + **无需校正** | 已识别批次列,但整合前 iLISI 已高,批次效应可忽略 | ok(`correction: "unnecessary"`) |
| ③ **无批次结构** | 单批次数据:候选列均不存在、为常量或结构性不合格,且证据充分 | ok(`batch: null`) |
| ④ 无法判定 | 证据不足,或候选耗尽且无合格者 | ok(`batch: null`),附结构化 warning 与全部证据 |

细胞类型判定:已识别 / 未识别(null + 结构化 warning);两者均不阻断流程。

## 3. obs 画像(确定性,eca-pp 内实现)

三层证据,整体写入 result.json;agent 与人工复核使用同一份证据:

1. **逐列统计与抽样值**:列名、dtype、基数、缺失率、是否常量、是否每细胞
   唯一、去重取值样本(≤20 个,附各取值的细胞数)。取值内容是主要信号;
   列名可能与内容不符,不能单独作为依据。空串、纯空白及明确的 NA 拼写
   在画像前统一归一为缺失值,不计为真实分组。
2. **分组规模统计**(候选分组列):组数、每组细胞数 min/median/max、
   归一化熵(均衡度)、微组(<25 细胞)占比。
3. **列间关系图**(候选分组列两两比较):**嵌套**(A 的每组恰落入 B 的
   一组,即 A 细分 B)与**等价**(相同划分、不同标签)。由
   groupby-nunique 计算,构成层级搜索所依据的偏序结构。

另做派生候选枚举:barcode 前/后缀拆分(基于 obs_names 的分隔符结构)、
双列复合(限基数乘积在阈值内的组合)。

## 4. 判定准则(写入 agent prompt,全文见附录 A)

1. **先归类,后选择**。分组列归类为:技术分组(lane/channel/library/run/
   pool/hash)、供体分组(donor/mouse/patient)、实验条件(disease/
   treatment/timepoint/genotype)、注释列、QC 数值列、标识符。
2. **候选分层**:PRIMARY = 技术分组与供体/样本分组,包括 barcode 派生的
   技术结构;FALLBACK = 实验条件及语义不明的分组。所有可行 primary 候选
   必须先试完,只有它们均不合格时才允许试验 fallback。采纳 fallback 时
   `warnings` 必须注明可能合并生物学差异的风险。
3. **结构性不合格**(不得作为批次):注释列、QC 数值列、每细胞唯一
   标识符、常量列。
4. **嵌套分组自底向上**:批次 = 最细的可行技术层。典型实验设计中 library
   嵌套于 condition,在最细技术层校正即同时实现跨条件对齐;细层病态
   (§5 预检)时上移一层。条件列直接作为批次,主要发生在数据未记录任何
   技术分组的情况下。
5. **正交分组**:v0.4 仅选择一列——取试验指标更优者,另一列的存在记录于
   证据中。多协变量校正为后续工作。
6. **无法判定时不猜测,以 `batch: null`、exit 0 和结构化 warning 完成。**

## 5. integration-probe · CLI 契约(确定性工具)

```bash
eca-pp-integration-probe SRC.h5ad -o OUTDIR --batch-col SPEC \
    [--cell-type-col SPEC] [--n-cells 5000] [--n-hvg 2000] [--seed 0]
```

- **抽样**:对每个候选做固定 seed 的分组保留抽样(总量 ≤ n-cells):先为
  每个 batch group 保留最多 5 个细胞,再从剩余细胞均匀补齐。同输入、候选
  与 seed 得到同一子集;避免稀有但有效的 barcode group 在子样本中消失。
  **规模自适应**:`n_cells = clamp(50 × 候选中最大批次数, 5000, 30000)`,
  保证被试批次在子样本中每批期望 ≥50 细胞,避免"全量数据健康、子样本
  呈现伪病态"的采样偏差;候选集确定后计算一次抽样规模。
  相应地,**病态预检(§6.1)仅依据全量画像的组规模**,子样本仅用于计算
  整合指标。
- 流程:抽样 → HVG → PCA → Harmony → 指标。
- **指标**(经 harmonypy `compute_lisi`,无新增依赖):
  - **iLISI**(批次混合度):整合前 / 后,原始值与归一值
    `(iLISI−1)/(n_batches−1)`。整合前 iLISI 已高 → 判定 ②(无需校正);
    整合后的提升量是校正收益的核心信号。
  - **cLISI**(细胞类型纯度):整合前 / 后,归一
    `(n_types−cLISI)/(n_types−1)`,衡量细胞类型结构的保持程度(同类
    聚集、异类分离)。无细胞类型候选列时,在**整合前的 Gaussian 加权
    exact-sklearn kNN 图上运行 Leiden**(固定 seed 与 resolution)
    生成伪标签——该构图路径避免每个 probe 子进程重复的 Numba 编译;
    `clisi_labels: "annotated" | "pseudo"` 与
    `pseudo_label_graph: "knn_gauss_sklearn"` 如实标注。
    作者注释存在缺失值时,cLISI 只在非缺失子集上计算;
    `cell_type_coverage_sampled` 与 `n_cells_clisi` 记录覆盖率和实际样本数。
    可用注释少于 50 个细胞或只剩一类时改用伪标签。
    作者注释使用 0.05 的下降容忍度;伪标签是较弱证据,使用 0.15,但明显
    破坏结构时仍否决候选。
  - `harmony_converged`:**不收敛或运行失败是有效观测结果**(status=ok,
    如实记录)——它是病态批次的最强信号。
  - 预检数据随附:n_batches、每组细胞数分布、微组占比、批次×细胞类型
    混杂度、批次对前若干主成分的回归方差(辅助判定 ②)。
- 退出码:0 = 试验完成(与整合效果无关);2 = 输入不合法;1 = 意外错误。
- 依赖经 extras `[probe]` 安装:scanpy、harmonypy、scikit-learn、
  leidenalg(伪标签 Leiden 使用)。

探测方法固定为 Harmony;正式 integration 环节自行选择方法(Harmony /
scVI / MrVI),批次列的身份判定与层级选择与整合方法无关。已知偏差:对
数百微批次但结构规整的数据(如 plate-based),probe 的结论偏保守(倾向
上移一层)——出现此情形时在 result.json 中记录说明,供正式整合阶段改判。

## 6. 病态批次:三层防护

1. **预检阈值**:微组(<25 细胞)细胞占比 <5% → 列可采纳,附注说明;
   **多数组为微组** → 判为结构病态,不进入试验(通常为过细的复合列或
   近似标识符的列);
2. **试验验证**:预检未能排除者将在 probe 中暴露(不收敛或 iLISI 异常),
   如实记录;
3. **层级回退**:自底向上搜索在细层失败后自然上移一层。

## 7. identify-columns · 环节契约

```bash
eca-pp-identify-columns SRC.h5ad -o OUTDIR \
    [--max-probes 6] [--n-cells 5000] [--no-probe] [--seed 0] [--model ID]
```

流程:

```
① obs 画像(确定性,§3)
② 候选生成:启发式枚举 ∪ agent 依画像补充(受判定准则 §4 约束)
③ 确定性预检:结构性不合格与病态候选排除(§4.3、§6.1)
④ 试验循环(≤ max-probes):agent 按自底向上顺序选择候选 → probe →
   评估指标 → 安全余量充足时由 metric fast path 采纳 /
   判定"无需校正";否则回到 agent 试验下一候选或复核
⑤ 细胞类型列:只把作者 annotation 作为正式输出;常量 annotation 仍是
   有效元数据,但不用于 cLISI。cluster 仅作为画像证据,不存在作者注释时
   输出 null,probe 自行生成 pseudo-label。
⑥ 判定写入 result.json;选中派生列时另写值文件(§8)
```

- **终止条件**:出现"收敛 + iLISI 提升显著 + cLISI 不劣化"的候选 →
  判定 ①;整合前 iLISI 已高 → 判定 ②;候选集为空且证据充分 → 判定 ③;
  候选耗尽或达 max-probes 上限而无合格者 → 判定 ④(exit 0 + warning)。
- `--no-probe`:仅输出画像、候选排名与细胞类型推断;batch=null、exit 0,
  并记录 `probe_disabled` warning。
- `--model ID`(或环境变量 `ECA_PP_AGENT_MODEL`):指定 agent 模型;
  缺省值随 `HARNESS` 后端选择。每轮实际使用的模型记录在
  `decisions[].usage.model`,汇总于 `metrics.llm.models`。
- `AGENT_WALL_MIN`:单次 agent run 墙钟上限,默认 2 分钟。超时后
  不重试同一请求,改用确定性策略继续,以保证无人值守流程有界。
- **无 API 凭据 / harness 不可用或中途失败**:自动切换确定性 policy
  继续试验,并记录 `agent_unavailable` / `agent_failed` warning。
- 退出码沿用全项目契约:0 / 3 / 2 / 1。

## 8. 产物

| 文件 | 内容 |
|---|---|
| `OUTDIR/result.json` | 画像、候选、每轮试验(指标、理由)、判定(§2) |
| `OUTDIR/batch.tsv` | 仅当选中**派生**批次列时产生(`cell_id<TAB>value`);选中真实 obs 列时 `columns.batch.value` 直接为列名,无此文件 |

下游工具(doublets、integration、probe 自身)的 `--batch-col` 统一接受
"obs 列名或 TSV 路径";输入 h5ad 不被修改。

```jsonc
{
  "schema_version": 2,
  "step": "identify_columns",
  "status": "ok | error",
  "exit_code": 0,
  "src": "…/standardized.h5ad",
  "params": { "max_probes": 6, "n_cells": 5000, "no_probe": false, "seed": 0 },
  "profile": {
    "columns": [ { "column": "", "dtype": "", "n_unique": 0, "entropy": 0.0,
                   "missing_frac": 0.0, "examples": { "<取值>": 0 },
                   "group_sizes": { "min": 0, "median": 0, "max": 0,
                                    "tiny_frac": 0.0 } } ],
    "relations": [ { "finer": "", "coarser": "", "kind": "nested | equivalent" } ],
    "derived": [ { "label": "", "kind": "barcode | composite", "n_groups": 0 } ]
  },
  "candidates": {
    "batch": [ { "label": "", "tier": "primary | fallback",
                   "equivalent_to": "先出现的 canonical 候选或缺省" } ],
    "cell_type": [ { "label": "", "output_eligible": true,
                     "usable_for_clisi": true } ]
  },
  "warnings": [ { "code": "cell_type_not_found",
                   "message": "...", "details": {} } ],
  // 决策审计:每轮 action/候选/理由 + agent 实际工具调用
  // + raw_reply:agent 该轮回复全文(含结构化决策块之外的分析叙述)
  "decisions": [ { "action": "", "candidate": "", "cell_type": null,
                   "source": "agent | metric_fast_path | deterministic | deterministic_fallback",
                   "reason": "",
                   "tools_used": [],
                   "raw_reply": "",
                   "usage": { "model": "", "cost_usd": 0.0, "input_tokens": 0,
                              "output_tokens": 0, "cache_creation_tokens": 0,
                              "cache_read_tokens": 0, "num_turns": 0 } } ],
  "trials": [ { "batch_col": "",
                "metrics": { "ilisi_pre": 0.0, "ilisi_post": 0.0,
                             "ilisi_norm_pre": 0.0, "ilisi_norm_post": 0.0,
                             "clisi_norm_pre": 0.0, "clisi_norm_post": 0.0,
                             "clisi_labels": "annotated | pseudo",
                             "pseudo_label_graph": "knn_gauss_sklearn | null",
                             "cell_type_coverage_sampled": 1.0,
                             "n_cells_clisi": 5000,
                             "harmony_converged": true,
                             "n_batches": 0, "pc_regression_r2": 0.0,
                             "timings": { "load": 0.0, "hvg_pca": 0.0,
                                          "labels": 0.0, "harmony": 0.0,
                                          "lisi": 0.0, "total": 0.0 } },
                "verdict": "adopted | rejected | correction_unnecessary",
                "reason": "" } ],
  "columns": {
    "batch": { "value": "obs 列名 或 batch.tsv 路径 或 null",
               "kind": "existing | derived | null",
               "correction": "recommended | unnecessary | null",
               "confidence": 0.0, "evidence": "" },
    "cell_type": { "value": "", "kind": "existing", "confidence": 0.0,
                   "evidence": "" } // 或 null;两者均为成功结果
  },
  "metrics": {
    "timings": {},
    // LLM 用量汇总:cost_complete=false 表示部分调用未上报单价(如订阅认证)
    "llm": { "calls": 0, "models": [], "input_tokens": 0, "output_tokens": 0,
             "cache_creation_tokens": 0, "cache_read_tokens": 0,
             "cost_usd": 0.0, "cost_complete": true,
             "billing_url": "随 HARNESS 后端选择" }
  }
}
```

## 9. 非功能要求

- **逐轮审计**:agent 的每次选择与理由记入 decisions/trials;同输入同 seed
  下 probe 指标完全可复现(agent 的选择顺序可能不同,但每条证据可独立
  复算验证)。
- **调用开销有界且可见**:agent 调用次数受 max-probes 与流程结构约束,
  不存在无上界的循环;每轮成功决策的 token 用量与费用记入
  `decisions[].usage`(OpenAI Responses 另记 reasoning token)。
  `metrics.llm` 中 `calls` 包含失败/超时尝试,
  并分列 `successful_calls` / `failed_calls` / `timeout_calls` /
  `failed_seconds` / `failures`(附账户级账单查询 URL)。
- 打包:extras `[probe]`(scanpy/harmonypy/scikit-learn/leidenalg)、
  `[agent]`(DSH + MCP)、`[openai]`(可选 OpenAI Agents SDK + Doubao
  对照后端)、`[claude]`(可选 Claude 后端);
  核心包不引入重依赖。
- OpenAI 对照后端直接注册 Python submit tool,不经过 MCP;
  默认不暴露本地文件工具,依赖每轮 prompt 中的完整 state。
- 验收沿用双环境流程(原生 + 容器);agent 相关测试以 mock 驱动,不发起
  真实 API 调用。

## 10. 验收测试

| 用例 | 期望 |
|---|---|
| 画像:常量列 / 每细胞唯一列 / 均衡分组列 | 字段正确;前两类被预检排除 |
| 嵌套列对(donor ⊃ lane,合成数据) | relations 正确;搜索先试验 lane |
| barcode 前缀批次(合成数据) | 派生候选被枚举,选中时物化为 TSV |
| probe:两批次合成数据,真实批次列 vs 打乱标签 | 真实列的 iLISI 提升显著高于打乱标签 |
| probe:注入批次效应 vs 无效应数据 | 后者整合前 iLISI 已高 → correction_unnecessary |
| probe:数百微批次(病态) | 程序正常完成,harmony_converged=false 如实记录 |
| 单批次数据(无分组列) | 判定 ③,batch=null,exit 0 |
| Agent 不可用 | 确定性 policy 接管并继续,exit 0,warning 完整 |
| `--no-probe` | batch=null,cell type 尽力识别,exit 0,`probe_disabled` warning |
| 仅有 cluster、没有作者 annotation | cell_type=null,exit 0,warning 说明原因 |
| 常量作者 annotation | 正常输出该列,但 probe 改用 pseudo-label |
| agent 流程(mock SDK) | 审计记录完整,终止条件逻辑正确 |
| Tabula Muris 实际数据 | 正确识别 channel/mouse.id 类批次列与细胞类型列 |

## 11. 范围之外(non-goals)

- organ / tissue 等其余角色——后续工作;
- 写入规范 obs 列、修改任何输入 h5ad——本环节仅产出结论与证据;
- 多协变量校正(Harmony 多 key)——后续工作;
- 整合方法选型(scVI/MrVI probe)——probe 固定使用 Harmony,
  正式整合是另一环节的职责;
- 细胞类型重注释——仅识别"哪一列是作者的标注",不生成新标注。

## 附录 A · agent prompt(与实现同步)

```
You identify two roles among the obs columns of a standardized scRNA-seq
dataset: the BATCH column (for integration) and the CELL TYPE column. You
work for an atlas-building pipeline: the purpose of integration is to align
cell identities across experiments.

Evidence provided: a three-layer profile (per-column stats with sampled
values and per-value cell counts; group-size health for grouping columns;
a nesting/equivalence graph among grouping columns), plus the metrics of every
probe trial run so far.

Doctrine:
1. Classify grouping columns first: technical (lane/channel/library/run/
   pool/hash), donor, experimental condition (disease/treatment/timepoint/
   genotype), annotation, QC numeric, identifier.
2. Batch candidates have two strict tiers. PRIMARY = technical and donor/
   sample factors, including technical structure derived from barcodes. Probe
   every viable primary candidate before considering FALLBACK = experimental
   condition or unknown grouping. Never jump to fallback while an untried
   primary candidate remains. A fallback is acceptable only when no primary
   candidate qualifies and its probe preserves biological structure; state
   this biological-risk consequence in the reason.
3. Never batch: annotation columns, QC numeric columns, per-cell-unique
   identifiers, constant columns.
4. Nested groupings: prefer the finest viable technical level (correcting
   at library level already aligns across conditions when libraries nest
   within conditions). If a level is pathological (mostly tiny groups, or
   a probe shows non-convergence / poor iLISI), move one level up.
5. Orthogonal groupings: choose ONE column - the one that probes better;
   record the other's existence.
6. Decide from probe metrics (iLISI gain, cLISI preservation and convergence).
   If pre-integration iLISI is already
   high, conclude "correction unnecessary". If nothing qualifies, stop and
   report undecidable rather than guessing.
7. Cell type column = the AUTHOR'S cell-type annotation (biological names
   such as "T cell", "hepatocyte", ontology terms), NOT an algorithmic
   clustering (leiden / louvain / seurat_clusters / numeric cluster IDs).
   candidates.cell_type is pre-ranked (class "annotation" before
   "cluster"; text labels before bare integers) and best_cell_type is
   the current default. Override it when the sampled values show the
   default is wrong. A constant author annotation is valid output but is not
   used for cLISI. If only cluster columns exist, set null; the probe creates
   its own pseudo-labels. Missing cell type is a successful unattended
   outcome. Set cell_type in EVERY reply, probe included.
8. Exhausted or ambiguous batch evidence is also a successful unattended
   outcome: use give_up so the pipeline reports batch=null with warnings.
Every action you take must name the candidate, the reason, and the
expected signal, in structured form.
```
