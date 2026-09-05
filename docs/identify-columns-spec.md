# identify-columns 环节 · 需求规格

> eca-pp · 2026-09-04 · v0.5 · 状态:✅ 回归验证通过(122 项)
> 前置:standardize 已交付;本环节从其范围中移出(原 F6),独立实现,不依赖 stanmetacols。

## 1. 目标与定位

从 standardized.h5ad 的 obs 中识别**批次列**(供 doublets 分组与 integration
校正)与**细胞类型列**,并以试验证据支撑结论:批次列选择不当往往到整合阶段
才暴露(批次过多或病态时 Harmony 不收敛或产出劣质嵌入),因此识别过程内置
**小规模整合试验**进行经验验证,而非仅依据列名与取值分布。

本环节服务于**图谱构建**:整合的目的是对齐不同实验中同类细胞的身份。由此
推出与差异分析场景不同的判定准则(§4)。

**这是本项目第一个 agent 型环节。** 对调用方它是普通 CLI(h5ad 输入 →
result.json 输出 → 退出码汇报);内部只有**一次**模型调用:模型通读每列的
取值计数表,一次给出批次列排序(≤3)与细胞类型列;随后程序按序用整合试验
验证批次候选。v0.5 起不再有"多轮决策循环"——v0.4 的名字分类器 + 梯队门禁 +
逐轮 agent 选择被证明是过度设计:模型只能在规则预筛后的列表里选,看到了
`ann0608` 的取值是细胞类型名也提交不了(abm-ilcp 测试,2026-09-04)。架构约束:

> **确定性环节内不引入 agent loop。agent 型环节须显式声明(本环节是第一个),
> 其所有工具调用必须是确定性 CLI——试验可复现、可比较、可审计;
> 模型只负责"读证据、做分类",试验与判定阈值由程序执行,分类连同理由记入审计记录。**

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
   唯一、取值计数表(基数 ≤50 时全量,否则最高频 50 个,附各取值的细胞数)。取值内容是主要信号;
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


模型在一次调用里读完整个画像后回答两个问题;程序只做存在性与可探测性校验,
不再按列名硬性过滤模型的答案。

1. **看值不看名**。程序给出的 `heuristic_class` 只是名字启发式的提示;列的
   取值计数表才是依据。
2. **批次列排序(≤3)**:技术因素(lane/channel/library/run/pool/well)优先,
   其次供体/样本/动物,再次实验条件;嵌套的技术层级取"组不以微组为主"的最细
   一层。嵌套在细胞类型样文本列之内、或取值形如 `<batch>-<celltype>` 的列是
   批次×注释的复合列,应改用较粗的技术列。注释列、QC 数值、每细胞标识符、
   常量、cluster ID 不能做批次;只允许排入程序标为可探测的列。空列表 = obs
   中不存在合理的批次结构。
3. **细胞类型列**:作者注释,由取值判断(谱系名、本体术语、proB/CDP/ILC2P 一类
   缩写),列名无关紧要(ann0608、ImmGen_refine、labels_v2 皆可);绝不是算法
   聚类(leiden/louvain/seurat_clusters/纯整数)。多列并存时取命名可辨且粒度
   可用者,其余在理由中提及;不存在则 null。
4. **程序验证**:按排序依次 probe(≤ `--max-probes`,默认 2),第一个满足
   "收敛 + iLISI 提升 ≥0.05 + cLISI 不劣化"或"整合前已混合"的候选即为结论;
   全部不合格 → `batch: null` + 结构化 warning,不猜测。
5. 采纳实验条件或语义不明的分组做批次时记 `biological_batch_fallback`。

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
  - 画像记录组规模与微组比例;probe 记录批次数、组大小摘要与批次对
    主成分的回归方差(辅助判定 ②)。当前未输出批次×细胞类型混杂度指标。
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
3. **排序回退**:模型给出的排序中前一候选不合格时试验下一个(≤ max-probes)。

## 7. identify-columns · 环节契约


```bash
eca-pp-identify-columns SRC.h5ad -o OUTDIR \
    [--max-probes 2] [--n-cells 5000] [--no-probe] [--seed 0] [--model ID]
```

流程:

```
① obs 画像(确定性,§3)→ 候选与可探测集(名字启发式 + 病态预检,§6.1)
② 证据表:逐列 取值计数表 / dtype / 基数 / 缺失率 / heuristic_class /
   组规模 / 嵌套父列 / 等价列,加派生候选与关系图
③ 一次分类调用(§4):batch_ranked(≤3)+ cell_type + 逐列类别 + 理由;
   校验:批次列须在可探测集内且不重复,cell_type 须是 obs 列、非每细胞唯一、
   非 cluster;不合法的提交连同原因退回模型在同一会话内改正
④ 细胞类型列:直接采用分类结果;名字启发式未识别为 annotation 的列被选中时
   记 `cell_type_identified_from_values`;常量注释仍是有效输出但不用于 cLISI
⑤ 按 batch_ranked 顺序 probe(≤ max-probes),首个合格者即结论;
   cell_type 作为 cLISI 标签列
⑥ 判定写入 result.json;选中派生列时另写值文件(§8)
```

- **终止条件**:某候选"收敛 + iLISI 提升显著 + cLISI 不劣化" → 判定 ①;
  整合前 iLISI 已高 → 判定 ②;batch_ranked 为空 → 判定 ③(`no_batch_candidate`);
  排序内候选均不合格或超出 max-probes → 判定 ④(`batch_evidence_insufficient`)。
- 细胞数 < probe 下限(300)时不做试验,记 `dataset_too_small_to_probe`,
  cell_type 照常产出。
- `--no-probe`:仅输出画像、分类与细胞类型;batch=null、exit 0,
  记录 `probe_disabled` warning。
- `--model ID`(或环境变量 `ECA_PP_AGENT_MODEL`):指定模型;缺省随 `HARNESS`。
  实际模型记录在 `classification.usage.model`,汇总于 `metrics.llm.models`。
- `AGENT_WALL_MIN`:单次 agent run 墙钟上限,默认 6 分钟(medium reasoning 需要)。
  超时不重试同一请求,改用名字启发式(`HeuristicClassifier`)继续。
- **无 API 凭据 / harness 不可用或调用失败**:名字启发式给出排序与 cell_type,
  记 `agent_unavailable` / `agent_failed` warning。
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
  "params": { "max_probes": 2, "n_cells": 5000, "no_probe": false, "seed": 0 },
  "profile": {
    "columns": [ { "column": "", "dtype": "", "n_unique": 0, "entropy": 0.0,
                   "missing_frac": 0.0, "examples": { "<取值>": 0 },
                   "group_sizes": { "min": 0, "median": 0, "max": 0,
                                    "tiny_frac": 0.0 } } ],
    "relations": [ { "finer": "", "coarser": "", "kind": "nested | equivalent" } ],
    "derived": [ { "label": "", "kind": "barcode | composite", "n_groups": 0 } ]
  },
  "candidates": {   // 名字启发式 + 病态预检:可探测集与兜底排序
    "batch": [ { "label": "", "class": "technical | donor | condition | other | derived",
                   "excluded": false, "note": "",
                   "nested_within": [ { "column": "", "class": "" } ],
                   "equivalent_to": "先出现的 canonical 候选或缺省" } ],
    "cell_type": [ { "label": "", "class": "annotation | other | cluster",
                     "usable_for_clisi": true, "note": "" } ]
  },
  "warnings": [ { "code": "cell_type_not_found",
                   "message": "...", "details": {} } ],
  // 一次分类调用的完整记录(raw_reply = 模型回复全文)
  "classification": { "source": "agent | deterministic | deterministic_fallback",
                      "batch_ranked": [ { "column": "", "class": "", "reason": "" } ],
                      "cell_type": null, "cell_type_reason": "",
                      "columns": { "<列名>": "<class>" }, "notes": "",
                      "tools_used": [], "raw_reply": "",
                      "usage": { "model": "", "cost_usd": 0.0, "input_tokens": 0,
                                 "output_tokens": 0, "reasoning_tokens": 0 } },
  // 兼容旧汇总脚本的审计形态:恒为一条 action=classify
  "decisions": [ { "action": "classify", "candidate": "", "cell_type": null,
                   "source": "", "reason": "", "tools_used": [], "raw_reply": "",
                   "usage": {} } ],
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
                "class": "模型给该候选的类别", "reason": "模型的排序理由",
                "probe_reasons": [] } ],
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

- **审计**:分类结果与理由记入 classification/trials;同输入同 seed
  下 probe 指标完全可复现(模型的排序可能不同,但每条证据可独立复算验证)。
- **调用开销有界且可见**:每次运行恰好一次模型调用(失败则零次 + 兜底);
  token 用量与费用记入 `classification.usage`(OpenAI Responses 另记
  reasoning token)。
  `metrics.llm` 中 `calls` 包含失败/超时尝试,
  并分列 `successful_calls` / `failed_calls` / `timeout_calls` /
  `failed_seconds` / `failures`(附账户级账单查询 URL)。
- 打包:extras `[probe]`(scanpy/harmonypy/scikit-learn/leidenalg)、
  `[agent]`(OpenAI + DSH/MCP)、`[openai]`(仅 OpenAI)、`[llm]`(DSH 兼容 extra)、
  `[claude]`(Claude 后端);核心包不引入模型 SDK。
- 默认 `HARNESS=openai`、`doubao-seed-2-1-turbo-260628`、reasoning=`medium`(`OPENAI_AGENTS_REASONING_EFFORT` 可调;minimal 曾漏看值计数表中的注释列)。
  `HARNESS=deepseek|claude` 显式选择其他后端;`--model` 或
  `ECA_PP_AGENT_MODEL` 独立选择模型。模型故障时使用确定性策略,不自动换后端。
- OpenAI 直接注册 Python submit tool,不经过 MCP;默认服务端会话续接、
  `parallel_tool_calls=False`,未提交时最多同会话 nudge 两次。
  三个后端在本步骤都只开放校验后的 submit tool,完整证据表放入 prompt。
- probe exit 2 是候选拒绝;exit 1 或异常退出会使主流程 exit 1、status=error。
  Harmony 不收敛仍作为试验观测记录。cLISI 排除缺失作者标签,不足时用
  pseudo labels;只有实际使用 annotated labels 才在输出 evidence 声明使用作者标签。
- 等价候选必须分区与缺失位置都一致才可合并;整数值浮点 group ID 也可
  参与分组,常见整数 QC 指标仍被排除。
- `batch.tsv` 原子发布。重跑时旧 result.json、batch.tsv、candidates/、trial_N/
  移入 `.history/identify_columns-*/`。需保留的 DSH session 日志使用唯一文件名,
  避免覆盖此前轮次;运行记录与正式结果文件的原子发布约定不同。
- 验收沿用双环境流程(原生 + 容器);agent 相关测试以 mock 驱动,不发起
  真实 API 调用。

## 10. 验收测试

| 用例 | 期望 |
|---|---|
| 画像:常量列 / 每细胞唯一列 / 均衡分组列 | 字段正确;前两类被预检排除 |
| 嵌套列对(donor ⊃ lane,合成数据) | relations 正确;候选带 nested_within |
| 排序验证(合成数据,首选不合格) | 依序 probe,第二候选被采纳;`--max-probes 0` 不试验 |
| 未按名识别的文本列被模型选为 cell_type | 采纳并用于 cLISI,记 `cell_type_identified_from_values` |
| barcode 前缀批次(合成数据) | 派生候选被枚举,选中时物化为 TSV |
| probe:两批次合成数据,真实批次列 vs 打乱标签 | 真实列的 iLISI 提升显著高于打乱标签 |
| probe:注入批次效应 vs 无效应数据 | 后者整合前 iLISI 已高 → correction_unnecessary |
| probe:数百微批次(病态) | 程序正常完成,harmony_converged=false 如实记录 |
| 单批次数据(无分组列) | 判定 ③,batch=null,exit 0 |
| Agent 不可用 | 确定性 policy 接管并继续,exit 0,warning 完整 |
| `--no-probe` | batch=null,cell type 尽力识别,exit 0,`probe_disabled` warning |
| 仅有 cluster、没有作者 annotation | cell_type=null,exit 0,warning 说明原因 |
| 常量作者 annotation | 正常输出该列,但 probe 改用 pseudo-label |
| 分类校验(mock SDK) | 排除列 / 不可探测列 / cluster 列 / 标识符被退回;合法答案通过 |
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
You are given the obs (cell metadata) profile of a standardized scRNA-seq
dataset: for every column its name, dtype, number of distinct values,
missing fraction, and a value -> cell-count table (all values when there are
at most 50, else the 50 most frequent), plus nesting/equivalence relations
between grouping columns and derived candidates (barcode prefix/suffix,
two-column composites). Read the VALUES of every column — names are hints,
values are the truth — and answer two questions in one submission.

1. BATCH column(s), ranked, at most 3. The program will run a small Harmony
   integration trial on each in order and keep the first one that qualifies
   (clear iLISI gain with cell-type structure preserved, or "already mixed").
   - Prefer technical factors (lane/channel/library/run/pool/10x well),
     then donor/sample/animal, then experimental condition. Among nested
     technical levels prefer the finest one whose groups are not mostly tiny.
   - A column nested inside a cell-type-like column, or whose values look like
     "<batch>-<cell type>" (e.g. "ABM2-ILC2P.4"), is batch x annotation:
     use the coarser technical column instead.
   - Never a batch: annotation columns, QC numbers, per-cell identifiers,
     constants, cluster IDs. Only columns listed as probeable are allowed.
   - An empty list means no plausible batch structure exists in obs.
2. CELL TYPE column: the AUTHOR'S cell-type annotation, judged from values
   (lineage names, ontology terms, abbreviations such as proB, CDP, ILC2P,
   "1:CDP-like"), whatever the column is called (ann0608, ImmGen_refine,
   labels_v2 ...). Never an algorithmic clustering (leiden/louvain/
   seurat_clusters/bare integers). When several exist prefer the one with
   recognizable names at usable granularity and mention the others. null if
   none exists.

Also classify each grouping column (technical/donor/condition/annotation/
cluster/qc_numeric/identifier/constant/other) in "columns".

Submit exactly this JSON through the provided tool:
{"batch_ranked": [{"column": "<name>", "class": "<class>", "reason": "<why, citing values>"}],
 "cell_type": "<column or null>", "cell_type_reason": "<quote 2-3 values>",
 "columns": {"<column>": "<class>"}, "notes": "<anything else worth recording>"}
```
