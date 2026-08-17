# identify-columns 环节 · 需求规格

> eca-pp · 2026-08-17 · v0.3 · 状态:设计中
> 前置:standardize 已交付;本环节从其范围中移出(原 F6),从零实现,不依赖 stanmetacols。

## 1. 目标与定位

从 standardized.h5ad 的 obs 中识别两个角色列——**批次列**(供 doublets 分组与
integration 校正)与**细胞类型列**——并以试验证据支撑结论:批次列选错往往要到
整合阶段才暴露(批次过多时 harmony 可能不收敛或产出劣质嵌入),因此识别过程
内置**小规模整合试验**做经验验证,而不是只看列名和值分布。

**这是本项目第一个 agent 型环节。** 对调用方它仍是普通 CLI(h5ad 进 →
result.json 出 → 退出码汇报);内部由 Agent SDK 驱动"提出候选 → 试验 →
评估 → 换候选"的循环。原红线相应改写:

> **确定性环节内不放 agent loop。agent 型环节须显式声明(本环节是第一个),
> 其所有工具调用必须是确定性 CLI——每轮试验可复现、可比较、可审计;
> agent 只承担候选选择与终止判断。**

## 2. 组件划分

```
ecasteps-identify-columns   agent 型环节:识别批次列 / 细胞类型列
    │ 内部调用(每轮均为确定性 CLI)
    ├─ obs 画像(库函数,确定性)
    └─ ecasteps-integration-probe   小规模整合试验(独立可用的确定性工具)
```

## 3. obs 画像(确定性,ecasteps 内实现)

对每个 obs 列计算:列名、dtype、基数(n_unique)、抽样值(≤20 个)、组间
均衡度(归一化熵)、缺失率、是否常量、是否每细胞唯一。另做两类**派生候选**
枚举:barcode 前/后缀拆分(检测 obs_names 的稳定分隔符结构)、双列复合
(限基数乘积合理的组合)。画像整体写入 result.json——它同时就是 agent 的
输入证据和人工复核的材料。

## 4. integration-probe · CLI 契约

```bash
ecasteps-integration-probe SRC.h5ad -o OUTDIR --batch-col SPEC \
    [--cell-type-col SPEC] [--n-cells 5000] [--n-hvg 2000] [--seed 0]
```

- `SPEC` = obs 列名,或派生列值文件路径(TSV:`cell_id<TAB>value`,见 §6)。
- 流程:分层抽样 ≤ n-cells → HVG → PCA → Harmony → 指标。全程固定 seed,
  同输入同参数 → 同指标。
- **指标**(result.json `metrics`):
  - 预检:`n_batches`、每批次细胞数(min/median)、批次×细胞类型混杂度;
  - `harmony_converged`(**不收敛/崩溃是合法观测结果**:status=ok、
    `harmony_converged=false`、异常信息入 reasons——这正是要探测的信号);
  - 邻域批次混合熵(整合前 / 后,提升量是核心信号);
  - 若给 cell_type:label silhouette(整合前 / 后,结构保持度)。
- 退出码:0 = 试验完成(无论整合好坏);2 = 输入不合法(列不存在、
  细胞过少);1 = 意外错误。
- 依赖走 extras `[probe]`:scanpy、harmonypy、scikit-learn。

## 5. identify-columns · 环节契约

```bash
ecasteps-identify-columns SRC.h5ad -o OUTDIR \
    [--max-probes 6] [--n-cells 5000] [--no-probe] [--seed 0]
```

流程:

```
① obs 画像(确定性)
② 候选生成:启发式枚举(列名线索 + 基数/均衡度合格)∪ agent 从画像补提
③ 确定性预检:淘汰明显不合格者(常量、每细胞唯一、基数过高、组内细胞过少)
④ 试验循环(≤ max-probes):agent 选下一个批次候选 → integration-probe
   → 读指标 → 采纳 / 换下一个
⑤ 细胞类型列:值词汇判断(agent)+ probe 的 silhouette 佐证(与 ④ 共用试验)
⑥ 定案写 result.json;选中派生列时另写值文件(§6)
```

- **终止条件**:出现"收敛 + 混合熵显著提升 + 细胞类型结构不劣化"的候选即
  定案;候选试尽或达 max-probes 仍无合格者 → exit 3,附全部轮次证据,
  由调用方决定(接受某个次优候选 / 判定无批次)。
- `--no-probe`:跳过试验,仅画像 + 启发式/agent 排名(降级模式,快而弱)。
- **无 API key / Agent SDK 不可用**:确定性降级——画像 + 启发式排名照常
  产出,status=needs_review,exit 3(识别结论需要调用方确认)。
- 退出码沿用全项目契约:0 = 定案;3 = 需要调用方决策;2 = 输入不合法;
  1 = 意外错误。

## 6. 产物

| 文件 | 内容 |
|---|---|
| `OUTDIR/result.json` | 画像、候选清单、每轮试验的指标与取舍、最终结论(见 schema) |
| `OUTDIR/batch.tsv` | 仅当选中**派生**批次列时:`cell_id<TAB>value`。真实 obs 列直接给列名,无此文件 |

下游工具(doublets、integration、probe 自身)的 `--batch-col` 统一接受
"obs 列名或 TSV 路径"——派生列不回写 h5ad,输入文件永不被修改。

```jsonc
{
  "schema_version": 1,
  "step": "identify_columns",
  "status": "ok | needs_review | rejected | error",
  "exit_code": 0,
  "src": "…/standardized.h5ad",
  "params": { "max_probes": 6, "n_cells": 5000, "no_probe": false, "seed": 0 },
  "profile": [ { "column": "", "dtype": "", "n_unique": 0, "entropy": 0.0,
                 "missing_frac": 0.0, "examples": [] } ],
  "candidates": { "batch": [], "cell_type": [] },
  "trials": [ { "batch_col": "", "metrics": {}, "verdict": "adopted | rejected",
                "reason": "" } ],
  "columns": {
    "batch":     { "value": "obs 列名 或 batch.tsv 路径", "kind": "existing | derived",
                   "confidence": 0.0, "evidence": "" },
    "cell_type": { "value": "", "kind": "existing", "confidence": 0.0,
                   "evidence": "" }        // 未识别到 → null
  },
  "metrics": { "timings": {} }
}
```

## 7. 非功能要求

- **每轮留痕**:agent 的每次候选选择与理由记入 trials;同输入同 seed 下,
  probe 指标逐位可复现(agent 的选择顺序可能不同,但每条证据可独立复算)。
- **费用有界**:LLM/agent 调用次数受 max-probes 与流程结构约束,无开放循环。
- 打包:`ecasteps` 新增 extras `[probe]`(scanpy/harmonypy/sklearn)与
  `[agent]`(claude-agent-sdk);核心包不引入重依赖。
- 双环境验收(原生 + 容器)照旧;agent 相关测试以 mock 驱动,不真调 API。

## 8. 验收测试(草案)

| 用例 | 期望 |
|---|---|
| 画像:常量列 / 每细胞唯一列 / 均衡分组列 | 字段正确,前两者被预检淘汰 |
| barcode 前缀批次(合成) | 派生候选被枚举并物化为 TSV |
| probe:两批次合成数据,真批次列 vs 打乱列 | 真列混合熵提升显著高于打乱列 |
| probe:数百微批次(harmony 病态场景) | 不崩,harmony_converged=false 如实记录 |
| identify-columns 降级(无 key,--no-probe) | exit 3,画像 + 排名完整 |
| agent 流程(mock SDK) | trials 留痕完整,采纳逻辑符合终止条件 |
| Tabula Muris 实跑 | 识别出 channel/mouse.id 类真实批次列,cell_type 列命中 |

## 9. 不做的事(non-goals)

- organ / tissue 等其余角色——后置,待批次/细胞类型闭环验证后再扩;
- 写规范 obs 列、修改任何输入 h5ad——识别只产出结论与证据;
- 整合方法选型(harmony vs MrVI 等)——probe 固定用 harmony 作为探测器,
  正式整合是另一个环节的职责;
- 细胞类型的重注释——只识别"哪一列是作者的细胞类型标注",不生产新标注。
