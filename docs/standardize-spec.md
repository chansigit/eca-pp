# standardize 步骤 · 需求规格(v2)

> eca-prefect-v3 · 2026-08-16
> 分期:**v0.1 已交付**(F1/F3/F2)→ **v0.2 本期**(F4a/F4/F5/F7)→ **v0.3 后置**(F6/F8)
> 原则:从 eca-prefect-v2 继承的是**需求**,不是代码;实现可以重写。

## 1. 目标与架构原则

把任意来源的单样本 `.h5ad` 变成下游可无条件信赖的**标准形**,以独立 CLI 交付:
无框架依赖、无常驻进程、默认不联网,可被任何驾驶员(脚本 / Snakemake / agent / 人)调用。

三条架构原则,贯穿所有功能块:

1. **两层分工**:步骤 = 确定性 CLI(h5ad in → files out + result.json);
   驾驶层消费 result.json 与退出码。Agent SDK / agent loop **只存在于驾驶层**。
2. **步骤内 LLM 上限 = 两处单次调用**(anthropic SDK,structured output),
   均显式开启(`--llm`)、默认关、有确定性回退:物种解析 T2(§5.2)、
   metacols 排名(§5.5,v0.3)。
3. **拍板权归驾驶层**:凡是"选哪个"的可逆决策(batch 列、物种存疑、counts 歧义),
   步骤只给证据与排名,不替驾驶员定;定不了就 exit 3 阻塞待决。

## 2. 输出契约(完整目标,验收以此为准)

产物:`OUTDIR/standardized.h5ad`(v0.2 起)+ `OUTDIR/result.json`(每次必写,含失败)。

| # | 不变量 | 落地期 | 下游依赖方 |
|---|---|---|---|
| I1 | `layers["counts"]` = 整数原始 counts | v0.2 | filter_cells、doublets、HVG |
| I2 | `X` = log1p(normalize_total(counts, 1e4)),float32,在**最终基因空间**上计算 | v0.2 | 嵌入、QC UMAP |
| I3 | `var_names` = 规范基因 symbol(make_unique 去重);**未映射基因默认丢弃**(可关);原名与 mapping 溯源在 `var` | v0.2 | mito/hb 检测、跨数据集合并 |
| I4 | 常规 QC 列存在且为本步**权威计算**:`pct_counts_mt`、`pct_counts_hb`、`total_counts`、`n_genes_by_counts`(scanpy 惯例名)——这是本步写入 obs 的**唯一常规内容** | v0.2 | QC 图、filter_cells、dissect 诊断 |
| I5 | 批次候选全部以真实 obs 列存在(barcode/复合派生列由本步物化)+ 完整排名溯源;**不写 `obs.sample`,不拍板** | v0.3 | doublets 分组、Harmony/MrVI 校正(经 `--batch-col`) |
| I6 | cell_type / organ / tissue 候选排名溯源(只记录,不写规范列) | v0.3 | annotation_qc、atlas、dissect(经显式列参数) |
| I7 | 已通过硬 QC 门,否则不产出 h5ad | v0.1 ✅ | 整条链 |
| I8 | 溯源完整:counts 来源、物种判定、基因映射、丢弃统计、角色排名全部可查(var / uns / result.json) | 全程 | 人 / agent 审查 |

obs 哲学:**原有列原样保留**(batch 候选活在里面),本步只增 I4 的 QC 列
(与 v0.3 的派生批次列);var 列只增不删,行(基因)按 I3 默认裁剪。

## 3. 功能块与分期

| 编号 | 功能 | 期 | 状态 |
|---|---|---|---|
| F1 | 输入校验(合法 h5ad) | v0.1 | ✅ 已交付 |
| F3 | 硬 QC 门(min_cells / min_genes),两道门设计,**先于 F2** | v0.1 | ✅ 已交付 |
| F2 | counts 定位与恢复(stancounts + 三层防线,§5.1) | v0.1 | ✅ 已交付 |
| F4a | 物种解析(F4/F5 的共同前置;四级阶梯,§5.2) | v0.2 | ✅ 已交付 |
| F4 | 基因名统一(stangene;改名 + **默认丢弃未映射**,§5.3) | v0.2 | ✅ 已交付 |
| F5 | QC 指标(物种感知 mt/hb;权威计算,§5.4) | v0.2 | ✅ 已交付 |
| F7 | 构建标准形 + 原子写盘(没有写盘,F4/F5 白算——随本期) | v0.2 | ✅ 已交付 |
| F6 | 元数据角色识别(stanmetacols;**只记录不拍板**,§5.5) | v0.3 | 后置 |
| F8 | 报告物(report.md + qc.png) | v0.3 | 后置 |

验收记录:v0.1 13 项 + v0.2 19 项 = **32 项测试全绿 × 双环境**(Sherlock dl2025
原生 155s;`python:3.12-slim` Apptainer 容器纯本地源码安装 384s)。
T1 推断另有 stangene 侧 9 项单元测试(test_infer_species.py)。

## 4. v0.2 执行流程

```
①  F1       不载入:文件存在 → is_hdf5 → 有 obs/var 结构
②  F3-pre   不载入:h5py 读 shape 取 n_cells;不过 → rejected(毫秒级快拒)
③  载入      read_h5ad
④  F3-pre   n_genes 预门控:X 非零结构可信时判;不可信 → 跳过
⑤  F2       counts 定位与恢复(§5.1);不可恢复 → 永久拒绝
⑥  F3-final 在真 counts 上复核两门 —— 权威判定
⑦  F4a      物种解析(§5.2);定不了 → exit 3 阻塞待决
⑧  F4       基因统一:stangene 映射 → 改名 → 默认丢弃未映射(§5.3);
            丢弃后复核 n_genes_detected 仍过 F3 门,不过 → rejected(final_gate)
⑨  F5       QC 指标:在最终基因空间上权威计算(§5.4)
⑩  F7       构建标准形(counts 层 + lognorm X + var/uns 溯源)
            → 原子写 OUTDIR/standardized.h5ad
⑪           原子写 result.json,按 §7 退出
```

原则:预门控只负责"快拒"(省算力),终门控负责"准"(保正确);同一套阈值。
①–⑥ 为 v0.1 已交付部分,行为不变。

## 5. 功能块细则

### 5.1 F2 · counts 定位:三层防线(步骤内零 LLM)

stancounts 白名单(`counts/count/raw_counts/umi/...`)接不住奇名 layer 时:

1. **确定性扫描 + 一致性校验**:普查所有非排除 layer 的整数性;对整数候选验证
   `log1p(normalize_total(candidate)) ≈ X`(抽样逐行非零集合 + 单一缩放比)。
   唯一通过 → 直接采用,记 `name_recognized: false`。数学确认,非猜测。
2. **显式覆盖**:`--counts-layer NAME`,驾驶员拍板后指定,跳过一切推断。
3. **needs_review 挂钩**:候选歧义 / 逆推但旁有未识别整数 layer / 不可恢复但有
   可疑 layer → result.json 附全部 layer 诊断表,供驾驶员决策后重跑。

counts 采纳后,若来源是奇名 layer,该层在 F7 **改名**为 `counts`(不存两份);
velocity 等其他 layer 保留(列随 F4 基因丢弃同步裁剪)。

### 5.2 F4a · 物种解析:四级阶梯

物种是 F4(基因统一)和 F5(mt/hb)的共同前提。解析顺序:

| 级 | 机制 | 说明 |
|---|---|---|
| T0 | `--species CODE` 显式声明 | 永远最高优先级;批量可复现跑法 |
| T1 | 确定性推断(实现进 stangene:`infer_species`) | 证据投票:Ensembl 前缀(ENSG/ENSMUSG/…,近乎判决)、命名惯例大小写、各物种 mt/hb 基因命中。只写主干规则,不追长尾 |
| T2 | LLM 单次调用(`--llm` 显式开启,默认关) | T1 证据矛盾时:抽样基因名 + T1 计票摘要,一次 structured-output 调用;有 T3 兜底 |
| T3 | 阻塞待决(exit 3) | result.json 给全部证据,驾驶员拍板后带 `--species` 重跑 |

物种确定后,mt/hb 识别是 stangene 精确基因集查询,**不涉及任何猜测**。
stangene 未覆盖的物种 → T3(补参考数据的问题,不是猜的问题)。

### 5.3 F4 · 基因统一:改名 + 默认丢弃未映射

stangene 五级映射瀑布(exact_id → 去版本 ID → approved symbol → alias/prev →
unmapped)逐行注释后,应用策略:

- **改名**:采纳 canonical symbol 为 `var_names`(重名 make_unique 保双列,
  不自动合并);原名存 `var["original_feature_name"]`,mapping 溯源列全部并入
  `var`(I8)。
- **丢弃(默认)**:`mapping_status ∈ {unmapped, ambiguous, non_gene_feature}`
  的基因从矩阵及所有 layers 移除;`--keep-unmapped` 关闭丢弃(保留者维持原名)。
- **统计入账**:result.json 记 per-status 丢弃数与总丢弃比例;
  丢弃比例 > 30% → `needs_review`(非阻塞,exit 0)——通常意味着物种判错或
  非基因型 feature 表。
- **顺序**:丢弃发生在 F5 / F7 之前——QC 分母与 normalize_total 都在最终基因
  空间上计算,全程自洽;丢弃后复核 F3 基因门。

### 5.4 F5 · QC 指标:权威计算

- 在**最终基因空间**(F4 丢弃后)的真 counts 上逐细胞计算:
  `pct_counts_mt`、`pct_counts_hb`(分子 = stangene 物种 mt/hb 基因集在
  harmonized `var_names` 上的精确命中)、`total_counts`、`n_genes_by_counts`。
- **权威性**:数据自带的同名列不作数——写入前同名原列改名保留
  (后缀 `__original`,信息无损),覆盖行为记入 result.json。
- mt/hb 基因集命中数(`n_mt_genes` / `n_hb_genes`)记入 result.json(metrics.qc);
  命中为 0 → 对应列全 0,**不触发 needs_review**——预过滤矩阵(上游已剔 mt 基因)
  与无 RBC 血红蛋白的物种(果蝇/线虫)都是正常情况,驾驶员按需查 metrics.qc。

### 5.5 F6 · 元数据角色识别:识别不拍板(v0.3)

v2 的 normalize_roles(top-1 过 0.5 分 → 拷贝 `sample`/`cell_type_*`/`organ`/
`tissue` 规范列)**整体废弃**——拍板权上移驾驶层。F6 只做两件事:

1. **全角色候选排名**(stanmetacols 十角色;LLM 单次可选、默认关、启发式回退),
   完整写入 `uns["metacols"]` 与 result.json。
2. **物化派生批次列**:barcode 前/后缀、复合列拼接这类"不存在于 obs 的候选"以
   确定性变换写成真实列(如 `batch_from_barcode`)——保证"每个候选都有列可指"。

约定:

- **不写任何规范角色列**;源列一概不动。
- **消费侧显式传列**:doublets、integration、dissect 一律带显式列参数
  (`--batch-col` 等),由驾驶层读排名拍板传入;无默认猜测。batch 的消费者含
  doublets(按批检测,链条上早于 integration)——driver 拍一次板,传两处。
- **iteration 归驾驶层**:integration 的 batch 选择由驾驶层(Agent SDK)按排名
  逐个试 `--batch-col`,以积分指标定案;standardize 不重跑,loop 不进步骤。
- 一个批次候选都没有 → `needs_review`(非阻塞):无批次时 doublets 整体跑、
  integration 不校正,是否接受由驾驶员定。

## 6. CLI 接口(签名即终态,后续只增不改)

```bash
python -m ecasteps.standardize SRC.h5ad -o OUTDIR \
    [--species CODE] [--llm] \
    [--min-cells 100] [--min-genes 5000] \
    [--counts-layer NAME] [--no-gate] [--keep-unmapped]
# 环境引导:bash run.sh standardize <同上参数>
```

- v0.1 产出:仅 `result.json`;**v0.2 起产出 `standardized.h5ad` + `result.json`**。
- 本期新增参数:`--species` / `--llm`(F4a)、`--keep-unmapped`(F4)。
- report.md / qc.png 随 F8(v0.3)加入。

## 7. 退出码与 result.json

| 退出码 | 语义 | 覆盖 | 驾驶员动作 |
|---|---|---|---|
| 0 | ok(含非阻塞 needs_review) | 门全过、counts 可得、物种已定 | 跑下一步;needs_review 另行复核 |
| 2 | 永久性数据问题 | QC 拒绝(含 F4 丢弃后)、counts 不可恢复、非法 h5ad | 跳过,记录,**不重试** |
| 3 | **阻塞待决** | 物种定不了;counts 歧义 | 读 result.json 证据 → 拍板 → 补参数重跑(`--species` / `--counts-layer`)。**agent 驾驶员的标准挂钩点** |
| 1 | 意外错误 | I/O、内存、bug | 重试或人工介入 |

区分标准:2 = 重试与补参数都救不回;3 = **补一个决策就能成**;1 = 原样重试可能就好。
`needs_review` 两种形态:非阻塞(产出成功但有存疑点,exit 0)与阻塞(exit 3,无产出)。

```jsonc
{
  "schema_version": 2,                       // v0.2 起
  "step": "standardize",
  "status": "ok | rejected | needs_review | error",
  "reasons": [],
  "rejected_at": "input | pre_gate | counts_recovery | final_gate | null",
  "src": "/abs/path.h5ad",
  "params": { "min_cells": 100, "min_genes": 5000, "counts_layer": null,
              "species": null, "llm": false, "keep_unmapped": false },
  "species": { "resolved": "human", "code": "hs",
               "source": "cli | inferred | llm | null",
               "confidence": 0.99,
               "evidence": { "ensembl_prefix_hits": {}, "naming_convention": "",
                             "mito_hits": {}, "hb_hits": {} } },
  "metrics": {
    "n_cells": 0, "n_vars": 0, "n_genes_detected": 0,
    "counts_source": "layer:<name> | X | raw | recovered",
    "x_normalization": { "is_log1p": true, "base": null, "is_integer": false },
    "counts_integer": true,
    // 键名带量纲:数的是基因(feature),永远不是细胞
    "harmonization": { "genes_kept": 0,
                       "genes_dropped": { "unmapped": 0, "ambiguous": 0,
                                          "non_gene_feature": 0 },
                       "genes_dropped_frac": 0.0,
                       "genes_unmappable": { "unmapped": 0, "ambiguous": 0,
                                             "non_gene_feature": 0 },
                       "keep_unmapped": false },
    "qc": { "n_mt_genes": 0, "n_hb_genes": 0, "overwritten_obs_cols": [] },
    "timings": { "f1_validate": 0.0, "load": 0.0, "f2_counts": 0.0,
                 "f3_gates": 0.0, "f4a_species": 0.0, "f4_harmonize": 0.0,
                 "f5_qc": 0.0, "f7_build_write": 0.0, "total": 0.0 }
  },
  "layers": [ { "name": "", "dtype": "", "is_integer": true,
                "sparsity": 0.9, "max": 0, "consistent_with_X": null } ],
  // v0.3(F6)加入:全角色排名 + 物化的派生批次列
  "metacols": { "method": "llm | heuristic", "ranking": {}, "derived_columns": [] }
}
```

## 8. 非功能需求

- **可复现**:同输入同参数 → 同输出;默认无网络、无 LLM。
- **幂等**:所有写盘原子化(临时文件 + rename),重跑安全。
- **失败三分**:rejected ≠ error ≠ needs_review,是给驾驶员的正式接口。
- **可移植 / 可发行**:ecasteps 是规范 pip 包(pyproject.toml,py≥3.10;依赖 =
  标准科学栈 + stancounts/stangene/stanmetacols,重依赖走 extras:`[llm]`
  `[doublets]` `[gpu]`)。代码不绑定任何集群。部署形态三选:
  ① Sherlock 原生:`run.sh` 引导 dl2025 venv(环境 fixup 只活在 run.sh);
  ② 任意机器:`pip install`(stan* 经 GitHub / PyPI);
  ③ 容器:Docker 镜像(CPU 与 GPU 分建),Sherlock 上以 Apptainer 消费同一镜像。
  **验收含双环境**:原生 + 容器(v0.1 已建立,`scripts/test-in-container.sh`)。
- 单进程 CPU,无 GPU,内存与数据同量级。

## 9. 验收测试

跑法(双环境,同一套测试):

```bash
bash run.sh test tests -q                    # Sherlock 原生(计算节点)
pip install .[test] && pytest tests -q       # 任意机器
bash scripts/test-in-container.sh            # python:3.12-slim 容器(Apptainer)
```

**v0.1(已过,13 项 × 双环境)**:整数 X / 白名单层 / 奇名层一致性采纳 /
仅 log1p 逆推 / scaled X + counts 层 / `--counts-layer` 覆盖 / `--no-gate` /
逆推旁有可疑层 needs_review / 双候选阻塞 / 指定层缺失阻塞 /
细胞不足预门快拒(未载入)/ 基因不足 / 非 h5ad。

**v0.2 新增**:

| 用例 | 期望 |
|---|---|
| ENSG 前缀基因名 | species=human,source=inferred |
| ENSMUSG 前缀 / 鼠式大小写 symbol | species=mouse,source=inferred |
| 证据矛盾、无 `--species` 无 `--llm` | exit 3,species.evidence 完整 |
| `--species hs` 显式声明 | source=cli,跳过推断 |
| `--llm` 开启且 T1 矛盾(mock LLM) | source=llm;LLM 失败 → 回退 exit 3 |
| 奇名/别名基因 | var_names 改为 canonical,原名在 var |
| 含 unmapped 基因(默认) | 被丢弃,harmonization 统计正确 |
| `--keep-unmapped` | 保留,维持原名 |
| 丢弃比例 > 30% | status=needs_review(exit 0) |
| 丢弃后基因不足 | exit 2,rejected_at=final_gate |
| 已知 MT-/HB- 基因构造 | pct_counts_mt / pct_counts_hb 数值正确 |
| 数据自带 pct_counts_mt 列 | 原列改名 `__original`,新列为权威值,入账 |
| 正常全流程 | standardized.h5ad 满足 I1–I4;counts 源层已改名;写盘原子(中断不留半成品) |
| 重跑同输入 | 输出逐字节可复现(I8 时间戳除外) |

## 10. 代码布局

```
src/ecasteps/
  standardize.py   CLI + 流程 ①–⑪
  countsloc.py     counts 定位:layer 普查 + 一致性证明(围绕 stancounts)
  species.py       物种解析四级阶梯(F4a)
  harmonize.py     基因名统一 + 默认丢弃(F4)
  qc.py            门指标 + 逐细胞 QC 列(F5)
  build.py         标准形装配 + 原子 h5ad 写盘(F7)
  result.py        result.json schema + 退出码(未来各步骤共用)
  atomic_io.py     原子写
tests/             验收测试(原生与容器同一套);dsets.py 构造真实基因名数据
run.sh             Sherlock 环境引导(唯一集群相关文件)
```

## 11. 不做的事(non-goals)

- 非 h5ad 输入(mtx / Seurat)——格式转换是上游步骤(stanobj 领域);
- watch-dir 扫描、样本发现、重试、调度、batch 列拍板、integration 迭代——驾驶层职责;
- 步骤内 Agent loop——**永不**。允许的 LLM 仅两处单次调用(§1 原则 2);
- 多样本批处理——一次调用一个样本,批量由驾驶员循环;
- obs 原始列剥离——如需"干净版",做成驾驶层定案后的可选清理,不是本步默认行为。
