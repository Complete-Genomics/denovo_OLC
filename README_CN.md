# denovo_OLC：面向 SE600 linked reads 的证据感知型逐 UMI 组装

## 摘要

我们提出 [denovo_OLC](https://github.com/Complete-Genomics/LFR_Pipeline/tree/main/modules/clfr/denovo)：一个面向稀疏 SE600 linked-read 池的确定性、逐 UMI OLC 组装器。边界 k-mer 索引、保守的重叠验证、内部锚点恢复和 collective pileup rescue 能在不全局放宽 mismatch 阈值的前提下改善延伸。学习得到的 read 质量预筛在保留图结构决策的同时，降低有界 conflict graph 的计算负担。

在 20,000-UMI benchmark 中，修正 component/merge 语义后，长度 >=1 kb 的 contig 从 16,908 增至 21,139；其中 99.3% 的变长 contig 具有 >=95% 的单层 read-back breadth。对 ML 选择的 draft 进行组装后 Racon polishing，在小规模 UMI 测试中带来 +1.77 个 identity 点的平均增益和 11.41% 的严重增益率，但仍是验证路径而非生产默认。另一个 candidate-chimera GBDT 仅以 shadow sidecar 接入：它记录规则已选择 contig 的风险分数，观测更多样品后决定是否上线。

## 引言

Linked-read barcode 提供了天然的组装单元：分配给一个 UMI 的 reads 构成局部、小规模证据池，而非整个宏基因组。这种结构使逐 UMI 组装具有吸引力，同时也暴露了在常规高深度组装中较不明显的失效模式。在 SE600 数据中，reads 足够长，因而可能遇到噪声末端和 read-through 序列；但单个 barcode 的 read 池也可能仅具有不完整的平铺覆盖证据。在 150–300 万 UMI 的规模下，为每个 barcode 启动通用组装器会被进程启动开销主导。

初始的贪心 OLC 实现快速且适用于低深度 read 池（包括单 read 情况），但其局限是根本性的：两条 read 的直接比较会叠加各自独立的错误。因此，即使完整 read 池支持一条连贯路径，真实 bridge 仍可能无法通过严格的两两 mismatch 测试。de Bruijn graph 在 reads 间汇总 k-mer 支持，对这一失效模式天然更稳健。这里的目标不是为每个 UMI 复刻全局 de Bruijn graph，而是在常规 read 池保留快速、保守 OLC，只在真正停滞的疑难情形部署共享的局部证据。

第二个挑战是 read 选择。一条 read 可能本身质量较低，但移除它是否改善组装，取决于它与*同一 UMI 内其他 reads*之间的不一致。这一区别促成图感知 ML 设计：ML 先移除有限数量的全局低质量 reads，再由 conflict graph 解决剩余的 read-pool 特异性冲突。

## 材料与方法

### 数据与评估队列

分析使用按 barcode 分组的 trimmed SE600 linked-read 数据、具有已知参考的 ZymoBIOMICS mock community，以及缺少完整真值参考的高多样性土壤（`hs6`）队列。主要的受控 OLC benchmark 对所有比较运行均使用同一 trimmed 输入中的前 20,000 个 barcodes。比较始终限定在两个输出中实际都存在的 barcode 交集；比较未匹配的完整数据集是无效的，因而被明确避免。

对于 Zymo，contig 和 reads 使用已知参考集的碱基级 identity 与比对覆盖度评估。对于没有完整参考的数据，使用 read-back breadth、经验证的双端支持以及无参考的 junction 信号。junction 信号统计能跨越 contig 内部位置、并在两侧都具有局部验证序列 identity 的 reads；其目的在于区分真实跨越证据与仅在疑似 chimera 拼接点上表面成立的 k-mer placement。

### Workflow 概览

图 1 概括完整 workflow。实线表示已实现的 SE workflow。ML prefilter 和 Racon 分支是可选的已评估路径，而非无条件生产默认。候选选择默认保持 `longest`；`gated_switch` 是 opt-in。candidate-chimera GBDT 是 shadow-only sidecar：在规则选择之后评分和报告，不存在回写到组装、选择或交付的边。组装之前会先运行异常检测 gate：只有相对 Zymo 呈现不利的深度和 read-length drift，且预测保留深度安全的 SE 样本才会触发 salvage；否则保留原始 SE read 池。

模型重训分支独立存在：`auto_retrain.smk` 不被主生产 workflow include，必须显式调用。调用时，首先运行轻量的 drift report；只有 verdict 为 `drift` 才会启动 candidate training、assembly/read-back A/B 测试和最终 promotion gate。

```mermaid
flowchart LR
    A[带 read quality 的输入 FASTQ] --> B[adapter trimming 与 barcode 分组]
    B --> C[按 barcode 分组的 read TSV]
    C --> D{样本异常检测 gate\n相对 Zymo 的深度与长度 drift\n以及保留深度安全性}
    D -->|正常或仅报告| E[原始逐 UMI read 池]
    D -->|SE salvage| F[移除短于 300 bp 的 reads\n每个 UMI 最多保留前 300 条合格 reads]
    F --> E

    E --> G{可选 ML prefilter\n已有训练模型}
    G -->|关闭或 graph-only 生产路径| H[Conflict graph read-quality filtering]
    G -->|已评估的 scheme B| I[对 reads 评分并移除 candidate pool 中\n评分最低的 K 条 reads]
    I --> H
    H --> J[逐 UMI OLC\n边界 k-mer、验证重叠、\n内部锚点与 collective rescue]
    J --> K[Draft contigs]
    K --> O[全部候选的 junction QC\n仅 gated switch 或 shadow 开启时执行]
    O --> P{候选选择\nlongest 默认；gated switch opt-in}
    P --> Q[规则选择的 contig]

    Q --> L{可选 Racon polishing\npilot / 验证路径}
    E -.->|pilot read pool，最多 50 条 reads| L
    L -->|关闭| M[交付 contigs 与 read-back / junction QC]
    L -->|开启| N[minimap2 read-to-contig 比对\nRacon partial-order consensus]
    N --> M

    O -.-> R[Shadow GBDT：P(chimera)\n仅评分规则选择的 contig]
    P -.-> R
    R -.-> S[Shadow reports：分数分布、\n规则分歧与人工复检队列]

    AA[Mock-control FASTQ 与已知参考<br/>auto_retrain.smk] -.->|显式调用| AB[轻量 feature-drift report]
    AB --> AC{Drift verdict}
    AC -->|no_drift 或 degradation| AD[保留 incumbent\n仅报告 / 调查该 run]
    AC -->|drift| AE[训练 candidate ML model]
    AE --> AF[Incumbent 与 candidate 的\nassembly 和 read-back A/B]
    AF --> AG{Promotion gate}
    AG -->|通过| AH[更新 models/in_use.lgb]
    AG -->|失败| AD
    AH -.->|供下一次生产 run 使用的 model handoff| G

    classDef production fill:#e8f5e9,stroke:#2e7d32,color:#1b5e20;
    classDef optional fill:#fff8e1,stroke:#ef6c00,color:#e65100;
    classDef shadow fill:#e8eaf6,stroke:#3949ab,color:#1a237e;
    class D,E,F,H,J,K,Q,M production;
    class G,I,L,N,O,P,AA,AB,AC,AD,AE,AF,AG,AH optional;
    class R,S shadow;
```

### OLC 组装

对于单个 barcode，先按长度降序和序列排序得到 unique reads，再通过重复的贪心 component 进行组装。候选延伸会在两个方向上评估。只有在共享边界 k-mer prefilter 与完整 suffix-prefix 比较同时满足配置的最小 overlap 和 mismatch-rate 限制时，才接受一个 overlap。因此，k-mer 仅用于发现候选，并非最终正确性判据。

实现采用以下有序证据阶梯：

1. 使用逐 UMI 边界 k-mer inverted index 的 boundary suffix-prefix extension。
2. 针对 literal end 含噪的 reads，使用双向 internal-anchor extension；锚点后必须存在已验证 overlap。
3. 前两条路径停滞后启动 collective rescue。使用多个共线 exact 17-mers 放置 reads，并排除重复锚点；仅当 draft contig 外的碱基获得至少两条独立 reads 支持时才追加。
4. 组装后多数投票 polishing：只有在满足 quorum 和 concordance 支持时才修改碱基。

所有 raw components 均保留到组装后 deduplication 和 merging 完成。用户可见的 `max_contigs` 限制仅在写出最终 contigs 时应用。这一分离至关重要：在 merge 前限制 raw attempts 会静默阻止后续 components 对原本有效的 bridge 作出贡献。

Collective step 可以将早先 raw-component attempt 消耗的 reads 作为证据，但每一个被接受的延伸必须至少消耗一条当前未使用的 read。这一进展不变量可防止在重复性输入中反复进行仅依赖证据、却不消耗 read 的延伸。

### Candidate indexing 与计算复杂度

朴素的贪心 OLC 会反复扫描每条剩余 read，在低质量 read 池中具有约 O(attempts × pool size) 的实际复杂度下限。`denovo_OLC` 每个 UMI 只构建一次 boundary-k-mer index，并只查询能够满足既有强制 k-mer prefilter 的 reads；后续 overlap validation 和 candidate order 保持不变。在一个已记录的病理 barcode 上，这使组装时间从 100.2 s 降至 5.07 s（11.7 倍），且相对未索引实现产生逐字节一致的输出；该优化改变的是搜索成本，不是接受语义。

### k-mer 策略

较短 k-mer 会提高对弱或短 overlap 证据的敏感度，但也增加偶然候选。较长 k-mer 降低 candidate traffic，因此在发现阶段更快、更具特异性；但它可能在真正 overlap 到达最终验证步骤之前就将其遗漏。当前生产 OLC boundary prefilter 使用 k=10；collective rescue 使用严格的 17-mer anchor。实验性的 `k=17 -> k=10` adaptive mode 先运行严格路径，只有当其不产生有效 contig 时才以 k=10 重跑整个 UMI。它不合并 contig pools，也不会仅因结果更长而选择某个结果。

在 2,000-barcode 测试中，k=10、k=17 和 adaptive k=17->10 的 Zymo mean identity 分别为 94.13%、94.15% 和 94.15%；k=10 和 k=17 的土壤 junction-suspect rate 分别为 33.45% 和 33.30%。这些差异处于观测噪声范围内。k=17 在 2,000-UMI 测试中约快 6%，但改变生产默认值仍需覆盖完整深度分布的更大规模评估。

### 图感知 ML read-quality filtering

组装前，将 ungapped overlap identity 低于 0.90 的 reads 在同一 barcode 内连接为 conflict edge。基线 filter 贪心地移除最高 degree vertices，直至不存在 conflict edge。该操作最适合被理解为 read-quality filter：它优先移除 indel/error-rich reads，而非已被证明的污染物或外源分子。

在 ML 实验中，使用 reads 相对 Zymo 已知参考的 identity 训练 LightGBM regressor。特征包括 quality summaries、read length、homopolymer statistics、minimizer spacing、pool depth 以及同一 UMI 内获得其他 reads 支持的 k-mers 比例。数据划分按 barcode 分组，以避免 read-pool level feature leakage。模型五折 MAE 为 2.10，而 global-mean baseline 为 8.48；但这不足以支持 standalone 部署。

评估了三种整合方式：

* **否决：global ML ranking。** 以相同 13.6% 移除率移除全局最低分 reads，mean contig identity 为 94.30%，低于 graph-only filter 的 95.05%。准确预测单条 read identity，并不等于准确选择损害特定 UMI 的 read。
* **有用但次要：图内 ML。** 保留 graph structure 和主要 degree rule。仅对 degree 相同的 vertices 按预测 identity 排序，优先移除预测分数较低的 read。在 3,000 个 held-out barcodes 上，该方法保留 13.58% 的移除率，并将 mean identity 提升至 95.22%。
* **首选：ML prefilter 后接 conflict graph（scheme B）。** 构图前，ML 从每个 candidate pool（最多 25 条 reads）移除评分最低的 8 条 reads；随后未改变的 graph 作出最终 conflict-aware 决策。该方案将 pairwise comparisons 降低 54.7%，并得到 95.28% 的 mean identity，高于同一评估中的 graph-only filtering（95.05%）和 ML tie-breaking（95.22%）。二者分工互补：ML 移除全局低质量 reads，graph 识别移除哪条 read 可以解决特定的 UMI 内冲突。

Scheme B 是当前首选的已评估 operating point，但尚非无条件生产默认：它移除了 24.48% 的 candidate reads，而 graph-only baseline 为 13.59%，且尚未在更多样本类型上复现。

### Racon polishing 与析因消融实验

在 polishing pilot 中，从原始逐 UMI read 池中最多取 50 条 reads，用 `minimap2 -x sr` 回贴到每条 draft contig，再交给 Racon 进行 partial-order consensus。若 contig 没有可用的 read-to-contig alignment，则原样透传。实验采用 2 × 2 析因设计：plain 或 ML tie-break draft，分别有或没有 Racon。因此，Racon(ML-draft)-versus-plain 的组合比较是一条端到端实验臂，而四个实验臂共同构成了可分离各自贡献及交互效应的消融实验。该 pilot 评估的是 ML tie-break draft，而非 scheme B。

### 候选选择与 shadow 监控

生产兼容的候选选择器默认保持历史 `longest` 行为。`candidate_select: gated_switch` 是 opt-in：它按 rank 选择第一个通过配置 junction gate（`span_cov_ratio >= 0.25` 且有足够 placed reads）的候选；若没有候选通过，则回退 primary candidate。该规则尚未成为默认，因为未 polish 的 severe-loss 在不同验证队列中并不一致。

配置 `shadow_score_model` 后，全部候选的 junction QC 与选择记录会输入冻结的 GBDT。模型只对规则实际交付的 contig 打分并写监控报告；没有任何下游 workflow rule 读取该分数。hs8 的 3,000-UMI integration run 已完成这份完整契约（2,964 个进入组装的 UMI）；GBDT 与 gated rule 在 89.71% 的交付上分歧。这个分歧是监控信号，不是模型应接管选择的证据。Shadow 的作用是在无参考真值的新批次上积累分数分布、规则分歧和可人工复检的高风险案例。

## 结果

### 正确的 component 语义恢复了有 read 支持的序列

在固定的 20,000-UMI benchmark 中，将 contig cap 移至最终输出、并让所有 raw components 都进入组装后 merging，使总 contig 数从 69,596 增至 124,240，长度 >=1 kb 的 contig 数从 16,908 增至 21,139。有 2,349 个 barcode 产生了更长的 primary contig。在这些变长 contig 中，2,333 个（99.3%）具有 >=95% 的单层 read-back breadth，2,312 个（98.4%）具有 >=80% 的双层 breadth。这些结果支持该语义修正，但不能证明每条长延伸都正确；仅单端支持和 breadth 回退队列仍是明确的审查对象。

### Collective rescue 缩小、但未消除 OLC 与图组装之间的差距

在 5,720 个具有 >=1 kb MEGAHIT 输出的匹配 barcode 集合中，collective rescue 后 OLC primary contig 至少与 MEGAHIT contig 同长的比例由 23.0% 增至 68.2%。这是结构比较而非真值主张：无论 OLC 更长还是 MEGAHIT 更长，都需要 raw-read evidence。代表性失败案例显示，MEGAHIT 路径可具有连续、全长的 raw-read coverage，而两两 OLC 仍接近 read 长度。Collective pileup 在不全局放宽 mismatch threshold 的前提下挽回了部分缺口，但不会使贪心 OLC 等同于 de Bruijn graph。

### ML prefilter 加 conflict graph 提供最佳的已测质量–速度权衡

Graph-only quality filter 将 Zymo mean contig identity 从未过滤的 94.16% 提升至 95.05%。尽管在 surrogate label 上具有较强的交叉验证表现，stand-alone ML filter 在相同移除率下明显更差。最佳测试策略为 scheme B：从 capped pool 中以 ML 移除评分最低的 8 条 reads，再应用未改变的 conflict graph。它在 24.48% 的 read-removal rate 下，将 pairwise comparisons 降低 54.7%，并将 mean identity 提升至 95.28%。在 2,997 个匹配 barcode 中，相对 graph-only filter，scheme B 在 805 个（26.9%）上更好、516 个（17.2%）上更差、1,676 个（55.9%）上相同。

图内 ML tie-breaking 同样为正向结果（13.58% 移除率下 identity 为 95.22%），但 scheme B 在该评估中更可取，因为它以略高的 identity 移除了超过一半的 pairwise graph 工作量。该结果支持一个具体结论：ML 在增强分子证据图时可以改善 read filtering，而不是用彼此独立的全局分数取代局部结构。

### 析因消融显示 ML 与 Racon 的增益大多可叠加

Racon(ML-draft)-versus-plain 是端到端比较，但单独来看不是消融实验：它同时改变了 draft selection 和 polishing。完整的 2 × 2 Zymo 实验使用 4,998 个具有已知参考评分的匹配 barcode，消除了这一歧义：

| Draft 选择 | Polishing | Mean identity | 相对 plain raw 的差异 |
|---|---|---:|---:|
| plain | 无 | 94.9698 | -- |
| ML tie-break | 无 | 95.1632 | +0.1934 |
| plain | Racon | 96.7406 | +1.7708 |
| ML tie-break | Racon | 96.9022 | **+1.9324** |

Racon(ML draft) 比 Racon(plain draft) 高 0.1616 个 identity 点（Wilcoxon p = 6.87e-13）。因此，约 84% 的原始 ML 增益在 polishing 后仍然保留；重叠收益仅为 0.0318 点。该结果支持二者具有互补机制：ML 改变形成 draft 的 read 集合，而 Racon 通过 read-to-contig realignment 和 partial-order consensus 纠正碱基。

同一组合臂还在 2,970 个匹配 hs8 barcode 上，以 fixed-reference local-truth 评估测试。相对 plain raw，Racon(ML draft) 将 mean identity 提升 1.7675 点；但其 1.01% 的 severe-loss rate 实际上正处于预先规定的 1% 安全边界。Racon 因此是有前景的组装后阶段，而非当前生产默认；启用该组合路径前仍需要更大规模 field-sample 验证。

## 讨论

`denovo_OLC` 被刻意设计为混合证据系统。常规 barcode read 池采用快速、确定性的 OLC 路径；困难 read 池获得保守的 internal-anchor 和 multi-read rescue，而非全局放宽 mismatch threshold。这既保持速度，也限制缺乏支持的拼接机会；但同时保留一项已知局限：局部证据聚合不是完整 de Bruijn graph，无法恢复每一条被错误遮蔽的路径。

ML 结果展示了类似的设计原则。模型可以准确预测 read identity，却仍可能在忽略局部 conflict structure 时作出不良组装决策。最佳的已测分工是 scheme B：ML 先移除固定数量的全局低质量 reads，随后由 graph 保留对局部冲突解决的决策权。这同时改善了测得的 identity 和 pairwise cost，但其较高的 read-removal rate 和有限的跨样本验证意味着应谨慎使用。

若干看似合理的干预已被测试并否决：整 UMI 移除、全局 mismatch relaxation、single-minimizer deduplication、homopolymer compression 和 standalone ML ranking。这些负面结果是方法的一部分，而非遗漏；它们将后续开发约束在保留证据的局部 graph/bridge rescue、更好的 indel-aware validation，以及跨 library type 的外部验证上。

## 代码与数据可用性

维护中的实现和 workflow 位于 [LFR pipeline denovo module](https://github.com/Complete-Genomics/LFR_Pipeline/tree/main/modules/clfr/denovo)。本仓库是该方法面向 release 的说明。Benchmark 细节、失败假设和决策记录将写入 [`tech_notes.md`](tech_notes.md)。
