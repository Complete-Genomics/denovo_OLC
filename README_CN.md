# denovo_OLC

[LFR pipeline denovo](https://github.com/Complete-Genomics/LFR_Pipeline/tree/main/modules/clfr/denovo) 模块的独立发布版本：一个进程内运行的Overlap-Layout-Consensus（OLC）组装器，专为**单端600bp（SE600）linked-read数据**设计，在百万分子级别的规模上为每个UMI/barcode重建出一条fragment。作为per-UMI de Bruijn图组装（例如每个barcode fork一次MEGAHIT子进程）的直接替代方案，本工具是围绕SE600数据——相对于短的双端read——实际呈现出的噪音特征和read结构特点专门设计的。

## 为什么SE600需要一个专门设计的组装器

Complete LFR（cLFR）以及类似的barcode分区测序技术，每个物理DNA分子会产生一小簇共享同一个barcode的短read。现有的大多数per-barcode组装工具都是围绕短的（~150bp）双端read设计的，在这种场景下mate-pair信息有助于消歧重复区域，adapter read-through也很少深到需要处理的程度。**SE600数据是完全不同的场景**：单端、约600bp长，没有mate可以交叉核对，read长度足以经常读穿比它短的物理插入片段，一路读进ligation adapter、样本index，甚至barcode本身；而且是在百万级barcode/次运行这样的规模上组装的——在低通量场景下算是罕见边缘情况的噪音，在这个规模下变成了常规的、高频出现的失败模式。由此直接产生两个问题：

- **De Bruijn图组装要求k-mer覆盖度≥2。** 低深度的barcode（在任何真实的barcode分区数据集里都占不小比例）根本无法组装——图会碎裂或者根本建不起来。
- **每个barcode fork一次完整的组装器子进程，规模上跑不动。** 在每次运行1.5-3百万个barcode的量级下，仅子进程fork本身的开销就会成为主导成本，跟实际的组装计算量无关。

`denovo_OLC`把这两个问题都解决了：直接在进程内实现经典的OLC范式（overlap检测 → layout → consensus），纯Python实现，不依赖任何外部组装器，也不需要每个barcode fork一次子进程——甚至能从**单条read**开始优雅地组装出结果——同时把SE600特有的噪音（边界噪音、多条独立read错误叠加、adapter/index/barcode read-through）当作需要检测和处理的一等问题，而不是可以忽略的边缘情况。

## 核心创新点

### 1. 基于边界k-mer倒排索引的输出敏感型候选发现

朴素的OLC做法是每次contig延伸都重新扫描整个read池（这在pairwise merge成功率低的barcode上会退化成O(尝试次数 × read池大小)），`denovo_OLC`则为每个UMI建一次"边界k-mer → 候选read"的倒排索引，每一步延伸都直接查询这个索引。在一个记录在案的真实最坏情况barcode上，这把单个UMI的组装时间从**100.2秒降到5.07秒（11.7倍）**，同时产出的组装结果跟不建索引的baseline逐字节完全一致（通过对隔离测试和500-barcode批量测试的摘要值相等来验证）——索引只是剪枝搜索空间，从不改变哪些overlap会被接受。

### 2. 双向内部k-mer锚定

真实read在自己的5′/3′端天然携带噪音（adapter read-through、read末端测序错误富集）。朴素的OLC只尝试匹配read的**边界**，即便再往里几个碱基就存在一个干净、无歧义的overlap，也会因为这类端部噪音而彻底延伸失败。`denovo_OLC`加入了内部锚点兜底机制——同时建了正向和反向索引——定位一个经过验证的内部锚点并从那里开始延伸，找回边界匹配完全会漏掉的合并。

### 3. 带跨尝试证据共享的集体pileup救援

Pairwise的单read overlap比较有一个结构性弱点：两条各自携带独立测序错误的read直接比较时，会把两边的错误率叠加在一起，经常把一个本来有效的overlap推过mismatch率阈值。与其全局放宽这个阈值（这样会给chimeric误组装开口子），`denovo_OLC`的做法是——只在常规pairwise延伸卡住时才触发——用长的、低碰撞率的k-mer同时锚定多条read，只延伸有pileup（多read）支持的区域。证据可以跨同一个barcode的不同贪心组装尝试共享，不局限于单次尝试内部，让被早前某次尝试消耗掉的read依然能为后一次尝试提供支持投票。一个"进度不变量"（每次被接受的延伸都必须消耗至少一条此前未使用过的read）保证了即便在对抗性重复输入下也能终止。

### 4. 轻量级多数投票consensus抛光

组装完成后，一个基于k-mer偏移量的投票步骤会把每条输入read重新对齐到草稿contig上（不做完整的动态规划重新对齐），只有在存在多数、达到法定票数支持的分歧时才翻转某个位置——纠正贪心延伸过程中引入的替换错误，而不需要完整多序列比对consensus步骤的成本和复杂度。

### 5. 构造上即确定性的输出

Python每个进程的字符串哈希随机化，意味着朴素地用`set()`做去重会让种子选择——进而整个贪心组装轨迹——在*相同*输入的不同运行之间悄悄变得不可复现。`denovo_OLC`里每一个对顺序敏感的步骤都使用显式的、基于内容的确定性tie-break，所以相同输入总是产出相同输出，跟解释器哈希种子、进程启动方式、worker数量都无关。

### 6. Adapter/index/barcode read-through检测

这正是SE600读长的双刃剑之处：一条600bp的单端read经常会读穿比它短的物理插入片段，一路读进ligation adapter和样本index，在（DNB风格的）滚环扩增文库制备中，甚至会读回barcode本身。如果不trim掉，这段技术性的read-through会被当作生物学序列一起组装进去，虚高contig的表观长度，把本应是单一连续的组装拆散成碎片。`denovo_OLC`的validation工具链包含了不依赖参考序列的取证方法（位置不变性测试、丰度匹配的保守区域对照、以及跟参考序列数据库的交叉核对），用来区分真正的adapter污染和真实的生物学保守区域——这对于正确归因这一SE600特有的长度/碎片化artifact类别至关重要，避免把它误诊为组装本身的bug。

### 7. 生产级别的多进程正确性

正确处理了`fork`（Linux默认）和`spawn`（macOS/Windows默认）多进程启动方式之间的语义差异——配置状态被显式传播给每个worker，而不是依赖`fork`的写时复制继承机制（这在`spawn`下会悄悄失效）。任务调度针对生物学数据固有的、高度不均匀的per-barcode成本分布做了调优（大多数barcode成本是毫秒级，一小部分尾部成本要高出几个数量级），避免了朴素固定chunksize调度器在这种分布下产生的worker饥饿问题。

### 8. 通过验证跨越深度崩塌检测chimera，不需要参考序列

一个UMI的read池可能混入了不止一个来源分子——这是组装之前发生的barcode碰撞或交叉污染，不是组装本身的bug——一旦发生，assembler自己的overlap逻辑可能会靠一段共享的短motif把两者桥接起来，产出一条看起来完全健康的chimeric contig：文件里每条read确实都属于这个barcode，所以朴素的read-back一致性检查根本看不出问题。`denovo_OLC`转而检查**经过验证的跨越深度**：对contig内部每个位置，统计有多少read能以真实的、两侧都成立的序列identity跨越它，而不只是k-mer能不能placement成功。真正的chimera拼接点会表现为验证过的跨越支持急剧崩塌，而普通覆盖度依然健康——因为没有任何一个物理分子能同时跨越两个不同物种的拼接点两侧。在一个参考序列完全已知的mock community对照上实测，这个信号对chimera的判别力达到AUC 0.827，而同一个read池上朴素的read-back一致性检查和一个现成的、不需要参考库的chimera检测工具，表现都接近随机。

### 9. 按样本多样性自适应调整QC强度，不需要用户输入样本类型

"可接受的fragment质量/产出权衡"很大程度上取决于样本本身：低多样性群落和高多样性环境样本，在同样的QC严格度下要付出的代价完全不同。与其让操作者声明样本类型，一个轻量级探针会直接测量这次run自己的cross-barcode序列identity分布——不同barcode的read按定义就是不同的来源分子，所以这个分布本身就是一个不需要外部参考库的内建阴性对照——并据此自动选择合适的QC预设，决策过程会被记录下来供审计。每个预设的fragment产出/准确度权衡都是**实测**出来的，不是假设的。

## 设计理念

- **不依赖任何外部组装器。** 纯标准库——没有编译扩展，没有子进程管理，没有需要考虑的非确定性第三方图算法库行为。
- **正确性优先于激进的启发式方法。** 每一个兜底机制（内部锚定、集体救援）都被刻意安排在更安全、更保守的选项成功或用尽之*后*才触发，每一项优化在被采纳前都要通过固定输出摘要值的回归验证。
- **在真实数据上验证过，不只是合成基准测试。** 合成测试套件（覆盖overlap检测、chimeric合并安全性、确定性tie-break、多进程配置传播等等的单元测试）跟针对真实测序数据的取证分析配合使用，用来捕捉合成测试在结构上无法覆盖的失败模式——不对称的read长度分布、平台特有的adapter构造、以及真实存在的病态深度离群值。
- **准确度是测出来的，不是假设出来的。** 用一个参考序列完全已知的mock community直接测量per-base identity和chimera率，而不是只靠群落组成"看起来合理"这种间接判断，每个QC预设的产出/准确度权衡都来自这个实测。好几个听起来很合理的改进方案——针对疑似read混杂的整UMI剔除、indel容忍的overlap打分、homopolymer-aware组装——都实现过，但在留出的ground truth上测不出实际收益后被否决，不是靠直觉留下来的。

## 用法

```bash
python3 src/denovo_seed_olc.py \
  --sequence_type se \
  --num_processes 30 \
  --min_ctg_len 400 \
  --r2 data_R2_sgrep.tsv
```

完整的配置项（overlap/mismatch阈值、内部锚定和集体救援开关、抛光参数、输出上限）请参见模块docstring。

## 状态

正在生产环境中用于百万UMI级别规模的per-UMI 16S rRNA amplicon和宏基因组fragment重建，生产环境的QC层（组装前read质量过滤、组装后chimera检测、按多样性自适应的QC预设）已在一个参考序列完全已知的mock community上完成验证。