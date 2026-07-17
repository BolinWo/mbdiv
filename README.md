# mbdiv — 宏基因组多样性分析流程

[![Python 3.9+](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![pyflakes clean](https://img.shields.io/badge/pyflakes-clean-green.svg)](https://github.com/PyCQA/pyflakes)

一条命令完成物种级宏基因组多样性分析。

```
raw_data  ──►  物种合并  ──►  归一化  ──►  Alpha 多样性（4 指数 + 箱线图）
                                  │
              Beta 多样性  ◄──────┘  （Bray-Curtis + PCoA + PERMANOVA）
                                  │
              过滤病毒/真菌  ──►  Top-N 物种堆叠图 + 热图
```

## 快速开始

```bash
# 安装
pip install mbdiv-1.0.0-py3-none-any.whl[full]

# 运行 — 就这么简单
mbdiv species.xlsx metadata.xlsx
```

所有结果（数据表、统计量、图表）自动输出到 `./result/` 目录。

## 安装方式

### 方式一：Wheel 安装（推荐）

```bash
# 仅核心功能（Alpha、Beta、Top-N，不含 scikit-bio）
pip install mbdiv-1.0.0-py3-none-any.whl

# 完整功能（scikit-bio PCoA + YAML 配置支持）
pip install mbdiv-1.0.0-py3-none-any.whl[full]
```

### 方式二：源码安装

```bash
git clone https://github.com/BolinWo/mbdiv.git
cd mbdiv
pip install .[full]
```

### 方式三：便携 ZIP（离线环境）

```bash
# 解压 mbdiv-1.0.0-portable.zip
# Windows: 双击 setup.bat
# Linux/Mac: bash setup.sh
```

## 使用方法

### 基本用法

```bash
# 最简模式 — 自动检测列名、样本前缀、分组
mbdiv species.xlsx metadata.xlsx

# 指定输出目录和 Top-N 数量
mbdiv species.xlsx metadata.xlsx -o my_results -t 15
```

### 全部参数（单字母旗标）

| 旗标 | 完整形式 | 默认值 | 说明 |
|------|----------|--------|------|
| *(位置参数)* | — | — | 物种丰度 Excel 文件 |
| *(位置参数)* | — | — | 样本元数据 Excel 文件 |
| `-r` | `--raw` | — | 物种丰度 Excel（位置参数的替代写法） |
| `-m` | `--meta` | — | 样本元数据 Excel（位置参数的替代写法） |
| `-o` | `--output` | `result` | 输出目录 |
| `-t` | `--top` | `10` | Top-N 物种数量（堆叠图） |
| `-c` | `--config` | — | YAML 配置文件（高级选项） |
| `-h` | `--help` | — | 显示帮助信息 |

### 两种调用风格

```bash
# 风格一：位置参数（简洁）
mbdiv species.xlsx metadata.xlsx -o results -t 10

# 风格二：显式旗标（清晰）
mbdiv -r species.xlsx -m metadata.xlsx -o results -t 10
```

### 高级：YAML 配置文件

如需可重复分析，可将设置保存为 `config.yaml`：

```bash
mbdiv species.xlsx metadata.xlsx -c config.yaml
```

参考 `assets/config_template.yaml` 了解全部可配置项，包括：
- 列名覆盖（物种列、分类学列、样本前缀）
- 分组顺序与颜色
- Alpha 多样性指标与统计检验方法
- Beta 距离度量与 PERMANOVA 置换次数
- 病毒/真菌过滤关键词
- 图表格式（PDF/PNG/两者）与 DPI

## 输入格式

### raw_data.xlsx — 物种丰度表

| #Species | Sample_01 | Sample_02 | ... | #Taxonomy |
|----------|-----------|-----------|-----|-----------|
| Escherichia coli [TAX_001] | 152.3 | 89.1 | ... | d__Bacteria;p__Proteobacteria;... |
| Bacteroides fragilis [TAX_002] | 203.5 | 0.0 | ... | d__Bacteria;p__Bacteroidota;... |

- **第一列**：物种名称（自动检测：`#Species`、`Species` 或任何以 `#` 开头的列）
- **中间列**：样本丰度值（数值型，自动检测）
- **最后一列**：分类学字符串（自动检测：包含 "tax" 的列）

### meta_data.xlsx — 样本元数据

| Sample | Group | Score |
|--------|-------|-------|
| Sample_01 | Control | 0 |
| Sample_02 | Treatment | 3 |

- **Sample** 列：样本 ID，需与 raw_data 列名匹配
- **Group** 列：分组信息，用于组间比较
- **数值列**（可选）：用于与 Alpha 多样性做 Spearman 相关性分析

**所有列名均自动检测。** 除非格式特殊，无需手动指定。

## 输出结构

```
<输出目录>/
├── pipeline_report.txt              # 完整流程审计报告
├── data/                            # 中间数据表
│   ├── merged_species_clean.xlsx    #   步骤1：合并 + 去零后
│   ├── zero_abundance_taxa.xlsx     #   步骤1：被移除的分类记录
│   ├── relative_abundance.xlsx      #   步骤2：归一化（小数）
│   ├── relative_abundance_percent.xlsx  # 步骤2：归一化（百分比）
│   └── sample_rpkm_summary.xlsx     #   步骤2：测序深度统计
└── result/
    ├── step3_alpha/                 # Alpha 多样性
    │   ├── alpha_diversity.xlsx     #   4 指数 × N 样本
    │   ├── alpha_statistics.xlsx    #   Kruskal-Wallis p 值
    │   ├── alpha_spearman.xlsx      #   Spearman 相关性（如有数值列）
    │   └── figures/
    │       ├── Observed_boxplot.{pdf,png}
    │       ├── Shannon_boxplot.{pdf,png}
    │       ├── Simpson_boxplot.{pdf,png}
    │       └── Chao1_boxplot.{pdf,png}
    ├── step4_beta/                  # Beta 多样性
    │   ├── distance/
    │   │   └── braycurtis_distance.xlsx
    │   ├── pcoa/
    │   │   ├── pcoa_coordinates.xlsx
    │   │   └── pcoa_variance.xlsx
    │   ├── statistics/
    │   │   └── PERMANOVA_result.txt
    │   └── figures/
    │       └── PCoA_BrayCurtis.{pdf,png}
    └── step5/                       # Top-N 物种组成
        ├── bacteria_species.xlsx
        ├── bacteria_relative_abundance.xlsx
        ├── top10_species.xlsx
        ├── top10_group_percentage.xlsx
        └── figures/
            ├── top10_group_bar.{pdf,png}
            └── top10_individual_heatmap.{pdf,png}
```

## 分析流程

### 步骤 1：合并与清洗
- 同名物种合并（丰度求和）
- 移除全零丰度分类
- 提取简洁物种名

### 步骤 2：归一化
- 按样本计算相对丰度（列和归一化为 1.0）
- 同时保存小数和百分比两个版本
- 测序深度统计

### 步骤 3：Alpha 多样性
- **Observed** — 物种丰富度
- **Shannon** — H' 指数（自然对数）
- **Simpson** — Gini-Simpson（1-D）
- **Chao1** — 丰富度估计量
- Kruskal-Wallis 组间检验
- 可选：与数值型元数据的 Spearman 相关性
- 箱线图 + 抖动点可视化（PDF + PNG）

### 步骤 4：Beta 多样性
- Bray-Curtis 距离矩阵
- PCoA 排序（scikit-bio，不可用时降级为 SVD）
- PERMANOVA 检验（999 次置换）
- PCoA 散点图 + 95% 置信椭圆
- 动态坐标轴标签（方差解释百分比）与 PERMANOVA p 值标注

### 步骤 5：Top-N 物种组成
- 根据分类学关键词过滤病毒和真菌
- 在仅细菌数据上重新计算相对丰度
- 按平均丰度选取 Top-N 物种
- 分组水平堆叠柱状图
- 个体样本热图

## 依赖

| 包名 | 是否必需 | 用途 |
|------|----------|------|
| pandas | 是 | 数据处理 |
| numpy | 是 | 数值计算 |
| scipy | 是 | 统计 + 距离计算 |
| matplotlib | 是 | 绘图 |
| seaborn | 是 | 统计可视化 |
| openpyxl | 是 | Excel 读写 |
| scikit-bio | 可选 | PCoA + PERMANOVA（不可用时降级为手动 SVD） |
| pyyaml | 可选 | YAML 配置文件支持 |

## 可重复性

流程在输出目录生成 `pipeline_report.txt`，包含：
- 输入文件路径
- 全部输出文件路径
- 关键统计结果（PERMANOVA p 值等）
- 总运行时间

如需精确复现，请保存 YAML 配置并使用相同的输入文件。

## 许可证

MIT — 详见 [LICENSE](LICENSE)

## 引用

如果在研究中使用了 mbdiv，请引用：

```
mbdiv: 单命令宏基因组多样性分析流程
https://github.com/BolinWo/mbdiv
```
