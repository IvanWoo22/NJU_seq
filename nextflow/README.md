# NJU-seq Nextflow Pipeline

## 📖 项目简介

NJU-seq Nextflow Pipeline 是一个基于 [Nextflow](https://www.nextflow.io/) 框架的自动化工作流，用于从 RNA-seq 数据中全面识别 2'-O-甲基化（Nm）修饰位点。

## 🌟 核心特性

- **模块化设计**: 每个分析步骤都是独立的模块，便于维护和扩展
- **并行执行**: 自动识别可并行化的任务，最大化利用计算资源
- **容器化支持**: 支持 Docker 和 Singularity 容器，确保可重复性
- **云就绪**: 原生支持 AWS Batch、Google Cloud Life Sciences 等云平台
- **详细日志**: 完整的日志记录和错误追踪
- **灵活配置**: 通过配置文件轻松调整参数和资源

## 📋 依赖要求

### 必需软件

- **Nextflow** (>= 21.04.0)
- **Java** (>= 11)
- **Perl** (>= 5.26)
- **R** (>= 4.0)
- **bowtie2** (>= 2.4)
- **samtools** (>= 1.12)
- **GNU parallel**

### R 包依赖

```R
install.packages(c(
    "ggplot2", "ggpubr", "gridExtra", "forcats", "dplyr",
    "VennDiagram", "splines", "RColorBrewer", "extrafont", "reshape2"
))
```

### Perl 包依赖

```bash
cpanm YAML::Syck AlignDB::IntSpan PerlIO::gzip Algorithm::Combinatorics
```

### 可选：容器

```bash
# Docker
docker pull ivanwoo/nju-seq:latest

# Singularity
singularity pull nju-seq.sif docker://ivanwoo/nju-seq:latest
```

## 🚀 快速开始

### 1. 安装 Nextflow

```bash
# 方法1: 使用 conda
conda install -c bioconda nextflow

# 方法2: 直接下载
curl -s https://get.nextflow.io | bash
mv nextflow /usr/local/bin/

# 方法3: 使用 brew
brew install nextflow
```

### 2. 准备参考文件

```bash
# 创建目录结构
mkdir -p NJU_seq_reads NJU_seq_output index

# 下载参考序列和注释文件
# (参考 README.md 获取详细下载说明)
```

### 3. 配置参数

编辑 `nextflow.config` 或创建自定义配置文件：

```groovy
// my_config.config
params {
    samples = ['SAMPLE1', 'SAMPLE2', 'SAMPLE3']
    output_dir = 'results'
    threads = 24
}
```

### 4. 运行流程

```bash
# 标准运行
nextflow run main.nf

# 使用自定义配置
nextflow run main.nf -c my_config.config

# 使用 Docker 容器
nextflow run main.nf -with-docker ivanwoo/nju-seq:latest

# 使用 Singularity 容器
nextflow run main.nf -with-singularity nju-seq.sif

# 恢复中断的运行
nextflow run main.nf -resume

# 查看执行计划（不实际运行）
nextflow run main.nf -dry-run
```

## 📁 项目结构

```
nextflow/
├── main.nf                 # 主流程文件
├── nextflow.config         # 配置文件
├── modules/                # 分析模块
│   ├── bacterial_alignment.nf
│   ├── bacterial_filter.nf
│   ├── remove_bacterial_reads.nf
│   ├── rrna_alignment.nf
│   ├── rrna_filter.nf
│   ├── rrna_count.nf
│   ├── rrna_scoring.nf
│   ├── remove_rrna_reads.nf
│   ├── mrna_alignment.nf
│   ├── mrna_filter.nf
│   ├── mrna_dedup.nf
│   ├── mrna_count.nf
│   ├── mrna_merge.nf
│   ├── mrna_scoring.nf
│   ├── regional_analysis.nf
│   ├── motif_analysis.nf
│   ├── visualization.nf
│   └── report.nf
├── conf/                   # 配置文件
├── results/               # 输出结果
└── reports/               # 分析报告
```

## 🔧 配置选项

### 参数配置

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `samples` | 见配置文件 | 样本ID列表 |
| `output_dir` | `NJU_seq_output` | 输出目录 |
| `threads` | 24 | CPU线程数 |
| `bowtie_threads` | 22 | Bowtie2线程数 |
| `parallel_jobs` | 20 | 并行任务数 |
| `min_score` | 0.5 | 最小评分阈值 |
| `min_coverage` | 5 | 最小覆盖度 |

### 资源限制

可以在 `nextflow.config` 中为不同类型的任务设置资源限制：

```groovy
process {
    withLabel: 'mrna_align' {
        cpus = 24
        memory = 12.GB
        time = 6.h
    }
}
```

## 📊 分析流程

```
┌─────────────────────────────────────────────────────────────┐
│                    Stage 1: 数据预处理                        │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐   │
│  │ 原始数据质控  │ -> │ 细菌RNA过滤  │ -> │ 数据清洗     │   │
│  └──────────────┘    └──────────────┘    └──────────────┘   │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│                    Stage 2: rRNA分析                         │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐   │
│  │ rRNA比对    │ -> │ rRNA过滤     │ -> │ rRNA计数评分 │   │
│  └──────────────┘    └──────────────┘    └──────────────┘   │
│                            ↓                                 │
│                    ┌──────────────┐                         │
│                    │ rRNA reads移除│                         │
│                    └──────────────┘                         │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│                    Stage 3: mRNA分析                         │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐   │
│  │ mRNA比对     │ -> │ 比对过滤     │ -> │ 去重复      │   │
│  └──────────────┘    └──────────────┘    └──────────────┘   │
│  ┌──────────────┐    ┌──────────────┐                       │
│  │ 计数合并     │ <- │ 唯一性判断   │                       │
│  └──────────────┘    └──────────────┘                       │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│                    Stage 4: Nm位点评分                       │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐   │
│  │ 样本评分     │ -> │ 区域分布分析 │ -> │ Motif分析   │   │
│  └──────────────┘    └──────────────┘    └──────────────┘   │
│                            ↓                                 │
│                    ┌──────────────┐                         │
│                    │ 可视化展示   │                         │
│                    └──────────────┘                         │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│                    Stage 5: 报告生成                         │
│  ┌──────────────┐    ┌──────────────┐                       │
│  │ HTML报告    │ <- │ 统计摘要     │                       │
│  └──────────────┘    └──────────────┘                       │
└─────────────────────────────────────────────────────────────┘
```

## 🎯 输出结果

### 样本级结果

每个样本在 `${output_dir}/${sample}/` 目录下包含：

```
${sample}/
├── bac_rna.raw.sam              # 细菌RNA比对原始结果
├── bac_rna.bowtie2.log          # 细菌RNA比对日志
├── bac_rna.out.list             # 细菌RNA过滤列表
├── rrna.filter.sam              # rRNA比对结果
├── rrna.filter.bowtie2.log     # rRNA比对日志
├── rrna.filter.list            # rRNA过滤列表
├── rrna.filter.{25s,18s,5-8s}.tsv  # rRNA计数结果
├── mrna.filter.sam             # mRNA比对结果
├── mrna.filter.bowtie2.log     # mRNA比对日志
├── mrna.dedup.filter.tmp       # 去重后结果
└── mrna.count.filter.tmp       # 计数结果
```

### 汇总结果

在 `${output_dir}/` 目录下：

```
NJU_seq_output/
├── mrna.filter.Ath_leaf_RF.tsv      # 叶组织Nm位点
├── mrna.filter.Ath_root_RF.tsv      # 根组织Nm位点
├── mrna.filter.Ath_stem_RF.tsv       # 茎组织Nm位点
├── mrna.filter.Ath_flower_RF.tsv     # 花组织Nm位点
├── rrna.filter.25s.Ath_leaf_RF.tsv
├── rrna.filter.18s.Ath_leaf_RF.tsv
├── rrna.filter.5-8s.Ath_leaf_RF.tsv
├── figures/                           # 可视化结果
│   ├── *.signature.pdf              # 签名计数图
│   └── *.region.pdf                 # 区域分布图
└── reports/                          # 分析报告
    ├── nju_seq_report.html          # HTML报告
    ├── nju_seq_report.md            # Markdown报告
    └── summary_statistics.json      # JSON统计摘要
```

## 🔍 监控和调试

### 查看运行状态

```bash
# 实时查看日志
nextflow log -f process,status,duration

# 查看已完成的任务
nextflow log -f process,hash,status,exit

# 查看特定运行的详细信息
nextflow log <run_name>
```

### 调试选项

```bash
# 显示详细日志
nextflow run main.nf -ansi-log false

# 生成执行图
nextflow run main.nf -with-dag flowchart.png

# 单步执行
nextflow run main.nf -N 1
```

## 🐛 故障排除

### 常见问题

#### 1. 内存不足

```
Error: Cannot allocate memory
```

**解决方案**: 减少并行任务数或增加内存限制

```groovy
process {
    withLabel: 'mrna_align' {
        memory = 24.GB
    }
}
```

#### 2. Java版本问题

```
Error: Unable to find Java version 11 or higher
```

**解决方案**: 确保安装 Java 11+

```bash
java -version
# OpenJDK 11.0.11 2021-04-20
```

#### 3. 模块找不到

```
Error: No such process: 'MODULE_NAME'
```

**解决方案**: 确保所有模块文件在 `modules/` 目录中

### 获取帮助

```bash
# 查看帮助信息
nextflow run main.nf --help

# 查看 Nextflow 帮助
nextflow -h

# 查看完整文档
nextflow run main.nf -with-help
```

## 📚 相关资源

- [Nextflow 文档](https://www.nextflow.io/docs/latest/)
- [NJU-seq 主页](https://github.com/IvanWoo22/NJU_seq)
- [bowtie2 文档](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools 文档](http://www.htslib.org/doc/samtools.html)

## 📄 许可证

本项目遵循 MIT 许可证。详见 [LICENSE](../LICENSE) 文件。

## 🤝 贡献

欢迎提交 Issue 和 Pull Request！

## 📞 联系方式

- GitHub: [IvanWoo22/NJU_seq](https://github.com/IvanWoo22/NJU_seq)
- 邮箱: contact@nju-seq.org

---

<div align="center">

**NJU-seq Pipeline** - 从测序数据中识别2'-O-甲基化位点的强大工具

*Built with Nextflow*

</div>
