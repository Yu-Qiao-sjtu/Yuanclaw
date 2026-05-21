# YuanSeq

<div align="center">

**An Integrated R/Shiny Platform for RNA-seq, Microarray & Single-Cell Downstream Analysis**

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%3E4.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Shiny-green.svg)](https://shiny.posit.co/)

**Developer:** [Yu Qiao (乔宇)](https://github.com/Passpoor) · PhD in Pharmacology · SJTU School of Pharmacy

**Supervisors:** [Prof. Feng Qian](https://pharm.sjtu.edu.cn/szdy/2862.html), [Prof. Lei Sun](https://pharm.sjtu.edu.cn/szdy/2870.html)

*Part of the **Yuanclaw** project at SJTU School of Pharmacy*

</div>

---

## 📖 Overview

**YuanSeq** is a modular, all-in-one bioinformatics analysis platform for comprehensive downstream analysis of bulk RNA-seq, microarray, and **single-cell cluster-derived pseudobulk data**.

### Key Features

| Module | Description | Methods |
|--------|-------------|---------|
| 🔬 **Differential Expression** | RNA-seq / microarray DE analysis | limma-voom, edgeR; 1v1 / nvn comparisons |
| 🧬 **Functional Enrichment** | KEGG, GO, GSEA analysis | clusterProfiler, local KEGG (no rate limit) |
| 🛤️ **Pathway Activity** | Pathway activity inference | ULM, WMEAN, AUCell, GSVA (decoupleR) |
| 🔬 **TF Activity v2.0** | Transcription factor activity (DoRothEA confidence levels + TF family enrichment) | CollecTRI, DoRothEA (A-D) + decoupleR; Fisher exact test |
| 🆕 **Single-Cell Support** | Downstream analysis for scRNA-seq cluster markers | Use cluster DEGs for enrichment & activity inference |
| 🤖 **AI Interpretation** | AI-powered biological interpretation | DeepSeek, OpenAI, Zhipu AI, Local models |
| 📊 **Visualization** | Interactive plots with sci-fi themed UI | ggplot2, plotly, pheatmap |

---


## 🗂️ Repository Layout

To reduce path confusion and keep runtime/debug flows predictable, the project now follows a documented layout baseline.

- Runtime app entrypoints: `app.R` and `inst/shiny/app.R`
- Runtime modules: `modules/` and `inst/shiny/modules/` (kept in sync)
- Operational scripts: `scripts/`
- Automated tests: `tests/`
- Historical materials: `archive/`

See detailed structure and cleanup conventions in [`docs/REPO_STRUCTURE.md`](docs/REPO_STRUCTURE.md).

---

## 🚀 Installation & Launch

### Quick Start

| Operation | Code |
|-----------|------|
| **Install** | `remotes::install_github("Passpoor/Yuanclaw", upgrade = "never")` |
| **Update** | `remotes::install_github("Passpoor/Yuanclaw", upgrade = "never", force = TRUE)` |
| **Launch** | `YuanSeq::run_app()` |

**First-time installation:**
```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Passpoor/Yuanclaw", upgrade = "never")
YuanSeq::run_app()
```

**Already installed? Just launch:**
```r
YuanSeq::run_app()
```

### Configure AI API (Optional)

```r
# Copy config template
config_dir <- file.path(Sys.getenv("USERPROFILE"), ".yuanseq")  # Windows
# config_dir <- "~/.yuanseq"  # Mac/Linux
dir.create(config_dir, showWarnings = FALSE)
file.copy(system.file("api_config.example.json", package = "YuanSeq"),
          file.path(config_dir, "api_config.json"))
```

Edit `api_config.json` with your API key:

```json
{
  "provider": "deepseek",
  "api_key": "YOUR_API_KEY",
  "model": "deepseek-chat"
}
```

**Supported AI Providers:**

| Provider | `provider` value | Recommended Model |
|----------|------------------|-------------------|
| [DeepSeek](https://www.deepseek.com/) | `deepseek` | deepseek-chat |
| [OpenAI](https://openai.com/) | `openai` | gpt-4o |
| [Zhipu AI](https://www.bigmodel.cn/) | `zhipu` | glm-4-flash |
| [Local Models](https://ollama.com/) | `local` | custom |

> ⚠️ **Data Security**: External API calls send data to third-party servers. For sensitive data, use local models.

---

### 方式一：本地安装

```bash
git clone https://github.com/Passpoor/Yuanclaw.git
cd Yuanclaw
```

```r
source("install.R")  # 安装依赖包
source("run.R")      # 启动应用
```

### 方式二：Pixi 管理（推荐，可复现环境）

[Pixi](https://pixi.sh/) 是一款跨平台的包管理工具，基于 conda-forge 和 bioconda 提供锁定版本的环境管理。

**前置条件：** [安装 pixi](https://pixi.sh/latest/#installation)

```bash
git clone https://github.com/Passpoor/Yuanclaw.git
cd Yuanclaw

# 安装所有依赖
pixi install

# 首次设置：配置 R 环境并安装额外包
pixi run setup

# 启动应用
pixi run app
```

**可用命令：**

| 命令 | 说明 |
|------|------|
| `pixi install` | 安装 conda/bioconda 依赖 |
| `pixi run setup` | 首次设置（配置 R 环境 + 安装额外包） |
| `pixi run app` | 启动 Shiny 应用（端口 3838） |
| `pixi run test` | 运行测试 |
| `pixi run check-r` | 查看 R 版本 |

> **提示：** `config/Rprofile.site` 会自动将系统已安装的 R 包加入 pixi 环境的搜索路径，因此 GO.db、org.Hs.eg.db、org.Mm.eg.db 等 bioconductor 注释数据库无需重复下载。

---

## 🤖 AI Interpretation

YuanSeq integrates AI-powered biological interpretation:

- **Context-aware analysis** based on sample metadata (organism, tissue, treatment)
- **Multi-format export**: Markdown, HTML (with embedded figures), PDF
- **Real-time progress indicator**
- **Conversation history tracking**

---

## 🆕 Single-Cell Downstream Analysis

YuanSeq supports downstream analysis of **single-cell RNA-seq data**:

- Import cluster marker genes from Seurat, Scanpy, or other scRNA-seq tools
- Perform KEGG/GO/GSEA enrichment on cluster-specific DEGs
- Infer pathway and TF activity for each cell cluster
- Compare biological programs across cell subpopulations

---


---

## 🔬 TF Activity v2.0

The TF activity module has been upgraded with major new features:

### 🆕 DoRothEA Database Support

| Feature | CollecTRI | DoRothEA |
|---------|-----------|----------|
| Coverage | ~2,500 TFs | ~1,500 TFs |
| Confidence levels | None | A/B/C/D four-level grading |
| Best for | Comprehensive exploration | High-confidence validation |

**Confidence Levels:**
- **A** (Highest): Orthogonal experimental validation (ChIP + KO)
- **B** (High): Multiple high-throughput experiments
- **C** (Medium): Expert-reviewed computational predictions
- **D** (Low): Single-source evidence

> 💡 **Recommendation**: Use A+B for publication-grade analysis; A-D for exploratory analysis.

### 🆕 TF Family Enrichment Analysis

Based on **Lambert et al. 2018 (Cell)**, the module now supports:

- **30+ TF families** (bZIP, bHLH, Homeobox, Nuclear Receptor, ETS, Forkhead, C2H2_ZF, etc.)
- **Fisher exact test** for family enrichment significance
- **BH correction** for multiple testing
- **Rich visualizations**: Enrichment bar chart, bubble chart, Lollipop plot

### Update Instructions

For users who already installed YuanSeq:

```r
# Simply reinstall to get the update
remotes::install_github("Passpoor/Yuanclaw", upgrade = "never", force = TRUE)

# Then launch as usual
YuanSeq::run_app()
```

**References:**
1. Lambert et al. (2018) The Human Transcription Factors. *Cell* 172(4):650-665.
2. Garcia-Alonso et al. (2019) Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Biology* 20:240.
3. Badia-i-Mompel et al. (2022) decoupleR. *Bioinformatics* 38(22):5174-5176.

## 🙏 Acknowledgments

YuanSeq integrates and acknowledges the following open-source projects:

| Category | Packages | Use |
|----------|----------|-----|
| **Framework** | [shiny](https://cran.r-project.org/package=shiny), [shinyjs](https://cran.r-project.org/package=shinyjs), [bslib](https://cran.r-project.org/package=bslib), [DT](https://cran.r-project.org/package=DT), [plotly](https://cran.r-project.org/package=plotly) | App framework & UI |
| **DE Analysis** | [edgeR](https://bioconductor.org/packages/edgeR/), [limma](https://bioconductor.org/packages/limma/) | Differential expression |
| **Enrichment** | [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/), [enrichplot](https://bioconductor.org/packages/enrichplot/), [GseaVis](https://bioconductor.org/packages/GseaVis/) | GO/KEGG/GSEA |
| **Local KEGG** | [biofree.qyKEGGtools](https://github.com/Passpoor/biofree.qyKEGGtools) | Offline KEGG enrichment |
| **Activity** | [decoupleR](https://bioconductor.org/packages/decoupleR/) | Pathway/TF activity |
| **Visualization** | [ggplot2](https://cran.r-project.org/package=ggplot2), [pheatmap](https://cran.r-project.org/package=pheatmap), [ggrepel](https://cran.r-project.org/package=ggrepel) | Plotting |
| **AI** | [httr](https://cran.r-project.org/package=httr), [jsonlite](https://cran.r-project.org/package=jsonlite) | API calls |

Thanks to R, Bioconductor, and all package developers!

---

## 📁 Project Structure

```
Yuanclaw/
├── app.R                      # Main entry
├── api_config.example.json    # API config template
├── modules/
│   ├── ui_theme.R             # Theme & layout
│   ├── data_input.R           # Data upload
│   ├── differential_analysis.R
│   ├── kegg_enrichment.R
│   ├── gsea_analysis.R
│   ├── pathway_activity.R
│   ├── tf_activity.R
│   ├── ai_interpretation.R    # AI interpretation
│   └── venn_diagram.R
├── R/                         # Package functions
├── inst/shiny/                # Packaged Shiny app
└── config/                    # Configuration
```

---

## 📄 License

[MIT License](LICENSE)
