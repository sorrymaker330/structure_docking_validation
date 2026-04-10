# 🔬 structure_docking_validation (scDock)

[![Framework: CellForge-Agent](https://img.shields.io/badge/Framework-CellForge--Agent-blue.svg)](https://github.com/xiaolab-xjtu/CellForge-Agent)
[![Method: AutoDock Vina](https://img.shields.io/badge/Algorithm-AutoDock%20Vina-success.svg)](#)

> **基于单细胞通讯与受体靶点的分子对接虚拟筛选技能 (Agent Skill)**

## 🧬 简介 (Overview)

`structure_docking_validation` 是为 **[CellForge Agent](https://github.com/xiaolab-xjtu/CellForge-Agent)** 框架开发的自动化 3D 结构药物验证技能。

在单细胞分析中，当我们通过细胞通讯分析（Cell-Cell Communication）或其他方法锁定了一个高度异常的**单一受体蛋白**时，本技能可以通过结构生物学（Structural Biology）手段，在微观物理层面上验证潜在药物与该受体的结合亲和力。它无缝集成了 AlphaFold 蛋白结构预测与 AutoDock Vina 对接引擎，实现了完全自动化的 *in silico* 靶向筛选。

## 🔗 在 CellForge 流水线中的定位

本技能属于 CellForge AI 药物发现管线中 **自下而上 (Bottom-up)** 的微观结构验证模块。
它通常作为 `scanpy_drug_prediction` (基于转录组特征逆转的宏观药筛) 技能的下游步骤。当 Agent 在转录组层面筛选出具有逆转潜力的药物名单后，可立即调用本技能，导入候选药物的 CAS 号与受体蛋白的 UniProt ID，进行 3D 物理构象的精确验证，从而完成**“转录组网络纠偏 + 物理结构验证”**的闭环。

## ⚙️ 核心机制 (Core Mechanism)

当大模型 Agent 触发此技能时，底层 Python 脚本会自动执行以下流水线：
1. **靶点准备**：利用 `download_alphafold.py` 通过输入的 `uniprot_id` 自动从 AlphaFold DB 下载高精度的靶点 3D 结构 (PDB格式)，并由 `prepare_receptor.py` 转化为对接所需的 PDBQT 格式。
2. **配体准备**：利用 `download_cas_pubchem.py` 根据药物的 `cas_number` 从 PubChem 抓取 SMILES 表达式，并使用 RDKit 和 Meina 生成 3D 构象的配体 PDBQT 文件。
3. **构象对接**：调用 AutoDock Vina 引擎，在指定的活性口袋（Grid Box）内进行构象搜索与热力学打分。
4. **结果提取**：解析 Vina 输出，将最强结合亲和力（Binding Affinity, kcal/mol）作为关键指标返回给大模型 Critic 层。

## 🛠️ 参数说明 (Parameters)

大模型 Agent 可自主推断或由用户配置以下参数：

* `uniprot_id`: 目标受体蛋白的 UniProt 唯一标识符（例如 `'P00533'` 代表 EGFR）。
* `cas_number`: 候选小分子药物的 CAS 注册号（例如 `'141430-65-1'` 代表 Imatinib）。
* `grid_center`: (可选) 对接中心坐标，格式为 `[x, y, z]`。若不提供，将尝试盲对接。
* `grid_size`: (可选) 对接盒子大小，格式为 `[x, y, z]`。默认为 `[20, 20, 20]`。

## 📊 输出与评估 (Outputs & Critic)

* **执行输出**：生成靶点与配体的 PDBQT 文件，以及 Vina 的完整对接日志与构象文件。
* **认知与反思 (Critic Layer)**：技能内部设有统计学防火墙，默认提取 `binding_affinity`。当亲和力 $\le -7.0$ kcal/mol 时判定为有效结合；若大于 $-5.0$ kcal/mol 则触发警告，提示大模型该配体与靶点结合较弱。

---
*Developed for CellForge Agent Pipeline.*