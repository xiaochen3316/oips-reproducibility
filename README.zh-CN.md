# OIPS 可重复性软件包

[English](README.md) | [复现指南](REPRODUCE.md)

本仓库是 OIPS 研究的**公开代码与处理后数据复现仓库**，提供以下论文配套的处理后数据、确定性分析代码、计算公式和冻结参考结果：

> *Choosing among Competing Pockets in Oligomeric Proteins: An
> OIPS-Assisted, Traceable Multi-Evidence Analysis*

该软件包可重建 cluster-v2 口袋候选、计算静态 OIPS-P 评分、映射后验
证据、计算评价指标，并重新生成公开的图表和报告。仓库的目标是让论文中的
处理后分析可检查、可追踪、可重复。

## 快速开始

请在仓库根目录使用干净的 Python 3.11 或 3.12 环境：

```console
python -m pip install --editable ".[test]"
python -m oips_repro reproduce --config config/manuscript.yaml --output results/reproduced
python -m oips_repro verify --config config/manuscript.yaml --bundle results/reproduced --snapshot tests/scientific/data/expected_summary.json
```

第二条命令是“一键复现”入口：先验证公开输入，再生成新的结果包。第三条命令
将结果包与冻结的科学快照逐项核对。环境设置、分步命令、预期文件和安全覆盖
规则见 [REPRODUCE.md](REPRODUCE.md)。

## 冻结的核心结果

公开快照包含 21 个靶点和 1,742 条来源工具记录。其中 1,564 条记录可映射，
178 条记录保留在排除审计表中；同工具合并后得到 1,417 个单元，跨工具聚类
得到 733 个 cluster-v2 候选。冻结结果中的最大簇直径为
11.914302329553335 A，低于含边界的 12 A 上限；任一簇内同一工具最多贡献
一个正式票数。

在 21 个靶点的参考配体评价中，冻结静态排序的 Top-1 为 0.5714285714，
Top-3 为 0.8571428571，Top-5 为 0.9047619048，MRR 为 0.7230158730。
机器可读的预期值和容差位于 `tests/scientific/data/expected_summary.json`；
`results/reference/` 中的文件是只读的规范参考输出。

## 静态排序与后验证据的边界

OIPS-P_static 只使用公开的静态口袋特征、跨工具一致性、几何、可配体性、
证据质量和寡聚体相关性。参考配体、分子动力学摘要和重对接摘要不参与候选
构建，也不进入静态评分。它们只在静态排序冻结后映射，用于审计、解释和评价。

精确公式、阈值包含关系、缺失值规则与排序约定见
[docs/oips-formula.md](docs/oips-formula.md)，方法流程见
[docs/methods.md](docs/methods.md)。

## 仓库结构

- `src/oips_repro/`：确定性实现和命令行入口。
- `config/`：冻结的论文配置、模式和作图契约。
- `data/`：经审核的公开输入、校验和、权利审计和外部归档清单。
- `results/reference/`：不可由命令行覆盖的规范参考表、作图源数据和报告。
- `results/reproduced/`：用户重新生成的输出，不作为参考真值。
- `figures/manuscript/`：仓库总结图；文件名不预设最终稿件图号，正式对应关系须在投稿前于 `FIGURE_SOURCE_DATA_INDEX.tsv` 中确认。
- `docs/`：方法、公式、数据字典、来源链、局限性和期刊声明。
- `tests/`：单元、集成、安全与冻结科学快照测试。

分析结果包还包含五个模块分别进行 ±20% 调整所得的 10 个权重扰动情景，以及
14 个五工具完整案例的原生单工具比较。对应的情景级、体系级和方法级 CSV 位于
`results/reference/analysis/`，均由标准复现流程重新生成。

论文审阅时可使用 [`MANUSCRIPT_RESULT_TRACEABILITY.tsv`](MANUSCRIPT_RESULT_TRACEABILITY.tsv)：
它将论文图 1–7、表 1–3 对应到仓库中的精确源表，并明确标出仍由作者管理的
结构渲染图或需要未公开轨迹的面板。

## 复现范围与限制

本仓库复现的是已经公开的处理后分析，**不**从头重新运行专有网络服务、商业
对接软件或完整分子动力学轨迹。原始服务包、商业软件二进制文件和日志、轨迹、
重启文件及检查点均未纳入 Git。其状态和再分发边界记录在
[THIRD_PARTY.md](THIRD_PARTY.md)、`data/metadata/asset-rights.tsv` 和
`data/external_archive_manifest.tsv` 中。

公开结构和标准化的服务来源测量仍受各自来源条款与署名要求约束。尚未核实的
上游版本、访问日期和归档标识符不会被推测，而是保留为明确的发布前检查项。

## 引用与发布状态

软件作者元数据见 [CITATION.cff](CITATION.cff)。当前为发布前软件包，论文、
代码归档、数据归档和 GitHub 仓库尚无已核实的公开标识符，因此本文不填入
DOI 或仓库网址。完成预留并核验落地页后，再按
[PRE_PUBLICATION_CHECKLIST.md](PRE_PUBLICATION_CHECKLIST.md) 补充。

计划版本为 `1.0.0`，拟使用标签 `v1.0.0-manuscript`。在发布清单完成、且
对应提交被归档之前，引用时请同时提供 `CITATION.cff` 中的软件作者信息和
材料所对应的具体 Git 提交。

## 许可证

代码采用 [BSD 3-Clause](LICENSE) 许可证。团队原创的文档、表格数据和图件
采用 [CC BY 4.0](LICENSES/CC-BY-4.0.txt)，但文件或记录明确标注第三方
条款的内容除外。本仓库不对第三方材料重新许可，详见
[THIRD_PARTY.md](THIRD_PARTY.md)。
