# 历史文件清理说明

原日期测试目录 `test_0924/`、`test_1001/`、`test_rust_0926/` 已清理。源码和可复用脚本迁至仓库根目录，完整参考迁至 `references/`；原始 FASTQ、原 `data/` 和样本说明保留。

删除的中间数据、重复测试包和大体积结果合计 **69.46 GiB 文件逻辑大小**，不等同于实际释放的磁盘空间。详情见 [清理汇总](cleanup.json)。[原 data 校验值](original_data.sha256.json) 覆盖 29 个原参考与样本说明文件，清理前后内容一致。

Git 中保留当前源码、使用文档、小份测试数据及 [验证汇总](../validation.md)。完整参考、原始 FASTQ、构建目录和运行结果由 `.gitignore` 排除。

旧源码文档压缩包 `dated_workflows_source.tar.gz`、逐文件清单 `removed_files.tsv.gz`、原链接记录 `original_raw_links.json`、详细历史汇总 `evidence/` 与 `evidence_manifest.json` 仍在原工作区本地保存，不随 Git 分发，也不是当前分析的依赖。它们可能引用已删除的结果路径。
