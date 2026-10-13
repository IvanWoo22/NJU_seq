# 旧 Nextflow 流程

此目录保留早期 Nextflow 编排源码，供查阅原分析步骤。配置、脚本参数和表格式尚未适配当前 Rust 程序，不作为当前的运行入口。

新分析请按 [主 README](../README.md) 安装环境并运行 `nju-seq plan/run`。旧步骤与当前命令的对应关系见 [迁移说明](../docs/migration.md)。

如果需要集群调度，可由调度脚本调用 Rust 命令，并显式指定样本表、参考、线程数和输出目录。
