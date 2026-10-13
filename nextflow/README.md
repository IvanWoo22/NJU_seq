# Nextflow 与集群运行

此目录的 Nextflow 脚本不支持当前 Rust 程序。使用 `nju-seq` 时，按 [主 README](../README.md) 安装环境，通过 `nju-seq plan/run` 运行分析。

集群任务可由调度脚本调用 Rust 命令，并指定样本清单、参考、线程数和输出目录。分批比对及资源参数见 [辅助脚本](../docs/scripts.md)。
