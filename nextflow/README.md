# Nextflow and cluster execution

**English** | [简体中文](README.zh-CN.md)

The Nextflow scripts in this directory do not support the current Rust program. To use `nju-seq`, install the environment described in the [main README](../README.md) and run the analysis with `nju-seq plan/run`.

A cluster job script can invoke the Rust commands with a sample manifest, references, thread count, and output directory. See [Helper scripts](../docs/scripts.md) for chunked alignment and resource options.
