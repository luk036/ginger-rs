# 🫚 ginger-rs - IFLOW Context

## 项目概述

`ginger-rs` 是一个用 Rust 编写的并行多项式求根算法库，主要实现了以下几种方法：

- **Aberth 方法**：用于寻找多项式的所有根
- **Bairstow 方法**：专门用于寻找实系数多项式的二次因子（从而获得复根）
- **Horner 方法**：用于多项式求值
- **Leja 排序**：数值稳定性优化
- **2x2 矩阵和 2D 向量**：底层数值计算支持

该库支持并行计算（使用 Rayon），通过多线程来加速求根过程。

## 构建和运行

### 依赖项

- Rust 工具链 (edition 2021)
- `rayon`：用于并行计算
- `num-complex`：用于复数运算
- `num-traits`：数值特征
- `lds-rs`：低差异序列（用于初始猜测）

### 构建命令

```bash
# 构建项目
cargo build

# 构建发布版本
cargo build --release

# 运行测试
cargo test

# 运行基准测试
cargo bench
```

### 安装命令

```bash
cargo install ginger-rs
```

## 主要模块和功能

### 核心 API

- `rootfinding` 模块：实现 Bairstow 方法（`pbairstow_even`, `pbairstow_even_mt`, `pbairstow_autocorr`, `pbairstow_autocorr_mt`）
- `aberth` 模块：实现 Aberth 方法（`aberth`, `aberth_mt`）
- `horner` 模块：实现 Horner 评估（`horner_eval_f`, `horner_eval_c`）
- `Options` 结构：控制算法参数（`max_iters`, `tolerance`, `tol_ind`）

### 关键类型

- `Vector2<f64>`：表示 2D 向量，用于 Bairstow 方法中的系数对
- `Matrix2<f64>`：表示 2x2 矩阵，用于数值计算
- `Complex<f64>`：用于 Aberth 方法

## 开发约定

### 代码风格

- 使用 Rust 标准编码规范
- 使用 `Vector2` 和 `Matrix2` 自定义类型进行数值计算
- 使用 `num-complex` 库处理复数

### 算法实现

- 实现了多种多项式求根算法，支持串行和并行版本
- 使用 `Options` 结构统一控制算法参数
- 包含专门针对自相关多项式的优化算法

### 测试

- 包含全面的单元测试
- 使用 `approx_eq` 进行浮点数比较
- 测试用例覆盖各种算法和边界情况

## 并行化

- 使用 Rayon 库实现多线程并行计算
- 提供多线程版本的 `pbairstow_even_mt` 和 `aberth_mt` 函数
- 通过并行处理同时更新多个根的估计值

## 应用场景

此库适用于需要高效多项式求根的科学计算、信号处理和控制系统等领域，尤其是当需要处理高次多项式或需要并行计算能力时。