"""Packaged for revise/; source: project/01this_study_dataset/sgb/sig_mag/sgb_da_and_plots.py
Figures: Fig4B
Note: MAG DA pipeline
"""
from pathlib import Path
import os, sys
PKG_ROOT = Path(__file__).resolve().parents[1]
SHARED = PKG_ROOT / 'shared_input'
OUTDIR = PKG_ROOT / 'output' / Path(__file__).stem
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from paths_config import PATHS, SHARED as _SHARED, PKG_ROOT as _PKG
    SHARED, PKG_ROOT = _SHARED, _PKG
except Exception:
    pass
# --- original script below; prefer SHARED/PATHS over hard-coded host paths ---

#!<USER_HOME>/software/miniconda3/bin/python
"""
SGB 差异丰度分析与多层级对比可视化

功能概述
--------
1. 使用与 MetaPhlAn 种水平完全相同的差异分析流程（复用 v2._compute_da），
   对 SGB 丰度表 `bin_abundance_table.tsv` 进行差异分析；
2. 读取 MetaPhlAn 种水平差异结果表 `Table_S_differential_species.csv`；
3. 通过提供的 SGB–物种映射文件，将 SGB 与 MetaPhlAn 物种对齐，生成双层级对比结果：
   - 两个层级的 log2FC、q 值（BH FDR）、方向（Normal / PoorPlaque）；
   - 是否双层级均显著（q < 0.05）且方向一致；
4. 输出：
   - SGB 差异分析完整结果表（含统计量、方向、q 值）；
   - 双层级对比汇总表；
   - Group A / Group B 富集 SGB 列表；
   - CNS 级三类图形（PDF + PNG，300 dpi）。

依赖
----
- Python 3.10+
- numpy, pandas, matplotlib, seaborn, statsmodels
- 本目录下的 `v2.py`（提供 `_compute_da`、`load_data`、`save_fig` 等函数）

输入文件说明
------------
- bin_abundance_table.tsv : SGB 丰度表，行为 SGB（或 MAG ID），列为样本
- this.study.group.tsv    : 两列 [sample, group]，group 与 v2.py 一致（Normal / PoorPlaque）
- Table_S_differential_species.csv : MetaPhlAn 种水平差异结果（由 v2.py 导出）
- sgb_map.tsv             : SGB 到 MetaPhlAn 物种名映射表（用户提供）
                            必须包含列：
                              - SGB      : 与 bin_abundance_table.tsv 的行名一致
                              - species  : 与 Table_S_differential_species.csv 的 species 列一致
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import v2  # 复用差异分析与绘图风格


def parse_args() -> argparse.Namespace:
    """解析命令行参数。

    Returns
    -------
    argparse.Namespace
        参数对象。
    """
    parser = argparse.ArgumentParser(
        description="SGB 差异丰度分析 + 与 MetaPhlAn 种水平结果对比 + CNS 级可视化",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--abund",
        "-a",
        type=Path,
        required=True,
        help="SGB 丰度表，制表符分隔，行为 SGB/MAG，列为样本（如 bin_abundance_table.tsv）",
    )
    parser.add_argument(
        "--group",
        "-g",
        type=Path,
        required=True,
        help="分组文件，两列 [sample, group]，与 v2.py 一致（如 this.study.group.tsv）",
    )
    parser.add_argument(
        "--species-da",
        "-s",
        type=Path,
        required=True,
        help="MetaPhlAn 种水平差异结果表（Table_S_differential_species.csv）",
    )
    parser.add_argument(
        "--sgb-map",
        "-m",
        type=Path,
        required=True,
        help="SGB 到 MetaPhlAn 物种名映射表（制表符分隔，至少包含列 SGB 和 species）",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        required=True,
        help="输出目录，将在其中创建 tables/ 与 figures/ 子目录",
    )
    parser.add_argument(
        "--group-a",
        type=str,
        default="Normal",
        help="定义 Group A 名称（用于富集方向拆分，默认为 Normal）",
    )
    parser.add_argument(
        "--group-b",
        type=str,
        default="PoorPlaque",
        help="定义 Group B 名称（用于富集方向拆分，默认为 PoorPlaque）",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="图像分辨率（DPI）",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="随机种子",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="输出调试级别日志",
    )
    return parser.parse_args()


def setup_logging(out_dir: Path, verbose: bool) -> None:
    """配置日志输出到终端和文件。

    Parameters
    ----------
    out_dir : Path
        输出目录，将在其中创建日志文件。
    verbose : bool
        是否使用 DEBUG 级别日志。
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "sgb_da_and_plots.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_path, mode="w"),
        ],
    )


def run_sgb_da(
    abund_path: Path,
    group_path: Path,
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, List[str], List[str]]:
    """对 SGB 丰度表执行与 v2.py 相同的差异分析流程。

    直接复用 v2.load_data 与 v2._compute_da，保证筛选条件、统计方法、
    多重校正方式与 MetaPhlAn 种水平分析完全一致。

    Parameters
    ----------
    abund_path : Path
        SGB 丰度表路径。
    group_path : Path
        分组文件路径。
    seed : int
        随机种子。

    Returns
    -------
    da_result : pd.DataFrame
        SGB 差异分析完整结果（含 log2FC、p 值、FDR、Robust_DA 等）。
    abund : pd.DataFrame
        过滤后用于差异分析的丰度矩阵。
    sn : list of str
        Group A（Normal）样本列表。
    sp : list of str
        Group B（PoorPlaque）样本列表。
    """
    logging.info("[SGB-DA] 读取丰度表与分组信息 ...")
    abund, groups, sn, sp = v2.load_data(str(abund_path), str(group_path))
    logging.info(
        "[SGB-DA] matrix=%d x %d | %s n=%d, %s n=%d",
        abund.shape[0],
        abund.shape[1],
        "Normal",
        len(sn),
        "PoorPlaque",
        len(sp),
    )

    logging.info("[SGB-DA] 运行 _compute_da（与 MetaPhlAn 种水平一致） ...")
    da_result, l2fc, fdr_wil, robust, robust_up, robust_dn, species = v2._compute_da(
        abund, groups, sn, sp, seed
    )

    # 为下游兼容性，显式添加列名别名
    da_result = da_result.copy()
    da_result.rename(columns={"species": "SGB"}, inplace=True)

    # 补充方向列（基于 log2FC，而不是 Robust_direction）
    da_result["Direction"] = np.where(
        da_result["log2FC"] > 0,
        "PoorPlaque",
        np.where(da_result["log2FC"] < 0, "Normal", "NS"),
    )

    return da_result, abund.loc[species], sn, sp


def load_species_da(species_da_path: Path) -> pd.DataFrame:
    """读取 MetaPhlAn 种水平差异结果表。

    Parameters
    ----------
    species_da_path : Path
        种水平差异结果 CSV 文件路径。

    Returns
    -------
    pd.DataFrame
        差异结果数据框，附加方向列。
    """
    logging.info("[META-DA] 读取种水平差异结果: %s", species_da_path)
    df = pd.read_csv(species_da_path)
    if "species" not in df.columns:
        raise ValueError("种水平差异结果表必须包含列 'species'")

    df = df.copy()
    if "log2FC" not in df.columns or "Wilcoxon_FDR" not in df.columns:
        raise ValueError("种水平差异结果表必须包含列 'log2FC' 与 'Wilcoxon_FDR'")

    df["Direction"] = np.where(
        df["log2FC"] > 0,
        "PoorPlaque",
        np.where(df["log2FC"] < 0, "Normal", "NS"),
    )
    return df


def load_sgb_map(sgb_map_path: Path) -> pd.DataFrame:
    """读取 SGB–物种映射表。

    Parameters
    ----------
    sgb_map_path : Path
        映射表路径，要求为制表符分隔，至少包含列 'SGB' 与 'species'。

    Returns
    -------
    pd.DataFrame
        映射表，包含列 ['SGB', 'species']。
    """
    logging.info("[MAP] 读取 SGB–species 映射表: %s", sgb_map_path)
    df = pd.read_csv(sgb_map_path, sep="\t")
    missing = {c for c in ("SGB", "species") if c not in df.columns}
    if missing:
        raise ValueError(f"SGB 映射表缺少列: {missing}")
    return df[["SGB", "species"]].drop_duplicates()


def build_cross_level_summary(
    sgb_da: pd.DataFrame,
    species_da: pd.DataFrame,
    sgb_map: pd.DataFrame,
    q_cutoff: float = 0.05,
) -> pd.DataFrame:
    """构建 SGB 与 MetaPhlAn 种水平的双层级对比结果表。

    Parameters
    ----------
    sgb_da : pd.DataFrame
        SGB 差异分析结果（含列 SGB, log2FC, Wilcoxon_FDR, Direction）。
    species_da : pd.DataFrame
        种水平差异分析结果（含列 species, log2FC, Wilcoxon_FDR, Direction）。
    sgb_map : pd.DataFrame
        SGB–物种映射表（列 SGB, species）。
    q_cutoff : float, default 0.05
        显著性阈值（BH FDR）。

    Returns
    -------
    pd.DataFrame
        双层级对比汇总表，每行对应一个 SGB。
    """
    logging.info("[SUMMARY] 构建双层级对比表 ...")

    sgb_da_sel = sgb_da[
        ["SGB", "log2FC", "Wilcoxon_FDR", "Direction"]
    ].rename(
        columns={
            "log2FC": "log2FC_SGB",
            "Wilcoxon_FDR": "q_SGB",
            "Direction": "Direction_SGB",
        }
    )

    species_da_sel = species_da[
        ["species", "log2FC", "Wilcoxon_FDR", "Direction"]
    ].rename(
        columns={
            "log2FC": "log2FC_Meta",
            "Wilcoxon_FDR": "q_Meta",
            "Direction": "Direction_Meta",
        }
    )

    merged = sgb_map.merge(sgb_da_sel, on="SGB", how="left").merge(
        species_da_sel, on="species", how="left"
    )

    # 显著性与方向一致性
    merged["sig_SGB"] = merged["q_SGB"] < q_cutoff
    merged["sig_Meta"] = merged["q_Meta"] < q_cutoff
    merged["both_sig"] = merged["sig_SGB"] & merged["sig_Meta"]
    merged["direction_consistent"] = (
        merged["Direction_SGB"].eq(merged["Direction_Meta"])
    ) & merged["Direction_SGB"].isin(["Normal", "PoorPlaque"])

    def classify(row: pd.Series) -> str:
        """根据显著性与方向关系进行分类标签。"""
        if pd.isna(row["q_SGB"]) or pd.isna(row["q_Meta"]):
            return "missing"
        if row["both_sig"] and row["direction_consistent"]:
            return "both_sig_same_dir"
        if row["both_sig"] and not row["direction_consistent"]:
            return "both_sig_opposite_dir"
        if row["sig_SGB"] and not row["sig_Meta"]:
            return "sgb_only_sig"
        if (not row["sig_SGB"]) and row["sig_Meta"]:
            return "meta_only_sig"
        return "both_nonsig"

    merged["consistency_category"] = merged.apply(classify, axis=1)
    logging.info(
        "[SUMMARY] 分类计数: %s",
        merged["consistency_category"].value_counts().to_dict(),
    )

    return merged


def export_enriched_lists(
    summary: pd.DataFrame,
    group_a: str,
    group_b: str,
    out_dir: Path,
) -> Tuple[Path, Path]:
    """基于双层级显著且方向一致的 SGB，导出 Group A / Group B 富集列表。

    Parameters
    ----------
    summary : pd.DataFrame
        双层级对比汇总表。
    group_a : str
        Group A 名称（通常为 Normal）。
    group_b : str
        Group B 名称（通常为 PoorPlaque）。
    out_dir : Path
        输出目录。

    Returns
    -------
    Tuple[Path, Path]
        (Group A 富集列表路径, Group B 富集列表路径)，
        列表中同时包含 SGB ID 与物种名。
    """
    mask_core = summary["consistency_category"].eq("both_sig_same_dir")
    core = summary.loc[mask_core].copy()

    cols_for_export = [
        "SGB",
        "species",
        "Direction_SGB",
        "log2FC_SGB",
        "q_SGB",
        "log2FC_Meta",
        "q_Meta",
        "Direction_Meta",
    ]

    a_df = (
        core.loc[core["Direction_SGB"].eq(group_a), cols_for_export]
        .dropna(subset=["SGB"])
        .drop_duplicates()
    )
    b_df = (
        core.loc[core["Direction_SGB"].eq(group_b), cols_for_export]
        .dropna(subset=["SGB"])
        .drop_duplicates()
    )

    tables_dir = out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    a_path = tables_dir / f"Group_{group_a}_enriched_SGB.txt"
    b_path = tables_dir / f"Group_{group_b}_enriched_SGB.txt"

    a_df.to_csv(a_path, sep="\t", index=False)
    b_df.to_csv(b_path, sep="\t", index=False)

    logging.info(
        "[OUTPUT] %s 富集 SGB: %d; %s 富集 SGB: %d",
        group_a,
        a_df.shape[0],
        group_b,
        b_df.shape[0],
    )
    return a_path, b_path


def plot_cross_level_scatter(
    summary: pd.DataFrame,
    out_dir: Path,
    dpi: int,
) -> None:
    """图 1：SGB vs MetaPhlAn 差异结果散点对比图。

    x 轴为 MetaPhlAn log2FC，y 轴为 SGB log2FC，颜色区分不同一致性类别，
    并高亮“双层级均显著且方向一致”的 SGB。
    """
    logging.info("[FIG1] 绘制双层级差异对比散点图 ...")
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    df = summary.dropna(subset=["log2FC_SGB", "log2FC_Meta"]).copy()

    # 调色：冷暖对比 + 低饱和度背景
    palette = {
        "both_sig_same_dir": "#B85450",  # 暖色，高亮
        "both_sig_opposite_dir": "#D095BF",  # 对照色
        "sgb_only_sig": "#2E8B8B",
        "meta_only_sig": "#92C5C5",
        "both_nonsig": "#DDDDDD",
        "missing": "#BBBBBB",
    }

    plt.close("all")
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(3.5, 3.0))  # 单栏宽度 ~ 89mm

    # 底层：非核心点
    for cat, sub in df.groupby("consistency_category"):
        if cat == "both_sig_same_dir":
            continue
        ax.scatter(
            sub["log2FC_Meta"],
            sub["log2FC_SGB"],
            s=20,
            c=palette.get(cat, "#CCCCCC"),
            alpha=0.6,
            linewidths=0,
            label=cat,
            rasterized=False,
        )

    # 高亮核心 SGB
    core = df[df["consistency_category"].eq("both_sig_same_dir")].copy()
    ax.scatter(
        core["log2FC_Meta"],
        core["log2FC_SGB"],
        s=40,
        c=palette["both_sig_same_dir"],
        alpha=0.9,
        edgecolor="black",
        linewidths=0.4,
        label="both_sig_same_dir",
        rasterized=False,
    )

    # 可选：为核心 SGB 添加标签（数量过多时可手动筛选）
    top_core = core.assign(abs_sgb=lambda x: x["log2FC_SGB"].abs()).sort_values(
        "abs_sgb", ascending=False
    ).head(20)
    for _, row in top_core.iterrows():
        ax.text(
            row["log2FC_Meta"],
            row["log2FC_SGB"],
            row["SGB"],
            fontsize=5.5,
            ha="left",
            va="center",
        )

    ax.axhline(0, color="#AAAAAA", lw=0.7, ls="--")
    ax.axvline(0, color="#AAAAAA", lw=0.7, ls="--")

    ax.set_xlabel(r"MetaPhlAn $\log_2$FC (Group B vs Group A)")
    ax.set_ylabel(r"SGB $\log_2$FC (Group B vs Group A)")

    # 自定义图例标签（中文友好）
    legend_labels = {
        "both_sig_same_dir": "两层级均显著且方向一致",
        "both_sig_opposite_dir": "两层级均显著但方向相反",
        "sgb_only_sig": "仅 SGB 显著",
        "meta_only_sig": "仅 MetaPhlAn 显著",
        "both_nonsig": "两层级均不显著",
        "missing": "缺少一侧结果",
    }
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [legend_labels.get(l, l) for l in labels]
    ax.legend(
        handles,
        new_labels,
        fontsize=6,
        frameon=False,
        loc="best",
        handlelength=1.2,
    )

    sns.despine(ax=ax)
    plt.tight_layout()
    v2.save_fig(fig, str(fig_dir / "Fig1_cross_level_scatter"), dpi=dpi)
    plt.close(fig)


def plot_enriched_waterfall(
    summary: pd.DataFrame,
    group_a: str,
    group_b: str,
    out_dir: Path,
    dpi: int,
) -> None:
    """图 2：筛选后 SGB 富集效果瀑布 / 水平条形图（背靠背）。

    仅使用“双层级均显著且方向一致”的 SGB，按 SGB log2FC 排序，
    左侧为 Group A 富集（负向），右侧为 Group B 富集（正向）。
    """
    logging.info("[FIG2] 绘制筛选后 SGB 富集瀑布图 ...")
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    core = summary[summary["consistency_category"].eq("both_sig_same_dir")].copy()
    if core.empty:
        logging.warning("[FIG2] 无双层级显著且方向一致的 SGB，跳过图 2 绘制")
        return

    core = core.sort_values("log2FC_SGB")
    y_labels = core["SGB"].tolist()
    x_vals = core["log2FC_SGB"].values
    directions = core["Direction_SGB"].values

    colors = np.where(
        np.array(directions) == group_b,
        "#B85450",  # Group B 富集：暖色
        "#2E8B8B",  # Group A 富集：冷色
    )

    plt.close("all")
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(4.0, 4.5))  # 稍接近双栏高度

    y_pos = np.arange(len(y_labels))
    ax.barh(
        y_pos,
        x_vals,
        color=colors,
        alpha=0.9,
        edgecolor="none",
    )
    ax.axvline(0, color="#333333", lw=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=6)
    ax.set_xlabel(r"SGB $\log_2$FC (Group B vs Group A)")
    ax.set_ylabel("SGB")

    # 自定义图例
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor="#2E8B8B", label=f"富集于 {group_a}"),
        Patch(facecolor="#B85450", label=f"富集于 {group_b}"),
    ]
    ax.legend(
        handles=legend_handles,
        fontsize=7,
        frameon=False,
        loc="lower right",
    )

    sns.despine(ax=ax)
    plt.tight_layout()
    v2.save_fig(fig, str(fig_dir / "Fig2_SGB_waterfall"), dpi=dpi)
    plt.close(fig)


def plot_sgb_boxplots(
    abund: pd.DataFrame,
    sn: List[str],
    sp: List[str],
    summary: pd.DataFrame,
    out_dir: Path,
    dpi: int,
) -> None:
    """图 3：筛选通过 SGB 的丰度箱线图（Group A vs Group B）。

    使用 v2._boxplot_on_ax 保持与 Figure 1 一致的风格。
    每个 SGB 单独一个子图，叠加散点与显著性符号。
    """
    logging.info("[FIG3] 绘制 SGB 丰度箱线图 ...")
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    core = summary[summary["consistency_category"].eq("both_sig_same_dir")].copy()
    sgbs = [s for s in core["SGB"].tolist() if s in abund.index]

    if not sgbs:
        logging.warning("[FIG3] 无可用 SGB 绘制箱线图，跳过图 3")
        return

    n = len(sgbs)
    ncols = 4
    nrows = int(np.ceil(n / ncols))

    plt.close("all")
    sns.set_style("ticks")
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(ncols * 2.0, nrows * 2.2), squeeze=False
    )

    for idx, sgb in enumerate(sgbs):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row, col]

        vals_n = abund.loc[sgb, sn].values.astype(float)
        vals_p = abund.loc[sgb, sp].values.astype(float)

        # 直接用 Mann-Whitney U 与 v2 保持一致
        from scipy import stats as _stats

        _, p_val = _stats.mannwhitneyu(vals_n, vals_p, alternative="two-sided")

        v2._boxplot_on_ax(
            ax=ax,
            vals_n=vals_n,
            vals_p=vals_p,
            ylabel="Relative abundance",
            p_val=p_val,
            seed=42,
            title=sgb,
        )

    # 隐藏多余子图
    for j in range(n, nrows * ncols):
        row = j // ncols
        col = j % ncols
        axes[row, col].set_visible(False)

    plt.tight_layout()
    v2.save_fig(fig, str(fig_dir / "Fig3_SGB_boxplots"), dpi=dpi)
    plt.close(fig)


def main() -> None:
    """主函数：串联 SGB 差异分析、双层级对比与可视化。"""
    args = parse_args()
    setup_logging(args.output, args.verbose)
    np.random.seed(args.seed)

    logging.info("[START] SGB 差异分析与双层级对比 | seed=%d | dpi=%d", args.seed, args.dpi)

    # 1. SGB 差异分析
    sgb_da, abund_da, sn, sp = run_sgb_da(args.abund, args.group, args.seed)

    # 2. 读取种水平 DA 与 SGB 映射
    species_da = load_species_da(args.species_da)
    sgb_map = load_sgb_map(args.sgb_map)

    # 3. 构建双层级汇总表
    summary = build_cross_level_summary(sgb_da, species_da, sgb_map, q_cutoff=0.05)

    # 4. 导出表格与列表
    tables_dir = args.output / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    # 为差异结果表补充物种注释，便于直接在物种名称维度进行比较
    sgb_da_annot = sgb_da.merge(sgb_map, on="SGB", how="left")
    # 调整列顺序: [SGB, species, ...其余统计量...]
    cols = sgb_da_annot.columns.tolist()
    if "species" in cols:
        cols = ["SGB", "species"] + [c for c in cols if c not in ("SGB", "species")]
        sgb_da_annot = sgb_da_annot[cols]

    sgb_da_path = tables_dir / "Table_SGB_differential_bins.tsv"
    sgb_da_annot.to_csv(sgb_da_path, sep="\t", index=False)

    summary_path = tables_dir / "Table_SGB_vs_MetaPhlAn_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    logging.info("[OUTPUT] SGB 差异结果: %s", sgb_da_path)
    logging.info("[OUTPUT] 双层级对比汇总表: %s", summary_path)

    export_enriched_lists(summary, args.group_a, args.group_b, args.output)

    # 5. 绘图（CNS 级）
    plot_cross_level_scatter(summary, args.output, args.dpi)
    plot_enriched_waterfall(summary, args.group_a, args.group_b, args.output, args.dpi)
    plot_sgb_boxplots(abund_da, sn, sp, summary, args.output, args.dpi)

    logging.info("[DONE] 所有结果已写入: %s", args.output)


if __name__ == "__main__":
    main()

