# S12 概念图 v2：数值从 SM Table S9 md 表解析，不手输
# 脚本要点：
#   1) 图内角落加 "Posterior medians only; uncertainty not shown" 提示
#   2) 画像标签加 (descriptive) 限定，与正文 descriptive profile 口径一致
#   3) 输出后同步复制 results/figures/（存档位）
import re, pathlib, shutil

root = pathlib.Path(__file__).resolve().parents[2]  # released: repo root
sm = (root / "submission/bmc_medicine/manuscript/06_Supplementary_Materials.md").read_text()
pat = r"\| (Noon|Evening) Main \| (mu|sigma) \| (\w+)_c \| ([-−0-9.]+) \| ([^|]+) \| ([0-9.]+|—) \|"
rows = {}
for m in re.finditer(pat, sm):
    sess, sub, item = m.group(1).lower(), m.group(2), m.group(3)
    med = float(m.group(4).replace('−', '-'))
    expv = float(m.group(6)) if m.group(6) != "—" else None
    rows[(sess, sub, item)] = (med, expv)
if len(rows) < 16:
    raise SystemExit(f"S9 解析行数不足: {len(rows)} (期望 16) — S9 表格式可能已变，停手核查")

def sigma_pct(sess, item):
    _, expv = rows[(sess, "sigma", item)]
    return (expv - 1) * 100

items = {
    ("noon", "aspiration"): "self-development (aspiration)", ("noon", "sleep_ontime"): "sleep on time",
    ("noon", "breakfast"): "breakfast", ("noon", "wake_early"): "early rising",
    ("evening", "exercise"): "exercise", ("evening", "breathing"): "breathing",
    ("evening", "walk"): "walking", ("evening", "dinner"): "dinner",
}

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5.6))
for (sess, item), label in items.items():
    mu = rows[(sess, "mu", item)][0]
    sp = sigma_pct(sess, item)
    c = "#2a7f62" if sess == "noon" else "#8a5a2b"
    ax.scatter(mu, sp, color=c, s=70, zorder=3)
    ax.annotate(label, (mu, sp), textcoords="offset points", xytext=(8, 6 if sp > 0 else -14), fontsize=9)
ax.axhline(0, color="gray", lw=0.8, ls="--")
ax.axvline(0, color="gray", lw=0.8, ls="--")
ax.set_xlabel("Mean-mood association (mu coefficient, points)")
ax.set_ylabel("Residual-SD association (%)  [− = compression, + = inflation]")
ax.set_title("Item-level location–scale profiles by session (conceptual overview)")
ax.text(0.55, -18, "descriptive stabilizer profile: σ↓ with modest μ", fontsize=9, color="#2a7f62")
ax.text(0.55, 9, "descriptive activator profile: μ↑ with σ↑", fontsize=9, color="#8a5a2b")
ax.annotate("aspiration (σ −30.5%, μ ≈ 0)", (0.021, -30.5), textcoords="offset points", xytext=(28, -4), fontsize=9, color="#2a7f62")
ax.annotate("exercise (μ +0.776, σ +12.9%)", (0.776, 12.9), textcoords="offset points", xytext=(-160, 12), fontsize=9, color="#8a5a2b")
ax.text(0.985, 0.02, "Posterior medians only; uncertainty not shown",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, style="italic", color="0.3")
fig.tight_layout()
out = root / "submission/bmc_medicine/figures/fig_s12_conceptual"
fig.savefig(str(out) + ".png", dpi=300)
fig.savefig(str(out) + ".svg")
# 注：submission/bmc_medicine/figures/ 与 results/figures/ 的 fig_s12_conceptual 为同一
# inode（硬链接），保存一处即两处同步，无需复制（SameFileError 实证）。
print("saved", out, "(hardlinked to results/figures/)")
