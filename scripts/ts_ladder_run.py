#!/usr/bin/env python3
"""TS 阶梯的串行执行器：逐级跑、逐级过门禁、按诊断路由。

配套 ts_ladder.py 生成的作业集。设计纪律：
  * 一次只跑一个 ORCA 作业（先等 ORCA 空闲）
  * 每级提交前过 preflight 门禁,跑完过 --check-output 交付物校验
  * L1 之后先判"最低模式是否牵动声明的反应对"——不牵动就**不进 L2**
  * L2 之后按 ts_verify 的判决等级路由：SADDLE-1 才继续；
    FLAT-VALLEY → 停并接受参考能垒；ORDER-n / MINIMUM → 停并报告
  * 每一步都打印"声明的反应对距离"与"相对起点的成键变化"——
    这是 L2-c 几何自洽性与 L3 连通性的程序化证据

用法：
  python ts_ladder_run.py --dir ladder_P1bHE --reactive 56,57 --reactive 21,57
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ts_common import (  # noqa: E402
    connectivity_delta, read_xyz_frames,
)

ORCA = os.environ.get("ORCA_EXE") or shutil.which("orca") or "orca"
GATE = os.environ.get(
    "TS_PREFLIGHT",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "preflight_check.py"),
)
HARTREE = 627.509474


def log(path, msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')}|LADDER|{msg}"
    print(line, flush=True)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def orca_running():
    if os.name == "nt":
        r = subprocess.run(["powershell", "-NoProfile", "-Command",
                            "@(Get-Process orca -ErrorAction SilentlyContinue).Count"],
                           capture_output=True)
        try:
            return int(r.stdout.decode("utf-8", "replace").strip()) > 0
        except ValueError:
            return True
    try:
        r = subprocess.run(["pgrep", "-x", "orca"], capture_output=True)
        return r.returncode == 0
    except FileNotFoundError:
        return False


def pair_dist(path, pairs):
    fr = read_xyz_frames(path, last_only=True)
    if not fr:
        return {}
    a = fr[0]
    out = {}
    for i, j in pairs:
        if 1 <= i <= len(a) and 1 <= j <= len(a):
            dx = a[i - 1][1][0] - a[j - 1][1][0]
            dy = a[i - 1][1][1] - a[j - 1][1][1]
            dz = a[i - 1][1][2] - a[j - 1][1][2]
            out[(i, j)] = (dx * dx + dy * dy + dz * dz) ** 0.5
    return out


def run_job(d, job, logf):
    inp = os.path.join(d, job + ".inp")
    if not os.path.exists(inp):
        log(logf, f"{job}|MISSING-INPUT")
        return None
    r = subprocess.run([sys.executable, GATE, "--job", job, "--dir", d],
                       capture_output=True)
    g = r.stdout.decode("utf-8", "replace").strip()
    log(logf, f"{job}|GATE|{g.splitlines()[0] if g else '(none)'}")
    if r.returncode != 0:
        log(logf, f"{job}|REJECT|{g.replace(chr(10), ' ; ')[:160]}")
        return None
    log(logf, f"{job}|START")
    with open(os.path.join(d, job + ".out"), "w", encoding="utf-8") as out:
        subprocess.run([ORCA, job + ".inp"], cwd=d, stdout=out,
                       stderr=subprocess.STDOUT)
    c = subprocess.run([sys.executable, GATE, "--check-output", job, "--dir", d],
                       capture_output=True)
    ct = c.stdout.decode("utf-8", "replace").strip()
    log(logf, f"{job}|{ct.replace(chr(10), ' ; ')[:200]}")
    return os.path.join(d, job + ".out")


def verify(out, reactive):
    cmd = [sys.executable,
           os.path.join(os.path.dirname(os.path.abspath(__file__)), "ts_verify.py"),
           "--out", out]
    for p in reactive:
        cmd += ["--reactive", f"{p[0]},{p[1]}"]
    r = subprocess.run(cmd, capture_output=True)
    return r.stdout.decode("utf-8", "replace")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="ts_ladder.py 的 --outdir")
    ap.add_argument("--reactive", action="append", required=True, metavar="I,J",
                    help="声明的反应对(1-based),可重复")
    ap.add_argument("--no-l3", action="store_true", help="跳过 L3 向下弛豫")
    ap.add_argument("--no-sp", action="store_true", help="跳过 L5 单点")
    args = ap.parse_args()

    root = args.dir
    pairs = [tuple(int(x) for x in s.replace(" ", ",").split(","))
             for s in args.reactive]
    logf = os.path.join(root, "ladder.log")
    report = os.path.join(root, "REPORT.txt")
    if os.path.exists(logf):
        os.remove(logf)

    # 从目录名反推作业名(前缀)
    subs = [d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]
    l1 = next((d for d in subs if d.endswith("-L1-constOptFreq")), None)
    if not l1:
        sys.exit(f"{root} 下找不到 *-L1-constOptFreq 目录")
    name = l1[: -len("-L1-constOptFreq")]
    L1, L2, L3, SP = (f"{name}-L1-constOptFreq", f"{name}-L2-TSopt",
                      f"{name}-L3-relax", f"sp-{name}-L2-TS")

    log(logf, f"阶梯启动 name={name} 反应对={pairs}")
    log(logf, "等待 ORCA 空闲(不并行)...")
    while orca_running():
        time.sleep(60)
    log(logf, "ORCA 空闲")

    R = []   # 报告行

    def emit(s=""):
        R.append(s)
        print(s)
        with open(report, "w", encoding="utf-8") as fh:
            fh.write("\n".join(R) + "\n")

    emit(f"TS 阶梯报告  name={name}")
    emit(f"声明的反应对 (1-based): {', '.join(f'{a}-{b}' for a, b in pairs)}")
    emit()

    # ================= L1 =================
    emit("=== L1 约束优化 + Freq（产出 Hessian 与模式谱）===")
    out1 = run_job(os.path.join(root, L1), L1, logf)
    if not out1:
        emit("L1 未通过门禁 —— 修输入后重跑。停止。")
        return 1
    d1 = pair_dist(os.path.join(root, L1, L1 + ".xyz"), pairs)
    emit("  L1 收敛点上各反应对距离: "
         + ", ".join(f"{a}-{b} = {v:.3f} A" for (a, b), v in d1.items()))
    v1 = verify(out1, pairs)
    v1f = os.path.join(root, "L1_verify.txt")
    open(v1f, "w", encoding="utf-8").write(v1)
    emit("  L1 判决: " + next((l.strip() for l in v1.splitlines()
                               if l.startswith("VERDICT")), "(无)"))
    # 门禁 L1：最低模式是否牵动声明的反应对
    involve = all(any(f"{a}-{b}" == ln.split()[0] or f"{b}-{a}" == ln.split()[0]
                      for ln in v1.splitlines())
                  for a, b in pairs) if pairs else False
    flat = "FLAT-VALLEY" in v1
    emit(f"  门禁 L1（最低模式是否牵动声明的反应对）: "
         f"{'PASS' if (involve or flat) else 'FAIL'}")
    if not involve and not flat:
        emit()
        emit(">> 出口：L1 的最低模式**不是**你的反应坐标 —— 禁止进 L2")
        emit(">> （InHess Read 会把优化器带上错模式，重演 §6.1 案例 1）")
        emit(">> 补救：改 --hold 的约束值 / 换方位重做 L1；或走 L4 换初猜来源")
        emit(">> （L4-a GOAT-TS / L4-b UniTS-Gen 多采样 / L4-c 弛豫扫描）")
        emit()
        emit(f"L1 已产出可用参考点：{os.path.join(root, L1)}（参考能垒从此处取）")
        return 0

    # ================= L2 =================
    emit()
    emit("=== L2 真鞍点（OptTS + Freq，读 L1 的 Hessian）===")
    # 把 L1 的 hess/gbw/xyz 让 L2 通过相对路径读到——同一 outdir 下已成立
    out2 = run_job(os.path.join(root, L2), L2, logf)
    if not out2:
        emit("L2 未通过门禁 —— 检查 L1 的 .hess/.gbw 是否产出。停止。")
        return 1
    v2 = verify(out2, pairs)
    open(os.path.join(root, "L2_verify.txt"), "w", encoding="utf-8").write(v2)
    for l in v2.splitlines():
        emit("  " + l)
    verdict = next((l.split()[1] for l in v2.splitlines()
                    if l.startswith("VERDICT:")), "?")
    emit()
    # L2-c：几何自洽性证据
    p2 = os.path.join(root, L2, L2 + ".xyz")
    d2 = pair_dist(p2, pairs)
    emit("=== 门禁 L2-c 几何自洽性（把几何读出来，逐个核对键长）===")
    emit("  L2 鞍点上各反应对距离: "
         + ", ".join(f"{a}-{b} = {v:.3f} A" for (a, b), v in d2.items()))
    st = read_xyz_frames(os.path.join(root, L1, L1 + ".xyz"), last_only=True)
    ts = read_xyz_frames(p2, last_only=True)
    if st and ts:
        gained, lost = connectivity_delta(st[0], ts[0])
        emit(f"  相对 L1，成键变化: 新成 {gained if gained else '无'}"
             f" / 断开 {lost if lost else '无'}")
        emit("  ⚠️ 这一步是必须人工读的：SADDLE-1 + 命中反应对**也可能不是你要的反应**"
             "（§6.1 案例 1）。把上面数值与预期的反应物/产物对照。")
    emit()
    if verdict == "FLAT-VALLEY":
        emit(">> 出口：目标通道本身是**宽平过渡谷**（真实物理，不是收敛失败）。")
        emit(">> 停。：用末态 Freq 作参考能垒，并在文中如实写明'参考能垒'。")
        emit(">> **不要再重启 OptTS** —— 换个初猜也变不出鞍点。")
        return 0
    if verdict != "SADDLE-1":
        emit(f">> 出口：判决 {verdict}，不是一阶鞍点。")
        emit(">> ORDER-n → 回 L1 调约束/换方位；MINIMUM → 约束值太靠产物或反应物。")
        emit(">> 或走 L4 换初猜来源。")
        return 0

    # ================= L3 =================
    if not args.no_l3:
        emit()
        emit("=== L3 连通性（从 TS 向下弛豫，看落到哪个极小点）===")
        out3 = run_job(os.path.join(root, L3), L3, logf)
        if out3:
            p3 = os.path.join(root, L3, L3 + ".xyz")
            d3 = pair_dist(p3, pairs)
            emit("  弛豫末态各反应对距离: "
                 + ", ".join(f"{a}-{b} = {v:.3f} A" for (a, b), v in d3.items()))
            if ts:
                en = read_xyz_frames(p3, last_only=True)
                if en:
                    g2, l2_ = connectivity_delta(ts[0], en[0])
                    emit(f"  相对 TS，成键变化: 新成 {g2 if g2 else '无'}"
                         f" / 断开 {l2_ if l2_ else '无'}")
            emit("  门禁 L3：这些变化必须与'预期产物'一致。要更严格就上 IRC")
            emit("  （`! IRC` + `%irc InitHess Read / Hess_Filename \"<L2>.hess\""
                 " / Direction both end`，IRC 不能与内坐标约束同用）")

    # ================= L5 =================
    if not args.no_sp:
        emit()
        emit("=== L5 高一级单点（保留 L2 的电子态设置）===")
        out5 = run_job(os.path.join(root, SP), SP, logf)
        if out5:
            t = open(out5, encoding="utf-8", errors="ignore").read()
            import re
            E = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
            es = re.findall(r"Total Energy after SMD CDS correction\s*=\s*(-?\d+\.\d+)", t)
            emit(f"  单点能: E = {E[-1] if E else '?'} Eh"
                 + (f"   (SMD 校正后 {es[-1]} Eh)" if es else ""))
            s2 = re.findall(r"Total S\*\*2\s+\.*\s+(\d+\.\d+)", t)
            emit(f"  ⟨S²⟩ = {s2[-1] if s2 else '未打印'}"
                 "   ← 必须与 L2 的 ⟨S²⟩ 一致（BS 单重态应 ~0.2–0.4，不是 0）")
    emit()
    emit("阶梯完成。把 L2 的能量与正确的参考态相减得能垒；")
    emit("跨'极小点 vs 鞍点'比较一律用电子能 ΔE（方法学铁律）。")
    log(logf, "LADDER-DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
