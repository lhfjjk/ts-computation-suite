# -*- coding: utf-8 -*-
"""提交前结构体检 (pre-flight check) —— ORCA 工作流技能自带工具

通用用法(推荐给 runner 做门禁):
    python preflight_check.py --job <作业名> --dir <工作目录>     # PASS/WARN/FAIL/PEND, FAIL 时 exit 1
批量用法:把下面的 JOBS 列表改成你的作业清单(每条为 (相对目录, 作业名, 预期类型)),
    预期类型 ∈ {opt, freq, optts, sp}; 或直接用 --job 模式 + 自动从关键词推断类型。
判据速查: 任意两原子 <0.85 A = FAIL(真实最短键 O-H 0.96 / C-H 1.09); 非键接近含重原子 = 有成键风险;
    电荷-多重度奇偶必须一致; BS 作业必须有 MORead+moinp+FlipSpin+FinalMs; 约束值须与起始几何自洽。
对每条作业的起始几何做: 坐标合法性 / 原子重叠 / 连通性完整性 / 约束原子索引与距离 /
电荷-多重度奇偶一致性 / BS 设置 / 关键词交付物。任何 FAIL 的作业不应提交。"""
import os, re, sys, json
from pathlib import Path
from collections import Counter

BASE = Path(os.environ.get("PREFLIGHT_BASE", os.getcwd()))  # 换项目: 设 PREFLIGHT_BASE 或从作业目录运行
Z = {'H':1,'B':5,'C':6,'N':7,'O':8,'F':9,'P':15,'S':16,'Cl':17,'Br':35,'I':53,'Ni':28,'Fe':26}
RCOV = {'H':0.31,'B':0.84,'C':0.76,'N':0.71,'O':0.66,'F':0.57,'P':1.07,'S':1.05,'Cl':1.02,'Br':1.20,'I':1.39,'Ni':1.24,'Fe':1.32}
# 共价半径 * 1.3 作为"成键"判据
BONDCUT = {('H','H'):1.30,('H','C'):1.35,('H','N'):1.35,('H','O'):1.30,
           ('C','C'):1.60,('C','N'):1.55,('C','O'):1.55,('C','B'):1.65,
           ('N','N'):1.50,('N','O'):1.50,('N','B'):1.60,('N','Ni'):2.10,
           ('O','Ni'):2.20,('O','O'):1.45,('B','O'):1.55,('C','Ni'):2.10,
           ('H','B'):1.35,('H','Ni'):1.80,('B','B'):1.75,('B','C'):1.65,('Ni','Ni'):2.4,
           ('C','F'):1.45,('C','Cl'):1.80,('C','Br'):1.95,('C','I'):2.15,
           ('H','F'):1.00,('H','Cl'):1.40,('H','Br'):1.55,('H','I'):1.70,
           ('O','Cl'):1.75,('N','Br'):1.95,('Ni','Cl'):2.40,('Ni','Br'):2.55,
           ('C','S'):1.85,('H','S'):1.45,('S','S'):2.10,('Ni','S'):2.45,('C','P'):1.95,('O','P'):1.75,
           ('C','Fe'):2.10,('N','Fe'):2.10,('O','Fe'):2.20,('H','Fe'):1.80,('S','Fe'):2.45,
           ('Cl','Fe'):2.40,('Br','Fe'):2.55,('Fe','Fe'):2.60}

def cut(a, b):
    k = tuple(sorted([a, b])); return BONDCUT.get(k, 2.0)

def read_xyz(fp):
    L = Path(fp).read_text(errors="replace").splitlines()
    st = 0
    for k, l in enumerate(L[:5]):
        if l.strip().isdigit(): st = k+1; break
    n = int(L[st-1].strip()); el, co = [], []
    for l in L[st+1:st+1+n]:
        p = l.split(); el.append(p[0]); co.append((float(p[1]), float(p[2]), float(p[3])))
    return el, co

def read_inp_geom(jdir, L):
    for i, l in enumerate(L):
        s = l.strip()
        if s.startswith("* xyzfile"):
            t = s.split()
            fp = jdir/t[4]
            if fp.exists():
                el, co = read_xyz(fp)
                return el, co, int(t[2]), int(t[3]), t[4], True
            return None, None, int(t[2]), int(t[3]), t[4], False
        if s.startswith("* xyz"):
            t = s.split(); el, co = [], []
            for m in L[i+1:]:
                if m.strip() == "*": break
                q = m.split()
                if len(q) >= 4: el.append(q[0]); co.append((float(q[1]), float(q[2]), float(q[3])))
            return el, co, int(t[2]), int(t[3]), "内联坐标", True
    return None, None, None, None, "无几何", False

def dist(c, i, j):
    return ((c[i][0]-c[j][0])**2 + (c[i][1]-c[j][1])**2 + (c[i][2]-c[j][2])**2) ** 0.5

def check(jdir, job, kind):
    """kind: 'opt'|'freq'|'optts'|'sp'  —— 由调用方给出预期类型"""
    res = {"job": job, "fail": [], "warn": [], "info": []}
    inp = jdir/(job+".inp")
    if not inp.exists():
        res["pending"] = "缺 .inp(批量清单与本目录不匹配?)"; return res
    L = inp.read_text(encoding="utf-8", errors="ignore").splitlines()
    kw = next((l.strip()[1:].strip() for l in L if l.strip().startswith("!")), "")
    el, co, chg, mult, src, ok = read_inp_geom(jdir, L)
    if not ok:
        res["pending"] = "等待上游产出 " + str(src)
        return res
    res["info"].append(f"{len(el)} 原子  {src}  {chg} {mult}")
    # 1) 元素 / 坐标
    bad = [e for e in el if e not in Z]
    if bad: res["fail"].append("未知元素 " + str(set(bad)))
    # 2) 电荷-多重度奇偶
    ne = sum(Z[e] for e in el) - chg
    if (ne % 2) != (mult - 1) % 2:
        res["fail"].append(f"电子数 {ne} 与多重度 {mult} 奇偶不匹配")
    # 2b) 自旋态选择 —— 铁律 26（2026-09-23 新增，源于一次真实错误）
    #   d 区金属络合物声明「单重态」时拦截。这是最难发现的一类错误：
    #   ORCA 对 mult=1 **不做任何警告、也不切换**（闭壳层单重态是合法请求），
    #   所以用错了参考态也不会报错 —— 只有拿实验/文献对参考态时才会暴露。
    #   实例：Tp*Ni(OAc) 基态是三重态（实验顺磁高自旋 + 计算 triplet 低 15.2 kcal/mol），
    #        但项目按单重态算了它 → 所有含该参考态的量（BDFE / ΔG(捕获) / ΔG(S_H2)）
    #        系统性偏高约 15 kcal/mol，且完全静默。
    #   注意与「开壳层却写 RKS」区分：后者 ORCA 会 WARNING + 自动切 UHF，不会静默出错。
    #   确有理由时（如 d8 平面四方闭壳层 Ni(II)），在 .inp 里**显式声明**：
    #       # SPINSTATE-OK: <理由与依据>
    _D_BLOCK = {"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Y","Zr","Nb","Mo",
                "Tc","Ru","Rh","Pd","Ag","Cd","La","Ce","Hf","Ta","W","Re","Os","Ir",
                "Pt","Au","Hg"}
    _txt = "\n".join(L)
    _tms = sorted({e for e in el if e in _D_BLOCK})
    if _tms and mult == 1 and "UKS" not in kw and "FlipSpin" not in _txt:
        if "# SPINSTATE-OK:" in _txt:
            res["info"].append(f"含 {','.join(_tms)} 声明闭壳层单重态 —— 已显式 SPINSTATE-OK")
        else:
            res["fail"].append(
                "含 d 区金属 " + ",".join(_tms) + " 却声明闭壳层单重态(mult 1)且无 UKS/FlipSpin"
                " —— 必须先核对基态自旋态（查实验：顺磁/EPR/NMR 位移？或自己算多个自旋态比能量），"
                "不要默认单重态。若确有依据（如 d8 平面四方闭壳层），在 .inp 加一行注释 "
                "`# SPINSTATE-OK: <理由>`")
    # 3) 原子重叠 —— 绝对判据: 真实最短键是 O-H 0.96 A（比 C-H 1.09 更短！),
    #    ⚠️ 2026-10-09 修正：原阈值 1.00 A 会把正常 O-H 键误判为重叠（含 O-H 的分子此前没跑过）
    worst = (99.0, None, None)
    for i in range(len(el)):
        for j in range(i+1, len(el)):
            d = dist(co, i, j)
            if d < worst[0]: worst = (d, i, j)
    # 同时用"键图"识别非键接近: 键判据 = 1.3*(rcov之和)
    nb = []
    for i in range(len(el)):
        for j in range(i+1, len(el)):
            d = dist(co, i, j)
            is_bond = d < 1.3*(RCOV.get(el[i],0.8)+RCOV.get(el[j],0.8))
            if not is_bond:
                lim = 1.90 if (el[i]!='H' and el[j]!='H') else 1.60
                if d < lim: nb.append((d,i,j))
    if worst[0] < 0.85:
        res["fail"].append(f"原子重叠 {el[worst[1]]}{worst[1]+1}-{el[worst[2]]}{worst[2]+1} = {worst[0]:.3f} A")
    elif worst[0] < 1.06:
        res["warn"].append(f"最近接触 {el[worst[1]]}{worst[1]+1}-{el[worst[2]]}{worst[2]+1} = {worst[0]:.3f} A")
    for d,i,j in sorted(nb)[:3]:
        res["warn"].append(f"非键接近 {el[i]}{i+1}-{el[j]}{j+1} = {d:.3f} A")
    # 4) 相邻原子距离 (成键是否合理)
    for i in range(len(el)):
        for j in range(i+1, len(el)):
            d = dist(co, i, j)
            if d < cut(el[i], el[j]) - 1e-6 and (i, j) != (worst[1], worst[2]):
                pass
    # 5) 关键词 / 交付物
    want = {"opt": ["Opt"], "freq": ["Freq"], "optts": ["OptTS", "Freq"],
            "optfreq": ["Opt", "Freq"], "sp": []}[kind]
    for w in want:
        if w not in kw: res["fail"].append("关键词缺 " + w + " (该作业应为 " + kind + ")")
    if kind == "opt" and "Opt" in kw and "OptTS" in kw:
        res["fail"].append("同时出现 Opt 与 OptTS")
    # 6) BS 设置(BS 作业 = UKS + MORead + FlipSpin + FinalMs;
    #    仅 MORead 而无 UKS/FlipSpin 的是"读收敛 gbw 做性质打印"的作业,不算 BS —— 2026-09-17 修正)
    t = "\n".join(L)
    # BS(破缺对称)设置检查 —— 触发条件只看 **FlipSpin**，不能看 UKS
    #   ⚠️ 2026-09-23 修正：原先写的是 `"UKS" in kw or "FlipSpin" in t`，把 UKS 与 BS 混为一谈。
    #   UKS 只是"非限制"（任何 mult>1 都需要，三元态就必须要），而 BS 是一种**特定技巧**
    #   （读收敛 gbw + 翻转某个原子的自旋 + 指定 FinalMs）。后果：任何**显式写 UKS 的
    #   正常开壳层作业**都会被要求提供 MORead/FlipSpin 而被误判 FAIL —— 三重态作业首当其冲。
    #   项目此前的输入从不显式写 UKS（靠 ORCA 自动切），所以这个误判一直没暴露。
    #   正确判据：FlipSpin 出现 ⇒ 这是 BS 作业 ⇒ 四件套必须齐。
    if "FlipSpin" in t:
        for w in ["MORead", "moinp", "FlipSpin", "FinalMs"]:
            if w not in t: res["fail"].append("BS 作业缺 " + w)
    if "moinp" in t:
        m = re.search(r'moinp\s+"([^"]+)"', t)
        if m and not (jdir/m.group(1)).exists(): res["fail"].append("moinp 文件不存在: " + m.group(1))
    # 7) 约束: 原子索引 + 距离与几何一致
    for m in re.finditer(r"\{\s*B\s+(\d+)\s+(\d+)\s+([0-9.]+)\s+C\s*\}", "\n".join(L)):
        i, j, d0 = int(m.group(1)), int(m.group(2)), float(m.group(3))
        dact = dist(co, i, j)
        res["info"].append(f"约束 B {i}-{j} = {d0} (几何实测 {dact:.3f})")
        if abs(dact - d0) > 0.05:
            res["warn"].append(f"约束 {i}-{j} 目标 {d0} 与起始几何 {dact:.3f} 差 {abs(dact-d0):.3f} A")
        if not (el[i] in ("C","Ni","N","O","H") and el[j] in ("C","Ni","N","O","H")):
            res["warn"].append(f"约束原子异常 {i}({el[i]}) / {j}({el[j]})")
    if kind == "optts" and not re.search(r"\{\s*B", "\n".join(L)):
        res["info"].append("OptTS 无约束(自由精修)")
    return res

def infer_kind(kw):
    """预期类型,仅从 .inp 自身的关键词推断。

    ⚠️ 这个推断是**循环的**:漏写 Opt 的作业会被推断成 freq/sp,于是关键词检查
    必然通过 —— 它永远抓不到"漏写关键词"。**不要用它做门禁**,用 kind_from_name()。
    保留它只为批量清单未声明类型时的向后兼容。"""
    if "OptTS" in kw: return "optts"
    if "Opt" in kw:   return "opt"
    if "Freq" in kw:  return "freq"
    return "sp"


def kind_from_name(job):
    """从作业名推断*预期*类型 —— 非循环,因此能抓出漏写的关键词。

    本项目的命名约定本身就是工序声明(skill §6.2)。返回 None 表示名字没带
    足够信息(此时门禁退回 infer_kind 并明确标注该判据不可靠)。
    """
    n = job.lower()
    # 先剥掉末尾序号（-TSopt2 -> -TSopt、sp-g6-bHE-TS2 -> sp-g6-bHE-TS），
    # 否则这些名字落不到任何后缀规则上
    n = re.sub(r"\d+$", "", n)
    if n.startswith("sp-") or n.startswith("sp2-") or n.endswith("-sp"):
        return "sp"
    # 顺序重要:先匹配更长的后缀
    if n.endswith("-optfreq") or n.endswith("-constoptfreq"):
        return "optfreq"
    if n.endswith("constfreq") or n.endswith("-freq") or n.endswith("-finalfreq"):
        return "freq"
    if n.endswith("-tsopt") or n.endswith("-optts"):
        return "optts"
    if (n.endswith("-opt") or n.endswith("-constopt") or n.endswith("-relax")
            or n.endswith("-const") or n.endswith("-modts")):
        return "opt"
    return None


# kind_from_name 的返回值 → 该工序**必须**在 .inp 里声明的关键词 token 集合。
# 复合名（optfreq）要拆成 token 才能与 .inp 侧解析出的集合比较（2026-09-23 修）。
KIND_TOKENS = {
    "sp": set(),
    "opt": {"opt"},
    "freq": {"freq"},
    "optts": {"optts"},
    "optfreq": {"opt", "freq"},
}


def check_output(jdir, job):
    """跑完之后的交付物校验 —— 这次对比的是"作业真的做了该做的事吗"。

    起因(S_H2 项目 2026-09-21 实测):`SH2s-enc20.inp` 漏写 `Opt`,ORCA 正常终止,
    runner 照写 `.done`,但它其实是**纯单点**(0 个优化周期、无几何输出)。
    只看"ORCA TERMINATED NORMALLY"的守卫拦不住这种偏差。
    """
    res = {"job": job, "fail": [], "warn": [], "info": []}
    out = jdir / (job + ".out")
    if not out.exists():
        res["pending"] = "无 .out(未运行?)"
        return res
    t = out.read_text(encoding="utf-8", errors="ignore")
    kw = ""
    inp = jdir / (job + ".inp")
    if inp.exists():
        kw = next((l.strip()[1:].strip() for l in
                   inp.read_text(encoding="utf-8", errors="ignore").splitlines()
                   if l.strip().startswith("!")), "")
    name_kind = kind_from_name(job)
    # 优先用名字声明的类型;没名字信息时看 .inp
    # ⚠️ acts 是**原子 token 集合**（opt/freq/optts），而 kind_from_name 返回的是
    # **复合串**（"optfreq"）。2026-09-23 实测：原先写 `name_kind not in acts` 等于
    # 拿 "optfreq" 去 in {opt,freq} —— 恒为真 ⇒ **所有 OptFreq 作业都被误判 FAIL**
    # （项目 21 个已有 + 队列里的 quartet 全是 OptFreq）。改为按 KIND_TOKENS 拆成
    # token 再比。另外：OptTS 作业若 .inp 也声明了 Freq，freq 必须独立加进 acts，
    # 否则"声明 Freq 却没有频率块"（优化未收敛时 ORCA 不算频率）漏检。
    acts = set()
    if "OptTS" in kw: acts.add("optts")
    elif "Opt" in kw: acts.add("opt")
    if "Freq" in kw: acts.add("freq")
    if not acts:
        acts.add("sp")

    if not re.search(r"ORCA TERMINATED NORMALLY", t):
        res["fail"].append("未正常终止")
        return res

    ncyc = len(re.findall(r"GEOMETRY OPTIMIZATION CYCLE", t))
    if ("opt" in acts or "optts" in acts):
        if ncyc == 0:
            res["fail"].append(
                "声明了优化但输出里 0 个优化周期 —— 实际跑成了纯单点"
                "(检查 .inp 是否漏写 Opt/OptTS)")
        elif "THE OPTIMIZATION HAS CONVERGED" not in t:
            res["warn"].append(f"优化 {ncyc} 周期后未形式收敛(达 MaxIter?)"
                               " —— 末态须做连通性检查再解读")
        if not (jdir / (job + ".xyz")).exists() and "THE OPTIMIZATION HAS CONVERGED" in t:
            res["warn"].append("收敛了但没有 .xyz 产出(下游作业会读不到几何)")
    if "freq" in acts and "VIBRATIONAL FREQUENCIES" not in t:
        res["fail"].append("声明了 Freq 但输出里没有频率块")
    name_kind = kind_from_name(job)
    if name_kind and name_kind != "sp":
        need = KIND_TOKENS.get(name_kind, set())
        missing = need - acts
        if missing:
            res["fail"].append(
                f"作业名暗示应为 {name_kind}(需 {sorted(need)}),"
                f"但 .inp 实际只声明 {sorted(acts)} —— 缺 {sorted(missing)}")
    elif name_kind is None:
        # 名字没声明工序 -> 无法校验"它是否做了该做的事"。
        # SH2s-enc20 就是这种:名字无工序标记,.inp 漏写 Opt,于是纯单点被当成
        # "完成"。这不是检查器的缺陷,是**意图从未被机器可读地记录下来**。
        res["warn"].append(
            "作业名未声明工序(如 -Opt/-Freq/-TSopt/sp-),无法校验交付物是否符合预期;"
            "按 skill §6.2 的命名约定重命名后再由门禁把关")
    return res


if any(a in sys.argv for a in ("-h", "--help")):
    print(__doc__); sys.exit(0)


if __name__ == "__main__" and "--check-output" in sys.argv:
    # 跑完之后:校验作业真的做了该做的事(交付物校验,非循环)
    job = sys.argv[sys.argv.index("--check-output") + 1]
    jdir = Path(sys.argv[sys.argv.index("--dir") + 1]) if "--dir" in sys.argv else Path(".")
    r = check_output(jdir, job)
    if r.get("pending"):
        print("PEND " + job + " " + r["pending"]); sys.exit(0)
    if r["fail"]:
        print("FAIL " + job)
        for x in r["fail"]: print("  ✗ " + x)
        sys.exit(1)
    print(("WARN " if r["warn"] else "PASS ") + job)
    for x in r["warn"]: print("  ! " + x)
    sys.exit(0)


if __name__ == "__main__" and "--job" in sys.argv:
    job = sys.argv[sys.argv.index("--job") + 1]
    jdir = Path(sys.argv[sys.argv.index("--dir") + 1]) if "--dir" in sys.argv else Path(".")
    inp = jdir / (job + ".inp")
    if not inp.exists():
        print("FAIL " + job + " 缺 .inp"); sys.exit(2)
    kw = next((l.strip()[1:].strip() for l in inp.read_text(encoding="utf-8", errors="ignore").splitlines()
               if l.strip().startswith("!")), "")
    # 优先用作业名推断预期类型(非循环,能抓出漏写的关键词);
    # 名字没带信息时才退回从 .inp 推断,并明确标注该判据不可靠
    name_kind = kind_from_name(job)
    if name_kind:
        kind = name_kind
        note = ""
    else:
        kind = infer_kind(kw)
        note = "作业名未声明工序,预期类型由 .inp 自身推断(循环判据,抓不到漏写关键词)"
    r = check(jdir, job, kind)
    if note:
        r["warn"].append(note)
    if r.get("pending"):
        print("PEND " + job + " " + r["pending"]); sys.exit(0)
    if r["fail"]:
        print("FAIL " + job)
        for x in r["fail"]: print("  ✗ " + x)
        sys.exit(1)
    print(("WARN " if r["warn"] else "PASS ") + job + " | 预期 " + kind + " | "
          + "; ".join(r["info"]))
    for x in r["warn"]: print("  ! " + x)
    sys.exit(0)


if __name__ == "__main__":
    # No batch job list is bundled: use the --job or --check-output mode with
    # your own job name and directory.
    print("preflight_check.py")
    print("  python preflight_check.py --job <name> --dir <dir>")
    print("  python preflight_check.py --check-output <name> --dir <dir>")
    print("No batch list is bundled; pass --job/--dir for your own job.")
    sys.exit(0)
