"""
Frame generator + visualizer, driven by the scheduler.
=====================================================
Two modes:
  * round -> one local syndrome-extraction round.
  * merge -> two rounds of a remote lattice-surgery merge. One communication ion
             per lane cycles herald -> deliver -> measure -> re-herald, so the
             two-round run shows the comm ion routing holds across QEC cycles.

Physical rules the frames obey:
  * a swap is a two-ion crystal rotated in a well, drawn merge -> rotate -> split.
  * junctions live only in the gate zone. The interface zone holds cavities and
    comm ions, no junction. An idle right-boundary ancilla that blocks a comm
    lane parks in the bottom cell's spare well (a gate-zone move) for the merge.
  * a comm ion shuttles out and back through the SPAM zone, one frame each way.
A live verifier checks every frame for overlap.
Run:  python3 qec_visualizer.py
"""
from __future__ import annotations
import json
from qec_scheduler import build_stabilizers, place, corner_step, stab_cell, num, seam_schedule, round_ops

D = 3
MERGE_ROUNDS = 2

# ----- human labels for the d=3 checks --------------------------------------
NAMES = {frozenset({1, 2, 4, 5}): "A1", frozenset({2, 3, 5, 6}): "A2",
         frozenset({4, 5, 7, 8}): "A3", frozenset({5, 6, 8, 9}): "A4",
         frozenset({2, 3}): "B1", frozenset({7, 8}): "B2",
         frozenset({1, 4}): "B3", frozenset({6, 9}): "B4"}
STABS = build_stabilizers(D)
LABEL = {s: NAMES[frozenset(num(rc, D) for rc in s.data)] for s in STABS}
KIND = {s: s.kind for s in STABS}
CELL = stab_cell(D)
SEAM = seam_schedule(D)                             # {s: {same, cross, junction}}

# ----- geometry -------------------------------------------------------------
CY = {0: 90, 1: 220, 2: 350}
XS = 150.0
def X(col): return 185 + (col + 0.5) * XS
def gap_y(a, b): return (CY[a] + CY[b]) / 2
def JX(c): return X(c) + XS / 4                     # gate-zone junction, in the gap
XIF = X(D - 0.5) + 140                              # interface zone (cavities), past the SPAM/readout end
def data_col(n):  return (n - 1) % D
def data_row(n):  return (n - 1) // D

ACOL = {}
CHAINS = place(D)
for i, cell in enumerate(CHAINS):
    for j, (kind, item) in enumerate(cell):
        if kind == "anc":
            lft = cell[j - 1] if j > 0 else None
            rgt = cell[j + 1] if j < len(cell) - 1 else None
            if lft and lft[0] == "data" and rgt and rgt[0] == "data":
                ACOL[item] = (lft[1][1] + rgt[1][1]) / 2
            elif lft is None:
                ACOL[item] = -0.5
            else:
                ACOL[item] = D - 0.5

ID_OF_STAB = {s: LABEL[s] for s in STABS}
ID_OF_DATA = {(data_row(n), data_col(n)): "d%d" % n for n in range(1, D * D + 1)}
IW = 11                                             # in-well half-offset


def is_right_boundary(s):
    return s.weight == 2 and s.kind == "Z" and all(c == D - 1 for r, c in s.data)


def base_ions(merge):
    """Ion table + home positions. In merge mode add one comm ion per lane,
    active lanes 0..d-2 (one per seam check) and the bottom lane held spare."""
    ions, home = {}, {}
    for n in range(1, D * D + 1):
        ions["d%d" % n] = ["d%d" % n, "data"]
        home["d%d" % n] = [X(data_col(n)), CY[data_row(n)]]
    for s in STABS:
        ions[LABEL[s]] = [LABEL[s], "X" if KIND[s] == "X" else "Z"]
        home[LABEL[s]] = [X(ACOL[s]), CY[CELL[s]]]
    if merge:
        for r in range(D):
            ions["C%d" % r] = ["C%d" % r, "comm" if r in SEAM else "spare"]
            home["C%d" % r] = [XIF, CY[r]]
    return ions, home


def build(merge=False, rounds=1):
    ions, home = base_ions(merge)
    pos = {i: home[i][:] for i in ions}
    FR = []
    parked = [(LABEL[s], CELL[s]) for s in STABS if is_right_boundary(s)] if merge else []
    parked_labels = {lab for lab, _ in parked}

    def snap(cap, hi=None, junc=None, badge="", merged=None):
        FR.append({"pos": {i: pos[i][:] for i in pos}, "cap": cap, "hi": hi or [],
                   "junc": junc or [], "badge": badge,
                   "merged": [sorted(p) for p in (merged or [])]})

    def setp(i, x, y): pos[i] = [x, y]

    def swap_pairs(pairs, cap, junc=None):
        """A swap = a two-ion crystal rotated in a well: merge, rotate, split."""
        if not pairs:
            return
        orig = {p: (pos[p[0]][:], pos[p[1]][:]) for p in pairs}
        hi = [a for a, b in pairs]
        for a, b in pairs:
            pa, pb = orig[(a, b)]; w = (pa[0] + pb[0]) / 2
            setp(a, w - IW, pa[1]); setp(b, w + IW, pb[1])
        snap(cap + " merge into one well.", merged=pairs, hi=hi, junc=junc)
        for a, b in pairs:
            pa, pb = orig[(a, b)]; w = (pa[0] + pb[0]) / 2
            setp(a, w + IW, pa[1]); setp(b, w - IW, pb[1])
        snap(cap + " the two-ion crystal rotates, swapping them.", merged=pairs, hi=hi, junc=junc)
        for a, b in pairs:
            pa, pb = orig[(a, b)]
            setp(a, pb[0], pb[1]); setp(b, pa[0], pa[1])
        snap(cap + " split; the two have exchanged wells.", hi=hi)

    # --- translate scheduler identifiers to visualizer ion ids -------------
    def A(s):   return LABEL[s]                       # ancilla stab -> label
    def Dt(rc): return ID_OF_DATA[rc]                 # data (r,c) -> "dN"
    def C(l):   return "C%d" % l                      # lane index -> comm ion id
    cx = X(D - 1)

    # --- one handler per operation the scheduler emits.  The sequence, the
    #     batching, and the beats all come from round_ops; this only places them.
    for op in round_ops(D, merge, rounds):
        v = op[0]
        if v == "prep":
            snap(f"Start of the {'merge' if op[1] else 'round'}: data and ancillas in memory wells, ancillas in their basis."
                 + (" Comm ions hold heralded Bell pairs at the cavities." if op[1] else ""))
        elif v == "round":
            snap(f"Round {op[1] + 1} of {op[2]}.", badge=f"round {op[1] + 1}")
        elif v == "park":
            jn = [(D - 1, cl) for _, cl in op[1]]; hi = [A(s) for s, _ in op[1]]
            for s, cl in op[1]:
                setp(A(s), JX(D - 1), gap_y(cl, cl + 1))
            snap("The idle right-boundary ancilla lifts into the gate-zone junction toward the bottom cell.", junc=jn, hi=hi)
            for s, cl in op[1]:
                setp(A(s), X(D - 0.5), CY[D - 1])
            snap("It settles in the bottom cell's spare well, clearing its comm lane for the whole merge.", junc=jn, hi=hi)
        elif v == "unpark":
            jn = [(D - 1, cl) for _, cl in op[1]]; hi = [A(s) for s, _ in op[1]]
            for s, cl in op[1]:
                setp(A(s), JX(D - 1), gap_y(cl, cl + 1))
            snap("The merge ends; the right-boundary ancilla lifts back into the gate-zone junction.", junc=jn, hi=hi)
            for s, cl in op[1]:
                setp(A(s), *home[A(s)])
            snap("It returns to its home well.", junc=jn, hi=hi)
        elif v == "inrow":
            L = op[1] + 1; hi = [A(s) for s, _ in op[2]]
            for s, rc in op[2]:
                d = Dt(rc); wx, wy = home[d]; setp(d, wx - IW, wy); setp(A(s), wx + IW, wy)
            snap(f"Step {L}: in-row ancillas step into their data's wells and gate.", hi=hi,
                 merged=[(A(s), Dt(rc)) for s, rc in op[2]])
            for s, rc in op[2]:
                d = Dt(rc); setp(d, *home[d]); setp(A(s), *home[A(s)])
            snap(f"Step {L}: they step back to their own wells.", hi=hi)
        elif v == "swap":
            L = op[1] + 1
            tag = "target-column data" if op[3] == "onto-column" else "own well"
            swap_pairs([(A(s), Dt(rc)) for s, rc in op[2]], f"Step {L}: the cross-row ancilla and its {tag}")
        elif v == "xlift":
            L = op[1] + 1
            for s, c, ac, tr in op[2]:
                setp(A(s), JX(c), gap_y(ac, tr))
            snap(f"Step {L}: the ancillas lift into the junctions in the gaps.", hi=[A(s) for s, *_ in op[2]],
                 junc=[(c, min(ac, tr)) for s, c, ac, tr in op[2]], badge="in transit")
        elif v == "xgate":
            L = op[1] + 1; lst = [(A(s), Dt(rc)) for s, rc in op[2]]
            for aid, did in lst:
                gx, gy = pos[did]; setp(did, gx - IW, gy); setp(aid, gx + IW, gy)
            snap(f"Step {L}: each lands in its data's well and the gate fires.", hi=[a for a, _ in lst], merged=lst)
        elif v == "xlower":
            L = op[1] + 1
            for s, c, ac, tr in op[2]:
                did = Dt((tr, c)); setp(did, pos[did][0] + IW, pos[did][1]); setp(A(s), JX(c), gap_y(ac, tr))
            snap(f"Step {L}: they lift back into the junctions.", hi=[A(s) for s, *_ in op[2]],
                 junc=[(c, min(ac, tr)) for s, c, ac, tr in op[2]])
        elif v == "xdrop":
            L = op[1] + 1
            for s, c, ac in op[2]:
                setp(A(s), X(c), CY[ac])
            snap(f"Step {L}: they drop back onto the target column.", hi=[A(s) for s, *_ in op[2]])
        elif v == "comm_out":
            L = op[1] + 1; lanes = op[2]; mode = op[3]
            tgt = (XIF + (cx if mode == "same" else JX(D - 1))) / 2
            for l in lanes:
                setp(C(l), tgt, CY[l])
            snap(f"Step {L}: comm ions shuttle out through the SPAM zone"
                 + (" to their same-cell boundary data." if mode == "same" else " toward their junctions."),
                 hi=[C(l) for l in lanes])
        elif v == "comm_lift":
            L = op[1] + 1
            for l in op[2]:
                setp(C(l), JX(D - 1), gap_y(l, l + 1))
            snap(f"Step {L}: they lift into the boundary junctions in the gap.", hi=[C(l) for l in op[2]],
                 junc=[(D - 1, l) for l in op[2]], badge="in transit")
        elif v == "comm_gate":
            L = op[1] + 1; mode = op[3]
            lst = [(C(l), Dt(rc), rc[0]) for l, rc in op[2]]
            for cid, did, row in lst:
                setp(cid, cx + IW, CY[row]); setp(did, cx - IW, CY[row])
            jn = [(D - 1, l) for l, _ in op[2]] if mode == "cross" else None
            snap(f"Step {L}: the seam gate fires in the boundary data's well.", hi=[c for c, _, _ in lst],
                 merged=[(c, d) for c, d, _ in lst], junc=jn)
        elif v == "comm_lower":
            L = op[1] + 1
            for l in op[2]:
                setp(Dt((l + 1, D - 1)), cx, CY[l + 1]); setp(C(l), JX(D - 1), gap_y(l, l + 1))
            snap(f"Step {L}: they lift back into the junctions.", hi=[C(l) for l in op[2]],
                 junc=[(D - 1, l) for l in op[2]])
        elif v == "comm_back":
            L = op[1] + 1; lanes = op[2]; mode = op[3]
            tgt = (XIF + (cx if mode == "same" else JX(D - 1))) / 2
            for l in lanes:
                if mode == "same":
                    setp(Dt((l, D - 1)), cx, CY[l])
                setp(C(l), tgt, CY[l])
            snap(f"Step {L}: they shuttle back through the SPAM zone.", hi=[C(l) for l in lanes])
        elif v == "comm_arrive":
            L = op[1] + 1
            for l in op[2]:
                setp(C(l), XIF, CY[l])
            snap(f"Step {L}: they are back at their cavities.", hi=[C(l) for l in op[2]])
        elif v == "measure":
            snap("The comm ions are measured at their cavities; this round's Bell pairs reach module B.",
                 hi=[C(l) for l in op[1]], badge=f"{len(op[1])} pairs -> B")
        elif v == "readout":
            swap_pairs([(A(s), Dt(rc)) for s, rc in op[1]], "Readout: an ancilla and its neighbour data")
        elif v == "syndromes":
            snap("Every ancilla is at the SPAM end; 493 nm readout gives the syndrome bits"
                 + (", the seam ones included." if op[1] else "."),
                 hi=[LABEL[s] for s in STABS if not (merge and is_right_boundary(s))], badge="syndromes")
        elif v == "reset":
            swap_pairs([(A(s), Dt(rc)) for s, rc in op[1]], "Reset: an ancilla and its data")
        elif v == "reset_done":
            snap(f"Ancillas are reset and back in their memory wells for round {op[1] + 2}.",
                 hi=[LABEL[s] for s in STABS if not (merge and is_right_boundary(s))])
        elif v == "reherald":
            snap("The comm ions re-herald a fresh Bell pair at their cavities.", hi=[C(l) for l in op[1]])
    return FR, ions, home


def bad_frame(f):
    mg = {frozenset(p) for p in f["merged"]}
    by = {}
    for i, (x, y) in f["pos"].items():
        by.setdefault(round(y), []).append((x, i))
    for y, arr in by.items():
        arr.sort()
        for k in range(1, len(arr)):
            if abs(arr[k][0] - arr[k - 1][0]) < 30 and frozenset((arr[k][1], arr[k - 1][1])) not in mg:
                return f"{arr[k-1][1]} & {arr[k][1]} at y={y}"
    return None


HTML = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>__TITLE__</title><style>
:root{--bg:#faf9f5;--panel:#fff;--ink:#2c2c2a;--mut:#6b6a64;--line:#d9d7cd;
--x:#3B7FD4;--xd:#1C4F8C;--z:#D2703A;--zd:#8F4620;--purple:#534AB7;--teal:#0F6E56;--amber:#BA7517;--red:#A32D2D;--redbg:#FCEBEB;}
@media(prefers-color-scheme:dark){:root{--bg:#22211e;--panel:#2c2b27;--ink:#eceae2;--mut:#a3a299;--line:#45453f;--redbg:#3a1414;}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:-apple-system,system-ui,"Segoe UI",sans-serif}
.wrap{max-width:1040px;margin:0 auto;padding:16px}h1{font-size:18px;margin:2px 0}
.sub{font-size:13px;color:var(--mut);margin-bottom:12px}
.bar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:8px}
button{border:1px solid var(--line);background:var(--panel);color:var(--ink);border-radius:8px;padding:6px 12px;font-size:14px;cursor:pointer}
button:disabled{opacity:.4}button.play{background:var(--purple);color:#fff;border-color:var(--purple)}
#slider{flex:1;min-width:150px}.sname{font-size:13px;font-weight:600;min-width:74px}
.badge{font-size:12px;font-weight:600;padding:3px 10px;border-radius:12px;background:#E1F5EE;color:var(--teal)}
.verify{font-size:12px;font-weight:600;padding:3px 10px;border-radius:12px;background:#E1F5EE;color:var(--teal)}
.verify.bad{background:var(--redbg);color:var(--red)}
.cap{font-size:13.5px;color:var(--mut);min-height:40px;margin:4px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:6px}svg{width:100%;height:auto;display:block}
.legend{display:flex;gap:16px;font-size:11.5px;color:var(--mut);margin-top:8px;flex-wrap:wrap}
.legend i{display:inline-block;width:12px;height:12px;vertical-align:-1px;margin-right:5px}
</style></head><body><div class="wrap">
<h1>__TITLE__</h1>
<div class="sub">__SUB__ A live verifier checks every frame for overlap.</div>
<div class="bar"><button id="prev">&lsaquo; Prev</button><button id="next">Next &rsaquo;</button>
<button id="play" class="play">&#9654; Play</button><input id="slider" type="range" min="0" max="0" value="0"/>
<span class="sname" id="sname"></span><span class="badge" id="badge"></span><span class="verify" id="verify"></span></div>
<div class="cap" id="cap"></div>
<div class="card"><svg id="stage" role="img"></svg></div>
<div class="legend"><span><i style="border-radius:50%;background:var(--panel);border:1px solid var(--line)"></i>data</span>
<span><i style="border-radius:2px;background:var(--x)"></i>X ancilla</span>
<span><i style="border-radius:2px;background:var(--z)"></i>Z ancilla</span>
__COMMLEG__</div>
</div><script>
const DATA=__DATA__;
const NS="http://www.w3.org/2000/svg",E=(t,a)=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;};
const $=i=>document.getElementById(i),svg=$("stage");
const W=DATA.xhi,H=DATA.celly[DATA.celly.length-1]+70;
svg.setAttribute("viewBox","0 0 "+W+" "+H);
const gZ=E("g",{}),gWell=E("g",{}),gJ=E("g",{}),gW=E("g",{}),gI=E("g",{});
[gZ,gWell,gJ,gW,gI].forEach(g=>svg.appendChild(g));
DATA.celly.forEach((y,ci)=>{gZ.appendChild(E("rect",{x:DATA.xlo,y:y-24,width:W-DATA.xlo-20,height:48,rx:10,fill:"none",stroke:"var(--line)","stroke-dasharray":"5 4",opacity:.55}));
 const t=E("text",{x:DATA.xlo+5,y:y-31,"font-size":11,fill:"var(--mut)"});t.textContent="cell "+ci;gZ.appendChild(t);});
if(DATA.xif){const t=E("text",{x:DATA.xif-6,y:DATA.celly[0]-31,"font-size":11,fill:"var(--mut)"});t.textContent="interface (cavities)";gZ.appendChild(t);
 DATA.celly.forEach((y,ci)=>gZ.appendChild(E("path",{d:"M "+(DATA.xif+30)+" "+(y-14)+" A 17 17 0 0 1 "+(DATA.xif+30)+" "+(y+14),fill:"none",stroke:"var(--teal)","stroke-width":2,opacity:.55})));}
(DATA.wells||[]).forEach(w=>{gWell.appendChild(E("rect",{x:w[0]-22,y:w[1]-19,width:44,height:38,rx:10,fill:"var(--panel)",stroke:"var(--line)","stroke-width":1.2}));});
const jel={};DATA.junctions.forEach(j=>{const lane=E("line",{x1:j.x,y1:j.y1+24,x2:j.x,y2:j.y2-24,stroke:"var(--line)","stroke-width":2.5,opacity:.5,"stroke-linecap":"round"});
 gJ.appendChild(lane);gJ.appendChild(E("circle",{cx:j.x,cy:j.y1+24,r:3.4,fill:"var(--mut)"}));gJ.appendChild(E("circle",{cx:j.x,cy:j.y2-24,r:3.4,fill:"var(--mut)"}));jel[j.c+"_"+j.b]=lane;});
const el={};for(const id in DATA.ions){const lab=DATA.ions[id][0],typ=DATA.ions[id][1],g=E("g",{});g.style.transition="transform .3s ease";
 if(typ==="data"){g.appendChild(E("circle",{cx:0,cy:0,r:15,fill:"var(--panel)",stroke:"var(--line)","stroke-width":1.4}));
  const t=E("text",{x:0,y:4,"text-anchor":"middle","font-size":10,fill:"var(--ink)"});t.textContent=lab;g.appendChild(t);el[id]={g};}
 else{const col=typ==="X"?"var(--x)":typ==="Z"?"var(--z)":typ==="comm"?"var(--teal)":"var(--line)";
  const dk=typ==="X"?"var(--xd)":typ==="Z"?"var(--zd)":typ==="comm"?"#0a3a2e":"var(--mut)";
  const r=E("rect",{x:-14,y:-13,width:28,height:26,rx:5,fill:col,stroke:dk,"stroke-width":1});g.appendChild(r);
  const t=E("text",{x:0,y:4,"text-anchor":"middle","font-size":9,fill:(typ==="spare")?"var(--mut)":"#fff","font-weight":600});t.textContent=lab;g.appendChild(t);el[id]={g,rc:r,dk:dk};}
 gI.appendChild(g);}
function verify(f){const mg=new Set((f.merged||[]).map(p=>p.join("|"))),by={};
 for(const id in f.pos){const x=f.pos[id][0],y=f.pos[id][1];(by[Math.round(y)]=by[Math.round(y)]||[]).push([x,id]);}
 for(const y in by){const a=by[y].sort((p,q)=>p[0]-q[0]);
  for(let i=1;i<a.length;i++){if(Math.abs(a[i][0]-a[i-1][0])<30){const k=[a[i-1][1],a[i][1]].sort().join("|");if(!mg.has(k))return a[i-1][1]+" & "+a[i][1];}}}return null;}
let step=0;function render(){const f=DATA.frames[step];
 $("cap").textContent=f.cap;$("badge").textContent=f.badge||"";$("sname").textContent=(step+1)+" / "+DATA.frames.length;$("slider").value=step;
 const v=verify(f),ve=$("verify");if(v){ve.textContent="⚠ overlap "+v;ve.classList.add("bad");}else{ve.textContent="✓ no overlap";ve.classList.remove("bad");}
 for(const k in jel){jel[k].setAttribute("stroke","var(--line)");jel[k].setAttribute("stroke-width",2.5);jel[k].setAttribute("opacity",.5);}
 (f.junc||[]).forEach(j=>{const k=j[0]+"_"+j[1];if(jel[k]){jel[k].setAttribute("stroke","var(--amber)");jel[k].setAttribute("stroke-width",4);jel[k].setAttribute("opacity",1);}});
 while(gW.firstChild)gW.removeChild(gW.firstChild);
 (f.merged||[]).forEach(pr=>{const a=f.pos[pr[0]],b=f.pos[pr[1]];if(!a||!b)return;
  const x=Math.min(a[0],b[0])-20,y=Math.min(a[1],b[1])-19,w=Math.abs(a[0]-b[0])+40,h=Math.abs(a[1]-b[1])+38;
  gW.appendChild(E("rect",{x:x,y:y,width:w,height:h,rx:11,fill:"var(--purple)","fill-opacity":.10,stroke:"var(--purple)","stroke-width":1.8,"stroke-dasharray":"5 3"}));
  const t=E("text",{x:x+w/2,y:y-5,"text-anchor":"middle","font-size":9.5,fill:"var(--purple)","font-weight":600});t.textContent="well";gW.appendChild(t);});
 for(const id in el){const p=f.pos[id];if(p)el[id].g.setAttribute("transform","translate("+p[0]+","+p[1]+")");
  if(el[id].rc){const on=f.hi.indexOf(id)>=0;el[id].rc.setAttribute("stroke-width",on?2.8:1);el[id].rc.setAttribute("stroke",on?"var(--amber)":el[id].dk);}}
 $("prev").disabled=step===0;$("next").disabled=step===DATA.frames.length-1;}
$("slider").max=DATA.frames.length-1;
$("next").onclick=()=>{if(step<DATA.frames.length-1){step++;render();}};
$("prev").onclick=()=>{if(step>0){step--;render();}};
$("slider").oninput=e=>{step=+e.target.value;render();};
let timer=null;$("play").onclick=()=>{if(timer){clearInterval(timer);timer=null;$("play").innerHTML="&#9654; Play";return;}
 $("play").innerHTML="&#10073;&#10073; Pause";timer=setInterval(()=>{if(step<DATA.frames.length-1){step++;render();}else{clearInterval(timer);timer=null;$("play").innerHTML="&#9654; Play";}},640);};
document.addEventListener("keydown",e=>{if(e.key==="ArrowRight")$("next").click();if(e.key==="ArrowLeft")$("prev").click();});
render();
</script></body></html>"""


def write_html(path, merge, rounds):
    FR, ions, home = build(merge, rounds)
    bad = [(k, bad_frame(f)) for k, f in enumerate(FR) if bad_frame(f)]
    wells = [home[i] for i in ions if ions[i][1] in ("data", "X", "Z")]
    data = {"ions": ions, "frames": FR, "celly": [CY[0], CY[1], CY[2]],
            "junctions": [{"c": c, "b": b, "x": JX(c), "y1": CY[b], "y2": CY[b + 1]}
                          for c in range(D) for b in range(D - 1)],
            "wells": wells, "xlo": X(-0.5) - 45,
            "xhi": (XIF + 70) if merge else (X(D - 0.5) + 60), "xif": XIF if merge else 0}
    title = (f"Distance-3 remote merge - {rounds} rounds, one comm ion per lane"
             if merge else "Distance-3 error-correction round - balanced {3,3,2}")
    sub = ("Two full merge rounds. One comm ion per lane cycles herald, deliver, measure, re-herald; "
           "junctions are gate-zone only and physical swaps happen in wells."
           if merge else "One code row per cell; ions sit in wells and swap in wells, junctions link neighbouring cells.")
    commleg = '<span><i style="border-radius:2px;background:var(--teal)"></i>comm ion</span>' if merge else ""
    html = (HTML.replace("__TITLE__", title).replace("__SUB__", sub)
            .replace("__COMMLEG__", commleg).replace("__DATA__", json.dumps(data)))
    with open(path, "w") as f:
        f.write(html)
    return len(FR), bad


if __name__ == "__main__":
    for merge, rounds, path in [(False, 1, "qec_round_sim_d3.html"),
                                (True, MERGE_ROUNDS, "qec_merge_full_sim_d3.html")]:
        n, bad = write_html(path, merge, rounds)
        tag = "merge" if merge else "round"
        print(f"{tag:6s}: {n} frames, {len(bad)} overlaps -> {path}")
        for k, b in bad[:8]:
            print(f"   frame {k}: {b}")