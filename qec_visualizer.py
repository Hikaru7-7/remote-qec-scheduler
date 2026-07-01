"""
Frame generator + visualizer, driven by the scheduler.
=====================================================
One program, two modes:
  * round  -> one local syndrome-extraction round (in-row gates, cross-row
              transits, readout).
  * merge  -> the same round plus the seam extraction of a remote lattice-surgery
              merge: each active comm ion comes out of the I/F zone, gates its
              same-cell boundary data, crosses a junction to gate its cross-cell
              data, and returns to be measured, while a herald ion holds each
              cavity. Steps come from qec_scheduler.seam_schedule.
Both write a self-contained HTML with an independent overlap check.
Run:  python3 qec_visualizer.py   (writes qec_round_sim_d3.html and qec_merge_full_sim_d3.html)
"""
from __future__ import annotations
import json
from qec_scheduler import build_stabilizers, place, corner_step, stab_cell, num, seam_schedule

D = 3

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
def JX(c): return X(c) + XS / 4
XIF = X(2 * D - 1.5) + 90                           # I/F zone, past the SPAM/readout end (a cell holds up to 2d ions)
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


def base_ions(merge):
    """The ion table and home positions. In merge mode add the comm ions: one
    delivering ion per lane at the cavity, plus a herald ion behind it, and the
    spare lane."""
    ions, home = {}, {}
    for n in range(1, D * D + 1):
        ions["d%d" % n] = ["d%d" % n, "data"]
        home["d%d" % n] = [X(data_col(n)), CY[data_row(n)]]
    for s in STABS:
        ions[LABEL[s]] = [LABEL[s], "X" if KIND[s] == "X" else "Z"]
        home[LABEL[s]] = [X(ACOL[s]), CY[CELL[s]]]
    if merge:
        for r in range(D):
            spare = (r == D - 1)
            ions["C%d" % r] = ["C%d" % r, "spare" if spare else "comm"]
            home["C%d" % r] = [XIF, CY[r]]
            if not spare:                            # herald ion behind the deliverer
                ions["C%d'" % r] = ["C%d'" % r, "herald"]
                home["C%d'" % r] = [XIF + 34, CY[r]]
    return ions, home


def build(merge=False):
    """Generate the frames for one round (merge=False) or one merge round."""
    ions, home = base_ions(merge)
    pos = {i: home[i][:] for i in ions}
    FR = []

    def snap(cap, hi=None, junc=None, badge="", merged=None):
        FR.append({"pos": {i: pos[i][:] for i in pos}, "cap": cap, "hi": hi or [],
                   "junc": junc or [], "badge": badge,
                   "merged": [sorted(p) for p in (merged or [])]})

    def setp(i, x, y): pos[i] = [x, y]

    def gates(step):
        inrow, cross = [], []
        for s in STABS:
            for (r, c) in s.data:
                if corner_step(s, (r, c)) != step:
                    continue
                g = (ID_OF_STAB[s], ID_OF_DATA[(r, c)], CELL[s], r, c)
                (inrow if r == CELL[s] else cross).append(g)
        return inrow, cross

    kind_txt = "merge round" if merge else "round"
    snap(f"Start of the {kind_txt}: data and ancillas in memory wells; ancillas pumped to their basis."
         + (" Comm ions wait at the cavities holding heralded Bell pairs." if merge else ""))

    for step in range(4):
        L = step + 1
        inrow, cross = gates(step)
        # --- local in-row gates: step into the well, gate fires, step back ---
        if inrow:
            for a, d, _, _, _ in inrow:
                wx, wy = home[d][0], home[d][1]
                setp(d, wx - IW, wy); setp(a, wx + IW, wy)
            mp = [(a, d) for a, d, _, _, _ in inrow]
            snap(f"Step {L}: in-row ancillas step into their data's wells.", hi=[a for a, *_ in inrow], merged=mp)
            snap(f"Step {L}: the in-row gates fire in the wells.", hi=[a for a, *_ in inrow], merged=mp)
            for a, d, _, _, _ in inrow:
                setp(d, home[d][0], home[d][1]); setp(a, home[a][0], home[a][1])
            snap(f"Step {L}: they step back to their own wells.", hi=[a for a, *_ in inrow])
        # --- local cross-row gates: swap on, lift through junction, land, gate ---
        if cross:
            for a, d, ac, tr, c in cross:
                dcol = ID_OF_DATA[(ac, c)]; ax = pos[a][0]
                setp(a, pos[dcol][0], pos[a][1]); setp(dcol, ax, pos[a][1])
            snap(f"Step {L}: cross-row ancillas swap onto the target column.", hi=[a for a, *_ in cross])
            for a, d, ac, tr, c in cross:
                setp(a, JX(c), gap_y(ac, tr))
            snap(f"Step {L}: they lift into the junction in the gap.", hi=[a for a, *_ in cross],
                 junc=[(c, min(ac, tr)) for a, d, ac, tr, c in cross], badge="in transit")
            for a, d, ac, tr, c in cross:
                tid = ID_OF_DATA[(tr, c)]; setp(tid, X(c) - IW, CY[tr]); setp(a, X(c) + IW, CY[tr])
            snap(f"Step {L}: each lands in its data's well and the gate fires.", hi=[a for a, *_ in cross],
                 merged=[(a, ID_OF_DATA[(tr, c)]) for a, d, ac, tr, c in cross])
            for a, d, ac, tr, c in cross:
                setp(ID_OF_DATA[(tr, c)], X(c), CY[tr]); setp(a, JX(c), gap_y(ac, tr))
            snap(f"Step {L}: split; the ancilla lifts back into the junction.", hi=[a for a, *_ in cross],
                 junc=[(c, min(ac, tr)) for a, d, ac, tr, c in cross])
            for a, d, ac, tr, c in cross:
                setp(a, X(c), CY[ac]); dcol = ID_OF_DATA[(ac, c)]
                setp(a, home[a][0], home[a][1]); setp(dcol, X(c), CY[ac])
            snap(f"Step {L}: swap back; the chain order is restored.", hi=[a for a, *_ in cross])
        # --- seam extraction on the same step (merge only): each comm ion comes
        #     out to the data, couples, and returns to its cavity so it never
        #     lingers in a cell. The herald ions stay put at the cavities.
        if merge:
            same = [s2 for s2, g in SEAM.items() if g["same"] == step]
            crs = [s2 for s2, g in SEAM.items() if g["cross"] == step]
            if same:
                for s2 in same:
                    setp("C%d" % s2, X(D - 1) + 2 * IW, CY[s2])
                snap(f"Step {L}: comm ions come out and gate their same-cell boundary data (herald ions hold the cavities).",
                     hi=["C%d" % s2 for s2 in same],
                     merged=[("C%d" % s2, ID_OF_DATA[(s2, D - 1)]) for s2 in same])
                for s2 in same:
                    setp("C%d" % s2, XIF, CY[s2])
            if crs:
                for s2 in crs:
                    setp("C%d" % s2, X(D - 1) + 2 * IW, CY[s2 + 1])
                snap(f"Step {L}: comm ions cross their junction and gate their cross-cell data.",
                     hi=["C%d" % s2 for s2 in crs], junc=[(D - 1, s2) for s2 in crs],
                     merged=[("C%d" % s2, ID_OF_DATA[(s2 + 1, D - 1)]) for s2 in crs])
                for s2 in crs:
                    setp("C%d" % s2, XIF, CY[s2])

    if merge:
        for s2 in SEAM:
            setp("C%d" % s2, XIF, CY[s2])                    # deliverers return to the cavity
        snap("The comm ions return and are measured; each Bell pair teleports its half-parity to module B.",
             hi=["C%d" % s2 for s2 in SEAM], badge=f"{len(SEAM)} pairs -> B")

    # --- readout: bubble ancillas to the SPAM end ---------------------------
    order = [[it for it in cell] for cell in CHAINS]

    def bubbled(c):
        seen = False
        for k, _ in c:
            if k == "anc": seen = True
            elif seen: return False
        return True

    while not all(bubbled(c) for c in order):
        for c in order:
            k = 0
            while k < len(c) - 1:
                if c[k][0] == "anc" and c[k + 1][0] == "data":
                    c[k], c[k + 1] = c[k + 1], c[k]; k += 2
                else:
                    k += 1
        for ci, c in enumerate(order):
            x = -0.5
            for kind, item in c:
                iid = ID_OF_STAB[item] if kind == "anc" else "d%d" % num(item, D)
                setp(iid, X(x), CY[ci]); x += 1
        snap("Readout: ancillas bubble to the SPAM end (adjacent swaps).",
             hi=[LABEL[s] for s in STABS], badge="swap-out")
    snap("Every ancilla is at the SPAM end; 493 nm readout gives all syndrome bits"
         + (", the seam ones included." if merge else "."),
         hi=[LABEL[s] for s in STABS], badge="syndromes")
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
#slider{flex:1;min-width:150px}.sname{font-size:13px;font-weight:600;min-width:64px}
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
if(DATA.xif){const t=E("text",{x:DATA.xif+4,y:DATA.celly[0]-31,"font-size":11,fill:"var(--mut)"});t.textContent="I/F (cavities)";gZ.appendChild(t);
 DATA.celly.forEach((y,ci)=>{if(ci<DATA.celly.length-1||DATA.celly.length===1)gZ.appendChild(E("path",{d:"M "+(DATA.xif+50)+" "+(y-14)+" A 17 17 0 0 1 "+(DATA.xif+50)+" "+(y+14),fill:"none",stroke:"var(--teal)","stroke-width":2,opacity:.6}));});}
(DATA.wells||[]).forEach(w=>{gWell.appendChild(E("rect",{x:w[0]-22,y:w[1]-19,width:44,height:38,rx:10,fill:"var(--panel)",stroke:"var(--line)","stroke-width":1.2}));});
const jel={};DATA.junctions.forEach(j=>{const lane=E("line",{x1:j.x,y1:j.y1+24,x2:j.x,y2:j.y2-24,stroke:"var(--line)","stroke-width":2.5,opacity:.5,"stroke-linecap":"round"});
 gJ.appendChild(lane);gJ.appendChild(E("circle",{cx:j.x,cy:j.y1+24,r:3.4,fill:"var(--mut)"}));gJ.appendChild(E("circle",{cx:j.x,cy:j.y2-24,r:3.4,fill:"var(--mut)"}));jel[j.c+"_"+j.b]=lane;});
const el={};for(const id in DATA.ions){const lab=DATA.ions[id][0],typ=DATA.ions[id][1],g=E("g",{});g.style.transition="transform .3s ease";
 if(typ==="data"){g.appendChild(E("circle",{cx:0,cy:0,r:15,fill:"var(--panel)",stroke:"var(--line)","stroke-width":1.4}));
  const t=E("text",{x:0,y:4,"text-anchor":"middle","font-size":10,fill:"var(--ink)"});t.textContent=lab;g.appendChild(t);el[id]={g};}
 else{const col=typ==="X"?"var(--x)":typ==="Z"?"var(--z)":typ==="comm"?"var(--teal)":typ==="herald"?"var(--panel)":"var(--line)";
  const dk=typ==="X"?"var(--xd)":typ==="Z"?"var(--zd)":typ==="comm"?"#0a3a2e":typ==="herald"?"var(--teal)":"var(--mut)";
  const r=E("rect",{x:-14,y:-13,width:28,height:26,rx:5,fill:col,stroke:dk,"stroke-width":typ==="herald"?2:1});g.appendChild(r);
  const t=E("text",{x:0,y:4,"text-anchor":"middle","font-size":9,fill:(typ==="herald"||typ==="spare")?"var(--mut)":"#fff","font-weight":600});t.textContent=lab;g.appendChild(t);el[id]={g,rc:r,dk:dk};}
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
  const t=E("text",{x:x+w/2,y:y-5,"text-anchor":"middle","font-size":9.5,fill:"var(--purple)","font-weight":600});t.textContent="gate well";gW.appendChild(t);});
 for(const id in el){const p=f.pos[id];if(p)el[id].g.setAttribute("transform","translate("+p[0]+","+p[1]+")");
  if(el[id].rc){const on=f.hi.indexOf(id)>=0;el[id].rc.setAttribute("stroke-width",on?2.8:1);el[id].rc.setAttribute("stroke",on?"var(--amber)":el[id].dk);}}
 $("prev").disabled=step===0;$("next").disabled=step===DATA.frames.length-1;}
$("slider").max=DATA.frames.length-1;
$("next").onclick=()=>{if(step<DATA.frames.length-1){step++;render();}};
$("prev").onclick=()=>{if(step>0){step--;render();}};
$("slider").oninput=e=>{step=+e.target.value;render();};
let timer=null;$("play").onclick=()=>{if(timer){clearInterval(timer);timer=null;$("play").innerHTML="&#9654; Play";return;}
 $("play").innerHTML="&#10073;&#10073; Pause";timer=setInterval(()=>{if(step<DATA.frames.length-1){step++;render();}else{clearInterval(timer);timer=null;$("play").innerHTML="&#9654; Play";}},820);};
document.addEventListener("keydown",e=>{if(e.key==="ArrowRight")$("next").click();if(e.key==="ArrowLeft")$("prev").click();});
render();
</script></body></html>"""


def write_html(path, merge):
    FR, ions, home = build(merge)
    bad = [(k, bad_frame(f)) for k, f in enumerate(FR) if bad_frame(f)]
    wells = [home[i] for i in ions if ions[i][1] in ("data", "X", "Z")]
    data = {"ions": ions, "frames": FR, "celly": [CY[0], CY[1], CY[2]],
            "junctions": [{"c": c, "b": b, "x": JX(c), "y1": CY[b], "y2": CY[b + 1]}
                          for c in range(D) for b in range(D - 1)],
            "wells": wells, "xlo": X(-0.5) - 45,
            "xhi": (XIF + 90) if merge else (X(D - 0.5) + 60), "xif": XIF if merge else 0}
    title = ("Distance-3 merge round - local syndrome extraction + seam"
             if merge else "Distance-3 error-correction round - balanced {3,3,2}")
    sub = ("The local round runs, and each active comm ion comes out of its cavity to gate its two boundary "
           "data (one across a junction) and carry the parity to the other module."
           if merge else "One code row per cell; ions sit in wells, junctions link neighbouring cells, gates fire off the junctions.")
    commleg = ('<span><i style="border-radius:2px;background:var(--teal)"></i>comm ion</span>'
               '<span><i style="border-radius:2px;background:var(--panel);border:2px solid var(--teal)"></i>herald ion</span>') if merge else ""
    html = (HTML.replace("__TITLE__", title).replace("__SUB__", sub)
            .replace("__COMMLEG__", commleg).replace("__DATA__", json.dumps(data)))
    with open(path, "w") as f:
        f.write(html)
    return len(FR), bad


if __name__ == "__main__":
    for merge, path in [(False, "qec_round_sim_d3.html"), (True, "qec_merge_full_sim_d3.html")]:
        n, bad = write_html(path, merge)
        tag = "merge" if merge else "round"
        print(f"{tag:6s}: {n} frames, {len(bad)} overlaps -> {path}")
        for k, b in bad[:6]:
            print(f"   frame {k}: {b}")