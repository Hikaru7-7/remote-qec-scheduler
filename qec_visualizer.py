"""
Frame generator + visualizer, driven by the scheduler.
=====================================================
Reads the placement and colouring from qec_scheduler.py, generates the ion-move
frames for one d=3 error-correction round on the balanced {3,3,2} layout with
the gate-zone geometry (data and ancillas in wells; junctions between cells;
gates off junctions), checks that no two ions ever overlap, and writes a
self-contained HTML that animates the frames with its own independent overlap
check.  Run:  python3 qec_visualizer.py
"""
from __future__ import annotations
import json
from qec_scheduler import build_stabilizers, place, corner_step, stab_cell, num

D = 3

# ----- human labels for the d=3 checks (by the data they read) --------------
NAMES = {
    frozenset({1, 2, 4, 5}): "A1", frozenset({2, 3, 5, 6}): "A2",
    frozenset({4, 5, 7, 8}): "A3", frozenset({5, 6, 8, 9}): "A4",
    frozenset({2, 3}): "B1", frozenset({7, 8}): "B2",
    frozenset({1, 4}): "B3", frozenset({6, 9}): "B4",
}
STABS = build_stabilizers(D)
LABEL = {s: NAMES[frozenset(num(rc, D) for rc in s.data)] for s in STABS}
KIND = {s: s.kind for s in STABS}
CELL = stab_cell(D)

# ----- geometry -------------------------------------------------------------
CY = {0: 90, 1: 220, 2: 350}                       # y-pixel of each cell axis
XS = 150.0                                         # column pitch (wide enough to leave gaps)
def X(col): return 185 + (col + 0.5) * XS          # x-pixel of a column
def gap_y(a, b): return (CY[a] + CY[b]) / 2        # gate-zone strip between cells
def JX(c): return X(c) + XS / 4                    # junction sits at the MIDPOINT between two
#                                                    wells, so it is on top of neither ion

def data_col(n):  return (n - 1) % D               # column of data qubit n (1..9)
def data_row(n):  return (n - 1) // D

# column of every ancilla, from the placement chain (gap = midpoint, ends offset)
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

# ----- ion table: id -> (label, type) ; home position -----------------------
ions = {}          # id -> {"lab","type"}
home = {}          # id -> [x,y]
for n in range(1, D * D + 1):
    i = "d%d" % n
    ions[i] = {"lab": "d%d" % n, "type": "data"}
    home[i] = [X(data_col(n)), CY[data_row(n)]]
for s in STABS:
    i = LABEL[s]
    ions[i] = {"lab": i, "type": "X" if KIND[s] == "X" else "Z"}
    home[i] = [X(ACOL[s]), CY[CELL[s]]]

ID_OF_STAB = {s: LABEL[s] for s in STABS}
ID_OF_DATA = {(data_row(n), data_col(n)): "d%d" % n for n in range(1, D * D + 1)}

# ----- frame engine ---------------------------------------------------------
pos = {i: home[i][:] for i in ions}
FR = []
IW = 11                                            # in-well half-offset: two ions side by side in one well
def snap(cap, hi=None, junc=None, badge="", merged=None):
    FR.append({"pos": {i: pos[i][:] for i in pos}, "cap": cap,
               "hi": hi or [], "junc": junc or [], "badge": badge,
               "merged": [sorted(p) for p in (merged or [])]})
def setp(i, x, y): pos[i] = [x, y]

# gates of a step, split into in-row and cross-row
def gates(step):
    inrow, cross = [], []
    for s in STABS:
        for (r, c) in s.data:
            if corner_step(s, (r, c)) != step:
                continue
            g = (ID_OF_STAB[s], ID_OF_DATA[(r, c)], CELL[s], r, c)
            (inrow if r == CELL[s] else cross).append(g)
    return inrow, cross

snap("Interleaved cells: data and ancillas in memory wells. Ancillas pumped to their check basis.")
snap("The X-check ancillas get a microwave pi/2 to |+>.",
     hi=[LABEL[s] for s in STABS if KIND[s] == "X"])

for step in range(D + 1) if False else range(4):
    L = step + 1
    inrow, cross = gates(step)
    # in-row gates: the ancilla steps INTO its data's well; the two sit side by
    # side in that one well while the gate fires, then both step back.
    if inrow:
        for a, d, _, _, _ in inrow:
            wx, wy = home[d][0], home[d][1]
            setp(d, wx - IW, wy); setp(a, wx + IW, wy)         # both ions inside the one well
        mp = [(a, d) for a, d, _, _, _ in inrow]
        snap(f"Step {L}: each in-row ancilla steps into its data's well.", hi=[a for a, *_ in inrow], merged=mp)
        snap(f"Step {L}: the two-qubit gate fires inside the well.", hi=[a for a, *_ in inrow], merged=mp)
        for a, d, _, _, _ in inrow:
            setp(d, home[d][0], home[d][1]); setp(a, home[a][0], home[a][1])
        snap(f"Step {L}: both step back to their own wells.", hi=[a for a, *_ in inrow])
    # cross-row gates: swap the ancilla into the target-column well, lift it
    # THROUGH the junction (in the gap), and land it IN the data's well of the
    # next cell, where the gate fires; then reverse every step.
    if cross:
        for a, d, ac, tr, c in cross:                          # swap ancilla into the target-column well
            dcol_id = ID_OF_DATA[(ac, c)]
            ax, dx = pos[a][0], pos[dcol_id][0]
            setp(a, dx, pos[a][1]); setp(dcol_id, ax, pos[a][1])
        snap(f"Step {L}: cross-row ancillas swap into the target-column well (well to well).",
             hi=[a for a, *_ in cross])
        for a, d, ac, tr, c in cross:                          # lift into the junction, in the gap
            setp(a, JX(c), gap_y(ac, tr))
        snap(f"Step {L}: they lift into the junction, which sits in the gap between wells.",
             hi=[a for a, *_ in cross], junc=[(c, min(ac, tr)) for a, d, ac, tr, c in cross],
             badge="in transit")
        for a, d, ac, tr, c in cross:                          # land IN the data's well of the next cell
            tid = ID_OF_DATA[(tr, c)]
            setp(tid, X(c) - IW, CY[tr]); setp(a, X(c) + IW, CY[tr])
        snap(f"Step {L}: each lands in its data's well; the gate fires there, on the junction side.",
             hi=[a for a, *_ in cross],
             merged=[(a, ID_OF_DATA[(tr, c)]) for a, d, ac, tr, c in cross])
        for a, d, ac, tr, c in cross:                          # data recentres, ancilla lifts back to the junction
            setp(ID_OF_DATA[(tr, c)], X(c), CY[tr]); setp(a, JX(c), gap_y(ac, tr))
        snap(f"Step {L}: gate done; the ancilla lifts back into the junction.",
             hi=[a for a, *_ in cross], junc=[(c, min(ac, tr)) for a, d, ac, tr, c in cross])
        for a, d, ac, tr, c in cross:                          # drop back into the home well, then swap back
            setp(a, X(c), CY[ac])
            dcol_id = ID_OF_DATA[(ac, c)]
            setp(a, home[a][0], home[a][1]); setp(dcol_id, X(c), CY[ac])
        snap(f"Step {L}: swap back; the chain order is restored.", hi=[a for a, *_ in cross])

# readout: bubble every ancilla to the SPAM (right) end of its cell
order = [[it for it in cell] for cell in CHAINS]
def bubbled(c):
    seen = False
    for k, _ in c:
        if k == "anc": seen = True
        elif seen: return False
    return True
layer = 0
while not all(bubbled(c) for c in order):
    for ci, c in enumerate(order):
        k = 0
        while k < len(c) - 1:
            if c[k][0] == "anc" and c[k + 1][0] == "data":
                c[k], c[k + 1] = c[k + 1], c[k]; k += 2
            else:
                k += 1
    # re-lay each cell left to right and snap
    for ci, c in enumerate(order):
        x = -0.5
        for kind, item in c:
            iid = ID_OF_STAB[item] if kind == "anc" else ID_OF_DATA[(item[0], item[1])] if False else None
            iid = ID_OF_STAB[item] if kind == "anc" else "d%d" % num((item), D)
            setp(iid, X(x), CY[ci]); x += 1
    layer += 1
    snap(f"Readout: ancillas bubble one step toward the SPAM end (adjacent swaps).",
         hi=[LABEL[s] for s in STABS], badge="swap-out")
snap("Every ancilla is contiguous at the SPAM end. 493 nm readout gives all 8 syndrome bits.",
     hi=[LABEL[s] for s in STABS], badge="8 syndromes")

# ----- independent overlap check (python side) ------------------------------
def bad_frame(f):
    mg = {frozenset(p) for p in f["merged"]}
    by = {}
    for i, (x, y) in f["pos"].items():
        by.setdefault(round(y), []).append((x, i))
    for y, arr in by.items():
        arr.sort()
        for k in range(1, len(arr)):
            if abs(arr[k][0] - arr[k - 1][0]) < 30 and \
               frozenset((arr[k][1], arr[k - 1][1])) not in mg:
                return f"{arr[k-1][1]} & {arr[k][1]} at y={y}"
    return None

# junction marks: OFFSET from the data columns (a routing lane between wells),
# in the gap between adjacent cells. Ions transit through them, never rest on them.
JUNC = [{"c": c, "b": b, "x": JX(c), "y1": CY[b], "y2": CY[b + 1]}
        for c in range(D) for b in range(D - 1)]

HTML = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>d=3 QEC round - balanced {3,3,2}</title><style>
:root{--bg:#faf9f5;--panel:#fff;--ink:#2c2c2a;--mut:#6b6a64;--line:#d9d7cd;
--x:#3B7FD4;--xd:#1C4F8C;--z:#D2703A;--zd:#8F4620;--purple:#534AB7;--teal:#0F6E56;--amber:#BA7517;--red:#A32D2D;--redbg:#FCEBEB;}
@media(prefers-color-scheme:dark){:root{--bg:#22211e;--panel:#2c2b27;--ink:#eceae2;--mut:#a3a299;--line:#45453f;--redbg:#3a1414;}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:-apple-system,system-ui,"Segoe UI",sans-serif}
.wrap{max-width:1000px;margin:0 auto;padding:16px}h1{font-size:18px;margin:2px 0}
.sub{font-size:13px;color:var(--mut);margin-bottom:12px}
.bar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:8px}
button{border:1px solid var(--line);background:var(--panel);color:var(--ink);border-radius:8px;padding:6px 12px;font-size:14px;cursor:pointer}
button:disabled{opacity:.4}button.play{background:var(--purple);color:#fff;border-color:var(--purple)}
#slider{flex:1;min-width:150px}.sname{font-size:13px;font-weight:600;min-width:64px}
.badge{font-size:12px;font-weight:600;padding:3px 10px;border-radius:12px;background:#E1F5EE;color:var(--teal)}
.verify{font-size:12px;font-weight:600;padding:3px 10px;border-radius:12px;background:#E1F5EE;color:var(--teal)}
.verify.bad{background:var(--redbg);color:var(--red)}
.cap{font-size:13.5px;color:var(--mut);min-height:40px;margin:4px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:6px}
svg{width:100%;height:auto;display:block}
.legend{display:flex;gap:16px;font-size:11.5px;color:var(--mut);margin-top:8px}
.legend i{display:inline-block;width:12px;height:12px;vertical-align:-1px;margin-right:5px}
.note{font-size:12px;color:var(--mut);border-top:1px solid var(--line);margin-top:12px;padding-top:10px;line-height:1.55}
</style></head><body><div class="wrap">
<h1>Distance-3 error-correction round &mdash; balanced {3,3,2} layout</h1>
<div class="sub">One code row per cell. Ions sit in wells; the junctions link neighbouring cells in the gate zone; gates happen off the junctions. A live verifier checks every frame for overlap.</div>
<div class="bar"><button id="prev">&lsaquo; Prev</button><button id="next">Next &rsaquo;</button>
<button id="play" class="play">&#9654; Play</button><input id="slider" type="range" min="0" max="0" value="0"/>
<span class="sname" id="sname"></span><span class="badge" id="badge"></span><span class="verify" id="verify"></span></div>
<div class="cap" id="cap"></div>
<div class="card"><svg id="stage" role="img" aria-label="QEC round"></svg></div>
<div class="legend"><span><i style="border-radius:50%;background:var(--panel);border:1px solid var(--line)"></i>Ba+ data</span>
<span><i style="border-radius:2px;background:var(--x)"></i>X-check ancilla</span>
<span><i style="border-radius:2px;background:var(--z)"></i>Z-check ancilla</span></div>
<div class="note"><b>Verifier.</b> Every frame is scanned along each cell axis; two non-merged ions closer than a spacing turn the badge red and name them, so green means nothing has collided or passed. <b>Junctions</b> (the crosses between cells) are gate-zone routing links, not wells: ions transit through them but gate at wells offset from them.</div>
</div><script>
const DATA=__DATA__;
const NS="http://www.w3.org/2000/svg",E=(t,a)=>{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;};
const $=i=>document.getElementById(i),svg=$("stage");
const W=DATA.xhi,H=DATA.celly[DATA.celly.length-1]+70;
svg.setAttribute("viewBox","0 0 "+W+" "+H);
const gZ=E("g",{}),gWell=E("g",{}),gJ=E("g",{}),gW=E("g",{}),gI=E("g",{});
svg.appendChild(gZ);svg.appendChild(gWell);svg.appendChild(gJ);svg.appendChild(gW);svg.appendChild(gI);
DATA.celly.forEach((y,ci)=>{
 gZ.appendChild(E("rect",{x:DATA.xlo,y:y-24,width:W-DATA.xlo-20,height:48,rx:10,fill:"none",stroke:"var(--line)","stroke-dasharray":"5 4",opacity:.55}));
 const t=E("text",{x:DATA.xlo+5,y:y-31,"font-size":11,fill:"var(--mut)"});t.textContent="cell "+ci;gZ.appendChild(t);});
(DATA.wells||[]).forEach(w=>{gWell.appendChild(E("rect",{x:w[0]-22,y:w[1]-19,width:44,height:38,rx:10,fill:"var(--panel)",stroke:"var(--line)","stroke-width":1.2}));});
const jel={};DATA.junctions.forEach(j=>{
 // the lane is the routing channel; the actual junctions are where it meets a cell.
 const lane=E("line",{x1:j.x,y1:j.y1+24,x2:j.x,y2:j.y2-24,stroke:"var(--line)","stroke-width":2.5,opacity:.5,"stroke-linecap":"round"});
 gJ.appendChild(lane);
 gJ.appendChild(E("circle",{cx:j.x,cy:j.y1+24,r:3.4,fill:"var(--mut)"}));   // junction: lane meets cell b
 gJ.appendChild(E("circle",{cx:j.x,cy:j.y2-24,r:3.4,fill:"var(--mut)"}));   // junction: lane meets cell b+1
 jel[j.c+"_"+j.b]=lane;});
const el={};for(const id in DATA.ions){const lab=DATA.ions[id][0],typ=DATA.ions[id][1],g=E("g",{});g.style.transition="transform .25s ease";
 if(typ==="data"){g.appendChild(E("circle",{cx:0,cy:0,r:15,fill:"var(--panel)",stroke:"var(--line)","stroke-width":1.4}));
  const t=E("text",{x:0,y:4,"text-anchor":"middle","font-size":10,fill:"var(--ink)"});t.textContent=lab;g.appendChild(t);el[id]={g};}
 else{const col=typ==="X"?"var(--x)":"var(--z)",dk=typ==="X"?"var(--xd)":"var(--zd)";
  const r=E("rect",{x:-15,y:-14,width:30,height:28,rx:5,fill:col,stroke:dk});g.appendChild(r);
  const t=E("text",{x:0,y:4,"text-anchor":"middle","font-size":10,fill:"#fff","font-weight":600});t.textContent=lab;g.appendChild(t);el[id]={g,rc:r,dk:dk};}
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
  if(el[id].rc){const on=f.hi.indexOf(id)>=0;el[id].rc.setAttribute("stroke-width",on?2.6:1);el[id].rc.setAttribute("stroke",on?"var(--amber)":el[id].dk);}}
 $("prev").disabled=step===0;$("next").disabled=step===DATA.frames.length-1;}
$("slider").max=DATA.frames.length-1;
$("next").onclick=()=>{if(step<DATA.frames.length-1){step++;render();}};
$("prev").onclick=()=>{if(step>0){step--;render();}};
$("slider").oninput=e=>{step=+e.target.value;render();};
let timer=null;$("play").onclick=()=>{if(timer){clearInterval(timer);timer=null;$("play").innerHTML="&#9654; Play";return;}
 $("play").innerHTML="&#10073;&#10073; Pause";timer=setInterval(()=>{if(step<DATA.frames.length-1){step++;render();}else{clearInterval(timer);timer=null;$("play").innerHTML="&#9654; Play";}},760);};
document.addEventListener("keydown",e=>{if(e.key==="ArrowRight")$("next").click();if(e.key==="ArrowLeft")$("prev").click();});
render();
</script></body></html>"""


def write_html(path):
    data = {"ions": {i: [ions[i]["lab"], ions[i]["type"]] for i in ions},
            "frames": FR, "celly": [CY[0], CY[1], CY[2]], "junctions": JUNC,
            "wells": [home[i] for i in ions],
            "xlo": X(-0.5) - 45, "xhi": X(D - 0.5) + 60}
    with open(path, "w") as f:
        f.write(HTML.replace("__DATA__", json.dumps(data)))


if __name__ == "__main__":
    bad = [(k, bad_frame(f)) for k, f in enumerate(FR) if bad_frame(f)]
    print(f"frames: {len(FR)}   overlaps: {len(bad)}")
    for k, b in bad[:8]:
        print(f"  frame {k}: {b}  [{FR[k]['cap'][:50]}]")
    write_html("qec_round_sim_d3.html")
    print("wrote qec_round_sim_d3.html")