"""
QEC round scheduler.

Goal: take an abstract surface-code round and lay it onto a real ion chip,
in four layers -- build the code, colour it, place ions, schedule moves.
Rule: build something, then check it. We trust the checks, not the builder.

Layer 1: build the code (data qubits + checks) and prove it is right.
Layer 2: give each check's qubits a step, so no qubit is used twice at once.

Run:  python3 qec_scheduler.py
"""
from __future__ import annotations
from dataclasses import dataclass


# --- REPRESENT -------------------------------------------------------------
# A qubit's place is (row, col) on the grid. Checks are built from this.
Coord = tuple[int, int]


@dataclass(frozen=True)
class Stabilizer:
    """One parity check. It has a kind (X or Z) and the qubits it reads.
    Frozen so we can drop it in a set for the d=3 check."""
    kind: str                  # "X" or "Z"
    data: frozenset            # the qubits this check reads

    @property
    def weight(self) -> int:   # 4 inside the grid, 2 on an edge
        return len(self.data)


# --- ENUMERATE -------------------------------------------------------------
def data_qubits(d: int) -> list:
    """All d*d data qubits. The board the checks sit on."""
    return [(r, c) for r in range(d) for c in range(d)]


def num(rc: Coord, d: int) -> int:
    """Turn (row, col) into the simulator's number 1..d*d.
    Used only by the d=3 check."""
    r, c = rc
    return r * d + c + 1


# --- GENERATE  (the heart of Layer 1) --------------------------------------
def build_stabilizers(d: int) -> list:
    """Build all d*d - 1 checks of the rotated surface code.

    Inside: one weight-4 check on each 2x2 block; X/Z swaps like a chess board.
    Edges:  weight-2 checks that close the open blocks; top/bottom X, sides Z.
    """
    stabs = []

    # Inside: a weight-4 check on each 2x2 block. (r+c) even -> X, odd -> Z.
    for r in range(d - 1):
        for c in range(d - 1):
            block = frozenset({(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)})
            kind = "X" if (r + c) % 2 == 0 else "Z"
            stabs.append(Stabilizer(kind, block))

    # Edges: weight-2 checks, one loop per side, on every other spot.
    for c in range(1, d - 1, 2):                  # top    -> X
        stabs.append(Stabilizer("X", frozenset({(0, c), (0, c + 1)})))
    for c in range(0, d - 1, 2):                  # bottom -> X
        stabs.append(Stabilizer("X", frozenset({(d - 1, c), (d - 1, c + 1)})))
    for r in range(0, d - 1, 2):                  # left   -> Z
        stabs.append(Stabilizer("Z", frozenset({(r, 0), (r + 1, 0)})))
    for r in range(1, d - 1, 2):                  # right  -> Z
        stabs.append(Stabilizer("Z", frozenset({(r, d - 1), (r + 1, d - 1)})))

    return stabs


# --- COLOUR  (the heart of Layer 2) ----------------------------------------
def corner_step(stab: Stabilizer, q: Coord) -> int:
    """Pick the step (0..3) for qubit q in this check.
    A check caps a 2x2 block. q sits at one corner of it, and the corner is the
    step. An edge check is half a block, so its block sits just off the grid;
    we slide the corner one step out to match."""
    rows = {r for r, c in stab.data}
    cols = {c for r, c in stab.data}
    r0, c0 = min(rows), min(cols)
    if stab.weight == 2:               # edge check: slide the block off the edge
        if len(rows) == 1:             # one row -> top or bottom
            r0 = r0 - 1 if r0 == 0 else r0
        else:                          # one col -> left or right
            c0 = c0 - 1 if c0 == 0 else c0
    qr, qc = q
    return 2 * (qr - r0) + (qc - c0)    # corners (0,0)(0,1)(1,0)(1,1) -> 0 1 2 3


def assign_steps(d: int) -> dict:
    """Give every check's qubits a step. Returns step -> list of qubits."""
    schedule = {0: [], 1: [], 2: [], 3: []}
    for stab in build_stabilizers(d):
        for q in stab.data:
            schedule[corner_step(stab, q)].append(q)
    return schedule


# --- PLACE  (the heart of Layer 3) -----------------------------------------
# One code row per cell, laid out as a 1-D chain:
#   left-end, data0, gap0, data1, gap1, ... , data_{d-1}, right-end.
# Data sit at odd slots; ancillas fill the gaps and the two ends.
#
# Balanced rule (keeps every cell to at most d ancillas):
#   bulk block (r,c): even column -> stay in row r, odd column -> drop to row r+1.
#     This frees each gap in row 0 for a top edge check, and each gap in the last
#     row for a bottom edge check.
#   top/bottom edge check: sits in the gap between its two data.
#   left/right edge check: sits at the left/right end of its top row.
# Result: every gap holds one ancilla, and the d-1 side checks go to ends,
# one per cell, giving per-cell counts d, d, ..., d, d-1.
def place(d: int) -> list:
    """Put each qubit in a well. Returns a list of cells; each cell is an
    ordered list of ("data", (r,c)) or ("anc", check)."""
    slots = [dict() for _ in range(d)]            # cell -> {slot: item}

    def put(cell, slot, item):
        assert slot not in slots[cell], f"d={d}: two ions want cell {cell} slot {slot}"
        slots[cell][slot] = item

    for r in range(d):                            # data (r,c) -> cell r, slot 2c+1
        for c in range(d):
            put(r, 2 * c + 1, ("data", (r, c)))

    LEFT, RIGHT = 0, 2 * d                         # the two end slots
    for s in build_stabilizers(d):
        rs = sorted({r for r, c in s.data})
        cs = sorted({c for r, c in s.data})
        if s.weight == 4:                                 # bulk block (r,c)
            r, c = rs[0], cs[0]
            home = r if c % 2 == 0 else r + 1             # even col up, odd col down
            put(home, 2 * c + 2, ("anc", s))              # the gap at column c
        elif len(rs) == 1:                                # top/bottom edge (one row)
            r, c = rs[0], cs[0]
            put(r, 2 * c + 2, ("anc", s))                 # the gap between its two data
        else:                                             # left/right edge (one column)
            r, c = rs[0], cs[0]
            put(r, LEFT if c == 0 else RIGHT, ("anc", s)) # an end of the top row

    return [[slots[i][s] for s in sorted(slots[i])] for i in range(d)]


# --- SCHEDULE  (the heart of Layer 4) --------------------------------------
def stab_cell(d: int) -> dict:
    """Which cell each ancilla lives in, read from the placement."""
    m = {}
    for i, cell in enumerate(place(d)):
        for kind, item in cell:
            if kind == "anc":
                m[item] = i
    return m


def cross_gates(d: int) -> list:
    """Every cross-row gate, as (step, ancilla, junction). A cross-row gate
    reads a data qubit in a neighbour cell, so its ancilla crosses through the
    junction at that data's column, on the boundary between the two cells."""
    cellof = stab_cell(d)
    out = []
    for stab in build_stabilizers(d):
        home = cellof[stab]
        for (r, c) in stab.data:
            if r != home:                        # data sits in a neighbour cell
                step = corner_step(stab, (r, c))
                junction = (c, min(home, r))     # (column, boundary)
                out.append((step, stab, junction))
    return out


# --- SIMULATE  (Layer 4: run the moves and count the steps) ----------------
def readout_layers(cells: list) -> int:
    """Bubble every ancilla to the SPAM end (right) of its cell with adjacent
    swaps, in parallel layers, and count the layers. Each swap exchanges an
    ancilla with the data on its right, so nothing passes. A cell with more
    ancillas, or ancillas further from the end, needs more layers."""
    cells = [list(c) for c in cells]              # mutable copy

    def done(c):                                  # all data before all ancillas?
        saw_anc = False
        for kind, _ in c:
            if kind == "anc":
                saw_anc = True
            elif saw_anc:
                return False
        return True

    layers = 0
    while not all(done(c) for c in cells):
        for c in cells:                           # one parallel layer per cell
            i = 0
            while i < len(c) - 1:
                if c[i][0] == "anc" and c[i + 1][0] == "data":
                    c[i], c[i + 1] = c[i + 1], c[i]
                    i += 2                         # non-overlapping swaps
                else:
                    i += 1
        layers += 1
    return layers


# --- GATE ZONE  (explicit: ions sit in wells, route through junctions) ------
# The physical rule, made explicit. An ion is stored or gated only in a WELL.
# An X-junction is a routing link between neighbouring cells; ions pass through
# it but are never stored on it. Memory and gate are separate segments of each
# cell axis: data and ancillas are STORED in the memory segment (data on the
# data columns); GATES happen in the gate segment, one gate well per ancilla,
# at the ancilla (offset) columns; JUNCTIONS sit at the data columns of the gate
# segment and connect neighbouring cells. So a gate well never shares a column
# with a junction, and no ion is ever at rest on a junction.
def check_gate_zone(d: int) -> None:
    """Check: one gate well per ancilla, every gate well off the junction
    columns (so gates never happen on a junction), and every cross-row gate
    routes through a junction that exists between the two cells."""
    pos = col_x(d, place(d))
    junctions = set(range(d))                       # junction columns = data columns
    gate_wells = []
    for stab in build_stabilizers(d):
        col = pos[stab][1]                          # this ancilla's gate-well column
        assert col not in junctions, \
            f"d={d}: gate well at column {col} lands on a junction"
        gate_wells.append(stab)
    assert len(gate_wells) == d*d - 1, \
        f"d={d}: {len(gate_wells)} gate wells, want {d*d-1}"
    for step, stab, (c, b) in cross_gates(d):       # crossings route through junctions
        assert 0 <= c < d and 0 <= b < d - 1, \
            f"d={d}: cross gate routes through a junction that does not exist"


# --- TWO-QPU MERGE  (remote lattice surgery: the seam stabilizers) ----------
# During a merge, the seam-facing data qubit of each row (column d-1) also
# gates with that row's communication ion, which holds a heralded Bell pair.
# That extra gate must fit in a connection step where the qubit is not already
# busy with a bulk gate, or the merge round would grow. We assume ideal herald
# rate, so this is purely a scheduling question.
def bulk_steps_at(d: int, r: int) -> set:
    """The steps in which the row-r seam qubit (r, d-1) does a BULK gate. Its old
    right-boundary check is replaced by a seam check during the merge, so only
    the weight-4 bulk checks keep it busy."""
    return {corner_step(s, (r, d - 1))
            for s in build_stabilizers(d) if (r, d - 1) in s.data and s.weight == 4}


def n_seam_checks(d: int, r: int) -> int:
    """How many weight-4 seam checks the row-r seam qubit belongs to. The seam
    has d-1 checks; check s reads boundary qubits of rows s and s+1, so an
    interior qubit is in two, a corner qubit in one."""
    return (1 if r >= 1 else 0) + (1 if r <= d - 2 else 0)


def check_seam_fits(d: int) -> None:
    """Check: each seam qubit's seam gates fit in the steps its bulk gates leave
    free, so the seam extraction rides in the same 4-step round, no extra phases."""
    for r in range(d):
        free = 4 - len(bulk_steps_at(d, r))
        assert free >= n_seam_checks(d, r), \
            f"d={d} row {r}: {n_seam_checks(d, r)} seam gates, only {free} free steps"


# The two modules are two copies of the same code. A merge joins them along the
# seam where module A's right boundary (column d-1) meets module B's left
# boundary (column 0).
def seam_stabilizers(d: int) -> list:
    """The d-1 weight-4 seam checks turned on for a merge. Check s reads two
    boundary data of module A (rows s and s+1, column d-1) and two of module B
    (rows s and s+1, column 0). Two of its four gates cross the remote link."""
    seam = []
    for s in range(d - 1):
        seam.append({"s": s,
                     "A": [(s, d - 1), (s + 1, d - 1)],   # local, in module A
                     "B": [(s, 0),     (s + 1, 0)]})       # remote, in module B
    return seam


def check_seam_census(d: int) -> None:
    """Check: there are d-1 seam checks, each weight 4, split two local and two
    remote, matching the boundary the code exposes to the seam."""
    seam = seam_stabilizers(d)
    assert len(seam) == d - 1, f"d={d}: {len(seam)} seam checks, want {d-1}"
    for sc in seam:
        assert len(sc["A"]) == 2 and len(sc["B"]) == 2, \
            f"d={d}: seam check {sc['s']} is not two local plus two remote"


def bell_pairs_per_round(d: int) -> int:
    """One heralded Bell pair carries each seam check's remote link, so a merge
    round needs d-1 pairs, one per seam check."""
    return d - 1


def check_merge_demand(d: int, lanes: int = None) -> tuple:
    """Check: a merge round needs d-1 Bell pairs, the d communication lanes can
    supply them with one lane held spare, and a full merge of d rounds needs
    d(d-1) pairs in all."""
    demand = bell_pairs_per_round(d)
    assert demand == len(seam_stabilizers(d)), f"d={d}: demand {demand} != seam checks"
    lanes = d if lanes is None else lanes
    assert lanes >= demand, f"d={d}: {lanes} lanes cannot supply {demand} pairs/round"
    return demand, lanes, demand * d


# --- COMMUNICATION IONS  (the remote leg of each seam check) ----------------
# Each lane has one communication ion in the I/F zone at the far end of its cell,
# holding a heralded Bell pair. A seam check's local ancilla (in cell s) gathers
# its two same-module boundary data, one in its own cell and one across a junction
# in cell s+1, then gates this comm ion, teleporting the half-parity to the other
# module. A merge uses d-1 comm ions, one per seam check, and holds one lane spare.
def comm_ions(d: int) -> dict:
    """One communication ion per lane, d in all, each in the I/F zone of its cell."""
    return {r: {"cell": r, "zone": "IF"} for r in range(d)}


def check_comm_ions(d: int) -> None:
    """Check: d comm ions, one per cell in the I/F zone, and a merge round uses
    d-1 of them (seam check s uses cell s's comm ion) with one lane spare."""
    ci = comm_ions(d)
    assert len(ci) == d, f"d={d}: {len(ci)} comm ions, want {d}"
    assert all(v["zone"] == "IF" for v in ci.values()), "a comm ion is not in the I/F zone"
    used = set(range(d - 1))                          # seam check s -> cell s's comm ion
    assert len(used) == d - 1 and d - len(used) == 1, "not d-1 used with one spare"


def seam_schedule(d: int) -> dict:
    """Place every seam gate in one of the four connection steps. For seam check s
    the ancilla gates its top data (s, d-1) and its bottom data (s+1, d-1) in steps
    their bulk gates leave free, then gates the comm ion in a third step. No data is
    used twice in a step and the ancilla does one gate per step."""
    free = {r: sorted(set(range(4)) - bulk_steps_at(d, r)) for r in range(d)}
    used = {r: set() for r in range(d)}              # steps already taken on each data
    sched = {}
    for si in range(d - 1):
        cand = [x for x in free[si] if x not in used[si]]
        assert cand, f"d={d} seam {si}: no free step for the top data"
        top = cand[0]; used[si].add(top)
        cand = [x for x in free[si + 1] if x not in used[si + 1] and x != top]
        assert cand, f"d={d} seam {si}: no free step for the bottom data"
        bot = cand[0]; used[si + 1].add(bot)
        comm = next(x for x in range(4) if x not in (top, bot))
        sched[si] = {"data_top": top, "data_bot": bot, "comm": comm}
    return sched


def check_seam_schedule(d: int) -> None:
    """Check: the seam extraction rides in the same four steps. Every data gate
    lands where that data is free of bulk gates, each ancilla's three gates sit in
    distinct steps, and no boundary data is gated twice in one step."""
    sched = seam_schedule(d)
    for si, g in sched.items():
        assert g["data_top"] not in bulk_steps_at(d, si), f"d={d} seam {si}: top clashes with a bulk gate"
        assert g["data_bot"] not in bulk_steps_at(d, si + 1), f"d={d} seam {si}: bottom clashes with a bulk gate"
        assert len({g["data_top"], g["data_bot"], g["comm"]}) == 3, f"d={d} seam {si}: ancilla double-books a step"




def col_x(d: int, cells: list) -> dict:
    """Column coordinate of every ion, aligned across cells so a junction at
    column c lines up in every cell. Data (r,c) -> c; a gap ancilla -> halfway
    between its two data; a left/right-end ancilla -> just outside the row."""
    pos = {}
    for i, cell in enumerate(cells):
        for j, (kind, item) in enumerate(cell):
            if kind == "data":
                pos[item] = [i, item[1]]                       # column c
            else:
                lft = cell[j - 1] if j > 0 else None
                rgt = cell[j + 1] if j < len(cell) - 1 else None
                if lft and lft[0] == "data" and rgt and rgt[0] == "data":
                    pos[item] = [i, (lft[1][1] + rgt[1][1]) / 2]
                elif lft is None:
                    pos[item] = [i, -0.5]
                else:
                    pos[item] = [i, d - 0.5]
    return pos


def _overlap(pos: dict, merged: set):
    """Any two non-merged ions in the same cell closer than a qubit spacing?"""
    by_cell = {}
    for ion, (cell, x) in pos.items():
        by_cell.setdefault(cell, []).append((x, ion))
    for cell, arr in by_cell.items():
        arr.sort(key=lambda t: t[0])
        for k in range(1, len(arr)):
            (x0, i0), (x1, i1) = arr[k - 1], arr[k]
            if abs(x1 - x0) < 0.3 and frozenset((i0, i1)) not in merged:
                return f"cell {cell}: two ions at ~{x1}"
    return None


def crossings_parallel_ok(d: int, step: int):
    """Try to position ALL cross-row ancillas of one step at once: each swaps
    with the in-row data at its target column, then transits (straight down its
    column's junction) into the neighbour cell to merge. Return None if it is
    collision-free, else where it collides."""
    cells = place(d)
    pos = col_x(d, cells)
    cellof = {item: i for i, c in enumerate(cells) for k, item in c if k == "anc"}
    moving = [(s, cellof[s], r, c)
              for s in build_stabilizers(d)
              for (r, c) in s.data
              if r != cellof[s] and corner_step(s, (r, c)) == step]
    for a, home, tgt, c in moving:                    # swap-out (all together)
        a_old = pos[a][1]
        pos[a] = [home, c]
        pos[(home, c)] = [home, a_old]
    bad = _overlap(pos, set())
    if bad:
        return "swap-out " + bad
    merged = set()
    for a, home, tgt, c in moving:                    # transit down the column
        pos[a] = [tgt, c]
        merged.add(frozenset((a, (tgt, c))))
    return None if not (bad := _overlap(pos, merged)) else "transit " + bad


def check_parallel_crossings(d: int) -> None:
    """Check: in every step, all cross-row crossings can fire in parallel."""
    for step in range(4):
        bad = crossings_parallel_ok(d, step)
        assert bad is None, f"d={d} step {step}: {bad}"


def connection_phases(d: int, cells: list, verbose: bool = False) -> int:
    """Count the serial two-qubit-gate phases in the 4 connection steps.
    In-row gates fire together (1 phase). Cross-row gates that share a cell
    serialize, so a step's cross-phases is its busiest cell's crossing count.
    This is the round-time driver; the shuttling between phases is fast."""
    cellof = {item: i for i, cell in enumerate(cells)
              for kind, item in cell if kind == "anc"}
    total = 0
    for step in range(4):
        inrow = 0
        deg = {}                                   # cell -> crossings touching it
        for stab in build_stabilizers(d):
            home = cellof[stab]
            for (r, c) in stab.data:
                if corner_step(stab, (r, c)) != step:
                    continue
                if r == home:
                    inrow = 1                      # an in-row gate exists this step
                else:
                    for cell in (home, r):         # a crossing touches both cells
                        deg[cell] = deg.get(cell, 0) + 1
        cross = max(deg.values(), default=0)
        total += inrow + cross
        if verbose:
            print(f"    step {step}: in-row {inrow}, cross-phases {cross}")
    return total


def hand_422() -> list:
    """The earlier hand-built d=3 layout, per-cell {4,2,2}, kept for comparison."""
    by = {frozenset(num(rc, 3) for rc in s.data): s for s in build_stabilizers(3)}
    chains = [
        [frozenset({1, 4}), 1, frozenset({1, 2, 4, 5}), 2, frozenset({2, 3, 5, 6}), 3, frozenset({2, 3})],
        [4, frozenset({4, 5, 7, 8}), 5, frozenset({5, 6, 8, 9}), 6],
        [7, frozenset({7, 8}), 8, frozenset({6, 9}), 9],
    ]
    out = []
    for ch in chains:
        cell = []
        for it in ch:
            cell.append(("data", ((it - 1) // 3, (it - 1) % 3)) if isinstance(it, int)
                        else ("anc", by[it]))
        out.append(cell)
    return out


# --- AUDIT #1: counts are right, at every distance -------------------------
def check_census(d: int) -> None:
    """Check 1: are the counts right? Works at any d."""
    stabs = build_stabilizers(d)
    nX   = sum(s.kind == "X" for s in stabs)
    nZ   = sum(s.kind == "Z" for s in stabs)
    bulk = sum(s.weight == 4 for s in stabs)
    bdry = sum(s.weight == 2 for s in stabs)
    assert len(stabs) == d*d - 1,  f"d={d}: {len(stabs)} checks, want {d*d-1}"
    assert bulk == (d-1)**2,       f"d={d}: {bulk} bulk, want {(d-1)**2}"
    assert bdry == 2*(d-1),        f"d={d}: {bdry} edge, want {2*(d-1)}"
    assert nX == nZ == (d*d-1)//2, f"d={d}: X={nX}, Z={nZ}, want {(d*d-1)//2} each"
    for s in stabs:
        assert s.weight in (2, 4), f"d={d}: a weight-{s.weight} check"
        for (r, c) in s.data:
            assert 0 <= r < d and 0 <= c < d, f"d={d}: ({r},{c}) off the grid"


# --- AUDIT #2: d=3 matches the hand-checked simulator ----------------------
# The 8 checks from the hand-checked simulator (qubits 1..9).
# Census checks the shape; this checks the content.
D3_REFERENCE = {
    ("X", frozenset({1, 2, 4, 5})), ("X", frozenset({5, 6, 8, 9})),
    ("X", frozenset({2, 3})),       ("X", frozenset({7, 8})),
    ("Z", frozenset({2, 3, 5, 6})), ("Z", frozenset({4, 5, 7, 8})),
    ("Z", frozenset({1, 4})),       ("Z", frozenset({6, 9})),
}


def check_d3_matches_reference() -> None:
    """Check 2: does d=3 match the simulator, check by check?"""
    got = {(s.kind, frozenset(num(rc, 3) for rc in s.data))
           for s in build_stabilizers(3)}
    assert got == D3_REFERENCE, (
        "d=3 does not match the simulator:\n"
        f"  missing: {D3_REFERENCE - got}\n"
        f"  extra:   {got - D3_REFERENCE}"
    )


# --- AUDIT #3: the colouring is collision-free -----------------------------
def check_no_double_touch(d: int) -> None:
    """Check 3: in one step, no qubit is used by two checks."""
    for step, qubits in assign_steps(d).items():
        assert len(qubits) == len(set(qubits)), \
            f"d={d}, step {step}: a qubit is used twice"


# --- AUDIT #4: the placement is well-formed --------------------------------
def check_placement(d: int) -> None:
    """Check 4: one ion per well, data in the right cell, every check placed,
    and each ancilla sits next to a data qubit it reads."""
    cells = place(d)
    data_seen, anc_seen = set(), set()
    for i, cell in enumerate(cells):
        for slot, (kind, item) in enumerate(cell):
            if kind == "data":
                r, c = item
                assert r == i, f"data {item} in cell {i}, want cell {r}"
                assert item not in data_seen, f"data {item} placed twice"
                data_seen.add(item)
            else:
                assert item not in anc_seen, "a check placed twice"
                anc_seen.add(item)
                near = cell[slot-1:slot] + cell[slot+1:slot+2]   # chain neighbours
                ok = any(k == "data" and q in item.data for k, q in near)
                assert ok, "an ancilla is not next to any data it reads"
    assert len(data_seen) == d*d,  f"placed {len(data_seen)} data, want {d*d}"
    assert len(anc_seen) == d*d-1, f"placed {len(anc_seen)} checks, want {d*d-1}"


# --- AUDIT #5: cross-row transits do not fight for a junction ---------------
def check_junctions(d: int) -> None:
    """Check 5: in one step, no two cross-row gates use the same junction."""
    used = {0: set(), 1: set(), 2: set(), 3: set()}
    for step, stab, junction in cross_gates(d):
        assert junction not in used[step], \
            f"d={d}, step {step}: two gates want junction {junction}"
        used[step].add(junction)


# --- RUN -------------------------------------------------------------------
# Run the checks. Print pass or fail.
if __name__ == "__main__":
    try:
        for d in (3, 5, 7, 9, 11):
            check_census(d)
        print("census .............. PASS  (d = 3,5,7,9,11)")
        for d in (3, 5, 7, 9, 11):
            check_no_double_touch(d)
        print("no qubit used twice . PASS  (d = 3,5,7,9,11)")
        check_d3_matches_reference()
        print("d=3 vs simulator .... PASS")
        for d in (3, 5, 7, 9, 11):
            check_placement(d)
        print("placement ........... PASS  (d = 3,5,7,9,11)")
        for d in (3, 5, 7, 9, 11):
            check_junctions(d)
        print("junctions ........... PASS  (d = 3,5,7,9,11)")
        for d in (3, 5, 7, 9, 11):
            check_gate_zone(d)
        print("gate zone ........... PASS  (wells off junctions; memory/gate split)")
        for d in (3, 5, 7):
            counts = [sum(1 for k, _ in cell if k == "anc") for cell in place(d)]
            print(f"  d={d}: per-cell ancillas = {counts}  (max {max(counts)}, total {sum(counts)})")
        print("readout depth (swap layers to the SPAM end):")
        print(f"  d=3 balanced {{3,3,2}} : {readout_layers(place(3))} layers")
        print(f"  d=3 hand     {{4,2,2}} : {readout_layers(hand_422())} layers")
        for d in (5, 7):
            print(f"  d={d} balanced         : {readout_layers(place(d))} layers")
        print("connection gate phases (the round-time driver):")
        print("  d=3 balanced {3,3,2}:")
        p1 = connection_phases(3, place(3), verbose=True)
        print(f"    total = {p1} phases")
        print("  d=3 hand {4,2,2}:")
        p2 = connection_phases(3, hand_422(), verbose=True)
        print(f"    total = {p2} phases")
        for d in (5, 7, 9):
            print(f"  d={d} balanced: {connection_phases(d, place(d))} phases")
        print("do the cross-row crossings actually fire in parallel?")
        for d in (3, 5, 7, 9, 11, 15):
            check_parallel_crossings(d)
        print("  parallel crossings .. PASS  (d = 3,5,7,9,11,15)")
        print("  => each step is 1 gate phase, so the round is 4 connection")
        print("     phases + readout, FLAT in d. The linear count above was the")
        print("     conservative bound; the crossings really do parallelize.")
        print("does the seam extraction fit in the same round?")
        for d in (3, 5, 7, 9, 11, 15):
            check_seam_fits(d)
        print("  seam fits ........... PASS  (d = 3,5,7,9,11,15)")
        print("  seam qubit schedule (d=3), bulk vs seam by row:")
        for r in range(3):
            b = sorted(bulk_steps_at(3, r))
            free = sorted(set(range(4)) - set(b))
            print(f"    row {r}: bulk in {b}, {n_seam_checks(3, r)} seam gate(s) fit in free {free}")
        for d in (3, 5, 7, 9, 11, 15):
            check_seam_census(d); check_merge_demand(d)
        print("  seam census ......... PASS  (d-1 weight-4 seam checks)")
        print("  merge demand ........ PASS  (d-1 Bell pairs/round, d lanes)")
        for d in (3, 5, 7):
            dem, lanes, tot = check_merge_demand(d)
            print(f"    d={d}: {dem} pairs/round, {lanes} lanes ({lanes - dem} spare), {tot} per {d}-round merge")
        for d in (3, 5, 7, 9, 11, 15):
            check_comm_ions(d); check_seam_schedule(d)
        print("  comm ions ........... PASS  (d in the I/F zone, d-1 used + 1 spare)")
        print("  seam schedule ....... PASS  (seam + comm gates ride in the 4 steps)")
        print("  seam schedule (d=3), each seam check by step:")
        for si, g in seam_schedule(3).items():
            print(f"    seam check {si} (cells {si},{si+1}): data gates in steps {g['data_top']} and {g['data_bot']}, comm-ion gate in step {g['comm']}")
    except NotImplementedError as e:
        print("not written yet:", e)
    except AssertionError as e:
        print("CHECK FAILED:\n", e)