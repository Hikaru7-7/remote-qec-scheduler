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
        for d in (3, 5, 7):
            counts = [sum(1 for k, _ in cell if k == "anc") for cell in place(d)]
            print(f"  d={d}: per-cell ancillas = {counts}  (max {max(counts)}, total {sum(counts)})")
        print("readout depth (swap layers to the SPAM end):")
        print(f"  d=3 balanced {{3,3,2}} : {readout_layers(place(3))} layers")
        print(f"  d=3 hand     {{4,2,2}} : {readout_layers(hand_422())} layers")
        for d in (5, 7):
            print(f"  d={d} balanced         : {readout_layers(place(d))} layers")
        print("\nLayer 4 in progress. Next: the connection-step moves (in-row + cross-row).")
    except NotImplementedError as e:
        print("not written yet:", e)
    except AssertionError as e:
        print("CHECK FAILED:\n", e)