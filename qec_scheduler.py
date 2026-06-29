"""
Yokomori QEC round scheduler
============================
Lowers an abstract surface-code error-correction round onto a concrete
trapped-ion chip, in four layers (construct -> colour -> place -> schedule+check).
Guiding principle: GENERATE-then-CHECK. Each layer produces something, and
separate, simpler routines certify it. We never trust the generator, only the
checks plus the hand-verified d=3 base case.

LAYER 1 (this file): build the abstract code itself -- its data qubits and its
stabilizer checks -- and prove we built the right one.

  Arc:  REPRESENT  (Coord, Stabilizer)
        ENUMERATE  the qubits        (data_qubits)
        GENERATE   the checks        (build_stabilizers)   <- the heart
        AUDIT      two ways          (check_census, check_d3_matches_reference)
        RUN+REPORT                   (__main__)

Run it with:  python3 qec_scheduler.py
"""
from __future__ import annotations
from dataclasses import dataclass


# --- REPRESENT -------------------------------------------------------------
# How we locate a qubit: its (row, col) on the code's 2D grid. Which qubits a
# stabilizer touches is decided by grid adjacency, so coordinates are the
# foundation every later step builds on.
Coord = tuple[int, int]


@dataclass(frozen=True)
class Stabilizer:
    """One parity check -- the atomic object the whole scheduler revolves around.

    Records the only two facts that matter: its kind (X or Z) and the set of
    data qubits it measures. frozen=True makes it immutable and hashable, so a
    Stabilizer can live inside a set -- which the d=3 audit relies on.
    """
    kind: str                  # "X" or "Z"
    data: frozenset            # the data-qubit Coords this check touches

    @property
    def weight(self) -> int:   # 4 for a bulk check, 2 for a boundary check
        return len(self.data)


# --- ENUMERATE -------------------------------------------------------------
def data_qubits(d: int) -> list:
    """Lay down the board: the d*d data qubits a distance-d code has.
    The stabilizers are placed on top of these positions."""
    return [(r, c) for r in range(d) for c in range(d)]


def num(rc: Coord, d: int) -> int:
    """Translator between our (row, col) and the simulator's 1..d^2 labels.
    Exists only so the d=3 audit can compare against the hand-verified instance."""
    r, c = rc
    return r * d + c + 1


# --- GENERATE  (the heart of Layer 1) --------------------------------------
def build_stabilizers(d: int) -> list:
    """Encode the rotated surface code's rules and emit all d^2 - 1 checks.
    This is the 'construct the code' stage -- the one piece that actually
    produces the code the chip will later have to run.

    BULK (weight-4): for every r, c in 0..d-2, a plaquette on the 2x2 block
        {(r,c),(r,c+1),(r+1,c),(r+1,c+1)}, kind "X" if (r+c) even else "Z".
    BOUNDARY (weight-2): caps the dangling plaquettes; each edge carries one
        kind on alternating positions (top X, bottom X, left Z, right Z).

    Returns a plain list of Stabilizer objects; their order does not matter.
    """
    stabs = []

    # BULK -- a weight-4 plaquette on every interior 2x2 block of data qubits.
    # The X/Z kind alternates like a checkerboard, set by the parity of (r + c).
    for r in range(d - 1):
        for c in range(d - 1):
            block = frozenset({(r, c), (r, c + 1), (r + 1, c), (r + 1, c + 1)})
            kind = "X" if (r + c) % 2 == 0 else "Z"
            stabs.append(Stabilizer(kind, block))

    # BOUNDARY -- weight-2 caps on the dangling plaquettes, one loop per edge.
    # Top and bottom carry X, left and right carry Z, each on alternating
    # positions so the boundary interlocks with the bulk checkerboard.
    for c in range(1, d - 1, 2):                  # top edge,    odd columns
        stabs.append(Stabilizer("X", frozenset({(0, c), (0, c + 1)})))
    for c in range(0, d - 1, 2):                  # bottom edge, even columns
        stabs.append(Stabilizer("X", frozenset({(d - 1, c), (d - 1, c + 1)})))
    for r in range(0, d - 1, 2):                  # left edge,   even rows
        stabs.append(Stabilizer("Z", frozenset({(r, 0), (r + 1, 0)})))
    for r in range(1, d - 1, 2):                  # right edge,  odd rows
        stabs.append(Stabilizer("Z", frozenset({(r, d - 1), (r + 1, d - 1)})))

    return stabs


# --- AUDIT #1: cheap, and holds at EVERY distance --------------------------
def check_census(d: int) -> None:
    """Does the output have the known surface-code census? Catches dropped,
    extra, or mistyped checks. We don't take the generator's word for it --
    we audit its counts against an established fact."""
    stabs = build_stabilizers(d)
    nX   = sum(s.kind == "X" for s in stabs)
    nZ   = sum(s.kind == "Z" for s in stabs)
    bulk = sum(s.weight == 4 for s in stabs)
    bdry = sum(s.weight == 2 for s in stabs)
    assert len(stabs) == d*d - 1,  f"d={d}: {len(stabs)} stabilizers, want {d*d-1}"
    assert bulk == (d-1)**2,       f"d={d}: {bulk} bulk, want {(d-1)**2}"
    assert bdry == 2*(d-1),        f"d={d}: {bdry} boundary, want {2*(d-1)}"
    assert nX == nZ == (d*d-1)//2, f"d={d}: X={nX}, Z={nZ}, want {(d*d-1)//2} each"
    for s in stabs:
        assert s.weight in (2, 4), f"d={d}: weight-{s.weight} check found"
        for (r, c) in s.data:
            assert 0 <= r < d and 0 <= c < d, f"d={d}: ({r},{c}) out of grid"


# --- AUDIT #2: strong, but only where we have ground truth ------------------
# The eight d=3 checks, copied from the hand-verified HTML simulator (data
# labelled 1..9). Census proves the SHAPE is right; this proves the CONTENT is.
D3_REFERENCE = {
    ("X", frozenset({1, 2, 4, 5})), ("X", frozenset({5, 6, 8, 9})),
    ("X", frozenset({2, 3})),       ("X", frozenset({7, 8})),
    ("Z", frozenset({2, 3, 5, 6})), ("Z", frozenset({4, 5, 7, 8})),
    ("Z", frozenset({1, 4})),       ("Z", frozenset({6, 9})),
}


def check_d3_matches_reference() -> None:
    """Does d=3 reproduce, check for check, the instance we verified by hand?"""
    got = {(s.kind, frozenset(num(rc, 3) for rc in s.data))
           for s in build_stabilizers(3)}
    assert got == D3_REFERENCE, (
        "d=3 does not match the simulator:\n"
        f"  missing: {D3_REFERENCE - got}\n"
        f"  extra:   {got - D3_REFERENCE}"
    )


# --- RUN + REPORT ----------------------------------------------------------
# Execute the audits and print a one-line PASS/FAIL, so every change gets an
# instant verdict on whether Layer 1 is still correct.
if __name__ == "__main__":
    try:
        for d in (3, 5, 7, 9, 11):
            check_census(d)
        print("census check ........ PASS  (d = 3, 5, 7, 9, 11)")
        check_d3_matches_reference()
        print("d=3 vs simulator .... PASS")
        print("\nLayer 1 done. Next: Layer 2 -- the 4-step colouring.")
    except NotImplementedError as e:
        print("build_stabilizers not written yet:", e)
    except AssertionError as e:
        print("CHECK FAILED:\n", e)