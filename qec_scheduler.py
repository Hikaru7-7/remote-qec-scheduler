"""
Yokomori QEC round scheduler — Layer 1: CONSTRUCT
=================================================
Build the distance-d rotated surface code and check its census.

This is the first of four layers. The guiding principle of the whole
scheduler is GENERATE-then-CHECK: each layer produces something, and a
separate, simpler routine asserts it is correct. We never have to trust the
generator — only the (much simpler) checker, plus the hand-verified d=3 case.

YOUR JOB for Layer 1: implement `build_stabilizers(d)` from the spec in its
docstring. Everything else is scaffolding. Then run:

    python3 qec_scheduler.py

and make both checks print PASS.
"""
from __future__ import annotations
from dataclasses import dataclass

# ----------------------------------------------------------------------
# Data structures + plumbing  (provided — you don't need to touch these)
# ----------------------------------------------------------------------

Coord = tuple[int, int]            # (row, col); row 0 = top, col 0 = left


@dataclass(frozen=True)
class Stabilizer:
    kind: str                      # "X" or "Z"
    data: frozenset                # the data-qubit Coords this check touches

    @property
    def weight(self) -> int:
        return len(self.data)


def data_qubits(d: int) -> list:
    """All d*d data-qubit coordinates."""
    return [(r, c) for r in range(d) for c in range(d)]


def num(rc: Coord, d: int) -> int:
    """Row-major 1..d^2 label, matching d1..d9 in the HTML simulator."""
    r, c = rc
    return r * d + c + 1


# ----------------------------------------------------------------------
# YOUR CORE FUNCTION — write this one.
# ----------------------------------------------------------------------

def build_stabilizers(d: int) -> list:
    """
    Return all (d*d - 1) stabilizers of the distance-d rotated surface code.

    Coordinates: data qubit at (r, c), with r, c in 0..d-1.

    BULK (weight-4):  for every r, c in 0..d-2, a plaquette on the 2x2 block
        {(r, c), (r, c+1), (r+1, c), (r+1, c+1)},
        kind = "X" if (r + c) is even, else "Z".              -> (d-1)^2 checks

    BOUNDARY (weight-2): caps the dangling plaquettes on the four edges.
        Each edge carries ONE kind, on alternating positions.
        You are GIVEN the top edge as a worked example:
            TOP    (row 0):   pairs {(0, c), (0, c+1)},  kind "X",  for c ODD.
        Work out the other three by symmetry — the d=3 check below will tell
        you immediately whether your parities are right:
            BOTTOM (row d-1): pairs {(d-1, c), (d-1, c+1)}, kind "X", for c ???
            LEFT   (col 0):   pairs {(r, 0), (r+1, 0)},     kind "Z", for r ???
            RIGHT  (col d-1): pairs {(r, d-1), (r+1, d-1)},  kind "Z", for r ???
                                                              -> 2(d-1) checks

    (Stuck on the three "???" parities? Ask me and I'll reveal them — but try
     the d=3 check first; it's a fast guess-and-confirm.)
    """
    raise NotImplementedError("write me — see the docstring above")


# ----------------------------------------------------------------------
# Independent CHECKS  (provided — don't edit; just make them pass)
# ----------------------------------------------------------------------

def check_census(d: int) -> None:
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


# The hand-verified d=3 stabilizers, exactly as in qec_round_simulator.html
# (data labelled 1..9 row-major). This is our base-case ground truth.
D3_REFERENCE = {
    ("X", frozenset({1, 2, 4, 5})), ("X", frozenset({5, 6, 8, 9})),  # A1, A4 bulk
    ("X", frozenset({2, 3})),       ("X", frozenset({7, 8})),         # B1, B2 boundary
    ("Z", frozenset({2, 3, 5, 6})), ("Z", frozenset({4, 5, 7, 8})),  # A2, A3 bulk
    ("Z", frozenset({1, 4})),       ("Z", frozenset({6, 9})),         # B3, B4 boundary
}


def check_d3_matches_reference() -> None:
    got = {(s.kind, frozenset(num(rc, 3) for rc in s.data))
           for s in build_stabilizers(3)}
    assert got == D3_REFERENCE, (
        "d=3 does not match the simulator:\n"
        f"  you are MISSING: {D3_REFERENCE - got}\n"
        f"  you have EXTRA : {got - D3_REFERENCE}"
    )


if __name__ == "__main__":
    try:
        for d in (3, 5, 7, 9, 11):
            check_census(d)
        print("census check ........ PASS  (d = 3, 5, 7, 9, 11)")
        check_d3_matches_reference()
        print("d=3 vs simulator .... PASS")
        print("\nLayer 1 done. Next: Layer 2 — the 4-step colouring.")
    except NotImplementedError as e:
        print("build_stabilizers not written yet:", e)
    except AssertionError as e:
        print("CHECK FAILED:\n", e)
