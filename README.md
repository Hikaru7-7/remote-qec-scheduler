# remote-qec-scheduler

A scheduling feasibility certifier for a modular trapped-ion quantum computer.
It builds one full round of surface-code error correction, and a complete
remote lattice-surgery merge across two modules, as an ordered list of physical
operations on a segmented ion trap. Then it checks, at every code distance,
that the schedule is legal. Ions cannot pass in their one-dimensional channels,
a well holds one ion except the pair merged for a gate, junctions are exclusive,
and every gate fires in a well, never on a junction.

This is the software behind Section 4.4 of the thesis *Trapped-Ion
Multi-Computer Design for Remote Lattice Surgery* (Keio University, 2026).
Every count and every claim in that section can be reproduced here.

**What it is not.** It moves ions, not quantum states. It proves the round can
be scheduled without conflict. It returns no logical error rate.

## Requirements

Python 3. Nothing else. Both programs use only the standard library.

## Run it

```
python3 qec_scheduler.py        # sweep d = 3..15, all structural checks (expect 15 PASS)
python3 qec_scheduler.py 7      # one-distance report: placement, seam, depth, tally
python3 qec_visualizer.py      # build the d=3 and d=5 HTML animations
python3 qec_visualizer.py 7    # build the d=7 HTML animations
```

Checks pass at every odd distance tried up to `python3 qec_scheduler.py 27`.

## Files

| File | What it is |
|---|---|
| `qec_scheduler.py` | The source of truth. Builds the distance-d rotated surface code, places every ion on the chip, emits the schedule (`round_ops`), packs it into parallel time-steps (`parallel_steps`), and runs the structural checks. `op_tally` counts every physical beat so a duration model can turn the schedule into a round time. |
| `qec_visualizer.py` | Renders the scheduler's operation list, unchanged, as an interactive HTML animation. It assigns pixel coordinates and nothing else. It also runs its own overlap check on every frame, written independently of the scheduler's checks, so two separate programs agree the motion is legal. |
| `qec_round_sim_d{3,5,7}.html` | One local error-correction round at that distance. |
| `qec_merge_full_sim_d{3,5,7}.html` | The full remote lattice-surgery merge, d rounds with the seam read by communication ions. |

## Viewing the animations

Open any HTML file in a browser. No server needed. Use Prev / Next, the
slider, or Play. Each frame shows a caption of the physical operation, the
parallel time-step it belongs to, and a live verifier badge that turns red if
any two ions in that frame overlap. It never does. The badge is computed in
the page, on the drawn positions, so the reader can see the check pass rather
than take it on trust.

## Headline numbers to reproduce

- Local round depth is exactly `27 + d` parallel time-steps at every distance.
  A flat 24-step connection band, then a readout tail of `d + 3`. This is an
  asserted check, not a printout.
- The full d-round merge packs into 81, 180, and 280 time-steps at d = 3, 5, 7.
- Resources are closed forms. `d^2` data ions, `d^2 - 1` ancillas, `d` comm
  lanes, `d - 1` Bell pairs per merge round, `4d(d - 1)` two-qubit gates per
  round, `floor((d-1)/2)` park wells.

## Design contract

The scheduler owns the schedule. What happens, in what order, with what motion,
is decided in `qec_scheduler.py` and certified there. The visualizer draws that
list and adds nothing. It rejects any operation it does not have a renderer
for, and its frame-by-frame overlap test is an independent second check.

## Cite

```bibtex
@misc{Yokomori2026QECScheduler,
  author       = {Yokomori, Hikaru},
  title        = {Remote-{QEC} feasibility scheduler},
  year         = {2026},
  howpublished = {\url{https://github.com/Hikaru7-7/remote-qec-scheduler}},
  note         = {Open-source software}
}
```

## License

MIT. See `LICENSE`.
