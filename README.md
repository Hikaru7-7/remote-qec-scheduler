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
python3 qec_scheduler.py        # named checks, then every check at every odd d = 3..27 (expect 17 PASS)
python3 qec_scheduler.py 7      # one-distance report: placement, seam, depth, tally
python3 qec_visualizer.py      # build the d=3 and d=5 HTML animations
python3 qec_visualizer.py 7    # build the d=7 HTML animations
python3 qec_visualizer.py all  # rebuild every odd d = 3..27 and print the sweep table
```

The default run ends by certifying every odd distance up to 27, so the claim
in the thesis is exactly what the command shows. One distance at a time is
`python3 qec_scheduler.py 27`.

## Files

| File | What it is |
|---|---|
| `qec_scheduler.py` | The source of truth. Builds the distance-d rotated surface code, places every ion on the chip, emits the schedule (`round_ops`), packs it into parallel time-steps (`parallel_steps`), and runs the structural checks. `op_tally` counts every physical beat so a duration model can turn the schedule into a round time. |
| `qec_visualizer.py` | Renders the scheduler's operation list, unchanged, as an interactive HTML animation. It assigns pixel coordinates and nothing else. It also runs its own overlap check on every frame, written independently of the scheduler's checks, so two separate programs agree the motion is legal. |
| `qec_round_sim_d{3,5,7}.html` | One local error-correction round at that distance. |
| `qec_merge_full_sim_d{3,5,7}.html` | Two rounds of the remote lattice-surgery merge, the seam read by communication ions. The full d-round merge is packed and certified in the scheduler. |

## Viewing the animations

Watch them in the browser, nothing to download:

| | d = 3 | d = 5 | d = 7 |
|---|---|---|---|
| One local round | [round d3](https://hikaru7-7.github.io/remote-qec-scheduler/qec_round_sim_d3.html) | [round d5](https://hikaru7-7.github.io/remote-qec-scheduler/qec_round_sim_d5.html) | [round d7](https://hikaru7-7.github.io/remote-qec-scheduler/qec_round_sim_d7.html) |
| Full merge | [merge d3](https://hikaru7-7.github.io/remote-qec-scheduler/qec_merge_full_sim_d3.html) | [merge d5](https://hikaru7-7.github.io/remote-qec-scheduler/qec_merge_full_sim_d5.html) | [merge d7](https://hikaru7-7.github.io/remote-qec-scheduler/qec_merge_full_sim_d7.html) |

Or open any HTML file locally in a browser. No server needed. Use Prev / Next, the
slider, or Play. Each frame shows a caption of the physical operation, the
parallel time-step it belongs to, and a live verifier badge that turns red if
any two ions in that frame overlap. It never does. The badge is computed in
the page, on the drawn positions, so the reader can see the check pass rather
than take it on trust.

## The full sweep

The thesis claim covers every odd distance from 3 to 27, and the three hosted
pairs above are samples of it. One command rebuilds the animation pages for
every swept distance locally (about 116 MB of HTML, well under a minute) and
reprints the table below:

```
python3 qec_visualizer.py all
```

Each row pairs the scheduler's own checks with the visualizer's independent
frame check at one distance. Local round steps are exactly `22 + 2d` at every
distance, a 2-round merge costs exactly twice its per-round share, and no
frame at any distance contains an overlap. The d = 3 row runs one extra
check, the hand-verified reference comparison.

| d | scheduler checks | local round steps | frames | overlaps | 2-round merge steps | frames | overlaps |
|--:|:----------------:|------------------:|-------:|---------:|--------------------:|-------:|---------:|
| 3 | 16 PASS | 28 | 59 | 0 | 52 | 149 | 0 |
| 5 | 15 PASS | 32 | 71 | 0 | 64 | 177 | 0 |
| 7 | 15 PASS | 36 | 83 | 0 | 72 | 205 | 0 |
| 9 | 15 PASS | 40 | 95 | 0 | 80 | 233 | 0 |
| 11 | 15 PASS | 44 | 107 | 0 | 88 | 261 | 0 |
| 13 | 15 PASS | 48 | 119 | 0 | 96 | 289 | 0 |
| 15 | 15 PASS | 52 | 131 | 0 | 104 | 317 | 0 |
| 17 | 15 PASS | 56 | 143 | 0 | 112 | 345 | 0 |
| 19 | 15 PASS | 60 | 155 | 0 | 120 | 373 | 0 |
| 21 | 15 PASS | 64 | 167 | 0 | 128 | 401 | 0 |
| 23 | 15 PASS | 68 | 179 | 0 | 136 | 429 | 0 |
| 25 | 15 PASS | 72 | 191 | 0 | 144 | 457 | 0 |
| 27 | 15 PASS | 76 | 203 | 0 | 152 | 485 | 0 |

## Headline numbers to reproduce

- Local round depth is exactly `22 + 2d` parallel time-steps at every distance.
  A flat 19-step connection band, a readout tail of `d + 3`, then `d` swap
  layers back to the rest state. This is an asserted check, not a printout.
- A cross-row ancilla swaps toward its junction only when a data ion blocks
  the way; the junction-adjacent half of the crossings lift directly, and
  those free lifts overlap the in-row gates. That is why the band is 19.
- Every schedule ends in the placement it started from, every swap replayed
  and undone. This is its own check, so rounds and merges chain.
- The full d-round merge packs into 78, 160, and 252 time-steps at d = 3, 5, 7,
  and every one of its d rounds takes exactly the same number of steps.
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
