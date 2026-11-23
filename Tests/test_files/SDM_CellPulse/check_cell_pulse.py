"""Check step-start SD births and integer doses using native AMReX particle records."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

# The neighboring SDM tests live under test_files, not on the installed Python path.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from SDM_Box3D_SampledMultiplicity.check_sampled_multiplicity import particle_records, read_inputs


INJ = "super_droplets_moisture.injection."


class Oracle:
    def __init__(self, deck):
        inputs = read_inputs(deck)
        self.dt = Fraction(inputs["erf.fixed_dt"][0])
        self.last = int(inputs["max_step"][0])
        self.period = Fraction(inputs[INJ + "period_s"][0])
        self.start = Fraction(inputs[INJ + "t_start"][0])
        self.stop = Fraction(inputs[INJ + "t_stop"][0])
        shape = tuple(map(int, inputs["amr.n_cell"]))
        self.cells = set(product(*(range(n) for n in shape)))
        self.lo = tuple(map(Fraction, inputs["geometry.prob_lo"]))
        extent = [Fraction(hi) - lo for hi, lo in zip(inputs["geometry.prob_hi"], self.lo)]
        self.dx = [length / n for length, n in zip(extent, shape)]
        # Round only the cumulative *domain* dose, not the dose per event/cell/box/rank.
        self.dose = Fraction(inputs[INJ + "rate"][0]) * self.period * math.prod(extent)
        self.low = math.floor(self.dose / len(self.cells))

    def events(self, time):
        # A state AT an event time still precedes the advance which injects it.
        return max(0, math.ceil((min(time, self.stop) - self.start) / self.period))

    def represented(self, events):
        return math.floor(events * self.dose + Fraction(1, 2))

    def snapshot(self, directory, step):
        plot = directory / f"plt{step:05d}"
        lines = (plot / "Header").read_text().splitlines()
        n = int(lines[1])
        assert int(lines[n+2]) == 3 and int(lines[n+4]) == 0, (plot, "not a single 3-D level")
        assert int(lines[n+9]) == step and Fraction(lines[n+3]) == step * self.dt, (
            plot, "wrong state step/time", lines[n+9], lines[n+3])
        records = {}
        for level, ident, active, pos, (weight,) in particle_records(plot, ["multiplicity"]):
            assert level == 0 and active > 0 and ident[0] > 0, (plot, level, ident, active)
            assert ident not in records, (plot, "duplicate particle id/cpu", ident)
            cell = tuple(math.floor((x - lo) / dx) for x, lo, dx in zip(pos, self.lo, self.dx))
            assert cell in self.cells, (plot, "outside domain", ident, pos)
            assert all(x == lo + (i + 0.5) * dx for x, lo, i, dx in zip(pos, self.lo, cell, self.dx)), (
                plot, "not at its cell center", ident, pos, cell)
            assert weight in (self.low, self.low + 1), (plot, ident, "unexpected multiplicity", weight)
            records[ident] = (cell, int(weight))
        return records

    def history(self, directory, ranks, first=0, previous=None, reference=None, last=None):
        """Track native id/cpu within a run. Compare cell weights, not newly assigned IDs, across runs."""
        states = {}
        for step in range(first, (self.last if last is None else last) + 1):
            state = self.snapshot(directory, step)
            where = (directory, step)
            n = self.events(step * self.dt)
            assert len(state) == n * len(self.cells), (where, "SD count", len(state), n * len(self.cells))
            total = sum(w for _, w in state.values())
            assert total == self.represented(n), (where, "cumulative represented dose", total, self.represented(n))
            doses = Counter()
            for cell, weight in state.values():
                doses[cell] += weight
            if doses:
                assert max(doses.values()) - min(doses.values()) <= 1, (where, "biased cumulative spatial dose")
            if step == first and previous is not None:
                # Checkpoint identities, positions and multiplicities must all be restored.
                assert state == previous, (where, "checkpoint particles changed on reload")
            else:
                old = previous if previous is not None else {}
                for ident, record in old.items():
                    assert state.get(ident) == record, (where, "particle lost/moved/reweighted", ident, record)
                born = {ident: state[ident] for ident in state.keys() - old.keys()}
                if n > self.events((step - 1) * self.dt):
                    by_cell = dict(born.values())
                    assert len(born) == len(by_cell) and by_cell.keys() == self.cells, (
                        where, "not exactly one birth per cell", len(born), len(by_cell))
                    assert {cpu for _, cpu in born} == set(range(ranks)), (where, "wrong birth ranks")
                    high = self.represented(n) - self.represented(n - 1) - len(self.cells) * self.low
                    assert 0 < high < len(self.cells), "fixture must exercise both integer weights"
                    assert Counter(by_cell.values()) == {self.low: len(self.cells) - high, self.low + 1: high}, (
                        where, "wrong floor/ceil mix", Counter(by_cell.values()), high)
                else:
                    assert not born, (where, "births between pulses or at the excluded upper endpoint", len(born))
            if reference is not None:
                assert Counter(state.values()) == Counter(reference[step].values()), (
                    where, "per-cell multiplicities differ from uninterrupted run")
            states[step] = state
            previous = state
        return states


def run(command, directory, *options, reject=None):
    directory.mkdir()
    log = directory / "simulation.log"
    with log.open("w") as stream:
        rc = subprocess.run([*command, *options], cwd=directory, stdout=stream,
                            stderr=subprocess.STDOUT, timeout=120, check=False).returncode
    if reject is None and rc == 0:
        return
    text = log.read_text(errors="replace")
    if reject is not None and rc != 0:
        # Do not accept echoed inputs, a launcher failure or an unrelated solver abort.
        diagnostics = [s for s in text.splitlines() if re.search(r"amrex::(abort|error)|assertion .*failed", s, re.I)]
        if any(re.search(reject, s, re.I) for s in diagnostics):
            return
    raise AssertionError(f"{log}: exit {rc}, expected {reject or 'success'}\n{text[-3000:]}")


def rejections(command, root):
    # SDInputs resolves injection.<key> and injection.0.<key> separately. Both
    # forbidden rate knobs must be rejected even if only a block supplies one.
    cases = [(tier + key, [INJ + tier + key + "=1"], key)
             for tier in ("", "0.") for key in ("sd_rate", "min_multiplicity")]
    cases += [
        ("period", [INJ + "period_s=0"], r"period"),
        ("rate", [INJ + "rate=0"], r"\brate\b|inj_rate"),
        ("zero_ppc", [INJ + "particles_per_cell=0"], r"particles_per_cell|\bppc\b"),
        ("empty_window", [INJ + "t_stop=5"], r"t_stop|time window"),
        ("too_little", [INJ + "rate=1e-13"], r"dose|multiplicity|physical particle|per cell"),
        ("sampled", [INJ + "multiplicity_type=sampled"], r"multiplicity|sampled"),
        ("moving", [INJ + "domain_velocity=1 0 0"], r"velocity|stationary|static"),
        ("partial", [INJ + "particle_box_hi=10000 10000 10000"], r"domain|box"),
        ("off_grid", [INJ + "period_s=11"], r"commensurate|align|multiple|time.?step|\bdt\b"),
        ("amr", ["amr.max_level=1", "amr.ref_ratio=2", "erf.coupling_type=TwoWay"],
         r"(?:period_s|pulse).*(?:level|AMR)|(?:level|AMR).*(?:period_s|pulse)"),
    ]
    for label, options, diagnostic in cases:
        run(command, root / label, "max_step=4", "erf.check_int=-1", "erf.plot_int_1=-1",
            *options, reject=diagnostic)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--ranks", type=int, choices=(1, 2), required=True)
    parser.add_argument("--reject", action="store_true")
    parser.add_argument("--reference-token", action="append", default=[],
                        help="one argv token of the one-rank launcher, for the MPI parity leg")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="-- [MPI launcher argv ...] ERF")
    args = parser.parse_args()
    if args.command[:1] != ["--"] or len(args.command) < 2:
        parser.error("supply -- followed by the solver command")
    solver = [args.command[-1], str(args.input.resolve()), "amrex.call_addr2line=0"]
    command = [*args.command[1:-1], *solver]
    root = Path(tempfile.mkdtemp(prefix="cell-pulse-", dir=Path.cwd()))
    print(f"Cell-pulse regression directory (retained on failure): {root}", flush=True)
    if args.reject:
        rejections(command, root)
    else:
        oracle = Oracle(args.input)
        straight = root / "straight"
        layout = "amr.max_grid_size_x=4" if args.ranks == 2 else "amr.max_grid_size_x=8"
        run(command, straight, layout)
        reference = oracle.history(straight, args.ranks)
        if args.ranks == 2:
            # Two ranks/two boxes vs one rank/one box, not just equal global totals.
            # IDs can differ between launches. Compare the per-cell weight of every pulse.
            serial = root / "one-box"
            last_birth = int((oracle.stop - oracle.period) / oracle.dt) + 1
            one_rank = [*args.reference_token, *solver]
            run(one_rank, serial, "amr.max_grid_size_x=8", "erf.check_int=-1", f"max_step={last_birth}")
            oracle.history(serial, 1, reference=reference, last=last_birth)
        # State at step 6 (15 s) is just BEFORE the second birth; step 9
        # (22.5 s) is BETWEEN births. Both checkpoints contain real particles.
        for start in (6, 9):
            directory = root / f"restart-{start}"
            checkpoint = straight / f"chk{start:05d}"
            run(command, directory, layout, f"erf.restart={checkpoint}", "erf.check_int=-1")
            oracle.history(directory, args.ranks, start, reference[start], reference)
    shutil.rmtree(root)
    print("PASS: cell-pulse " + ("input rejections" if args.reject else f"births and restarts ({args.ranks} rank(s))"))


if __name__ == "__main__":
    main()
