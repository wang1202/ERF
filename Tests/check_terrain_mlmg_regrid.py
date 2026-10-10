#!/usr/bin/env python3
"""Verify a restart actually changes the level-1 BoxArray and recreates its grids."""
import argparse
import pathlib
import re
import sys

BOX_ARRAY = re.compile(r"\(\d+ 0\n(?:\(\([^\n]+\)\)\n)+\)")


def parse_box_array(path: pathlib.Path) -> str:
    match = BOX_ARRAY.search(path.read_text())
    if match is None:
        raise ValueError(f"could not find BoxArray in {path}")
    return match.group(0)


def compare(checkpoint: pathlib.Path, final: pathlib.Path, log: pathlib.Path) -> None:
    old = parse_box_array(checkpoint)
    new = parse_box_array(final)
    if old == new:
        raise ValueError("level-1 BoxArray is unchanged across restart/regrid")
    if "REMAKING WITH NEW BA AT LEVEL 1" not in log.read_text():
        raise ValueError("restart log does not prove level-1 grids were remade")
    print("checkpoint level-1 BoxArray:", " ".join(old.splitlines()))
    print("post-regrid level-1 BoxArray:", " ".join(new.splitlines()))
    print("Level-1 BoxArray changed across the restart regrid.")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint-cell-h", type=pathlib.Path)
    parser.add_argument("--final-cell-h", type=pathlib.Path)
    parser.add_argument("--restart-log", type=pathlib.Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        try:
            if parse_box_array_text("(1 0\n((0,0,0) (7,7,7) (0,0,0))\n)\n") == parse_box_array_text(
                "(1 0\n((8,0,0) (15,7,7) (0,0,0))\n)\n"
            ):
                raise AssertionError("distinct layouts compared equal")
            try:
                parse_box_array_text("malformed")
            except ValueError:
                print("Terrain MLMG regrid checker self-test passed.")
                return 0
            raise AssertionError("malformed header was accepted")
        except AssertionError as exc:
            print(f"self-test failed: {exc}", file=sys.stderr)
            return 1
    if args.checkpoint_cell_h is None or args.final_cell_h is None or args.restart_log is None:
        parser.error("provide all three file arguments or use --self-test")
    try:
        compare(args.checkpoint_cell_h, args.final_cell_h, args.restart_log)
    except (OSError, ValueError) as exc:
        print(f"Terrain MLMG regrid check failed: {exc}", file=sys.stderr)
        return 1
    return 0


def parse_box_array_text(text: str) -> str:
    match = BOX_ARRAY.search(text)
    if match is None:
        raise ValueError("could not find BoxArray")
    return match.group(0)


if __name__ == "__main__":
    raise SystemExit(main())
