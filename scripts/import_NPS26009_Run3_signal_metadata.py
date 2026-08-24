#!/usr/bin/env python3
"""Import the private NPS-26-009 Run-3 Z' signal metadata."""

from __future__ import annotations

import copy
import json
import math
import os
from pathlib import Path
import re
import tempfile


SOURCE_BASE = Path(
    "/data9/Users/taehee_public/HerwigSample/RKZp_13p6TeV/RS/"
    "samples_nEvt-100000/Info"
)
ERAS = ("2022", "2022EE", "2023", "2023BPix")
ALLOWED_MASSES = (12, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)
ALIAS_PATTERN = re.compile(r"^Zp_M-(\d+)_Pt-(\d+)to(\d+)_hw7$")
XSEC_SCALE = 4.0 * math.pi


def load_json(path: Path) -> dict:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Cannot read valid JSON from {path}: {error}") from error
    if not isinstance(value, dict):
        raise RuntimeError(f"Expected a JSON object in {path}")
    return value


def sort_key(alias: str) -> tuple[int, int, int]:
    match = ALIAS_PATTERN.fullmatch(alias)
    if match is None:
        raise RuntimeError(f"Invalid signal alias: {alias}")
    return tuple(int(value) for value in match.groups())


def write_json_atomic(path: Path, value: dict) -> None:
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        json.dump(value, handle, indent=4)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.chmod(path.stat().st_mode if path.exists() else 0o644)
    temporary.replace(path)


def write_lines_atomic(path: Path, lines: list[str]) -> None:
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write("".join(f"{line}\n" for line in lines))
        temporary = Path(handle.name)
    if path.exists():
        temporary.chmod(path.stat().st_mode)
    temporary.replace(path)


def main() -> None:
    destination_env = os.environ.get("SKNANO_DATA")
    if not destination_env:
        raise RuntimeError("SKNANO_DATA is not set; source setup.sh before running")

    destination_base = Path(destination_env).expanduser().resolve()
    repository = Path(__file__).resolve().parents[1]
    signal_list = repository / "SampleLists" / "Run3signal.txt"
    allowed_masses = set(ALLOWED_MASSES)

    aliases_by_era: dict[str, list[str]] = {}
    source_common_by_era: dict[str, dict] = {}
    destination_common_by_era: dict[str, dict] = {}
    source_files_by_era: dict[str, dict[str, Path]] = {}

    # Complete all consistency checks before writing any destination file.
    for era in ERAS:
        source_sample = SOURCE_BASE / era / "Sample"
        source_for_snu = source_sample / "ForSNU"
        source_common_path = source_sample / "CommonSampleInfo.json"
        destination_common_path = (
            destination_base / era / "Sample" / "CommonSampleInfo.json"
        )
        destination_for_snu = destination_base / era / "Sample" / "ForSNU"

        if not source_for_snu.is_dir():
            raise RuntimeError(f"Missing source directory: {source_for_snu}")
        if not destination_for_snu.is_dir():
            raise RuntimeError(f"Missing destination directory: {destination_for_snu}")

        source_common = load_json(source_common_path)
        destination_common = load_json(destination_common_path)
        source_files: dict[str, Path] = {}

        for path in source_for_snu.iterdir():
            if not path.is_file() or path.suffix != ".json":
                continue
            alias = path.stem
            match = ALIAS_PATTERN.fullmatch(alias)
            if match is None or int(match.group(1)) not in allowed_masses:
                continue
            if alias in source_files:
                raise RuntimeError(f"Duplicate alias in {era}: {alias}")
            source_files[alias] = path

        aliases = sorted(source_files, key=sort_key)
        masses = {sort_key(alias)[0] for alias in aliases}
        if masses != allowed_masses:
            missing = sorted(allowed_masses - masses)
            unexpected = sorted(masses - allowed_masses)
            raise RuntimeError(
                f"Mass mismatch in {era}: missing={missing}, unexpected={unexpected}"
            )
        if len(aliases) != len(set(aliases)):
            raise RuntimeError(f"Duplicate aliases discovered in {era}")

        for alias in aliases:
            if alias not in source_common:
                raise RuntimeError(
                    f"{alias} has no entry in source {source_common_path}"
                )
            metadata = source_common[alias]
            if not isinstance(metadata, dict):
                raise RuntimeError(f"Malformed metadata for {alias} in {era}")
            if metadata.get("isMC") != 1:
                raise RuntimeError(f"{alias} in {era} is not marked as MC")
            source_xsec = metadata.get("xsec")
            if (
                isinstance(source_xsec, bool)
                or not isinstance(source_xsec, (int, float))
                or not math.isfinite(source_xsec)
                or source_xsec == 0
            ):
                raise RuntimeError(f"Invalid xsec for {alias} in {era}: {source_xsec}")

            for_snu = load_json(source_files[alias])
            if for_snu.get("name") != alias:
                raise RuntimeError(f"ForSNU name mismatch for {alias} in {era}")
            if for_snu.get("isMC") != 1:
                raise RuntimeError(f"ForSNU file is not MC for {alias} in {era}")
            if not isinstance(for_snu.get("path"), list):
                raise RuntimeError(f"ForSNU path is not a list for {alias} in {era}")
            for_snu_xsec = for_snu.get("xsec")
            if (
                isinstance(for_snu_xsec, bool)
                or not isinstance(for_snu_xsec, (int, float))
                or not math.isfinite(for_snu_xsec)
                or for_snu_xsec == 0
            ):
                raise RuntimeError(
                    f"Invalid ForSNU xsec for {alias} in {era}: {for_snu_xsec}"
                )
            if for_snu_xsec != source_xsec:
                raise RuntimeError(
                    f"Source Common/ForSNU xsec mismatch for {alias} in {era}: "
                    f"{source_xsec} != {for_snu_xsec}"
                )

        aliases_by_era[era] = aliases
        source_common_by_era[era] = source_common
        destination_common_by_era[era] = destination_common
        source_files_by_era[era] = source_files

    reference_aliases = aliases_by_era[ERAS[0]]
    reference_set = set(reference_aliases)
    for era in ERAS[1:]:
        era_set = set(aliases_by_era[era])
        if era_set != reference_set:
            raise RuntimeError(
                f"Signal alias mismatch for {era}: "
                f"only_{ERAS[0]}={sorted(reference_set - era_set)}, "
                f"only_{era}={sorted(era_set - reference_set)}"
            )

    examples: dict[str, list[tuple[str, float, float, float]]] = {}
    for era in ERAS:
        aliases = aliases_by_era[era]
        source_common = source_common_by_era[era]
        destination_common = destination_common_by_era[era]
        original_non_signal = {
            alias: copy.deepcopy(metadata)
            for alias, metadata in destination_common.items()
            if alias not in reference_set
        }

        for alias in aliases:
            imported = copy.deepcopy(source_common[alias])
            source_xsec = imported["xsec"]
            imported["xsec"] = source_xsec * XSEC_SCALE
            destination_common[alias] = imported

            other_source = {k: v for k, v in source_common[alias].items() if k != "xsec"}
            other_destination = {k: v for k, v in imported.items() if k != "xsec"}
            if other_destination != other_source:
                raise RuntimeError(f"Non-xsec metadata changed for {alias} in {era}")
            ratio = imported["xsec"] / source_xsec
            if not math.isclose(ratio, XSEC_SCALE, rel_tol=1e-15, abs_tol=0.0):
                raise RuntimeError(f"xsec validation failed for {alias} in {era}")

        after_non_signal = {
            alias: metadata
            for alias, metadata in destination_common.items()
            if alias not in reference_set
        }
        if after_non_signal != original_non_signal:
            raise RuntimeError(f"Non-signal destination metadata changed in {era}")

        common_path = destination_base / era / "Sample" / "CommonSampleInfo.json"
        write_json_atomic(common_path, destination_common)

        destination_for_snu = destination_base / era / "Sample" / "ForSNU"
        for alias in aliases:
            source_path = source_files_by_era[era][alias]
            destination_path = destination_for_snu / source_path.name
            source_for_snu = load_json(source_path)
            imported_for_snu = copy.deepcopy(source_for_snu)
            source_for_snu_xsec = imported_for_snu["xsec"]
            imported_for_snu["xsec"] = source_for_snu_xsec * XSEC_SCALE

            other_source_for_snu = {
                key: value for key, value in source_for_snu.items() if key != "xsec"
            }
            other_destination_for_snu = {
                key: value for key, value in imported_for_snu.items() if key != "xsec"
            }
            if other_destination_for_snu != other_source_for_snu:
                raise RuntimeError(
                    f"Non-xsec ForSNU metadata changed for {alias} in {era}"
                )
            ratio = imported_for_snu["xsec"] / source_for_snu_xsec
            if not math.isclose(ratio, XSEC_SCALE, rel_tol=1e-15, abs_tol=0.0):
                raise RuntimeError(f"ForSNU xsec validation failed for {alias} in {era}")

            write_json_atomic(destination_path, imported_for_snu)
            if load_json(destination_path) != imported_for_snu:
                raise RuntimeError(f"ForSNU write validation failed for {alias} in {era}")

        example_indexes = (0, len(aliases) // 2, len(aliases) - 1)
        examples[era] = []
        for index in example_indexes:
            alias = aliases[index]
            source_xsec = source_common[alias]["xsec"]
            destination_xsec = destination_common[alias]["xsec"]
            examples[era].append(
                (alias, source_xsec, destination_xsec, destination_xsec / source_xsec)
            )

    write_lines_atomic(signal_list, reference_aliases)

    bins_by_mass = {
        mass: sum(sort_key(alias)[0] == mass for alias in reference_aliases)
        for mass in ALLOWED_MASSES
    }
    print(f"Source: {SOURCE_BASE}")
    print(f"Destination SKNANO_DATA: {destination_base}")
    print(f"Imported {len(reference_aliases)} samples in each era: {', '.join(ERAS)}")
    print("pT bins by mass:")
    for mass in ALLOWED_MASSES:
        print(f"  M={mass}: {bins_by_mass[mass]}")
    print("First 10 entries:")
    for alias in reference_aliases[:10]:
        print(f"  {alias}")
    print("Last 10 entries:")
    for alias in reference_aliases[-10:]:
        print(f"  {alias}")
    print(f"xsec scale (4*pi): {XSEC_SCALE:.15g}")
    for era in ERAS:
        print(f"Representative xsec checks for {era}:")
        for alias, source_xsec, destination_xsec, ratio in examples[era]:
            print(
                f"  {alias}: source={source_xsec:.17g}, "
                f"destination={destination_xsec:.17g}, ratio={ratio:.17g}"
            )
    print("Validated Common and ForSNU xsec ratio 4*pi for every imported entry.")
    print("Validated that all pre-existing non-signal metadata is unchanged.")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as error:
        raise SystemExit(f"ERROR: {error}") from error
