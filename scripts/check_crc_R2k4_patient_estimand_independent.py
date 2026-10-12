#!/usr/bin/env python
"""Independent arithmetic audit of the frozen broad-pair patient estimand.

Date: 2026-10-12
Random seed: 42 (PCG64 patient bootstrap, for archived-CI comparison only)
Input: pair-by-IS-group edge counts and E4 tables newly reconstructed in this run.
This independent aggregation does not itself rebuild the graph, labels or scores;
those steps are performed upstream by the public-input runner.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


PAIRS = ("tumor->tumor", "tumor->stroma", "stroma->tumor", "stroma->stroma")
TOL = 1e-12
# Internal `archived_*` names come from the original audit. In this release the
# CLI passes E4 tables reconstructed earlier in the same public-input run.


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_bool(value: str) -> bool:
    if value not in ("True", "False"):
        raise ValueError(f"Unexpected boolean value: {value!r}")
    return value == "True"


def close(actual: float, archived: str, label: str, errors: list[str]) -> None:
    if archived == "":
        errors.append(f"{label}: archived value missing")
    elif not np.isclose(actual, float(archived), rtol=0, atol=TOL):
        errors.append(f"{label}: independent {actual:.17g}, archived {archived}")


def pair_effect(rows: dict[bool, dict[str, str]], label: str, errors: list[str]) -> tuple[float, int] | None:
    if set(rows) != {False, True}:
        errors.append(f"{label}: IS-group rows are not exactly low/high")
        return None
    low, high = rows[False], rows[True]
    low_ok, high_ok = as_bool(low["pair_passed_minimum"]), as_bool(high["pair_passed_minimum"])
    if low_ok != high_ok:
        errors.append(f"{label}: pair gate disagrees between IS groups")
    if not (low_ok and high_ok):
        if int(low["n_pair_retained"]) or int(high["n_pair_retained"]):
            errors.append(f"{label}: failed pair has retained edges")
        return None
    counts = [int(low["n_pair_retained"]), int(high["n_pair_retained"])]
    if min(counts) < 30:
        errors.append(f"{label}: passed pair has fewer than 30 edges in one IS group")
        return None
    contrasts = []
    for outcome in ("observed", "shifted"):
        field = f"n_{outcome}_SS_high_on_retained"
        numerators = [int(low[field]), int(high[field])]
        if any(n < 0 or n > d for n, d in zip(numerators, counts)):
            errors.append(f"{label}: {outcome} successes outside edge denominator")
        contrasts.append(numerators[1] / counts[1] - numerators[0] / counts[0])
    return contrasts[0] - contrasts[1], min(counts)


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> int:
    errors: list[str] = []
    input_paths = {
        "pair_ledger": args.pair_ledger,
        "e4_shifts": args.archived_shifts,
        "e4_sections": args.archived_sections,
        "e4_patients": args.archived_patients,
        "e4_summary": args.archived_summary,
    }
    ledger = read_rows(args.pair_ledger)
    shifts = read_rows(args.archived_shifts)
    sections = read_rows(args.archived_sections)
    patients = read_rows(args.archived_patients)
    summaries = read_rows(args.archived_summary)

    # Keys are dataset, section, row offset, column offset. Keep technical and
    # biological denominators distinct; a failed pair never contributes a zero.
    cells: dict[tuple[str, str, str, str], dict[str, dict[bool, dict[str, str]]]] = defaultdict(lambda: defaultdict(dict))
    for row in ledger:
        key = tuple(row[k] for k in ("dataset", "sample", "shift_array_row", "shift_array_col"))
        pair = row["pathology_pair"]
        group = as_bool(row["IS_high_source"])
        if pair not in PAIRS or group in cells[key][pair]:
            errors.append(f"Invalid or repeated pair/group: {key}/{pair}/{group}")
        cells[key][pair][group] = row

    shift_lookup = {tuple(r[k] for k in ("dataset", "sample", "shift_array_row", "shift_array_col")): r for r in shifts}
    if len(shift_lookup) != len(shifts):
        errors.append("Archived shifts contain duplicate keys")
    offset_rows: list[dict] = []
    passing_by_section: dict[tuple[str, str], list[float]] = defaultdict(list)
    for key in sorted(cells):
        pair_groups = cells[key]
        if set(pair_groups) != set(PAIRS):
            errors.append(f"{key}: expected all four pathology pairs")
        pair_values: list[tuple[float, int]] = []
        pair_names: list[str] = []
        for pair in PAIRS:
            if pair not in pair_groups:
                continue
            finding = pair_effect(pair_groups[pair], f"{key}/{pair}", errors)
            if finding is not None:
                pair_values.append(finding)
                pair_names.append(pair)
        archived = shift_lookup.get(key)
        if archived is None:
            errors.append(f"{key}: missing archived shift")
            continue
        if len(pair_values) != int(archived["n_matched_pathology_pairs"]):
            errors.append(f"{key}: qualifying pair count differs")
        ledger_statuses = {row["status"] for groups in pair_groups.values() for row in groups.values()}
        if len(ledger_statuses) != 1:
            errors.append(f"{key}: pair/group rows disagree on offset status")
        ledger_status = next(iter(ledger_statuses))
        if ledger_status != archived["status"]:
            errors.append(f"{key}: ledger status {ledger_status} != archived {archived['status']}")
        if ledger_status == "PASS":
            total_weight = sum(weight for _, weight in pair_values)
            if total_weight <= 0:
                errors.append(f"{key}: PASS with no qualifying pair")
                continue
            effect = sum(delta * weight for delta, weight in pair_values) / total_weight
            close(effect, archived["observed_minus_shifted_IS_contrast"], f"{key} offset effect", errors)
            passing_by_section[key[:2]].append(effect)
        else:
            effect = None
            total_weight = 0
            if pair_values or archived["observed_minus_shifted_IS_contrast"]:
                errors.append(f"{key}: failed shift has pair/effect")
        offset_rows.append({"dataset": key[0], "sample": key[1], "shift_array_row": key[2],
                            "shift_array_col": key[3], "status": ledger_status,
                            "qualifying_pairs": "|".join(pair_names), "n_qualifying_pairs": len(pair_values),
                            "total_pair_weight": total_weight, "effect_fraction": effect if effect is not None else ""})
    if len(cells) != 1920 or len(ledger) != 15360:
        errors.append(f"Ledger size {len(cells)} offsets/{len(ledger)} pair-group rows, expected 1920/15360")

    section_rows: list[dict] = []
    eligible_by_patient: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in sections:
        dataset, sample, patient = row["dataset"], row["sample"], row["patient_id"]
        if dataset == "GSE294385" and row["tissue_context"] != "Primary_CRC":
            continue
        if dataset not in ("GSE294385", "Valdeolivas_2024"):
            continue
        values = passing_by_section.get((dataset, sample), [])
        if row["n_accepted_shifts"] and int(float(row["n_accepted_shifts"])) != len(values):
            errors.append(f"{dataset}/{sample}: accepted shift count differs")
        passed = row["status"] == "PASS"
        if passed != (len(values) >= 24):
            errors.append(f"{dataset}/{sample}: section gate disagrees with 24-offset rule")
        if passed:
            section_effect = sum(values) / len(values)
            close(section_effect, row["pathology_matched_shift_effect_fraction"],
                  f"{dataset}/{sample} section effect", errors)
            context = "all" if dataset == "Valdeolivas_2024" else "Primary_CRC"
            eligible_by_patient[(dataset, context, patient)].append(section_effect)
        else:
            section_effect = None
            if row["pathology_matched_shift_effect_fraction"]:
                errors.append(f"{dataset}/{sample}: failed section has effect")
        section_rows.append({"dataset": dataset, "sample": sample, "patient_id": patient,
                             "context": "all" if dataset == "Valdeolivas_2024" else "Primary_CRC",
                             "status": row["status"], "n_accepted_offsets": len(values),
                             "effect_fraction": section_effect if section_effect is not None else ""})

    patient_rows: list[dict] = []
    patient_lookup = {(r["dataset"], r["context"], r["patient_id"]): r for r in patients}
    for key in sorted(eligible_by_patient):
        values = eligible_by_patient[key]
        effect = sum(values) / len(values)
        archived = patient_lookup.get(key)
        if archived is None:
            errors.append(f"{key}: archived patient absent")
        else:
            close(effect, archived["pathology_matched_shift_effect_fraction"], f"{key} patient effect", errors)
        patient_rows.append({"dataset": key[0], "context": key[1], "patient_id": key[2],
                             "n_passing_sections": len(values), "effect_fraction": effect})

    cohort_rows: list[dict] = []
    summary_lookup = {(r["dataset"], r["context"]): r for r in summaries}
    for dataset, context in (("Valdeolivas_2024", "all"), ("GSE294385", "Primary_CRC")):
        cohort_patients = [r for r in patient_rows if r["dataset"] == dataset and r["context"] == context]
        values = np.array([r["effect_fraction"] for r in cohort_patients], dtype=float)
        if not len(values):
            errors.append(f"{dataset}/{context}: no eligible patients")
            continue
        estimate = float(values.mean())
        archive = summary_lookup[(dataset, context)]
        close(estimate, archive["estimate_fraction"], f"{dataset}/{context} cohort effect", errors)
        if len(values) != int(archive["n_eligible_patients"]):
            errors.append(f"{dataset}/{context}: patient count differs")
        n_sections = sum(r["n_passing_sections"] for r in cohort_patients)
        if n_sections != int(archive["n_eligible_sections"]):
            errors.append(f"{dataset}/{context}: passing section count differs")
        loo = [(sum(values) - value) / (len(values) - 1) for value in values]
        close(min(loo), archive["leave_one_patient_out_min_fraction"], f"{dataset}/{context} LOO min", errors)
        close(max(loo), archive["leave_one_patient_out_max_fraction"], f"{dataset}/{context} LOO max", errors)
        rng = np.random.default_rng(42)
        draws = rng.choice(values, size=(10000, len(values)), replace=True).mean(axis=1)
        ci_low, ci_high = (float(x) for x in np.quantile(draws, [0.025, 0.975]))
        close(ci_low, archive["ci_lower_fraction"], f"{dataset}/{context} bootstrap low", errors)
        close(ci_high, archive["ci_upper_fraction"], f"{dataset}/{context} bootstrap high", errors)
        cohort_rows.append({"dataset": dataset, "context": context, "n_patients": len(values),
                            "n_passing_sections": n_sections, "n_accepted_section_offsets": sum(
                                r["n_accepted_offsets"] for r in section_rows
                                if r["dataset"] == dataset and r["status"] == "PASS"),
                            "estimate_fraction": estimate, "ci_lower_fraction": ci_low,
                            "ci_upper_fraction": ci_high, "loo_min_fraction": min(loo),
                            "loo_max_fraction": max(loo)})

    args.out.mkdir(parents=True, exist_ok=True)
    write_csv(args.out / "independent_offset_effects.csv", offset_rows, list(offset_rows[0]))
    write_csv(args.out / "independent_section_effects.csv", section_rows, list(section_rows[0]))
    write_csv(args.out / "independent_patient_effects.csv", patient_rows, list(patient_rows[0]))
    write_csv(args.out / "independent_cohort_effects.csv", cohort_rows, list(cohort_rows[0]))
    report = {"status": "PASS" if not errors else "FAIL", "errors": errors,
              "scope": "independent arithmetic from newly rebuilt pair/group counts and E4 tables; graph and scores are rebuilt upstream",
              "input_sha256": {name: file_sha256(path) for name, path in input_paths.items()},
              "n_pair_group_rows": len(ledger), "n_offsets": len(cells),
              "cohorts": cohort_rows}
    (args.out / "audit_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if not errors else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-ledger", type=Path, required=True)
    parser.add_argument("--e4-shifts", dest="archived_shifts", type=Path, required=True)
    parser.add_argument("--e4-sections", dest="archived_sections", type=Path, required=True)
    parser.add_argument("--e4-patients", dest="archived_patients", type=Path, required=True)
    parser.add_argument("--e4-summary", dest="archived_summary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
