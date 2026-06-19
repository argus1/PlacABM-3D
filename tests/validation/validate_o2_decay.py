#!/usr/bin/env python3
"""
Phase 1 oxygen validation for PlacABM-3D.

Analytical check used:
  For spatially uniform oxygen with no flux/Dirichlet forcing and no cellular uptake,
  the diffusion-reaction PDE reduces to:
      dC/dt = -lambda * C
  so that:
      C(t1) = C(t0) * exp(-lambda * (t1 - t0))

This script compares simulated oxygen at two snapshots against that analytical solution
and reports relative L2 error.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import scipy.io as sio


EPS = 1e-16


def parse_multicellds_snapshot(xml_path: Path, substrate_name: str):
    tree = ET.parse(xml_path)
    root = tree.getroot()

    metadata = root.find("metadata")
    if metadata is None:
        raise RuntimeError(f"No <metadata> in {xml_path}")

    current_time_node = metadata.find("current_time")
    if current_time_node is None or current_time_node.text is None:
        raise RuntimeError(f"No <current_time> in {xml_path}")

    current_time = float(current_time_node.text)

    domain = root.find("./microenvironment/domain")
    if domain is None:
        raise RuntimeError(f"No <microenvironment>/<domain> in {xml_path}")

    variables = domain.find("variables")
    if variables is None:
        raise RuntimeError(f"No <variables> section in {xml_path}")

    substrate_index = None
    decay_rate = None
    substrate_names = []

    for i, var_node in enumerate(variables.findall("variable")):
        name = var_node.get("name", "")
        substrate_names.append(name)
        if name == substrate_name:
            substrate_index = i
            decay_node = var_node.find("./physical_parameter_set/decay_rate")
            if decay_node is None or decay_node.text is None:
                raise RuntimeError(f"Missing decay_rate for substrate '{substrate_name}' in {xml_path}")
            decay_rate = float(decay_node.text)

    if substrate_index is None:
        raise RuntimeError(
            f"Substrate '{substrate_name}' not found in {xml_path}. Found: {substrate_names}"
        )

    data_node = domain.find("./data/filename")
    if data_node is None or data_node.text is None:
        raise RuntimeError(f"No microenvironment MATLAB filename found in {xml_path}")

    mat_path = xml_path.parent / data_node.text.strip()
    if not mat_path.exists():
        raise FileNotFoundError(f"Referenced MATLAB file not found: {mat_path}")

    mat = sio.loadmat(mat_path)
    if "multiscale_microenvironment" not in mat:
        raise RuntimeError(f"Expected 'multiscale_microenvironment' in {mat_path}")

    m = mat["multiscale_microenvironment"]
    # Rows 0..2: xyz, row 3: volume, rows 4+: substrate values
    row = 4 + substrate_index
    if row >= m.shape[0]:
        raise RuntimeError(
            f"Substrate row out of bounds in {mat_path}. row={row}, shape={m.shape}"
        )

    concentration = np.asarray(m[row, :], dtype=float)

    return {
        "time": current_time,
        "decay_rate": float(decay_rate),
        "concentration": concentration,
        "xml_path": str(xml_path),
        "mat_path": str(mat_path),
    }


def compute_decay_validation(c0: np.ndarray, c1: np.ndarray, decay_rate: float, dt: float):
    expected = c0 * math.exp(-decay_rate * dt)

    diff = c1 - expected
    rel_l2 = float(np.linalg.norm(diff) / max(np.linalg.norm(expected), EPS))
    max_abs = float(np.max(np.abs(diff)))
    mean_abs = float(np.mean(np.abs(diff)))

    return {
        "rel_l2_error": rel_l2,
        "max_abs_error": max_abs,
        "mean_abs_error": mean_abs,
        "expected_mean": float(np.mean(expected)),
        "simulated_mean": float(np.mean(c1)),
    }


def main():
    parser = argparse.ArgumentParser(description="Validate oxygen decay against analytical solution.")
    parser.add_argument("--output-dir", default="output", help="Directory with MultiCellDS XML/MAT files")
    parser.add_argument("--start-xml", default="initial.xml", help="Start snapshot XML filename")
    parser.add_argument("--end-xml", default="final.xml", help="End snapshot XML filename")
    parser.add_argument("--substrate", default="oxygen", help="Substrate name to validate")
    parser.add_argument(
        "--max-rel-l2",
        type=float,
        default=1e-2,
        help="Pass threshold for relative L2 error (default: 1e-2)",
    )
    parser.add_argument("--json", action="store_true", help="Print JSON report")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    start_xml = output_dir / args.start_xml
    end_xml = output_dir / args.end_xml

    if not start_xml.exists():
        raise FileNotFoundError(f"Start snapshot not found: {start_xml}")
    if not end_xml.exists():
        raise FileNotFoundError(f"End snapshot not found: {end_xml}")

    s0 = parse_multicellds_snapshot(start_xml, args.substrate)
    s1 = parse_multicellds_snapshot(end_xml, args.substrate)

    dt = s1["time"] - s0["time"]
    if dt < 0:
        raise RuntimeError(f"Negative dt detected: t_end={s1['time']} < t_start={s0['time']}")

    metrics = compute_decay_validation(
        s0["concentration"],
        s1["concentration"],
        s1["decay_rate"],
        dt,
    )

    report = {
        "substrate": args.substrate,
        "start_xml": str(start_xml),
        "end_xml": str(end_xml),
        "t_start": s0["time"],
        "t_end": s1["time"],
        "dt": dt,
        "decay_rate": s1["decay_rate"],
        **metrics,
        "threshold_rel_l2": args.max_rel_l2,
        "pass": metrics["rel_l2_error"] <= args.max_rel_l2,
    }

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("[Phase1 O2 Validation]")
        print(f"  substrate         : {report['substrate']}")
        print(f"  snapshots         : {report['start_xml']} -> {report['end_xml']}")
        print(f"  dt                : {report['dt']:.6g} min")
        print(f"  decay_rate        : {report['decay_rate']:.6g} 1/min")
        print(f"  expected_mean     : {report['expected_mean']:.6g}")
        print(f"  simulated_mean    : {report['simulated_mean']:.6g}")
        print(f"  rel_l2_error      : {report['rel_l2_error']:.6g}")
        print(f"  max_abs_error     : {report['max_abs_error']:.6g}")
        print(f"  mean_abs_error    : {report['mean_abs_error']:.6g}")
        print(f"  threshold_rel_l2  : {report['threshold_rel_l2']:.6g}")
        print(f"  PASS              : {report['pass']}")

    if report["pass"]:
        return 0

    return 1


if __name__ == "__main__":
    sys.exit(main())
