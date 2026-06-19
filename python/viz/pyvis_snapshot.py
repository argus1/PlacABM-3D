#!/usr/bin/env python3
"""Render PlacABM tissue masks as a quick 3D snapshot with PyVista.

This script expects mask CSVs in ``config/tissue_masks`` where each line begins
with a voxel index (additional columns are ignored).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List
import xml.etree.ElementTree as ET


MASK_FILES = {
    "decidua": "decidua.csv",
    "syncytiotrophoblast": "syncytiotrophoblast.csv",
    "villous_stroma": "villous_stroma.csv",
    "fetal_capillary": "fetal_capillary.csv",
}

MASK_COLORS = {
    "decidua": "gold",
    "syncytiotrophoblast": "tomato",
    "villous_stroma": "cornflowerblue",
    "fetal_capillary": "limegreen",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render PlacABM tissue masks with PyVista.")
    parser.add_argument(
        "--xml",
        type=Path,
        default=Path("config/placabm_default.xml"),
        help="Path to PhysiCell XML config.",
    )
    parser.add_argument(
        "--mask-dir",
        type=Path,
        default=Path("config/tissue_masks"),
        help="Directory containing mask CSV files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/mask_snapshot.png"),
        help="Output PNG path.",
    )
    parser.add_argument("--point-size", type=float, default=7.0, help="Mask point size.")
    return parser.parse_args()


def _get_float(root: ET.Element, path: str) -> float:
    node = root.find(path)
    if node is None or node.text is None:
        raise ValueError(f"Missing XML value at path: {path}")
    return float(node.text)


def parse_domain(xml_path: Path):
    root = ET.parse(xml_path).getroot()

    x_min = _get_float(root, "domain/x_min")
    x_max = _get_float(root, "domain/x_max")
    y_min = _get_float(root, "domain/y_min")
    y_max = _get_float(root, "domain/y_max")
    z_min = _get_float(root, "domain/z_min")
    z_max = _get_float(root, "domain/z_max")
    dx = _get_float(root, "domain/dx")
    dy = _get_float(root, "domain/dy")
    dz = _get_float(root, "domain/dz")

    nx = int(round((x_max - x_min) / dx))
    ny = int(round((y_max - y_min) / dy))
    nz = int(round((z_max - z_min) / dz))

    if nx <= 0 or ny <= 0 or nz <= 0:
        raise ValueError("Computed invalid mesh dimensions from XML domain settings.")

    return (x_min, y_min, z_min), (dx, dy, dz), (nx, ny, nz)


def read_mask_indices(mask_csv: Path) -> List[int]:
    indices: List[int] = []
    if not mask_csv.exists():
        return indices

    with mask_csv.open("r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            first = row[0].strip()
            if not first or first.startswith("#"):
                continue
            try:
                indices.append(int(first))
            except ValueError:
                continue
    return indices


def voxel_index_to_xyz(index: int, origin, spacing, shape):
    x0, y0, z0 = origin
    dx, dy, dz = spacing
    nx, ny, nz = shape

    if index < 0 or index >= nx * ny * nz:
        raise ValueError("Voxel index out of range")

    i = index % nx
    j = (index // nx) % ny
    k = index // (nx * ny)

    x = x0 + (i + 0.5) * dx
    y = y0 + (j + 0.5) * dy
    z = z0 + (k + 0.5) * dz
    return x, y, z


def build_points(mask_indices: List[int], origin, spacing, shape):
    points = []
    for idx in mask_indices:
        try:
            points.append(voxel_index_to_xyz(idx, origin, spacing, shape))
        except ValueError:
            continue
    return points


def render_masks(points_by_layer: Dict[str, List[tuple]], output_path: Path, point_size: float):
    try:
        import pyvista as pv
    except ImportError as exc:
        raise RuntimeError(
            "PyVista is required. Install with `pip install pyvista` in your Python environment."
        ) from exc

    output_path.parent.mkdir(parents=True, exist_ok=True)

    plotter = pv.Plotter(off_screen=True)
    plotter.set_background("white")

    for layer, points in points_by_layer.items():
        if not points:
            continue
        poly = pv.PolyData(points)
        plotter.add_mesh(
            poly,
            color=MASK_COLORS[layer],
            render_points_as_spheres=True,
            point_size=point_size,
            label=f"{layer} (n={len(points)})",
        )

    plotter.add_legend(size=(0.25, 0.25), bcolor="white")
    plotter.add_axes()
    plotter.camera_position = "iso"
    plotter.screenshot(str(output_path))


def main() -> None:
    args = parse_args()

    origin, spacing, shape = parse_domain(args.xml)

    points_by_layer: Dict[str, List[tuple]] = {}
    for layer, filename in MASK_FILES.items():
        mask_path = args.mask_dir / filename
        indices = read_mask_indices(mask_path)
        points_by_layer[layer] = build_points(indices, origin, spacing, shape)

    if all(len(points) == 0 for points in points_by_layer.values()):
        raise RuntimeError(
            f"No mask points found in {args.mask_dir}. "
            "Generate mask CSVs first or run the simulation once with exported masks."
        )

    render_masks(points_by_layer, args.output, args.point_size)
    print(f"Saved mask snapshot: {args.output}")


if __name__ == "__main__":
    main()
