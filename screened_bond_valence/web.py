"""Gradio web interface for CIF-based bond valence fitting."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

import gradio as gr

from .api import ScreenedBondValenceService, build_material_from_cif


def calculate_bv_params(cation: str, anion: str, cif_file: Any) -> tuple[str | None, str]:
    try:
        cif_path = Path(getattr(cif_file, "name", cif_file))
        material = build_material_from_cif(cif_path, cation=cation, anion=anion)
        service = ScreenedBondValenceService()
        result = service.fit_material(material)

        if result.theoretical is None or not result.theoretical.has_solution:
            return None, "No network solution found"

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as handle:
            json.dump(result.theoretical.bond_valences, handle, indent=2)
            sij_path = handle.name

        summary = result.aggregate()
        if summary is None:
            reason = "; ".join(result.failure_reasons) if result.failure_reasons else "No consistent R0/B values found"
            return sij_path, reason

        return sij_path, (
            f"R0: {summary.r0:.4f} A\n"
            f"B: {summary.b:.4f}\n"
            f"Algorithms used: {summary.n_algos}"
        )
    except Exception as exc:  # pragma: no cover - interactive surface
        return None, f"Error: {exc}"


def create_interface() -> gr.Interface:
    return gr.Interface(
        fn=calculate_bv_params,
        inputs=[
            gr.Textbox(label="Cation", placeholder="Enter cation (for example Li)"),
            gr.Textbox(label="Anion", placeholder="Enter anion (for example O)"),
            gr.File(label="CIF File", file_types=[".cif"]),
        ],
        outputs=[
            gr.File(label="Download Sij Results"),
            gr.Text(label="Bond Valence Parameters"),
        ],
        title="Bond Valence Parameter Search",
        description="Upload a CIF file to calculate bond valence parameters (R0 and B).",
    )


def launch_app(**launch_kwargs: Any) -> None:
    create_interface().launch(**launch_kwargs)
