"""Gradio web interface for bond valence parameter calculation from CIF files."""

import gradio as gr
import json
import tempfile
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN
from BVparams_search import TheoreticalBondValenceSolver, BVParamSolver


def calculate_bv_params(cation, anion, cif_file):
    try:
        structure = Structure.from_file(cif_file.name)

        # Get oxidation states from CIF
        element2charge = {}
        for site in structure:
            if not hasattr(site.specie, "oxi_state"):
                return None, "Oxidation states missing in CIF!"
            element2charge[site.specie.symbol] = site.specie.oxi_state

        if cation not in element2charge or anion not in element2charge:
            return None, "Cation/Anion not found in structure!"

        # Generate bonded structure
        cnn = CrystalNN(cation_anion=True, distance_cutoffs=[0, 1])
        bonded_structure = cnn.get_bonded_structure(structure)

        # Calculate Sij
        sij_solver = TheoreticalBondValenceSolver(element2charge)
        sij_results = sij_solver.get_sij("temp", structure, bonded_structure)
        dict_sij = sij_results[0]

        # Save Sij to temp file
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
            json.dump(dict_sij, f)
            sij_path = f.name

        # Calculate R0/B for each algorithm
        algos = ["shgo", "brute", "diff", "dual_annealing", "direct"]
        R0_vals, B_vals = [], []

        for algo in algos:
            with tempfile.TemporaryDirectory() as tmpdir:
                solver = BVParamSolver(tmpdir, algo)
                try:
                    R0, B = solver.solve_R0Bs(
                        cation=cation,
                        anion=anion,
                        bond_type_list=sij_results[1],
                        networkValence_dict=sij_results[0],
                        materID="mp-6340",
                        bondLen_dict=sij_results[2],
                        chem_formula="LiAlSi2O6",
                        R0_bounds=(0, 5),
                    )
                    R0_vals.append(R0)
                    B_vals.append(B)
                except Exception:
                    continue

        # Check consistency
        if len(R0_vals) == 5 and len(B_vals) == 5:
            R0_avg = sum(R0_vals) / 5
            B_avg = sum(B_vals) / 5
            if all(abs(r - R0_avg) < 1e-2 for r in R0_vals) and all(
                abs(b - B_avg) < 1e-2 for b in B_vals
            ):
                return sij_path, f"R\u2080: {R0_avg:.4f} \u00c5\nB: {B_avg:.4f}"

        return sij_path, "No consistent R\u2080/B values found"

    except Exception as e:
        return None, f"Error: {str(e)}"


interface = gr.Interface(
    fn=calculate_bv_params,
    inputs=[
        gr.Textbox(label="Cation", placeholder="Enter Cation (e.g., Li)"),
        gr.Textbox(label="Anion", placeholder="Enter Anion (e.g., O)"),
        gr.File(label="CIF File", file_types=[".cif"]),
    ],
    outputs=[
        gr.File(label="Download Sij Results"),
        gr.Text(label="Bond Valence Parameters"),
    ],
    title="Bond Valence Parameter Search",
    description="Upload a CIF file to calculate bond valence parameters (R\u2080 and B).",
)

interface.launch()
