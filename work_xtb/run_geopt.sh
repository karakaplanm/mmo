#!/bin/bash
# ============================================================
#  xTB Geometry Optimization Example — Ethanol (C₂H₅OH)
# ============================================================
#
# This script demonstrates several geometry optimization
# scenarios using GFN2-xTB, from a basic gas-phase run to
# solvated and frequency-verified optimizations.
#
# Output files produced by xtb:
#   xtbopt.xyz   – optimized geometry (Cartesian coordinates)
#   xtbopt.log   – optimization trajectory
#   xtbout       – main output (energy, gradient, etc.)
#   wbo          – Wiberg bond orders
#   charges      – atomic partial charges
#   xtbrestart   – restart information
# ============================================================

XTB="/home/muka/Downloads/xtb-dist/bin/xtb"
MOL="ethanol.xyz"

echo "========================================"
echo " 1) Basic GFN2-xTB Geometry Optimization"
echo "========================================"
# --opt         : run geometry optimization (default level = normal)
# --gfn 2       : use GFN2-xTB parametrization (default, shown explicitly)
# -c 0          : neutral molecule
# --uhf 0       : closed-shell (no unpaired electrons)
$XTB "$MOL" --opt --gfn 2 -c 0 --uhf 0 2>&1 | tee opt_basic.log

# Save optimized structure
cp xtbopt.xyz ethanol_opt_basic.xyz
echo ""
echo ">>> Optimized geometry saved to ethanol_opt_basic.xyz"
echo ""

echo "========================================"
echo " 2) Tight Optimization + Frequency Check"
echo "========================================"
# --ohess tight : optimize at 'tight' level, then compute Hessian
#                 to verify the structure is a true minimum
#                 (no imaginary frequencies)
$XTB "$MOL" --ohess tight --gfn 2 2>&1 | tee opt_ohess.log

cp xtbopt.xyz ethanol_opt_tight.xyz
echo ""
echo ">>> Optimized geometry saved to ethanol_opt_tight.xyz"
echo ">>> Check opt_ohess.log for vibrational frequencies."
echo ""

echo "========================================"
echo " 3) Optimization in Water Solvent (ALPB)"
echo "========================================"
# --alpb water  : use the Analytical Linearized Poisson-Boltzmann
#                 implicit solvation model with water as solvent
$XTB "$MOL" --opt --gfn 2 --alpb water 2>&1 | tee opt_solvent.log

cp xtbopt.xyz ethanol_opt_water.xyz
echo ""
echo ">>> Optimized geometry (in water) saved to ethanol_opt_water.xyz"
echo ""

echo "========================================"
echo " All done!"
echo "========================================"
echo " Key output files:"
echo "   ethanol_opt_basic.xyz  — gas-phase optimized geometry"
echo "   ethanol_opt_tight.xyz  — tight optimization + frequencies"
echo "   ethanol_opt_water.xyz  — ALPB water-solvated optimization"
echo ""
echo " Inspect any .log file for full energies and convergence details."
