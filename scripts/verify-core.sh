#!/bin/sh
# Focused regression certificates for the current paper and exported geometry.
set -eu
unset PYTHONOPTIMIZE
cd "$(dirname "$0")/.."
for certificate in \
  exact_slice \
  exact_field \
  exact_core_descent \
  exact_reach \
  exact_kauffman2 \
  exact_trunc_nbhd \
  exact_kkt_dagger \
  exact_second \
  exact_branch \
  exact_opening \
  exact_unique_h \
  exact_common_overlap \
  exact_global_core \
  exact_census_roots \
  exact_outline \
  exact_notch \
  n4_exact_kkt \
  n4_exact_cell \
  n4_cell_crosscheck \
  n4_exact_valid \
  n4_exact_opening \
  n4_exact_local \
  n4_krawczyk \
  n4_certify \
  check_n4_geometry \
  check_n4_rerun \
  test_guards
do
  uv run python "verify/$certificate.py"
done
uv run python verify/export_figure_geometry.py --check
uv run python verify/tref_render.py --check
uv run python verify/render_families.py --check
uv run python verify/make_fig_notation.py --check
uv run python verify/render_junction.py --check
uv run python verify/n4_pair.py --check
uv run python verify/test_n4_figure.py
printf "%s\n" "PASS: core certificates, mutation controls and figure checks"
