#!/bin/zsh
# Regression test: a fixed-seed run must reproduce the golden .dat files byte for byte.
# Regenerate the golden (only after an intended change to the physics):
#   ./build/wormCP2Mu2dAFMplaq test/golden closed 1000 2.50 4 4 200 probabilities.txt 1 1 42

set -e
cd "${0:A:h}/.."    # algorithm/

seed=${1:-42}
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

cp test/golden/probabilities.txt "$scratch/"
./build/wormCP2Mu2dAFMplaq "$scratch" closed 1000 2.50 4 4 200 probabilities.txt 1 1 $seed > /dev/null

for kind in C M; do
  f="${kind}_wormCP2AFMmu2dplaq_BC_closed_nMeas_1000_beta_2.50_XYT_4_4_200_mu3_0.000_mu8_0.000_nSim_1.dat"
  if ! cmp -s "test/golden/$f" "$scratch/$f"; then
    echo "FAIL: $f differs from golden"
    diff "test/golden/$f" "$scratch/$f" | head -5
    exit 1
  fi
done
echo "PASS: all .dat files match golden"
