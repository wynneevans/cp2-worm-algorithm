#!/bin/zsh
# Usage: scripts/run_worm.zsh   (override with CP2_BIN=... CP2_OUT=...)
root=${0:A:h:h}
BIN=${CP2_BIN:-$root/build/release/src/wormCP2Mu2dAFMplaq}
OUT=${CP2_OUT:-$root/data}

 #"$BIN" "$OUT" closed 100000 10.0 20 4 800 probabilities.txt 1 1

 "$BIN" "$OUT" closed 100000 2.50 4 4 200 probabilities.txt 1 1 $((RANDOM + 1))
 "$BIN" "$OUT" closed 100000 5.0 20 4 400 probabilities.txt 1 1 $((RANDOM + 1))
 "$BIN" "$OUT" closed 100000 5.0 8 8 400 probabilities.txt 1 1 $((RANDOM + 1))
  "$BIN" "$OUT" closed 100000 5.0 8 4 400 probabilities.txt 1 1 $((RANDOM + 1))
 "$BIN" "$OUT" closed 100000 10.0 8 8 800 probabilities.txt 1 1 $((RANDOM + 1))
 "$BIN" "$OUT" closed 100000 20.0 4 4 1600 probabilities.txt 1 1 $((RANDOM + 1))

 #"$BIN" "$OUT" y-open 100000 100.0 10 2 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" y-open 100000 100.0 20 4 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" y-open 100000 100.0 32 6 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" y-open 100000 100.0 56 8 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" y-open 100000 100.0 90 10 8000 probabilities.txt 1 1

 #"$BIN" "$OUT" closed 100000 100.0 10 2 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" closed 100000 100.0 20 4 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" closed 100000 100.0 32 6 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" closed 100000 100.0 56 8 8000 probabilities.txt 1 1
 #"$BIN" "$OUT" closed 100000 100.0 90 10 8000 probabilities.txt 1 1
 