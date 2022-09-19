#!/usr/bin/env bash

set -euo pipefail

bash test/treeknit_tests.sh
OUT=$?
if [ "$OUT" != 0 ]; then
  exit 1
fi