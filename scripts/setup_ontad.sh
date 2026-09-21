#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
commit=dc36df93598bd2dce92fe6e1e97425a60c73178b
mkdir -p tools
if [ ! -d tools/OnTAD/.git ]; then
  git clone https://github.com/anlin00007/OnTAD.git tools/OnTAD
fi
git -C tools/OnTAD checkout --detach "$commit"
cd tools/OnTAD/src
"${CXX:-c++}" -std=c++11 main.cpp step1.cpp step2.cpp step3.cpp step4.cpp common.cpp straw.cpp -lm -lcurl -lz -o OnTAD
