#!/bin/bash


# Check existence of aliBuild
hash aliBuild && FOUND_ALIBUILD="1" || FOUND_ALIBUILD=""

# Check existence of alidist in this directory
FOUND_ALIDIST="$(find . -name alidist -type -d)"

# Install aliBuild if necessary
[[ "$FOUND_ALIBUILD" == "" ]] && pip install --user alibuild

[[ "$FOUND_ALIDIST" == "" ]] && git clone https://github.com/benedikt-voelkel/alidist.git
cd alidist
git checkout concurrentEngines
cd ..

echo -e "Configuring done...\n"
