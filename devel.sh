#!/usr/bin/env bash

SCRIPT="
using AuChimiste;
using ChemicalComponents;
using ChemicalKinetics;
using ChemicalReactors;
using PhysicalChemistry;
using CombustionChemistry;
"

julia --project -e "${SCRIPT}" -i
