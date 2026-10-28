#!/usr/bin/env python3
"""Compile the production scale-weight readers against small ROOT-free stubs."""

from pathlib import Path
import os
import shlex
import subprocess
import sys
import tempfile


def function_source(path, signature):
    text = path.read_text()
    start = text.index(signature)
    opening = text.index("{", start)
    depth = 1
    end = opening + 1
    while depth:
        depth += (text[end] == "{") - (text[end] == "}")
        end += 1
    return text[start:end]


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
    core = function_source(root / "Analyzers/src/AnalyzerCore.cc", "float AnalyzerCore::GetScaleVariation(")
    niso = function_source(root / "Analyzers/src/NIsoMuon.cc", "float NIsoMuon::ReadTheoryWeight(")
    harness = r'''
#include <array>
#include <cassert>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using std::cout;
using std::endl;
using TString = std::string;
namespace MyCorrection { enum class variation { nom, up, down }; }
using V = MyCorrection::variation;
struct AnalyzerCore {
    unsigned int nLHEScaleWeight = 0;
    std::vector<float> LHEScaleWeight;
    float GetScaleVariation(const V &, const V &);
};
struct NIsoMuon : AnalyzerCore {
    int nLHEPdfWeight = 0;
    std::vector<float> LHEPdfWeight;
    float ReadTheoryWeight(const TString &, int) const;
};
'''
    harness += core + "\n" + niso + r'''
int main() {
    // Independent physical grid: rows are muR=(0.5,1,2), columns muF=(0.5,1,2).
    // Distinct marker values make wrong source indices observable.
    const std::vector<float> nine = {11,12,13,21,1,23,31,32,33};
    const std::vector<float> eight = {11,12,13,21,23,31,32,33};
    const std::array<float,9> expected = {1,23,21,32,33,31,12,13,11};
    const std::array<V,9> muF = {V::nom,V::up,V::down,V::nom,V::up,V::down,V::nom,V::up,V::down};
    const std::array<V,9> muR = {V::nom,V::nom,V::nom,V::up,V::up,V::up,V::down,V::down,V::down};
    for (const auto &grid : {eight, nine}) {
        NIsoMuon reader;
        reader.nLHEScaleWeight = grid.size();
        reader.LHEScaleWeight = grid;
        for (int i=0; i<9; ++i) {
            const float a = reader.GetScaleVariation(muF[i],muR[i]);
            const float b = reader.ReadTheoryWeight("PDFScale",i);
            if (a != expected[i] || b != expected[i]) {
                std::cerr << grid.size() << " entries, PDFScale" << i
                          << ": expected " << expected[i] << ", core=" << a
                          << ", NIsoMuon=" << b << std::endl;
                return 1;
            }
        }
    }
    NIsoMuon missing;
    for (int i=0; i<9; ++i) {
        assert(missing.GetScaleVariation(muF[i],muR[i]) == 1);
        assert(missing.ReadTheoryWeight("PDFScale",i) == 1);
    }
    assert(std::isnan(missing.ReadTheoryWeight("PDFScale",-1)));
    assert(std::isnan(missing.ReadTheoryWeight("PDFScale",9)));
    NIsoMuon malformed;
    malformed.nLHEScaleWeight = 10;
    malformed.LHEScaleWeight.resize(10,1);
    bool coreThrows=false, nisoThrows=false;
    try { malformed.GetScaleVariation(V::up,V::nom); }
    catch (const std::runtime_error &) { coreThrows=true; }
    try { malformed.ReadTheoryWeight("PDFScale",1); }
    catch (const std::runtime_error &) { nisoThrows=true; }
    assert(coreThrows && nisoThrows);
    NIsoMuon nonfinite;
    nonfinite.nLHEScaleWeight=9;
    nonfinite.LHEScaleWeight=nine;
    nonfinite.LHEScaleWeight[5]=std::numeric_limits<float>::quiet_NaN();
    assert(nonfinite.ReadTheoryWeight("PDFScale",1)==1);
    std::cout << "PASS: both production readers, all nine points in 8/9-entry grids; "
                 "missing, invalid-size, invalid-index and non-finite handling." << std::endl;
}
'''
    with tempfile.TemporaryDirectory(prefix="lhe_scale_test_") as directory:
        source = Path(directory) / "test.cc"
        binary = Path(directory) / "test"
        source.write_text(harness)
        compiler = shlex.split(os.environ.get("CXX", "c++"))
        subprocess.run(compiler + ["-std=c++17", "-Wall", "-Wextra", "-Werror", str(source), "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True)


if __name__ == "__main__":
    main()
