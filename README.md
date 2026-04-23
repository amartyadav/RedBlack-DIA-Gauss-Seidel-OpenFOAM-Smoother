[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.19706030.svg)](https://doi.org/10.5281/zenodo.19706030)

# RBDIAGaussSeidel

A drop-in Red-Black Gauss-Seidel smoother plugin for **OpenFOAM Foundation v13** that operates
on matrices stored in **Diagonal (DIA) sparse format**, targeting structured hexahedral
meshes for improved memory access performance and data-race-free parallelism.

Selectable at runtime via `fvSolution` — no source modification required.

```
solver          p
{
    solver          PCG;
    preconditioner  DIC;
    smoother        RBDIAGaussSeidel;   // <-- drop this in
    ...
}
```

---

## Background

OpenFOAM's default Gauss-Seidel smoother operates over the indirect LDU (Lower-Upper
Diagonal) sparse format, which involves pointer indirection at every matrix-vector
multiply. On structured hexahedral meshes, the non-zero sparsity pattern is fixed and
predictable — the DIA format exploits this by encoding off-diagonal bands as contiguous
arrays, eliminating indirection and dramatically improving cache line utilisation.

Red-Black ordering partitions the mesh cells into two independent sets — **red** cells
and **black** cells — such that no cell in one set has a neighbour in the same set.
This means all red cells can be updated simultaneously without data races, followed by
all black cells. The result is a smoother that is both cache-friendly and safe for
shared-memory parallelism (OpenMP), without sacrificing convergence behaviour relative
to the standard lexicographic Gauss-Seidel.

**Key characteristics of the target problem:**
- Structured hex meshes (e.g. cavity, duct, channel flow)
- Symmetric matrices with uniform or variable coefficients per diagonal band
- Fixed diagonal offsets: `{±1, ±N, ±N²}` for a mesh of stride `N`
- Red-Black cell ordering derived from mesh index parity

---

## Scope and Limitations

This plugin is intentionally scoped. It is **not** a general-purpose smoother.

- ✅ Structured hexahedral meshes
- ✅ Symmetric matrices
- ✅ Uniform and variable coefficients per diagonal band
- ✅ Data-race-free parallel sweep (OpenMP-ready)
- ✅ OpenFOAM Foundation v13
- ❌ Unstructured or polyhedral meshes (falls back gracefully — see below)
- ❌ ESI/OpenCFD (foam-extend) variants — untested

The plugin includes a structured-mesh detection check at initialisation. If the mesh
does not satisfy the DIA assumptions, it will either warn and fall back or abort
cleanly, depending on configuration.

---

## Relation to DIAGaussSeidel

This plugin is the Red-Black extension of
[DIAGaussSeidel](https://github.com/amartyadav/DIAGaussSeidel-Smoother-OpenFOAM)
(DOI: [10.5281/zenodo.19706068](https://doi.org/10.5281/zenodo.19706068)).

Both share the same DIA matrix format and structured-mesh assumptions. The difference
is the sweep order: `DIAGaussSeidel` uses standard lexicographic ordering (sequential,
cache-efficient); `RBDIAGaussSeidel` uses Red-Black ordering (parallel-safe, identical
convergence on symmetric problems). For single-threaded use the two are comparable;
Red-Black becomes advantageous when OpenMP threading is enabled.

---

## Build

Standard OpenFOAM wmake workflow:

```bash
source /opt/openfoam13/etc/bashrc   # or your OF installation path
wmake libso
```

This produces `libRBDIAGaussSeidel.so` in `$FOAM_USER_LIBBIN`.

---

## Usage

**1. Load the library in `system/controlDict`:**

```cpp
libs
(
    "libRBDIAGaussSeidel.so"
);
```

**2. Select the smoother in `system/fvSolution`:**

```cpp
solvers
{
    p
    {
        solver          PCG;
        preconditioner  DIC;
        smoother        RBDIAGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }
}
```

The smoother name `RBDIAGaussSeidel` is registered in OpenFOAM's runtime selection
table and will be picked up automatically.

---

## Repository Structure

```
RBDIAGaussSeidel/
├── Make/
│   ├── files
│   └── options
├── include/
│   └── RBDIAGaussSeidelSmoother.H
├── src/
│   └── RBDIAGaussSeidelSmoother.C
├── CITATION.cff
├── LICENSE
└── README.md
```

---

## Citation

If this work contributes to a publication, conference paper, or technical report,
please cite it. GitHub will surface a **Cite this repository** button using the
included `CITATION.cff`, or you can cite manually:

```
Amartya Yadav. RBDIAGaussSeidel: A Red-Black DIA-format Gauss-Seidel Smoother Plugin
for OpenFOAM. 2026. https://github.com/amartyadav/RBDIAGaussSeidel-Smoother-OpenFOAM
```

BibTeX:
```bibtex
@software{yadav2026rbdiagaussseidel,
  author  = {Yadav, Amartya},
  title   = {{RBDIAGaussSeidel}: A Red-Black {DIA}-format {Gauss-Seidel} Smoother Plugin for {OpenFOAM}},
  year    = {2026},
  url     = {https://github.com/amartyadav/RBDIAGaussSeidel-Smoother-OpenFOAM},
  license = {MIT}
}
```

---

## License

MIT License — © 2026 Amartya Yadav.

You are free to use, modify, and distribute this software for any purpose. The
copyright notice and this permission notice must be preserved in all copies or
substantial portions of the software. See [LICENSE](LICENSE) for full terms.

---

## Author

**Amartya Yadav**  
HPC Software Engineer, Calligo Technologies  
MSc High Performance Computing, EPCC, University of Edinburgh  
[github.com/amartyadav](https://github.com/amartyadav)  
[amartyadav.com](https://amartyadav.com)
