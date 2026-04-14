# DIAGaussSeidel

A drop-in Gauss-Seidel smoother plugin for **OpenFOAM Foundation v13** that operates
on matrices stored in **Diagonal (DIA) sparse format**, targeting structured hexahedral
meshes for improved memory access performance.

Selectable at runtime via `fvSolution` — no source modification required.

```
solver          p
{
    solver          PCG;
    preconditioner  DIC;
    smoother        DIAGaussSeidel;   // <-- drop this in
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

This plugin replaces the LDU smoother with a DIA-format implementation for the specific
case where the mesh topology permits it. The result is fewer DRAM transactions, lower
instruction counts, and measurable wall-clock improvements on modern out-of-order CPUs.

**Key characteristics of the target problem:**
- Structured hex meshes (e.g. cavity, duct, channel flow)
- Symmetric matrices ~~with uniform coefficients per diagonal band~~ (now supports uniform or variable coefficients per diagonal band - should cover more use cases now)
- Fixed diagonal offsets: `{±1, ±N, ±N²}` for a mesh of stride `N`

---

## Scope and Limitations

This plugin is intentionally scoped. It is **not** a general-purpose smoother.

- ✅ Structured hexahedral meshes
- ✅ Symmetric matrices
- ✅ ~~Uniform coefficients per diagonal band~~ Now supports both uniform and variable coefficients
- ✅ OpenFOAM Foundation v13
- ❌ Unstructured or polyhedral meshes (falls back gracefully — see below)
- ❌ ~~Non-uniform or anisotropic coefficient distributions~~
- ❌ ESI/OpenCFD (foam-extend) variants — untested

The plugin includes a structured-mesh detection check at initialisation. If the mesh
does not satisfy the DIA assumptions, it will either warn and fall back or abort
cleanly, depending on configuration.

---

## Build

Standard OpenFOAM wmake workflow:

```bash
source /opt/openfoam13/etc/bashrc   # or your OF installation path
wmake libso
```

This produces `libDIAGaussSeidel.so` in `$FOAM_USER_LIBBIN`.

---

## Usage

**1. Load the library in `system/controlDict`:**

```cpp
libs
(
    "libDIAGaussSeidel.so"
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
        smoother        DIAGaussSeidel;
        tolerance       1e-6;
        relTol          0.1;
    }
}
```

The smoother name `DIAGaussSeidel` is registered in OpenFOAM's runtime selection
table and will be picked up automatically.

---

## Repository Structure

```
DIAGaussSeidel/
├── Make/
│   ├── files
│   └── options
├── include/
│   └── DIAGaussSeidelSmoother.H
├── src/
│   └── DIAGaussSeidelSmoother.C
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
Amartya Yadav. DIAGaussSeidel: A DIA-format Gauss-Seidel Smoother Plugin
for OpenFOAM. 2026. https://github.com/amartyadav/DIAGaussSeidel-Smoother-OpenFOAM
```

BibTeX:
```bibtex
@software{yadav2026diagaussseidel,
  author  = {Yadav, Amartya},
  title   = {{DIAGaussSeidel}: A {DIA}-format {Gauss-Seidel} Smoother Plugin for {OpenFOAM}},
  year    = {2026},
  url     = {https://github.com/amartyadav/DIAGaussSeidel-Smoother-OpenFOAM},
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
HPC Software Engineer, Calligo Technologies <br/>
MSc High Performance Computing, EPCC, University of Edinburgh  
[github.com/\amartyadav](https://github.com/amartyadav) <br/>
[amartyadav.com/](https://amartyadav.com)
