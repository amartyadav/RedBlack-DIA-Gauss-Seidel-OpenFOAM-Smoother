#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_START(a)
#define LIKWID_MARKER_STOP(a)
#define LIKWID_MARKER_CLOSE
#endif

#ifdef LIKWID_PERFMON
__attribute__((constructor))
static void likwid_init_wrapper()
{
    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;
}

__attribute__((destructor))
static void likwid_close_wrapper()
{
    LIKWID_MARKER_CLOSE;
}
#endif

#include "RBDIAGaussSeidelSmoother.H"
#include "className.H"
#include "label.H"
#include "lduInterfaceFieldPtrsList.H"
#include "lduMatrix.H"
#include "scalar.H"
#include "scalarField.H"
#include <set>
#include <vector>
#include "GaussSeidelSmoother.H"
#include <chrono>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(RBDIAGaussSeidelSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<RBDIAGaussSeidelSmoother>
        addRBDIAGaussSeidelSmootherSymMatrixConstructorToTable_;

    lduMatrix::smoother::addasymMatrixConstructorToTable<RBDIAGaussSeidelSmoother>
        addRBDIAGaussSeidelSmootherAsymConstructorToTable_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBDIAGaussSeidelSmoother::RBDIAGaussSeidelSmoother(
    const word &fieldName,
    const lduMatrix &matrix,
    const FieldField<Field, scalar> &interfaceBouCoeffs,
    const FieldField<Field, scalar> &interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList &interfaces)
    : lduMatrix::smoother(
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    structuredMesh_(false),
    Nx_(0),
    Ny_(0),
    Nz_(0),
    useDIA_(false),
    upperCoeffI_(0.0),
    upperCoeffJ_(0.0),
    upperCoeffK_(0.0),
    lowerCoeffI_(0.0),
    lowerCoeffJ_(0.0),
    lowerCoeffK_(0.0),
    upperCoeffFieldI_(),
    upperCoeffFieldJ_(),
    upperCoeffFieldK_()
{
    const label nCells = matrix_.diag().size();

    const labelUList upperAddr = matrix_.lduAddr().upperAddr();
    const labelUList lowerAddr = matrix_.lduAddr().lowerAddr();
    label nFaces = upperAddr.size();
    const scalarField& upper = matrix_.upper();



    // unique offsets
    std::set<label> offsets;
    for(label f = 0; f < lowerAddr.size(); f++)
    {
        offsets.insert(upperAddr[f] - lowerAddr[f]);
    }
    std::vector<label> sortedOffsets(offsets.begin(), offsets.end());

    // checking if the size of the set (offsets) is 2 (for 2d) or 3 (for 3d)
    if(sortedOffsets.size() == 2)
    {
        if(sortedOffsets[0] == 1 && nCells % sortedOffsets[1] == 0)
        {
            Nx_ = sortedOffsets[1];
            Ny_ = nCells / Nx_;
            Nz_ = 1;
            structuredMesh_ = true;
        }
    }
    else if(sortedOffsets.size() == 3)
    {
        if(sortedOffsets[0] == 1 && sortedOffsets[2] % sortedOffsets[1] == 0
            && nCells % sortedOffsets[2] == 0)
        {
            Nx_ = sortedOffsets[1];
            Ny_ = sortedOffsets[2] / sortedOffsets[1];
            Nz_ = nCells / sortedOffsets[2];
            structuredMesh_ = true;
        }
    }

    // Logging the outcome
    // DEBUG: dump first 5 faces per direction to see actual coefficient values
    if (structuredMesh_)
    {
        label countI = 0, countJ = 0, countK = 0;
        // Info<< "DEBUG: sampling upper coefficients by direction" << endl;
        for (label f = 0; f < upperAddr.size() && (countI < 5 || countJ < 5 || countK < 5); f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            if (offset == 1 && countI < 5)
            {
                countI++;
            }
            else if (offset == Nx_ && countJ < 5)
            {
                countJ++;
            }
            else if (offset == Nx_ * Ny_ && countK < 5)
            {
                countK++;
            }
        }
    }
    // checking coefficient uniformity
    bool allUniform = false;
    if(structuredMesh_ && !matrix_.asymmetric())
    {
        bool firstI = true, firstJ = true, firstK = true;

        allUniform = true;

        for(label f = 0; f < nFaces; f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            const scalar coeff = upper[f];

            if(offset == 1) // I
            {
                if(firstI)
                {
                    upperCoeffI_ = coeff;
                    firstI = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffI_), SMALL);
                    if(std::abs(coeff - upperCoeffI_) > tol)
                    {

                        allUniform = false;
                        break;
                    }
                }
            }
            else if(offset == Nx_) // J
            {
                if(firstJ)
                {
                    upperCoeffJ_ = coeff;
                    firstJ = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffJ_), SMALL);
                    if(std::abs(coeff - upperCoeffJ_) > tol)
                    {

                        allUniform = false;
                        break;
                    }
                }
            }
            else if(offset == Nx_*Ny_) // K
            {
                if(firstK)
                {
                    upperCoeffK_ = coeff;
                    firstK = false;
                }
                else
                {
                    const scalar tol = 1e-12 * std::max(std::abs(upperCoeffK_), SMALL);
                    if(std::abs(coeff - upperCoeffK_) > tol)
                    {

                        allUniform = false;
                        break;
                    }
                }
            }
        }
        if(allUniform)
        {
            lowerCoeffI_ = upperCoeffI_;

            lowerCoeffJ_ = upperCoeffJ_;
            lowerCoeffK_ = upperCoeffK_;
        }
    }

    useDIA_ = structuredMesh_ && !matrix.asymmetric();

    if (useDIA_)
    {
        if (allUniform)
        {
            Info<< "RBDIAGaussSeidel: fast path enabled (uniform coefficients),  Nx=" << Nx_
                << " Ny=" << Ny_ << " Nz=" << Nz_
                << " upperCoeffs=(" << upperCoeffI_
                << ", " << upperCoeffJ_
                << ", " << upperCoeffK_ << ")"
                << " (nCells=" << nCells << ")" << endl;
        }
        else
        {
            Info<< "RBDIAGaussSeidel: fast path enabled (variable coefficients),  Nx=" << Nx_
                << " Ny=" << Ny_ << " Nz=" << Nz_
                << " upperCoeffs=(" << upperCoeffI_
                << ", " << upperCoeffJ_
                << ", " << upperCoeffK_ << ")"
                << " (nCells=" << nCells << ")" << endl;
        }

        // initialising coefficient field
        upperCoeffFieldI_.setSize(nCells);
        upperCoeffFieldJ_.setSize(nCells);
        upperCoeffFieldK_.setSize(nCells);

        upperCoeffFieldI_ = 0.0;
        upperCoeffFieldJ_ = 0.0;
        upperCoeffFieldK_ = 0.0;

        for(label f = 0; f < nFaces; f++)
        {
            const label offset = upperAddr[f] - lowerAddr[f];
            const label owner = lowerAddr[f];

            if(offset == 1) // I
            {
                upperCoeffFieldI_[owner] = upper[f];
            }
            else if (offset == Nx_) // J
            {
                upperCoeffFieldJ_[owner] = upper[f];
            }
            else if(offset == Nx_*Ny_) // K
            {
                upperCoeffFieldK_[owner] = upper[f];
            }
            else{}
        }

    }
    else if (!structuredMesh_)
    {
        Info<< "RBDIAGaussSeidel: mesh not structured, falling back (nCells="
            << nCells << ")" << endl;
    }
    else if (matrix_.asymmetric())
    {
        Info<< "RBDIAGaussSeidel: asymmetric matrix, falling back" << endl;
    }

    if (useDIA_ && nCells == 125000)  // only log on finest level of your test case
    {
        Info<< "DEBUG: first 5 entries of upperCoeffFieldI_: ";
        for (label i = 0; i < 5; i++) Info<< upperCoeffFieldI_[i] << " ";
        Info<< endl;
        // same for J and K
        Info<< "DEBUG: first 5 entries of upperCoeffFieldJ_: ";
        for (label i = 0; i < 5; i++) Info<< upperCoeffFieldJ_[i] << " ";
        Info<< endl;
        Info<< "DEBUG: first 5 entries of upperCoeffFieldK_: ";
        for (label i = 0; i < 5; i++) Info<< upperCoeffFieldK_[i] << " ";
        Info<< endl;

        Info<< "DEBUG: boundary slots (should be 0): "
            << "upperCoeffFieldI_[" << (Nx_-1) << "]=" << upperCoeffFieldI_[Nx_-1] << " "
            << "upperCoeffFieldJ_[" << (Nx_*(Ny_-1)) << "]=" << upperCoeffFieldJ_[Nx_*(Ny_-1)] << " "
            << "upperCoeffFieldK_[" << (nCells-1) << "]=" << upperCoeffFieldK_[nCells-1]
            << endl;
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Define the virtual smooth method (the 4-arg const one). Its body should just call Foam::GaussSeidelSmoother::smooth(...) — the stock one — passing fieldName_, psi, matrix_, source, interfaceBouCoeffs_, interfaces_, cmpt, nSweeps. All the underscore-suffixed members are inherited from the protected base class.
void Foam::RBDIAGaussSeidelSmoother::smooth(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    if (useDIA_)
    {

        scalar* __restrict__ psiPtr = psi.begin();
        const scalar* __restrict__ diagPtr = matrix_.diag().begin();
        const label nCells = psi.size();

        scalarField bPrime(nCells);
        scalar* bPrimePtr = bPrime.begin();

        // Cache Nx, Ny, Nz into local labels (avoids repeated member access in hot loop)
        label Nx = Nx_;
        label Ny = Ny_;
        label Nz = Nz_;

        const scalar* __restrict__ upperIPtr = upperCoeffFieldI_.begin();
        const scalar* __restrict__ upperJPtr = upperCoeffFieldJ_.begin();
        const scalar* __restrict__ upperKPtr = upperCoeffFieldK_.begin();

        const scalar* __restrict__ lowerIPtr = upperIPtr;
        const scalar* __restrict__ lowerJPtr = upperJPtr;
        const scalar* __restrict__ lowerKPtr = upperKPtr;


        label jStride = Nx;
        label kStride = Nx * Ny;

        // Parallel boundary initialisation.  The parallel boundary is treated
            // as an effective jacobi interface in the boundary.
            // Note: there is a change of sign in the coupled
            // interface update.  The reason for this is that the
            // internal coefficients are all located at the l.h.s. of
            // the matrix whereas the "implicit" coefficients on the
            // coupled boundaries are all created as if the
            // coefficient contribution is of a source-kind (i.e. they
            // have a sign as if they are on the r.h.s. of the matrix.
            // To compensate for this, it is necessary to turn the
            // sign of the contribution.

            FieldField<Field, scalar> &mBouCoeffs =
                const_cast<FieldField<Field, scalar> &>(
                    interfaceBouCoeffs_);

            forAll(mBouCoeffs, patchi)
            {
                if (interfaces_.set(patchi))
                {
                    mBouCoeffs[patchi].negate();
                }
            }

        // ---- Sweep loop ----
        static double totalSweepTime = 0.0;
        static label totalCalls = 0;
        auto t0 = std::chrono::high_resolution_clock::now();
		LIKWID_MARKER_START("RBDIA_sweep");
        for (label sweep = 0; sweep < nSweeps; sweep++)
        {
            bPrime = source;

            matrix_.initMatrixInterfaces(
                        mBouCoeffs,
                        interfaces_,
                        psi,
                        bPrime,
                        cmpt);

                    matrix_.updateMatrixInterfaces(
                        mBouCoeffs,
                        interfaces_,
                        psi,
                        bPrime,
                        cmpt);

            // red cells sweep
            #pragma omp parallel for
            for(label idx = 0; idx < nCells; idx++)
            {
                label j = idx % Nx;
                label i = (idx / Nx) % Ny;
                label k = idx / kStride;

                if((i + j + k) % 2 != 0) continue; // skip half the cells
                scalar psii = bPrimePtr[idx];

                // forward neighbours (coefficient is at upperXPtr[idx])
                if(j < Nx - 1)
                {
                    psii -= upperIPtr[idx] * psiPtr[idx + 1];

                }
                if(i < Ny - 1)
                {
                    psii -= upperJPtr[idx] * psiPtr[idx + jStride];
                }
                if(k < Nz - 1)
                {
                    psii -= upperKPtr[idx] * psiPtr[idx + kStride];
                }

                // backward neighbours (coefficient is at upperXPtr[idx - 1])
                if(j > 0)
                {
                    psii -= upperIPtr[idx - 1] * psiPtr[idx - 1];
                }
                if(i > 0)
                {
                    psii -= upperJPtr[idx - jStride] * psiPtr[idx - jStride];
                }
                if(k > 0)
                {
                    psii -= upperKPtr[idx - kStride] * psiPtr[idx - kStride];
                }

                psiPtr[idx] = psii / diagPtr[idx];

            }
            // black cells sweep
            #pragma omp parallel for
            for(label idx = 0; idx < nCells; idx++)
            {
                label j = idx % Nx;
                label i = (idx / Nx) % Ny;
                label k = idx / kStride;

                if((i + j + k) % 2 != 1) continue;
                scalar psii = bPrimePtr[idx];

                // forward neighbours (coefficient is at upperXPtr[idx])
                if(j < Nx - 1)
                {
                    psii -= upperIPtr[idx] * psiPtr[idx + 1];
                }
                if(i < Ny - 1)
                {
                    psii -= upperJPtr[idx] * psiPtr[idx + jStride];
                }
                if(k < Nz - 1)
                {
                    psii -= upperKPtr[idx] * psiPtr[idx + kStride];
                }

                // backward neighbours (coefficient is at upperXPtr[idx - 1])
                if(j > 0)
                {
                    psii -= upperIPtr[idx - 1] * psiPtr[idx - 1];
                }
                if(i > 0)
                {
                    psii -= upperJPtr[idx - jStride] * psiPtr[idx - jStride];
                }
                if(k > 0)
                {
                    psii -= upperKPtr[idx - kStride] * psiPtr[idx - kStride];
                }

                psiPtr[idx] = psii / diagPtr[idx];

            }
        }
        LIKWID_MARKER_STOP("RBDIA_sweep");
        auto t1 = std::chrono::high_resolution_clock::now();
        totalSweepTime += std::chrono::duration<double>(t1 - t0).count();
        totalCalls++;

        if (totalCalls % 1000 == 0)
        {
            Info<< "RBDIA_sweep: totalTime=" << totalSweepTime
                << "s calls=" << totalCalls << endl;
        }

        // ---- Restore interface coefficients ----
        forAll(mBouCoeffs, patchi)
            {
                if (interfaces_.set(patchi))
                {
                    mBouCoeffs[patchi].negate();
                }
            }
    }
    else
    {
        // fallback - existing forwarding call
        GaussSeidelSmoother::smooth(
            fieldName_,
            psi,
            matrix_,
            source,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt,
            nSweeps
        );
    }
}
