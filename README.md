# ARF optimization artifact

This source-only artifact provides exact and certified-gap anytime ARF
optimization. Both modes use the same executable.
It includes full rank-offset accounting. Benchmark orchestration, results,
historical variants and experimental datasets are not included.

## Build requirements

Linux x86-64, GNU Make, a C++17 compiler and OpenMP are required, together with:

- IBM ILOG CPLEX Optimization Studio 22.1.1 (CPLEX and Concert development files).
  Obtain and install it separately under an appropriate IBM license.
- Google OR-Tools C++ 9.11, Protobuf and compatible Abseil libraries.
- Reentrant Qhull development library (`libqhull_r`).
- GMP C++ development library (`libgmpxx`).
- lp_solve 5.5.2.0 library, matching the bundled third-party headers.

Dependency binaries, proprietary installers and license keys are not supplied.
Set library search paths for your installation when necessary.

```bash
make -j2 CPLEX_DIR=/path/to/CPLEX_Studio2211 PREFIX=/usr/local
```

If lp_solve is not installed in the linker search path:

```bash
make -j2 CPLEX_DIR=/path/to/CPLEX_Studio2211 \
  LP_SOLVE_LIB=/path/to/liblpsolve55.a
```

The executable is `bin/arf`. `make clean` removes local build outputs only.
Optimization flags are configurable with `CXXFLAGS`; timings depend on the
compiler flags, dependencies and hardware used.

## Run the synthetic example

From the repository root, exact mode:

```bash
OMP_NUM_THREADS=1 ANYTIME_DISABLE=1 ./bin/arf \
  -d 2 -pf examples/products.dat -uf examples/users.dat \
  -wca generated -wr 1.0 -wn 4 -ws 31
```

Anytime mode with a 0.5% certified relative-gap target:

```bash
env -u ANYTIME_DISABLE -u ANYTIME_STOP_GAP_PERCENT OMP_NUM_THREADS=1 \
  ./bin/arf -d 2 -pf examples/products.dat -uf examples/users.dat \
  -wca generated -wr 1.0 -wn 4 -ws 31 -agap 0.005
```

Use `-agap 0.01` for 1% and `-agap 0.02` for 2%. In this interface,
`-agap 1` means 100%, not 1%. Clear inherited anytime environment variables.
The relative gap is `(incumbent - certified_lower_bound) / certified_lower_bound`
when the certified lower bound is positive. This is a certificate against a
bound, not a measurement against an independently known optimum.

The program prints `ARF:` and `TIME_SEC:`. Set `ANYTIME_INTERVAL_SEC=5` for
additional snapshots at solver checkpoints; this is not a hard real-time timer.

## Input format

Whitespace-separated text with no header:

- Products: integer ID, then `d` lower coordinates, then `d` upper coordinates.
- Users: integer ID followed by `d` preference coordinates.

Dimensions must match `-d`. The bundled example is synthetic demonstration
data, not an experimental dataset. `-ws` selects the generated WCA seed;
the user seed is determined when creating the user input file.


## Third-party notices

See `THIRD_PARTY_NOTICES.md`. Preserve the notices in the bundled headers.
