# Designing Products that Rise in the Rankings

Anonymous source artifact for Aggregate Rank Minimization (ARM).

## Paper overview

The paper studies how to design a new product that ranks favorably across users
with different preferences. Products have numerical attributes, and each user
scores a product using a non-negative linear preference vector. A product's rank
for a user is the number of existing products with a strictly higher score.
The aggregate rank function (ARF) sums those ranks across users.

ARM seeks a feasible product minimizing ARF, subject to linear constraints on
its attributes. The implementation uses geometric preprocessing, upper-boundary
trimming, slab filtering, LP lower bounds and tightening, branch-and-bound, and
reduced MILP subproblems. The same executable supports exact search and anytime
search with a certified relative-gap stopping target.

This repository contains ARM (Exact and anytime), the EAC and MILP competitor
implementations, build instructions, and small synthetic example inputs.
It does not contain benchmark runners, ablation
variants, experimental results, proprietary libraries, or the paper's large
experimental datasets. The examples demonstrate the input format and commands;
they are not a reproduction of the paper's experimental evaluation.

## Requirements

Use Linux x86-64 with GNU Make, a C++17 compiler, OpenMP, and these development
dependencies:

- IBM ILOG CPLEX Optimization Studio 22.1.1, including CPLEX and Concert.
  Obtain and install it separately under an appropriate IBM license.
- Google OR-Tools C++ 9.11, Protobuf, and compatible Abseil libraries.
- Reentrant Qhull (`libqhull_r`).
- GMP C++ (`libgmpxx`).
- lp_solve 5.5.2.0, matching the bundled third-party headers.

Dependency binaries, proprietary installers, and license keys are not included.

## Build

Run from the repository root. Replace the dependency placeholders with your own
installation locations:

```bash
make -j2 CPLEX_DIR=/path/to/CPLEX_Studio2211 PREFIX=/path/to/dependency-prefix
```

`PREFIX` defaults to `/usr/local`. If lp_solve is not in the linker search path,
provide its library explicitly:

```bash
make -j2 CPLEX_DIR=/path/to/CPLEX_Studio2211 \
  PREFIX=/path/to/dependency-prefix \
  LP_SOLVE_LIB=/path/to/liblpsolve55.a
```

The resulting executable is `bin/arf`. Configure your dynamic-library search
path if your shared dependencies are outside the system search paths:

```bash
export LD_LIBRARY_PATH=/path/to/dependency-prefix/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
```

`make clean` removes only the local `build/` and `bin/` output directories.
Compiler options can be overridden with `CXXFLAGS`.

## Run the IND example

The example has 64 synthetic products with four pseudo-random uniform attributes
in the unit interval and five illustrative user preference vectors. Point
products are represented by equal lower and upper coordinates. The generated
feasible region uses radius 1, four constraint hyperplanes, and WCA seed 31.
There is no application timeout in these commands.

The sample coordinates are written with enough decimal digits to round-trip the
32-bit values used internally. Both example commands return `ARF: 71`; elapsed
time depends on the machine.

### Exact mode

```bash
OMP_NUM_THREADS=1 ANYTIME_DISABLE=1 ./bin/arf \
  -d 4 \
  -pf examples/ind_products.dat \
  -uf examples/ind_users.dat \
  -wca generated -wr 1.0 -wn 4 -ws 31
```

### Anytime mode: 2% certified gap

```bash
env -u ANYTIME_DISABLE -u ANYTIME_STOP_GAP_PERCENT OMP_NUM_THREADS=1 \
  ./bin/arf -d 4 \
  -pf examples/ind_products.dat \
  -uf examples/ind_users.dat \
  -wca generated -wr 1.0 -wn 4 -ws 31 -agap 0.02
```

Use `-agap 0.005` for 0.5%, `-agap 0.01` for 1%, or `-agap 0.04` for 4%.
For values at most 1, `-agap` is a fraction: `-agap 1` means 100%, not 1%.
Clear inherited `ANYTIME_DISABLE` and `ANYTIME_STOP_GAP_PERCENT` settings when
running anytime mode. `ANYTIME_DISABLE=1` takes precedence over a gap target.

The relative-gap certificate is
`(incumbent - certified_lower_bound) / certified_lower_bound` when the certified
lower bound is positive. It is a guarantee based on a lower bound, not a
comparison against an independently known optimum. A requested gap is an upper
limit on the certificate; the returned solution can be better than that limit.

## Output

The program prints the selected slab combination, a feasible-point candidate,
`ARF:`, and `TIME_SEC:`. `TIME_SEC` is elapsed solver time in seconds. Anytime
mode also prints incumbent and lower-bound information at search checkpoints.
For periodic checkpoint snapshots, set `ANYTIME_INTERVAL_SEC=5`; this setting
does not impose a hard real-time timer or a runtime cutoff.

The example is a smoke test, not a general correctness proof. Independently
reevaluate a returned point against the original inputs when using it as a
solution certificate.

## Run other inputs

Use the same commands with your own product and user files and set `-d` to their
dimension. Both formats are whitespace-separated text with no header.

Products contain an integer ID, then `d` lower coordinates, then `d` upper
coordinates. Non-point products are scored using their bounding-box centers.
A four-dimensional point product can be written as:

```text
1 0.2 0.4 0.6 0.8 0.2 0.4 0.6 0.8
```

Users contain an integer ID followed by `d` non-negative preference weights
summing to 1:

```text
1 0.25 0.25 0.25 0.25
```

Input precision matters near score ties: product coordinates are stored as
32-bit floating-point values, and user weights are truncated downward to four
decimal places on loading. Independent rank checks must apply the same input
conversion; using the original decimal values can change a near-tie comparison.
For non-point products, the lower-plus-upper coordinate addition used to compute
each center is also performed in 32-bit floating-point arithmetic.
The supplied IND inputs avoid additional conversion differences.

All product attributes in these examples are normalized to the unit interval,
with higher values representing more desirable attributes. The option `-ws`
selects the generated WCA seed, not a user seed; the user preferences come from
the supplied user file. The original two-dimensional smoke-test inputs remain
available as `examples/products.dat` and `examples/users.dat`, using `-d 2`.

## Competitors: EAC and MILP

The two competitors are separate executables; they do not change the ARM build
or its anytime mode. EAC enumerates slab combinations and checks feasibility
with lp_solve. MILP solves a whole-problem mixed-integer model with CPLEX.
Both use Qhull for feasible-region vertices and the same input-file formats.
They are exact baselines, not anytime implementations.

Build from the repository root:

```bash
make -C competitors/eac -j2 PREFIX=/path/to/dependency-prefix \
  LP_SOLVE_LIB=/path/to/liblpsolve55.a
make -C competitors/milp -j2 PREFIX=/path/to/dependency-prefix \
  CPLEX_DIR=/path/to/CPLEX_Studio2211
```

EAC needs lp_solve and Qhull. MILP needs CPLEX/Concert and Qhull; it does not
link lp_solve. Neither competitor needs OR-Tools or GMP. The same shared-library
search-path setup described above applies when dependencies are not on the
system search path. No compiled libraries or proprietary dependencies are bundled.

Run the original small 2D example:

```bash
./competitors/eac/bin/eac -d 2 \
  -pf examples/products.dat -uf examples/users.dat -N 4 -wr 1.0 -ws 31
./competitors/milp/bin/milp -d 2 \
  -pf examples/products.dat -uf examples/users.dat -N 4 -wr 1.0 -ws 31
```

Both commands return `ARF: 2`. For the shared 4D IND example, use `-d 4`,
`-pf examples/ind_products.dat`, and `-uf examples/ind_users.dat`, keeping
`-N 4 -wr 1.0 -ws 31`; both return `ARF: 71`. EAC's enumeration can be much
slower than ARM or MILP, even on this small example.

For competitors, `-N` selects the number of generated constraint hyperplanes
(ARM uses `-wn`). `-wr` is the radius and `-ws` is the WCA seed.
Keep these parameters identical when comparing methods. Both print the selected
slabs, a solution point, `ARF:`, and total elapsed `TIME_SEC:` including input
loading and preprocessing. Neither command has an application timeout.

The packaged versions retain constant rank contributions from products pruned
above all feasible scores. Their reported ARF therefore includes those products,
not just the remaining slab indices. This reporting correction does not change
the minimizing slab search. They also reject missing or malformed inputs, print
solution coordinates at full double precision, and build without generated
benchmark or plotting files. The MILP solver requires an optimal status and sets
its relative and absolute MIP gap targets to zero (subject to solver tolerances).

## Third-party notices

See `THIRD_PARTY_NOTICES.md`. Preserve the original notices and license text in
the bundled third-party headers. External dependencies retain their own licenses.
