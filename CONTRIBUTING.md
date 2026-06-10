# Contributing to DistributionTool

Thanks for your interest in improving DistributionTool. This is a research tool;
contributions that improve correctness, performance, documentation or test
coverage are all welcome.

## Development setup

```bash
git clone https://github.com/TiagoBe0/DistributionTool.git
cd DistributionTool
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
```

Dependencies: a C++17 compiler (GCC ≥ 7 / Clang ≥ 5), CMake ≥ 3.14, Eigen3 and
(optionally) OpenMP. See the [README](README.md#build) for the full list.

## Guidelines

- **Build clean.** Code must compile with `-Wall -Wextra` and no new warnings.
- **Add tests.** New numerical kernels or parsers should come with a test in
  [`tests/`](tests/) using the lightweight `test_harness.h` macros. Run the
  suite before opening a PR.
- **Keep the hot paths fast.** The SOAP / neighbour-list code is performance
  critical (see [docs/architecture.md](docs/architecture.md)); avoid per-query
  heap allocations and prefer the existing `CellList` index.
- **Match the surrounding style.** Header-only utilities live in
  [`include/`](include/), implementations in [`src/`](src/).
- **Document physics changes.** If you change a method or default, update the
  relevant section of [docs/methodology.md](docs/methodology.md).

## Pull requests

1. Branch off `main`.
2. Keep commits focused and write descriptive messages.
3. Ensure CI (build + `ctest`) is green.
4. Describe *what* changed and *why* in the PR body; link any related issue.

## Reporting issues

Open a GitHub issue with: the command you ran, the input (or a minimal LAMMPS
dump that reproduces it), the expected vs. actual output, and your OS / compiler
versions.
