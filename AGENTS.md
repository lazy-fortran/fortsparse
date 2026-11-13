# Contributor guidance

Read [CONTRIBUTING.md](CONTRIBUTING.md), [the user guide](docs/user-guide.md)
and [implementation status](docs/status.md) before changing the library.

- Keep the backend abstraction, matrix-class contracts and status codes.
- Generate project-owned pointwise mathematics through `gen/`; keep CPU and
  CUDA output reproducible. Preserve independent residual and derivative tests.
- Keep the permissive library separate from the optional GPL helper process.
- Use native CMake/CTest checks and affected fpm tests; do not weaken gates.
- Preserve unrelated edits and state unsupported paths and verification limits.
