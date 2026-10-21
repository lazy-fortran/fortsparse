#!/bin/bash
# Execute the repository's original native CI checks in isolated Slurm scratch.
set -euo pipefail
: "${SLURM_JOB_ID:?ci/run.sh must run in Slurm}"
tier=$(cat "$(dirname "$(dirname "$(realpath "$0")")")/../tier.txt")
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$(dirname "$src")
cd "$src"
. "$root/task.env"
. /home/ert/opt/kin6d-env.sh
export PATH=/home/ert/.local/bin:$PATH
export GIT_CONFIG_GLOBAL=$root/gitconfig
: > "$GIT_CONFIG_GLOBAL"
git config --global protocol.file.allow always
shopt -s nullglob
for mirror in "$root"/mirrors/*/*.git; do
    slug=${mirror#"$root/mirrors/"}; slug=${slug%.git}
    git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug.git"
    git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug"
    git config --global --add "url.file://$mirror.insteadOf" "git@github.com:$slug.git"
done
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export CMAKE_BUILD_PARALLEL_LEVEL=${SCI_CPUS:-4} CTEST_PARALLEL_LEVEL=${SCI_CPUS:-4}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
mkdir -p ci-out
{ hostname; /usr/bin/gfortran --version | head -1; cat "$root/revision.txt"; } | tee ci-out/toolchain.txt
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP"
load_job_environment() {
    while IFS='=' read -r key value; do
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] && export "$key=$value"
    done < "$GITHUB_ENV"
    while IFS= read -r item; do [ -z "$item" ] || export PATH="$item:$PATH"; done < "$GITHUB_PATH"
}
case "$tier" in
  ci-native-analysis)
    export OMP_NUM_THREADS="8"
    export OPENBLAS_NUM_THREADS="1"
    export OMP_DYNAMIC="FALSE"
    export OMP_WAIT_POLICY="passive"
    echo 'Configure native Release'
    (
cmake -S . -B native -G Ninja -DCMAKE_BUILD_TYPE=Release -DFORTSPARSE_BUILD_TESTING=ON -DFORTSPARSE_ENABLE_SUPERLU=OFF -DFORTSPARSE_ENABLE_UMFPACK_IPC=OFF -DFORTSPARSE_ENABLE_MUMPS=OFF

    )
    load_job_environment
    echo 'Build affected native gates'
    (
cmake --build native -j2 --target test_fortsparse_slu_factor test_fortsparse_slu_z_api test_fortsparse_slu_threads

    )
    load_job_environment
    echo 'Test values and parallel reuse'
    (
ctest --test-dir native -j1 --output-on-failure --no-tests=error -R '^test_fortsparse_(slu_factor|slu_z_api|slu_threads_uncapped|slu_threads_repeat)$'

    )
    load_job_environment
    ;;
  ci-cmake)
    export OMP_NUM_THREADS="1"
    export OPENBLAS_NUM_THREADS="1"
    echo Configure
    (
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DFORTSPARSE_BUILD_TESTING=ON \
  -DFORTSPARSE_USE_SYSTEM_DEPS=ON


    )
    load_job_environment
    echo 'Build focused solver gates'
    (
cmake --build build -j 2 --target \
  test_fortsparse_status \
  test_fortsparse_solve_real \
  test_fortsparse_solve_complex \
  test_fortsparse_solve_derivatives \
  test_fortsparse_complex_solve_derivatives \
  test_fortsparse_superlu_symmetric \
  test_fortsparse_errors \
  test_fortsparse_superlu_status \
  test_fortsparse_umfpack_ipc \
  test_fortsparse_umfpack_controls \
  fortsparse_umfpack_helper \
  test_fortsparse_precision \
  test_fortsparse_mixed_lu \
  test_fortsparse_stability \
  test_mixed_allocations \
  test_fortsparse_native_general_api \
  test_fortsparse_slu_z_api \
  test_fortsparse_ldlt_blocked \
  test_fortsparse_native_ldlt \
  test_fortsparse_native_ldlt_z \
  test_fortsparse_inertia


    )
    load_job_environment
    echo 'Test focused solver gates'
    (
ctest --test-dir build -j 2 --output-on-failure --no-tests=error -R '^test_fortsparse_(status|solve_real|solve_complex|solve_derivatives|complex_solve_derivatives|superlu_symmetric|errors|superlu_status|umfpack_ipc|umfpack_controls|precision|mixed_lu|stability|native_general_api|slu_z_api|ldlt_blocked|native_ldlt|native_ldlt_z|inertia)$|^test_mixed_allocations$'

    )
    load_job_environment
    ;;
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
