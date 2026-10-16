// MIT adapter. FreeFem SDK headers retain their own LGPL terms; see README.md.
#include "ff++.hpp"
#include "fortsparse_native.h"
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>

namespace {
static_assert(std::is_same<int,int32_t>::value, "FreeFem I4 indices required");
static_assert(std::is_same<Complex,fs_native_complex>::value,
              "FreeFem/C ABI complex types must match");

template<class R> struct NativeABI;
template<> struct NativeABI<double> {
  static constexpr int scalar_kind = 0;
  static int32_t factor(void *h, int32_t n, int32_t nz, const int32_t *p,
                        const int32_t *r, const double *v, bool reuse) {
    return reuse ? fs_native_refactor_d(h,n,nz,p,r,v) : fs_native_factor_d(h,n,nz,p,r,v);
  }
  static int32_t solve(void *h, int32_t trans, int32_t n, int32_t nrhs,
                       const double *b, double *x) {
    return fs_native_solve_d(h,trans,n,nrhs,b,x);
  }
};
template<> struct NativeABI<Complex> {
  static constexpr int scalar_kind = 1;
  static int32_t factor(void *h, int32_t n, int32_t nz, const int32_t *p,
                        const int32_t *r, const Complex *v, bool reuse) {
    return reuse ? fs_native_refactor_z(h,n,nz,p,r,v) : fs_native_factor_z(h,n,nz,p,r,v);
  }
  static int32_t solve(void *h, int32_t trans, int32_t n, int32_t nrhs,
                       const Complex *b, Complex *x) {
    // FreeFem's complex trans=1 is conjugate transpose (UMFPACK_At),
    // not plain transpose (UMFPACK_Aat). No complex scalar crosses by value.
    return fs_native_solve_z(h,trans ? 2 : 0,n,nrhs,b,x);
  }
};

template<class R> class NativeSolver : public VirtualSolver<int,R> {
public:
  static const int orTypeSol = 1|16|32; // SDK general serial real/complex solver flags.
  NativeSolver(HashMatrix<int,R> &input, const Data_Sparse_Solver &, Stack)
    : matrix_(&input), handle_(nullptr), pointer_(nullptr), row_(nullptr),
      values_(nullptr), factored_(false) {
    if (mpisize > 1) ExecError("FortSparse native FreeFem solver requires one process");
    check(fs_native_create(NativeABI<R>::scalar_kind, &handle_));
    try {
      // Explicit ordering is required by both concrete Maxwell consumers.
      check(fs_native_set_ordering(handle_, FS_NATIVE_ORDER_AMD));
    } catch (...) {
      fs_native_free(&handle_);
      throw;
    }
  }
  ~NativeSolver() override { fs_native_free(&handle_); }
  void UpdateState() override {
    const bool numeric = matrix_->GetReDoNumerics();
    const bool symbolic = matrix_->GetReDoSymbolic();
    if (symbolic) this->state = std::min(this->state, 1);
    else if (numeric) this->state = std::min(this->state, 2);
  }
  void fac_symbolic() override {
    if (matrix_->n <= 0 || matrix_->n != matrix_->m ||
        matrix_->n >= std::numeric_limits<int32_t>::max() ||
        matrix_->nnz >= static_cast<size_t>(std::numeric_limits<int32_t>::max()))
      ExecError("FortSparse FreeFem requires square 32-bit sparse input");
    // Preserve SDK real/complex symmetric/Hermitian half-storage semantics.
    matrix_->UnHalf();
    if (matrix_->nnz >= static_cast<size_t>(std::numeric_limits<int32_t>::max()))
      ExecError("FortSparse FreeFem expanded matrix exceeds 32-bit nnz");
    matrix_->setfortran(0);
    matrix_->CSC(pointer_, row_, values_);
  }
  void fac_numeric() override {
    // SDK updates can relocate the arrays. The C ABI copies them during this call.
    matrix_->CSC(pointer_, row_, values_);
    check(NativeABI<R>::factor(handle_, matrix_->n, static_cast<int32_t>(matrix_->nnz),
                              pointer_, row_, values_, factored_));
    factored_ = true;
    // UnHalf/CSC may themselves set SDK change flags. This factor includes
    // those transformations; do not refactor an unchanged second RHS.
    matrix_->GetReDoNumerics();
    matrix_->GetReDoSymbolic();
  }
  void dosolver(R *x, R *b, int nrhs, int transpose) override {
    if (nrhs == 0) return;
    if (nrhs < 0 || (transpose != 0 && transpose != 1))
      ExecError("Invalid FortSparse FreeFem solve operation");
    check(NativeABI<R>::solve(handle_, transpose, matrix_->n, nrhs, b, x));
  }
private:
  NativeSolver(const NativeSolver &) = delete;
  NativeSolver &operator=(const NativeSolver &) = delete;
  HashMatrix<int,R> *matrix_;
  void *handle_;
  int *pointer_, *row_;
  R *values_;
  bool factored_;
  void check(int32_t status) const {
    if (status == 0) return;
    char detail[160] = {};
    fs_native_last_error(handle_, detail, sizeof detail);
    const std::string message = "FortSparse native FreeFem failure " +
      std::to_string(status) + ": " + detail;
    ExecError(message.c_str());
  }
};

void initialize() {
  if (mpisize > 1) ExecError("FortSparse native FreeFem plugin requires one process");
  addsolver<NativeSolver<double>>("FORTSPARSE", 1000, 3);
  addsolver<NativeSolver<Complex>>("FORTSPARSE", 1000, 3);
  if (verbosity) cout << "FreeFem sparse solver: FortSparse native real/complex (no fallback)" << endl;
}
}
LOADFUNC(initialize)
