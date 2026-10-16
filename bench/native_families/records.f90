module native_family_records
    use native_family_fixtures, only: dp
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private
    public :: family_record_t, process_metrics, parent_cpu_seconds
    interface
        function parent_cpu_seconds() result(seconds) &
            bind(C, name='fsbench_parent_cpu_seconds')
            import :: c_double
            real(c_double) :: seconds
        end function parent_cpu_seconds
        subroutine process_metrics(rss, child_rss, child_cpu) &
            bind(C, name='fsbench_process_metrics')
            import :: c_double
            real(c_double), intent(out) :: rss, child_rss, child_cpu
        end subroutine
    end interface
    type :: family_record_t
        character(32) :: backend = '', family = '', workload = ''
        character(32) :: variant = '', operation = 'N'
        character(32) :: native_factor_family = ''
        integer :: level = 0, n = 0, nnz = 0, threads = 1, nrhs = 1, status_code = 0
        integer :: memory_policy = -1
        integer(int64) :: host_memory_limit_bytes = -1
        integer :: active_threads = -1, native_factor_serial = 0
        logical :: reused = .false.
        real(dp) :: total_s = -1, parent_cpu_s = -1, helper_cpu_s = -1
        real(dp) :: parent_main_cpu_s = -1
        real(dp) :: parent_rss_kib = -1, helper_rss_kib = -1
        real(dp) :: factor_s = -1, solve_s = -1, refactor_s = -1, prime_factor_s = -1
    real(dp) :: analysis_s = -1, numeric_s = -1, backward_error = -1, forward_error = -1
        real(dp) :: factor_mem_mb = -1
        integer(int64) :: factor_nnz = -1, host_current_bytes = -1, host_peak_bytes = -1
    contains
        procedure :: write_csv
    end type
contains
    subroutine write_csv(self)
        class(family_record_t), intent(in) :: self
        character(24) :: outcome
        outcome = 'failure'
        if (self%status_code == 0) outcome = 'ok'
        if (self%status_code == 3) outcome = 'backend_unavailable'
        write (*, '(a)') 'backend,family,level,threads,workload,variant,operation,'// &
            'n,nnz,nrhs,status,status_code,factor_s,solve_s,refactor_s,reused,'// &
            'total_s,parent_cpu_s,parent_main_cpu_s,helper_cpu_s,'// &
            'parent_rss_kib,helper_rss_kib,'// &
            'analysis_s,numeric_s,prime_factor_s,factor_nnz,factor_mem_mb,'// &
            'host_current_bytes,host_peak_bytes,active_threads,'// &
            'native_factor_serial,native_factor_family,backward_error,forward_error,'// &
            'memory_policy,host_memory_limit_bytes'
        write (*, '(*(g0,:,","))') trim(self%backend), trim(self%family), &
            self%level, self%threads, trim(self%workload), trim(self%variant), &
            trim(self%operation), self%n, self%nnz, self%nrhs, trim(outcome), &
          self%status_code, self%factor_s, self%solve_s, self%refactor_s, self%reused, &
            self%total_s, self%parent_cpu_s, self%parent_main_cpu_s, &
            self%helper_cpu_s, self%parent_rss_kib, &
            self%helper_rss_kib, self%analysis_s, self%numeric_s, self%prime_factor_s, &
            self%factor_nnz, self%factor_mem_mb, self%host_current_bytes, &
            self%host_peak_bytes, self%active_threads, self%native_factor_serial, &
            trim(self%native_factor_family), self%backward_error, self%forward_error, &
            self%memory_policy, self%host_memory_limit_bytes
    end subroutine
end module native_family_records
