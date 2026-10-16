! Thin Fortran interface to the optional retained CUDA execution context.
! This module compiles without CUDA: the build capability remains available
! while C bindings are exposed only when their implementation is linked.
module fortsparse_native_cuda
    use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_size_t
    implicit none
    private
    public :: native_cuda_available
#ifdef FORTSPARSE_HAVE_NATIVE_CUDA
    logical, parameter :: native_cuda_available = .true.
    public :: spd_cuda_create, spd_cuda_destroy, spd_cuda_upload, spd_cuda_download
    public :: spd_cuda_synchronize, spd_cuda_panel_factor, spd_cuda_factor_info
    public :: spd_cuda_panel_solve, spd_cuda_schur_update, spd_cuda_schur_product
    public :: spd_cuda_allocated_bytes
    public :: spd_cuda_set_maps, spd_cuda_begin_factor, spd_cuda_scatter_update
    public :: spd_cuda_schur_scatter
    public :: spd_cuda_schur_task, spd_cuda_set_schur_batches, spd_cuda_schur_batch
    public :: spd_cuda_host_plan_bytes
    public :: spd_cuda_level_panel, spd_cuda_level_task, spd_cuda_set_levels
    public :: spd_cuda_factor_level, spd_cuda_update_level

    type, bind(C) :: spd_cuda_schur_task
        integer(c_size_t) :: left_offset, right_offset, destination_offset
        integer(c_size_t) :: row_map_offset, column_map_offset
        integer(c_int) :: ldl, ldr, m, n, k, destination_rows, destination_columns
    end type spd_cuda_schur_task
    type, bind(C) :: spd_cuda_level_panel
        integer(c_size_t) :: offset
        integer(c_int) :: width, rows, lda, level
    end type spd_cuda_level_panel
    type, bind(C) :: spd_cuda_level_task
        type(spd_cuda_schur_task) :: task
        integer(c_size_t) :: source, destination
    end type spd_cuda_level_task
    public :: spd_cuda_create_deferred, spd_cuda_prepare_blas, spd_cuda_front_factor

    interface
        integer(c_int) function spd_cuda_set_levels(context, panels, panel_count, &
            tasks, task_count, level_count) bind(C, name='fsparse_spd_cuda_set_levels')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context, panels, tasks
            integer(c_size_t), value :: panel_count, task_count, level_count
        end function spd_cuda_set_levels

        integer(c_int) function spd_cuda_factor_level(context, level) &
            bind(C, name='fsparse_spd_cuda_factor_level')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: level
        end function spd_cuda_factor_level

        integer(c_int) function spd_cuda_update_level(context, level) &
            bind(C, name='fsparse_spd_cuda_update_level')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: level
        end function spd_cuda_update_level
        integer(c_int) function spd_cuda_set_schur_batches(context, tasks, task_count, &
            batch_ptr, batch_count) bind(C, name='fsparse_spd_cuda_set_schur_batches')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context, tasks, batch_ptr
            integer(c_size_t), value :: task_count, batch_count
        end function spd_cuda_set_schur_batches

        integer(c_int) function spd_cuda_schur_batch(context, batch) &
            bind(C, name='fsparse_spd_cuda_schur_batch')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: batch
        end function spd_cuda_schur_batch

        integer(c_size_t) function spd_cuda_host_plan_bytes(context) &
            bind(C, name='fsparse_spd_cuda_host_plan_bytes')
            import :: c_ptr, c_size_t
            type(c_ptr), value :: context
        end function spd_cuda_host_plan_bytes

        integer(c_int) function spd_cuda_create(context, device, capacity) &
            bind(C, name='fsparse_spd_cuda_create')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), intent(out) :: context
            integer(c_int), value :: device
            integer(c_size_t), value :: capacity
        end function spd_cuda_create

        integer(c_int) function spd_cuda_create_deferred(context, device, capacity) &
            bind(C, name='fsparse_spd_cuda_create_deferred')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), intent(out) :: context
            integer(c_int), value :: device
            integer(c_size_t), value :: capacity
        end function spd_cuda_create_deferred

        integer(c_int) function spd_cuda_prepare_blas(context) &
            bind(C, name='fsparse_spd_cuda_prepare_blas')
            import :: c_ptr, c_int
            type(c_ptr), value :: context
        end function spd_cuda_prepare_blas

        integer(c_int) function spd_cuda_front_factor(context, offset, width, rows, lda) &
            bind(C, name='fsparse_spd_cuda_front_factor')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: offset
            integer(c_int), value :: width, rows, lda
        end function spd_cuda_front_factor

        integer(c_int) function spd_cuda_destroy(context) &
            bind(C, name='fsparse_spd_cuda_destroy')
            import :: c_ptr, c_int
            type(c_ptr), value :: context
        end function spd_cuda_destroy

        integer(c_int) function spd_cuda_upload(context, offset, source, count) &
            bind(C, name='fsparse_spd_cuda_upload')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context, source
            integer(c_size_t), value :: offset, count
        end function spd_cuda_upload

        integer(c_int) function spd_cuda_download(context, offset, destination, count) &
            bind(C, name='fsparse_spd_cuda_download')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context, destination
            integer(c_size_t), value :: offset, count
        end function spd_cuda_download

        integer(c_int) function spd_cuda_synchronize(context) &
            bind(C, name='fsparse_spd_cuda_synchronize')
            import :: c_ptr, c_int
            type(c_ptr), value :: context
        end function spd_cuda_synchronize

        integer(c_int) function spd_cuda_panel_factor(context, offset, n, lda) &
            bind(C, name='fsparse_spd_cuda_panel_factor')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: offset
            integer(c_int), value :: n, lda
        end function spd_cuda_panel_factor

        integer(c_int) function spd_cuda_factor_info(context, info) &
            bind(C, name='fsparse_spd_cuda_factor_info')
            import :: c_ptr, c_int
            type(c_ptr), value :: context
            integer(c_int), intent(out) :: info
        end function spd_cuda_factor_info

        integer(c_int) function spd_cuda_panel_solve(context, l_offset, n, ldl, &
                                                   b_offset, m, ldb) &
            bind(C, name='fsparse_spd_cuda_panel_solve')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: l_offset, b_offset
            integer(c_int), value :: n, ldl, m, ldb
        end function spd_cuda_panel_solve

        integer(c_int) function spd_cuda_schur_update(context, c_offset, ldc, &
                                                    left_offset, ldl, &
                                                    right_offset, ldr, m, n, k) &
            bind(C, name='fsparse_spd_cuda_schur_update')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: c_offset, left_offset, right_offset
            integer(c_int), value :: ldc, ldl, ldr, m, n, k
        end function spd_cuda_schur_update

        integer(c_int) function spd_cuda_schur_product(context, c_offset, ldc, &
                                                     left_offset, ldl, &
                                                     right_offset, ldr, m, n, k) &
            bind(C, name='fsparse_spd_cuda_schur_product')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: c_offset, left_offset, right_offset
            integer(c_int), value :: ldc, ldl, ldr, m, n, k
        end function spd_cuda_schur_product

        integer(c_size_t) function spd_cuda_allocated_bytes(context) &
            bind(C, name='fsparse_spd_cuda_allocated_bytes')
            import :: c_ptr, c_size_t
            type(c_ptr), value :: context
        end function spd_cuda_allocated_bytes

        integer(c_int) function spd_cuda_set_maps(context, maps, count) &
            bind(C, name='fsparse_spd_cuda_set_maps')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context, maps
            integer(c_size_t), value :: count
        end function spd_cuda_set_maps

        integer(c_int) function spd_cuda_begin_factor(context) &
            bind(C, name='fsparse_spd_cuda_begin_factor')
            import :: c_ptr, c_int
            type(c_ptr), value :: context
        end function spd_cuda_begin_factor

        integer(c_int) function spd_cuda_scatter_update(context, update_offset, m, n, &
            destination_offset, destination_rows, destination_columns, row_map, column_map) &
            bind(C, name='fsparse_spd_cuda_scatter_update')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: update_offset, destination_offset, row_map, column_map
            integer(c_int), value :: m, n, destination_rows, destination_columns
        end function spd_cuda_scatter_update

        integer(c_int) function spd_cuda_schur_scatter(context, left_offset, ldl, &
            right_offset, ldr, m, n, k, destination_offset, destination_rows, &
            destination_columns, row_map, column_map) &
            bind(C, name='fsparse_spd_cuda_schur_scatter')
            import :: c_ptr, c_int, c_size_t
            type(c_ptr), value :: context
            integer(c_size_t), value :: left_offset, right_offset, destination_offset
            integer(c_size_t), value :: row_map, column_map
            integer(c_int), value :: ldl, ldr, m, n, k
            integer(c_int), value :: destination_rows, destination_columns
        end function spd_cuda_schur_scatter
    end interface
#else
    logical, parameter :: native_cuda_available = .false.
#endif
end module fortsparse_native_cuda
