submodule (fortsparse_native_slu) fortsparse_native_slu_z_factor
    ! Complex numeric factorization of the static-structure supernodal LU.
    ! The same left-looking scheme, OpenMP subtree tasks, block-local
    ! pivoting with counted perturbation and statistics as the real factor;
    ! the arithmetic is the generated complex leaves (fortsparse_slu_z_leaves)
    ! and ZGEMM/ZTRSM. U12 is stored as a plain transpose (Ut), never
    ! conjugated: the adjoint solve conjugates the right-hand side instead.
    ! Nothing is allocated (the OpenMP runtime owns its task descriptors).
    use fortsparse_slu_z_leaves, only: slu_z_update, slu_z_scale, &
        slu_z_block_factor, slu_z_trsm_upper, slu_z_trsm_unit_lower_t, &
        slu_z_schur_update, slu_z_scatter_mapped_update
    implicit none

    interface
        subroutine zgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import dp
            character, intent(in) :: ta, tb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(dp), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            complex(dp), intent(inout) :: c(ldc,*)
        end subroutine zgemm
        subroutine ztrsm(side, uplo, trans, diag, m, n, alpha, a, lda, b, ldb)
            import dp
            character, intent(in) :: side, uplo, trans, diag
            integer, intent(in) :: m, n, lda, ldb
            complex(dp), intent(in) :: alpha, a(lda,*)
            complex(dp), intent(inout) :: b(ldb,*)
        end subroutine ztrsm
    end interface

    complex(dp), parameter :: zone = (1.0_dp, 0.0_dp), zmone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zzero = (0.0_dp, 0.0_dp)
contains
#ifdef FORTSPARSE_HAVE_OMP_DEPEND_ITERATOR
    module subroutine factor_tasks_z(self)
        ! One task per elimination-subtree group, depending on its
        ! predecessor groups (as the real factor).
        class(slu_factor_t), target, intent(inout) :: self
        integer, pointer :: token(:), sources(:)
        integer :: s, incoming_start, incoming_end
        token => self%sched%group_tokens
        sources => self%sched%group_sources
        !$omp parallel num_threads(self%sched%active_threads) default(shared)
        !$omp single
        !$ self%sched%active_threads = omp_get_num_threads()
        do s = 1, self%sched%task_groups
            incoming_start = self%sched%group_ptr(s)
            ! OpenMP iterator end is exclusive, including in Fortran.
            incoming_end = self%sched%group_ptr(s+1)
            !$omp task firstprivate(s) &
            !$omp depend(iterator(incoming=incoming_start:incoming_end), &
            !$omp in:token(sources(incoming))) depend(out:token(s))
            call factor_group_z(self, s)
            !$omp end task
        end do
        !$omp end single
        !$omp end parallel
    end subroutine factor_tasks_z

    subroutine factor_group_z(self, group)
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: group
        integer :: member, thread
        thread = self%sched%workspace_thread()
        do member = self%sched%group_member_ptr(group), &
                self%sched%group_member_ptr(group+1)-1
            call factor_destination_z(self, self%sched%group_members(member), thread)
        end do
    end subroutine factor_group_z
#else
    module subroutine factor_tasks_z(self)
        ! Without iterator dependencies: one task per panel and per update,
        ! ordered by the per-panel failure flags as dependence tokens.
        class(slu_factor_t), target, intent(inout) :: self
        integer, pointer :: token(:)
        integer :: s, group, first, q, dest
        token => self%failures
        !$omp parallel num_threads(self%sched%active_threads) default(shared)
        !$omp single
        !$ self%sched%active_threads = omp_get_num_threads()
        do s = 1, self%node_count
            !$omp task firstprivate(s) depend(inout:token(s))
            call factor_panel_z(self, s, self%sched%workspace_thread())
            !$omp end task
            first = self%panels(s)%width+1
            do while (first <= self%panels(s)%m)
                call self%sched%update_group(self%panels(s)%supernode_t, first, group, &
                    dest)
                q = group-first
                !$omp task firstprivate(s,dest,first,q) &
                !$omp depend(in:token(s)) depend(inout:token(dest))
                call update_panel_z(self, s, dest, first, q)
                !$omp end task
                first = group
            end do
        end do
        !$omp end single
        !$omp end parallel
    end subroutine factor_tasks_z

    subroutine update_panel_z(self, source, dest, first, q)
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: source, dest, first, q
        integer :: thread
        thread = self%sched%workspace_thread()
        call set_rowmap(self, dest, thread, .true.)
        call accumulate_update_z(self, source, dest, first, q, thread)
        call set_rowmap(self, dest, thread, .false.)
    end subroutine update_panel_z
#endif

    subroutine set_rowmap(self, dest, thread, on)
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: dest, thread
        logical, intent(in) :: on
        integer :: p
        if (on) then
            do p = 1, self%panels(dest)%m
                self%sched%rowmap(self%panels(dest)%rows(p),thread) = p
            end do
        else
            do p = 1, self%panels(dest)%m
                self%sched%rowmap(self%panels(dest)%rows(p),thread) = 0
            end do
        end if
    end subroutine set_rowmap

    module subroutine factor_destination_z(self, dest, thread)
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: dest, thread
        integer :: p
        call set_rowmap(self, dest, thread, .true.)
        do p = self%sched%incoming_ptr(dest), self%sched%incoming_ptr(dest+1)-1
            call accumulate_update_z(self, self%sched%incoming_source(p), dest, &
                self%sched%incoming_first(p), self%sched%incoming_q(p), thread)
        end do
        call factor_panel_z(self, dest, thread)
        call set_rowmap(self, dest, thread, .false.)
    end subroutine factor_destination_z

    module subroutine assemble_z(self, ptr, rows, values, status)
        ! M(k,l) = drp(k)*A(i,c)*dcp(l) into the owning L or Ut panel;
        ! ||M||_inf (moduli, row sums in the real part of the solve scratch)
        ! and max |M|.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: ptr(:), rows(:)
        complex(dp), intent(in) :: values(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: l, c, e, k
        integer(i8) :: p
        complex(dp) :: v
        self%zvalues = zzero
        self%zsolve_work(1:self%n) = zzero
        self%max_abs_m = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'supernodal LU: pattern outside the analyzed structure')
        do l = 1, self%n
            c = self%col_old(l)
            do e = ptr(c), ptr(c+1)-1
                k = self%row_new(rows(e))
                v = slu_z_scale(self%drp(k), values(e), self%dcp(l))
                self%zsolve_work(k) = slu_z_update(self%zsolve_work(k), zmone, &
                    cmplx(abs(v), 0.0_dp, dp))
                self%max_abs_m = max(self%max_abs_m, abs(v))
                p = entry_offset(self, k, l)
                if (p == 0_i8) then
                    self%zsolve_work(1:self%n) = zzero
                    return
                end if
                self%zvalues(p) = v
            end do
        end do
        self%anorm = 0.0_dp
        do k = 1, self%n
            self%anorm = max(self%anorm, real(self%zsolve_work(k),dp))
        end do
        self%tau_abs = self%tau*self%anorm
        self%zsolve_work(1:self%n) = zzero
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine assemble_z

    integer(i8) function entry_offset(self, k, l) result(p)
        ! Panel offset of M(k,l) in the static structure, or 0 outside it.
        class(slu_factor_t), intent(in) :: self
        integer, intent(in) :: k, l
        integer :: dk, dl, pos
        p = 0_i8
        dl = self%sched%colnode(l)
        dk = self%sched%colnode(k)
        if (dk == dl) then
            p = self%panels(dl)%offset_l &
                +int(l-self%panels(dl)%first,i8)*self%panels(dl)%m &
                +(k-self%panels(dl)%first+1)
        else if (k > l) then
            pos = supernode_find_row(self%panels(dl)%supernode_t, k)
            if (pos == 0) return
            p = self%panels(dl)%offset_l &
                +int(l-self%panels(dl)%first,i8)*self%panels(dl)%m+pos
        else
            pos = supernode_find_row(self%panels(dk)%supernode_t, l)
            if (pos == 0) return
            p = self%panels(dk)%offset_u+int(k-self%panels(dk)%first,i8) &
                *(self%panels(dk)%m-self%panels(dk)%width) &
                +(pos-self%panels(dk)%width)
        end if
    end function entry_offset

    subroutine accumulate_update_z(self, s, d, f, q, thread)
        ! GEMM1 (destination L panel, full diagonal block included) and
        ! GEMM2 (destination Ut panel) of the real accumulate_update.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: s, d, f, q, thread
        integer :: m, w, p, position
        logical :: contiguous_rows
        if (self%failures(s) /= 0) self%failures(d) = self%failures(s)
        if (self%failures(d) /= 0) return
        m = self%panels(s)%m
        w = self%panels(s)%width
        if (w == 1) then
            call rank_one_update_z(self, s, d, f, q, thread)
            return
        end if
        contiguous_rows = .false.
        position = self%sched%rowmap(self%panels(s)%rows(f),thread)
        if (self%use_blas) then
            contiguous_rows = .true.
            do p = f+1, m
                if (self%sched%rowmap(self%panels(s)%rows(p),thread) == &
                    position+p-f) cycle
                contiguous_rows = .false.
                exit
            end do
        end if
        if (contiguous_rows) then
            call direct_update_z(self, s, d, f, q, position)
        else
            call mapped_update_z(self, s, d, f, q, thread)
        end if
    end subroutine accumulate_update_z

    subroutine rank_one_update_z(self, s, d, f, q, thread)
        ! Width-one source: scatter each column directly through the map.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: s, d, f, q, thread
        integer :: m, w, m2, md, mud, wd, j, col
        integer(i8) :: ls, us, ld, ud
        m = self%panels(s)%m
        w = self%panels(s)%width
        m2 = m-f+1
        md = self%panels(d)%m
        wd = self%panels(d)%width
        mud = md-wd
        ls = self%panels(s)%offset_l
        us = self%panels(s)%offset_u
        ld = self%panels(d)%offset_l
        ud = self%panels(d)%offset_u
        do j = 1, q
            col = self%panels(s)%rows(f+j-1)-self%panels(d)%first+1
            call slu_z_scatter_mapped_update(self%zvalues(ld+int(col-1,i8)*md+1), &
                self%panels(s)%rows(f:), self%sched%rowmap(:,thread), &
                self%zvalues(ls+f), m2, self%zvalues(us+f-w+j-1))
            if (m2 > q) call slu_z_scatter_mapped_update( &
                self%zvalues(ud+int(col-1,i8)*mud+1-wd), self%panels(s)%rows(f+q:), &
                self%sched%rowmap(:,thread), self%zvalues(us+f+q-w), m2-q, &
                self%zvalues(ls+f+j-1))
        end do
    end subroutine rank_one_update_z

    subroutine direct_update_z(self, s, d, f, q, position)
        ! Destination rows and columns are both consecutive: ZGEMM in place.
        ! Ut_d is addressed from w_d entries before its column start, so the
        ! row map positions (> w_d) index it directly.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: s, d, f, q, position
        integer :: m, w, mu, m2, md, mud, wd, col
        integer(i8) :: ls, us, ld, ud
        m = self%panels(s)%m
        w = self%panels(s)%width
        mu = m-w
        m2 = m-f+1
        md = self%panels(d)%m
        wd = self%panels(d)%width
        mud = md-wd
        ls = self%panels(s)%offset_l
        us = self%panels(s)%offset_u
        ld = self%panels(d)%offset_l
        ud = self%panels(d)%offset_u
        col = self%panels(s)%rows(f)-self%panels(d)%first+1
        call zgemm('N', 'T', m2, q, w, zmone, self%zvalues(ls+f), m, &
            self%zvalues(us+f-w), mu, zone, &
            self%zvalues(ld+int(col-1,i8)*md+position), md)
        if (m2 > q) call zgemm('N', 'T', m2-q, q, w, zmone, &
            self%zvalues(us+f+q-w), mu, self%zvalues(ls+f), m, zone, &
            self%zvalues(ud+int(col-1,i8)*mud+position+q-wd), mud)
    end subroutine direct_update_z

    subroutine mapped_update_z(self, s, d, f, q, thread)
        ! Products into the thread's update scratch (ZGEMM or the generated
        ! Schur leaf), then the generated mapped scatter.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: s, d, f, q, thread
        integer :: m, w, mu, m2, md, mud, wd, j, col
        integer(i8) :: ls, us, ld, ud
        m = self%panels(s)%m
        w = self%panels(s)%width
        mu = m-w
        m2 = m-f+1
        md = self%panels(d)%m
        wd = self%panels(d)%width
        mud = md-wd
        ls = self%panels(s)%offset_l
        us = self%panels(s)%offset_u
        ld = self%panels(d)%offset_l
        ud = self%panels(d)%offset_u
        call update_product_z(self, thread, m2, q, w, ls+f, m, us+f-w, mu)
        ! zupdate_work holds -L_s Ut_s^T; the mapped scatter subtracts it times -1.
        do j = 1, q
            col = self%panels(s)%rows(f+j-1)-self%panels(d)%first+1
            call slu_z_scatter_mapped_update(self%zvalues(ld+int(col-1,i8)*md+1), &
                self%panels(s)%rows(f:), self%sched%rowmap(:,thread), &
                self%zupdate_work((j-1)*m2+1,thread), m2, zmone)
        end do
        if (m2 == q) return
        call update_product_z(self, thread, m2-q, q, w, us+f+q-w, mu, ls+f, m)
        do j = 1, q
            col = self%panels(s)%rows(f+j-1)-self%panels(d)%first+1
            call slu_z_scatter_mapped_update( &
                self%zvalues(ud+int(col-1,i8)*mud+1-wd), self%panels(s)%rows(f+q:), &
                self%sched%rowmap(:,thread), &
                self%zupdate_work((j-1)*(m2-q)+1,thread), m2-q, zmone)
        end do
    end subroutine mapped_update_z

    subroutine update_product_z(self, thread, rows, cols, k, left, ldl, right, ldr)
        ! zupdate_work(1:rows,1:cols) := -left(1:rows,1:k) right(1:cols,1:k)^T.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: thread, rows, cols, k, ldl, ldr
        integer(i8), intent(in) :: left, right
        if (self%use_blas) then
            call zgemm('N', 'T', rows, cols, k, zmone, self%zvalues(left), ldl, &
                self%zvalues(right), ldr, zzero, self%zupdate_work(1,thread), rows)
        else
            self%zupdate_work(1:rows*cols,thread) = zzero
            call slu_z_schur_update(self%zupdate_work(1,thread), rows, &
                self%zvalues(left), ldl, self%zvalues(right), ldr, rows, cols, k)
        end if
    end subroutine update_product_z

    subroutine factor_panel_z(self, d, thread)
        ! Restricted-pivot block LU, exchanges applied to the Ut columns,
        ! then L21 := L21 U^-1 and Ut := Ut L^-T (plain transposes).
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: d, thread
        integer :: m, w, mu, first, perturbed, swaps
        integer(i8) :: l0, u0, e
        real(dp) :: min_pivot, max_u
        if (self%failures(d) /= 0) return
        m = self%panels(d)%m
        w = self%panels(d)%width
        mu = m-w
        first = self%panels(d)%first
        l0 = self%panels(d)%offset_l
        u0 = self%panels(d)%offset_u
        perturbed = 0
        swaps = 0
        min_pivot = huge(1.0_dp)
        max_u = 0.0_dp
        call slu_z_block_factor(self%zvalues(l0+1), m, w, self%tau_abs, &
            self%pivots(first), perturbed, swaps, min_pivot, max_u)
        if (mu > 0) call panel_solves_z(self, d, max_u)
        self%thread_stats(1,thread) = self%thread_stats(1,thread)+perturbed
        self%thread_stats(2,thread) = self%thread_stats(2,thread)+swaps
        self%thread_stats(3,thread) = min(self%thread_stats(3,thread), min_pivot)
        self%thread_stats(4,thread) = max(self%thread_stats(4,thread), max_u)
        do e = l0+1, u0+int(mu,i8)*w
            if (ieee_is_finite(real(self%zvalues(e),dp))) then
                if (ieee_is_finite(aimag(self%zvalues(e)))) cycle
            end if
            self%failures(d) = 1
            return
        end do
    end subroutine factor_panel_z

    subroutine panel_solves_z(self, d, max_u)
        ! Block exchanges as Ut column swaps, then both triangular solves.
        class(slu_factor_t), intent(inout) :: self
        integer, intent(in) :: d
        real(dp), intent(inout) :: max_u
        integer :: m, w, mu, k, r, i, first
        integer(i8) :: l0, u0, e
        complex(dp) :: saved
        m = self%panels(d)%m
        w = self%panels(d)%width
        mu = m-w
        first = self%panels(d)%first
        l0 = self%panels(d)%offset_l
        u0 = self%panels(d)%offset_u
        do k = 1, w
            r = self%pivots(first+k-1)
            if (r == k) cycle
            do i = 1, mu
                saved = self%zvalues(u0+int(k-1,i8)*mu+i)
                self%zvalues(u0+int(k-1,i8)*mu+i) = self%zvalues(u0+int(r-1,i8)*mu+i)
                self%zvalues(u0+int(r-1,i8)*mu+i) = saved
            end do
        end do
        if (self%use_blas .and. w >= 96 .and. mu >= 32) then
            call ztrsm('R', 'U', 'N', 'N', mu, w, zone, self%zvalues(l0+1), m, &
                self%zvalues(l0+w+1), m)
            call ztrsm('R', 'L', 'T', 'U', mu, w, zone, self%zvalues(l0+1), m, &
                self%zvalues(u0+1), mu)
        else
            call slu_z_trsm_upper(self%zvalues(l0+1), m, w, self%zvalues(l0+w+1), m, mu)
            call slu_z_trsm_unit_lower_t(self%zvalues(l0+1), m, w, &
                self%zvalues(u0+1), mu, mu)
        end if
        do e = u0+1, u0+int(mu,i8)*w
            max_u = max(max_u, abs(self%zvalues(e)))
        end do
    end subroutine panel_solves_z
end submodule fortsparse_native_slu_z_factor
