program test_fortsparse_compiler_self_assignment
    use fortsparse_compiler, only: require_safe_self_assignment
    implicit none

    call require_safe_self_assignment()
    print '(a)', 'PASS'
end program test_fortsparse_compiler_self_assignment
