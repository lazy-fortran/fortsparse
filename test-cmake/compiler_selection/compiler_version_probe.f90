program compiler_version_probe
    use, intrinsic :: iso_fortran_env, only: compiler_version
    implicit none

    print '(a)', compiler_version()
end program compiler_version_probe
