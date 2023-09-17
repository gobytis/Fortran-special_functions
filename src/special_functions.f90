!> 特殊関数を格納するクラスです。
module special_functions
    use, intrinsic :: iso_fortran_env
    implicit none

    private

    ! 変数の kind
    integer, parameter :: sp = real32
    integer, parameter :: dp = real64

    ! public 関数一覧
    public :: beta
    public :: digamma, trigamma, polygamma

    !> ベータ関数です。
    interface beta
        pure module function int_beta(a, b) result(retval)
            integer, intent(in) :: a
            integer, intent(in) :: b
            real(dp) :: retval
        end function

        pure module function sp_beta(a, b) result(retval)
            real(sp), intent(in) :: a
            real(sp), intent(in) :: b
            real(sp) :: retval
        end function

        pure module function dp_beta(a, b) result(retval)
            real(dp), intent(in) :: a
            real(dp), intent(in) :: b
            real(dp) :: retval
        end function
    end interface


    !> ディガンマ関数用のラッパーです。
    interface digamma
        impure module function dp_digamma(x) result(retval)
            real(dp), intent(in) :: x
            real(dp) :: retval
        end function
    end interface


    !> トリガンマ関数用のラッパーです。
    interface trigamma
        impure module function dp_trigamma(x) result(retval)
            real(dp), intent(in) :: x
            real(dp) :: retval
        end function
    end interface


    !> ポリガンマ関数です。
    interface polygamma
        impure recursive module function dp_polygamma(k, x) result(retval)
            integer, intent(in) :: k
            real(dp), intent(in) :: x
            real(dp) :: retval
        end function
    end interface

end module
