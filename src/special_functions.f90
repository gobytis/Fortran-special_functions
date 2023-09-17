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
    !> `beta(a, b)` = $int_0^1 t^{a-1} (1 - t)^{a - 1} dt$
    interface beta
        pure module function int_beta(a, b) result(retval)
            !> argument 1
            integer, intent(in) :: a
            !> argument 2
            integer, intent(in) :: b
            !> 戻り値
            real(dp) :: retval
        end function

        pure module function sp_beta(a, b) result(retval)
            !> argument 1
            real(sp), intent(in) :: a
            !> argument 2
            real(sp), intent(in) :: b
            !> 戻り値
            real(sp) :: retval
        end function

        pure module function dp_beta(a, b) result(retval)
            !> argument 1
            real(dp), intent(in) :: a
            !> argument 2
            real(dp), intent(in) :: b
            !> 戻り値
            real(dp) :: retval
        end function
    end interface


    !> ディガンマ関数用のラッパーです。
    interface digamma
        impure module function dp_digamma(x) result(retval)
            !> 入力
            real(dp), intent(in) :: x
            !> 戻り値
            real(dp) :: retval
        end function
    end interface


    !> トリガンマ関数用のラッパーです。
    interface trigamma
        impure module function dp_trigamma(x) result(retval)
            !> 入力
            real(dp), intent(in) :: x
            !> 戻り値
            real(dp) :: retval
        end function
    end interface


    !> ポリガンマ関数です。
    interface polygamma
        impure recursive module function dp_polygamma(k, x) result(retval)
            !> 関数パラメータ
            !> k = 0, 1, 2, corresponds to digamma, trigamma, tetragamma.... respectively.
            integer, intent(in) :: k
            !> 入力
            real(dp), intent(in) :: x
            !> 戻り値
            real(dp) :: retval
        end function
    end interface

end module
