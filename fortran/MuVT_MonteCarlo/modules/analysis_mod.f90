! -----------------------------------------------------------------
! Module for statistical analysis: computes averages and standard errors.
! Uses Kahan compensated summation to minimize numerical errors.
!
! Módulo para análisis estadístico: calcula promedios y errores estándar.
! Utiliza suma compensada de Kahan para minimizar errores numéricos.
! -----------------------------------------------------------------
module analysis_mod
    use iso_fortran_env, only: real64
    implicit none

    private
    public :: compute_average
    public :: compute_std_error

contains

    !---------------------------------------------------------------
    ! Compute the average of a 1D array using Kahan compensated summation.
    ! Provides higher numerical accuracy for large datasets.
    !
    ! Calcula el promedio de un array 1D usando suma compensada de Kahan.
    ! Proporciona mayor precisión numérica para datasets grandes.
    !---------------------------------------------------------------
    function compute_average(data) result(avg)
        real(real64), intent(in) :: data(:)
        real(real64) :: avg, c, y, t
        integer :: i, n

        n = size(data)
        if (n == 0) then
            error stop "compute_average: el array de entrada está vacío"
        end if

        avg = 0.0_real64
        c = 0.0_real64
        do i = 1, n
            y = data(i) - c
            t = avg + y
            c = (t - avg) - y
            avg = t
        end do

        avg = avg / real(n, real64)
    end function compute_average

    !---------------------------------------------------------------
    ! Compute the standard error of the mean for a 1D array using Kahan summation.
    ! Requires at least 2 elements. Uses compensated summation for accurate variance.
    !
    ! Calcula el error estándar del promedio para un array 1D usando suma compensada de Kahan.
    ! Requiere al menos 2 elementos. Usa suma compensada para mayor precisión de la varianza.
    !---------------------------------------------------------------
    function compute_std_error(data) result(err)
        real(real64), intent(in) :: data(:)
        real(real64) :: err, avg, diff, sum_sq, c, y, t
        integer :: i, n

        n = size(data)
        if (n < 2) then
            error stop "compute_std_error: se requieren al menos 2 elementos"
        end if

        ! Promedio con Kahan
        avg = compute_average(data)

        ! Suma de cuadrados con Kahan
        sum_sq = 0.0_real64
        c = 0.0_real64
        do i = 1, n
            diff = data(i) - avg
            y = diff*diff - c
            t = sum_sq + y
            c = (t - sum_sq) - y
            sum_sq = t
        end do

        err = sqrt(sum_sq / real(n*(n-1), real64))
    end function compute_std_error

end module analysis_mod
