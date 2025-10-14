!===============================================================
! Handles the global system state for Monte Carlo μVT simulations
!
! Maneja el estado global del sistema para simulaciones Monte Carlo μVT
!===============================================================
module system_mod
    use iso_fortran_env, only: real64
    use rng_mod, only: random_uniform, random_vector
    implicit none

    private
    public :: init_system, update_system
    public :: n_particles, volume, temperature, chemical_potential
    public :: positions, max_disp, energy_total, density

    !----------------------------------------
    ! Variables globales del sistema
    !----------------------------------------
    integer :: n_particles                   ! Number of particles / Número de partículas
    real(real64) :: volume                   ! System volume / Volumen del sistema
    real(real64) :: temperature              ! System temperature / Temperatura del sistema
    real(real64) :: chemical_potential       ! Chemical potential / Potencial químico
    real(real64) :: max_disp                 ! Maximum MC displacement / Desplazamiento máximo de MC
    real(real64) :: energy_total             ! Total system energy / Energía total del sistema
    real(real64) :: density                  ! System density / Densidad del sistema
    real(real64), allocatable :: positions(:,:)  ! Particle positions (3 x N) / Posiciones de partículas (3 x N)

contains

    !----------------------------------------
    ! Initializes the physical system for μVT Monte Carlo
    !
    ! Inicializa el sistema físico para Monte Carlo μVT
    !----------------------------------------
    subroutine init_system(N, V, T, mu)
        integer, intent(in) :: N
        real(real64), intent(in) :: V, T, mu
        integer :: i
        real(real64) :: box_length

        !---------------------------------------------------
        ! Input validation / Validación de parámetros de entrada
        !---------------------------------------------------
        if (N <= 0) then
            error stop "init_system: El número de partículas debe ser mayor que 0"
        end if
        if (V <= 0.0_real64) then
            error stop "init_system: El volumen debe ser positivo"
        end if
        if (T <= 0.0_real64) then
            error stop "init_system: La temperatura debe ser positiva"
        end if

        !---------------------------------------------------
        ! Assign global variables / Asignación de variables globales
        !---------------------------------------------------
        n_particles = N
        volume = V
        temperature = T
        chemical_potential = mu
        max_disp = 0.1_real64  ! Desplazamiento máximo por defecto

        !---------------------------------------------------
        ! Allocate particle positions / Asignación de memoria para posiciones
        !---------------------------------------------------
        if (allocated(positions)) deallocate(positions)
        allocate(positions(3, n_particles))

        !---------------------------------------------------
        ! Initialize positions uniformly in the cubic box / Inicialización de posiciones dentro del cubo
        !---------------------------------------------------
        box_length = volume**(1.0_real64/3.0_real64)
        do i = 1, n_particles
            positions(:,i) = random_vector(3) * box_length
        end do

        !---------------------------------------------------
        ! Initialize total energy and density / Inicialización de energía total y densidad
        !---------------------------------------------------
        energy_total = 0.0_real64
        call update_system()
    end subroutine init_system

    !----------------------------------------
    ! Updates system-dependent properties (e.g., density)
    !
    ! Actualiza propiedades dependientes del sistema (por ejemplo, densidad)
    !----------------------------------------
    subroutine update_system()
        !---------------------------------------------------
        ! Update density safely / Actualiza densidad de manera segura
        !---------------------------------------------------
        if (volume > 0.0_real64) then
            density = real(n_particles, real64) / volume
        else
            error stop "update_system: Volumen no puede ser cero"
        end if
    end subroutine update_system

    !----------------------------------------
    ! Optional debug/monitoring function to display system properties
    !
    ! Función adicional opcional para depuración/monitoreo
    ! Muestra las propiedades actuales del sistema
    !----------------------------------------
    subroutine print_system_state()
        print '(A,I8)', "Número de partículas: ", n_particles
        print '(A,F12.6)', "Volumen: ", volume
        print '(A,F12.6)', "Temperatura: ", temperature
        print '(A,F12.6)', "Potencial químico: ", chemical_potential
        print '(A,F12.6)', "Desplazamiento máximo: ", max_disp
        print '(A,F12.6)', "Energía total: ", energy_total
        print '(A,F12.6)', "Densidad: ", density
    end subroutine print_system_state

end module system_mod
