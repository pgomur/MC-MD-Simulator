! SIMULACIÓN DE DINÁMICA MOLECULAR 3D / 3D MOLECULAR DYNAMICS SIMULATION
! Autor: Pedro Gomez Urdiales / Author: Pedro Gomez Urdiales
! Fecha: 2025-09-12 / Date: 2025-09-12
! ===================================================================

! -----------------------------------------------------------------
! Monte Carlo simulation in the μVT ensemble with optional Metadynamics/Umbrella bias.
! Reads input parameters from the command line, initializes the system and RNG,
! applies Monte Carlo moves, adjusts step sizes periodically, and writes binary output.
!
! Simulación Monte Carlo en el ensamble μVT con posible sesgo de Metadynamics/Umbrella.
! Lee los parámetros de entrada desde la línea de comandos, inicializa el sistema y el RNG,
! aplica movimientos Monte Carlo, ajusta los pasos periódicamente y guarda la salida binaria.
! -----------------------------------------------------------------
program mc_muvt
    use iso_fortran_env, only: real64
    use system_mod
    use mc_mod
    use io_mod
    use rng_mod     
    implicit none

    !----------------------------------------
    ! Input parameters and simulation configuration
    ! Parámetros de entrada y configuración de la simulación
    !----------------------------------------
    integer :: nsteps, n, seed
    real(real64) :: T, V, mu
    real(real64) :: xi_target, k_bias
    character(len=1024) :: output_file
    character(len=1024) :: arg
    integer :: i, ios
    integer, parameter :: ADJUST_INTERVAL = 100
    integer, parameter :: LOG_INTERVAL = 500

    !----------------------------------------
    ! Default values
    ! Valores por defecto
    !----------------------------------------
    nsteps      = 1000
    xi_target   = 0.0_real64
    k_bias      = 0.0_real64
    sigma_meta  = 0.0_real64
    height_meta = 0.0_real64
    output_file = "output.bin"
    seed        = 0   ! Por defecto => reloj del sistema

    !----------------------------------------
    ! Check minimum required arguments
    ! Verifica argumentos mínimos requeridos
    !----------------------------------------
    if (command_argument_count() < 4) then
        print *, "Uso: mc_muvt N V T μ [nsteps output_file xi_target k_bias sigma_meta height_meta seed]"
        stop
    end if

    !----------------------------------------
    ! Read main parameters from command line
    ! Lee parámetros principales desde la línea de comandos
    !----------------------------------------
    call get_command_argument(1,arg); read(arg,*) n
    call get_command_argument(2,arg); read(arg,*) V
    call get_command_argument(3,arg); read(arg,*) T
    call get_command_argument(4,arg); read(arg,*) mu

    !----------------------------------------
    ! Read optional arguments safely
    ! Lee argumentos opcionales de manera segura
    !----------------------------------------
    if (command_argument_count() >= 5) then
        call get_command_argument(5,arg)
        read(arg,*,iostat=ios) nsteps
        if (ios /= 0) nsteps = 1000
    end if

    if (command_argument_count() >= 6) then
        call get_command_argument(6,arg)
        output_file = trim(arg)
    end if

    if (command_argument_count() >= 7) then
        call get_command_argument(7,arg)
        read(arg,*,iostat=ios) xi_target
        if (ios /= 0) xi_target = 0.0_real64
    end if

    if (command_argument_count() >= 8) then
        call get_command_argument(8,arg)
        read(arg,*,iostat=ios) k_bias
        if (ios /= 0) k_bias = 0.0_real64
    end if

    if (command_argument_count() >= 9) then
        call get_command_argument(9,arg)
        read(arg,*,iostat=ios) sigma_meta
        if (ios /= 0) sigma_meta = 0.0_real64
    end if

    if (command_argument_count() >= 10) then
        call get_command_argument(10,arg)
        read(arg,*,iostat=ios) height_meta
        if (ios /= 0) height_meta = 0.0_real64
    end if

    !----------------------------------------
    ! Initialize random number generator
    ! Inicializa generador de números aleatorios
    !----------------------------------------
    call set_seed()

    !----------------------------------------
    ! Initialize physical system
    ! Inicializa el sistema físico
    !----------------------------------------
    call init_system(n, V, T, mu)

    !----------------------------------------
    ! Setup Metadynamics/Umbrella bias if requested
    ! Configuración de Metadynamics/Umbrella si está activado
    !----------------------------------------
    if (k_bias > 0.0_real64) then
        if (sigma_meta <= 0.0_real64) sigma_meta = 0.1_real64
        if (height_meta <= 0.0_real64) height_meta = k_bias
        call set_metadynamics(sigma_meta, height_meta)
        print *, "Umbrella/Metadinámica activada"
        print *, " xi_target = ", xi_target, ", k_bias = ", k_bias
        print *, " Metadynamics sigma = ", sigma_meta, ", height = ", height_meta
    end if

    !----------------------------------------
    ! Initial simulation summary
    ! Resumen inicial de la simulación
    !----------------------------------------
    print *, "Simulación μVT Monte Carlo iniciada"
    print *, "N = ", n_particles
    print *, "V = ", volume
    print *, "T = ", temperature
    print *, "μ = ", chemical_potential
    print *, "Número de pasos = ", nsteps
    if (meta_active) print *, "Metadynamics activo"
    print *, "Seed RNG = ", seed

    !----------------------------------------
    ! Main Monte Carlo loop
    ! Bucle principal Monte Carlo
    !----------------------------------------
    do i = 1, nsteps
        !---------------------------------------------------------
        ! Perform one Monte Carlo step
        ! Ejecuta un paso Monte Carlo
        !---------------------------------------------------------
        call monte_carlo_step()

        !---------------------------------------------------------
        ! Adjust maximum displacement periodically
        ! Ajusta desplazamiento máximo periódicamente
        !---------------------------------------------------------
        if (mod(i, ADJUST_INTERVAL) == 0 .and. i < nsteps) then
            call adjust_max_disp()
        end if

        !---------------------------------------------------------
        ! Periodic logging of acceptance statistics
        ! Registro periódico de tasas de aceptación
        !---------------------------------------------------------
        if (mod(i, LOG_INTERVAL) == 0) then
            print *, "Paso ", i, " de ", nsteps
            call print_acceptance()
        end if
    end do

    !----------------------------------------
    ! Final acceptance statistics
    ! Estadísticas finales de aceptación
    !----------------------------------------
    print *, "Acceptance rates (final window):"
    call print_acceptance()

    !----------------------------------------
    ! Write simulation results to binary file
    ! Guarda resultados de la simulación en archivo binario
    !----------------------------------------
    if (meta_active) then
        call write_binary(output_file, positions, n_particles, energy_total, density, max_disp, &
                          .true., xi_target, k_bias, .true., sigma_meta, height_meta, meta_bias)
    else
        call write_binary(output_file, positions, n_particles, energy_total, density, max_disp)
    end if

    print *, "Simulación finalizada. Datos guardados en ", trim(output_file)

end program mc_muvt
