! ===================================================================
! SIMULACIÓN DE DINÁMICA MOLECULAR 3D / 3D MOLECULAR DYNAMICS SIMULATION
! Potencial: Lennard-Jones 12-6 con shift / Potential: Lennard-Jones 12-6 with shift
! Integrador: Langevin BAOAB completo / Integrator: Full Langevin BAOAB
! Cell lists O(n) / O(n) cell lists
! Paralelizado con OpenMP / Parallelized with OpenMP
! Salida binaria / Binary output
! Autor: Pedro Gomez Urdiales / Author: Pedro Gomez Urdiales
! Fecha: 2025-09-12 / Date: 2025-09-12
! ===================================================================

! -----------------------------------------------------------------
! Main program: Molecular Dynamics simulation using BAOAB integrator 
! and cell lists with OpenMP parallelization.
! Reads parameters from the command line, initializes the system, 
! computes forces, and writes binary output frames at each step.
!
! Programa principal: Simulación de Dinámica Molecular con integrador BAOAB 
! y listas de celdas paralelizadas con OpenMP.
! Lee los parámetros desde la línea de comandos, inicializa el sistema, 
! calcula fuerzas y escribe frames binarios en cada paso.
! -----------------------------------------------------------------
program md_baoab_cell_omp
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    include 'omp_lib.h'

    ! -----------------------------------------------------------------
    ! Default simulation parameters
    ! Parámetros por defecto de la simulación
    ! -----------------------------------------------------------------
    integer(int32), parameter :: DEFAULT_N = 10000
    integer(int32), parameter :: DEFAULT_STEPS = 1000
    real(real64),   parameter :: DEFAULT_DT = 0.005_real64
    real(real64),   parameter :: DEFAULT_BOX = 50.0_real64
    real(real64),   parameter :: DEFAULT_CUTOFF = 2.5_real64
    real(real64),   parameter :: DEFAULT_SKIN = 0.3_real64
    integer(int32), parameter :: DEFAULT_SEED = 12345
    integer(int32), parameter :: DEFAULT_NTYPES = 2

    ! -----------------------------------------------------------------
    ! Main simulation variables
    ! Variables principales de simulación
    ! -----------------------------------------------------------------
    integer(int32) :: n, steps, seed, ntypes
    real(real64) :: dt, box, cutoff, skin, kT, gamma
    integer(int32) :: i

    ! -----------------------------------------------------------------
    ! Dynamic arrays for particle properties
    ! Arreglos dinámicos para las propiedades de las partículas
    ! -----------------------------------------------------------------
    real(real64), allocatable :: x(:,:), v(:,:), f(:,:)
    real(real64), allocatable :: mass(:), sigma(:), epsilon(:)
    real(real64), allocatable :: sigma6(:,:), epsilon_mat(:,:)

    ! -----------------------------------------------------------------
    ! Cell list data structures for neighbor searching
    ! Estructuras de listas de celdas para búsqueda de vecinos
    ! -----------------------------------------------------------------
    integer, allocatable :: head(:), linked_list(:)
    integer :: ncell, nx, ny, nz
    real(real64) :: cell_size

    ! -----------------------------------------------------------------
    ! Energy variables (potential, kinetic, total)
    ! Variables de energía (potencial, cinética, total)
    ! -----------------------------------------------------------------
    real(real64) :: pe, ke, total

    ! -----------------------------------------------------------------
    ! Binary output file unit
    ! Unidad de archivo binario de salida
    ! -----------------------------------------------------------------
    integer :: bin_unit = 10

    ! -----------------------------------------------------------------
    ! Read command-line parameters or use defaults
    ! Lee parámetros desde la línea de comandos o usa valores por defecto
    ! -----------------------------------------------------------------
    call read_command_line_arguments(n, steps, dt, box, cutoff, skin, seed, kT, gamma)

    ntypes = DEFAULT_NTYPES

    ! -----------------------------------------------------------------
    ! System initialization: memory allocation, setup, and interactions
    ! Inicialización del sistema: reserva de memoria, configuración e interacciones
    ! -----------------------------------------------------------------
    call allocate_arrays(n, ntypes)
    call initialize_system(n, ntypes, box, kT, seed)
    call compute_interaction_matrices(ntypes, sigma, epsilon, sigma6, epsilon_mat)

    ! -----------------------------------------------------------------
    ! Cell list setup based on cutoff and skin distance
    ! Configuración de la lista de celdas según el cutoff y el skin
    ! -----------------------------------------------------------------
    cell_size = cutoff + skin
    nx = max(1,int(box/cell_size))
    ny = nx
    nz = nx
    ncell = nx*ny*nz
    allocate(head(ncell))
    allocate(linked_list(n))

    ! -----------------------------------------------------------------
    ! Open binary output file for trajectory storage
    ! Abre archivo binario para almacenar la trayectoria
    ! -----------------------------------------------------------------
    open(unit=bin_unit, file='simulation_NVT_BAOAB.bin', status='replace', form='unformatted')

    ! -----------------------------------------------------------------
    ! Main simulation loop
    ! Bucle principal de simulación
    ! -----------------------------------------------------------------
    do i = 1, steps
        ! -------------------------------------------------------------
        ! Build cell list for neighbor interactions
        ! Construye lista de celdas para las interacciones vecinas
        ! -------------------------------------------------------------
        call build_cell_list(n, x, nx, ny, nz, box, head, linked_list)

        ! -------------------------------------------------------------
        ! Integrate particle motion using BAOAB scheme (thread-safe RNG)
        ! Integra el movimiento de las partículas con el esquema BAOAB (RNG seguro por hilo)
        ! -------------------------------------------------------------
        call baoab_step_omp(n, dt, gamma, kT, mass, x, v, f, box, seed)

        ! -------------------------------------------------------------
        ! Compute interparticle forces using cell list (parallel)
        ! Calcula las fuerzas entre partículas usando listas de celdas (paralelo)
        ! -------------------------------------------------------------
        call compute_forces_cell_omp(n, ntypes, x, f, pe, sigma6, epsilon_mat, nx, ny, nz, head, linked_list, box, cutoff)

        ! -------------------------------------------------------------
        ! Compute kinetic energy and total energy
        ! Calcula la energía cinética y la energía total
        ! -------------------------------------------------------------
        call compute_kinetic_energy_omp(n, mass, v, ke)
        total = ke + pe

        ! -------------------------------------------------------------
        ! Write current simulation frame to binary file
        ! Escribe el frame actual de la simulación en el archivo binario
        ! -------------------------------------------------------------
        call write_frame_binary(bin_unit, n, x, v, f, ke, pe, total)
    end do
    ! -----------------------------------------------------------------
    ! Close output file and deallocate all arrays
    ! Cierra el archivo de salida y libera toda la memoria
    ! -----------------------------------------------------------------
    close(bin_unit)
    deallocate(x, v, f, mass, sigma, epsilon, sigma6, epsilon_mat, head, linked_list)

contains

! -----------------------------------------------------------------
! Read simulation parameters from command-line arguments.
! If an argument is missing or invalid, default values are assigned.
! Supports basic and thermodynamic parameters (kT, gamma).
!
! Lee los parámetros de simulación desde los argumentos de línea de comandos.
! Si falta un argumento o es inválido, se asignan valores por defecto.
! Soporta parámetros básicos y termodinámicos (kT, gamma).
! -----------------------------------------------------------------
subroutine read_command_line_arguments(n, steps, dt, box, cutoff, skin, seed, kT, gamma)
    integer(int32), intent(out) :: n, steps, seed
    real(real64), intent(out) :: dt, box, cutoff, skin, kT, gamma
    character(len=32) :: arg
    integer :: ios

    ! Lectura de parámetros básicos
    call get_command_argument(1, arg); read(arg,*,iostat=ios) n; if(ios/=0) n=DEFAULT_N
    call get_command_argument(2, arg); read(arg,*,iostat=ios) steps; if(ios/=0) steps=DEFAULT_STEPS
    call get_command_argument(3, arg); read(arg,*,iostat=ios) dt; if(ios/=0) dt=DEFAULT_DT
    call get_command_argument(4, arg); read(arg,*,iostat=ios) box; if(ios/=0) box=DEFAULT_BOX
    call get_command_argument(5, arg); read(arg,*,iostat=ios) cutoff; if(ios/=0) cutoff=DEFAULT_CUTOFF
    call get_command_argument(6, arg); read(arg,*,iostat=ios) skin; if(ios/=0) skin=DEFAULT_SKIN
    call get_command_argument(7, arg); read(arg,*,iostat=ios) seed; if(ios/=0) seed=DEFAULT_SEED

    ! Lectura de parámetros termodinámicos (opcional)
    call get_command_argument(8, arg); read(arg,*,iostat=ios) kT; if(ios/=0) kT=1.0_real64
    call get_command_argument(9, arg); read(arg,*,iostat=ios) gamma; if(ios/=0) gamma=0.1_real64
end subroutine read_command_line_arguments

! -----------------------------------------------------------------
! Allocate all main simulation arrays for particles and interaction parameters.
! Includes positions, velocities, forces, masses, and Lennard-Jones matrices.
!
! Reserva todos los arreglos principales de la simulación para partículas y parámetros de interacción.
! Incluye posiciones, velocidades, fuerzas, masas y matrices de Lennard-Jones.
! -----------------------------------------------------------------
subroutine allocate_arrays(n, ntypes)
    integer, intent(in) :: n, ntypes
    allocate(x(3,n), v(3,n), f(3,n))
    allocate(mass(n), sigma(n), epsilon(n))
    allocate(sigma6(ntypes,ntypes), epsilon_mat(ntypes,ntypes))
end subroutine allocate_arrays

! -----------------------------------------------------------------
! Initialize the particle system with positions, velocities, masses, and types.
! Positions are placed on a lattice with small random perturbations,
! ensuring a minimum interparticle distance. Velocities follow Maxwell-Boltzmann.
! Random seed initialization ensures reproducibility.
!
! Inicializa el sistema de partículas con posiciones, velocidades, masas y tipos.
! Las posiciones se colocan en una rejilla con pequeñas perturbaciones aleatorias,
! asegurando una distancia mínima entre partículas. Las velocidades siguen Maxwell-Boltzmann.
! La inicialización de la semilla aleatoria garantiza reproducibilidad.
! -----------------------------------------------------------------
subroutine initialize_system(n, ntypes, box, kT, seed)
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none
    integer, intent(in) :: n, ntypes
    real(real64), intent(in) :: box, kT
    integer, intent(in) :: seed

    integer :: i, j, nx, ny, nz, ix, iy, iz, rsize
    integer, allocatable :: seed_array(:)
    real(real64) :: rnd(3), dx(3), min_dist, scale
    real(real64) :: spacing

    ! -----------------------------------------------------------------
    ! Initialize random seed for reproducible stochastic behavior
    ! Inicializa la semilla aleatoria para garantizar reproducibilidad
    ! -----------------------------------------------------------------
    call random_seed(size=rsize)
    allocate(seed_array(rsize))
    seed_array = seed
    call random_seed(put=seed_array)
    deallocate(seed_array)

    ! -----------------------------------------------------------------
    ! Assign particle masses and types
    ! Asigna masas y tipos de partícula
    ! -----------------------------------------------------------------
    do i = 1, n
        mass(i) = 1.0_real64 + 0.5_real64*mod(i-1, ntypes)
        if (mod(i-1, ntypes)==0) then
            sigma(i) = 1.0_real64
            epsilon(i) = 1.0_real64
        else
            sigma(i) = 1.2_real64
            epsilon(i) = 1.2_real64
        end if
    end do

    ! -----------------------------------------------------------------
    ! Generate initial lattice configuration with small random perturbations
    ! Genera configuración inicial en rejilla con pequeñas perturbaciones aleatorias
    ! -----------------------------------------------------------------
    nx = int(n**(1.0_real64/3.0_real64)) + 1
    ny = nx
    nz = nx
    spacing = box / real(nx, real64)  ! base spacing between particles / separación base entre partículas
    min_dist = 0.8_real64             ! minimum allowed distance / distancia mínima permitida

    i = 0
    do ix = 0, nx-1
        do iy = 0, ny-1
            do iz = 0, nz-1
                i = i + 1
                if (i > n) exit
                ! Colocar en rejilla con pequeña perturbación
                call random_number(rnd)
                scale = 0.1_real64 * spacing
                x(1,i) = (ix + 0.5_real64 + (rnd(1)-0.5_real64)*0.2_real64)*spacing
                x(2,i) = (iy + 0.5_real64 + (rnd(2)-0.5_real64)*0.2_real64)*spacing
                x(3,i) = (iz + 0.5_real64 + (rnd(3)-0.5_real64)*0.2_real64)*spacing
            end do
        end do
    end do

    ! -----------------------------------------------------------------
    ! Enforce minimum interparticle distance by repositioning overlaps
    ! Asegura distancia mínima entre partículas reposicionando solapamientos
    ! -----------------------------------------------------------------
    do i = 2, n
        do j = 1, i-1
            dx = x(:,i) - x(:,j)
            dx = dx - box*nint(dx/box)
            if (sqrt(sum(dx**2)) < min_dist) then
                ! Mover ligeramente la partícula i
                call random_number(rnd)
                x(:,i) = x(:,i) + (rnd - 0.5_real64)*min_dist
            end if
        end do
    end do

    ! -----------------------------------------------------------------
    ! Initialize velocities from Maxwell-Boltzmann distribution
    ! Inicializa velocidades según la distribución de Maxwell-Boltzmann
    ! -----------------------------------------------------------------
    do i = 1, n
        call random_number(rnd)
        v(:,i) = (rnd - 0.5_real64)*sqrt(2.0_real64*kT/mass(i))
    end do

end subroutine initialize_system

! -----------------------------------------------------------------
! Compute interaction matrices for all particle type pairs.
! Calculates sigma^6 for Lennard-Jones and the geometric mean of epsilon.
!
! Calcula las matrices de interacción para todos los pares de tipos de partículas.
! Calcula sigma^6 para Lennard-Jones y la media geométrica de epsilon.
! -----------------------------------------------------------------
subroutine compute_interaction_matrices(ntypes, sigma, epsilon, sigma6, epsilon_mat)
    integer, intent(in) :: ntypes
    real(real64), intent(in) :: sigma(:), epsilon(:)
    real(real64), intent(out) :: sigma6(:,:), epsilon_mat(:,:)
    integer :: i,j
    do i=1,ntypes
        do j=1,ntypes
            sigma6(i,j) = ((sigma(i)+sigma(j))*0.5_real64)**6
            epsilon_mat(i,j) = sqrt(epsilon(i)*epsilon(j))
        end do
    end do
end subroutine compute_interaction_matrices

! -----------------------------------------------------------------
! Compute the total kinetic energy in parallel using OpenMP.
! Each particle's contribution is summed with a reduction for thread safety.
!
! Calcula la energía cinética total en paralelo usando OpenMP.
! La contribución de cada partícula se suma mediante una reducción para seguridad entre hilos.
! -----------------------------------------------------------------
subroutine compute_kinetic_energy_omp(n, mass, v, ke)
    integer, intent(in) :: n
    real(real64), intent(in) :: mass(:)
    real(real64), intent(in) :: v(3,n)
    real(real64), intent(out) :: ke
    integer :: i
    ke = 0.0_real64
    !$omp parallel do reduction(+:ke)
    do i=1,n
        ke = ke + 0.5_real64*mass(i)*sum(v(:,i)**2)
    end do
    !$omp end parallel do
end subroutine compute_kinetic_energy_omp

! -----------------------------------------------------------------
! Write a simulation frame to a binary file.
! Outputs positions, velocities, forces, kinetic, potential, and total energy.
!
! Escribe un frame de la simulación en un archivo binario.
! Salidas: posiciones, velocidades, fuerzas, energía cinética, potencial y total.
! -----------------------------------------------------------------
subroutine write_frame_binary(unit, n, x, v, f, ke, pe, total)
    integer, intent(in) :: unit, n
    real(real64), intent(in) :: x(3,n), v(3,n), f(3,n), ke, pe, total
    write(unit) x, v, f, ke, pe, total
end subroutine write_frame_binary

! -----------------------------------------------------------------
! Generate a normally distributed random number using the Box-Muller transform.
! Inputs u1 and u2 must be uniform random numbers in (0,1).
!
! Genera un número aleatorio con distribución normal usando la transformación Box-Muller.
! Las entradas u1 y u2 deben ser números aleatorios uniformes en (0,1).
! -----------------------------------------------------------------
function box_muller(u1, u2) result(z)
    real(real64), intent(in) :: u1, u2
    real(real64) :: z
    real(real64), parameter :: twopi = 6.283185307179586_real64
    z = sqrt(-2.0_real64*log(u1)) * cos(twopi*u2)
end function box_muller

! -----------------------------------------------------------------
! Perform one BAOAB Langevin integration step in parallel.
! Thread-safe random numbers ensure reproducibility across threads.
! Applies periodic boundary conditions and velocity updates in BAOAB order.
!
! Ejecuta un paso de integración Langevin BAOAB en paralelo.
! Números aleatorios thread-safe garantizan reproducibilidad entre hilos.
! Aplica condiciones de contorno periódicas y actualizaciones de velocidad en el orden BAOAB.
! -----------------------------------------------------------------
subroutine baoab_step_omp(n, dt, gamma, kT, mass, x, v, f, box, seed)
    use iso_fortran_env, only: real64
    use omp_lib
    implicit none

    integer, intent(in) :: n, seed
    real(real64), intent(in) :: dt, gamma, kT, mass(:), box
    real(real64), intent(inout) :: x(3,n), v(3,n), f(3,n)

    integer :: i, tid, nthreads
    real(real64) :: expg, c1, c2, z(3), rnd(2)
    integer :: rng_size
    integer, allocatable :: seed_local(:)
    
    ! Parámetros BAOAB
    expg = exp(-gamma*dt)
    c1   = expg
    c2   = sqrt(1.0_real64 - expg**2)

    ! Obtener número de hilos
    !$omp parallel
    nthreads = omp_get_num_threads()
    !$omp end parallel

    ! Inicializar semillas locales para cada hilo (una vez por paso)
    !$omp parallel private(tid,rng_size,seed_local)
    tid = omp_get_thread_num() + 1

    ! Inicializamos el RNG del hilo de forma reproducible
    call random_seed(size=rng_size)
    allocate(seed_local(rng_size))
    seed_local = seed + tid*1234  ! offset estable por hilo
    call random_seed(put=seed_local)

    !$omp do schedule(static)
    do i = 1, n
        ! ----- Paso B: Langevin aleatorio -----
        call random_number(rnd); z(1)=box_muller(rnd(1), rnd(2))
        call random_number(rnd); z(2)=box_muller(rnd(1), rnd(2))
        call random_number(rnd); z(3)=box_muller(rnd(1), rnd(2))
        z = z * sqrt(kT/mass(i))
        v(:,i) = c1*v(:,i) + c2*z

        ! ----- Paso A: mitad de kick -----
        v(:,i) = v(:,i) + 0.5_real64*f(:,i)/mass(i)*dt

        ! ----- Paso O: drift completo -----
        x(:,i) = x(:,i) + v(:,i)*dt
        x(:,i) = modulo(x(:,i), box)  ! PBC

        ! ----- Kick final (A) -----
        v(:,i) = v(:,i) + 0.5_real64*f(:,i)/mass(i)*dt
    end do
    !$omp end do

    deallocate(seed_local)
    !$omp end parallel

end subroutine baoab_step_omp

! -----------------------------------------------------------------
! Build cell lists for efficient neighbor searching.
! Uses head and linked_list arrays to store particles in each cell.
! Safe for parallel force computation since no forces are computed here.
!
! Construye listas de celdas para búsqueda eficiente de vecinos.
! Usa los arrays head y linked_list para almacenar partículas en cada celda.
! Seguro para paralelización en el cálculo de fuerzas, ya que aquí no se calculan fuerzas.
! -----------------------------------------------------------------
subroutine build_cell_list(n, x, nx, ny, nz, box, head, linked_list)
    integer, intent(in) :: n, nx, ny, nz
    real(real64), intent(in) :: x(3,n), box
    integer, intent(out) :: head(:), linked_list(:)
    integer :: i, ix, iy, iz, cell

    head = 0
    linked_list = 0

    do i=1,n
        ix = modulo(int(x(1,i)/box*nx), nx)
        iy = modulo(int(x(2,i)/box*ny), ny)
        iz = modulo(int(x(3,i)/box*nz), nz)
        cell = ix + nx*iy + nx*ny*iz + 1
        linked_list(i) = head(cell)
        head(cell) = i
    end do
end subroutine build_cell_list

! -----------------------------------------------------------------
! Compute pairwise forces and potential energy using cell lists.
! Thread-local arrays are used to ensure safe OpenMP parallelization.
! Energy shift applied to smoothly truncate the potential at cutoff.
!
! Calcula fuerzas y energía potencial entre pares usando listas de celdas.
! Se usan arrays locales por hilo para garantizar seguridad en paralelización con OpenMP.
! Se aplica un shift de energía para truncar suavemente el potencial en el cutoff.
! -----------------------------------------------------------------
subroutine compute_forces_cell_omp(n, ntypes, x, f, pe, sigma6, epsilon_mat, &
                                   nx, ny, nz, head, linked_list, box, cutoff)
    use iso_fortran_env, only: real64
    use omp_lib
    implicit none

    integer, intent(in)    :: n, ntypes, nx, ny, nz
    real(real64), intent(in) :: x(3,n)
    real(real64), intent(in) :: sigma6(ntypes,ntypes), epsilon_mat(ntypes,ntypes)
    real(real64), intent(in) :: box, cutoff
    integer, intent(in)    :: head(:), linked_list(:)

    real(real64), intent(out) :: f(3,n), pe
    integer :: i, ti, tj, ix, iy, iz
    integer :: neigh_ix, neigh_iy, neigh_iz
    integer :: neighbor_cell, neighbor
    real(real64) :: dx(3), r2, r6, r12, inv_r2, fij(3), e_shift, rcut2
    integer :: tid, nthreads

    ! ---- Arrays locales por hilo ----
    real(real64), allocatable :: f_locals(:,:,:)
    real(real64), allocatable :: pe_locals(:)

    rcut2 = cutoff**2
    f  = 0.0_real64
    pe = 0.0_real64

    ! Número de hilos
    !$omp parallel
    nthreads = omp_get_num_threads()
    !$omp end parallel

    allocate(f_locals(3,n,nthreads))
    allocate(pe_locals(nthreads))
    f_locals = 0.0_real64
    pe_locals = 0.0_real64

    ! ---- Bucle paralelo con arrays locales por hilo ----
    !$omp parallel default(none) &
    !$omp shared(n,ntypes,x,sigma6,epsilon_mat,nx,ny,nz,head,linked_list,box,cutoff,rcut2,f_locals,pe_locals) &
    !$omp private(i,ti,tj,ix,iy,iz,neigh_ix,neigh_iy,neigh_iz,neighbor_cell,neighbor,dx,r2,r6,r12,inv_r2,fij,e_shift,tid)

    tid = omp_get_thread_num() + 1

    !$omp do schedule(dynamic)
    do i = 1, n
        ti = mod(i-1, ntypes) + 1
        ix = int(x(1,i)/box*nx)
        iy = int(x(2,i)/box*ny)
        iz = int(x(3,i)/box*nz)

        ! ---- Preparar arrays temporales vectorizables ----
        do neigh_ix = ix-1, ix+1
            do neigh_iy = iy-1, iy+1
                do neigh_iz = iz-1, iz+1
                    neighbor_cell = modulo(neigh_ix,nx) + nx*modulo(neigh_iy,ny) + nx*ny*modulo(neigh_iz,nz) + 1
                    neighbor = head(neighbor_cell)
                    ! ---- Vectorización posible: acumulamos fuerzas por bloque ----
                    do while(neighbor /= 0)
                        if (neighbor > i) then
                            tj = mod(neighbor-1, ntypes) + 1

                            ! Diferencia de posiciones (array-wide, vectorizable)
                            dx = x(:,i) - x(:,neighbor)
                            dx = dx - box*nint(dx/box)
                            r2 = sum(dx**2)

                            if (r2 < rcut2 .and. r2 > 0.0_real64) then
                                inv_r2 = 1.0_real64 / r2
                                r6  = sigma6(ti,tj) * inv_r2**3
                                r12 = r6**2
                                ! Fuerza vectorizada
                                fij = 24.0_real64*epsilon_mat(ti,tj)*(2.0_real64*r12 - r6)*inv_r2*dx

                                ! Shift de energía
                                e_shift = 4.0_real64*epsilon_mat(ti,tj) * &
                                          ((sigma6(ti,tj)/rcut2**3)**2 - sigma6(ti,tj)/rcut2**3)

                                ! Guardar en arrays locales por hilo
                                f_locals(:,i,tid)        = f_locals(:,i,tid) + fij
                                f_locals(:,neighbor,tid) = f_locals(:,neighbor,tid) - fij
                                pe_locals(tid) = pe_locals(tid) + 4.0_real64*epsilon_mat(ti,tj)*(r12 - r6) - e_shift
                            end if
                        end if
                        neighbor = linked_list(neighbor)
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel

    ! ---- Reducción final ----
    do tid = 1, nthreads
        f = f + f_locals(:,:,tid)
        pe = pe + pe_locals(tid)
    end do

    deallocate(f_locals)
    deallocate(pe_locals)

end subroutine compute_forces_cell_omp

end program md_baoab_cell_omp
