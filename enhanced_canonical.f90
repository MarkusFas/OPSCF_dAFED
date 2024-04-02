module Physical_constants
      use, intrinsic :: iso_fortran_env
      implicit none
      real(real64), parameter :: pi = 3.1415926535897931 
      real(real64), parameter :: mass = 1.0
      real(real64), parameter :: w = 1.0
      real(real64), parameter :: hbar = 1.0
end module Physical_constants


module Simulation_parameters 
      use, intrinsic :: iso_fortran_env
      implicit none
      
      ! MD parameters
      real(real64) :: delta_t 
      integer :: eqsteps, numsteps, skipsteps 
      
      ! OPSCF parameters, time, beta, P,...
      real(real64) :: time, tau_sq, tau_mod
      integer :: P 
      integer, parameter :: M = 2
      real(real64) :: beta
      real(real64) :: wp 
      real(real64) :: gamma_
      character(len=128) :: potential_name

      ! d-AFED parameters
      real(real64) :: k_ext(2)
      real(real64) :: mass_ext(2)
      real(real64) :: beta_ext
      
      ! Parameters for data I/O
      character(len=128) :: OUT_FILE = "_trj_output.txt" 
      character(len=128) :: ENE_FILE = "_energy_output.txt" 
      character(len=128) :: FOR_FILE = "_force_output.txt" 
      integer :: index_

contains
      ! subroutine create output name
      !     concatenates the given index and output file names
      !     3 different output files are created: 
      !     OUT_FILE (containing trajactories)
      !     ENE_FILE (containing energies of the system and subsystems)
      !     FOR_FILE (containing forces for beads and extended system)
      subroutine create_output_name
            implicit none
            character(len=8) :: index_inplace
            if (index_>=100) then
                  write(index_inplace, '(i2)') int(index_)
            else if (index_>10) then
                  write(index_inplace, '(i3)') int(index_)
            else
                  write(index_inplace, '(i1)') int(index_)
            end if
            OUT_FILE = trim(index_inplace) // trim(OUT_FILE)
            OUT_FILE = trim(OUT_FILE)
            ENE_FILE = trim(index_inplace) // trim(ENE_FILE)
            ENE_FILE = trim(ENE_FILE)
            FOR_FILE = trim(index_inplace) // trim(FOR_FILE)
            FOR_FILE = trim(FOR_FILE)
      end subroutine create_output_name

      ! subroutine read_parameters:
      !     reads in parameters for:
      !     MD timesteps, OSPCF, d-AFED, I/O
      subroutine read_parameters()
        implicit none
        character(len=200) :: line
        character(len=20) :: key
        integer :: iostat

        ! Open the file
            do
                  read(*, '(A)', iostat=iostat) line ! Reading from standard input
                  if (iostat /= 0) exit ! Exit on end of file or error

                  ! Parse the line into key and value
                  read(line, *, iostat=iostat) key
                  if (iostat /= 0) then
                      print *, 'Error reading key from line: ', trim(line)
                      cycle
                  endif

                  ! Update module variables based on the key
                  select case (trim(key))

                  case ('timesteps')
                      read(line, *, iostat=iostat) key, delta_t
                      if (iostat /= 0) print *, 'Error reading timesteps from line: ', trim(line)
                  case ('eqsteps')
                      read(line, *, iostat=iostat) key, eqsteps
                      if (iostat /= 0) print *, 'Error reading eqsteps from line: ', trim(line)
                  case ('numsteps')
                      read(line, *, iostat=iostat) key, numsteps
                      if (iostat /= 0) print *, 'Error reading numsteps from line: ', trim(line)
                  case ('index')
                      read(line, *, iostat=iostat) key, index_
                      if (iostat /= 0) print *, 'Error reading index from line: ', trim(line)
                  case ('time')
                      read(line, *, iostat=iostat) key, time
                      if (iostat /= 0) print *, 'Error reading time from line: ', trim(line)
                  case ('P')
                      read(line, *, iostat=iostat) key, P
                      if (iostat /= 0) print *, 'Error reading P from line: ', trim(line)
                  case ('beta')
                      read(line, *, iostat=iostat) key, beta
                      if (iostat /= 0) print *, 'Error reading temperature from line: ', trim(line)
                  case ('potential')
                      read(line, *, iostat=iostat) key, potential_name
                      if (iostat /= 0) print *, 'Error reading potential from line: ', trim(line)
                  case ('dAFED_k1')
                      read(line, *, iostat=iostat) key, k_ext(1)
                      if (iostat /= 0) print *, 'Error reading temperature from line: ', trim(line)
                  case ('dAFED_k2')
                      read(line, *, iostat=iostat) key, k_ext(2)
                      if (iostat /= 0) print *, 'Error reading temperature from line: ', trim(line)
                  case ('dAFED_mass1')
                      read(line, *, iostat=iostat) key, mass_ext(1)
                      if (iostat /= 0) print *, 'Error reading temperature from line: ', trim(line)
                  case ('dAFED_mass2')
                      read(line, *, iostat=iostat) key, mass_ext(2)
                      if (iostat /= 0) print *, 'Error reading temperature from line: ', trim(line)
                  case ('dAFED_beta')
                      read(line, *, iostat=iostat) key, beta_ext
                      if (iostat /= 0) print *, 'Error reading temperature from line: ', trim(line)
                  case default
                      print *, 'Unknown parameter: ', trim(key)
                  end select
            end do
      end subroutine read_parameters

      ! subroutine assign_data:
      !     calls other methods and updates constants for OPSCF given the input
      subroutine assign_data
            use Physical_constants, only : hbar
            implicit none
            call read_parameters
            call create_output_name
            tau_sq = time*time + beta*beta*hbar*hbar/4
            tau_mod = sqrt(tau_sq)
            wp = sqrt(real(P))/tau_mod
            gamma_ = wp
      end subroutine assign_data

end module Simulation_parameters

! The TAMD extension module contains several functions related to properties unique to the extended
! d-AFED/TAMD system, such as energies or forces derived for it. The current implementation is fixed
! to 2 CVs.
module TAMD_extension
      use, intrinsic :: iso_fortran_env
      use Simulation_parameters, only : beta_ext, k_ext
      implicit none
      integer, parameter :: num_coll_var = 2
contains
      
      ! get_energy_extended:
      !     Takes bead positions of the open chain x(P+1) and exended variables z(num_coll_var) 
      !     and returns the energy from the harmonic coupling of CV to z.       
      function get_energy_extended(x,z)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P, mass_ext
            implicit none
            real(real64), dimension(P+1) :: x
            real(real64), dimension(num_coll_var) ::  z, cv
            real(real64) :: get_energy_extended
            cv = coll_var(x) 
            get_energy_extended = 0.5 * sum(k_ext * ((cv - z) * (cv - z)))
      end function get_energy_extended 

      ! get_forces_extendend:
      !     Takes bead positions of the open chain x(P+1) and exended variables z(num_coll_var) 
      !     and returns the forces on the extended variabels as array with dim = num_coll_var
      function get_forces_extended(x,z)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P
            implicit none
            real(real64), dimension(P+1) :: x
            real(real64), dimension(num_coll_var) :: get_forces_extended, z
            get_forces_extended = k_ext * (coll_var(x) - z)
      end function get_forces_extended


      ! coll_var(x):
      !      takes the open chain coordinates as input and returns an array containing 
      !      the collective variable of size num_coll_vars
      function coll_var(x)
            use, intrinsic :: iso_fortran_env
            use Physical_constants
            use Simulation_parameters, only : P
            implicit none
            real(real64), dimension(P+1) :: x
            real(real64), dimension(num_coll_var) :: coll_var
            coll_var(1) = 0.5 * (x(1) + x(P+1))
            coll_var(2) = x(1) - x(P+1)
            coll_var = coll_var
      end function coll_var
end module TAMD_extension

! The GGMT extension deals provides both a routine to calculate the required thermostat masses 
! for the GGMT and the energy of the GGMT variables.

module GGMT_extension 
      use, intrinsic :: iso_fortran_env
      use Simulation_parameters, only : M, beta_ext, k_ext, mass_ext, delta_t
      use Physical_constants, only : pi
      implicit none
      integer,dimension(M+1) :: C
      real(real64), dimension(M) :: A, Q_

contains
      ! The thermostat masses Q depend on the chosen timescale (tau_GGMT). 
      ! It is chosen depending on a reasonable timescale of the system. As this thermostat is only
      ! involved in the extended system, the timescale is solely depending on the coupling strength of the systems.
      ! Experience showed, that it should not be chosen too low. Parameter names follow from (Liu, Tuckerman , 2000)
      subroutine get_thermostat_masses
            implicit none
            integer :: i
            real(real64) :: tau_GGMT 
            tau_GGMT = 2*pi*sqrt(mass_ext(1)/k_ext(1))/5
            C = [1, 3, 15]
            A(1) = ((1/beta_ext) / C(1)) * (1)
            A(2) = (((1/beta_ext)**3) /C(2)) * (C(2)/C(1) + C(3)/C(2))
            Q_(1) = A(1) * (tau_GGMT**2)
            Q_(2) = A(2) * (tau_GGMT**2) 
      end subroutine get_thermostat_masses
      ! GGMT energy returns the energy for both thermostats on CV1 and CV2
      ! The same masses Q are chosen for the thermostats for both CVs.
      function GGMT_energy(n1,n2,pn1,pn2)
            use TAMD_extension, only : num_coll_var
            implicit none
            real(real64), dimension(num_coll_var) :: n1, n2, pn1, pn2
            real(real64) :: kinetic, potential, GGMT_energy
            kinetic = (sum(pn1*pn1)/Q_(1) + sum(pn2*pn2)/Q_(2))/2
            potential = (sum(n1) + sum(n2))/beta_ext
            GGMT_energy = kinetic + potential
      end function GGMT_energy
end module GGMT_extension

! The module Staging_transformation contains transformations from and to staging coordinates.
! The implementation follows the equations in Perez and Tuckerman, 2011
module Staging_transformation
      use Simulation_parameters, only : P
      implicit none
contains
      ! q_trafor(x) returns the staging coordinates q, obtained from open chain coordinates x
      function q_trafo(x)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64), dimension(P+1) :: q_trafo
            real(real64), dimension(P+1) :: x
            integer :: b
            q_trafo(1) = 0.5*(x(1) + x(P+1))
            do b=2,P
                  q_trafo(b) = x(b) - ((b-1)*x(b+1)+x(1))/b 
            end do
            q_trafo(P+1) = x(1) - x(P+1)
      end function q_trafo
      
      ! x_trafo(q) returns the open chain coordinates x from the staging coordinates q.
      function x_trafo(q)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64), dimension(P+1) :: q
            real(real64), dimension(P+1) :: x_trafo
            real(real64) :: temp
            integer :: t,s 
            x_trafo(1) = q(1) + 0.5*q(P+1)
            do s=2,P 
                  temp = 0
                  do t=s,P
                        temp = temp + q(t)/(t-1)
                  end do
                  x_trafo(s) = q(1) + q(P+1)*((P/2) - s + 1)/P + temp*(s-1) 
            end do
            x_trafo(P+1) = q(1) - 0.5*q(P+1)
      end function x_trafo
      
      ! force_trafo(Force_x) transforms the forces on the open chain coord x 
      ! to the forces acting on the staging variables
      function force_trafo(Force_x)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64), dimension(P+1) :: Force_x
            real(real64), dimension(P+1) :: force_trafo
            real(real64) :: F_temp, F_sum, force_exq
            integer :: b, s
            
            ! First adapt forces: add harmonic term or adapt prefactors:  
            F_temp = 0
            F_sum = 0
            do b = 2,P
                  F_sum = F_sum + Force_x(b)
                  F_temp = F_temp + Force_x(b)*((P/2) - b + 1)
            end do

            ! Force on first and last staging variable
            force_trafo(1) = (Force_x(1) + Force_x(P+1)) + F_sum
            force_trafo(P+1) = 0.5*(Force_x(1) - Force_x(P+1)) + F_temp/P

            ! For the other variables we calculate the forces iteratively
            do b=2,P
                  force_trafo(b) = Force_x(b) + (b-2)*force_trafo(b-1)/(b-1)
            end do
      end function force_trafo
      
      ! mass_trafo() returns the masses of the staging coordinates (dim = P+1) 
      ! given the mass of the particle
      function mass_trafo()
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : mass
            implicit none

            real(real64), dimension(P+1) :: mass_trafo
            integer :: s
            mass_trafo(1) = 0
            mass_trafo(P+1) = mass/P
            do s=2,P
                  mass_trafo(s) = (mass*s)/(s-1)
            end do
      end function mass_trafo
end module staging_transformation


module masses
      use, intrinsic :: iso_fortran_env
      use Physical_constants, only : mass
      use Simulation_parameters, only : P
      use Staging_transformation, only : mass_trafo
      implicit none
      real(real64), dimension(:), allocatable :: mass_staging, mass_staging_prime

contains
      ! the subroutine get_mass() updates the staging masses and the arbitrary masses for the staging 
      ! open chain simulation, denoted as mass_staging_prime (see Perez & Tuckerman, 2011)
      subroutine get_mass()
            implicit none
            allocate(mass_staging(P+1), mass_staging_prime(P+1))
            mass_staging = mass_trafo()
            mass_staging_prime = mass_staging
            mass_staging_prime(1) = mass
      end subroutine get_mass
end module masses

! The module Random_Variables provides a function to sample from a Gaussian distribution.
! For that purpose a Box-Müller transform is used.
module Random_Variables
      implicit none
contains
      ! Box-Müller transforms to sampe a RV from N(0,1): returns real
      function Gaussian_noise()
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : pi
            implicit none
            real(real64) :: n1,n2
            real(real64):: Gaussian_noise
            
            call random_number(n1)
            call random_number(n2)
            Gaussian_noise = sqrt(-2*log(n1))*cos(2*pi*n2)
      end function Gaussian_noise

      ! sample from a Gaussian with 0 mean and standard deviation (stdev): returns real RV
      function sample_Gaussian(mu, stdev)
            use, intrinsic :: iso_fortran_env
            implicit none      
            real(real64) :: mu
            real(real64) :: stdev
            real(real64) :: sample_Gaussian
            sample_Gaussian = mu + Gaussian_noise()*stdev
      end function sample_Gaussian
end module Random_Variables



module Forces_Potentials
      use, intrinsic :: iso_fortran_env
      !use Simulation_parameters, only : P
      use Random_Variables
      implicit none
      
      interface
            function potential_inter(x1) 
                  use, intrinsic :: iso_fortran_env
                  implicit none
                  real(real64) :: x1
                  real(real64) :: potential_inter
            end function potential_inter
      end interface
      
      procedure(potential_inter), pointer :: potential => null()
      procedure(potential_inter), pointer :: potential_df_1 => null()
      procedure(potential_inter), pointer :: potential_df_2 => null()
      procedure(potential_inter), pointer :: potential_df_3 => null()

contains
      
      ! Harmonic Potential and derivatives
      function harm_pot(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: harm_pot
            harm_pot = 0.5*x1*x1
      end function harm_pot

      function harm_pot_df1(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: harm_pot_df1
            harm_pot_df1 = x1
      end function harm_pot_df1

      function harm_pot_df2(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: harm_pot_df2
            harm_pot_df2 = 1
      end function harm_pot_df2

      function harm_pot_df3(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64):: harm_pot_df3
            harm_pot_df3 = 0
      end function harm_pot_df3
      
      ! Anharmonic Potential and derivatives
      function anharm_pot(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: anharm_pot
            anharm_pot = 0.01*x1*x1*x1*x1 + 0.1*x1*x1*x1 + 0.5*x1*x1
      end function anharm_pot

      function anharm_pot_df1(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: anharm_pot_df1
            anharm_pot_df1 = 0.04*x1*x1*x1 + 0.3*x1*x1 + x1
      end function anharm_pot_df1

      function anharm_pot_df2(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: anharm_pot_df2
            anharm_pot_df2 = 0.12*x1*x1 + 0.6*x1 + 1
      end function anharm_pot_df2

      function anharm_pot_df3(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64):: anharm_pot_df3
            anharm_pot_df3 = 0.24*x1 + 0.6
      end function anharm_pot_df3
      
      ! Quartic Potential and derivatives
      function quart_pot(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: quart_pot
            quart_pot = 0.25*x1*x1*x1*x1
      end function quart_pot

      function quart_pot_df1(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: quart_pot_df1
            quart_pot_df1 = x1*x1*x1
      end function quart_pot_df1

      function quart_pot_df2(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64) :: quart_pot_df2
            quart_pot_df2 = 3*x1*x1
      end function quart_pot_df2

      function quart_pot_df3(x1)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64) :: x1
            real(real64):: quart_pot_df3
            quart_pot_df3 = 6*x1
      end function quart_pot_df3
      
      ! To select the functions that should be used by the program
      ! This avoids 
      subroutine select_potential
            use, intrinsic :: iso_fortran_env
            use Simulation_Parameters, only : potential_name
            implicit none
            select case (potential_name)
            case ('harmonic')
                  potential => harm_pot
                  potential_df_1 => harm_pot_df1
                  potential_df_2 => harm_pot_df2
                  potential_df_3 => harm_pot_df3
            case ('anharmonic')
                  potential => anharm_pot
                  potential_df_1 => anharm_pot_df1
                  potential_df_2 => anharm_pot_df2
                  potential_df_3 => anharm_pot_df3
            case ('quartic')
                  potential => quart_pot
                  potential_df_1 => quart_pot_df1
                  potential_df_2 => quart_pot_df2
                  potential_df_3 => quart_pot_df3
            case default
                  print *, "Invalid choice, no function selected."
                  stop
            end select
      end subroutine select_potential

      ! get_forces(x,x_ext) returns the forces on the open chain coordinates. 
      ! It includes forces from the potential U, the effective OPSCF potential V, 
      ! and the coupling to the extended system. (the spring interbead forces will be computed
      ! in the staging coordinates)
      function get_forces(x,x_ext)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P, beta
            use TAMD_extension, only : num_coll_var, k_ext, coll_var
            implicit none
            real(real64) :: coll_var_term_1, coll_var_term_2
            real(real64), dimension(P+1) :: x
            real(real64), dimension(num_coll_var) :: x_ext
            real(real64), dimension(P+1) :: OPSCF_der, OPSCF_der__
            real(real64), dimension(P+1) :: get_forces
            real(real64), dimension(num_coll_var) :: collective_variables
            integer :: b

            !OPSCF_der: derivative of V
            OPSCF_der = - get_V_forces(x)/(2*beta)
            
            get_forces(1) = - 0.5 * potential_df_1(x(1)) / P
            get_forces(P+1) = - 0.5 * potential_df_1(x(P+1)) / P
            
            collective_variables = coll_var(x)
            ! ADDITIONAL TERMS DUE TO EXTENDED VARIABLES COUPLED TO COLLECTIVE VARIABLES
            coll_var_term_1 = collective_variables(1) - x_ext(1)
            coll_var_term_2 = collective_variables(2) - x_ext(2)
            
            get_forces(1) = get_forces(1) - k_ext(1) * ((coll_var_term_1 / 2) + coll_var_term_2)
            get_forces(P+1) = get_forces(P+1) - k_ext(2) * ((coll_var_term_1 / 2) - coll_var_term_2)
            
            ! REST OF THE COORDINATES (NOT ASSOCIATED WITH COLLECTIVE VARIABLES)
            do b=2,P
                  get_forces(b) = - potential_df_1(x(b)) / P
            end do
            get_forces = get_forces + OPSCF_der
      end function get_forces

      ! get_V_potential(x) returns the energy from the effective potential V
      ! It requires to compute the inverse and determinant of M. 
      function get_V_potential(r)
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : mass, hbar
            use simulation_parameters, only : P, tau_sq, time, beta
            implicit none
            
            real(real64) :: gamma_phasefactor, A, common_factor
            real(real64), dimension(P+1) :: r
            real(real64), dimension(P-1) :: diagonal
            real(real64), dimension(P-2) :: superdiagonal
            real(real64), dimension(P-1) :: K_vec
            real(real64), dimension(P-1, P-1) :: M_inv
            real(real64) :: M_det, KMK
            real(real64) :: first_term(P+1), testing(P+1), testing2(P+1), testing3(P+1)
            real(real64) :: get_V_potential
            integer :: b,i,j
      
            ! calculate parameters A, gammma
            common_factor = mass*P/tau_sq
            A = common_factor*beta/4
            gamma_phasefactor = common_factor*time/hbar  
            
            ! calculate the M matrix and K vector: 
            ! M can be described by a diagonal and sub or superdiagonal, since it is tridiagonal and symmetric 
            do b=2,P
                  K_vec(b-1) = gamma_phasefactor*(2*r(b) - r(b-1) - r(b+1)) - time*potential_df_1(r(b))/(P*hbar)
                  diagonal(b-1) = 2*A + beta*potential_df_2(r(b))/(4*P)
            end do
            superdiagonal = - A
           
            ! get the inverse M_inv and determinant M_det: 
            call Inv_Det(diagonal, superdiagonal, M_inv, M_det)
            
            ! compute the matrix product
            KMK = 0.0d0
            do i=1,P-1
                  do j=1,P-1
                        KMK = KMK + K_vec(j)*M_inv(j,i)*K_vec(i) 
                  end do
            end do
            get_V_potential = (KMK + log(M_det))/(2*beta)
      end function get_V_potential
     
      ! numerical_forces(r) computed the forces on the open chain coordinates numerically
      ! for that it used the finite difference between the potential energies at differetn times 
      function numerical_forces(r)
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : mass, hbar
            use simulation_parameters, only : P, tau_sq, time, beta, wp
            implicit none
            real(real64), dimension(P+1) :: r, r_var, r_var2
            real(real64), parameter :: dr = 0.0000000001
            real(real64) :: numerical_forces(P+1) 
            integer :: k
            do k=1,P+1
                  r_var = r
                  r_var2 = r
                  r_var(k) = r_var(k) + dr
                  r_var2(k) = r_var2(k) - dr
                  numerical_forces(k) = -(pot_energy(r_var,r) - pot_energy(r_var2,r))/(2*dr)
            end do
      end function numerical_forces

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!       V for OPSCF       !!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ! Evaluate the derivative of V wrt to one variable rk: except the first term ! 
      function evaluate_der(M_inv, K_vec, r, k)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P, beta
            implicit none
            real(real64) :: M_inv(P-1, P-1)
            real(real64) :: K_vec(P-1)
            real(real64) :: r
            integer :: j,i,k
            
            real(real64) :: evaluate_der
            evaluate_der = 0
            do j=1,P-1
                  do i=1,P-1
                        evaluate_der = evaluate_der - K_vec(i)*M_inv(i,k)*M_inv(k,j)*K_vec(j)
                  end do
            end do
            evaluate_der = evaluate_der * potential_df_3(r)*beta/(4*P)
            ! the detrminant atually rops out of the equation!! Only invrse
            evaluate_der = evaluate_der + M_inv(k,k) * beta*potential_df_3(r)/(4*P)
            !print *, M_inv(k,k) * beta*potential_df_3(r)/(4*P)

      end function evaluate_der
      
      function get_V_forces(r)
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : mass,  hbar
            use simulation_parameters, only : P, tau_sq, time, beta 
            implicit none
            real(real64) :: x
            
            real(real64) :: gamma_phasefactor, A, common_factor
            real(real64), dimension(P+1) :: r
            real(real64), dimension(P-1) :: diagonal
            real(real64), dimension(P-2) :: superdiagonal
            real(real64), dimension(P-1) :: K_vec
            real(real64), dimension(P-1, P-1) :: M_inv
            real(real64) :: M_det
            real(real64) :: first_term(P+1), testing(P+1), testing2(P+1), testing3(P+1)
            real(real64) :: get_V_forces(P+1)
            integer :: b,k,j
      
            ! calculate parameters A, gammma
            
            common_factor = mass*P/tau_sq
            A = common_factor*beta/4
            gamma_phasefactor = common_factor*time/hbar  
            
            ! calculate the M matrix and K vector
            do b=2,P
                  K_vec(b-1) = gamma_phasefactor*(2*r(b) - r(b-1) - r(b+1)) - time*potential_df_1(r(b))/(P*hbar)
                  diagonal(b-1) = 2*A + beta*potential_df_2(r(b))/(4*P)
            end do
            superdiagonal = - A
            !print *,gamma_phasefactor*(2*r(3) - r(3-1) - r(3+1))
            
            call Inv_Det(diagonal, superdiagonal, M_inv, M_det)
      
            !still need to do for k=1,P-1: only dependenc ethrough term in K!
            first_term = 0
            do j=1,P-1
                  ! only the term with k+1 = i survives which only leaves first matrix element
                  first_term(1) = first_term(1) - gamma_phasefactor*M_inv(1,j)*K_vec(j)
                  ! only the term with k-1 = i survives, which only leaves last matrix element
                  first_term(P+1) = first_term(P+1) - gamma_phasefactor*M_inv(P-1,j)*K_vec(j)

                  ! also the beads closest to ends, since different expression for the first term
                  first_term(2) = first_term(2) + gamma_phasefactor*(2*M_inv(1,j)*K_vec(j) - M_inv(2,j)*K_vec(j)) &
                                    - potential_df_2(r(2)) * M_inv(1,j) * K_vec(j) * time / (P*hbar) 
                  first_term(P) = first_term(P) + gamma_phasefactor*(2*M_inv(P-1,j)*K_vec(j) - M_inv(P-2,j)*K_vec(j)) &
                                    - potential_df_2(r(P)) * M_inv(P-1,j) * K_vec(j) * time / (P*hbar) 
            end do
           
            testing = 0
            testing2 = 0
            testing3 = 0
            get_V_forces = 0
            get_V_forces(2) = evaluate_der(M_inv, K_vec, r(2), 1) 
            get_V_forces(P) = evaluate_der(M_inv, K_vec, r(P), P-1)
            !$OMP PARALLEL DO
            do k=3,P-1 !was P-1 before?!?!?!? index of gradient: adapted to not run out of bounds with the off diagonal terms
                  do j=1,P-1
                        first_term(k) = first_term(k) +&
                        gamma_phasefactor*(2*M_inv(k-1,j)*K_vec(j))
                        testing(k) = testing(k) -  gamma_phasefactor * (M_inv(k-2,j)*K_vec(j))
                        testing3(k) = testing3(k) - gamma_phasefactor * (M_inv(k,j)*K_vec(j)) 
                        testing2(k) = testing2(k) - potential_df_2(r(k)) * M_inv(k-1,j) * K_vec(j) * time / (P*hbar)                                        
                  end do
            
                  get_V_forces(k) = evaluate_der(M_inv, K_vec, r(k), k-1)  
            end do
            !$OMP END PARALLEL DO
            !print *
            
            ! since the inverse has to be symmetric we just multiply this result by 2
            first_term = first_term + testing+ testing2 + testing3
            get_V_forces = get_V_forces + 2*first_term

      end function get_V_forces

      ! subroutine Inv_Det calculates the inverse and determinant for a given symmetric tridiagonal
      ! matrix, which is defined by its diagonal and superdiag. For performing a LU factorization, the Lapack
      ! method dgttrf optimized for tridiagonal matrices, is used. 
      subroutine Inv_Det(diagonal, superdiagonal, inv, determinant) 
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P  
            implicit none

            real(real64), intent(in) :: diagonal(P-1), superdiagonal(P-2)    ! Diagonal elements
            real(real64) :: subdiagonal(P-2) ! Subdiagonal elements
            real(real64) :: superdiagonal2(P-3) ! for the superdiagonal values of U
            integer :: i, j
            integer :: ipiv(P-1)
            real(real64), intent(out) :: inv(P-1,P-1) ! Symmetric tridiagonal matrix
            real(real64), intent(out) :: determinant
            integer :: info

            ! symmetric property:
            subdiagonal = superdiagonal 

            ! Perform LU factorization
            call dgttrf(P-1, subdiagonal, diagonal, superdiagonal, superdiagonal2, ipiv, info)
            if (info /= 0) stop "LU factorization failed."

            ! Calculate determinant using LU decomposition
            determinant = 1.0
            do i = 1, P-1
                  determinant = determinant * diagonal(i)
            end do
            
            ! Calculate the inverse using LU factors
            inv  = 0.0
            do j = 1,P-1 
                  inv(j,:) = 0.0
                  inv(j,j) = 1.0      
                  call dgttrs("N", P-1, 1, subdiagonal, diagonal, superdiagonal, superdiagonal2,  ipiv, inv(j,:), P-1, info)
                  if (info /= 0) stop "Solution of linear system failed."
            end do    
      end subroutine Inv_Det 

      ! get_HTR returns the Heat Transfer Rate for each CV as an array with dim=num_coll_var.
      ! the HTR follows from the work performed along the given CV
      function get_HTR(x,x_ext,v_ext)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P, k_ext
            use TAMD_extension, only : num_coll_var, coll_var
            implicit none
            real(real64), dimension(num_coll_var) :: x_ext, v_ext, get_HTR
            real(real64) :: x(P+1)
            get_HTR = k_ext * (coll_var(x) - x_ext) * v_ext 
      end function get_HTR

      ! tot_energy(x,v,q) returns the total energy of the system. 
      function tot_energy(x,v,q)
            use, intrinsic :: iso_fortran_env
            use simulation_parameters, only : P, time
            use masses, only : mass_staging_prime
            implicit none
            real(real64) :: x(P+1), v(P+1), q(P+1)
            real(real64) :: tot_energy, kinetic
            integer :: i
            kinetic = 0
            do i=1,P+1
                  kinetic = kinetic + v(i) * v(i) * mass_staging_prime(i)
            end do
            kinetic = kinetic  * 0.5

            tot_energy = kinetic + pot_energy(x,q)
      end function tot_energy
     
      function tot_energy_ext(x, x_ext, v_ext)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P, k_ext, mass_ext
            use TAMD_extension, only : num_coll_var, coll_var
            use masses, only : mass_staging_prime
            implicit none
            real(real64), dimension(num_coll_var):: x_ext, v_ext, cv
            real(real64) :: tot_energy_ext, kinetic, pot, x(P+1)
            integer :: i
            cv = coll_var(x)
            kinetic = 0
            pot = 0
            do i=1,num_coll_var
                  kinetic = kinetic + v_ext(i) * v_ext(i)*mass_ext(i) 
                  pot = pot + 0.5 * k_ext(i) * (cv(i) - x_ext(i))*(cv(i) - x_ext(i))
            end do
            kinetic = kinetic  * 0.5 

            tot_energy_ext = kinetic + pot
            !tot_energy_ext=pot
      end function tot_energy_ext
       
      function force_terms(x,x_ext)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : P, beta
            use TAMD_extension, only : num_coll_var, k_ext, coll_var
            implicit none
            real(real64) :: coll_var_term_1, coll_var_term_2
            real(real64), dimension(P+1) :: x
            real(real64), dimension(num_coll_var) :: x_ext, couple_term
            real(real64), dimension(P+1) :: U_term, Vtilde_term
            real(real64), dimension(12) :: force_terms
            real(real64), dimension(num_coll_var) :: collective_variables
            integer :: b
            
            !!!!!!!!! U POTENTIAL TERM !!!!!!!!!!!!! 
            U_term(1) = - 0.5 * potential_df_1(x(1)) / P
            U_term(P+1) = - 0.5 * potential_df_1(x(P+1)) / P
            ! REST OF THE COORDINATES (NOT ASSOCIATED WITH COLLECTIVE VARIABLES)
            do b=2,P
                  U_term(b) = - potential_df_1(x(b)) / P
            end do
            
            !!!!!!!!! COUPLED TERMS !!!!!!!!!!!!
            collective_variables = coll_var(x)
            ! ADDITIONAL TERMS DUE TO EXTENDED VARIABLES COUPLED TO COLLECTIVE VARIABLES
            coll_var_term_1 = collective_variables(1) - x_ext(1)
            coll_var_term_2 = collective_variables(2) - x_ext(2)
            
            couple_term(1) = - k_ext(1) * ((coll_var_term_1 / 2) + coll_var_term_2)
            couple_term(2) = - k_ext(2) * ((coll_var_term_1 / 2) - coll_var_term_2)
            
            
            Vtilde_term = -get_V_forces(x)/(2*beta)
            force_terms = [U_term(1), U_term(int((P+1)/4)), U_term(int((P+1)/2)), &
                               U_term(3*int((P+1)/4)), U_term(P+1), &
                         Vtilde_term(1), Vtilde_term((int(P+1)/4)), Vtilde_term(int((P+1)/2)), &
                              Vtilde_term(3*int(((P+1)/4))), Vtilde_term(P+1), couple_term(1), &
                              couple_term(2)]
      end function force_terms

      function energy_terms(x,q,v,x_ext,v_ext)
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : mass, hbar
            use TAMD_extension, only : num_coll_var, get_energy_extended
            use simulation_parameters, only : P, tau_sq, time, beta, wp, mass_ext
            use masses, only : mass_staging_prime
            implicit none
            real(real64) :: x(P+1), q(P+1),v(P+1), x_ext(num_coll_var), v_ext(num_coll_var)
            real(real64) :: sp_potential, spring_potential, energy_terms(6), kinetic
            integer :: i
            ! GET SINGLE PARTICLE POTENTIAL (sp_potential)
            sp_potential = 0.5 * (potential(x(1)) + potential(x(P+1)))
            do i=2,P
                  sp_potential = sp_potential + potential(x(i))
            end do
            sp_potential = sp_potential/P

            ! GET SPRING TERM
            spring_potential = 0

            spring_potential = mass * wp * wp * q(P+1) * q(P+1) / (2 * P)
            do i=2,P
                  spring_potential = spring_potential + mass_staging_prime(i) * wp * wp * q(i) * q(i) /2
            end do
            
            kinetic = 0
            do i=1,num_coll_var
                  kinetic = kinetic + v_ext(i) * v_ext(i) * mass_ext(i) 
            end do
            kinetic = kinetic  * 0.5 
            kinetic =0
            do i=1,P+1
                  kinetic = kinetic + v(i) * v(i) * mass_staging_prime(i)*0.5
            end do

            ! GET extended coupling
            energy_terms = [sp_potential, spring_potential, get_V_potential(x), get_energy_extended(x,x_ext),kinetic, 0.d0]
            energy_terms(6) = sum(energy_terms)

      end function energy_terms

      
      function pot_energy(x,q)
            use, intrinsic :: iso_fortran_env
            use Physical_constants, only : mass
            use simulation_parameters, only : P, wp
            use masses, only : mass_staging_prime
            implicit none
            real(real64) :: x(P+1), q(P+1)
            real(real64) :: sp_potential, spring_potential, pot_energy, V_tilde
            integer :: i
            ! GET SINGLE PARTICLE POTENTIAL (sp_potential)
            sp_potential = 0.5 * (potential(x(1)) + potential(x(P+1)))
            do i=2,P
                  sp_potential = sp_potential + potential(x(i))
            end do
            sp_potential = sp_potential/P

            ! GET SPRING TERM
            spring_potential = 0

            spring_potential = mass * wp * wp * q(P+1) * q(P+1) / (2 * P)
            do i=2,P
                  spring_potential = spring_potential + mass_staging_prime(i) * wp * wp * q(i) * q(i) /2
            end do
            
            ! GET Vtilde
            V_tilde = get_V_potential(x)

            pot_energy = sp_potential + spring_potential + V_tilde
            !pot_energy = sp_potential + V_tilde
      end function pot_energy

end module Forces_Potentials
 

module MD_routines     
      use Physical_constants, only :  mass, w
      use Simulation_parameters, only : P, wp, delta_t, beta
      use Random_Variables, only : sample_Gaussian, Gaussian_noise
      use TAMD_extension, only : num_coll_var
      use masses
      implicit none
contains
      ! get_total_force(x,q,x_ext) returns the total force acting on the staging variables, dim = P+1
      ! for that it first gets forces on the open chain cordinates, transforms them (chain rule)
      ! and adds the interbeads spring terms from the path integral.
      function get_total_force(x,q,x_ext)
            use, intrinsic :: iso_fortran_env
            use Forces_Potentials
            use Staging_transformation, only : force_trafo
            implicit none
            real(real64) ::x(P+1), q(P+1)
            real(real64) :: x_ext(num_coll_var)
            real(real64) :: get_total_force(P+1)
            real(real64) :: F(P+1), F_num(P+1)
            real(real64) :: F_q(P+1)
            integer :: b
            
            ! Forces on the beads coodinates (including coupling to extended sytems)
            F = get_forces(x,x_ext)
            F_q = force_trafo(F)
    
            ! Add interbead potential
            get_total_force(1) = F_q(1) 
            get_total_force(P+1) = F_q(P+1) - mass*wp*wp*q(P+1)/P 
            do b=2,P
                  get_total_force(b) = F_q(b) - mass_staging(b)*wp*wp*q(b)
            end do  
      end function get_total_force
      
      ! The subroutine GGMT_step performs one update of the GGMt variables
      ! and scaling of the velocities of the variables it is coupled to (here the extended system)  
      ! It performs the updates for both thermostats. n1(2), n2(2) contain the GGMT variables,
      ! pn1(2), pn2(2) their conj. momenta. del_t is the timestep for propagation of the dynamics.
      subroutine GGMT_step(v_ext, n1, n2, pn1, pn2, del_t)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : beta_ext, mass_ext
            use TAMD_extension, only : num_coll_var, coll_var, get_forces_extended
            use GGMT_extension, only : Q_
            implicit none
            real(real64), intent(in) :: del_t
            real(real64), dimension(num_coll_var), intent(inout) :: n1, n2, pn1, pn2
            real(real64), dimension(num_coll_var), intent(inout) :: v_ext
            real(real64), dimension(num_coll_var) :: kin_term, lambda, alpha
            
            ! Translation of thermostat momenta
            kin_term = mass_ext * v_ext * v_ext
            pn1 = pn1 + 0.5 * del_t * (kin_term - 1/beta_ext)
            pn2 = pn2 + 0.5 * del_t * (kin_term*kin_term/3 - 1/(beta_ext*beta_ext))
            
            ! Scaling of extendend velocities
            lambda = pn1/Q_(1) + pn2/(Q_(2)*beta_ext)
            alpha = pn2*mass_ext / (3*Q_(2))
            v_ext = v_ext * exp(- 0.25 * del_t * lambda)
            v_ext = v_ext * sqrt(1/(1 + v_ext*v_ext*alpha*del_t)) 
            v_ext = v_ext * exp(- 0.25 * del_t * lambda)
            
            ! Translation of thermostat positions 
            n2 = n2 + del_t * ((1/beta_ext) + mass_ext*v_ext*v_ext) * pn2 / Q_(2)
            n1 = n1 + del_t * pn1 / Q_(1)    
            
            ! Scaling of extendend velocities
            v_ext = v_ext * exp(- 0.25 * del_t * lambda)
            v_ext = v_ext * sqrt(1/(1 + v_ext*v_ext*alpha*del_t)) 
            v_ext = v_ext * exp(- 0.25 * del_t * lambda)
            
            ! Tranlation of thermostat momenta
            kin_term = mass_ext * v_ext * v_ext
            pn1 = pn1 + 0.5 * del_t * (kin_term - 1/beta_ext)
            pn2 = pn2 + 0.5 * del_t * (kin_term*kin_term/3 - 1/(beta_ext*beta_ext))
      end subroutine GGMT_step
      
      ! Langevin_dynamics is the integration scheme for the system. 
      ! It runs Langevin dynamics in the physical variables in a BAOAB scheme,
      !and a velocity verlet with GGMT in a middle scheme for the extended d-AFED variables. 
      subroutine MD_step(x,v,q,x_ext, v_ext, total_force, total_force_ext,n1, n2, pn1, pn2) 
            use, intrinsic :: iso_fortran_env
            use TAMD_extension, only : num_coll_var, beta_ext, k_ext, coll_var, get_forces_extended 
            use Simulation_parameters, only : mass_ext, gamma_
            use masses
            use GGMT_extension, only : Q_
            use Staging_transformation, only : x_trafo, q_trafo 
            implicit none
            integer :: b, i,j,n_bath
            real(real64) :: exp_term, w(3), delta_t_bath
            real(real64), dimension(P+1) :: F_q
            real(real64), dimension(num_coll_var) :: F_ext
            real(real64), dimension(P+1) :: F
            real(real64), dimension(P+1), intent(inout) :: x,q,v
            real(real64), dimension(P+1), intent(inout) :: total_force
            real(real64), dimension(num_coll_var), intent(inout) :: x_ext, v_ext
            real(real64), dimension(num_coll_var), intent(inout) :: n1, n2, pn1, pn2
            real(real64), dimension(num_coll_var), intent(inout) :: total_force_ext
            
            ! the weights for the 4th order Suzuki-Yoshida scheme for integrating the GGMT propagation
            w = [1.3512071919596578,-1.7024143839193155,1.3512071919596578]

            ! B: 
            ! get intermediate velocities by 1st velocity update
            v = v + delta_t * total_force /(2*mass_staging_prime) 
            v_ext = v_ext + delta_t * total_force_ext / (2 * mass_ext)
            
            ! A:
            ! update positions with a half step
            q = q + delta_t_bath*v/2 !contiguous data transfer
            x_ext = x_ext + delta_t_bath * v_ext/2
                  
            ! O:
            ! Langevin Thermostat
            exp_term = exp(-gamma_*delta_t_bath)
            do b=1,P+1
                  v(b) = exp_term*v(b) &
                  + Gaussian_noise() * sqrt((1 - exp_term * exp_term) / (mass_staging_prime(b) * beta))
            end do
            
            ! GGMT steps
            do i=1,3
                  call GGMT_step(v_ext, n1, n2, pn1, pn2, delta_t_bath*w(i))
            end do
            
            ! A: 
            ! the second position update with thermostatted velocities
            q = q + delta_t_bath*v/2 !OK since it moves contiguous chunks of data
            x_ext = x_ext + delta_t_bath * v_ext/2
            
            ! To recalculate the forces, we update first the primitive variables 
            x = x_trafo(q)

            total_force = get_total_force(x,q,x_ext)
            total_force_ext = get_forces_extended(x,x_ext)
        
            ! B: 
            ! get final velocities by 2nd velocity update
            v = v + delta_t * total_force /(2*mass_staging_prime) 
            v_ext = v_ext + delta_t * total_force_ext / (2 * mass_ext)
      end subroutine MD_step

end module MD_routines

module Trajectories
      use Physical_constants, only :  mass, w
      use Simulation_parameters, only : P, wp, delta_t, beta
      use Random_Variables, only : sample_Gaussian, Gaussian_noise
      use TAMD_extension, only : num_coll_var
      use masses
      implicit none
contains
      subroutine Initialization(x, v, x_ext, v_ext) 
            use, intrinsic :: iso_fortran_env
            use Staging_transformation, only : q_trafo
            use TAMD_extension, only : coll_var
            implicit none
            real(real64) :: x_dist_harm_osc_stdev, x_mu
            real(real64) :: v_dist_harm_osc_stdev 
            real(real64), dimension(P+1), intent(inout) :: x, v
            real(real64), dimension(P+1) :: q
            real(real64), dimension(num_coll_var), intent(inout) :: x_ext, v_ext
            integer :: b 
            ! FOr initialization: We sample bead positions from a harmonic oscillator potential
            ! Boltzmann distribution. The corresponding Gaussian dist. is sampled via simple
            ! transformation from the standard normal distribution (sample_Gaussian)
            
            x_dist_harm_osc_stdev = sqrt(1/(beta*mass*(w**2)))
            v_dist_harm_osc_stdev = sqrt(1/(beta*mass)) 
            x_mu = 0
            
            do b = 1, P+1
                  x(b) = sample_Gaussian(x_mu, x_dist_harm_osc_stdev)
                  v(b) = sample_Gaussian(x_mu, v_dist_harm_osc_stdev)
            end do
            q = q_trafo(x)
            x_ext = coll_var(x) 
            v_ext = [0.0d0, 0.0d0]
      end subroutine Initialization 

      ! The subroutine Read_in will initialize the trajectories for the MD simulation.
      ! It will check if a trj file exists. If yes it will read the last entry and continue the previous run
      ! If no, the subroutine Initialization is called.
      subroutine Read_in(x, v, x_ext, v_ext)
            use, intrinsic :: iso_fortran_env
            use Simulation_parameters, only : OUT_FILE, P
            use TAMD_extension, only : num_coll_var
            implicit none
            logical :: file_exists
            real(real64), dimension(2*P + 2 + 2*num_coll_var) :: temp_array
            integer :: iounit, num_elements, i, status
            real(real64), dimension(P+1), intent(inout) :: x, v
            real(real64), dimension(num_coll_var), intent(inout) :: x_ext, v_ext

            INQUIRE(FILE=OUT_FILE, EXIST=file_exists)
            if (file_exists) then
                  ! Open the file
                  open(unit=iounit, file=OUT_FILE, status='old', action='read')

                  ! Read to the end of the file to get the last line
                  
                  do
                        read(iounit, *, IOSTAT=status) temp_array
                        if (status /= 0) exit
                  end do
                  
                  ! Save the entries to thw trajectories of position and vellocities of the open chain
                  ! and the extended d-AFED system
                  x = temp_array(1:P+1)
                  x_ext = temp_array(P+2:P+1+num_coll_var)
                  v = temp_array(P+2+num_coll_var:2*P+2+num_coll_var)
                  v_ext = temp_array(2*P+3+num_coll_var:2*P+2+2*num_coll_var)
                  close(iounit)
            else
                  ! File doesn't exist, call Initialization subroutine
                  call Initialization(x, v, x_ext, v_ext)
                   
                  open(iounit, file=OUT_FILE, status='NEW')
                  close(iounit)
            end if

      end subroutine Read_in
end module Trajectories 


program main
      use, intrinsic :: iso_fortran_env
      use Physical_constants, only : mass, w
      use Simulation_parameters, only : time, P, numsteps, eqsteps, skipsteps, OUT_FILE, ENE_FILE,&
      FOR_FILE, assign_data, beta, mass_ext
      use Forces_Potentials, only : tot_energy, pot_energy, energy_terms, force_terms, get_HTR
      use masses
      use TAMD_extension
      use Staging_transformation
      use Forces_Potentials
      use GGMT_extension
      use Trajectories
      use MD_routines
      implicit none
      ! type declaration statements
      real(real64), dimension(:), allocatable :: x, q, v, force, force_ext   
      real(real64), dimension(num_coll_var) :: n1, n2, pn1, pn2
      real(real64), dimension(num_coll_var) :: x_ext, v_ext, collective_var, HTR
      integer :: t, i, l
      real(real64), dimension(:,:), allocatable :: output, force_output, energy_output
      real(real64) :: kin_E, pot_E, tot_E_ext, tot_E_GGMT
      integer :: values(8)
      
      print *,'start job at'
      call date_and_time(values=values)
      
      ! read in config file with given t and P for this simulation. Compute staging masses and select potential:
      call assign_data
      call get_mass
      call select_potential

      ! now knowing the size, we can allocate memory for the trajectory arrays
      allocate(x(P+1), q(P+1), v(P+1), force(P+1), force_ext(num_coll_var),&
      output(numsteps,2*P+2+2*num_coll_var), force_output(numsteps,P+1+num_coll_var),&
      energy_output(numsteps,6))
      x = 0
      v = 0
      call get_trj(x, v, x_ext, v_ext)

      force = get_total_force(x,q,x_ext)
      force_ext = get_forces_extended(x,x_ext)
      n1 = [0,0]
      n2 = [0,0]
      pn1 = [1,1]
      pn2 = [1,1]
      call get_thermostat_masses
      
      ! MD simulation
      
      ! EQUILIBRATION
      do t = 1, eqsteps
            call MD_step(x, v, q, x_ext, v_ext, force, force_ext, n1, n2, pn1, pn2) 
      end do

      ! DATA COLLECTION
      do t = 1, numsteps
            do l=1,skipsteps
                  call MD_step(x, v, q, x_ext, v_ext, force, force_ext, n1, n2, pn1, pn2)
            end do
            HTR = get_HTR(x,x_ext,v_ext)
            tot_E_GGMT = GGMT_energy(n1,n2,pn1,pn2)
            tot_E_ext = tot_energy_ext(x,x_ext,v_ext) 
            pot_E = tot_energy(x,v,q)
            kin_E = sum(mass_staging_prime*v*v)*0.5
            output(t,1:P+1) = x
            output(t,P+2:P+1+num_coll_var) = x_ext
            output(t,P+2+num_coll_var:2*P+2+num_coll_var) = v
            output(t,2*P+3+num_coll_var:2*P+2+2*num_coll_var) = v_ext

            force_output(t,:P+1) = force
            force_output(t,P+2:) = force_ext
            energy_output(t,:) = [HTR(1),HTR(2), kin_E, pot_E ,tot_E_ext ,tot_E_GGMT]
      end do
      open(1, file=OUT_FILE, status='old')
            do i=1,size(output(:,1))
                  write(1,*) output(i,:)
            end do
      close(1)
      open(2, file=ENE_FILE, status='NEW')
            do i=1,size(energy_output(:,1))
                  write(2,*) energy_output(i,:)
            end do
      close(2)
      open(3, file=FOR_FILE, status='NEW')
            do i=1,size(force_output(:,1))
                  write(3,*) force_output(i,:)
            end do
      close(3)
      print *,'finished at'
      call date_and_time(values=values)
end program main
