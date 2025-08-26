      ! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib
      use gyre_lib
      
      implicit none
      
      real(dp), save :: initial_mass_change, initial_star_age, acc_duration
      real(dp), save ::  Lacc, Ladd
      character(LEN=strlen) :: logs_dir
      integer, save :: steps_interval
      ! these routines are called by the standard run_star check_model
      contains


      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         character(len = 19) :: name_base
         character(len = 5) :: mass_name
         character(len = 24) :: final_name
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  
         s% other_energy => energy_routine

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
         s% other_adjust_mdot => other_adjust_mdot



         s% accretion_h2 = s% x_ctrl(1)*1d-6
         s% accretion_he3 = s% x_ctrl(2)*1d-6
         s% accretion_he4 = 0.224 + 2*s% initial_z - s% x_ctrl(2)*1d-6
         s% accretion_h1 = 1 - s% accretion_he3 - s% accretion_he3 -s% accretion_he4 - s % initial_z

         s% xa_central_lower_limit(1) = s% accretion_h1  - 0.001

         name_base = "../../models/model_"
         write(mass_name, '(f5.3)') s% initial_mass
         final_name = name_base // mass_name
         s% job% saved_model_name = final_name
         initial_mass_change = s% mass_change
         s% job% years_for_initial_dt = 2d-3/initial_mass_change/500
         write(*,*) s% job% years_for_initial_dt
         acc_duration = 3d-3/initial_mass_change
         s% max_timestep = acc_duration/500*secyer
         s% job% years_for_initial_dt = s% max_timestep/secyer
      end subroutine extras_controls
      
      

      subroutine energy_routine(id, ierr)
        use const_def
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        type (star_info), pointer :: s

        real(dp) :: mass, mdot, mdot_msun, alpha, radius, mcumul

        real(dp) :: xposk, dif, xpos, thingie, thingie2, Mstar
        integer :: k, xposk_orig
        ierr = 0
        call star_ptr(id, s, ierr)

                 

        alpha = 0.1

        mass = s% mstar
        mdot = s% mstar_dot
        if (s% mstar_dot < 1d-10) then
          s% newton_itermin_until_reduce_min_corr_coeff = 15
          s% corr_coeff_limit = 0.05
          s% use_gold_tolerances = .true.
        end if
        if (mdot == 0.) return

        radius = s% r(1)

        mcumul = 0.
        k = 0
        Lacc = 0.5 * (standard_cgrav * (s% mstar_dot) * (mass)) / radius

        Ladd = alpha/2*standard_cgrav * (s% star_mass *Msun) * (s% mstar_dot) / (s% r(1)) ![erg/s]
        Lacc =  safe_log10((1-alpha)/2*standard_cgrav * (s% star_mass *Msun) *&
                (s% mstar_dot) / (s% r(1))  /Lsun) ![erg/s]
        Mstar=(s% star_mass *Msun) ![g]


       
        xpos = s% x_ctrl(5)
        do k = 1, s% nz
          if (s% star_mdot > 0.0d0 .and. s% m(k)>=(1.0d0-xpos)*(s% star_mass)*msol ) then
            s% extra_heat(k) = 2.0d0*Ladd/(Mstar*xpos**2)*(s% m(k)/Mstar-(1.0d0-xpos))
          end if
        end do

       end subroutine energy_routine


       subroutine other_adjust_mdot(id, ierr)
          use star_def
          implicit none
          integer, intent(in) :: id
          integer, intent(out) :: ierr
          type(star_info), pointer :: s
          ierr = 0
          call star_ptr(id, s, ierr) ! retrieve star id
          if (s% star_age< initial_star_age + acc_duration) then 
            s% mstar_dot = initial_mass_change*(1-(s%star_age-initial_star_age)/acc_duration)*(1-(s%star_age-initial_star_age)/acc_duration)*Msun/secyer
            steps_interval = 20 *s% x_integer_ctrl(1)
          else
            s% max_timestep = 1d9*secyer
            s% mstar_dot = 0
            steps_interval = s% x_integer_ctrl(1)
          end if

        end subroutine other_adjust_mdot

      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         real(dp) :: mass_dif
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         initial_star_age = s% star_age


         ! Initialize GYRE

         call gyre_init('gyre.in')

         ! Set constants

         call gyre_set_constant('G_GRAVITY', standard_cgrav)
         call gyre_set_constant('C_LIGHT', clight)
         call gyre_set_constant('A_RADIATION', crad)

         call gyre_set_constant('M_SUN', msol)
         call gyre_set_constant('R_SUN', rsol)
         call gyre_set_constant('L_SUN', lsol)
         logs_dir = s% log_directory
      end subroutine extras_startup
      
!      example lines from inlist

!      x_integer_ctrl(1) = 10 ! output GYRE info at this step interval
!      x_logical_ctrl(1) = .false. ! save GYRE info whenever save profile

!      x_integer_ctrl(2) = 2 ! max number of modes to output per call
!      x_logical_ctrl(2) = .false. ! output eigenfunction files

!      x_integer_ctrl(3) = 1 ! mode l (e.g. 0 for p modes, 1 for g modes)
         ! should match gyre.in mode l
!      x_integer_ctrl(4) = 1 ! order
!      x_ctrl(1) = 0.6d-4 ! freq < this (Hz)
!      x_ctrl(2) = 0.3d4 ! growth < this (days)

  integer function gyre_in_mesa_extras_finish_step(id)
    integer, intent(in) :: id

    integer                   :: ierr
    logical                   :: call_gyre
    type (star_info), pointer :: s
    real(dp), allocatable     :: global_data(:)
    real(dp), allocatable     :: point_data(:,:)
    integer                   :: ipar(5), step_interval, mode_l
    real(dp)                  :: rpar(2), growth_lim
    character(LEN=strlen)     :: summary_filename		!from me
    
    ierr = 0
    call star_ptr(id, s, ierr)
    if (ierr /= 0) return
    
    gyre_in_mesa_extras_finish_step = terminate

    ! Extract parameters

    step_interval = steps_interval!s% x_integer_ctrl(1)
    mode_l = s% x_integer_ctrl(3)
    growth_lim = s% x_ctrl(5)

    if (s% need_to_save_profiles_now .and. s% x_logical_ctrl(1)) then
       call_gyre = .TRUE.
    elseif (step_interval > 0) then
       call_gyre = MOD(s% model_number, step_interval) == 0
    else
       call_gyre = .FALSE.
    endif

    ! If necessary, call GYRE

    if (call_gyre) then

       call star_get_pulse_data(s%id, 'GYRE', .FALSE., .FALSE., .FALSE., global_data, point_data, ierr)
       if (ierr /= 0) then
          print *,'Failed when calling star_get_pulse_data'
          return
       end if

       call gyre_set_model(global_data, point_data, 101)



       write(summary_filename, 120) TRIM(logs_dir) , '/summary.', s% model_number, '.dat'    	!from me
120      format(A, A,I0,A,I0,A)									!from me
       open(UNIT=10, FILE=summary_filename, STATUS='REPLACE')					!from me


       if (growth_lim > 0d0) then
          write(*, 100) 'model', 'order', 'freq (Hz)', &
             'P (sec)', 'P (min)', 'P (day)', 'growth (day)', '(4pi*im/re)'
100       format(2A8,6A20)


          write(10, 200) 'model_number', 'n_pg', 'n_g', 'n_p','j', 'freq(Hz)', 'freq(day)',  &	!from me
             'P(sec)', 'P(min)', 'P(day)', 'growth(day)', '(4pi*im/re)'	, 'W', 'eta'		!from me
200       format(5A8,9A20)									!from me


       else
          write(*, 100) 'model', 'order', 'freq (Hz)', &
             'P (sec)', 'P (min)', 'P (day)', '(4pi*im/re)'
101       format(2A8,5A20)
       end if

       ipar(1) = s% model_number
       ipar(2) = s% x_integer_ctrl(4) ! order_target
       rpar(1) = s% x_ctrl(4)         ! freq_target
       rpar(2) = growth_lim
       ! write_flag
       if (s% x_logical_ctrl(2)) then
          ipar(3) = 1
       else
          ipar(3) = 0
       endif
       ipar(4) = s% x_integer_ctrl(2) ! max_to_write
       ipar(5) = 0 ! num_written

       call gyre_get_modes(mode_l, process_mode_, ipar, rpar)

       close(10)										!from me
    end if

    gyre_in_mesa_extras_finish_step = keep_going

  contains

    subroutine process_mode_ (md, ipar, rpar, retcode)

      type(mode_t), intent(in) :: md
      integer, intent(inout)   :: ipar(:)
      real(dp), intent(inout)  :: rpar(:)
      integer, intent(out)     :: retcode

      character(LEN=strlen) :: filename
      integer               :: unit, k, model_number, num_written, max_to_write !, order_target
      complex(dp)           :: cfreq
      real(dp)              :: freq, growth, freq_target, growth_lim
      logical               :: write_flag, freq_okay
      
      max_to_write = ipar(4)
      num_written = ipar(5)
      if (num_written >= max_to_write) return
      ipar(5) = num_written + 1
      
      model_number = ipar(1)
      !order_target = ipar(2)
      freq_target = rpar(1)
      growth_lim = rpar(2)
      write_flag = (ipar(3) == 1)

      cfreq = md% freq('HZ')
      write(*,*) md% eta(), md% W()
      freq = REAL(cfreq)
      freq_okay = (abs(freq - freq_target) < freq_target*3d-2)

      if (growth_lim > 0d0) then ! report growth
         if (AIMAG(cfreq) > 0._dp) then ! unstable
            growth = 1d0/(2*pi*24*3600*AIMAG(cfreq))
            write(*, 100) ipar(1), md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), growth, &
               4*pi*AIMAG(cfreq)/freq, md% W(), md% eta()
100         format(2I8,E20.4,3F20.4,4E20.4)


            if (freq_okay .and. growth < growth_lim) & ! .and. md%n_pg == order_target) &
               write(*,*) 'matched target frequency'
         else ! stable
            write(*, 110) ipar(1), md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), 'stable'
110         format(2I8,E20.4,3F20.4,A20)
         end if


       growth = 1d0/(2*pi*24*3600*AIMAG(cfreq))							!from me
       write(10, 200) ipar(1), md%n_pg, md%n_g, md%n_p, md%j, freq, freq*24*3600, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600), growth, &
            4*pi*AIMAG(cfreq)/freq, md% W(), md% eta()						!from me
200         format(5I8,2E20.8,3F20.8,4E20.8)							!from me


      else ! growth_lim <= 0 means ignore it
         write(*, 111) ipar(1), md%n_pg, freq, 1d0/freq, 1d0/(freq*60), 1d0/(freq*24*3600)
111      format(2I8,E20.4,3F20.4)
         if (freq_okay) & ! .and. md%n_pg == order_target) &
            write(*,*) 'matched target frequency'
      endif

      if (write_flag) then

      ! Write the mode radial & horizontal eigenfunctions, together with the differential work
       
         write(filename, 120) TRIM(logs_dir) , '/eigfunc.', model_number, '.', md%n_pg,'.', md%j, '.dat'
120      format(A, A,I0,A,I0,A,I0,A)

         print *,'Writing eigenfunction to file:', TRIM(filename)
         write(*,*)

         open(NEWUNIT=unit, FILE=filename, STATUS='REPLACE')

         write(unit, 130) 'x=r/R', 'Real(xi_r/R)', 'Imag(xi_r/R)', 'Real(xi_h/R)', 'Imag(xi_h/R)', 'dW/dx', 'dE/dx'
130      format(7(1X,A24))

         do k = 1, md%n_k
            write(unit, 140) md%gr%pt(k)%x, md%xi_r(k), md%xi_h(k), md%dW_dx(k), md%dE_dx(k)
140         format(7(1X,E24.16))
         end do

         close(unit)

      end if

      retcode = 0

    end subroutine process_mode_

  end function gyre_in_mesa_extras_finish_step
         integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         real(dp) :: mass_dif 
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0

  end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
        integer function extras_check_model(id)
        integer, intent(in) :: id
        integer :: ierr
        character(len = 13) :: name_base
        character(len = 5) :: mass_name
        character(len = 18) :: final_name
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
  

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 0
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.
         

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 3
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         names(1) = 'log_opacity'
         names(2) = 'grad_temperature'
         names(3) = 'grad_L'
         do k = 1, s% nz
            vals(k,1) = safe_log10(s% opacity(k))
            vals(k,2) = safe_log10(s% grad_temperature(k))
            vals(k,3) = safe_log10(s% gradL(k))
         end do
 
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         extras_finish_step = gyre_in_mesa_extras_finish_step(id)

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! Finalize GYRE

         call gyre_final()
      end subroutine extras_after_evolve


      
      end module run_star_extras
      
