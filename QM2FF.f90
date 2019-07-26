        program fitsnap
        use fit_lammps_class
        use common_var
        implicit none
        integer                        :: l
        character (len=100)            :: command,input,datas,output,new_datas
       
         if(iargc().eq.0)then
          write(*,*) 'FitSnap Usage:'
          write(*,*) '-datas : list of files to be fitted'
          write(*,*) '-inp : Lammps like input to set put calculation'
          write(*,*) '-pot_fit : Lammps like input for the potentials to be &
                       fitted '
          write(*,*) '-pot_fix : Lammps like input to be appended, not &
                       touched by the optimization'
          write(*,*) '-ener   : switch on the fitting of energies'
          write(*,*) '-forces : switch on the fitting of forces'
          write(*,*) '-tensor : switch on the fitting of a tensor'
          write(*,*) '-compress <val> : activates ridge-regression'
          write(*,*) '-pca            : activate principal components &
                                        analysys.'
          write(*,*) '-print_bi       : print bispectrum components'
          write(*,*) '-out'
         stop
         endif   

         do l=1,iargc()

          call getarg(l,command)

          if(trim(command).eq.'-datas')then        
           call getarg(l+1,command)
           read(command,*) datas
          endif
          if(trim(command).eq.'-inp')then        
           call getarg(l+1,command)
           read(command,*) sys%inp
          endif
          if(trim(command).eq.'-pot_fit')then        
           call getarg(l+1,command)
           read(command,*) sys%inp_fit
          endif
          if(trim(command).eq.'-pot_fix')then        
           call getarg(l+1,command)
           read(command,*) sys%inp_fix
          endif
          if(trim(command).eq.'-forces')then        
           fit_forces=.true.
          endif
          if(trim(command).eq.'-ener')then        
           fit_ener=.true.
          endif
          if(trim(command).eq.'-skip_fit')then        
           skip_fit=.true.
          endif
          if(trim(command).eq.'-print_bi')then        
           print_bi=.true.
          endif
          if(trim(command).eq.'-pca')then        
           pca=.true.
          endif
          if(trim(command).eq.'-refine')then        
           refine=.true.
           call getarg(l+1,command)
           read(command,*) refine_maxiter
           call getarg(l+2,command)
           read(command,*) refine_temp
           call getarg(l+3,command)
           read(command,*) thr_kernel
          endif
          if(trim(command).eq.'-compress')then        
           cs=.true.
           call getarg(l+1,command)
           read(command,*) cm_val
          endif
          if(trim(command).eq.'-E0')then        
            !!!! -E0 0 -E0 1:
           call getarg(l+1,command)
           read(command,*) e0cs
          endif

         enddo

         iter=0
         new_datas=datas

         do while ( (iter.le.refine_maxiter .and. refine) .or. iter.eq.0 )

          if(iter.gt.0)then

           call execute_command_line('tail -n $( head -n 1 '&
                //trim(datas)//' ) '//trim(datas)//' > new_datas')
           call execute_command_line('sed -i "1i $(( 1+$( head -n 1 '&
                //trim(datas)//') ))" new_datas'  )

           new_datas='new_datas'
           open(9,file=new_datas,access='append')

           if(fit_ener .and. (.not. fit_forces))then
            write(9,*) trim(sys%data(1)%inp_data)//' ',&
                        iter,'new_geo.xyz',' new_geo.ener 1.0'
            close(9)
           endif
           if(fit_ener .and. fit_forces)then
            write(9,*) trim(sys%data(1)%inp_data)//' ',&
                        iter,'new_geo.xyz',&
                        ' new_geo.ener new_geo.force 1.0'
            close(9)
           endif
           if((.not. fit_ener) .and. fit_forces)then
            write(9,*) trim(sys%data(1)%inp_data)//' ',&
                        iter,'new_geo.xyz',&
                        ' new_geo.force 1.0'
            close(9)
           endif

          endif

          call sys%read_sys(new_datas,fit_ener,fit_forces)
          call get_lsmf_snap
          call get_chi2
          if(refine) call refine_snap
          iter=iter+1

         enddo

        return
        end program fitsnap


        subroutine get_chi2 
        use fit_lammps_class
        use common_var
        use LAMMPS
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        implicit none
        double precision     :: chi_val_ener,chi_val_fx,chi_val_fy,chi_val_fz
        integer              :: l,i,k,n,v
        real (C_double), pointer :: ener => null()
        real (C_double), pointer :: fx(:) => null()
        real (C_double), pointer :: fy(:) => null()
        real (C_double), pointer :: fz(:) => null()
        real (C_double), pointer   :: id_dbl(:)=> null()
        integer, allocatable       :: map(:),id(:)
       
        chi_val_ener=0.0d0
        chi_val_fx=0.0d0
        chi_val_fy=0.0d0
        chi_val_fz=0.0d0

        open(222,file='energy_rms.dat')
        open(333,file='force_rms.dat')
        write(222,*) 'RMS Energies'
        write(333,*) 'RMS Forces'

!!      calcola chi2
        do i=1,sys%ndata

         do l=1,sys%data(i)%frames

           call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
           call lammps_command (lmp(i), 'run 0')

           if(fit_ener)then

            call lammps_extract_compute (ener, lmp(i), 'pe_ener', 0, 0)

            chi_val_ener=chi_val_ener+( (ener-sys%data(i)%ener(l))/sys%data(i)%nats )**2

!            write(222,*) i,l,ener/sys%data(i)%nats,sys%data(i)%ener(l)/sys%data(i)%nats,&
!                 (ener-sys%data(i)%ener(l))/sys%data(i)%nats
            write(222,*) i,l,ener,sys%data(i)%ener(l),&
                 (ener-sys%data(i)%ener(l))

            ener=>null()

           endif

           if(fit_forces)then

            call lammps_extract_compute (fx, lmp(i), 'f_x', 1, 1)
            call lammps_extract_compute (fy, lmp(i), 'f_y', 1, 1)
            call lammps_extract_compute (fz, lmp(i), 'f_z', 1, 1)
            if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
            if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
            call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)
          
            id=INT(id_dbl)
            id_dbl=>null()

            do k=1,sys%data(i)%nats          
             map(id(k))=k
            enddo

            do k=1,sys%data(i)%nats
             chi_val_fx=chi_val_fx+( fx(map(k))-sys%data(i)%fx(l,k) )**2
             chi_val_fy=chi_val_fy+( fy(map(k))-sys%data(i)%fy(l,k) )**2
             chi_val_fz=chi_val_fz+( fz(map(k))-sys%data(i)%fz(l,k) )**2
             write(333,*) k,fx(map(k)),fy(map(k)),fz(map(k)),sys%data(i)%fx(l,k),sys%data(i)%fy(l,k),sys%data(i)%fz(l,k)
            enddo

            fx=>null()
            fy=>null()
            fz=>null()

           endif
                  
         enddo      ! ciclo su frames

        call lammps_close (lmp(i))
        if(allocated(map)) deallocate(map)
        if(allocated(id)) deallocate(id)

        enddo   ! ciclo su data

        if(allocated(lmp)) deallocate(lmp)

        write(222,*) 'Total RMS Energies (Kcal/mol/atom):',sqrt(chi_val_ener/sys%tot_frames)
        write(333,*) 'Total RMS Forces (Kcal/mol/Ang): ',sqrt(chi_val_fx/sys%tot_frames),&
                                          sqrt(chi_val_fy/sys%tot_frames),&
                                          sqrt(chi_val_fz/sys%tot_frames)
        close(222)
        close(333)

        return
        end subroutine get_chi2


