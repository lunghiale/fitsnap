        subroutine refine_snap
        use fit_lammps_class 
        use lapack_diag_simm
        use lapack_inverse
        use random_numbers_class
        use common_var
        use LAMMPS
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        implicit none
        integer                         :: l,i,j,v,k
        double precision, allocatable   :: r(:)

        character (len=10000)           :: snap_string
        double precision, allocatable   :: cutoff(:),radii(:)
        integer, allocatable            :: kind_count(:),type2fit(:)
        character (len=2),allocatable   :: label(:)
        integer                         :: nkinds,bi_order,npar,npar2fit,tot_frames
        integer                         :: nats2fit,quadflag
        double precision                :: gen_cutoff

        real (C_double), pointer        :: id_dbl(:)=> null()
        real (c_double), pointer        :: kind_nat(:) => null()
        integer, allocatable            :: map(:),id(:)

        double precision, allocatable   :: max_kernel(:)
        double precision                :: dft_ener,val,rand_seed

        real (C_double), pointer   :: B(:,:) => null()
        real (c_double), pointer   :: ff_ener => null()

         allocate(lmp(1))

         call lammps_open_no_mpi ('lmp -log none', lmp(1))
         call lammps_file (lmp(1),sys%inp)
         call lammps_command (lmp(1),'read_data '//trim(sys%data(1)%inp_data))
         call lammps_file (lmp(1),sys%inp_fix)
         call lammps_command (lmp(1), 'compute pe_ener all pe') 

         open(16,file=sys%inp_fit)
         read(16,*) gen_cutoff,bi_order,npar,quadflag
         read(16,*) nkinds
         allocate(label(nkinds))
         allocate(type2fit(nkinds))
         allocate(cutoff(nkinds))
         allocate(radii(nkinds))
         allocate(kind_count(nkinds))
         kind_count=0
         do i=1,nkinds
          read(16,*) label(i),type2fit(i),radii(i),cutoff(i)
          write(*,*) label(i),type2fit(i),radii(i),cutoff(i)
         enddo         
         close(16)
         write(snap_string,*) gen_cutoff,'1.0000',bi_order,(radii(i),i=1,nkinds),(cutoff(i),i=1,nkinds),&
                 'quadraticflag ',quadflag

         call lammps_command (lmp(1), &
                  'compute sna_e all sna/atom '//trim(snap_string)//&
                  ' diagonal 3 rmin0 0 switchflag 1')
         call lammps_command (lmp(1),&
                'compute type all property/atom type')
         call lammps_command (lmp(1),&
                'compute id all property/atom id')

        ! minimize energy

         call lammps_command (lmp(1),'thermo 1')
         call lammps_command (lmp(1),&
                'thermo_style custom step time temp pe etotal ') 
         call lammps_command (lmp(1), &
                'minimize 1.0e-8 1.0e-8 1000 100000')


         write(snap_string,*) (label(i)//' ',i=1,size(label))

         call lammps_command (lmp(1),&
                'dump xyz_dump all xyz 5 geo_opt.xyz')

         call lammps_command (lmp(1),&
                'dump_modify            xyz_dump element '//trim(snap_string))

        ! set velocities

         write(snap_string,*) refine_temp
         
         call lammps_command (lmp(1),'timestep 0.5')
         call lammps_command (lmp(1),'variable t equal'//trim(snap_string))

         call random_number(rand_seed)
         write(snap_string,*) nint(rand_seed*10000.0d0)

         call lammps_command (lmp(1), &
                'velocity all create $t '//trim(snap_string)//&
                ' dist gaussian')
         call lammps_command (lmp(1),'velocity all zero linear')
         call lammps_command (lmp(1),'velocity all zero angular')
         
         call lammps_command (lmp(1), &
                'fix 1 all nvt temp $t $t 100.0 tchain 3')

         open(13,file='new_geo.xyz',access='append')
         open(14,file='kernel_max.dat',access='append')
         open(15,file='new_geo.ener',access='append') 

         if(.not.allocated(max_kernel)) &
               allocate(max_kernel(sys%data(1)%nats))
!         if(.not.allocated(max_kernel)) &
!                allocate(max_kernel(kernel%nkinds))

         write(14,*) '## nkinds: ',kernel%nkinds,'nenvs: ',(kernel%K(i)%nenvs,i=1,kernel%nkinds)
         flush(14)

         do i=1,20000

          call lammps_command (lmp(1), 'run 5')
          call lammps_extract_compute (B, lmp(1), 'sna_e', 1, 2)
          call lammps_extract_compute (ff_ener, lmp(1), 'pe_ener', 0, 0)
          call lammps_extract_compute (kind_nat, lmp(1), 'type', 1, 1)
          if(.not. allocated(id)) allocate(id(sys%data(1)%nats))
          if(.not. allocated(map)) allocate(map(sys%data(1)%nats))
          call lammps_extract_compute (id_dbl, lmp(1), 'id', 1, 1)   

          id=INT(id_dbl)
          id_dbl=>null()

          do k=1,sys%data(1)%nats          
           map(id(k))=k
          enddo

          max_kernel=0.0d0

          do k=1,size(B,2)
           l=nint(kind_nat(map(k)))
           do v=1,kernel%K(l)%nenvs
            val=0.0d0
            do j=1,size(B,1)
             val=val-(B(j,map(k))-kernel%K(l)%B(v,j))**2  
            enddo
            val=exp(val/2*kernel%K(l)%sigma**2)
            if(val.gt.max_kernel(map(k))) max_kernel(map(k))=val
           enddo
          enddo

          write(14,*) i,max_kernel
          flush(14)

          if(any(max_kernel.lt.thr_kernel))then
           if(allocated(r)) deallocate(r)
           call lammps_gather_atoms (lmp(1),'x', 3, r)
           write(13,*) sys%data(1)%nats
           write(13,*)
           l=1
           do j=1,sys%data(1)%nats
            write(13,*) sys%data(1)%label(j),r(l),r(l+1),r(l+2)
            l=l+3
           enddo
           flush(13)
           call execute_command_line('./run_DFT_scf.x')
           open(16,file='last_dft_ener.dat')
           read(16,*) dft_ener
           close(16)
           write(15,*) dft_ener,ff_ener
           flush(15)
           close(13)
           close(14)
           close(15)
           call lammps_close (lmp(1))
           deallocate(lmp)
           return
           write(*,*) 'PUPPA'
           flush(6)
          endif 

         enddo

         call lammps_close (lmp(1))
         deallocate(lmp)
         refine=.false.
         close(13)
         close(14)
         close(15)

        return   
        end subroutine refine_snap
