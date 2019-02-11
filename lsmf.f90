        subroutine get_lsmf_snap
        use fit_lammps_class 
        use lapack_diag_simm
        use lapack_inverse
        use common_var
        use LAMMPS
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        implicit none
        integer                         :: aaa,bbb,l,i,k,n,v,s,t,ll,tt
        integer                         :: ss,vv,inf,LWORK
        double precision, allocatable   :: x(:),B(:),A(:,:),work(:),AA(:,:)
        double precision, allocatable   :: W(:),Tval(:,:),YY(:,:),Y(:,:),A2(:,:)
        double precision, allocatable   :: ACM(:,:),BCM(:),BB(:)
        double precision                :: dump,ave

        real (C_double), pointer   :: ener(:,:) => null()
        real (C_double), pointer   :: f(:,:) => null()
        real (c_double), pointer   :: kind_nat(:) => null()
        real (c_double), pointer   :: Eref => null()
        real (C_double), pointer   :: fx_ref(:) => null()
        real (C_double), pointer   :: fy_ref(:) => null()
        real (C_double), pointer   :: fz_ref(:) => null()
        real (C_double), pointer   :: id_dbl(:)=> null()
        integer, allocatable       :: map(:),id(:)

        integer                         :: nkinds,bi_order,npar,npar2fit,tot_frames
        integer                         :: nats2fit,quadflag
        double precision                :: gen_cutoff
        character (len=10000)           :: snap_string,snap_string2
        double precision, allocatable   :: cutoff(:),radii(:),sigma(:)
        integer, allocatable            :: kind_count(:),type2fit(:)
        character (len=2),allocatable   :: label(:)


         tot_frames=sys%tot_frames         
         npar2fit=sys%npar2fit
      
         if(print_bi) open(121,file='Bi_compoenents.dat')

         open(16,file=sys%inp_fit)
         read(16,*) gen_cutoff,bi_order,npar,quadflag
         read(16,*) nkinds
         allocate(label(nkinds))
         allocate(type2fit(nkinds))
         allocate(cutoff(nkinds))
         allocate(sigma(nkinds))
         allocate(radii(nkinds))
         allocate(kind_count(nkinds))
         kind_count=0
         do i=1,nkinds
          read(16,*) label(i),type2fit(i),radii(i),cutoff(i),sigma(i)
          write(*,*) label(i),type2fit(i),radii(i),cutoff(i),sigma(i)
         enddo         
         close(16)
         write(snap_string,*) gen_cutoff,'1.0000',bi_order,(radii(i),i=1,nkinds),(cutoff(i),i=1,nkinds),&
                 'quadraticflag ',quadflag
         write(snap_string2,*) (type2fit(i),i=1,nkinds)
         write(*,*) trim(snap_string)
         write(*,*) 'Fitting types: ',trim(snap_string2)

         if(.not.skip_fit)then

          open(222,file='snapcoeff')

          write(222,*) nkinds,npar
          do i =1,nkinds
           write(222,*) label(i),radii(i),cutoff(i)
           do n=1,npar
            write(222,*) 0.0000000
           enddo
          enddo
          flush(222)

          close(222)

         endif

         
         allocate(lmp(sys%ndata))

         do i=1,sys%ndata
           
           call lammps_open_no_mpi ('lmp -log none', lmp(i))
           call lammps_file (lmp(i),sys%inp)
           call lammps_command (lmp(i),'read_data '//trim(sys%data(i)%inp_data))
           call lammps_command (lmp(i),'group fitsnap type '//trim(snap_string2))
           call lammps_file (lmp(i),sys%inp_fix)
           call lammps_command (lmp(i), &
                  'compute sna_e all sna/atom '//trim(snap_string)//&
                  ' diagonal 3 rmin0 0 switchflag 1')
           call lammps_command (lmp(i),&
                'compute type all property/atom type')
           call lammps_command (lmp(i),&
                'compute id all property/atom id')
           call lammps_command (lmp(i), 'compute pe_ener all pe') 

           if(fit_forces)then
            call lammps_command (lmp(i), &
                  'compute sna_f all snad/atom '//trim(snap_string)//&
                  ' diagonal 3 rmin0 0 switchflag 1')
            call lammps_command (lmp(i), &
                  'compute f_x all property/atom fx')
            call lammps_command (lmp(i), &
                  'compute f_y all property/atom fy')
            call lammps_command (lmp(i), &
                  'compute f_z all property/atom fz')
           endif

         enddo


         !!! do kernel

         if(refine)then

          call kernel%dealloc()

          kernel%nkinds=nkinds
          allocate(kernel%K(nkinds))
          do i=1,nkinds
           kernel%K(i)%sigma=sigma(i)
           kernel%K(i)%nenvs=0
          enddo

          do i=1,sys%ndata
           call lammps_command (lmp(i), 'run 0')
           call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
           call lammps_extract_compute (ener, lmp(i), 'sna_e', 1, 2)
           kind_count=0
           do t=1,sys%data(i)%nats
            kind_count(nint(kind_nat(t)))=kind_count(nint(kind_nat(t)))+1
           enddo
           do t=1,nkinds
            kernel%K(t)%nenvs=kernel%K(t)%nenvs+kind_count(t)*sys%data(i)%frames
           enddo
          enddo
          do t=1,nkinds
           allocate(kernel%K(t)%B(kernel%K(t)%nenvs,size(ener,1)))
          enddo

          do i=1,nkinds
           kernel%K(i)%nenvs=0
          enddo

          do i=1,sys%ndata
           do l=1,sys%data(i)%frames

            call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
            call lammps_command (lmp(i), 'run 0')
            call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
            call lammps_extract_compute (ener, lmp(i), 'sna_e', 1, 2)
            call lammps_extract_compute (Eref, lmp(i), 'pe_ener', 0, 0) 
            if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
            if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
            call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)
          
            id=INT(id_dbl)
            id_dbl=>null()

            do k=1,sys%data(i)%nats          
             map(id(k))=k
            enddo

            do t=1,sys%data(i)%nats
             v=kernel%K(nint(kind_nat(map(t))))%nenvs+1
             do k=1,size(ener,1)
              kernel%K(nint(kind_nat(map(t))))%B(v,k)=ener(k,map(t))
             enddo
             kernel%K(nint(kind_nat(map(t))))%nenvs=&
                 kernel%K(nint(kind_nat(map(t))))%nenvs+1
            enddo

            id_dbl=>null()
            ener=>null()
            Eref=>null()
            kind_nat=>null()

           enddo

           deallocate(id)
           deallocate(map)

          enddo

         endif ! refine
        
         !!! do fitting

        if(.not.skip_fit)then

         v=1
         vv=1

         do i=1,sys%ndata

          do l=1,sys%data(i)%frames

           call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
           call lammps_command (lmp(i), 'run 0')
           call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
           call lammps_extract_compute (ener, lmp(i), 'sna_e', 1, 2)
           call lammps_extract_compute (Eref, lmp(i), 'pe_ener', 0, 0) 
           if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
           if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
           call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)
          
           id=INT(id_dbl)
           id_dbl=>null()

           do k=1,sys%data(i)%nats          
            map(id(k))=k
           enddo

           if(l.eq.1)then

            npar=size(ener,1)+1
            kind_count=0
            
            do t=1,sys%data(i)%nats
             kind_count(nint(kind_nat(t)))=kind_count(nint(kind_nat(t)))+1
            enddo

            if(.not.cs)then
             if(.not.allocated(B)) allocate(B(npar2fit+nkinds-1))
!             if(.not.allocated(B)) allocate(B(npar2fit+nkinds))
             if(.not.allocated(A)) then
              allocate(A(npar2fit+nkinds-1,npar*nkinds))
!              allocate(A(npar2fit+nkinds,npar*nkinds))
              A=0.0d0
             endif
            else    ! compressive sensing
             if(.not.allocated(B)) allocate(B(npar2fit+nkinds-1+(npar-1)*nkinds))
!             if(.not.allocated(B)) allocate(B(npar2fit+nkinds+(npar-1)*nkinds))
             if(.not.allocated(A)) then
              allocate(A(npar2fit+(npar-1)*nkinds+nkinds-1,npar*nkinds))
!              allocate(A(npar2fit+(npar-1)*nkinds+nkinds,npar*nkinds))
              A=0.0d0
             endif
            endif

            if(pca)then
             open(1313,file='PCA.dat')
             if(.not.allocated(A2)) allocate(A2(tot_frames,(npar-1)*nkinds))
             if(.not.allocated(BB)) allocate(BB(tot_frames))
             if(.not.allocated(W)) allocate(W((npar-1)*nkinds))
             if(.not.allocated(Y)) allocate(Y(tot_frames,(npar-1)*nkinds))
             if(.not.allocated(YY)) allocate(YY((npar-1)*nkinds,(npar-1)*nkinds))
             if(.not.allocated(Tval)) allocate(Tval(tot_frames,(npar-1)*nkinds))
            endif

           endif
       
           if(fit_ener)then 

            s=1
            do k=1,nkinds
             A(v,s)=kind_count(k)
             s=s+npar
            enddo
            
            do k=2,npar
             do t=1,sys%data(i)%nats

              if(print_bi) write(121,*) ener(k-1,map(t))

               s=k+(nint(kind_nat(t))-1)*npar
               A(v,s)=A(v,s)+ener(k-1,map(t))*sys%data(i)%weight
               if(pca)then
                ss=k-1+(nint(kind_nat(t))-1)*(npar-1)
                A2(vv,ss)=A(vv,ss)+ener(k-1,map(t))
               endif

             enddo
            enddo

            B(v)=(sys%data(i)%ener(l)-Eref)*sys%data(i)%weight
            v=v+1
            vv=vv+1
            
            ener=>null()
            Eref=>null()

           endif ! fitener
        
           if(fit_forces )then

            call lammps_extract_compute (f, lmp(i), 'sna_f', 1, 2)
            call lammps_extract_compute (fx_ref, lmp(i), 'f_x', 1, 1)
            call lammps_extract_compute (fy_ref, lmp(i), 'f_y', 1, 1)
            call lammps_extract_compute (fz_ref, lmp(i), 'f_z', 1, 1)

            do t=1,sys%data(i)%nats

             s=1
             do k=1,nkinds
              A(v,s)=0.0
              s=s+npar
             enddo
             
             do n=1,nkinds
              s=(n-1)*(3*(npar-1))+1  ! npar-1?
              do k=2,npar
               A(v,((n-1)*npar+k))=f(s,map(t))*sys%data(i)%weight
               s=s+1
              enddo
             enddo
             B(v)=(sys%data(i)%fx(l,t)-fx_ref(map(t)))*sys%data(i)%weight
             v=v+1

             s=1
             do k=1,nkinds
              A(v,s)=0.0
              s=s+npar
             enddo

             do n=1,nkinds
              s=(n-1)*(3*(npar-1))+1+npar-1
              do k=2,npar
               A(v,((n-1)*npar+k))=f(s,map(t))*sys%data(i)%weight
               s=s+1
              enddo
             enddo
             B(v)=(sys%data(i)%fy(l,t)-fy_ref(map(t)))*sys%data(i)%weight
             v=v+1

             s=1
             do k=1,nkinds
              A(v,s)=0.0
              s=s+npar
             enddo

             do n=1,nkinds
             s=(n-1)*(3*(npar-1))+1+2*(npar-1)
              do k=2,npar
               A(v,((n-1)*npar+k))=f(s,map(t))*sys%data(i)%weight
               s=s+1
              enddo
             enddo
             B(v)=(sys%data(i)%fz(l,t)-fz_ref(map(t)))*sys%data(i)%weight
             v=v+1

            enddo

            f=>null()
            fx_ref=>null()
            fy_ref=>null()
            fz_ref=>null()


           endif ! end if on forse

          enddo      ! ciclo on frames

          if(allocated(map)) deallocate(map)
          if(allocated(id)) deallocate(id)
    
         enddo   ! ciclo su data

         s=1
         do k=1,nkinds-1
!         do k=1,nkinds!-1
          B(v)=0.0d0
          A(v,s)=1.0d0
          s=s+npar
          v=v+1
         enddo

! compressive sensing

         if(cs)then
          do k=1,nkinds
           do l=2,npar
            B(v)=0.0d0
            s=(k-1)*npar+l
            A(v,s)=cm_val
            v=v+1
           enddo
          enddo
         endif


! principal component analysis

        if(pca)then

         Y=A2
         do l=1,size(Y,2)
          ave=0.0d0
          do s=1,size(Y,1)
           ave=ave+Y(s,l)
          enddo
          Y(1:size(Y,1),l)=Y(1:size(Y,1),l)-ave/size(Y,1)
         enddo
         YY=matmul(transpose(Y),Y)
         call new_diag(size(YY,1),YY,W)

         if(.not.allocated(AA)) allocate(AA(tot_frames,70))

         AA=matmul(A2,YY(:,1:70))
       
         W=sqrt(abs(W))
         W=W/sum(W)
         write(1313,*) 'Principal Components Values'
         do k=1,size(W)
          write(1313,*) W(k)
         enddo

!!!
!         aaa=70
!         bbb=tot_frames
!         lwork=bbb+64*bbb+1000
!         allocate(work(lwork))
!         call dgels('N',bbb,aaa,1,AA,bbb,BB,bbb,WORK,LWORK,inf)
!         deallocate(work)
!         if(inf.ne.0)then
!          write(*,*) 'zgels failed',inf
!           stop
!         else

!         YY=transpose(YY)
!         BB=matmul(YY(:,1:70),BB(1:70))

!         open(222,file='snapcoeff')

!         l=1
!         write(222,*) nkinds,npar
!         do i =1,nkinds
!          write(222,*) label(i),radii(i),cutoff(i)
!          write(222,*) 0.000000000
!          do n=1,npar-1
!           write(222,*) BB(l)
!           l=l+1
!          enddo
!         enddo
!         close(222)

!        endif
!!!

        endif
        
        aaa=npar*nkinds
        if(cs)then
!         bbb=npar2fit+(nkinds*npar)!-1
         bbb=npar2fit+(nkinds*npar)-1
        else
!         bbb=npar2fit+nkinds!-1
         bbb=npar2fit+nkinds-1
        endif
        lwork=bbb+64*bbb+1000
        allocate(work(lwork))
        call dgels('N',bbb,aaa,1,A,bbb,B,bbb,WORK,LWORK,inf)
        deallocate(work)
        if(inf.ne.0)then
         write(*,*) 'zgels failed',inf
         stop
        else
        

        open(222,file='snapcoeff')

        l=1
        write(222,*) nkinds,npar
        do i =1,nkinds
         write(222,*) label(i),radii(i),cutoff(i)
         do n=1,npar
          write(222,*) B(l)
          l=l+1
         enddo
        enddo

!        dump=0.0d0       
!        do i=npar*nkinds+1,size(B)
!         dump=dump+B(i)**2
!        enddo

!        write(222,*) 'dgels residulal: ',dump 
        close(222)

        open(333,file='snapparam')        
        write(333,*) 'rcutfac ',gen_cutoff
        write(333,*) 'twojmax ',bi_order
        write(333,*) 'quadraticflag ',quadflag
        write(333,*) 'rfac0 1.00000'
        write(333,*) 'rmin0 0'
        write(333,*) 'diagonalstyle 3'
        write(333,*) 'switchflag 1'

        close(333)

        do i=1,sys%ndata
         call lammps_file (lmp(i),sys%inp_fix)          
        enddo

        endif

        endif ! skip_fit 
             
        return
        end subroutine get_lsmf_snap
