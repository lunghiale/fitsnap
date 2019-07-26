        module common_var
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use fit_lammps_class

         TYPE(system)                   :: sys
         type (C_ptr), allocatable      :: lmp(:)
         type(kernel_global)            :: kernel
         logical                        :: fit_forces=.false.,pca=.false.,cs=.false.,fit_ener=.false.
         logical                        :: fit_tensors=.true.
         logical                        :: skip_fit=.false.,print_bi=.false.,refine=.false.,metrix=.true.
         double precision               :: cm_val=1.0d0,refine_temp,thr_kernel=0.5d0
         integer                        :: refine_maxiter,iter,tens_order,e0cs=0

                                                
        end module common_var

