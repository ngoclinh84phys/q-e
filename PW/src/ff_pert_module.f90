MODULE ff_module
      !
      USE kinds,                      ONLY : dp
      !
      SAVE
      !
      CHARACTER(LEN=256) :: pert_pot_file
      LOGICAL            :: finite_field_pert = .FALSE.
      COMPLEX(DP), ALLOCATABLE :: dpotff_r (:,:)
      REAL(DP), ALLOCATABLE :: rhoin_r (:,:)
      REAL(DP) :: ff_scale_factor=0.001
      !
ENDMODULE ff_module

SUBROUTINE init_ff_calc()
      !
      USE kinds,                ONLY : DP
      USE io_global,            ONLY : stdout
      USE control_flags,        ONLY : restart, gamma_only
      USE fft_base,             ONLY : dffts, dfftp
      USE fft_interfaces,       ONLY : fwfft, invfft
      USE gvecs,                ONLY : doublegrid
      USE scf,                  ONLY : rho, rho_core, rhog_core, &
                                       v, vltot, vrs, kedtau, vnew
      USE lsda_mod,             ONLY : nspin
      USE gvect,                ONLY : ngm
      USE io_base,              ONLY : read_rhog
      USE ff_module,            ONLY : dpotff_r, rhoin_r, ff_scale_factor, &
                                       pert_pot_file, finite_field_pert
      USE input_parameters,     ONLY : ffield_response_calc, input_ffpot_file 
      USE fft_rho,              ONLY : rho_g2r, rho_r2g
      !
      ! local variable
      !
      INTEGER :: is
      CHARACTER :: slabel
      CHARACTER(LEN=6) :: ilabel, jlabel
      CHARACTER(LEN=:), ALLOCATABLE :: fname
      COMPLEX(DP), ALLOCATABLE :: pot_tmp(:,:)!, aux(:), aux_g(:)
      REAL(DP), ALLOCATABLE :: aux_r(:,:)
      COMPLEX(DP), ALLOCATABLE :: aux_g(:,:)
      ! 
      IF (.not. restart) CALL errore( 'init_ff_calc', 'Finite field calculation need a restart calculation ' ,1 )
      IF (.not. gamma_only) CALL errore( 'init_ff_calc', 'Finite field calculation need a gamma trick method ' ,1 )
      !
      ! Set global var
      finite_field_pert = ffield_response_calc
      pert_pot_file = input_ffpot_file
      !
      ! Set global var
      ALLOCATE (dpotff_r(dffts%nnr, nspin))
      !
      ! ... Read input pertubing potential in G space in file
      !
      WRITE(stdout,*) "... Reading potential from file: ", pert_pot_file  
      !
      ALLOCATE (pot_tmp(ngm, nspin))
      pot_tmp(:,:) = (0.0_DP, 0.0_DP)
      CALL read_rho_xml_ff(pot_tmp, pert_pot_file)
      !
      ! ... Bring to R space
      ALLOCATE (aux_r(dffts%nnr, nspin))
      ALLOCATE (aux_g(ngm, nspin))
      aux_r(:,:) = (0.0_DP, 0.0_DP)
      aux_g(:,:) = (0.0_DP, 0.0_DP)
      aux_g(:,:) = pot_tmp(:,:) 
      CALL rho_g2r (dfftp, aux_g(:,:), aux_r(:,:))
      DO is = 1, nspin
         dpotff_r(:, is) = aux_r(:,is)
      ENDDO
      DEALLOCATE (aux_g)
      DEALLOCATE (aux_r)
      DEALLOCATE (pot_tmp)
      !
      WRITE(stdout,*) "... Done reading potential..."
      !
      ! ... Allocate density rhoin_r to save gs density
      !
      ! Set global var
      ALLOCATE (rhoin_r(dffts%nnr, nspin))
      !
      rhoin_r(:,:) = 0.0_DP
      rhoin_r(:,:) = rho%of_r(:,:)
      !
      ! Add perturbing potential to vHxc local, and reset vrs
      !
      CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
      !
      WRITE(stdout,*) SIZE(vrs)
      DO ir = 1, 50
      WRITE(stdout,*) ir, vrs(ir,1), dpotff_r(ir,1)
      ENDDO
      !
ENDSUBROUTINE

SUBROUTINE ends_ff_calc()
      !
      USE kinds,              ONLY : DP
      USE io_global,          ONLY : stdout
      USE control_flags,      ONLY : restart, gamma_only
      USE cell_base,          ONLY : omega
      USE scf,                ONLY : rho
      USE fft_base,           ONLY : dffts, dfftp
      USE fft_interfaces,     ONLY : fwfft,invfft
      USE lsda_mod,           ONLY : nspin
      USE gvect,              ONLY : ngm
      USE mp,                 ONLY : mp_sum
      USE mp_bands,           ONLY : intra_bgrp_comm
      USE ff_module,          ONLY : dpotff_r, rhoin_r, ff_scale_factor, &
                                     pert_pot_file
      !
      ! local variable
      !
      INTEGER  :: is
      REAL(DP) :: rPi
      CHARACTER :: slabel
      CHARACTER(LEN=6) :: ilabel, jlabel
      CHARACTER(LEN=:), ALLOCATABLE :: fname
      REAL(DP), ALLOCATABLE :: drho_r(:,:)
      COMPLEX(DP), ALLOCATABLE :: pot_tmp(:,:), aux(:)
      !
      ! calculate variational density
      !
      ALLOCATE (drho_r(dffts%nnr, nspin))
      !
      drho_r(:,:) = 0.0_DP
      !
      IF ( nspin==2 ) THEN
         !
         rhoin_r(:,:) = rho%of_r(:,:) - rhoin_r(:,:)
         drho_r(:,1)  = (rhoin_r(:,1) + rhoin_r(:,2))*0.5_DP
         drho_r(:,2)  = (rhoin_r(:,1) - rhoin_r(:,2))*0.5_DP
         !
      ELSE
         !
         drho_r(:,1) = rho%of_r(:,1) - rhoin_r(:,1)
         !
      ENDIF
      !
      drho_r(:,:) = drho_r(:,:) / ff_scale_factor
      !
      rPi = SUM(DBLE(dpotff_r(:, 1:nspin)) * drho_r(:, 1:nspin))
      rPi = omega * rPi / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
      CALL mp_sum( rPi, intra_bgrp_comm )
      WRITE(stdout,*) "Relaxed Koopmans rPi=", rPi
      !
      ALLOCATE (pot_tmp(ngm, nspin))
      ALLOCATE (aux(dffts%nnr))
      DO is = 1, nspin
         aux(:) = CMPLX(drho_r(:,is), KIND=DP)
         CALL fwfft('Rho', aux, dfftp)
         pot_tmp(:,is) = aux(dfftp%nl(:))
      ENDDO
      DEALLOCATE (aux)
      !
      ! Write G-space response density to hdf5  format
      !
      fname = 'Resdens'//pert_pot_file
      CALL write_rho_xml_ff(pot_tmp, fname)
      !
      DEALLOCATE (pot_tmp)
      DEALLOCATE (drho_r)
      ! Reset global var
      DEALLOCATE (dpotff_r)
      DEALLOCATE (rhoin_r)
      !
ENDSUBROUTINE
!
SUBROUTINE read_rho_xml_ff(dpot, filename)
      !
      USE kinds,                ONLY : DP
      USE lsda_mod,             ONLY : lsda, nspin
      USE fft_base,             ONLY : dfftp
      USE gvect,                ONLY : ig_l2g
      USE mp_pools,             ONLY : my_pool_id
      USE mp_bands,             ONLY : my_bgrp_id, root_bgrp_id, root_bgrp, &
                                       intra_bgrp_comm, inter_bgrp_comm, nbgrp
      USE mp_images,            ONLY : intra_image_comm
      USE mp,                   ONLY : mp_bcast
      USE io_base,              ONLY : read_rhog, write_rhog
      !
      CHARACTER(LEN=*), INTENT(IN) :: filename
      COMPLEX(DP), INTENT(INOUT) :: dpot(dfftp%ngm, nspin)
      !
      ! ... Read input pertubing potential in G space
      !
      IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
      CALL read_rhog(filename, root_bgrp, intra_bgrp_comm, &
                     ig_l2g, nspin, dpot(:,1:nspin) )
      IF( nbgrp > 1 ) CALL mp_bcast( dpot, root_bgrp_id, inter_bgrp_comm )
      !
ENDSUBROUTINE
!
SUBROUTINE write_rho_xml_ff(drho, filename)
      !
      USE kinds,                ONLY : DP
      USE lsda_mod,             ONLY : lsda, nspin
      USE fft_base,             ONLY : dfftp
      USE control_flags,        ONLY : gamma_only
      USE gvect,                ONLY : ig_l2g, mill
      USE mp_pools,             ONLY : my_pool_id
      USE mp_bands,             ONLY : my_bgrp_id, root_bgrp_id, root_bgrp, &
                                       intra_bgrp_comm, inter_bgrp_comm, nbgrp
      USE mp_images,            ONLY : intra_image_comm
      USE io_base,              ONLY : read_rhog, write_rhog
      !
      CHARACTER(LEN=*), INTENT(IN) :: filename
      COMPLEX(DP), INTENT(IN) :: drho(dfftp%ngm, nspin)
      REAL(DP) :: dum(3) = [0_dp, 0_dp, 0_dp]
      !
      IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
      CALL write_rhog(filename, &
                 root_bgrp, intra_bgrp_comm, dum, dum, dum, &
                 gamma_only, mill, ig_l2g, drho(:,1:nspin) )
      !
ENDSUBROUTINE
