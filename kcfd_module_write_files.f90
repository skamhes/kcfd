module module_write_files

        implicit none
        
        public :: write_tecplot_file_b
        public :: write_tecplot_file_v
    contains
        subroutine write_tecplot_file_b
            use module_common_data, only : nnodes      , & !# of nodes
                                    x, y, z     , & !nodal coords
                                    ntria, tria , & !# of triangles and triangle list
                                    nquad, quad , & !# of quads and quad list
                                    p2, zero, bc_type, nb
            use module_common_data, only : filename_tecplot_b
            !To access the solution data.
            use module_ccfv_data_grid, only : cell, ncells, bound, bound_export
            use module_ccfv_data_soln, only : w, ir, iu, iv, iw, ip, Temp, mu

            use module_input_parameter       , only : project_name, navier_stokes
            use module_ccfv_gradient         , only : my_alloc_int_ptr
            
            implicit none 
            type bnode_type
                integer                               :: nbnodes 
                integer, dimension(:),  pointer       :: bnodes  
            end type bnode_type
            integer :: i, os, ibn
            
            integer                           :: j, k, ib, bcell_i, candidate_node, nk, j_count
            real(p2), dimension(:,:), pointer :: wn
            integer , dimension(:  ), pointer :: nc
            real(p2), dimension(:  ), allocatable :: tn, mun, Mn
            type(bnode_type), dimension(nb)   :: bnode_data
            logical                           :: already_added

            allocate(wn(5, nnodes))
            allocate(nc(   nnodes))
            if (navier_stokes) allocate(tn(   nnodes))
            if (navier_stokes) allocate(mun(  nnodes))
            allocate(mn(   nnodes))
            
            nc = 0
            wn = zero
            if (navier_stokes) tn = zero
            if (navier_stokes) mun = zero

            bound_loop : do ib = 1,nb
                bface_loop : do i = 1,bound(ib)%nbfaces
                    bcell_i = bound(ib)%bcell(i)
                    do k = 1,cell(bcell_i)%nvtx
                        wn(:,cell(bcell_i)%vtx(k)) = wn(:,cell(bcell_i)%vtx(k))  + w(:,bcell_i)
                        if (navier_stokes) then
                            tn(  cell(bcell_i)%vtx(k)) = tn(  cell(bcell_i)%vtx(k))  + temp(bcell_i)
                            mun(  cell(bcell_i)%vtx(k)) = mun(  cell(bcell_i)%vtx(k))  + mu(bcell_i)
                        end if
                        nc(cell(bcell_i)%vtx(k))   = nc(cell(bcell_i)%vtx(k)) + 1
                    end do
                end do bface_loop
            end do bound_loop
            do j = 1,nnodes
                wn(:,j) = wn(:,j) / nc(j) ! copmute an average
                if (navier_stokes) then
                    tn(j) = tn(j) / nc(j)
                    mun(j) = mun(j) / nc(j)
                end if
                Mn(j) = sqrt(wn(2,j)**2 + wn(3,j)**2 + wn(4,j)**2) ! mach number
            end do

            allocate(bound_export(nb))
            boundary_loop : do ib = 1,nb
                bnode_data(ib)%nbnodes = zero
                allocate(bnode_data(ib)%bnodes(1))
                bnode_data(ib)%bnodes(1) = zero
                bound_export(ib)%nbfaces = bound(ib)%nbfaces
                allocate( bound_export(ib)%bfaces( 5,bound(ib)%nbfaces ) )  
                bound_export(ib)%bfaces = zero
                bfaces_loop : do i = 1,bound(ib)%nbfaces
                    bound_export(ib)%bfaces(1,i) = bound(ib)%bfaces(1,i)
                    bface_vertex_loop : do k = 2,(bound(ib)%bfaces(1,i) + 1) ! loop through number of vertices for face
                        if (bnode_data(ib)%nbnodes == 0 ) then
                            bnode_data(ib)%nbnodes = bnode_data(ib)%nbnodes + 1
                            bnode_data(ib)%bnodes  = bound(ib)%bfaces(k,i)
                            bound_export(ib)%bfaces(k,i) = 1
                            cycle bface_vertex_loop 
                        end if 
                        candidate_node = bound(ib)%bfaces(k,i)
                        already_added = .false.
                        bnodes_loop : do nk = 1,bnode_data(ib)%nbnodes
                            if (candidate_node == bnode_data(ib)%bnodes(nk)) then
                                already_added = .true.
                                bound_export(ib)%bfaces(k,i) = nk
                                exit bnodes_loop  
                            end if 
                        end do bnodes_loop 
                        if (.not.already_added) then
                            bnode_data(ib)%nbnodes = bnode_data(ib)%nbnodes + 1
                            call my_alloc_int_ptr(bnode_data(ib)%bnodes,bnode_data(ib)%nbnodes)
                            bnode_data(ib)%bnodes(bnode_data(ib)%nbnodes) = candidate_node
                            bound_export(ib)%bfaces(k,i) = bnode_data(ib)%nbnodes
                        end if
                    end do bface_vertex_loop 
                end do bfaces_loop
            end do boundary_loop

            write(*,*)
            write(*,*) "-------------------------------------------------------"
            write(*,*) ' Writing Tecplot file = ', trim(filename_tecplot_b)
            write(*,*)
        
            !Open the output file.
            open(unit=8, file=filename_tecplot_b, status="unknown", iostat=os)   

            !---------------------------------------------------------------------------

            !(0)Header information

            write(8,*) 'TITLE = "GRID"'
            if (navier_stokes) then
                write(8,*) 'VARIABLES = "x","y","z","rho","u","v","w","p","T","mu","M"'
            else
                write(8,*) 'VARIABLES = "x","y","z","rho","u","v","w","p","M"'
            end if
            ! (1) Nodal Information
            do ib = 1,nb
                write(8,*) 'ZONE T = "',ib,'"  n=', bnode_data(ib)%nbnodes, &
                                ',e=', bound(ib)%nbfaces,' , zonetype=fequadrilateral, datapacking=point'
                do j_count = 1,bnode_data(ib)%nbnodes
                    j = bnode_data(ib)%bnodes(j_count)
                    if (navier_stokes) then
                        write(8,'(11es25.15)') x(j), y(j), z(j), wn(ir,j), wn(iu,j), wn(iv,j), wn(iw,j), wn(ip,j),tn(j),mun(j),Mn(j)
                    else 
                        write(8,'(11es25.15)') x(j), y(j), z(j), wn(ir,j), wn(iu,j), wn(iv,j), wn(iw,j), wn(ip,j), Mn(j)
                    end if
                end do
                ! Loop through faces
                do i = 1,bound_export(ib)%nbfaces
                    if (bound_export(ib)%bfaces(1,i) == 3) then ! write tri as a degenerate quad
                        write(8,'(4i10)') bound_export(ib)%bfaces(2,i), bound_export(ib)%bfaces(3,i), & 
                                            bound_export(ib)%bfaces(4,i), bound_export(ib)%bfaces(4,i)
                    else if (bound_export(ib)%bfaces(1,i) == 4) then ! write quad
                        write(8,'(4i10)') bound_export(ib)%bfaces(2,i), bound_export(ib)%bfaces(3,i), & 
                        bound_export(ib)%bfaces(4,i), bound_export(ib)%bfaces(5,i)
                    end if

                end do
            end do
            close(8)
            write(*,*)
            write(*,*) ' End of Writing Tecplot file = ', trim(filename_tecplot_b)
            write(*,*) "-------------------------------------------------------"
            write(*,*)
        end subroutine write_tecplot_file_b

        subroutine write_tecplot_file_v
            use module_common_data, only : nnodes      , & !# of nodes
                                    x, y, z     , & !nodal coords
                                    ntria, tria , & !# of triangles and triangle list
                                    nquad, quad , & !# of quads and quad list
                                    p2, zero
            use module_common_data, only : filename_tecplot_v
            !To access the solution data.
            use module_ccfv_data_grid, only : cell, ncells, bound, bound_export
            use module_ccfv_data_soln, only : w, ir, iu, iv, iw, ip

            use module_input_parameter       , only : project_name
            integer :: i, os, ibn
            integer                           :: j, k, ib, bcell_i, candidate_node, nk, j_count
            real(p2), dimension(:,:), pointer :: wn
            integer , dimension(:  ), pointer :: nc
            logical                           :: already_added

            allocate(wn(5,nnodes))
            allocate(nc(  nnodes))

            nc = 0
            wn = zero

            do i = 1, ncells
                !Loop over vertices of the cell i
                do k = 1, cell(i)%nvtx
                    wn( :,cell(i)%vtx(k) ) = wn( :,cell(i)%vtx(k) ) + w(:,i) !<- Add up solutions
                    nc( cell(i)%vtx(k)  ) = nc( cell(i)%vtx(k)  ) + 1      !<- Count # of contributing cells
                end do
            end do
            do j = 1,nnodes
                wn(:,j) = wn(:,j) / nc(j) ! copmute an average
            end do
            write(*,*)
            write(*,*) "-------------------------------------------------------"
            write(*,*) ' Writing Tecplot file = ', trim(filename_tecplot_v)
            write(*,*)
        
            !Open the output file.
            open(unit=8, file=filename_tecplot_v, status="unknown", iostat=os)   

            !---------------------------------------------------------------------------

            !(0)Header information

            write(8,*) 'TITLE = "GRID"'
            write(8,*) 'VARIABLES = "x","y","z","rho","u","v","w","p"'
            ! (1) Nodal Information
            write(8,*) ' ZONE T="Volume" n=', nnodes,', e=',ncells,', zonetype=FEBRICK, datapacking=point'
            do j = 1,nnodes
                write(8,'(10es25.15)') x(j), y(j), z(j), wn(ir,j), wn(iu,j), wn(iv,j), wn(iw,j), wn(ip,j)
            end do
            ! Write nodes based on tecplot guidlines.  Tets and prisms are written as
            ! degenerate hex cells.  For tets nodes the pattern is 1,2,3,3,4,4,4,4.  
            ! For Prisms: 1,2,3,3,4,5,6,6.
            ! Fer: C:/Users/Karsten%20Hendrickson/Dropbox/Karsten%20Storage/Karsten%20Schoolwork/PhD/I%20Do%20Like%20CFD%20Source%20Code/Reference%20Papers/Tecplot%20360_data_format_guide.pdf
            ! 
            do i=1,ncells
                if ( cell(i)%nvtx == 4 ) then ! tet
                    write(8,'(8i10)') cell(i)%vtx(1),cell(i)%vtx(2),cell(i)%vtx(3),cell(i)%vtx(3), &
                        cell(i)%vtx(4),cell(i)%vtx(4),cell(i)%vtx(4),cell(i)%vtx(4)
                else if ( cell(i)%nvtx == 6 ) then
                    write(8,'(8i10)') cell(i)%vtx(1),cell(i)%vtx(2),cell(i)%vtx(3),cell(i)%vtx(3), &
                        cell(i)%vtx(4),cell(i)%vtx(5),cell(i)%vtx(6),cell(i)%vtx(6)                    
                end if
            end do
            !---------------------------------------------------------------------------

            !Close the output file.
            close(8)

            write(*,*)
            write(*,*) ' End of Writing Tecplot file = ', trim(filename_tecplot_v)
            write(*,*) "-------------------------------------------------------"
            write(*,*)
        end subroutine write_tecplot_file_v

        subroutine write_data_file
            use module_ccfv_data_grid , only : ncells
            use module_common_data, only : filename_data
            use module_ccfv_data_soln, only : u, res_norm, res_norm_initial
            
            implicit none
            integer :: j, os

            !Open the output file.
            open(unit=8, file=filename_data, status="unknown", iostat=os)   

            write(*,*)
            write(*,*) "-------------------------------------------------------"
            write(*,*) ' Writing Data file = ', trim(filename_data)
            write(*,*)

            write(8,'(5es25.15)') res_norm(1),res_norm(2),res_norm(3),res_norm(4),res_norm(5)
            write(8,'(5es25.15)') res_norm_initial(1),res_norm_initial(2), &
                        res_norm_initial(3),res_norm_initial(4),res_norm_initial(5)

            do j = 1,ncells
                write(8,'(5es25.15)')  u(1,j), u(2,j), u(3,j), u(4,j), u(5,j)
            end do
            close(8)

            
        end subroutine write_data_file
end module module_write_files
