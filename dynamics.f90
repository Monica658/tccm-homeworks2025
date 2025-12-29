      program dynamics
!Allocating and freeing 2D arrays
 
      double precision, allocatable :: a(:,:)
      integer :: i_stat

      allocate(a(m,n), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed!"
                stop
        end if

        deallocate(a)

      subroutine read_molecule(inp.txt, Natoms, coord, mass)
      integer, intent(in) :: inp.txt
      integer, intent(in) :: Natoms
      double precision, intent(out) :: coord(Natoms,3)
      double precision, intent(out) :: mass(Natoms)
      !reads the atomic coordinates and masses and stores the data in
      !the arrays coord and mass provided as parameters. coord is a
      !two-dimensinal array distance of size (Natoms*Natoms)
      end subroutine read_molecule

      subroutine compute_distances(Natoms, coord, distance)
      implicit none 
      integer, intent(in) :: Natoms
      double precision, intent(in) :: coord(Natoms,3)
      double precision, intent(out) :: coord(Natoms,Natoms)
      end subroutine compute_distances

