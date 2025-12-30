      program dynamics
!Allocating and freeing 2D arrays
      implicit none
      integer :: input_file, Natoms, atoms
      double precision, allocatable :: coord(:,:)
      double precision, allocatable :: mass(:)
      integer :: i_stat
      input_file = 9
      open(input_file, file='inp.txt')
      open(10,file='dynamics.xyz')
      atoms = read_Natoms(input_file,Natoms)
      write(6,*) atoms
      call read_molecule(input_file, Natoms, coord, mass)
!      DO i=1,Natoms
!        write(6, *) (coord(i,j), j=1,3), mass(i)
!      END DO
contains
      integer function read_Natoms(input_file,Natoms) 
      implicit none
      integer, intent(in) :: input_file
      integer :: Natoms
      read(input_file, *) Natoms
      end function read_Natoms

      subroutine read_molecule(input_file, Natoms, coord, mass)
      implicit none
      integer, intent(in) :: input_file
      integer, intent(in) :: Natoms
      double precision, allocatable, intent(out) :: coord(:,:)
      double precision, allocatable, intent(out) :: mass(:)
      integer :: i_stat
      allocate(coord(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for coord!"
                stop
        end if
      allocate(mass(Natoms), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for mass!"
                stop
        end if
      DO i=1,Natoms
        read(input_file, *) (coord(i,j), j=1,3), mass(i)
      END DO 
      DO i=1,Natoms
        write(6, *) (coord(i,j), j=1,3), mass(i)
      END DO
      end subroutine read_molecule
end program dynamics
