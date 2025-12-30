      program dynamics
!Allocating and freeing 2D arrays
      implicit none
      integer :: input_file, Natoms,i,j
      double precision, allocatable :: coord(:,:)
      double precision, allocatable :: mass(:)
      double precision, allocatable :: distance(:,:)
      integer :: i_stat
      input_file = 9
      open(input_file, file='inp.txt')
      open(10,file='dynamics.xyz')
      Natoms = read_Natoms(input_file)
      write(6,*) Natoms
      call read_molecule(input_file, Natoms, coord, mass)
      DO i=1,Natoms
        write(6, *) (coord(i,j), j=1,3), mass(i)
      END DO

      call compute_distances(Natoms, coord, distance)


contains
      integer function read_Natoms(input_file) result(atoms)
      implicit none
      integer, intent(in) :: input_file
      read(input_file, *) atoms
      end function read_Natoms

      subroutine read_molecule(input_file, Natoms, coord, mass)
      implicit none
      integer, intent(in) :: input_file
      integer, intent(in) :: Natoms
      double precision, allocatable, intent(out) :: coord(:,:)
      double precision, allocatable, intent(out) :: mass(:)
      integer :: i_stat,i,j
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
      end subroutine read_molecule

!subroutine that calculates the distances between each possible pair 
      subroutine compute_distances(Natoms, coord, distance)
      implicit none
      integer :: i,j,k
      double precision :: diff
      integer, intent(in) :: Natoms
      double precision, intent(in) :: coord(Natoms,3)
      double precision, allocatable, intent(out) :: distance(:,:)
      allocate(distance(Natoms,Natoms), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for distance!"
                stop
      end if
      DO j=1,Natoms
        DO i=1,Natoms
          if(i<j) then
            distance(j,i)=distance(i,j)
          end if
          diff=0.d0
          DO k=1,3
            diff=diff+(coord(j,k)-coord(i,k))**2
          end do
          distance(j,i)=sqrt(diff)
        end do
      end do
      write(6,*) 
      DO i=1,Natoms
        write(6, *) (distance(i,j), j=1,Natoms)
      END DO
      end subroutine compute_distances

end program dynamics
