      program dynamics
!Allocating and freeing 2D arrays
      implicit none
      integer :: input_file, Natoms,i,j
      double precision, allocatable :: coord(:,:)
      double precision, allocatable :: mass(:)
      double precision, allocatable :: distance(:,:)
      double precision, allocatable :: velocity(:,:)
      double precision :: epsilon, sigma, Vljtot, Ttot, Etot
! In this program will be considered 5 Xenon atoms, the epsilon value is 1.77 kJ/mol
! the sigma value is 4.10 Angstroms
! time must be expressed in picoseconds and distance in nanometers
      input_file = 9
      epsilon = 1.77
      sigma = 4.10
      sigma=sigma*0.1d0
      open(input_file, file='inp.txt')
      open(10,file='dynamics.xyz')
      Natoms = read_Natoms(input_file)
      write(6,*) Natoms
      call read_molecule(input_file, Natoms, coord, mass)
      DO i=1,Natoms
        write(6, *) (coord(i,j), j=1,3), mass(i)
      END DO
      call compute_distances(Natoms, coord, distance)
      DO i=1,Natoms
        write(6, *) (distance(i,j), j=1,Natoms)
      END DO
      Vljtot=V(epsilon, sigma, Natoms, distance)
      write(6,*) Vljtot
      Ttot=T(Natoms, velocity, mass)
      write(6,*) Ttot
      Etot=E(Vljtot,Ttot)
      write(6,*) Etot
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
      integer :: i,j,k,i_stat
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
          !converts Angstrom in nm
          distance(j,i)= distance(j,i)*0.1d0
        end do
      end do
      write(6,*) 
      end subroutine compute_distances

      double precision function V(epsilon, sigma, Natoms, distance) result(Vlj)
      implicit none
      double precision, intent(in) :: epsilon, sigma
      integer, intent(in) :: Natoms
      double precision, intent(in) :: distance(Natoms,Natoms)
      double precision :: Vljexp6
      integer :: i,j
      write(6,*) epsilon, sigma

      Vlj=0.d0
      do j=1,Natoms-1
        do i=j+1,Natoms
          Vljexp6=(sigma/distance(j,i))**6
          Vlj=Vlj+Vljexp6*Vljexp6-Vljexp6
        end do
      end do
      Vlj=4.d0*epsilon*Vlj
      end function V

      double precision function T(Natoms, velocity, mass) result(K)
      implicit none
      integer, intent(in) :: Natoms
      double precision, allocatable, intent(out) :: velocity(:,:)
      double precision, intent(in) :: mass(Natoms)
      integer :: i_stat

      allocate(velocity(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for velocity!"
                stop
      end if
      DO j=1,3
        DO i=1,Natoms
          velocity(i,j)=0.d0 
        END DO
      END DO
      K=0.d0
      DO i=1,Natoms
        K=K+mass(i)*(velocity(i,1)**2+velocity(i,2)**2+velocity(i,3)**2)
      END DO
      K=0.5d0*K
      end function T

      double precision function E(Vljtot,Ttot) result(Energy)
      implicit none
      double precision, intent(in) :: Vljtot,Ttot
      Energy=Vljtot+Ttot
      end function E

end program dynamics
