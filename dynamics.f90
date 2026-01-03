      program dynamics
!Allocating and freeing 2D arrays
      implicit none
      integer :: input_file, Natoms,i,j, i_stat
      double precision, allocatable :: coord(:,:)
      double precision, allocatable :: mass(:)
      double precision, allocatable :: distance(:,:)
      double precision, allocatable :: velocity(:,:)
      double precision, allocatable :: acceleration(:,:)
      double precision, allocatable :: acc(:,:)
      double precision :: epsilon, sigma, Vljtot, Ttot, Etot
! In this program will be considered 5 Xenon atoms, the epsilon value is 1.77 kJ/mol
! the sigma value is 4.10 Angstroms
! time must be expressed in picoseconds and distance in nanometers
      input_file = 9
      epsilon = 1.77d0
      sigma = 4.10d0
      sigma=sigma*0.1d0
      open(input_file, file='inp.txt')
      open(10,file='dynamics.xyz')
      Natoms = read_Natoms(input_file)
      write(6,*) Natoms
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
      call read_molecule(input_file, Natoms, coord, mass)
      DO i=1,Natoms
        write(6, *) (coord(i,j), j=1,3), mass(i)
      END DO
      allocate(distance(Natoms,Natoms), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for distance!"
                stop
      end if
      call compute_distances(Natoms, coord, distance)
      DO i=1,Natoms
        write(6, *) (distance(i,j), j=1,Natoms)
      END DO
      Vljtot=V(epsilon, sigma, Natoms, distance)
      write(6,*) Vljtot
      allocate(velocity(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for velocity!"
                stop
      end if
      DO j=1,3
        DO i=1,Natoms
          velocity(i,j)=0.d0
        END DO
      ENDDO
      Ttot=T(Natoms, velocity, mass)
      write(6,*) Ttot
      Etot=E(Vljtot,Ttot)
      write(6,*) Etot
      allocate(acceleration(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for acceleration!"
                stop
      end if
      call compute_acc(Natoms, mass, acceleration, coord, distance, sigma, epsilon)
      DO i=1,Natoms
        write(6, *) (acceleration(i,j), j=1,3)
      END DO
      allocate(acc(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for acc!"
                stop
      end if
      call Verlet(Natoms, coord, mass, distance, acceleration, acc, sigma, epsilon)     


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
      double precision, intent(out) :: coord(Natoms,3)
      double precision, intent(out) :: mass(Natoms)
      DO i=1,Natoms
        read(input_file, *) (coord(i,j), j=1,3), mass(i)
      END DO
      DO j=1,3
       DO i=1,Natoms 
        coord(i,j)=0.1d0*coord(i,j)
       enddo
      enddo

      end subroutine read_molecule

      subroutine Verlet(Natoms, coord, mass, distance, acceleration,acc, sigma, epsilon)
      implicit none
      integer, intent(in) :: Natoms
      integer i, j , k, g
      double precision deltat, Vljtot, Ttot, Etot
      double precision, intent(in) :: epsilon, sigma
      double precision, intent(inout) :: coord(Natoms,3)
      double precision, intent(in) :: mass(Natoms)
      double precision, intent(inout) :: distance(Natoms,Natoms)
      double precision, intent(inout) :: acceleration(Natoms,3)
      double precision :: acc(Natoms,3)
      deltat=0.05d0
      g=0
      DO k=0,2999      
        DO j=1,Natoms
          DO i=1,3 
           acc(j,i)=acceleration(j,i)
           coord(j,i)=coord(j,i)+velocity(j,i)*deltat+acc(j,i)*(deltat**2/2.d0)
           enddo
        enddo
        call compute_acc(Natoms, mass, acceleration, coord, distance, sigma, epsilon)
        DO j=1,Natoms
          DO i=1,3
           velocity(j,i)=velocity(j,i)+0.5d0*deltat*(acc(j,i)+acceleration(j,i))
          enddo
          write(6,*) (coord(j,i), i=1,3)
        enddo
        Vljtot=V(epsilon, sigma, Natoms, distance)
        Ttot=T(Natoms, velocity, mass)
        Etot=E(Vljtot,Ttot)
        if(k.eq.g) then
          g=g+10
          write(10,*) Natoms
          write(10,*)  'E=', Etot, 'V=', Vljtot, 'T=',Ttot
          DO j=1,Natoms
            write(10,*) 'Xe', coord(j,1)*10d0, coord(j,2)*10d0, coord(j,3)*10d0
          end do
        endif
      enddo
      end subroutine Verlet

      subroutine compute_acc(Natoms, mass, acceleration, coord, distance, sigma, epsilon)
      implicit none
      integer, intent(in) :: Natoms
      double precision acc, Uexp6, accith
      integer i, j , k
      double precision, intent(in) :: epsilon, sigma
      double precision :: coord(Natoms,3)
      double precision, intent(in) :: mass(Natoms)
      double precision :: distance(Natoms,Natoms)
      double precision, intent(out) :: acceleration(Natoms,3)

      call  compute_distances(Natoms, coord, distance)
      DO k=1,Natoms
        DO j=1,3
          acc=0.d0
          DO i=1,Natoms
            if(i.ne.k) then
              Uexp6=(sigma/distance(k,i))**6
              accith=(Uexp6-2.d0*Uexp6*Uexp6)*((coord(k,j)-coord(i,j))/distance(k,i))*(1.d0/distance(k,i))
              acc=acc+accith
            endif
          enddo
            acceleration(k,j)=-24.d0*epsilon*acc/mass(k)
        enddo
      enddo
      end subroutine compute_acc

!subroutine that calculates the distances between each possible pair 
      subroutine compute_distances(Natoms, coord, distance)
      implicit none
      integer :: i,j,k,i_stat
      double precision :: diff
      integer, intent(in) :: Natoms
      double precision, intent(in) :: coord(Natoms,3)
      double precision, intent(out) :: distance(Natoms,Natoms)

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
      double precision, intent(in) :: velocity(Natoms,3)
      double precision, intent(in) :: mass(Natoms)
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
