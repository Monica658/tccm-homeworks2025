      program dynamics
!dynamic program that allows to write the trajectory of Natoms, with Nsteps!
!given the initial coordinates and with the starting velocities set to zero!
      implicit none
      integer :: input_file, Natoms,i,j, i_stat, Nstep
      character(len=3) :: symbol
      double precision, allocatable :: coord(:,:)
      double precision, allocatable :: mass(:)
      double precision, allocatable :: distance(:,:)
      double precision, allocatable :: velocity(:,:)
      double precision, allocatable :: acceleration(:,:)
      double precision, allocatable :: acc(:,:)
      double precision :: epsilon, sigma, Vljtot, Ttot, Etot, deltat
      input_file = 9
!-----input and output files------------------------------------------------------------------------------------------!
      open(input_file, file='inp.txt')
      open(8, file='param.inp')
      open(10,file='dynamics.xyz')
!---------------------------------------------------------------------------------------------------------------------!
      read(8,*) epsilon
      read(8,*) sigma
      read(8,*) deltat
      read(8,*) Nstep
      read(8,*) symbol
      sigma=sigma*0.1d0 !unit conversion from angstrom to nanometers
      Natoms = read_Natoms(input_file) 

!-----Matrices are now allocated after reading the number of atoms----------------------------------------------------!
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
      allocate(distance(Natoms,Natoms), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for distance!"
                stop
      end if
      allocate(velocity(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for velocity!"
                stop
      end if
      allocate(acceleration(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for acceleration!"
                stop
      end if
      allocate(acc(Natoms,3), stat=i_stat)
      if (i_stat /= 0) then
                print *,  "Memory allocation failed, for acc!"
                stop
      end if
!-----end of the allocation part--------------------------------------------------------------------------------------!

      call read_molecule(input_file, Natoms, coord, mass)
      call compute_distances(Natoms, coord, distance)
      Vljtot=V(epsilon, sigma, Natoms, distance)
      DO j=1,3
        DO i=1,Natoms
          velocity(i,j)=0.d0!velocities are initially set to zero
        END DO
      ENDDO
      Ttot=T(Natoms, velocity, mass)
      Etot=E(Vljtot,Ttot)
      call compute_acc(Natoms, mass, acceleration, coord, distance, sigma, epsilon)
      call Verlet(Natoms, coord, mass, distance, velocity, acceleration, acc, sigma, epsilon, deltat,Nstep, symbol) 
      
!----starts the part of the program where that builds the trajectory-----------------------------------------------!
contains
      integer function read_Natoms(input_file) result(atoms)
      implicit none
      integer, intent(in) :: input_file
      read(input_file, *) atoms
      end function read_Natoms

      subroutine read_molecule(input_file, Natoms, coord, mass)
      implicit none
      integer :: i,j
      integer, intent(in) :: input_file
      integer, intent(in) :: Natoms
      double precision, intent(out) :: coord(Natoms,3)
      double precision, intent(out) :: mass(Natoms)
      do i=1,Natoms
        read(input_file, *) (coord(i,j), j=1,3), mass(i)
      end do
      do i=1,3
       do j=1,Natoms
        coord(j,i)=0.1d0*coord(j,i) !unit conversion from angstrom to nanometers
       enddo
      enddo
      end subroutine read_molecule
!---------------------------------------------------------------------------------------------------------------------!
!-----subroutine that calculates the trajectory of the Natoms, updating each step:------------------------------------!
!-----the coordinates with the velocities calculated in the previous step---------------------------------------------!
!-----the accelerations by calling the subroutine compute_acc---------------------------------------------------------!
!-----the velocities with the acceleration values of the previous and actual step-------------------------------------!
!-----The deltat is defined in the input file Nstep iterationsa are done----------------------------------------------!
!-----The subroutine writes every ten step the calculated values in the output file-----------------------------------!
      subroutine Verlet(Natoms, coord, mass, distance, velocity, acceleration,acc, sigma, epsilon,deltat,Nstep, symbol)
      implicit none
      integer :: i, j , k, g
      character(len=3) :: symbol
      double precision :: Vljtot, Ttot, Etot
      integer, intent(in) :: Natoms, Nstep
      double precision, intent(in) :: epsilon, sigma, deltat
      double precision, intent(inout) :: coord(Natoms,3)
      double precision, intent(in) :: mass(Natoms)
      double precision, intent(inout) :: distance(Natoms,Natoms)
      double precision, intent(inout) :: acceleration(Natoms,3)
      double precision, intent(inout) :: velocity(Natoms,3)
      double precision :: acc(Natoms,3)
      g=0
      DO k=0,Nstep-1     
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
        enddo
        Vljtot=V(epsilon, sigma, Natoms, distance)
        Ttot=T(Natoms, velocity, mass)
        Etot=E(Vljtot,Ttot)
        if(k.eq.g) then
          g=g+10
          write(10,*) Natoms
          write(10,*)  'E=', Etot, 'V=', Vljtot, 'T=',Ttot
          DO j=1,Natoms
            write(10,*) symbol, coord(j,1)*10d0, coord(j,2)*10d0, coord(j,3)*10d0 ! from nanometers to angstrom
          end do
        endif
      enddo
      end subroutine Verlet
!---------------------------------------------------------------------------------------------------------------------!

!-----subroutine that calculates the accelerations--------------------------------------------------------------------!
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
      do k=1,Natoms
        do j=1,3
          acc=0.d0
          do i=1,Natoms
            if(i.ne.k) then !condition to exclude the distance diagonal elements that are equal to zero
              Uexp6=(sigma/distance(k,i))**6
              accith=(Uexp6-2.d0*Uexp6*Uexp6)*((coord(k,j)-coord(i,j))/distance(k,i))*(1.d0/distance(k,i))
              acc=acc+accith
            endif
          enddo
            acceleration(k,j)=-24.d0*epsilon*acc/mass(k)
        enddo
      enddo
      end subroutine compute_acc
!---------------------------------------------------------------------------------------------------------------------!

!-----subroutine that calculates the distances between each possible pair of atoms------------------------------------!
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
            distance(j,i)=distance(i,j) !symmetric matrix
          end if
          diff=0.d0
          DO k=1,3
            diff=diff+(coord(j,k)-coord(i,k))**2
          end do
          distance(j,i)=sqrt(diff)
        end do
      end do
      end subroutine compute_distances
!---------------------------------------------------------------------------------------------------------------------!

!-----double precison function that calculate the Lennard Jones potential based on the previously calculated distances!
      double precision function V(epsilon, sigma, Natoms, distance) result(Vlj)
      implicit none
      integer :: i,j
      double precision, intent(in) :: epsilon, sigma
      integer, intent(in) :: Natoms
      double precision, intent(in) :: distance(Natoms,Natoms)
      double precision :: Vljexp6
      Vlj=0.d0
      do j=1,Natoms-1
        do i=j+1,Natoms
          Vljexp6=(sigma/distance(j,i))**6
          Vlj=Vlj+Vljexp6*Vljexp6-Vljexp6
        end do
      end do
      Vlj=4.d0*epsilon*Vlj
      end function V
!---------------------------------------------------------------------------------------------------------------------!

!-----double precision function that calculates the kinetic energy based on the atoms mass and velocities-------------!
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
!---------------------------------------------------------------------------------------------------------------------!

!-----double precision function that calculates the total energy------------------------------------------------------!
      double precision function E(Vljtot,Ttot) result(Energy)
      implicit none
      double precision, intent(in) :: Vljtot,Ttot
      Energy=Vljtot+Ttot
      end function E
!---------------------------------------------------------------------------------------------------------------------!

end program dynamics
