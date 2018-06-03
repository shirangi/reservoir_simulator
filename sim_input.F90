Subroutine Sim_Input
!!
Use globalvar
!-------------------------------------Variable Declaration Part-------------------------------------------------------
Implicit None
!

Character*80 chline
Character*13 ch1
Character*10 ch2
Character*20 ch3
	

Integer :: i,j , state




!-------------------------------------Execution Part------------------------------------------------------------------
Open(300,DEFAULTFILE='.\Code Input Files\',file='Sim_Input.dat',status='old')

!!***********************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='DIMENS') then 
      read(300,'(A80)')chline
      exit ! find the location
   endif
Enddo ! 
!!
read(300,*) Nx,Ny,Nz
!
rewind 300
!!
NTOTAL = Nx*Ny*Nz   !! number of model parameters to be calculated
!!***********************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='SWOF') then 
      exit ! find the location
   endif
Enddo ! 
!!

i = 0
! obtaining number of sw and krw and kro which we are reading
do while (chline(1:1)/='/')
   i = i+1
    read(300,'(A80)')chline
!
end do
n_Kr = i - 1
Allocate(Sw_mat(n_Kr))
Allocate(Kro(n_Kr))
Allocate(Krw(n_Kr))
Allocate(Pcow(n_Kr))
rewind 300

Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='SWOF') then 
      exit ! find the location, again!
   endif
Enddo ! 

j = 1
do while (j < (n_kr + 1))

   Read(300, FMT=*) Sw_mat(j), Krw(j), Kro(j),Pcow(j)
   j = j + 1   
!
end do
rewind 300

!                    Do i=1,n_Kr
!                       write(*,*) Kr_Sw(i,3)
!                    end do
!                    pause



!                                                      write(*,*) 'kar mikone!'
!                                                      pause 
!!***********************************************

! Initial Reservoir Pressure and the reference depth of that
Do 
   read(300,'(A80)')chline
   if(chline(1:5)=='EQUIL') then 
      exit ! find the location
    endif
Enddo ! 
!!
Allocate(initialp(4))
read(300,*) initialp(1) , initialp(2) , initialp(3), initialp(4)
Datum_D = initialp(1)
Datum_P = initialp(2)
OWC_depth = initialp(3)
OWC_PCOW = initialp(4)
!
rewind 300

!!*************************************************

! density of Oil and Water
Do 
   read(300,'(A80)')chline
   if(chline(1:7)=='DENSITY') then 
      exit ! find the location
    endif
Enddo ! 
!!
read(300,*)den_o,den_w 
!
rewind 300

!                                                      write(*,*) 'kar mikone!3'
!                                                      pause 
!!*************************************************
! well locations
! at first finding number of data
Do 
   read(300,'(A80)')chline
   if(chline(1:8)=='WELSPECS') then 
      exit ! find the location
    endif
Enddo ! 
i=0
! obtaining number of Wells
do while (chline(1:1)/='/')
   i = i+1
    read(300,'(A80)')chline
!
end do
Nwell = i - 1
Allocate(Well_name(Nwell))
Allocate(Well_loc(Nwell,2))
! 
rewind 300
! reading well data
Do 
   read(300,'(A80)')chline
   if(chline(1:8)=='WELSPECS') then 
      exit ! find the location
    endif
Enddo ! 
Allocate (Well_k(Nwell,2))
j=1

do while (j .LE. Nwell )

   !Read(300, FMT=*) Well_name(j), ch2 , well_loc(j,1), well_loc(j,2),well_k(j,1),well_k(j,2)
   Read(300, FMT=*) Well_name(j) , well_loc(j,1), well_loc(j,2),well_k(j,1),well_k(j,2)
   j = j + 1   
!
end do

!                    do i = 1, Nwell
!                       Write(*,*) 'i=',well_loc(i,1) ,'    j=' ,well_loc(i,2)
!                    enddo
!                    pause
rewind 300
! *******************************************************
! Radius of wells
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='RADIUS') then 
      exit ! find the location
    endif
Enddo ! 
!!
allocate(r_w(Nwell))
read(300,*) r_w(1:Nwell)  
rewind 300

!!*************************************************
! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:3)=='SWI') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) SWI
 
rewind 300
! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:8)=='PRESSURE') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) Initial_P
 
rewind 300
! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:5)=='USEIP') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) USEIP
 
rewind 300

!!************************************************* 
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='USEISW') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) USEISW
 
rewind 300

!!************************************************* 
! Minimum and Maximum pressure for constant rate injection or production
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='MINBHP') then 
      exit ! find the location
    endif
Enddo ! 
!!
read(300,*) MinBHP
rewind 300

Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='MAXBHP') then 
      exit ! find the location
    endif
Enddo ! 
!!
read(300,*) MaxBHP
rewind 300

!!*************************************************
!!*************************************************

! 
Do 
   read(300,'(A80)')chline
   if(chline(1:9)=='time_step') then 
      exit ! find the location
    endif
Enddo ! 
!!
read(300,*) time_step
rewind 300
!!*************************************************
! reading Well type 
Do 
   read(300,'(A80)')chline
   if(chline(1:8)=='WellType') then 
      exit ! find the location
   endif
Enddo ! 
!!
allocate(W_type(Nwell))
Do i=1, Nwell
   read(300,*) j, W_type(i)  
EndDo
rewind 300
!!*************************************************

Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='PVDO') then 
      exit ! find the location
    endif
Enddo ! 
i=0
! obtaining number of Wells
do while (chline(1:1)/='/')
   i = i+1
    read(300,'(A80)')chline
!
end do
N_mup = i - 1
N_Bo = N_mup

Allocate(mu_p(N_mup,2))
Allocate(Bo_p(N_Bo,2))
! 
rewind 300

Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='PVDO') then 
      exit ! find the location
    endif
Enddo ! 
j=1

do while (j .LE. N_mup )

   Read(300, FMT=*) Bo_p(j,1),Bo_p(j,2), mu_p(j,2)  ! the first Column is pressure for both Viscosity and FVF
   mu_p(j,1) = Bo_p(j,1)
   j = j + 1 
!
end do

rewind 300

!!*************************************************

! PVT data of WATER, including C_w 
Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='PVTW') then 
      exit ! find the location
    endif
Enddo ! 
!!
allocate(pvtw(4))
read(300,*) pvtw(1:4)   ! consists of 4 entries :
                        ! The reference pressure, The water formation volume factor at the reference pressure , The water compressibility, The water viscosity
                            !
                            !write(*,*) pvtw(2), '       ',pvtw(4)
                            !pause
rewind 300
! *******************************************************

Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='ROCK') then 
      exit ! find the location
    endif
Enddo ! 
!!
read(300,*) P_ref  , C_ROCK
                            !
                            !write(*,*) pvtw(2), '       ',pvtw(4)
                            !pause
rewind 300


! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='MAXWOR') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) max_WOR
 
rewind 300
! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:3)=='OMB') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) OMB
 
rewind 300
! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='MRHS') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) MRHS
 
rewind 300
!!*************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:3)=='SOR') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) SOR
SW_AtSor = 1 - SOR
 
rewind 300


! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:4)=='MOMB') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) MOMB
 
rewind 300

! *******************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='M_ITER') then 
      exit ! find the location
    endif
Enddo ! 

read(300,*) M_ITER
 
rewind 300
! ******************************************************* 
 ! reading time steps
Do 
   read(300,'(A80)')chline
   if(chline(1:6)=='NTSTEP') then 
      exit ! find the location
    endif
Enddo ! 



read(300,*) N_Dt

Allocate(Deltat(N_Dt))
Allocate(Cu_t(N_Dt))
! 
rewind 300

Do 
   read(300,'(A80)')chline
   if(chline(1:5)=='TSTEP') then  !TSTEP
      exit ! find the location
    endif
Enddo ! 
Read(300, FMT=*) Deltat(1:N_Dt)

Cu_t(1)= Deltat(1)
Do i=2,N_Dt
   Cu_t(i)  = Cu_t(i-1) + Deltat(i)
End Do

rewind 300
!************************************************************
Do 
   read(300,'(A80)')chline
   if(chline(1:9)=='NSchedule') then 
      exit ! find the location
    endif
Enddo ! 
backspace 300
!!
read(300,*)ch1,NSchedule
!
read(300,'(A80)')chline 
read(300,'(A80)')chline 
read(300,'(A80)')chline 
!!
Allocate(Schedule(Nwell,5,NSchedule))
Allocate(t_Schedule(NSchedule+1))
!! First Column is Well number
!! Second Column is Constraint Types: 1: Oil flow Rate(qo) or Water injection Rate (qw_inj)
!!                                    2: Bottomhole Pressure (Pwf)
!! Third Column is Constraint
!! Fourth and Fifth Columns are Constraint Duration
!!t_Schedule is a vector which containes the begin and end times of all Schedules
Do i=1,NSchedule
   Do j=1, Nwell
      Read(300, FMT=*)Schedule(j,1,i),Schedule(j,2,i),Schedule(j,3,i),Schedule(j,4,i),Schedule(j,5,i)
   Enddo
Enddo
t_Schedule(1) = Schedule(1,4,1)
Do i=1,NSchedule
   t_Schedule(i+1) = Schedule(1,5,i)
endDo

ALLOCATE(Schedule_T(Nwell,5,NSchedule))
Schedule_T = Schedule
!!
rewind 300
!!***********************************************************************
 
Close(300)



! reading Delta x, Delta y  and Delta z
Open(15,DEFAULTFILE='.\Code Input Files\',file='geo.Dat',status='old')

ALLOCATE(Deltax(Nx))

Do 
   read(15,'(A80)')chline
   if(chline(1:3)=='DXV') then 
      exit ! find the location
    endif
Enddo ! 

Do i = 1 , Nx
    Read(15, FMT=*) Deltax(i)
end DO

rewind 15

! @         @          @          @

ALLOCATE(Deltay(Ny))

Do 
   read(15,'(A80)')chline
   if(chline(1:3)=='DYV') then 
      exit ! find the location
    endif
Enddo ! 

Do i = 1 , Ny
    Read(15, FMT=*) Deltay(i)
end DO

rewind 15

! @         @          @          @
Do 
   read(15,'(A80)')chline
   if(chline(1:3)=='DZ') then 
      exit ! find the location
    endif
Enddo ! 

ALLOCATE(Deltaz(Ntotal))

Read(15,*) Deltaz(1:Ntotal)
!				Open(55,file='DZ.Dat')
!				do i =1 , Ntotal
!					write(55,*) Deltaz(i)
!				EndDo
!				close(55)

rewind 15

! @         @          @          @

Do 
   read(15,'(A80)')chline
   if(chline(1:4)=='TOPS') then 
      exit ! find the location
    endif
Enddo ! 

ALLOCATE(tops(Ntotal))

Read(15,*) tops(1:Ntotal)

!				Open(55,file='tops.Dat')
!				do i =1 , Ntotal
!					write(55,*) tops(i)
!				EndDo
!				close(55)
! ********************************************

close(15)






EndSubroutine Sim_Input