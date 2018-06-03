
 Module Grd_subs
!************************************************************************************************************************************************************    
 use globalvar 
!************************************************************************************************************************************************************   
integer :: nz_num2
 !ia1, ja1, a1, x1, rhs1, itr_max1, mr1, tol_abs1, tol_rel1
Real*8, Allocatable, Dimension(:) :: a2

integer*4 ,Allocatable, Dimension(:) :: ia2, ja2

 
 contains 
!!!	DO WHILE (.NOT. EOF(2)) 

!************************************************************************************************************************************************************   
!************************************************************************************************************************************************************   
!************************************************************************************************************************************************************   
 Subroutine Read_Mtrue  !! This Subroutine is used Only when True Model is Avaiable and We want to generate Observed Data from The True Model

	 Implicit None

	integer :: i
! ***************************************	
	if(Nz > 1) then	 
		ALLOCATE( m_true( 2 * Ntotal ,1 ) )	
	Else		
		ALLOCATE( m_true( Ntotal ,1 ) )	
	EndIf
		
	 !! True Model is already in the folder of Prior Ensembles and The 
					   !! Program will read the true model and generate Observed data from it
	Open(2,DEFAULTFILE='.\Code Input Files\' , FILE = 'true_lnkh.Dat', status='old') 
	Do  i = 1, Ntotal
 		 read(2,*) m_true( i , [1])
	EndDO
	close(2)	    

! reading Porosity
	Open(12,DEFAULTFILE='.\Code Input Files\',file='true_poro.Dat',status='old')

	ALLOCATE(phi(Ntotal))
	Do i = 1 , Ntotal
		Read(12, FMT=*) phi(i)
	end DO
	close(12)

!!#################################################	
	
	If(Nz > 1) then
		Open(2,DEFAULTFILE='.\Code Input Files\' , FILE = 'True_lnkz.Dat', status='old') 
		Do  i = 1 + Ntotal ,  2 * Ntotal 
 			 read(2,*) m_true( i , [1])
		EndDO
		close(2)	    		 
	EndIf			 
		 
EndSubroutine Read_Mtrue	 
	 
!************************************************************************************************************************************************************   
!************************************************************************************************************************************************************   

	 Subroutine Gen_dobs
		 
		 ! ***************************************
		implicit none
		 ! ***************************************
		 integer :: i  , j   , alfa , i_1 , i_2 , j_1, j_2 , iwell, i_well, j_well , i_D  , j_D 
		 integer :: i_data , j_data 
		 Real*8 ,Allocatable, Dimension(:) :: CD , LCD, Z
		 character(len=30)::charac,charac1
		 
		 
		 ALLOCATE( CD( N_d ) )
		 ALLOCATE( d_obs( N_d , 1) )
		 
		 	 ! CD construction

		 Allocate(z(N_d))
		 	 
		 CD = 0
		 d_obs = 0	 

		 ! CD construction
		 j = 1
		 Do i = 1 , N_d
			If ( diagonals(i, 1) .NE. -1 )  then
				If ((diagonals(i, 1) == 0) .OR. (diagonals(i, 1) .LT. MINERROR_WOR) ) diagonals(i, 1) = MINERROR_WOR
				CD(i) = diagonals(i, 1)
			Else 				
				CD(i) = diagonals(i, 1)				
			EndIf
		 EndDo	 
		 
		 ! ***************************************	 	 	 
		 ! ***************************************	 	 	 
 		 call random1(N_d)
		 Open(1,FILE='randomdev.DAT', status='old')
		 read(1,*) z(1:N_d)
		 close(1)
		 Do i = 1 , N_d
			If ( diagonals(i, 1) .NE. -1 )  then			
				d_obs(i,1) = d_pre(i,1) + sqrt(CD(i)) * Z(i)
			Else 
				d_obs(i,1) = -1
			EndIf
		 EndDO
		 
		 Open(2,DEFAULTFILE='.\sim_output\' , FILE='Observed_Data.DAT', status='unknown' , access = 'append')
		 Do i = 1 , N_d
			If ( diagonals(i, 1) .NE. -1 )  then			
				write(2,'( < 3 >(F15.1 ,2x), < 2 >(F16.9 ,2x)  )') time, well_No(i,1) , Data_type(i,1) , diagonals(i,1),d_obs(i,1)
			EndIf
		 EndDo
		 close(2)	
		 		 
		 DeALLOCATE(d_obs)

	 EndSubroutine Gen_dobs

!************************************************************************************************************************************************************   

!*******************************************************************************************
!*******************************************************************************************
!*******************************************************************************************
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	 


EndModule Grd_subs
