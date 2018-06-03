!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! 3D oil reservoir Simulator project 
!! by Mehrdad Gharib Shirangi From Spring 2010 until Fall 2010
!! mehrdad.ghsh@gmail.com
!! Copy Right by Mehrdad Gharib Shirangi
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		   
    subroutine Simulate

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
	use globalvar
	Use Grd_subs
	Use Trans
	USE Well
	Use Mat_Calc    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	
		Real*8, Dimension(Ntotal) :: permx, permy, permz , S_w     	! PS_w is Past Sw vector 
		Real*8, Dimension(Ntotal + Nwell) :: p					 ! P_p is past pressure vector, while P is the vector we want to put the pressure into
		integer :: alfa , B_index , i , j , k 
			    
		  ! putting pressure and Saturation values into Pressure and Saturation Matrices
		  Do k = 1 , Nz
			 Do j = 1 , Ny
				Do i = 1 , Nx      
					B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
					alfa = 2 * (B_index - 1 ) + 1                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
					YY(alfa) = past_p(B_index) 
					YY(alfa + 1) = past_Sw(B_index) 
				EndDO
			 EndDO
		  EndDO	   
		  Do i=1, Nwell
			 B_index = Ntotal + i ! B_index is the block number, while 
			 alfa = 2 * Ntotal + i                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
			 YY(alfa) = pressure(B_index) 
		  EndDo
		        		   
	    ! 3) solving the equation and updating the vector of pressures and saturations   
	    call wellIndex_cons
	    Call wellIndex
		Converged = .False.
	    i_iter = 0
	    Do while ( .NOT. Converged ) 
		   i_iter = i_iter  + 1
		   
		   call solve_Equ		   			
		   call con_chek ! convergence chek
			 ! call Calc_rate is included in con_check		   
 	    Enddo
 	      		        			
   Endsubroutine Simulate

!*******************************************************************************************
!*******************************************************************************************   
!*******************************************************************************************   

   subroutine Solve_Equ
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
	use globalvar
	Use Grd_subs
	Use Trans
	USE Well
	Use Mat_Calc    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   
   integer :: B_index, alfa , i
   real :: RUN_TIME1 , RUN_TIME2 , RUN_TIME3
   !! Solver Variables 
   integer :: mr1 , itr_max1  
   Real*8 ::  tol_abs1 , tol_rel1


!   Parameters used in Solver 
    itr_max1 = 50
    mr1 = 30
    tol_abs1 = 1.0D-05
    tol_rel1 = 1.0D-05
        				

      Allocate(del_y(N_equ))
      del_y = 0
		   Call CPU_TIME(RUN_TIME1)
	  
        call matrix
        
		DeAllocate( k_rw , k_ro , B_o , B_w , mu_o , mu_w , DkrodS , DkrwdS )   
		Call CPU_TIME(RUN_TIME2)      
		
		call mgmres_st ( N_equ , nz_num1, ia1, ja1, a1, del_y, - FCN, itr_max1, mr1, tol_abs1, tol_rel1 )   
		
!		Call SOLVERN(NX,NY,NZ,NTOTAL,Nwell,MGS,A,X,B,m1,isol,maxiter, eps, iter)
		
		Call CPU_TIME(RUN_TIME3)            
		YY = YY + del_y
      
		maxFCN = maxval(abs(del_y))
		maxentry = maxloc(abs(del_y))


      DeAllocate ( ia1)
      DeAllocate ( ja1)
      DeAllocate ( a1)
      DeAllocate(del_y)
      
	  Allocate( k_rw(Ntotal), k_ro(Ntotal) , B_o(Ntotal) , B_w(Ntotal) , mu_o(Ntotal) , mu_w(Ntotal) , DkrodS(Ntotal) , DkrwdS(Ntotal) )
            
      
      ! putting pressure and Saturation values into Pressure and Saturation Matrices
      Do k = 1 , Nz
         Do j = 1 , Ny
            Do i = 1 , Nx      
				B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
				alfa = 2 * (B_index - 1 ) + 1                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
				pressure(B_index) = YY(alfa)
				Sw(B_index) = YY(alfa + 1)
				if (Sw(B_index) > 1) then 
					!write(*,*) 'Oh! , Sw =' , Sw(B_index) , '   Block#' , B_index
					Sw(B_index) = 1.0
				elseif ( Sw(B_index) < 0.1)	then 
					!write(*,*) 'Oh! , Sw =' , Sw(B_index) , '   Block#' , B_index
					Sw(B_index) = 0.1
				EndIf
				W_press(B_index) = Pressure(B_index) - PCO(B_index)
				
				k_ro(B_index) = k_r( B_index , 'o')
				k_rw(B_index) = k_r( B_index , 'w')
				
	            B_o(B_index) = B_f( pressure(B_index) , 'o')
				B_w(B_index) = B_f( W_press(B_index) , 'w') 	              			  
				
				mu_o(B_index) = mu_f ( pressure(B_index),'o')				
				mu_w(B_index) = mu_f ( W_press(B_index),'w')				
				
				DkrodS(B_index) = Dkr_dS(B_index,'o')
				DkrwdS(B_index) = Dkr_dS(B_index,'w')					
						
            EndDO
         EndDO
      EndDO
      
      Do i=1, Nwell
         B_index = Ntotal + i ! B_index is the block number, while 
         alfa = 2 * Ntotal + i                ! alfa is the location in the Jacobian Matrix , I will use B_index to point the Block Number and pressure of that Block
         pressure(B_index) = YY(alfa)
      EndDo
      
      Call wellIndex
	  call check_BHP 

!!      Open(3,DEFAULTFILE='.\sim_output\' , FILE='RunTime.DAT', status='unknown' , access = 'append')
!!      write(3,'(F15.3 , 2x , F15.3 , 2x  , F15.3)') RUN_TIME1 , RUN_TIME2 , RUN_TIME3
!!      close(3) 

     
   Endsubroutine Solve_Equ
!*******************************************************************************************   
!*******************************************************************************************   
!*******************************************************************************************   
    
    subroutine Simulator(permx , permy , permz , porosity) !!  , pred_data
    !**************************************************************
    !
    ! Written by Mehrdad Gharib Shirangi , 
    ! 4th July 2010
    !
    !**************************************************************    
    ! Important:
    !			 When this sub is called, the Result of Rerun which is Pressure and Sw, is in Global Vectors of "Pressure"  and "Sw"
    !			 and all P and Sw would be written into files 
    !**************************************************************    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    
	use globalvar
!!	use EnKFglobalvar
	Use Grd_subs
	Use Trans
	USE Well
	Use Mat_Calc    
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    
    use Globalvar
    use Trans
    USE Well
    USE Derivative
    !**************************************************************        
!!!	Real*8,intent(out) :: pred_data(N_obs) 
	Real*8,intent(in) :: permx(Ntotal), permy(Ntotal), permz(Ntotal) , porosity(Ntotal)
		    
!!    REAL*8 :: DEPTH  ,&!! Updated Woc depth of the realization
!!			  final_time ! the final time we want to have rerun form time zero up to
	integer :: alfa , B_index , i , j , k 
	integer :: dt_Number
!!	character(len=50) :: STR , charac
	character*50 :: STR , charac
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	N_sim = N_sim + 1		
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	write(charac,'(I<7>)')j_ens
	charac='d_pre'//charac(6:)
	charac=charac(1:7)//'.dat'
	Open(212,DEFAULTFILE='.\sim_output\',FILE=charac, status='unknown')	!  , access = 'append'

	write(charac,'(I<11>)')j_ens
	charac='d_pre_Obs'//charac(10:)
	charac=charac(1:11)//'.dat'
	Open(215,DEFAULTFILE='.\sim_output\',FILE=charac, status='unknown')	! , access = 'append'
	Write(215 , '(a, I5)') 'N_lanczos = ' , N_lanczos
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^		
!!	kx = permx
!!	ky = permy
!!	kz = permz
!!	phi = porosity

	Schedule = Schedule_T   ! It is possible that the model can not produce with the specified BHP, which it has
							! to produce with constant BHP
    
    ! This subroutine runs the Simulator from time zero 
	if (USEIP == 'NOT') then
!!		We can add Updated Depth to our Model Parameters 	
!!		OWC_depth = depth 
		call ref_press_table
	EndIF

	call initialize


	
	Do k = 1 , Nz
		Do j = 1 , Ny
			Do i = 1 , Nx      
				B_index = (k-1) * Nx * Ny + (j-1) * Nx + i ! B_index is the block number, while 
				DkrodS(B_index) = Dkr_dS(B_index,'o')
				DkrwdS(B_index) = Dkr_dS(B_index,'w')				
			EndDO
		EndDO
	EndDO
	call wellIndex_cons
	Call wellIndex
		
	call initialize_BHP !! This Subroutine is required just because a value should be specified for initial BHP, however it may be just arbitrary as Cal_Delpwf will take care of it
	If (Nz .NE. 1) then
		call Calc_DeltaPwf
	Else
		Del_Pwf = 0
	EndIf	

	past_p = pressure
	

	
	past_Sw = Sw
	Past_pw = W_press
		
	time = 0.0
	time_index=1
	sche_N = 1
	dt_Number = 0
	datasche_N = 1
	Obs_Index = 1 ! The next Observation that we are going to have is Obs_Index ( 1st one or ...), this is necessary for the matrix Nt_dobs		
	
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
!	WRITING  intial P and Sw into file
	IF (j_ens == j_base ) Then
	
		Call System('del .\P_SW\File_Name.DAT')			
		
		str='P_SW000'
		str = str(1:7)//'.DAT'		
		Open(10,DEFAULTFILE='.\P_SW\' , FILE=str, status='unknown',form='UNFORMATTED')		  		  
		write(10) Pressure , Sw
		CLOSE(10) 	
		
!!!	! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         						
!!!		Open(11,FILE='BASED_SVD.DAT', access = 'append')
!!!		write(11,*) '--------------------------------------------'
!!!		Write(11,'(<2>(a,I5))' ) 'j_base =' , j_base , ', j_Ens =' , j_ens		
!!!		write(11,*) '--------------------------------------------'		
!!!		Close(11) 					
!!!	! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^         		
	EndIf
	! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!	Open(11,FILE='Sim_J_Ens.DAT', access = 'append')
!!!	write(11,*) '--------------------------------------------'
!!!	Write(11,'(<2>(a,I5))' ) 'j_base =' , j_base , ', j_Ens =' , j_ens		
!!!	write(11,*) '--------------------------------------------'		
!!!	Close(11) 					
	
	
	Do while ( (time+0.001) .LE. Cu_t(N_Dt))
	   ! 1) Time step considerations
		if ((time + time_step)>= Cu_t(time_index)) Then
			d_time = Cu_t(time_index) - time
			time = Cu_t(time_index)
			time_index = time_index + 1
		Else 
			time = time + time_step
			d_time = time_step
		EndIf	
		dt_Number = dt_Number + 1
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^		
		! finding the appropriate Schedule
		i= sche_N +1
		if (time > t_Schedule(i)) then
			sche_N = i ! 
		Endif   
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        i= datasche_N
        if (time > DATATIME(i)) then
            datasche_N = i + 1 !
        Endif  
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
		! 2) putting updated vectors into vectors of previous time_step
		past_p = pressure
		past_Sw = Sw
		Past_pw = W_press

		Call Simulate
		
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		! 5) writing data + P + Sw to files 
		If ( j_ens == j_base ) Then !! for the MAP estimate
			write(str,'(I<int(log10(real(dt_Number)))+1>)')dt_Number
			If(dt_Number < 10)then
				str='P_SW00'//str(1:)
			Elseif( dt_Number < 100)then
				str='P_SW0'//str(1:)
			Else
				str='P_SW'//str(1:)
			Endif
			str = str(1:7)//'.DAT'		
			Open(10,DEFAULTFILE='.\P_SW\' , FILE=str, status='unknown',form='UNFORMATTED')		  		  
			write(10) Pressure , Sw
			CLOSE(10) 	
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		!! We save the d_step and the time in the file , so that in Adjoint we can read that time and the Number to open P_SW file
			Open(10,DEFAULTFILE='.\P_SW\' , FILE='File_Name.DAT', status='unknown', access='append')	! ,form='UNFORMATTED'
			Write(10, '(F10.1,I<5>)')	time , dt_Number
		EndIf		
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		call output							
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
!!		if    ( ( mod(time, T_OBS) == 0.0 ) .AND. (time .LE. T_End) .AND. (Gen_Obs)) then !Generate Observed Data from True Data   
		If (Gen_Obs) Then
			if    ( ( time == time_obs(obs_index) )) then !Generate Observed Data from True Data 		
				call Gen_dobs
				IF ( Obs_Index < N_time_obs ) then
					Obs_Index = Obs_Index + 1						
				EndIF
			EndIF
		EndIf
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		DeAllocate(d_pre)
		DeAllocate( Well_No)
		DeAllocate( Data_type)				
		DeAllocate( diagonals)   
		! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%						
	EndDo
	! Saving the Satuation MAP
	If (Gen_Obs) Then
		Open(2,DEFAULTFILE='.\sim_output\',FILE='True_Sw.DAT', status='unknown')	! , access = 'append'
		Do i =1 , Ntotal
			Write(2 , '(F15.7)') Sw(i)
		EndDo
		close(2)

	EndIF
	
	close(212) ; close(215)
	        
	Endsubroutine Simulator
!*******************************************************************************************  
