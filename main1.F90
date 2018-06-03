	

program main

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! This is a 3D Simulator
!! The code was written in 4 month of work
!! by "Mehrdad Gharib Shirangi", From Spring 2010 until Fall 2010
!! mehrdad.ghsh@gmail.com
!! Copy Right by Mehrdad Gharib Shirangi, 2011
!! This Program Shall not be used without the Permission of either TUPREP director, Professor Albert Reynolds, or Mehrdad Gharib Shirangi
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
use globalvar
Use Grd_subs
Use Trans
USE Well
Use Mat_Calc
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
implicit none
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
integer :: i , j , k , B_index , ii
character*15 ::str
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^  SVD Input Parameters   ^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Integer:: L , n , i2,j1,j2,l1,l2,l3,l4,state 

REAL(8) :: time_end , gama_max

!! SVD variables
Integer::itmax, nsing, nsing_min, nsing_max, itmax_bidig, info, ounit,nlanc        
REAL(8):: sing_cut, min_sing_cut, sing_tol, toloff, toldig
REAL(8), Allocatable, Dimension(:,:):: usvec, vsvec
REAL(8), Allocatable, Dimension(:):: yvec, wvec, singold,alfa, beta, vaux, singval 
Integer, Allocatable, Dimension(:):: perm
Integer , Allocatable, Dimension(:):: Sec_Iter
Integer:: I_SV_Change,n_i_sv,n_SV_max  
!!------------------------------------------------------------------------------

Real(8), Allocatable, Dimension(:)::m_t,RHS,del_m,udotRHS,vT_dotm_t,delm_tilda , Om_ml, Om_m , D_scale , del_m_scaled , C_MAP
REAL(8)::Norm1, Norm2,dnorm2, d_m, d_Om ,O_m
REAL(8)::&
	      step0,&      ! Initial step size
          Epm,Epmf,&   ! Initial and final tolerance for criterion on model change
	      EpOm,EpOmf,& ! Initial and final tolerance for criterion on function change
		  !OmL,&        ! Objective function value
		  O_ml,O_mll,dir_max, &
		  Std_Lnk      ! Standard Deviation of lnk
INTEGER::IERR 
Real(8)::step , sum
!!-----------------------------------------------------------------------------
!! Varible declaration (sampling from posterior pdf)
Integer, Allocatable, Dimension(:)::Number
Integer:: &    !Number of Realizations
         I_SVD&    !Index for choosing a proper model to obtain SVD
         , max_LM_iter
!!
real*8, Allocatable, Dimension(:,:)::m_r,g_re , m_tempt
real*8, Allocatable, Dimension(:)::Z,Om,Oml,Od_ml,On_m,On_ml,gama,m_ave,m_avel,g_ave,m_aEnKF , Od_m , alpha
character(len=30)::charac
real*8::d_ml,d_Oml,em,eOm,Pre_Time

real*8, Allocatable, Dimension(:):: u, Om_re , step_j
real*8::Od_damp,Od_dampl , Sum1
real*8, Allocatable, Dimension(:):: gama_tempt
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! Added by Mehrdad:
real*8 :: gama_i , s_cut_M , gama_inc_f , gama_dec_f , step_c_f , step_i_f ,  min_step , max_step
Logical :: GN , N_SING  
Real*8	:: Od_mll , Om_mll , Min_GAMMA , dObj_ISVD , Factor , coef
Integer :: LM_counter  , &
		   nsing_Limit , n_Inc_Sing , PreI_SVD
Integer :: k_index
! j_base is the realization which the Lanczos is based on 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! Solver :
Integer ::	  xcounter,  ycounter
!------------------------------------------MAIN CODE----------------------------------------------------------
		
!-------------------------------------
	
	CU_TIME = 0.0 ; RUN_TIME = 0.0 ; time_begin= 0.0 ;
	call CPU_TIME(time_begin)
	
	!! Reading Input Parameters of the simulator			
	call sim_input

	Call Read_OBSERVED_Input	
	
	!===================
	Gen_obs = .True.
	!===================

	N_equ = 2 * Ntotal + Nwell
	sche_N = 1

	FOPR = 0 
	Allocate (Zblock(Ntotal))
	Zblock = tops + 1.0/2.0 * Deltaz ! all these 3 matrices have one dimention of N_total


	ALLOCATE (YY(N_equ))
	Allocate( k_rw(Ntotal), k_ro(Ntotal) , B_o(Ntotal) , B_w(Ntotal) , mu_o(Ntotal) , mu_w(Ntotal) , DkrodS(Ntotal) , DkrwdS(Ntotal) )
	Allocate(Well_index_cons(Nwell , Nz))
	Allocate(R_equ(Nwell , Nz))
	Allocate( KH(Nwell , Nz) )
	Allocate( W_index_o(Nwell , Nz) )
	Allocate( W_index_w(Nwell , Nz) )
	Allocate (oilrate(Nwell))
	Allocate (waterrate(Nwell))
	Allocate (WOR(Nwell))
	Allocate(Pressure(Ntotal + Nwell), Past_p(Ntotal + Nwell), New_p(Ntotal + Nwell))
	Allocate(W_press(Ntotal) , Past_pw(Ntotal) )
	Allocate(Sw(Ntotal))
	Allocate(past_Sw(Ntotal))
	Allocate(qoL(Nwell, Nz))
	Allocate(qwL(Nwell, Nz))
	ALLOCATE (FCN(N_equ))
	ALLOCATE(Kz(Ntotal),Kx(Ntotal),Ky(Ntotal) )
	Allocate( Del_pwf(Nwell, Nz) , maxentry(1),perm_x(Ntotal),perm_y(Ntotal),perm_z(Ntotal))                        

	
	Open(211,DEFAULTFILE='.\sim_output\' , FILE='RunTime_Step.DAT', status='unknown')
	
		
	!	We open and close the following file, so that it won't have data written in it from previous Runs
	Call System('del .\P_SW\File_Name.DAT')		
	Call System('del .\sim_output\Predicted_Data.DAT')		
	Call System('del .\sim_output\Predicted_Data_Obs.DAT')		
	Call System('del .\sim_output\Observed_Data.DAT')		
	Call System('del .\sim_output\True_Data.DAT')		
	Call System('del .\sim_output\True_FOPR.DAT')		
	
		
		
	if (USEIP == 'NOT') then
		! If USEIP is .True. we use a constant P and Sw for all gridblocks
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


	call  Read_Mtrue
	Do i = 1 , Ntotal
		perm_x(i) = exp( m_true( i , 1 ) )
	EndDo
	If (Nz > 1)	then 
		Do i =  1 , Ntotal
			perm_z(i ) = exp( m_true( i + Ntotal , 1 ) )
		EndDo
	Else
		Perm_z = 0.25 * perm_x
	EndIf		
	
	Schedule = Schedule_T			
	Perm_y = perm_x
	
	kx = perm_x		
	ky = perm_y
	kz = perm_z		
	Call Simulator(kx , ky , kz , phi) 



	close(211) 	   	     
	


!*******************************************************************************************
!*******************************************************************************************
contains 
	Subroutine Read_OBSERVED_Input

	implicit none
	!
	Character*80 chline
	Character*13 ch1
	!
	!read the parameter data from inverse.par 
	!
	Open(1000,DEFAULTFILE='.\Code Input Files\',file='Observ_Input.par',status='old')
	!!
! *******************************************************
	Do 
	   read(1000,'(A80)')chline
	   if(chline(1:8)=='SIGMA_QO') then 
		  exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, SIGMA_QO

	rewind 1000
! *******************************************************
	Do 
	   read(1000,'(A80)')chline
	   if(chline(1:9)=='SIGMA_WOR') then 
		  exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, SIGMA_QO

	rewind 1000
! *******************************************************
	Do 
		read(1000,'(A80)')chline	   
		if(chline(1:8)=='SIGMA_QW') then 	   
		  exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, SIGMA_QW

	rewind 1000

! *******************************************************
	Do 
		read(1000,'(A80)')chline	   
		if(chline(1:6)=='VAR_P ') then 
		  exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, VAR_P

	rewind 1000

! *******************************************************
	Do 
		read(1000,'(A80)')chline	   
		if(chline(1:5)=='T_OBS') then 
		  exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, T_OBS

	rewind 1000

! *******************************************************
Do 
   read(1000,'(A80)')chline
   if(chline(1:8)=='DATATIME') then 
      exit ! find the location
   endif    
Enddo ! 
i = 0
Do 
   read(1000,'(A80)')chline
   if(chline(1:1)=='/') then 
      exit ! find the location
   endif    
   i = i + 1
Enddo ! Finding number of times available 
N_DataType = i 
ALLOCATE(DATATIME(i))
rewind 1000
Do 
   read(1000,'(A80)')chline
   if(chline(1:8)=='DATATIME') then 
      exit ! find the location
   endif    
Enddo ! 
Read(1000,*)DATATIME(1:i)
rewind 1000 

 ! *******************************************************
Do 
   read(1000,'(A80)') chline
   if(chline(1:12)=='N_DATA_INDEX') then 
      exit ! find the location
    endif
Enddo ! 
!!
read(1000,*)  N_DATA_INDEX
Allocate(DATA_INDEX(Nwell, 5 , N_DATA_INDEX))
read(1000,*) 
Do j = 1 , N_DATA_INDEX
	Do i = 1 , Nwell
		read(1000,*) DATA_INDEX( i , 1 : 5 , j)
	EndDo
EndDo
rewind 1000
 ! ******************************************************* 
 	Do 
		read(1000,'(A80)')chline	   
		if(chline(1:12)=='MINERROR_WOR') then 
			exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, MINERROR_WOR

	rewind 1000
 ! ******************************************************* 
 	Do 
		read(1000,'(A80)')chline	   
		if(chline(1:12)=='MINERROR_RAT') then 
			exit ! find the location
		endif
	Enddo ! 
	backspace 1000	
	!!
	read(1000,*) ch1, MINERROR_RAT

	rewind 1000
 ! ******************************************************* 	


	close(1000)

!!***********************************************
! Reading observation Times

	Open(400,DEFAULTFILE='.\Code Input Files\',file='t_obs.dat',status='old')
	j = 0
	Do 
		read(400 , FMT = * , iostat = state ) i
		If(state/=0) Exit		
		j = j + 1
	Enddo
	rewind 400
	N_time_obs =  j 
	Allocate(time_obs( j))
	Do i = 1 , j
		read(400 , FMT = * , iostat = state ) time_obs(i)
	Enddo
	Close(400)






	
	Endsubroutine Read_OBSERVED_Input
	
	
End program main