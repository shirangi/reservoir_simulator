!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!! 3D oil reservoir Simulator project 
!! by Mehrdad Gharib Shirangi From Spring 2010 until Fall 2010
!! mehrdad.ghsh@gmail.com
!! Copy Right by Mehrdad Gharib Shirangi
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

module globalvar

Integer:: Nx , Ny , Nz, Ntotal, Nwell, &
         n_Kr  !Number of rows in rel. perm table
         
Integer::N_mup , & ! simulation time
          N_Bo ,&!N_mup is the number of data of oil viscosity as a function of pressure
         N_Dt ,& ! number of time intervals at which we want to have data
         NSchedule , time_index , &
         N_equ , &! number of equations
         size_table , i_iter

integer:: sche_N , &! Sche_N is the schedule number    
          M_ITER ! =1 if converged , and = 0 if not converged       
          
Logical :: Converged 

real*8 :: den_o , den_w  , max_WOR ! oil density and water density 

real :: time , OWC_depth , OWC_PCOW , Datum_D , Datum_P
real*8, parameter :: betac = 0.001127	 , pi =3.14159265
real*8, parameter :: alfac = 5.615 , c_fluid = 4.255d-5 , Epsilon = 0.0001 ! the value of this Epsilon is used for the saturation of the grid blocks in Aquifer
                                                                         ! it is used in the fucntion of PCO_inv

real*8 :: C_t , d_time , C_rock , P_ref , MinBHP , MaxBHP , OMB , MOMB , MRHS , SOR , SW_AtSor, kr_inj	, Initial_P , SWI  , time_step

Character*3 :: USEIP
Logical :: USEISW

!Real*8,Allocatable, Dimension(:,:) :: Schedule
Real*8,Allocatable, Dimension(:,:) ::  mu_p, Bo_p, KH , Del_pwf, qoL, qwL ! qoL is oil rate of each layer 
Real*8 ,Allocatable, Dimension(:) :: Sw_mat , Op_table, Wp_table , Sw_table , Z_table  , WOR
Real*8 ,Allocatable, Dimension(:) :: kx, ky , kz , phi, Deltax, Deltay, Deltaz, tops ,& ! PVTW  has 4 data: Ref. Pressure, Water FVF , Water Compressibility, Water viscosity at P ref 
                                    pvtw , & !C_R is the rock compressibility which has two values, first reference Pressure and 2nd is compressibility   
                                    Deltat , Cu_t, & ! Cu_t is culmulative time, which we want to get data
                                    Zblock, & ! Zblock is block center depth , used to calculate potential difference
                                    krw , kro , Pcow ! the entry in each of the vectors Krw , Kro and Pcow refer to the same entry in Sw vector
Real*8 , Allocatable, Dimension(:) :: initialp, &  ! initialp is a 2 entry vector which the first entry is the reference pressure and the second entry is the reference depth
                                    pressure , W_press, Past_pw, new_p, past_p, & ! past_p is the pressure vector at the previous time step
                                    YY,YY_Past, Sw, past_Sw, FCN , del_y, Sw_exact ! pressure is the pressure vector at the previous iteration
                                      ! and new p is the pressure vector at the new time step
                                    ! YY is the vector including all the unknowns ( Sw, Po, Pwf )
                                    ! del_y is the solution to J.del_y =-FCN
                                    ! Sw_exact is the actual solution, in order to make the solver work, we must replace zero values of saturation by 0.001 or less in YY

Real*8,Allocatable, Dimension(:) :: r_w , &! well radius
                                    t_Schedule, &
                                    Oilrate , waterRate , k_ro, k_rw , B_o , B_w , mu_o , mu_w , DkrodS , DkrwdS
Real*8 :: TOILP ! TOILP is the total Oil production
                                   
Real*8 ,Allocatable, Dimension(:) :: perm_x , perm_y , perm_z   

Real*8 ,Allocatable, Dimension(:,:) :: Well_index_cons  , &! permeability of each block, perm (Ntotal,1:3), which 1 refers to x, 2 refers to y, and 3 refers to z
                                      W_index_o , W_index_w , & ! well index for oil flow; well index for water flow
                                      R_equ !  R_equ is a vector which contains  the equivalent radius of each of the wells
                                      
Real,Allocatable, Dimension(:,:,:) :: Schedule  , Schedule_T

Integer , Allocatable, Dimension(:,:) ::well_loc, well_k


Character*6 , Allocatable, Dimension(:) :: well_name
Character*4 , Allocatable, Dimension(:) :: W_type ! W_type(i) is either 'inje' or 'prod'

integer :: nz_num1
 
Real*8, Allocatable, Dimension(:) :: a1

integer*4 ,Allocatable, Dimension(:) :: ia1, ja1

integer ,Allocatable, Dimension(:) :: maxentry

real*8 :: maxFCN

real*8 :: CU_TIME , RUN_TIME,  time_begin

real*8 :: SIGMA_QO , SIGMA_QW  , SIGMA_WOR , VAR_P , T_OBS

integer :: N_OBS , Obs_Index
Logical :: Gen_Ensem , COV_LOC   , GEN_OBS
Real*8 :: MINERROR_RAT  

Real*8 :: FOPR 



Real*8, Allocatable, Dimension(:) :: DATATIME , pred_data
Real*8, Allocatable, Dimension(:,:) :: Phi_prior 
integer :: datasche_N

Integer :: N_DATA_INDEX
Real*8, Allocatable, Dimension(:,:,:)	:: DATA_INDEX
Real*8 , Allocatable, Dimension(:,:) ::  LCM 
Integer :: Ne

!! Adjoint  Varibales

Logical :: Adjoint_Sol
Integer :: N_L , & ! N_L is the number of simulation time_steps
			N_m    ! N_m is the number of model parameters 

Integer :: l_adj , N_adjoint , N_sim , N_Forward ,N_lanczos
Real*8 :: d_tn_Plus !! Is the time_step from n+1 to n  while d_time is the time_step from n to n-1
Real*8 , Allocatable, Dimension(:) :: sim_time , aCM  , Mean_Dpr !!model_par ! df_dyn 
Integer , Allocatable, Dimension(:) :: iCM, jCM

!!From EnKFglobalvar

integer ::  N_d , j_ens
		!		   N_d is the number of data ( observed data)		
real :: T_End ! T_end is the end of assimulation

Real*8 ,Allocatable, Dimension(:,:) :: Ln_K  , m_uc  , m_true , m_p , d_pre , d_obs , diagonals ,OWC_Ens ! predicted data vector

Real*8 ,Allocatable, Dimension(:,:) :: well_No, Data_type , d_observ

Integer :: N_t_Obs 
Integer ,Allocatable, Dimension(:,:) :: Nd_tobs	! in Nd_tobs 1st column is the time of observation and 2nd column is the number of data in that time

Logical ::EST_PHI , EST_OWC , EST_RP 

Real*8 ,Allocatable, Dimension(:,:) :: true_perm , true_por , phi_uc , phi_true , phi_p , FOPR_prior

integer :: n_oil , n_water ,n_BHP , n_wor
integer :: N_par ! Number of model parameters , depending if we wanna update porosity and k_z too or not
				 ! The nest plan is to put all model parameters into a state vector and update them all together
Character*3 :: DATAT !! can be RAT (== rate) or WOR

Integer ,Allocatable, Dimension(:,:) :: WELLDATA 
real*8 :: MINERROR_WOR


Logical ::  Forw_RUN
Integer ::  N_DataType
Integer :: N_time_obs
Real*8 ,Allocatable, Dimension(:) ::  time_obs

Real*8 ,Allocatable, Dimension(:,:) ::  Ens_FOPR ,   Ens_FWPR

Real*8 :: FWPR


Real*8 ,Allocatable, Dimension(:) ::D_CM


endmodule globalvar