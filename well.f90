Module Well

Use globalvar
Use Trans

!************************************************************************************************************************************
!                                                    Explanations
! This Madul contains some Subroutines, the common point of these subroutines is that they are all about Wells, or Well rates

!Subroutine 'wellIndex_cons' calculates the matrix of constant part of Well Index matrix, 

!   subroutine wellIndex , multiplies the constant term of Well Index matrix by (kro)/(muo)(Bo) and (krw)/(muw)(Bw) , 
! Needless to say, the Well Index is a matrix with Nwell columns and Nz rows, so for each layer that is not ptoducing 'oil' or 'water ' , 
! the entry of W_Index_O or W_Index_ w will be zero respectively. It is possible that a layer just produces oil or just Water, then kro or
! krw would be zero

! subroutine calc_rate, computes the oil rate and water rate of each well, and puts the value into a vector, the vectors Oilrate and Waterrate are
! Global vectors and are used in other subroutines like 'con_check' and Output

! In subroutine Calc_DeltaPwf, the value of DeltaPwf is calculated for each layer in the well and these values are in a matrix with Nwell columns and Nz rows
! The algorithm used in this sub is from Dr. Li's course Notes . There are 4 cases considered here : Production well, with q oil sp or BHP sp , and 
! Injector Well with q water sp, or BHP sp. 
! in case of constant rate production or injection, the value of BHP is determined using a loop. This loop is very exact and after BHP is determined, there 
! is no need to improve it, as Fwell (would become)= 1d-12 
! In case of BHP sp, the iterations aim to precisely calculate del_pwf and the convergence criteria is that oil rate or water rate is not changing with
! recalculating del_pwf

!************************************************************************************************************************************
! contains transmissibility functions
use globalvar
implicit none

Contains 

    subroutine wellIndex_cons ! calculates the matrix of constant part of Well Index
    implicit none
    integer:: i_well, i_k,  alfa_well , counter	, m , n



        Do i_well = 1 , Nwell
            Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)
                alfa_well = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1 ) * Nx + (i_k - 1 ) * Nx * Ny
                R_equ(i_well , i_k ) = 0.28 * sqrt( ky(alfa_well) * Deltax( Well_loc(i_well,1))**2 + kx(alfa_well) * Deltay( Well_loc(i_well,2))**2 ) / ( sqrt(kx(alfa_well)) + sqrt(ky(alfa_well)) )
                KH(i_well , i_k) =  ( kx(alfa_well) * ky(alfa_well) ) ** (0.5)
                ! the above terms are not a function of pressure or Saturation, but the below one is, so I can calculate the above terms for once and then just
                ! calculate the pressure and saturation dependent term during program Run
                Well_index_cons(i_well , i_k ) = 2 *  3.14159265 * betac * KH(i_well , i_k) * Deltaz(alfa_well) / log( R_equ(i_well,i_k) / r_w(i_well) )
            EndDo					
        EndDo
            n = n_Kr
            m = 1
            do while (abs(n-m)>1)
                if (SW_AtSor < Sw_mat(floor( m + (n-m+1)/2.0) ) )  then
                    n = floor(m + ( n - m + 1) / 2.0);
                else
                    m = floor(m + ( n - m + 1) / 2.0);
                end if
            end do
        
		kr_inj = krw(m)+ (SW_AtSor -Sw_mat(m))/(Sw_mat(n)-Sw_mat(m))*(krw(n)-krw(m))        
    Endsubroutine wellIndex_cons
    
    !********************************************************************************************************
    subroutine wellIndex ! this sub would be called at each time step from solver subroutine
    integer:: i_well, i_k,  alfa_well , counter
        !counter = 0

        W_index_o = 0.0
        W_index_w = 0.0
        Do i_well = 1 , Nwell
			if ( W_type(i_well) == 'inje') then
				Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)         
					alfa_well = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1 ) * Nx +  ( i_k - 1) * Nx * Ny
					W_index_o(i_well,i_k) = 0.0 
					if ( k_r(alfa_well , 'w' ) .GT. kr_inj )  then 
		                W_index_w(i_well,i_k) = Well_index_cons(i_well,i_k) * k_r(alfa_well , 'w' ) / (mu_f( pressure(alfa_well),'w' ) * B_f( W_press(alfa_well) , 'w' ) )					
					Else 
						W_index_w(i_well,i_k) = Well_index_cons(i_well,i_k) * kr_inj / (mu_f( pressure(alfa_well),'w' ) * B_f( W_press(alfa_well) , 'w' ) )
					EndIf					
				EndDo			
			Else			
				Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)         
					alfa_well = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1 ) * Nx +  ( i_k - 1) * Nx * Ny
					W_index_o(i_well,i_k) = Well_index_cons(i_well,i_k) * k_r(alfa_well , 'o' ) / (mu_f( pressure(alfa_well),'o' ) * B_f( pressure(alfa_well) , 'o' ) )
					W_index_w(i_well,i_k) = Well_index_cons(i_well,i_k) * k_r(alfa_well , 'w' ) / (mu_f( pressure(alfa_well),'w' ) * B_f( W_press(alfa_well) , 'w' ) )
				EndDo
			EndIF
        EndDO
        
    Endsubroutine wellIndex
    
    !********************************************************************************************************
    subroutine calc_rate
    
    Implicit none
    
    integer :: i_well , alfa_w  , alfa_wprim , alfa_wpp , i_k ,i , W_index
    integer :: counter
    real*8 :: qo, qw 

    ! writing oil rates to a matrix 
    oilrate = 0.0
    waterrate = 0.0
    TOILP = 0
    WOR = 0
        
    If (Nz .NE. 1) then
        call Calc_DeltaPwf
	EndIf        
        counter = 0 
        WOR = 0
        Do i_well = 1, Nwell
            if ( ( W_type(i_well) == 'prod') .AND. ( Schedule(i_well , 2 , sche_N) .NE. 0)  )then
				qo = 0
				qw = 0
				W_index = Ntotal + i_well  ! Pwf position in pressure vector
				Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2) 
					 alfa_w = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx +  ( i_k - 1) * Nx * Ny
					 qo = qo + w_index_o(i_well, i_k ) * (pressure(alfa_w) - pressure(W_index ) - Del_Pwf(i_well, i_k))
					 qoL(i_well, i_k) = w_index_o(i_well, i_k ) * (pressure(alfa_w) - pressure(W_index ) - Del_Pwf(i_well, i_k))
				EndDO

				Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2) 				
					alfa_w = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx +  ( i_k - 1) * Nx * Ny
					qw = qw + w_index_w(i_well, i_k ) * (pressure(alfa_w) - pressure(W_index ) - Del_Pwf(i_well, i_k))                 
					qwL(i_well, i_k) = w_index_w(i_well, i_k ) * (pressure(alfa_w) - pressure(W_index ) - Del_Pwf(i_well, i_k))
			    EndDo

				oilrate(i_well) = qo
				waterRate(i_well) = qw
				if ( abs(qo) .GE. 1 ) then
					WOR(i_well) = qw / qo
				EndIf
				TOILP = TOILP + oilrate(i_well) * d_time !! Total OIL Production in that Time_step				
!				if ( Schedule(i_well , 3 ,sche_N) == 0.0 ) then
!				   oilrate (i_well) = 0
!				   waterRate(i_well) = 0 
!				EndIf
			ElseIf ( Schedule(i_well , 2 , sche_N) .NE. 0) Then 
				oilrate (i_well) = 0
				qw = 0
				W_index = Ntotal + i_well  ! Pwf position in pressure vector
				Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2) 
					 alfa_w = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx +  ( i_k - 1) * Nx * Ny
					 qw = qw + w_index_w(i_well, i_k ) * (pressure(alfa_w) - pressure(W_index ) - Del_Pwf(i_well, i_k))                 
				EndDO				
				waterRate(i_well) = qw				
			EndIf
        EndDo
    
    Endsubroutine calc_rate
    
    !********************************************************************************************************
        
    subroutine Calc_DeltaPwf ! calculates the matrix of Del_Pwf
    implicit none
    
    integer :: i_well, i_k,  alfa_well , alfa_w , W_index , B_index , N_index , ii , counter
    integer :: con ! convergence , when == 0 , not converged 
    Real, Allocatable, Dimension(:) :: O_Volume, W_Volume , mix_den ! O_Volume is the vector of Volume of oil in each completion ( for dt = 1 day , and dp = 1 psi )
    Real, Allocatable, Dimension(:) :: qoN , qoP , qwN , qwP ! q oil New, q water New ... for convergence check in case of BHPsp
    real*8 :: denom , eps_pwf , Fwell ,max_dif , c1 , c2
    
    
    Del_Pwf = 0
    Allocate ( O_Volume(Nz) )
    Allocate ( W_Volume(Nz ) ) ! Well_k(i_well ,2) - Well_k(i_well ,1)
    Allocate ( mix_den(Nz ) )                              
    
    Allocate ( qoN(Nz) )
    Allocate ( qwN(Nz) )
    Allocate ( qoP(Nz) )
    Allocate ( qwP(Nz) )


        Do i_well = 1 , Nwell       
            W_index = Ntotal + i_well  ! Pwf position in pressure vector
            alfa_w = 2 * Ntotal + i_well 
            if ( W_type(i_well) == 'prod') then
               if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! oil rate is given
                       counter  = 0 
                       con = 0 ! convergence criteria
                       Do while (con == 0) 
                           if ( counter .NE. 0 ) then
                              ! checking convergence   
                              max_dif = 0
                              Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                       
                                 ! finding the maximum difference
                                 B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                   
                                 max_dif = max ( max_dif , abs( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) ) )
                              EndDo
                              if  ( (max_dif < 1.0D-7 ) .OR.  ( abs(eps_pwf / max_dif) < 0.01 ) )  then 
                                 con = 1
                                 exit
                              EndIf
                           EndIf
                           counter = counter + 1 
                           O_Volume = 0
                           W_Volume = 0
                           mix_den = 0
                           i_K =  Well_k(i_well ,2)  
                           ! stage 1 : obtaining Oil Volume and Water Volume , 
                           Do while (i_k .GE. Well_k(i_well , 1) )              
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                              O_Volume(i_k) = w_index_o(i_well, i_k ) * 1.0 * ( pressure(B_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
                              W_Volume(i_k) = w_index_w(i_well, i_k ) * 1.0 * ( pressure(B_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
                              if (i_k .NE. Well_k(i_well ,2) ) then
                                 O_Volume(i_k) = O_Volume(i_k) + O_Volume(i_k + 1)
                                 W_Volume(i_k) = W_Volume(i_k) + W_Volume(i_k + 1)
                              EndIf
                              mix_den(i_k) = ( W_Volume(i_k) * den_w + O_Volume(i_k) * den_o ) / ( W_Volume(i_k) * den_w / density( pressure(W_index) + Del_Pwf(i_well, i_k),'w') + O_Volume(i_k) * den_o / density( pressure(W_index) + Del_Pwf(i_well, i_k),'o' ) ) 
                              ! BHP or Pwf for oil and water are the same , so pressure( W_index) works for both
                              i_K = i_K - 1 
                           EndDo
                           ! stage 2 : obtaining Del_pwf
                           Do i_k = Well_k(i_well , 1) + 1 , Well_k(i_well ,2)                             
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                              N_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 - 1) * Nx * Ny
                              Del_Pwf(i_well , i_k ) = mix_den(i_k - 1) / 144.0 *( Deltaz(B_index) +  Deltaz(N_index) ) / 2.0
                              !if (i_k .NE. Well_k(i_well , 1) + 1 ) then ! if this is not the first layer which has del_pwf
                              Del_Pwf(i_well , i_k ) = Del_Pwf(i_well , i_k ) + Del_Pwf(i_well , i_k - 1)
                           EndDO                  
                           ! stage 3 : modifying pwf
                           Fwell = - Schedule(i_well , 3 , sche_N)
                           denom = 0
                           Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                             
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                
                              Fwell = Fwell + w_index_o(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k))
                              denom = denom + w_index_o(i_well, i_k )
                           EndDO
!                           if ( Schedule(i_well , 3 , sche_N) .EQ. 0.0 ) then ! total liquid rate == 0
!						       Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                        
!                                  B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                
!								  Fwell = Fwell + w_index_w(i_well, i_k ) * ( W_press(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k))
!                                  denom = denom + w_index_w(i_well, i_k )
!                               EndDo
!                           EndIf
                           eps_pwf =  Fwell / denom
                           pressure(W_index ) = pressure(W_index ) + eps_pwf
                           YY(alfa_w) = YY(alfa_w) + eps_pwf
                       EndDO
               Elseif ( Schedule(i_well , 2 , sche_N) .EQ. 3 ) then ! Total Liquid Rate Specified
                       counter  = 0 
                       con = 0 ! convergence criteria
                       Do while (con == 0) 
                           if ( counter .NE. 0 ) then
                              ! checking convergence   
                              max_dif = 0
                              Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                       
                                 ! finding the maximum difference
                                 B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                   
                                 max_dif = max ( max_dif , abs( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) ) )
                              EndDo
                              if  ( (max_dif < 1.0D-7 ) .OR. abs(eps_pwf / max_dif) < 0.01)  then 
                                 con = 1
                                 exit
                              EndIf
                           EndIf
                           counter = counter + 1 
                           O_Volume = 0
                           W_Volume = 0
                           mix_den = 0
                           i_K =  Well_k(i_well ,2)  
                           ! stage 1 : obtaining Oil Volume and Water Volume , 
                           Do while (i_k .GE. Well_k(i_well , 1) )              
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                              O_Volume(i_k) = w_index_o(i_well, i_k ) * 1.0 * ( pressure(B_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
                              W_Volume(i_k) = w_index_w(i_well, i_k ) * 1.0 * ( pressure(B_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
                              if (i_k .NE. Well_k(i_well ,2) ) then
                                 O_Volume(i_k) = O_Volume(i_k) + O_Volume(i_k + 1)
                                 W_Volume(i_k) = W_Volume(i_k) + W_Volume(i_k + 1)
                              EndIf
                              mix_den(i_k) = ( W_Volume(i_k) * den_w + O_Volume(i_k) * den_o ) / ( W_Volume(i_k) * den_w / density( pressure(W_index) + Del_Pwf(i_well, i_k),'w') + O_Volume(i_k) * den_o / density( pressure(W_index) + Del_Pwf(i_well, i_k),'o' ) ) 
                              ! BHP or Pwf for oil and water are the same , so pressure( W_index) works for both
                              i_K = i_K - 1 
                           EndDo
                           ! stage 2 : obtaining Del_pwf
                           Do i_k = Well_k(i_well , 1) + 1 , Well_k(i_well ,2)                             
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                              N_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 - 1) * Nx * Ny
                              Del_Pwf(i_well , i_k ) = mix_den(i_k - 1) / 144.0 *( Deltaz(B_index) +  Deltaz(N_index) ) / 2.0
                              !if (i_k .NE. Well_k(i_well , 1) + 1 ) then ! if this is not the first layer which has del_pwf
                              Del_Pwf(i_well , i_k ) = Del_Pwf(i_well , i_k ) + Del_Pwf(i_well , i_k - 1)
                           EndDO                  
                           ! stage 3 : modifying pwf
                           Fwell = - Schedule(i_well , 3 , sche_N)
                           denom = 0
                           Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                             
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                
                              Fwell = Fwell + w_index_o(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k)) & 
											+ w_index_w(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k))
                              denom = denom + w_index_o(i_well, i_k )  + w_index_w(i_well, i_k ) 
                           EndDO
                           eps_pwf =  Fwell / denom
                           pressure(W_index ) = pressure(W_index ) + eps_pwf
                           YY(alfa_w) = YY(alfa_w) + eps_pwf
                       EndDO
                       
               ElseIf ( Schedule(i_well , 2 , sche_N) .EQ. 2 ) then ! BHP specified case 
                       qoN = 0
                       qoP = 0
                       qwN = 0
                       qwP = 0
                       counter  = 0 
                       con = 0 ! convergence criteria
                       Do while (con == 0) 
                           counter = counter + 1 
                           O_Volume = 0
                           W_Volume = 0
                           mix_den = 0
                           i_K =  Well_k(i_well ,2)  
                           ! stage 1 : obtaining Oil Volume and Water Volume , 
                           Do while (i_k .GE. Well_k(i_well , 1) )              
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                              O_Volume(i_k) = w_index_o(i_well, i_k ) * 1.0 * ( pressure(B_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
                              W_Volume(i_k) = w_index_w(i_well, i_k ) * 1.0 * ( pressure(B_index) - pressure(W_index) - Del_Pwf(i_well, i_k) )
                              if (i_k .NE. Well_k(i_well ,2) ) then
                                 O_Volume(i_k) = O_Volume(i_k) + O_Volume(i_k + 1)
                                 W_Volume(i_k) = W_Volume(i_k) + W_Volume(i_k + 1)
                              EndIf
                              mix_den(i_k) = ( W_Volume(i_k) * den_w + O_Volume(i_k) * den_o ) / ( W_Volume(i_k) * den_w / density( pressure(W_index) + Del_Pwf(i_well, i_k),'w') + O_Volume(i_k) * den_o / density( pressure(W_index) + Del_Pwf(i_well, i_k),'o' ) ) 
                              ! BHP or Pwf for oil and water are the same , so pressure( W_index) works for both
                              i_K = i_K - 1 
                           EndDo
                           ! stage 2 : obtaining Del_pwf
                           Do i_k = Well_k(i_well , 1) + 1 , Well_k(i_well ,2)                             
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                              N_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 - 1) * Nx * Ny
                              Del_Pwf(i_well , i_k ) = mix_den(i_k - 1) / 144.0 *( Deltaz(B_index) +  Deltaz(N_index) ) / 2.0
                              !if (i_k .NE. Well_k(i_well , 1) + 1 ) then ! if this is not the first layer which has del_pwf
                              Del_Pwf(i_well , i_k ) = Del_Pwf(i_well , i_k ) + Del_Pwf(i_well , i_k - 1)
                           EndDO                  
                           ! stage 3 : convergence 
                            Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                       
                               B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                                             
                               qoN(i_k) = w_index_o(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) )
                               qwN(i_k) = w_index_w(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) )
                            EndDo
                           if ( counter .NE. 1 ) then
                              ! checking convergence 
                              max_dif = 0
                              Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                                                        
                                 if  ( abs(qwN(i_k)) > 1   ) then
                                    max_dif = max ( max_dif , abs(qwN(i_k) - qwP(i_k) ) /qwN(i_k) )
                                 EndIf
                                 if  ( abs(qoN(i_k)) > 1   ) then
									max_dif = max ( max_dif , abs(qoN(i_k) - qoP(i_k) ) /qoN(i_k) )                             
                                 EndIf
							  EndDo
                              if  (max_dif < 0.000001)  then 
                                 con = 1
                                   !exit
                              EndIf
                           EndIf
                          qoP = qoN
                          qwP = qwN                        

                       EndDO               
                   EndIf
            Else ! Injection Well
               if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! injection rate is given
                       counter  = 0 
                       con = 0 ! convergence criteria
                       Do while (con == 0) 
                           if ( counter .NE. 0 ) then
                              ! checking convergence   
                              max_dif = 0
                              Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                       
                                 ! finding the maximum difference
                                 B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                   
                                 max_dif = max ( max_dif , abs( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) ) )
                              EndDo
                              if  (abs(eps_pwf / max_dif) < 0.01)  then 
                                 con = 1
                                 exit
                              EndIf
                           EndIf
                           counter = counter + 1 

                           ! obtaining Del_pwf
                           Do i_k = Well_k(i_well , 1) + 1 , Well_k(i_well ,2)                             
                              c1 = 0
                              c2 = 2
                              Do while ( abs(c1-c2) > 0.1) 
                                  c1 = Del_Pwf(i_well , i_k )
                                  B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                                  N_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 - 1) * Nx * Ny
                                  Del_Pwf(i_well , i_k ) = density( pressure(W_index) + Del_Pwf(i_well, i_k),'w') / 144.0 *( Deltaz(B_index) +  Deltaz(N_index) ) / 2.0
                                  !if (i_k .NE. Well_k(i_well , 1) + 1 ) then ! if this is not the first layer which has del_pwf
                                  Del_Pwf(i_well , i_k ) = Del_Pwf(i_well , i_k ) + Del_Pwf(i_well , i_k - 1)
                                  c2 = Del_Pwf(i_well , i_k )
                              EndDo
                           EndDO
                           ! stage 3 : modifying pwf
                           Fwell = Schedule(i_well , 3 , sche_N)
                           denom = 0
                           Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                             
                              B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                
                              Fwell = Fwell + w_index_w(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k))
                              denom = denom + w_index_w(i_well, i_k )
                           EndDO
                           
                           eps_pwf =  Fwell / denom
                           pressure(W_index ) = pressure(W_index ) + eps_pwf
                           YY(alfa_w) = YY(alfa_w) + eps_pwf
                       EndDO
               Elseif ( Schedule(i_well , 2 , sche_N) .EQ. 2 ) then  ! BHP of Injection specified case 
                       qwN = 0
                       qwP = 0
                       counter  = 0 
                       con = 0 ! convergence criteria
                       Do while (con == 0) 
                           counter = counter + 1 
                           O_Volume = 0
                           W_Volume = 0
                           mix_den = 0
                           i_K =  Well_k(i_well ,2)  
                           ! stage 2 : obtaining Del_pwf
                           Do i_k = Well_k(i_well , 1) + 1 , Well_k(i_well ,2)                             
                              c1 = 0
                              c2 = 2
                              Do while ( abs(c1-c2) > 0.1) 
                                  c1 = Del_Pwf(i_well , i_k )
                                  B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 ) * Nx * Ny
                                  N_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 - 1) * Nx * Ny
                                  Del_Pwf(i_well , i_k ) = density( pressure(W_index) + Del_Pwf(i_well, i_k),'w') / 144.0 *( Deltaz(B_index) +  Deltaz(N_index) ) / 2.0
                                  !if (i_k .NE. Well_k(i_well , 1) + 1 ) then ! if this is not the first layer which has del_pwf
                                  Del_Pwf(i_well , i_k ) = Del_Pwf(i_well , i_k ) + Del_Pwf(i_well , i_k - 1)
                                  c2 = Del_Pwf(i_well , i_k )
                              EndDo
                           EndDO                           ! stage 3 : convergence 
                           if ( counter .NE. 1 ) then
                              ! checking convergence 
                              max_dif = 0
                              Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                       
                                 B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                                             
                                 qwN(i_k) = w_index_w(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) )
                                 max_dif = max ( max_dif , abs(qwN(i_k) - qwP(i_k) ) /qwN(i_k) )
                              EndDo
                              if  (max_dif < 0.01)  then 
                                 con = 1
                                 !exit
                              EndIf
                           Else
                               Do i_k = Well_k(i_well , 1) , Well_k(i_well ,2)                                                       
 								   B_index = Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx + ( i_k - 1 )* Nx * Ny                                                             
								   qoN(i_k) = w_index_o(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) )
								   qwN(i_k) = w_index_w(i_well, i_k ) * ( pressure(B_index) - pressure(W_index ) - Del_Pwf(i_well, i_k) )
							   EndDo
                           EndIf
                           qoP = qoN
                           qwP = qwN                                                       
                       EndDO                             
            EndIf
       EndIf
          
   EndDo
   DeAllocate (O_Volume)        
   DeAllocate (W_Volume)        
   DeAllocate (mix_den)        
   
  Endsubroutine Calc_DeltaPwf 
    !********************************************************************************************************
!    subroutine calc_fwell
!    
!    Endsubroutine calc_fwell
    !********************************************************************************************************
subroutine check_BHP

integer :: i_well , W_index, alfa_w

integer :: iwell_sche ! Number of well that its schedule has changed
integer :: sche_type , b_index ! type of new schedule
real*8 :: sche_val ! value of new schdule
Logical :: change_Sche ! If a schdule has changed , it would be .True.


change_Sche = .False.
!	If (time == 500) then 
!			time = 500
!	endIf
	Do i_well = 1 , Nwell       
		b_index =  Well_loc(i_well,1) + ( Well_loc(i_well,2) - 1) * Nx 
		change_Sche = .False.
		W_index = Ntotal + i_well  ! Pwf position in pressure vector
		alfa_w = 2 * Ntotal + i_well 
		if ( W_type(i_well) == 'prod') then
            if ( WOR(i_well) > max_WOR) then ! shut in the well
				change_Sche = .True.
			    Schedule(i_well , 2 , sche_N) = 1 ! switching to zero rate
		        Schedule(i_well , 3 , sche_N) = 0.0
		        iwell_sche = i_well ; sche_type = 1  ; sche_val = 0.0 ; change_Sche = .True. ; exit
		    EndIf		  		
		    if ((( Schedule(i_well , 2 , sche_N) .EQ. 1 ) .OR. ( Schedule(i_well , 2 , sche_N) .EQ. 3 ) ) .AND. ( Schedule(i_well , 3 , sche_N) .NE. 0 )  )then ! oil rate is given
			   if ( pressure(W_index) < minBHP) then
					change_Sche = .True.
					Schedule(i_well , 2 , sche_N) = 2 ! switching to constant BHP
					Schedule(i_well , 3 , sche_N) = minBHP
					iwell_sche = i_well ; sche_type = 2  ; sche_val = minBHP ; change_Sche = .True. ; exit
			   EndIF
			EndIf		   
		    
		Else
			if ( Schedule(i_well , 2 , sche_N) .EQ. 1 ) then ! injection rate is given		    
				if ( pressure(W_index) > maxBHP) then
					change_Sche = .True.
					Schedule(i_well , 2 , sche_N) = 2 ! switching to constant BHP
					Schedule(i_well , 3 , sche_N) = maxBHP
					iwell_sche = i_well ; sche_type = 2  ; sche_val = minBHP ; change_Sche = .True. ; exit			     			     
				EndIF
			EndIf
		EndIF
		
		If (change_Sche) then
			Open(2,DEFAULTFILE='.\sim_output\' , FILE='WARNING_Change of Schedule.DAT', status='unknown' , access = 'append')		        
			write(2,*) 'Time =' , Time
			write(2,*) 'well #  ', iwell_sche
			write(2,*) 'New Schedule Type', Schedule(iwell_sche , 2 , sche_N)
			write(2,*) 'New Schedule Value', Schedule(iwell_sche, 3 , sche_N)				
			close(2)
		EndIf		
	EndDo


Endsubroutine check_BHP    
    
    !********************************************************************************************************     
    
EndModule Well