subroutine output

use globalvar

implicit none

integer :: iWell, alfawell , j 
character(len=50)::charac,charac1
Real*8 ,Allocatable, Dimension(:,:) :: BHP_mat , flow_mat  , Oil_mat , Water_mat , wor_mat
integer ::  n_Flow , i1, i2   , i  , i_water , i_oil , i_wor , ii
integer :: N_report
Real*8 :: acc
!integer :: n_oil , n_water ,n_BHP 
!! The value of sigmad^2 of pressure is usually 2 psi ( from Input File)
!! The value of sigmad of Oil rate is 2.0 percent	  ( Form Input File)
!! The Value of sigmad of Water Rate is 3.0 percent ( Form Input File)

!! To save the predicted data, 

!! D_obs has 5 columns 
! 1st Column is Time Vector, 2nd column is Well # that the data is came from 
! 3rd Column is Data Type , 4th column is Variance  !! IF the data is not used in assimilation it would be saved as -1 (fake data )
! 5th column is the Data

! D_Pre has 1 columns which contains the Predicted Data
! The time and Well # and Variance and Type of data is the same as Observations

! Data type can be : 0 for water rate , 1 for Oil rate , 2 for BHP , 3 for Total Liquid Rate , and 4 for WOR ( or Water cut)

! Notice that the vector Diagonals referes to the diagonal entries of CD, Unless we run the true model to generate observed data
! We don't need this vector, however in case of generating predicted data with prior model or any model , this vector (diagonals ) shows 
! that which entry is a fake data ( we dont have any observation of that data ) and which entry is really a data

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		DATAT = 'RAT'
		!! If you want to have WOR as data instead of water rate,
		!! set DATAT='WOR'


		n_BHP = 0
		n_water = 0
		n_oil = 0
		n_wor = 0
		Do iWell=1,Nwell
			if ( W_type(iwell) == 'prod') then
				if (( Schedule_T(iwell , 2 , sche_N) .NE. 0 ) )then 
					n_BHP = n_BHP + 1	        					
					n_Oil = n_Oil + 1					
					n_water = n_water + 1
					n_wor = n_wor + 1
				EndIf
			Elseif ( W_type(iwell) == 'inje') then		
				if (( Schedule_T(iwell , 2 , sche_N) .NE. 0 ) )then 				
					n_BHP = n_BHP + 1
					n_water = n_water + 1 
				ENDif
			EndIF
		EndDO

		
		Allocate( BHP_mat(n_BHP,2) )		    
		Allocate( water_mat(n_water,2) )
		Allocate( oil_mat(n_oil,2) )
		Allocate( WOR_mat(n_wor,2) )		
		
		! 1st Entry of BHP_mat is the well Number ( that the data is came from )
		! 2nd Entry of BHP_mat is the value of the data
		
		i1 = 1				 
		i_Oil = 1
		i_water = 1
		i_wor = 1   
		Do iWell=1,Nwell
			if ( W_type(iwell) == 'prod') then
				if (Schedule_T(iwell , 2 , sche_N) .NE. 0 ) then ! water rate is a data							
					BHP_mat(i1,2) = pressure(Ntotal + iWell)
					BHP_mat(i1,1) = iWell
					i1 = i1 + 1
				
					oil_mat(i_oil,2) = oilrate(iWell)
					oil_mat(i_oil,1) = iWell
					i_oil = i_oil + 1					

					water_mat(i_water,1) = iWell		        										
					water_mat(i_water,2) = waterRate(iWell)		        					
					i_water = i_water + 1
					
					wor_mat(i_wor,1) = iWell
					wor_mat(i_wor,2) = wor(iWell)						
					i_wor = i_wor + 1					
				EndIf			
			Elseif ( W_type(iwell) == 'inje') then		
				if (Schedule_T(iwell , 2 , sche_N) .NE. 0 ) then ! water rate is a data										
					BHP_mat(i1,2) = pressure(Ntotal + iWell)
					BHP_mat(i1,1) = iWell
					i1 = i1 + 1				

					water_mat(i_water,2) = waterRate(iWell)		        
					water_mat(i_water,1) = iWell
					i_water = i_water + 1					        
				EndIf
			EndIF
		EndDO

			
		N_d = 0 
		If (DATAT == 'RAT') then
			N_d = N_oil + N_BHP + N_water
		Else
			N_d = N_oil + N_BHP + N_WOR
		EndIf
		
		Allocate( diagonals(N_d , 1) ) 
		Allocate( Well_No(n_d,1) )
		Allocate( Data_type(n_d,1))							
		Allocate(d_pre( N_d , 1 ) )


!!		FOPR = TOILP + FOPR ! FOPR is the total amount of oil it produced from time Zero
		if ( Gen_obs )Then  
			acc = 0
			do i = 1 , Nwell
				if ( W_type(i) .NE. 'inje') then
					acc = acc + waterrate(i)
				EndIf
			EndDo						
			FWPR = FWPR + acc*d_time
			FOPR = FOPR + SUM(oilrate)*d_time
			Open(2,DEFAULTFILE='.\sim_output\' , FILE='True_FOPR.DAT', status='unknown' , access = 'append')
			write(2,'(F10.1,2x , F15.1)') time , FOPR
			close(2)	

			Open(2,DEFAULTFILE='.\sim_output\' , FILE='True_FWPR.DAT', status='unknown' , access = 'append')
			write(2,'(F10.1,2x , F15.1)') time , FWPR
			close(2)	
		EndIf		

   		Do i = 1 , n_BHP
			d_pre(i ,1) = BHP_mat(i,2) 
			Well_No(i,1) = BHP_mat(i,1)  ! the well that the came from
			Data_type(i,1) = 2.0
			If ( DATA_INDEX( Well_No(i,1) , 2 , datasche_N) == -1 ) Then ! fake data 
				diagonals(i , 1) = -1				  !!! FOR CD 
			Else			
				diagonals(i , 1) = VAR_P			  !!! FOR CD 
			EndIf
		EndDo
		Do i =  1 ,n_oil
			j = i + n_BHP
			d_pre(j ,1) = oil_mat(i,2)
			Well_No(j,1) = oil_mat(i,1)  ! the well that the came from
			Data_type(j,1) = 1.0			
			If ( DATA_INDEX( Well_No(j,1) , 3 , datasche_N) == -1 ) Then ! fake data 
				diagonals(j , 1) = -1		
			Else	
				IF 	(d_pre(j ,1) .LE. 150) then								
					diagonals(j , 1) = ( SIGMA_QO * d_pre(j ,1)) ** 2 !! FOR CD
				Else
					diagonals(j , 1) = ( 3.0 ) ** 2	
				EndIF
				IF (diagonals(j , 1) < MINERROR_RAT) diagonals(j , 1) = MINERROR_RAT				
			EndIf			
		EndDo		
		
		If (DATAT == 'RAT') then !! We have the rate data
			Do i = 1 ,  n_water
				 j = i +  n_BHP + n_oil
				d_pre(j ,1) = Water_mat(i,2)
				Well_No(j,1) = water_mat(i,1)  ! the well that the came from
				Data_type(j,1) = 0.0						
				If ( DATA_INDEX( Well_No(j,1) , 4 , datasche_N) == -1 ) Then ! fake data 
					diagonals(j , 1) = -1		
				Else					
					IF 	(abs(d_pre(j ,1)) .LE. 150) then
						diagonals(j , 1) = ( SIGMA_QW * d_pre(j ,1)) ** 2
					Else
						diagonals(j , 1) = ( 3.0 ) ** 2
					EndIF
					IF (diagonals(j , 1) < MINERROR_RAT) diagonals(j , 1) = MINERROR_RAT								
				EndIf							
			EndDo
		Else						!! We have the WOR data
			Do i = 1 ,  n_wor
				j = i +  n_BHP + n_oil
				d_pre(j ,1) = WOR_mat(i,2)
				Well_No(j,1) = WOR_mat(i,1)  ! the well that the came from
				Data_type(j,1) = 4						
				If ( DATA_INDEX( Well_No(j,1) , 5 , datasche_N) == -1 )  Then ! fake data 
					diagonals(j , 1) = -1		
				Else						
					diagonals(j , 1) = ( SIGMA_WOR * d_pre(j ,1)) ** 2			
				EndIf											
			EndDo		
		EndIF
				
 		if (( abs(time) .GT. 1d-4 ) .and. (Gen_obs) )then   
			Open(2,DEFAULTFILE='.\sim_output\' , FILE='True_Data.DAT', status='unknown' , access = 'append')
			Do i = 1 , N_d
				If ( diagonals(i, 1) .NE. -1 )  then								
					write(2,'( (F15.1 ,2x), < 4 >(F15.5 ,2x))') time, well_No(i,1) , Data_type(i,1) , diagonals(i,1),d_pre(i,1)
				EndIF
			EndDo
			close(2)		
		Else 
			Do i = 1 , N_d
			   write(212,'( (F15.1 ,2x), < 4 >(F15.5 ,2x))') time, well_No(i,1) , Data_type(i,1) ,d_pre(i,1)
			EndDo
		EndIf	
		
		if ( ( .NOT.(Gen_obs) ) )Then  
			IF ( time == Nd_tobs(Obs_Index,1) ) then
						
				j = 0
				Do i = 1 , N_d
					If ( diagonals(i, 1) .NE. -1 )  then					
						write(215,'( (F15.1 ,2x), < 4 >(F15.5 ,2x))') time, well_No(i,1) , Data_type(i,1) ,d_pre(i,1)							
						pred_data( Nd_tobs(Obs_Index,3) + j ) = d_pre(i,1)
						j = j + 1					
					EndIF
				EndDo
				
				If (Forw_RUN ) then
					If (Obs_Index == 1) then
						Ens_FOPR(Obs_Index , j_ens) = Ens_FOPR(Obs_Index , j_ens) + sum(oilrate)*d_time
						acc = 0
						do i = 1 , Nwell
							if ( W_type(i) .NE. 'inje') then
								acc = acc + waterrate(i)
							EndIf
						EndDo
						Ens_FWPR(Obs_Index , j_ens) = Ens_FWPR(Obs_Index , j_ens) + acc*d_time
					Else
						Ens_FOPR(Obs_Index , j_ens) = Ens_FOPR(Obs_Index , j_ens)+ Ens_FOPR(Obs_Index - 1 , j_ens) + sum(oilrate) * d_time
						acc = 0
						do i = 1 , Nwell
							if ( W_type(i) .NE. 'inje') then
								acc = acc + waterrate(i)
							EndIf
						EndDo						
						Ens_FWPR(Obs_Index , j_ens) = Ens_FWPR(Obs_Index , j_ens) + Ens_FWPR(Obs_Index - 1 , j_ens)+ acc * d_time
					EndIF
				EndIF
				
				Obs_Index = Obs_Index + 1			
			Else
				If (Forw_RUN ) then
					Ens_FOPR(Obs_Index , j_ens) = Ens_FOPR(Obs_Index , j_ens) + sum(oilrate) * d_time
					acc = 0
					do i = 1 , Nwell
						if ( W_type(i) .NE. 'inje') then
							acc = acc + waterrate(i)
						EndIf
					EndDo						
					Ens_FWPR(Obs_Index , j_ens) = Ens_FWPR(Obs_Index , j_ens)+ acc * d_time
!!					EndIF
				EndIF			
			EndIF
		EndIf

	DeAllocate(BHP_mat , Oil_mat , water_mat, WOR_mat)		    
!********************************************************************************************************
!)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
!#########################################################################################################
	
	

Endsubroutine output