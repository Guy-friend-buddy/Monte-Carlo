
 PROGRAM MAIN

 USE nrtype  ! the first three are modules from book nr recipies modified 
		! the modification includes swaping single precision to double precision
		! therefore in the code SP mean double precision and dp means single
		! this was done since most of the code in the book is written for SP
 USE nr
 USE moduleglobal
 USE inversecdf  ! the module is used to define the cdf functions for the variable y and k

 IMPLICIT NONE

    integer, parameter :: ndim = 3,no_blocks = itermain/1000.0     ! ndim = dimension of integral
    integer(kind =8) :: temp,ind,catch_y,catch_yprime ! integers used to find index 
    REAL(kind =8), DIMENSION(ndim) :: rand ! rand is the ndimensional ( ndim)  random number
    real(kind=8) :: t1, t2, elapsed_time ! hunt is temp variable to find the index of the wallnormal/delst
    integer(kind=8) :: tclock1, tclock2, clock_rate   ! These variables are used to estimate the time program takes for completion
    real(kind =8), dimension(no_blocks,number_points) :: sum_1000,error
    real(kind =8), dimension(number_points) :: ARGX1   
    real(kind=8),dimension(number_points) :: SUM_total, SQR_DEVIATION,SUMT,SUMIN,SUMMD,SUMOT
    real(kind = 8),dimension(number_points) :: SUMI1_ININ,SUMI2_ININ,SUMI1_INMD,SUMI2_INMD, &
    & SUMI1_INOT,SUMI2_INOT,SUMI1_MDIN,SUMI2_MDIN,SUMI1_MDMD, &
    & SUMI2_MDMD, SUMI1_MDOT, SUMI2_MDOT,SUMI1_OTIN,SUMI2_OTIN, &
    & SUMI1_OTMD,SUMI2_OTMD,SUMI1_OTOT,SUMI2_OTOT    ! these are the values of integrals broken into different stratas across BL, however due to complexity are not used, but will be included in future. 
    real(kind =8) :: BRKSQ,k_bar,x1,x1a,C1,C2,C_mean    ! different variables like k^2, random number for k3, correlation length etc
    real(kind =8) :: CIN,EXPIN,VOLIN,X2AI,XF2I,ARG2I,X2I ! different varibles used for variance reduction
    integer :: I,valid_evaluations,J,step,erc     ! iteration variables
    real(kind =8) :: CMD,XPMD,EXPMD,VOLMD,X2AM,X2M,XF2M,ARG2M   ! variables of program used for variance reduction in the mid zone of BL
    real(kind=8) :: COT,XPOT,EXPOT,VOLOT,X2AT,X2T,XF2T,ARG2T   ! varibales used for variance reduction in the outer zone of BL
    real(kind=8) :: X3AI,X3I,XF3I,ARG3I,X3AM,X3M,XF3M,ARG3M         
    real(kind=8) :: X3AT,X3T,XF3T,ARG3T,X4,c1i,c1m,c1t
    real(kind=8) :: X5I,XF5I,X5M,XF5M,X5T,XF5T, &
    &  X6I,XF6I,X6M,XF6M,X6T,XF6T, MSinner_I1,&
    & MSinner_I2,MSmid_I1,MSmid_I2,MSouter_I1,MSouter_I2,LHS,RHS,MeanU_LHS,MeanU_RHS  ! variables used for inner BL 
    real(kind=8) :: tiinner_I1,tiinner_I2,timid_I1,timid_I2, &
    & tiouter_I1, tiouter_I2,a,b
    real(kind=8) :: I_termK,I_termY ! The contribution of the integral
    real(kind=8),dimension(number_points) :: I1total_inin,I2total_inin,I1total_inmd,I2total_inmd, & 
    &  I1total_inot,I2total_inot,I1total_mdin,I2total_mdin,I1total_mdmd, &
    &  I2total_mdmd,I1total_mdot,I2total_mdot,I1total_otin,I2total_otin, &
    &  I1total_otmd,I2total_otmd,I1total_otot,I2total_otot ! sum of contribution of each layer due to different parts of integral
    real(kind = 8),dimension(number_points) :: I1_total,I2_total,I_total,k1_tilda,pi_tilda,frequency ! sum of integral  
    real(kind = 8),parameter :: xacc = 0.00001,left_bound = 0.d0, &
                &       right_bound = 100.d0 ! variables used for the numerical inversion of function 
    real(kind = 8) :: pdf_k3,block_size,anisotropy,Uc,fact,lambda_22_interp   

      INTEGER(kind=8), PARAMETER :: NP=500,images = 2051 ! number of samples used by data handed over
    INTEGER(kind =8), PARAMETER :: points_intp = 500 ! number of points used for interpolation for Lambda R22, input of gabriele

    REAL(kind = 8), DIMENSION(NP) :: input_variable,target_variable,second_derivative_est_Umean,second_derivative_est_vrms, &
                                 target_vrms,hunt_index,target_lambda_norm
    Real(kind = 8),dimension(NP) :: y_wallnormal,Meanshear_input,vertical_velocity_std,Lambda_norm ! inputs from a text file 
    real(kind = 8), parameter :: yprime_1 = 2.*(10.)**30 ,yprime_n = 0.d0  ! the value of deravative at end points
                                        ! The value 2.*(10.)**30 means algorithm will calculate natural spline see NR for details.
    real(kind =8), dimension(images,NP) :: vertical_velocity
    real(kind =8) :: c1y,c1yprime,c1k ! variables used for variance reduction in y & k 
    real(kind = 8), dimension(points_intp) ::  lambda_22,y_wallnormal_fitted,input_y,second_derivative_est_lambda ! lambda R22 as a function of wall normal distance at location pp26 
   real(kind = 8) , external :: erf_s_sher ! external function to evaluate the error function
   real(kind = 8) :: numer,denom,cons1,cons2,cons3 ! temporary variables used to estimate the value of pdf and contant chosen for variance reduction


!!!!!!!!!!!!!!!!!!!!!Start the clock!!!
        call system_clock(tclock1)
        call cpu_time(t1)

        ! Read the file which contains the mean data
        open(unit=45, file="input_gabriele.txt", status='old')
        do i = 1,NP
                read(45,*) y_wallnormal(i),Meanshear_input(i),vertical_velocity_std(i),Lambda_norm(i)

        end do

        close(45)
	
	! calculated the second derivative using spline function to be used in future 
        
        input_variable = y_wallnormal ! input  y /del
        target_variable = Meanshear_input ! Meanshear input rans
        call spline(input_variable,target_variable,yprime_1,yprime_n,second_derivative_est_Umean)


        target_vrms = vertical_velocity_std ! Vrms input rans
        call spline(input_variable,target_vrms,yprime_1,yprime_n,second_derivative_est_vrms)

	target_lambda_norm = lambda_norm ! Lambda input/del 
	call spline(input_variable,target_lambda_norm,yprime_1,yprime_n,second_derivative_est_lambda)



		! Start the calculation
do i =1,25!number_points ! here number_points = variation of k1
	SUM_total(i) = 0.d0
	SQR_DEVIATION(i) = 0.d0
	SUMIN(i) =0.d0
	SUMMD(i) =0.d0
	SUMOT(i) = 0.d0
	valid_evaluations = 0.
	J = 1 ! intialize index of sum_1000
	step = 100.d0
	! Intialize the Integrals  
	SUMI1_ININ(i)=0.d0
	SUMI2_ININ(i)=0.d0
	block_size = 0.d0 ! intialize block size
	call sobseq(rand,1) ! intialize the random number generator
	CALL init_variables ! call general variables for the calculation

	open(unit=25, file='wall_spectra_calculated.txt', status='unknown') ! save the spectra
	open(unit=15, file= 'PDF.txt', status = 'unknown') ! Save the pdf for comparision and testing
	50  do while(valid_evaluations < itermain) ! Calculate spectra for individual values of K1
                        
                        
		temp = 1
		ind = NP ! intialize the temp variables

			

		call sobseq(rand) ! generate a three dimensional random number
			
		! generate random numbers based on pdf of y and y prime
		! for this we need an estimate of the correlation length/del 
		   ! and for correlation length we need y and y prime so it becomes a chicken egg problem.
		
		cons1 = 0.4d0 ! corresponds to a, see the PDF of the document attached
		cons2 = 5.d0  ! corresponds to b
		
		cons3 = 0.3d0 ! corresponds to h
		! here the random number generated  acording to pdf is X2I for var y
   
		call brent_root(func1,left_bound,right_bound,xacc,rand(1),1.d0,X2I,cons1,cons2,cons3,streamwise_wavenumber(i))


		!print *, 'randout',x2i
		! here the random number generated  acording to pdf is X3I for var yprime
		! here the sixth input is useless for the determination of func1 and hence arbitary value is chosen
		
		
		call brent_root(func1,left_bound,right_bound,xacc,rand(2),1.d0,X3I,cons1,cons2,cons3,streamwise_wavenumber(i))
		
		! here the random number generated acording to pdf is X3I for var yprime
		

		! calculate the pdf of y and y prime
		! pdf of y
		numer = exp(-cons2*streamwise_wavenumber(i)*(x2i-((cons1*cons3)/streamwise_wavenumber(i)))**2.d0)* &
		&           2.d0*sqrt(cons2)*sqrt(streamwise_wavenumber(i))
		denom = sqrt(math_pi)*(erf_s_sher(( cons1*sqrt(cons2)*cons3 )/sqrt(streamwise_wavenumber(i)))+   & 
		 	erf_s_sher(( ((-cons1*cons3)+streamwise_wavenumber(i)) *sqrt(cons2) )/sqrt(streamwise_wavenumber(i)) ) )
		xf2I = numer/denom
		! pdf of yprime
		numer = exp(-cons2*streamwise_wavenumber(i)*(x3i-((cons1*cons3)/streamwise_wavenumber(i)))**2.d0)* &
                &           2.d0*sqrt(cons2)*sqrt(streamwise_wavenumber(i))

		xf3I =  numer/denom

		 
    		c1  = splint(input_variable,target_lambda_norm,second_derivative_est_lambda,X2I)
    		c2 =  splint(input_variable,target_lambda_norm,second_derivative_est_lambda,X3I)

   		c1k = (c1+c2)/2.d0

	!print *, 'c1', c1
	!print *, 'randin' , rand(1)
	!print *, 'c1y', c1y
	!print *, 'x2i',x2i
  ! the value taken for c1k for the varaince reduction taken for variable k3.
		! the impact of this variable is to be seen. 

		! Calculate the random number generated for the variable k3
		
		call brent_root(func2,left_bound,right_bound,xacc,rand(3),c1k,x1,1.d0,1.d0,1.d0,1.d0)
		! Last three inputs are useless for the evaluation of func2 and hence are given an aribatary value 
		! x1 or k3 above is the random variable generated by the pdf function in subroutine
		! inversecdf. 
			
		! The wavenumber k^2 
		BRKSQ = streamwise_wavenumber(i)**2+X1**2
		! The wave number k
		k_bar = sqrt(BRKSQ)
		!print *, 'k1', streamwise_wavenumber(i)
		!Calculate the pdf of k3 
		pdf_k3 = (c1k*( exp(-((x1**2) * c1k**2  )/4.0d0 ) ) ) / ( ( (math_pi)**0.5d0 )* & 
		& erf_s_sher( (c1k*spanwisewavenumber_upperlimit )/2.0d0 ))




		! save the pdf generated to be compared in the matlab script called (pdf_check.m)
	!	write(15,'(e17.10,e17.10,e17.10)') x1,X2I,X3I 

	

! Evaluation of different terms of the integral to be used for the calculation of spectra 

		MSinner_I1 = splint(input_variable,target_variable,second_derivative_est_Umean,x2I) ! Meanshear
		MSinner_I2 = splint(input_variable,target_variable,second_derivative_est_Umean,x3I)




		tiinner_I1 = splint(input_variable,target_vrms,second_derivative_est_Vrms,X2I) ! Vrms
		tiinner_I2 = splint(input_variable,target_vrms,second_derivative_est_Vrms,X3I)
		
		
		c_mean = c1k ! Take the mean correlation length
		
		! Calculate the anisotropy based on eqn 11 of AIAA-51663-AIAAJ-Remmler
		if (streamwise_wavenumber(i) .gt. 5.d0) then
			anisotropy = 1.d0
		else if ( streamwise_wavenumber(i) .le. 1.d0) then
			anisotropy = 3.d0 
		else
			anisotropy = -0.5d0*streamwise_wavenumber(i)+3.5d0
		end if

!print *, 'MSInner_I2', MSinner_I2
!print *, 'x3i',x3I
!print *, 'vrms',tiinner_I2
!print *, 'lambda',c2




		!print *, 'alpha', anisotropy
		!print *, 'k1', streamwise_wavenumber(i)

! The contribution of the term y and the contribution of the term yprime 

		call contributiony(I_termY,X2i,X3i,K_bar,MSinner_I2,tiinner_I2,MSinner_I1,tiinner_I1,xf2I,xf3I,c_mean)



! The contribution of the term K3

		call contributionk(I_termK,anisotropy,BRKSQ,streamwise_wavenumber(i),x1,pdf_k3,c_mean)


! Having generated the pdf, we can verify that the integral is as  flat as possible. To do this we evaluate the integral


! Constants of Integration 




! The complete Integral


I1total_inin = anisotropy*(I_termY*I_termK) ! the complete contribution to the integral for a specific k1


! Sum up the contributions

SUMI1_ININ = I1total_inin+SUMI1_ININ


! Close the loop 

			valid_evaluations = valid_evaluations+1.
                       	block_size = block_size+1.d0
			if (block_size .eq. 1000.d0) then
                               sum_1000(j,I) = (sumi1_inin(i))/valid_evaluations
                                block_size = 0.d0! reset
                                j = j+1.
			endif



	end do ! end of the loop for spectra calculation for a single k1 
	I1_total(i) = (SUMI1_ININ(i))
	fact = 17961314271180.262d0/(0.1356d0**3.d0) ! Factor taken from the code of gabriele, since no information on the dimensions 

!fact = 330261151577.32904d0

	!of input of rans were avialable 
	I_total(I) = fact*(I1_total(I))/valid_evaluations
	K1_tilda(i) = (streamwise_wavenumber(I)/del) ! dimensional streamwise wavenumber 
	Uc = 0.7289094489762658d0 ! value of convection velocity taken by Gabriele from his code 
	frequency(i) = ((K1_tilda(i)*Uc)/(2.d0*math_pi))*(16.0d0/0.1356d0) ! re-dimensionalization 
	PI_tilda(i) = 10.d0*log10( I_total(i)  ) ! Ptilda in db
	

	print * , "freq ", frequency(i)

	print * , "PI ", PI_tilda(i)

	write(25,'(e17.10,e17.10)') frequency(i), PI_tilda(i)

	


end do

close(25)
close(15) ! close the files which were opened
! Now calculat the time for calculation
        call cpu_time(t2)
        call system_clock(tclock2, clock_rate)
        elapsed_time = float(tclock2 - tclock1) / float(clock_rate)
        print *, " "
        print 11, elapsed_time
        11 format("Elapsed time = ",f12.8, " seconds")
        print 12, t2-t1
        12 format("CPU time = ",f12.8, " seconds")



end program main
