module moduleglobal

    implicit none
    real (kind=8) :: velocity_ratio,thickness_ratio, &
       &   reynolds_utau,lowerlimit_meanshear,SVK1,SVK2,SVK3,&
       &  SVK4,SVK5, SVK6,SVK7,SVK8,ALIM,&
       &  AK,BK,B1, &
       &  spanwisewavenumber_upperlimit, YSTR 
    integer, parameter :: K=100,itermain=200000,number_points =25
    real(kind =8), parameter :: anisotropy_lowwavenumber = 3.d0, & 
      & anisotropy_midwavenumber = 2.d0, anisotropy_highwavenumber =1.d0 

    real(kind=8), parameter ::VKC=0.41,NU=1.5433e-5,C1IN=4.5154, &
      &  C1MD=1.5012,C1OT=1.1110,UTAU=0.544734710702843d0,math_pi=4.*atan(1.), &
      & UINF=16.,ALFAC=0.837,PI=0.6,DEL=0.042849683d0,DELST=0.002425d0,EM=1.67, &
      & SHEAR_RATIO = 27.d0
	real(kind =8), dimension(number_points), parameter ::&
	& streamwise_wavenumber =  (/0.313d0 ,                 &
	& 3.792511888453973246e-01*1.d0,4.594733048539725306e-01*1.d0, &
	& 5.566646172320730557e-01*1.d0,6.744145803565532171e-01*1.d0, &
	& 8.170719174843563692e-01*1.d0,9.899052271209940468e-01*1.d0, &   
	& 1.199297562078097901e+00*1.d0,1.452982167383452872e+00*1.d0, &
	& 1.760328083279166567e+00*1.d0,2.132686161153358295e+00*1.d0, &
	& 2.583808271411721424e+00*1.d0,3.130355185408624674e+00*1.d0, & 
	& 3.792511888453969249e+00*1.d0,4.594733048539725750e+00*1.d0, &
	& 5.566646172320732333e+00*1.d0,6.744145803565526620e+00*1.d0, &
	& 8.170719174843563692e+00*1.d0,9.899052271209940912e+00*1.d0, &
	& 1.199297562078097990e+01*1.d0,1.452982167383452783e+01*1.d0, &
	& 1.760328083279166833e+01*1.d0,2.132686161153358384e+01*1.d0, &
	& 2.583808271411722401e+01*1.d0,3.130355185408624408e+01*1.d0  /)
      ! k  number of iterations between sucessive display of results
      ! iter  Total number of iteration for a given wavenumber BRK1
      ! BRK1  Wave numeber K for which the calculations are to be performed
      ! ALHPALOW,MID,HIGH are the value of anistropicity factor which should depend on K (BRK1) 
      ! YSTR = value of Y2* at the outer boundary of viscous sublayer
      ! ULIMK=Limit of K3 desired
      ! A and B are the values taken from mean shear calculations taken from Bull 1969 paper.
      ! VKC is the von karman constant taken to be 0.41. NU = kinatic viscosity constant taken (present value based on PS Klebanoff 1954. C1IN,C1MD,C1OT are the values used for variance reduction in inner,middle and  upper layer respectively.
     ! UTAU and UINF are friction and free stream velocity taken from the paper of PS Klebanoff 1954 paper
     ! ALFAC is an emperical constant taken for the estimation of mean shear from Bull 1969 paper
     ! PI is the coles wake parameter. This is equal to 0.6 for a flate plate
     ! DEl and DELST are boundary layer thickness and displacement thickness
     ! EM is another emperical constant taken from Bull 1969 paper.

  contains

    subroutine init_variables
 
      ! Value of wavenumbers for which the calculations are run
      spanwisewavenumber_upperlimit=30.d0 

      ! Calculation of Dependent boundary layer parameter
       YSTR = 8.
       velocity_ratio=UTAU/UINF
       thickness_ratio=DELST/DEL
       reynolds_utau=UTAU*DEL/NU ! Reynolds number based on inner boundary layer parameter
 !      print *, "reynolds number " , reynolds_utau
       ! Calculation of mean shear parameters
       lowerlimit_meanshear=33.2/reynolds_utau ! lower limit of mean shear
       SVK1 = 1./VKC    ! all the constants SVK1, SVK2,.... are the shear velocity constant
       SVK2=ALFAC*VKC
       SVK3=(math_pi*PI)/ALFAC
 !       print *, "svk3 " , svk3
    

       SVK4=math_pi/ALFAC
       SVK5=1./SVK2
       SVK6=1.-ALFAC
       SVK7=EM-1.
       SVK8=reynolds_utau/1.
       
      
       !Sublayer parameter
         ALIM=YSTR/reynolds_utau
         AK=0.75/ALIM**1.5-16.3 ! AK and BK are  constants taken for velocity intensity in inner layer
         BK=-0.45/ALIM**2.5 
        

    end subroutine init_variables    
end module moduleglobal
