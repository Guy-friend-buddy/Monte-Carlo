module inversecdf
USE moduleglobal
contains

real(kind =8) function func2(x,rand_in,c1,a1,b1,h1,k3in1) ! pdf function for the variable k3
    implicit none
    real(kind =8), intent(in) :: x,c1,rand_in,a1,b1,k3in1,h1
    real(kind = 8), external  :: erf_s_sher ! the error function used for variance reduction. 
 	      
    func2 = (rand_in*erf_s_sher( ( c1*spanwisewavenumber_upperlimit)/2.d0 )  ) - erf_s_sher( (c1*x)/2.d0 )
end function func2

	real(kind=8) function func1(x,rand_in,c2,a,b,h,k3in) ! pdf function for the variable y and yprime
    implicit none
    real(kind =8), intent(in) :: x,a,b,c2,h,rand_in,k3in ! here c2 is the integral length scale
    real(kind =8), external  :: erf_s_sher ! the error function used for variance reduction. 
   real (kind=8) ::  temp3,temp4,temp5,temp6 ! temprory variable
	temp3 = (sqrt(b)*a*h)/(sqrt(k3in))
	temp4 = (sqrt(b)*(k3in-(a*h)))/(sqrt(k3in))
	temp5 = (sqrt(b)*((k3in*x)-(a*h)))/(sqrt(k3in))
	temp6 = (sqrt(b)*a*h)/(sqrt(k3in))
!	print *, 'k1', k3in
!        print *, 'a', a
!        print *, 'b', b
!        print *, 'h', h
!	print *, 'temp3', temp3
!	print *, 'temp4', temp4
!	print *, 'temp5', temp5
!	print *, 'temp6', temp6
func1 = (rand_in*(erf_s_sher(temp3)+erf_s_sher(temp4))) - (erf_s_sher(temp5)+erf_s_sher((temp6)))

end function func1

end module inversecdf
