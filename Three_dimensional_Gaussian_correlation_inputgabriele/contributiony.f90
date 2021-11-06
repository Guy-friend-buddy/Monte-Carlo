! This subroutine calculates the contribution of y and y prime
subroutine contributiony(I_termY,X2,X3,k_t,ms_yprime,vin_yprime,ms_y,vin_y,pdf_y,pdf_yprime,correlationlength)

USE moduleglobal
      implicit none
      real(kind=8), intent(in) ::  X2,X3,k_t,ms_yprime,vin_yprime,ms_y,vin_y,pdf_y,pdf_yprime,correlationlength
      real(kind=8), intent(out) :: I_termY ! the contribution of Y term
      real(kind=8) :: tempr
     CALL init_variables ! call general variables for the calculation

tempr = ((ms_yprime*vin_yprime*ms_y*vin_y)/(math_pi*2.d0)) 

I_termY =( ( (exp(-K_t*(X2+X3)))*exp(- ((X2-X3)**2 /correlationlength**2)) * tempr  )/((pdf_yprime) * (pdf_y)))

end subroutine contributiony

