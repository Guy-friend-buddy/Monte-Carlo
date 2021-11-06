! This subroutine calculates the contribution of k 

subroutine contributionk(I_termK,alpha,K_sqr,k1,k3,pdf_k,correlationlength)

use moduleglobal

implicit none
      real(kind=8), intent(in) ::  alpha,K_sqr,k1,k3,pdf_k,correlationlength
      real(kind=8), intent(out) :: I_termK ! the contribution of K terms
      real(kind =8) :: tempr,tempr2,tempr3 ! temp variable used to make the multiplication look neat

	CALL init_variables ! call general variables for the calculation

tempr = alpha**2*(k1**2)+k3**2
tempr2 = exp(-(tempr)*(correlationlength**2/4.d0))
tempr3 = (k1**2.d0/k_sqr)*(tempr)*(correlationlength**4*(del**3))
I_termK=(tempr3*(tempr2)) /(pdf_k)

end subroutine contributionk






