
subroutine subx3(I_termY2prime,X3,MS,VIN,K_BAR,pdf_ybar)

implicit none
      real(kind=8), intent(in) ::  X3,MS,VIN,K_BAR,pdf_ybar
      real(kind=8), intent(out) :: I_termY2prime

   
    ! MS is the mean shear term
    ! VIN is the viscous intensity term
    !X3 is  Y2 prime bar 
   I_termY2prime =((exp(-(K_BAR*X3))*MS*VIN))/(1.*pdf_ybar) 
      
end subroutine subx3
