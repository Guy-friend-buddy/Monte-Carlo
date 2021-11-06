
subroutine subx2(I_termY2,X2,MS,VIN,K_BAR,pdf_y)

implicit none
      real(kind=8), intent(in) ::  X2,MS,VIN,K_BAR,pdf_y
      real(kind=8), intent(out) :: I_termY2



   I_termY2 =((exp(-(K_BAR*X2))*MS*VIN))/(pdf_y) 
      
end subroutine subx2
