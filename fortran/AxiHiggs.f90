    module AxiHiggs
    use AxionStandard
    use constants
    use results
    use classes
    implicit none

    private
    
    real(dl), parameter :: Mpl_GeV = 2.435415204336538e18_dl
    real(dl), parameter :: phi_to_GeV = 2.437532984704441e18_dl

    type, extends(TAxionStandard) :: TAxiHiggs
        
    contains
    
    procedure :: HiggsVEV => TAxiHiggs_HiggsVEV

    end type TAxiHiggs

    public TAxiHiggs
    contains


    subroutine TAxiHiggs_HiggsVEV(this, a, ini_vev, vev)
    !Get delta_v
    class(TAxiHiggs), intent(inout) :: this
    real(dl), intent(in) :: a, ini_vev
    real(dl), intent(out) :: vev
    real(dl) phi, phidot, grhoa_t, gpres, phi2
    
    if (this%omah2 > this%omah2_min) then
        if (a < this%a_osc) then
            call this%BackgroundDensityAndPressure(a, grhoa_t, gpres, phi, phidot)
            phi2 = phi**2
        else !DM approximation
            phi2 = this%phi2_osc*(this%a_osc/a)**3
        end if
        vev = 1._dl + (ini_vev-1._dl)*phi2/(this%initial_phi**2)
    else
        vev = 1._dl   !This can be either ini_vev or 1 depending on limit definition
    end if
    
    end subroutine TAxiHiggs_HiggsVEV


    end module AxiHiggs
