    module AxionInterface
    use precision
    use interpolation
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TAxionModel
        integer :: num_perturb_equations = 0
    contains
    procedure :: Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy !Get density perturbation and heat flux for sources
    procedure :: PerturbationInitial
    procedure :: PerturbationEvolve
    procedure :: HiggsVEV

    end type TAxionModel

    public TAxionModel
    contains


    subroutine Init(this, State)
    use classes
    class(TAxionModel), intent(inout) :: this
    class(TCAMBdata), intent(inout), target :: State

    end subroutine Init
	

    subroutine BackgroundDensityAndPressure(this, a, grhoa_t, gpres_ax, phi_ax, phidot_ax)
    class(TAxionModel), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhoa_t
    real(dl), optional, intent(out) :: gpres_ax, phi_ax, phidot_ax
    real(dl) phi, phidot

    grhoa_t=0._dl
    gpres_ax=0._dl

    end subroutine BackgroundDensityAndPressure


    subroutine PerturbedStressEnergy(this, dgrho_a, dgq_a, a, grhoa_t, ay, ayprime, ax_ix)
    class(TAxionModel), intent(inout) :: this
    real(dl), intent(out) :: dgrho_a, dgq_a
    real(dl), intent(in) ::  a, grhoa_t
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: ax_ix

    dgrho_a=0
    dgq_a=0

    end subroutine PerturbedStressEnergy


    subroutine PerturbationEvolve(this, ayprime, ax_ix, a, adotoa, k, z, ay, phi, phidot)
    class(TAxionModel), intent(inout) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa, k, z, ay(:), phi, phidot
    integer, intent(in) :: ax_ix
    end subroutine PerturbationEvolve


    subroutine PerturbationInitial(this, y, a, tau, k)
    class(TAxionModel), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    !Get intinitial values for perturbations at a (or tau)
    !For standard adiabatic perturbations can usually just set to zero to good accuracy

    y = 0._dl

    end subroutine PerturbationInitial


    subroutine HiggsVEV(this, a, ini_vev, vev)
    !Get delta_v
    class(TAxionModel), intent(inout) :: this
    real(dl), intent(in) :: a, ini_vev
    real(dl), intent(out) :: vev
    
    vev = ini_vev
        
    end subroutine HiggsVEV


    end module AxionInterface
