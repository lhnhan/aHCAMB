    ! Equations module allowing for fairly general quintessence models
    !
    ! by Antony Lewis (http://cosmologist.info/)

    !!FIX March 2005: corrected update to next treatment of tight coupling
    !!Fix Oct 2011: a^2 factor in ayprime(EV%w_ix+1) [thanks Raphael Flauger]
    ! Oct 2013: update for latest CAMB, thanks Nelson Lima, Martina Schwind
    ! May 2020: updated for CAMB 1.x+

    ! Notes at http://antonylewis.com/notes/CAMB.pdf

    !This module is not well tested, use at your own risk!

    !Need to specify Vofphi function, and also initial_phi
    !You may also need to change other things to get it to work with different types of quintessence model

    !It works backwards, in that it assumes Omega_de is Omega_Q today, then does a binary search on the
    !initial conditions to find what is required to give that Omega_Q today after evolution.

    module Quintessence
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none
    private

    real(dl), parameter :: Tpl= sqrt(kappa*hbar/c**5)  ! sqrt(8 pi G hbar/c^5), reduced planck time

    ! General base class. Specific implemenetations should inherit, defining Vofphi and setting up
    ! initial conditions and interpolation tables
    type, extends(TDarkEnergyModel) :: TQuintessence
        integer :: DebugLevel = 0 !higher then zero for some debug output to console
        real(dl) :: astart = 1e-7_dl
        real(dl) :: integrate_tol = 1e-6_dl
        real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a

        !!Nhan!!
        real(dl) :: grhocc = 0._dl
        real(dl) :: hubble_osc = 1._dl
        real(dl) :: osc_factor = 20._dl
        real(dl) :: grhoa_osc = 0._dl
        real(dl) a_osc

        real(dl) :: grho_to_omega = 2.993431338748136e6_dl
        real(dl) :: GeV_to_Mpl = 4.10964744299298_dl
        real(dl) :: isec_to_eV = 6.582119e-16_dl
        real(dl) :: Mpc_to_m = 3.086e22_dl
        real(dl) :: lightspeed = 2.99792458e8_dl
        !!!!!!!!

        ! Steps for log a and linear spacing, switching at max_a_log (set by Init)
        integer, private :: npoints_linear, npoints_log
        real(dl), private :: dloga, da, log_astart, max_a_log
        real(dl), private, dimension(:), allocatable :: ddphi_a, ddphidot_a
        class(CAMBdata), pointer, private :: State
    contains
    procedure :: Vofphi !V(phi) potential [+ any cosmological constant]
    procedure :: ValsAta_pre
    procedure :: ValsAta !get phi and phi' at scale factor a, e.g. by interpolation in precomputed table
    procedure :: Hubble_func !get Hubble function
    procedure :: Init => TQuintessence_Init
    procedure :: PerturbedStressEnergy => TQuintessence_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TQuintessence_PerturbationEvolve
    procedure :: BackgroundDensityAndPressure => TQuintessence_BackgroundDensityAndPressure
    procedure :: EvolveBackground
    procedure :: EvolveBackgroundLog
    procedure, private :: phidot_start => TQuintessence_phidot_start
    end type TQuintessence

    ! Specific implementation for early quintessence + cosmologial constant, assuming the early component
    ! energy density fraction is negligible at z=0.
    ! The specific parameterization of the potential implemented is the axion model of arXiv:1908.06995
    type, extends(TQuintessence) :: TEarlyQuintessence
        !real(dl) :: n = 3._dl
        !real(dl) :: f = 0.05 ! sqrt(8*pi*G)*f
        !real(dl) :: m = 5d-54 !m in reduced Planck mass units
        !real(dl) :: theta_i = 3.1_dl !initial value of phi/f
		
		!!Nhan!!
        ! real(dl) :: f = 0.411 !sqrt(8*pi*G)*f
        real(dl) :: m = 8.22d-58 !m in reduced Planck mass units
        real(dl) :: omah2 = 0._dl !initial value of phi/f
        ! real(dl) :: theta_i = 0._dl !initial value of phi/f

		! real(dl) :: f_a = 1.e18_dl !just for displaying
		real(dl) :: m_a = 2.e-30_dl !just for displaying
        !!!!!!!!
		
        integer :: npoints = 5000 !baseline number of log a steps; will be increased if needed when there are oscillations
        integer :: min_steps_per_osc = 10

    contains
	!!Nhan!!
	procedure :: Copy => TEarlyQuintessence_Copy
    !!!!!!!!
	procedure :: Vofphi => TEarlyQuintessence_VofPhi
    procedure :: Init => TEarlyQuintessence_Init
    procedure :: ReadParams =>  TEarlyQuintessence_ReadParams
    procedure, nopass :: PythonClass => TEarlyQuintessence_PythonClass
    procedure, nopass :: SelfPointer => TEarlyQuintessence_SelfPointer

    end type TEarlyQuintessence

    procedure(TClassDverk) :: dverk

    public TQuintessence, TEarlyQuintessence
    contains

    function VofPhi(this, phi, deriv)
    !Get the quintessence potential as function of phi
    !The input variable phi is sqrt(8*Pi*G)*psi, where psi is the field
    !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
    !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
    class(TQuintessence) :: this
    real(dl) phi,Vofphi
    integer deriv

    call MpiStop('Quintessence classes must override to provide VofPhi')
    VofPhi = 0

    end function VofPhi


    subroutine TQuintessence_Init(this, State)
    class(TQuintessence), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    !Make interpolation table, etc,
    !At this point massive neutrinos have been initialized
    !so grho_no_de can be used to get density and pressure of other components at scale factor a

    select type(State)
    class is (CAMBdata)
        this%State => State
    end select

    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2

    this%log_astart = log(this%astart)

    end subroutine  TQuintessence_Init

    subroutine TQuintessence_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    !Get grhov_t = 8*pi*rho_de*a**2 and (optionally) equation of state at scale factor a
    class(TQuintessence), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) V, a2, grhov_lambda, phi, phidot

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1_dl
    elseif ((a >= this%astart) .and. (a < this%a_osc)) then
        a2 = a**2
        call this%ValsAta(a,phi,phidot)
        ! V = this%Vofphi(phi,0)
        V = this%Vofphi(phi,0) + this%grhocc
        grhov_t = phidot**2/2 + a2*V
        if (present(w)) then
            w = (phidot**2/2 - a2*V)/grhov_t    !!!Check whether w oscillating around -1
        end if
    elseif (a >= this%a_osc) then
        grhov_t = this%grhoa_osc*this%a_osc**3/a + a**2*this%grhocc !DM approximation
        if (present(w)) then
            w = -1._dl
        end if
    else
        grhov_t=0
        if (present(w)) w = -1_dl
    end if

    end subroutine TQuintessence_BackgroundDensityAndPressure

    real(dl) function Hubble_func(this, a, phi, phidot)
    class(TQuintessence) :: this
    real(dl) a
    real(dl), optional :: phi, phidot
    real(dl) grhonoax, grhotot
    real(dl) :: grhoax = 0._dl

    ! grhoax = 0.5d0*phidot**2/a**2 + this%Vofphi(phi,0) !uncomment to include backreaction
    grhonoax = this%state%grho_no_de(a)/a**4 + this%grhocc
    grhotot = grhonoax + grhoax
    Hubble_func = sqrt(grhotot/3.0d0)
    
    end function Hubble_func

    subroutine EvolveBackgroundLog(this,num,loga,y,yprime)
    ! Evolve the background equation in terms of loga.
    ! Variables are phi=y(1), a^2 phi' = y(2)
    ! Assume otherwise standard background components
    class(TQuintessence) :: this
    integer num
    real(dl) y(num),yprime(num)
    real(dl) loga, a

    a = exp(loga)
    call this%EvolveBackground(num, a, y, yprime)
    yprime = yprime*a !dy/dloga = dy/da * a

    end subroutine EvolveBackgroundLog

    subroutine EvolveBackground(this,num,a,y,yprime)
    ! Evolve the background equation in terms of a.
    ! Variables are phi=y(1), a^2 phi' = y(2)
    ! Assume otherwise standard background components
    class(TQuintessence) :: this
    integer num
    real(dl) y(num),yprime(num)
    real(dl) a, a2, hubble
    real(dl) phi, phidot, adot
    
    a2=a**2
    phi = y(1)
    phidot = y(2)/a2

    hubble = this%Hubble_func(a,phi,phidot)
    this%hubble_osc = hubble    !!!!!!!!!!!!!

    adot=hubble*a2 !a'
    yprime(1)=phidot/adot !dy(1)/da = dphi/da
    yprime(2)= -a2**2*this%Vofphi(phi,1)/adot !dy(2)/da

    end subroutine EvolveBackground


    real(dl) function TQuintessence_phidot_start(this,phi)
    class(TQuintessence) :: this
    real(dl) :: phi

    TQuintessence_phidot_start = 0

    end function TQuintessence_phidot_start
	
	!!Nhan!!
    subroutine ValsAta(this,a,aphi,aphidot)
    class(TQuintessence) :: this
    !Do interpolation for background phi and phidot at a (precomputed in Init)
    real(dl) a, aphi, aphidot
    real(dl) a0,b0,ho2o6,delta,da
    integer ix

    if (a >= 0.9999999d0) then !!scale factor bigger than nowadays(not necessary)!!
        aphi= this%State%phi_a(this%npoints_linear+this%npoints_log)
        aphidot= this%State%phidot_a(this%npoints_linear+this%npoints_log)
        return
    elseif (a < this%astart) then
        aphi = this%State%phi_a(1)
        aphidot = 0
        return
    elseif (a > this%max_a_log) then    !!Linear interpolation(not necessary)!!
        delta= a-this%max_a_log
        ix = this%npoints_log + int(delta/this%da)
    else
        delta= log(a)-this%log_astart
        ix = int(delta/this%dloga)+1
    end if

    if (ix + 1 > this%npoints_log) then
        call MpiStop("Not enough data for interpolation")
    end if

    da = this%State%sampled_a(ix+1) - this%State%sampled_a(ix)
    a0 = (this%State%sampled_a(ix+1) - a)/da
    b0 = 1 - a0
    ho2o6 = da**2/6._dl
    aphi=b0*this%State%phi_a(ix+1) + a0*(this%State%phi_a(ix)-b0*((a0+1)*this%State%ddphi_a(ix)+(2-a0)*this%State%ddphi_a(ix+1))*ho2o6)
    aphidot=b0*this%State%phidot_a(ix+1) + a0*(this%State%phidot_a(ix)-b0*((a0+1)*this%State%ddphidot_a(ix)+(2-a0)*this%State%ddphidot_a(ix+1))*ho2o6)

    end subroutine ValsAta


    subroutine ValsAta_pre(this,a,aphi,aphidot)
    class(TQuintessence) :: this
    !Do interpolation for background phi and phidot at a (precomputed in Init)
    real(dl) a, aphi, aphidot
    real(dl) a0,b0,ho2o6,delta,da
    integer ix
    
    if (a < this%astart) then
        aphi = this%State%phi_a(1)
        aphidot = 0
        return
    else
        delta= log(a)-this%log_astart
        ix = int(delta/this%dloga)+1
    end if

    if (ix + 1 > this%npoints_log) then
        call MpiStop("Not enough data for interpolation")
    end if

    da = this%sampled_a(ix+1) - this%sampled_a(ix)
    a0 = (this%sampled_a(ix+1) - a)/da
    b0 = 1 - a0
    ho2o6 = da**2/6._dl
    aphi=b0*this%phi_a(ix+1) + a0*(this%phi_a(ix)-b0*((a0+1)*this%ddphi_a(ix)+(2-a0)*this%ddphi_a(ix+1))*ho2o6)
    aphidot=b0*this%phidot_a(ix+1) + a0*(this%phidot_a(ix)-b0*((a0+1)*this%ddphidot_a(ix)+(2-a0)*this%ddphidot_a(ix+1))*ho2o6)

    end subroutine ValsAta_pre
	!!!!!!!!
	
    subroutine TQuintessence_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    !Get density perturbation and heat flux
    class(TQuintessence), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) ::  a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) phi, phidot, clxq, vq

    call this%ValsAta(a,phi,phidot)
    clxq=ay(w_ix)
    vq=ay(w_ix+1)
    dgrhoe= phidot*vq +clxq*a**2*this%Vofphi(phi,1)     !!! density source
    dgqe= k*phidot*clxq     !!! heat-flux source

    end subroutine TQuintessence_PerturbedStressEnergy


    subroutine TQuintessence_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    !Get conformal time derivatives of the density perturbation and velocity
    class(TQuintessence), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) clxq, vq, phi, phidot


    call this%ValsAta(a,phi,phidot) !wasting time calling this again..
    clxq=y(w_ix)
    vq=y(w_ix+1)
    ayprime(w_ix)= vq
    ayprime(w_ix+1) = - 2*adotoa*vq - k*z*phidot - k**2*clxq - a**2*this%Vofphi(phi,2)*clxq

    end subroutine TQuintessence_PerturbationEvolve

    ! Early Quintessence example, axion potential from e.g. arXiv: 1908.06995

    function TEarlyQuintessence_VofPhi(this, phi, deriv) result(VofPhi)
    !The input variable phi is sqrt(8*Pi*G)*psi
    !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
    !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
    class(TEarlyQuintessence) :: this
    real(dl) phi,Vofphi
    integer deriv
    real(dl) theta
    real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2

    ! Assume f = sqrt(kappa)*f_theory = f_theory/M_pl
    ! m = m_theory/M_Pl
    ! theta = phi/this%f
    if (deriv==0) then
        ! Vofphi = units*this%m**2*this%f**2*(1 - cos(theta)) + this%State%grhov
        Vofphi = units*this%m**2*phi**2/2
    else if (deriv ==1) then
        ! Vofphi = units*this%m**2*this%f*sin(theta)
        Vofphi = units*this%m**2*phi
    else if (deriv ==2) then
        ! Vofphi = units*this%m**2*cos(theta)
        Vofphi = units*this%m**2
    end if

    end function TEarlyQuintessence_VofPhi


    subroutine TEarlyQuintessence_Init(this, State)
    class(TEarlyQuintessence), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) aend, afrom
    integer, parameter ::  NumEqs=2
    real(dl) c(24),w(NumEqs,9), y(NumEqs)
    integer ind, i, ix
    real(dl), parameter :: splZero = 0._dl
    real(dl) lastsign, da_osc, last_a
    real(dl) initial_phi, initial_phidot, a2
    real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a
    real(dl) scaling
    integer npoints, tot_points
    
    real(dl) ma_hubble
    real(dl) phi_test, phidot_test, hubble_test, hubble_1
    real(dl) loga_1, loga_2, loga_test

    !Make interpolation table, etc,
    !At this point massive neutrinos have been initialized
    !so grho_no_de can be used to get density and pressure of other components at scale factor a

    call this%TQuintessence%Init(State)

    this%dloga = (-this%log_astart)/(this%npoints-1)

    !use log spacing in a up to max_a_log, then linear. Switch where step matches
    this%max_a_log = 1.d0/this%npoints/(exp(this%dloga)-1)
    npoints = (log(this%max_a_log)-this%log_astart)/this%dloga + 1

    if (allocated(this%phi_a)) then
        deallocate(this%phi_a,this%phidot_a)
        deallocate(this%ddphi_a,this%ddphidot_a, this%sampled_a)
    end if
    allocate(phi_a(npoints),phidot_a(npoints), sampled_a(npoints))
	
    this%m = this%m_a * this%GeV_to_Mpl * (1.e-28_dl)
    ma_hubble = this%m_a*this%Mpc_to_m/this%lightspeed/this%isec_to_eV / this%osc_factor
    this%grhocc = this%State%grhov - (this%omah2/this%grho_to_omega)

    initial_phi = 1._dl                             !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initial_phi = this%theta_i*this%f


    y(1)=initial_phi
    initial_phidot =  this%astart*this%phidot_start(initial_phi)
    y(2)= initial_phidot*this%astart**2

    phi_a(1)=y(1)
    phidot_a(1)=y(2)/this%astart**2
    sampled_a(1)=this%astart
    da_osc = 1
    last_a = this%astart

    ind=1
    afrom=this%log_astart
    do i=1, npoints-1
        aend = this%log_astart + this%dloga*i
        ix = i+1
        sampled_a(ix)=exp(aend)
        a2 = sampled_a(ix)**2
        call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
        call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
        phi_a(ix)=y(1)
        phidot_a(ix)=y(2)/a2


        if (i==1) then
            lastsign = y(2)
        elseif (y(2)*lastsign < 0) then
            !derivative has changed sign. Use to probe any oscillation scale:
            da_osc = min(da_osc, exp(aend) - last_a)
            last_a = exp(aend)
            lastsign= y(2)
        end if

        if (sampled_a(ix)*(exp(this%dloga)-1)*this%min_steps_per_osc > da_osc) then
            !Step size getting too big to sample oscillations well
            print*, "Break because of too big step size"
            exit
        end if

        if (this%hubble_osc < ma_hubble) then 
            !Break loop when oscillation starts
            print*, "Break because of oscillation condition"
            exit 
        end if

        if (i == npoints-1) then
            print*, "Finish loop without breaking"
        end if
    
    end do

    this%npoints_log = ix

    allocate(this%phi_a(ix),this%phidot_a(ix))
    allocate(this%sampled_a(ix), this%ddphi_a(ix),this%ddphidot_a(ix))

    this%sampled_a(1:ix) = sampled_a(1:ix)
    this%phi_a(1:ix) = phi_a(1:ix)
    this%phidot_a(1:ix) = phidot_a(1:ix)

    call spline(this%sampled_a,this%phi_a,ix,splZero,splZero,this%ddphi_a)          !!Calculate second derivatives of phi!!
    call spline(this%sampled_a,this%phidot_a,ix,splZero,splZero,this%ddphidot_a)

    hubble_1 = this%hubble_osc
    this%a_osc = sampled_a(ix)
    loga_1 = log(this%a_osc)
    loga_2 = log(sampled_a(ix-1))

    ! print*, "Current log scale factor 1 ", loga_1
    ! print*, "Current log scale factor 2 ", loga_2
    ! print*, "Current Hubble 1", hubble_1/this%Mpc_to_m*this%lightspeed*this%isec_to_eV

    do while (abs(loga_1-loga_2)>1e-5)
        loga_test = (loga_1+loga_2)/(2.0d0)
        call this%ValsAta_pre(exp(loga_test),phi_test,phidot_test)
        hubble_test = this%Hubble_func(exp(loga_test),phi_test,phidot_test)
                
        if ((hubble_test-ma_hubble)*(hubble_1-ma_hubble) < 0) then
            loga_2 = loga_test
        else
            loga_1 = loga_test
            hubble_1 = hubble_test
        end if
    end do
    this%a_osc = exp(loga_test)
    this%grhoa_osc = (phidot_test/this%a_osc)**2/(2.0d0) + this%Vofphi(phi_test,0)

    !!Rescale axion!!
    scaling = (this%grhoa_osc*this%a_osc**3)*this%grho_to_omega/this%omah2
    do i=1, ix
        this%phi_a(i) = this%phi_a(i) / sqrt(scaling)
        this%phidot_a(i) = this%phidot_a(i) / sqrt(scaling)
    end do
    this%grhoa_osc = this%grhoa_osc / scaling
    !!!!!!!!!!!!!!!!! 

    end subroutine TEarlyQuintessence_Init

	
	!!Nhan!!
	subroutine TEarlyQuintessence_Copy(this, State)
    class(TEarlyQuintessence), intent(inout) :: this
    class(TCAMBdata), intent(inout), target :: State
    real(dl) a, ghroa_t, delta_a, w
    integer i, j

    call this%TQuintessence%Init(State)

    this%State%phi_a(1:this%npoints_log) = this%phi_a(1:this%npoints_log)
    this%State%phidot_a(1:this%npoints_log) = this%phidot_a(1:this%npoints_log)
    this%State%sampled_a(1:this%npoints_log) = this%sampled_a(1:this%npoints_log)
    this%State%ddphi_a(1:this%npoints_log) = this%ddphi_a(1:this%npoints_log)
    this%State%ddphidot_a(1:this%npoints_log) = this%ddphidot_a(1:this%npoints_log)

    deallocate(this%phi_a, this%phidot_a)
    deallocate(this%sampled_a, this%ddphi_a, this%ddphidot_a)
	
	! this%State%omah2 = this%omah2
	this%State%m_a = this%m_a

    print*, "Break at ix = ", this%npoints_log

    delta_a = (1._dl - this%State%sampled_a(this%npoints_log))/(8832-this%npoints_log)
    do i=this%npoints_log+1, 8832
        j = i - this%npoints_log
        this%State%sampled_a(i) = this%State%sampled_a(this%npoints_log) + j*delta_a
    end do

    open(1, file='omah2_approx.dat', status='REPLACE')
    do i=1, 8832
        a = this%State%sampled_a(i)
        call this%BackgroundDensityAndPressure(this%State%grhov, a, ghroa_t, w)
        ghroa_t = ghroa_t - a**2*this%grhocc
        write(1,*) a, ghroa_t/a**2*this%grho_to_omega, w
    end do
    close(1)

    call MpiStop('Finish writing')

    end subroutine  TEarlyQuintessence_Copy
	!!!!!!!!
    

    subroutine TEarlyQuintessence_ReadParams(this, Ini)
    use IniObjects
    class(TEarlyQuintessence) :: this
    class(TIniFile), intent(in) :: Ini
    
    call this%TDarkEnergyModel%ReadParams(Ini)
	
	this%m_a = Ini%Read_Double('EarlyQuintessence_m')
    
	! this%theta_i = Ini%Read_Double('EarlyQuintessence_theta_i')
    this%omah2 = Ini%Read_Double('omah2')
	
    end subroutine TEarlyQuintessence_ReadParams


    function TEarlyQuintessence_PythonClass()
    character(LEN=:), allocatable :: TEarlyQuintessence_PythonClass

    TEarlyQuintessence_PythonClass = 'EarlyQuintessence'

    end function TEarlyQuintessence_PythonClass

    subroutine TEarlyQuintessence_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TEarlyQuintessence), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TEarlyQuintessence_SelfPointer


    end module Quintessence
