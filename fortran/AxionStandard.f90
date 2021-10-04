    module AxionStandard
    use AxionInterface
    use constants
    use results
    use classes
    implicit none

    private

    real(dl), parameter :: Tpl= sqrt(kappa*hbar/c**5)  ! sqrt(8 pi G hbar/c^5), reduced planck time
    real(dl), parameter :: units = MPC_in_sec**2 /Tpl**2  !convert to units of 1/Mpc^2
    real(dl), parameter :: grho_to_omega = 2.993431338748136e6_dl
    real(dl), parameter :: GeV_to_Mpl = 4.10964744299298_dl
    real(dl), parameter :: isec_to_eV = 6.582119e-16_dl
    real(dl), parameter :: imeter_to_eV = 1.97327e-7_dl
    real(dl), parameter :: kg_to_eV = 1.782662e-36_dl
    real(dl), parameter :: kappa_eV = kappa*kg_to_eV/(imeter_to_eV**3)*(isec_to_eV**2)
    real(dl), parameter :: Mpl_GeV = (1.e-9_dl)/sqrt(kappa_eV)
    real(dl), parameter :: phi_to_GeV = (GeV_to_Mpl*(1.e-37_dl))/Tpl*isec_to_eV/sqrt(kappa_eV)

    type, extends(TAxionModel) :: TAxionStandard

        ! real(dl) :: theta_i = 0._dl !initial value of phi/f
        ! real(dl) :: f_a = 1.e18_dl !just for displaying
        real(dl) :: m = 8.e-58_dl !m in reduced Planck mass units
        real(dl) :: initial_phi = 0._dl
        real(dl) :: omah2 = 0._dl
        real(dl) :: omah2_min = 1.e-5_dl !minimum value of axion density
        real(dl) :: m_a = 1.e-29_dl !just for displaying

        real(dl) :: hubble_osc = 1._dl
        real(dl) :: osc_factor = 3._dl
        real(dl) :: hubble_ma = 1._dl
        real(dl) :: grhoa_osc = 0._dl
        real(dl) :: phi2_osc = 0._dl
        real(dl) a_osc
        logical  :: first_iter = .true.

        integer :: npoints = 10000 !baseline number of log a steps;
        integer :: npoints_coarse = 1000 !number of coarse log a steps;
        integer :: max_npoints = 30000
        integer :: npoints_log = 0 !final number of log a steps
        real(dl) :: astart = 1e-9_dl
        real(dl) :: a_coarse = 1e-5_dl
        real(dl) :: integrate_tol = 1e-8_dl

        ! Steps for log a spacing, switching at max_a_log (set by Init)
        real(dl), private :: dloga, dloga_coarse, da, log_astart, log_acoarse
        class(CAMBdata), pointer, private :: State
        
    contains
    
    procedure :: phi_interp !get phi and phi' at scale factor a, e.g. by interpolation in precomputed table
    ! procedure :: ValsAta !get phi and phi' at scale factor a, e.g. by interpolation in precomputed table
    procedure :: Hubble_func !get Hubble function
    procedure :: EvolveBackground
    procedure :: EvolveBackgroundLog

    procedure :: ReadParams => TAxionStandard_ReadParams
    procedure :: VofPhi => TAxionStandard_VofPhi
    procedure :: Init => TAxionStandard_Init
    procedure :: PerturbedStressEnergy => TAxionStandard_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TAxionStandard_PerturbationEvolve
    procedure :: PerturbationInitial => TAxionStandard_PerturbationInitial
    procedure :: BackgroundDensityAndPressure => TAxionStandard_BackgroundDensityAndPressure
    procedure :: EquationOfState => TAxionStandard_EquationOfState
    procedure :: AdiabaticSoundSpeed => TAxionStandard_AdiabaticSoundSpeed
    procedure, private :: AxionDynamicSolver => TAxionStandard_AxionDynamicSolver
    procedure, private :: phidot_start => TAxionStandard_phidot_start
    
    end type TAxionStandard

    procedure(TClassDverk) :: dverk

    public TAxionStandard
    contains

    real(dl) function Hubble_func(this, a, phi, phidot)
    class(TAxionStandard) :: this
    real(dl) a
    real(dl), optional :: phi, phidot
    real(dl) grhonoax, grhotot
    real(dl) :: grhoax = 0._dl

    grhoax = phidot**2/2/a**2 + this%Vofphi(phi,0) !uncomment to include backreaction
    grhonoax = this%state%grho_no_de(a)/a**4 + this%State%grhov
    grhotot = grhonoax + grhoax
    Hubble_func = sqrt(grhotot/3.0d0)
    
    end function Hubble_func


    subroutine EvolveBackgroundLog(this,num,loga,y,yprime)
    ! Evolve the background equation in terms of loga.
    ! Variables are phi=y(1), a^2 phi' = y(2)
    ! Assume otherwise standard background components
    class(TAxionStandard) :: this
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
    class(TAxionStandard) :: this
    integer num
    real(dl) y(num),yprime(num)
    real(dl) a, a2, hubble
    real(dl) phi, phidot, adot
    
    a2=a**2
    phi = y(1)
    phidot = y(2)/a2

    hubble = this%Hubble_func(a,phi,phidot)
    this%hubble_osc = hubble

    adot=hubble*a2 !a'
    yprime(1)=phidot/adot !dy(1)/da = dphi/da
    yprime(2)= -a2**2*this%Vofphi(phi,1)/adot !dy(2)/da

    end subroutine EvolveBackground


    real(dl) function TAxionStandard_phidot_start(this,phi)
    class(TAxionStandard) :: this
    real(dl) :: phi

    TAxionStandard_phidot_start = 0

    end function TAxionStandard_phidot_start


    function TAxionStandard_VofPhi(this, phi, deriv) result(VofPhi)
    !The input variable phi is sqrt(8*Pi*G)*psi
    !Returns (8*Pi*G)^(1-deriv/2)*d^{deriv}V(psi)/d^{deriv}psi evaluated at psi
    !return result is in 1/Mpc^2 units [so times (Mpc/c)^2 to get units in 1/Mpc^2]
    class(TAxionStandard) :: this
    real(dl) phi,Vofphi
    integer deriv
    real(dl) theta
    
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
    
    end function TAxionStandard_VofPhi


    subroutine TAxionStandard_AdiabaticSoundSpeed(this, a, adotoa, phi, phidot, csquared_ad)
    class(TAxionStandard), intent(inout) :: this
    real(dl), intent(in) :: a, adotoa, phi, phidot
    real(dl), intent(out) :: csquared_ad
    real(dl) a2, con_hubble, csquared_ad_ini
        
    if (a < this%a_osc) then
        csquared_ad_ini = -(7.d0/3.d0)
        a2 = a**2
        con_hubble = adotoa
        csquared_ad = 1.d0 + (2.d0/3.d0)*a2*this%Vofphi(phi,1)/(con_hubble*phidot)
        if ((csquared_ad < csquared_ad_ini) .or. (abs(phidot) < 1.e-50_dl)) then
            csquared_ad = csquared_ad_ini
        end if
    else
        csquared_ad = 0._dl !DM approximation
    end if

    end subroutine TAxionStandard_AdiabaticSoundSpeed


    subroutine TAxionStandard_EquationOfState(this, a, phi, phidot, wa_eff)
    class(TAxionStandard), intent(inout) :: this
    real(dl), intent(in) :: a, phi, phidot
    real(dl), intent(out) :: wa_eff
    real(dl) a2, V
    
    if (a < this%a_osc) then
        a2 = a**2
        V = this%Vofphi(phi,0)
        wa_eff = (phidot**2/2/a2-V)/(phidot**2/2/a2+V)
    else
        wa_eff = 0._dl !DM approximation
    end if
        
    end subroutine TAxionStandard_EquationOfState


    subroutine TAxionStandard_BackgroundDensityAndPressure(this, a, grhoa_t, gpres_ax, phi_ax, phidot_ax)
    !Get grhoa_t = 8*pi*rho_a*a**2
    class(TAxionStandard), intent(inout) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: grhoa_t
    real(dl), optional, intent(out) :: gpres_ax, phi_ax, phidot_ax
    real(dl) phi, phidot
    
    phi = 0._dl
    phidot = 0._dl

    if (this%omah2 > this%omah2_min) then
        if (a < this%a_osc) then
            ! call this%ValsAta(a,phi,phidot)
            call this%phi_interp(a,phi,phidot)
            grhoa_t = phidot**2/2 + a**2*this%Vofphi(phi,0)
            if (present(gpres_ax)) then
                gpres_ax = phidot**2/2 - a**2*this%Vofphi(phi,0)
            end if
        else !DM approximation
            grhoa_t = this%grhoa_osc*this%a_osc**3/a
            if (present(gpres_ax)) gpres_ax = 0._dl
        end if
        if (present(phi_ax)) phi_ax = phi
        if (present(phidot_ax)) phidot_ax = phidot
    else
        grhoa_t = 0._dl
        if (present(gpres_ax)) gpres_ax = 0._dl
        if (present(phi_ax)) phi_ax = 0._dl
        if (present(phidot_ax)) phidot_ax = 0._dl
    end if
    
    end subroutine TAxionStandard_BackgroundDensityAndPressure


    subroutine TAxionStandard_PerturbedStressEnergy(this, dgrho_a, dgq_a, a, grhoa_t, ay, ayprime, ax_ix)
    !Get density perturbation and heat flux multiplying a^2
    class(TAxionStandard), intent(inout) :: this
    real(dl), intent(out) :: dgrho_a, dgq_a
    real(dl), intent(in) ::  a, grhoa_t
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: ax_ix
    real(dl) phi, phidot, clxa, ua

    clxa=ay(ax_ix)
    ua=ay(ax_ix+1)
    
    dgrho_a= clxa*grhoa_t
    dgq_a= ua*grhoa_t

    end subroutine TAxionStandard_PerturbedStressEnergy


    subroutine TAxionStandard_PerturbationEvolve(this, ayprime, ax_ix, a, adotoa, k, z, ay, phi, phidot)
    !Get conformal time derivatives of the density perturbation and velocity
    class(TAxionStandard), intent(inout) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, ay(:), phi, phidot
    integer, intent(in) :: ax_ix
    real(dl) clxa, ua, wa_eff, csquared_ad, csquared_a, ka_const

    clxa=ay(ax_ix)
    ua=ay(ax_ix+1)
    
    if (this%omah2 > this%omah2_min) then
        if (a < this%a_osc) then
            call this%EquationOfState(a, phi, phidot, wa_eff)
            call this%AdiabaticSoundSpeed(a, adotoa, phi, phidot, csquared_ad)
            ayprime(ax_ix)= -k*ua - (1.d0+wa_eff)*k*z - 3*adotoa*(1.d0-wa_eff)*clxa - 9*adotoa**2*(1.d0-csquared_ad)*ua/k
            ayprime(ax_ix+1)= 2*adotoa*ua + k*clxa + 3*adotoa*(wa_eff-csquared_ad)*ua 
        else !DM approximation
            ka_const = k**2/(4*this%Vofphi(phi,2)*a**2)  !!Potential dependence!!
            csquared_a = ka_const/(1.d0+ka_const)
            ayprime(ax_ix)= -k*ua - k*z - 3*adotoa*csquared_a*clxa - 9*adotoa**2*csquared_a*ua/k
            ayprime(ax_ix+1)= -adotoa*ua + k*csquared_a*clxa + 3*csquared_a*adotoa*ua
        end if
    else
        ayprime(ax_ix)=0.d0
        ayprime(ax_ix+1)=0.d0
    end if

    end subroutine TAxionStandard_PerturbationEvolve


    subroutine TAxionStandard_PerturbationInitial(this, y, a, tau, k)
    class(TAxionStandard), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k
    
    y = 0._dl !Adiabatic initial condition

    !Isocurvature initial condition can be added in here
    
    end subroutine TAxionStandard_PerturbationInitial


    subroutine TAxionStandard_AxionDynamicSolver(this, State, initial_phi, omah2, npoints_log)

    class(TAxionStandard), intent(inout) :: this
    class(TCAMBdata), intent(inout), target :: State
    real(dl) aend, afrom
    integer, parameter ::  NumEqs=2
    real(dl) c(24),w(NumEqs,9), y(NumEqs)
    integer ind, i, ix
    real(dl), parameter :: splZero = 0._dl
    real(dl) initial_phidot, a2
    real(dl), dimension(:), allocatable :: sampled_a, phi_a, phidot_a
    
    real(dl), intent(in) :: initial_phi
    real(dl), intent(out) :: omah2
    integer, intent(out), optional :: npoints_log
    real(dl) phi_test, phidot_test, hubble_test
    real(dl) loga_1, loga_test


    y(1)=initial_phi
    initial_phidot =  this%astart*this%phidot_start(initial_phi)
    y(2)= initial_phidot*this%astart**2
    
    allocate(phi_a(this%max_npoints),phidot_a(this%max_npoints), sampled_a(this%max_npoints))

    phi_a(1)=y(1)
    phidot_a(1)=y(2)/this%astart**2
    sampled_a(1)=this%astart

    ind=1
    afrom=this%log_astart
    do i=1, this%npoints_coarse-1

        aend = this%log_astart + this%dloga_coarse*i
        ix = i+1
        sampled_a(ix)=exp(aend)
        a2 = sampled_a(ix)**2
        call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
        call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
        phi_a(ix)=y(1)
        phidot_a(ix)=y(2)/a2
        
        if (this%first_iter) then

            if (this%hubble_osc < this%hubble_ma*this%osc_factor) then
                this%a_coarse = sampled_a(ix)
                this%log_acoarse = log(this%a_coarse)
                this%first_iter = .false.
                ! print*, "acoarse is: ", this%a_coarse
            end if

        end if

        if (this%hubble_osc < this%hubble_ma) exit

    end do
    
    
    if (ix >= this%npoints_coarse) then

        ind=1
        afrom=this%log_acoarse
        do i=1, this%npoints
            aend = this%log_acoarse + this%dloga*i
            ix = ix+1
            if (ix > this%max_npoints) call MpiStop("High oscillation factor requires max_npoints!")
            sampled_a(ix)=exp(aend)
            a2 = sampled_a(ix)**2
            call dverk(this,NumEqs,EvolveBackgroundLog,afrom,y,aend,this%integrate_tol,ind,c,NumEqs,w)
            call EvolveBackgroundLog(this,NumEqs,aend,y,w(:,1))
            phi_a(ix)=y(1)
            phidot_a(ix)=y(2)/a2
            
            if (this%hubble_osc < this%hubble_ma) exit  !Break loop when effective regime starts

        end do

    end if

    this%npoints_log = ix
    if (present(npoints_log)) npoints_log = ix
                
    this%State%sampled_a(1:ix) = sampled_a(1:ix)
    this%State%phi_a(1:ix) = phi_a(1:ix)
    this%State%phidot_a(1:ix) = phidot_a(1:ix)
    
    ! Allocate enough memory for this%State arrays
        
    this%a_osc = sampled_a(ix)
    loga_test = log(sampled_a(ix))
    phi_test = this%State%phi_a(ix)
    phidot_test = this%State%phidot_a(ix)    
    
    this%grhoa_osc = (phidot_test/this%a_osc)**2/(2.0d0) + this%Vofphi(phi_test,0)
    this%phi2_osc = phi_test**2
    omah2 = (this%grhoa_osc*this%a_osc**3)*grho_to_omega
        
    end subroutine TAxionStandard_AxionDynamicSolver


    subroutine TAxionStandard_Init(this, State)
    class(TAxionStandard), intent(inout) :: this
    class(TCAMBdata), intent(inout), target :: State
    real(dl) ini_phi_1, ini_phi_2, ini_phitest, omah2_1, omah2_2, omah2_test, ommh2

    integer i,j,ix
    real(dl) scaling, delta_a, cH
    integer n_shoot, npoints_log_test, npoints_log_1, npoints_log_2, num_point_converge
    real(dl) a, gpres_ax, phi_ax, phidot_ax, ghroa_t, wa_eff, csquared_ad, delta_phi

    
    select type(State)
    class is (CAMBdata)
        this%State => State
    end select
    this%num_perturb_equations = 2

    this%npoints_coarse = this%max_npoints
    this%log_astart = log(this%astart)
    this%dloga_coarse = (-this%log_astart)/(this%npoints_coarse-1)
            
    this%m = this%m_a*GeV_to_Mpl*(1.e-28_dl)
    this%hubble_ma = this%m_a/isec_to_eV*MPC_in_sec / this%osc_factor
    this%State%m_a = this%m_a
    this%omah2 = this%State%CP%omah2


    if (this%omah2>this%omah2_min) then

        ! Guessing step
        ommh2 = this%State%CP%ombh2 + this%State%CP%omch2
        ini_phi_1 = sqrt((6.d0/9.d0)*(this%omah2/ommh2))*Mpl_GeV/phi_to_GeV !Guess ini_phi
        call this%AxionDynamicSolver(State, ini_phi_1, omah2_1, npoints_log_1)
        
        ! Shooting for correct initial phi
        this%npoints_coarse = int(this%log_acoarse-this%log_astart+1)*300
        this%npoints = 10000 + int(this%osc_factor/10.d0)*1000 !Adapt to oscillation factor
        this%dloga_coarse = (this%log_acoarse-this%log_astart)/(this%npoints_coarse-1)
        this%dloga = (log(2.d0*this%a_osc)-this%log_acoarse)/(this%npoints-1)

        ini_phi_1 = ini_phi_1*sqrt(this%omah2/omah2_1)
        call this%AxionDynamicSolver(State, ini_phi_1, omah2_1, npoints_log_1)

        delta_phi = abs(ini_phi_1*sqrt(this%omah2/omah2_1)-ini_phi_1)

        n_shoot = 0
        do
            if (omah2_1<this%omah2) then
                ini_phi_2 = ini_phi_1 + 2.d0*delta_phi
            else
                ini_phi_2 = ini_phi_1 - 2.d0*delta_phi
            end if
            call this%AxionDynamicSolver(State, ini_phi_2, omah2_2, npoints_log_2)

            n_shoot=n_shoot+1
            if (n_shoot == 10) then
                print*, "ombh2: ", this%State%CP%ombh2
                print*, "omch2: ", this%State%CP%omch2
                print*, "omah2: ", this%omah2
                print*, "me: ", this%State%CP%meR
                print*, "H0: ", this%State%CP%H0
                call MpiStop("Getting stuck at finding two bounds of phi!")
            end if

            if ((omah2_2-this%omah2)*(omah2_1-this%omah2) < 0) then
                ini_phitest = (ini_phi_1+ini_phi_2)/(2.0d0)
                call this%AxionDynamicSolver(State, ini_phitest, omah2_test, npoints_log_test)
                exit
            end if

        end do

        do while (abs(omah2_test-this%omah2)>this%omah2_min)

            ini_phitest = (ini_phi_1+ini_phi_2)/(2.0d0)

            call this%AxionDynamicSolver(State, ini_phitest, omah2_test, npoints_log_test)
            if ((omah2_test-this%omah2)*(omah2_1-this%omah2) < 0) then
                ini_phi_2 = ini_phitest
                omah2_2 = omah2_test
                npoints_log_2 = npoints_log_test
            else
                ini_phi_1 = ini_phitest
                omah2_1 = omah2_test
                npoints_log_1 = npoints_log_test
            end if

            if (npoints_log_1 == npoints_log_2) then
                ini_phitest = ini_phi_1 + (ini_phi_2 - ini_phi_1) * (this%omah2 - omah2_1) / (omah2_2 - omah2_1)
                call this%AxionDynamicSolver(State, ini_phitest, omah2_test, npoints_log_test)
            end if


            if (abs(ini_phi_2-ini_phi_1)<1e-15) then
                
                if (npoints_log_1 /= npoints_log_2) then
                    print*, "Discontinuing function"
                else
                    print*, "Unknown reason"
                end if

                print*, "ombh2: ", this%State%CP%ombh2
                print*, "omch2: ", this%State%CP%omch2
                print*, "omah2: ", this%omah2
                print*, "me: ", this%State%CP%meR
                print*, "H0: ", this%State%CP%H0
                print*, "Two index of phi: ", npoints_log_1, npoints_log_2
                print*, "Two bounds of phi: ", ini_phi_1, ini_phi_2
                print*, "Two bounds of omah2: ", omah2_1, omah2_2
                call MpiStop("Getting stuck at fine-tuning phi!")
            end if

        end do

        this%initial_phi = ini_phitest
        this%State%phi_i = ini_phitest*phi_to_GeV/(1.e17_dl)

    end if


    !!Fill a array!!
    ! ix = npoints_log_test
    ! print*, "Break at ix = ", ix

    ! delta_a = (1._dl - this%State%sampled_a(ix))/(this%npoints-ix)
    ! do i=ix+1, this%npoints
        ! j = i - ix
        ! this%State%sampled_a(i) = this%State%sampled_a(ix) + j*delta_a
    ! end do

    ! open(1, file='out_soundspeed.dat', status='REPLACE')
    ! do i=1, 30000
    !     a = this%State%sampled_a(i)
    !     call this%EquationOfState(a, this%State%phi_a(i), this%State%phidot_a(i), wa_eff)
    !     cH = a * this%Hubble_func(a,this%State%phi_a(i),this%State%phidot_a(i))
    !     call this%AdiabaticSoundSpeed(a, cH, this%State%phi_a(i), this%State%phidot_a(i), csquared_ad)
    !     write(1,*) a, wa_eff, csquared_ad
    ! end do
    ! close(1)

    !!Debugging!!
    ! a = 1.e-12
    ! call this%BackgroundDensityAndPressure(a, ghroa_t, gpres_ax, phi_ax, phidot_ax)
    ! call this%EquationOfState(a, phi_ax, phidot_ax, wa_eff)
    ! call this%AdiabaticSoundSpeed(a, 0._dl, phi_ax, phidot_ax, csquared_ad)

    ! print*, phi_ax, phidot_ax, ghroa_t
    ! print*, wa_eff, csquared_ad

    ! call MpiStop('Finish writing')
    
    end subroutine TAxionStandard_Init


    subroutine TAxionStandard_ReadParams(this, Ini)
    use IniObjects
    class(TAxionStandard) :: this
    class(TIniFile), intent(in) :: Ini
                        
    this%m_a = Ini%Read_Double('m_a')
    this%osc_factor = Ini%Read_Double('osc_factor')

    end subroutine TAxionStandard_ReadParams


    subroutine phi_interp(this,a,aphi,aphidot)
    class(TAxionStandard) :: this
    !Do linear interpolation for background phi and phidot at a (precomputed in Init)
    real(dl) a, aphi, aphidot
    real(dl) delta,a_ix,a_ix_1
    integer ix

    if (a < this%astart) then
        aphi = this%State%phi_a(1)
        aphidot = 0._dl
        return
    elseif ((a >= this%astart) .and. (a < this%a_coarse)) then
        delta = log(a)-this%log_astart
        ix = int(delta/this%dloga_coarse) + 1
    elseif ((a >= this%a_coarse) .and. (a < this%a_osc)) then
        delta = log(a)-this%log_acoarse
        ix = int(delta/this%dloga) + this%npoints_coarse
    else
        call MpiStop("Cannot interpolate after oscillation")
    end if

    a_ix = this%State%sampled_a(ix)
    a_ix_1 = this%State%sampled_a(ix+1)

    if (a >= a_ix .and. a <= a_ix_1) then
        aphi = this%State%phi_a(ix) + (this%State%phi_a(ix+1) - this%State%phi_a(ix))*(a - a_ix)/(a_ix_1 - a_ix);
        aphidot = this%State%phidot_a(ix) + (this%State%phidot_a(ix+1) - this%State%phidot_a(ix))*(a - a_ix)/(a_ix_1 - a_ix);
    else
        print*, "a, a_ix, a_ix_1 ", a, a_ix, a_ix_1
        print*, "ombh2: ", this%State%CP%ombh2
        print*, "omch2: ", this%State%CP%omch2
        print*, "omah2: ", this%omah2
        print*, "me: ", this%State%CP%meR
        print*, "H0: ", this%State%CP%H0
        call MpiStop("Cannot interpolate phi and phidot!")
    end if

    end subroutine phi_interp



    ! subroutine ValsAta(this,a,aphi,aphidot)
    ! class(TAxionStandard) :: this
    ! !Do interpolation for background phi and phidot at a (precomputed in Init)
    ! real(dl) a, aphi, aphidot
    ! real(dl) a0,b0,ho2o6,delta,da
    ! integer ix

    ! if (a < this%astart) then
    !     aphi = this%State%phi_a(1)
    !     aphidot = 0._dl
    !     return
    ! elseif ((a >= this%astart) .and. (a < this%a_osc)) then
    !     delta= log(a)-this%log_astart
    !     ix = int(delta/this%dloga)+1
    ! else
    !     call MpiStop("Cannot interpolate after oscillation")
    ! end if
    
    ! if (ix + 1 > this%npoints_log) then
    !     call MpiStop("Not enough data for interpolation")
    ! end if
    
    ! da = this%State%sampled_a(ix+1) - this%State%sampled_a(ix)
    ! a0 = (this%State%sampled_a(ix+1) - a)/da
    ! b0 = 1 - a0
    ! ho2o6 = da**2/6._dl
    
    ! aphi=b0*this%State%phi_a(ix+1) + a0*(this%State%phi_a(ix)-b0*((a0+1)*this%State%ddphi_a(ix)+(2-a0)*this%State%ddphi_a(ix+1))*ho2o6)
    ! aphidot=b0*this%State%phidot_a(ix+1) + a0*(this%State%phidot_a(ix)-b0*((a0+1)*this%State%ddphidot_a(ix)+(2-a0)*this%State%ddphidot_a(ix+1))*ho2o6)
    
    ! end subroutine ValsAta

    end module AxionStandard
