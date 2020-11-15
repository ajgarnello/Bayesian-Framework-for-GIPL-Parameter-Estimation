! File: main.f90
!
! Author: anthony john garnello
! Attribution: A. Garnello, Projecting permafrost thaw of sub-Arctic tundra with a thermodynamic model calibrated to site measurements JGR Biogeosciences, 2021.

! Original MCMC model script written for Xu, T., White, L., Hui, D., & Luo, Y. (2006). Probabilistic inversion of a terrestrial ecosystem model: Analysis of uncertainty in parameter estimation and model prediction. Global Biogeochemical Cycles, 20(2), 1–15. https://doi.org/10.1029/2005GB002468

! Outline:
! This code was written to allow for an MCMC chain construction using the Metropolis-Hastings Algorithm which was used to construct posterior distributions of soil thermodynamic parameters included in the Geophysical Institute Permafrost Laboratory Model (GIPL 2.0)
! This code is written in Fortran 90, and exists in a series of subroutines each called inside the primary function. Each subroutine contains a brief explanation, but for a full understanding please review the manuscript.

! <-------------------> Begin Code
subroutine GetObsData(O5_cm, O10_cm, O20_cm, O40_cm, O50_cm ,O60_cm, O100_cm, O200_cm, O400_cm, O470_cm,      &
    &       O800_cm, O1500_cm, O2400_cm, OALT)
!   ***            
    implicit none
    character(len=120) observedfile,Day
    integer, parameter :: ilines = 5000
    real O5_cm(ilines),O10_cm(ilines),O20_cm(ilines),O40_cm(ilines),O50_cm(ilines),O60_cm(ilines), O100_cm(ilines),      &
    &    O200_cm(ilines),O400_cm(ilines),O470_cm(ilines),O800_cm(ilines),O1500_cm(ilines),O2400_cm(ilines),OALT(ilines)
    integer istat1, m
    observedfile ='./in/Observed.txt'
    observedfile = trim(observedfile)
    observedfile = adjustl(observedfile)
    open(unit =3, file = observedfile, status = 'old', action = 'read')
    read(3, *)
    m=0
    do
        m = m+1
        read(3, *, IOSTAT = istat1)Day,O5_cm(m),O10_cm(m),O20_cm(m),O40_cm(m),O50_cm(m),O60_cm(m),O100_cm(m),      &
        &   O200_cm(m),O400_cm(m),O470_cm(m),O800_cm(m),O1500_cm(m),O2400_cm(m),OALT(m)
!
        if (istat1 <0)exit
    enddo
    close(3)

end subroutine GetObsData
! subroutine GetObsData opens and reads the file "./in/Obsserved.txt", which is the file containing the measured soil temperature and Active Layer Thickness data from the Eight Mile Lake Experimenal Site. These data strings are assigned to matrices. These are the data that the MCMC method uses for calibration.

subroutine GetParams(parval)
    implicit none
    character(len=140) paramfile
    integer, parameter :: partotal = 110
    real,dimension(partotal):: parval

    paramfile ='./in/SoilProperties.txt'
    paramfile = trim(paramfile)
    paramfile = adjustl(paramfile)
    open(unit =15, file = paramfile, status = 'old', action = 'read')
    read(15, *) ! skip the first line, since comments
    read(15, *) ! skip the second line, since header

    !   this is requiring very particular formatting.
    !           thick,     wvc,      ac,      bc,        cc,    cond_th,  cond_fr,  cap_th,   cap_fr    FIT ! all needed for gipl
    read(15,*)parval(1),parval(2),parval(3),parval(4),parval(5),parval(6),parval(7),parval(8),parval(9),parval(10) ! layer 1 3cm
    read(15,*)parval(11),parval(12),parval(13),parval(14),parval(15),parval(16),parval(17),parval(18),parval(19),parval(20) !2 10cm
    read(15,*)parval(21),parval(22),parval(23),parval(24),parval(25),parval(26),parval(27),parval(28),parval(29),parval(30) !3 19cm
    read(15,*)parval(31),parval(32),parval(33),parval(34),parval(35),parval(36),parval(37),parval(38),parval(39),parval(40) !4 41cm
    read(15,*)parval(41),parval(42),parval(43),parval(44),parval(45),parval(46),parval(47),parval(48),parval(49),parval(50) !5 56cm
    read(15,*)parval(51),parval(52),parval(53),parval(54),parval(55),parval(56),parval(57),parval(58),parval(59),parval(60) !6 76cm
    read(15,*)parval(61),parval(62),parval(63),parval(64),parval(65),parval(66),parval(67),parval(68),parval(69),parval(70) !7 99cm
    read(15,*)parval(71),parval(72),parval(73),parval(74),parval(75),parval(76),parval(77),parval(78),parval(79),parval(80) !8 305cm
    read(15,*)parval(81),parval(82),parval(83),parval(84),parval(85),parval(86),parval(87),parval(88),parval(89),parval(90) !9 513cm
    read(15,*)parval(91),parval(92),parval(93),parval(94),parval(95),parval(96),parval(97),parval(98),parval(99),parval(100) !10 15m 
    read(15,*)parval(101),parval(102),parval(103),parval(104),parval(105),parval(106),      &   
    &   parval(107),parval(108),parval(109),parval(110) !11 100m

    close(15)
    
end subroutine GetParams
! subroutine GetParams is responsible for opening and reading the soil layer parameter file ('./in/SoilProperties.txt') that the GIPL requires. Each soil parameter is then assigned to a specific spot in the parval matrix.

subroutine GetResData(R0_cm, R5_cm, R10_cm, R20_cm, R40_cm, R50_cm, R60_cm, R100_cm, R200_cm,R400_cm,R470_cm,R800_cm,         &
    &       R1500_cm,R2400_cm,sFrost,talik,pfTab,pfBase,pfBase2,talik2,pfTab2,pfTab3,Air_t,Sn_Dep,Sn_Den, nObs)      

    implicit none
    character(len=120) modeloutputfile,Day
    integer, parameter :: ilines = 5000
    real R0_cm(ilines),R5_cm(ilines),R10_cm(ilines),R20_cm(ilines),R40_cm(ilines),R50_cm(ilines),R60_cm(ilines),      &
    &    R100_cm(ilines),R200_cm(ilines),R400_cm(ilines),R470_cm(ilines),R800_cm(ilines),R1500_cm(ilines),R2400_cm(ilines),        &
    &       sFrost(ilines), talik(ilines), pfTab(ilines),pfBase(ilines),pfBase2(ilines),talik2(ilines),pfTab2(ilines),      &
    &   pfTab3(ilines),Air_T(ilines),Sn_Dep(ilines),Sn_Den(ilines)
    integer istat1, m, nObs, i
    integer istat2
    modeloutputfile='./out/Result_daily.txt'
    modeloutputfile = trim(modeloutputfile)
    modeloutputfile = adjustl(modeloutputfile)
    open(unit =2, file = modeloutputfile, status = 'old', action = 'read')
    
    read(2, *) ! skip the first line.
    m=0
    do ! read result file
        m = m+1
        read(2, *, IOSTAT = istat1)Day,Air_T(m),Sn_Dep(m),Sn_Den(m),sFrost(m),talik(m),pfTab(m),     &
        &   pfBase(m),pfBase2(m),talik2(m),pfTab2(m),pfTab3(m),R0_cm(m),R5_cm(m),R10_cm(m),R20_cm(m),R40_cm(m),R50_cm(m),       &
        &   R60_cm(m),R100_cm(m),R200_cm(m),R400_cm(m),R470_cm(m),R800_cm(m),R1500_cm(m),R2400_cm(m)
        nObs = m
        if (istat1 <0)exit
    enddo
     i=0       ! ALT Here: I think this works, Nov 2019
    ! Examine the 5cm temperature series output by the GIPL and assign Active Layer Thickness. Active Layer Thickness is not a direct output by the GIPL, and instead needs to be inferred from the output variables talik and pfTab.
    do i=1,m ! iterate through to turn pfTab3 into ALT column
        if( R5_cm(i) .le. 0) then         
            pfTab3(i) = 0
        endif
        
        if( ( R5_cm(i) .gt. 0) .AND. (talik(i) .le. pfTab(i) )) then         
            pfTab3(i) = talik(i)
        endif
        
         if( (R5_cm(i) .gt. 0) .AND. (talik(i) .gt. pfTab(i))) then
             pfTab3(i) = pfTab(i)
         endif
         
         if( (R5_cm(i) .gt. 0) .AND. (talik2(i) .ne. 0)) then
             pfTab3(i) = talik2(i)
         endif
    enddo
            
    close(2)

end subroutine GetResData
! subroutine GetResData grabs the output of a single GIPL model run, contained in the text file './out/Result_daily.txt'. This file holds daily temperature timeseries for pre-defineds soil layers, including calculated soil ice layer depths. Each data series is assigned into a matrix. Also, Active Layer Thickness is assigned here.

subroutine GetDAcheckbox(DApar,parmin,parmax,DAparfile)
    implicit none
    integer,dimension(110):: DApar
    real,dimension(110):: parmax,parmin
    character(len=50) DAparfile
    DAparfile = './in/EML_da_pars.txt'
    DAparfile=TRIM(DAparfile)

    open(unit =15, file = DAparfile, status = 'old')
    read(15, *) ! skip the first line.
    read(15,*)  ! this line is the header.
    !           thick,  wvc,    ac,         bc,     cc,     cond_th, cond_fr, cap_th, cap_fr    FIT
    read(15,*)DApar(1),DApar(2),DApar(3),DApar(4),DApar(5),DApar(6),DApar(7),DApar(8),DApar(9),DApar(10) ! layer 1
    read(15,*)DApar(11),DApar(12),DApar(13),DApar(14),DApar(15),DApar(16),DApar(17),DApar(18),DApar(19),DApar(20) !2
    read(15,*)DApar(21),DApar(22),DApar(23),DApar(24),DApar(25),DApar(26),DApar(27),DApar(28),DApar(29),DApar(30) !3
    read(15,*)DApar(31),DApar(32),DApar(33),DApar(34),DApar(35),DApar(36),DApar(37),DApar(38),DApar(39),DApar(40) !4
    read(15,*)DApar(41),DApar(42),DApar(43),DApar(44),DApar(45),DApar(46),DApar(47),DApar(48),DApar(49),DApar(50) !5
    read(15,*)DApar(51),DApar(52),DApar(53),DApar(54),DApar(55),DApar(56),DApar(57),DApar(58),DApar(59),DApar(60) !6
    read(15,*)DApar(61),DApar(62),DApar(63),DApar(64),DApar(65),DApar(66),DApar(67),DApar(68),DApar(69),DApar(70) !7
    read(15,*)DApar(71),DApar(72),DApar(73),DApar(74),DApar(75),DApar(76),DApar(77),DApar(78),DApar(79),DApar(80) !8
    read(15,*)DApar(81),DApar(82),DApar(83),DApar(84),DApar(85),DApar(86),DApar(87),DApar(88),DApar(89),DApar(90) !9
    read(15,*)DApar(91),DApar(92),DApar(93),DApar(94),DApar(95),DApar(96),DApar(97),DApar(98),DApar(99),DApar(100) !10
    read(15,*)DApar(101),DApar(102),DApar(103),DApar(104),DApar(105),DApar(106),        &
    &   DApar(107),DApar(108),DApar(109),DApar(110) !11

    
    read(15, *) ! skip the first line
    read(15,*)  ! this line is the header.
    !           thick,  wvc,    ac,         bc,     cc,     cond_th, cond_fr, cap_th, cap_fr
    read(15,*)parmin(1),parmin(2),parmin(3),parmin(4),parmin(5),parmin(6),parmin(7),parmin(8),parmin(9),parmin(10) ! layer 1
    read(15,*)parmin(11),parmin(12),parmin(13),parmin(14),parmin(15),parmin(16),parmin(17),parmin(18),parmin(19),parmin(20) !2
    read(15,*)parmin(21),parmin(22),parmin(23),parmin(24),parmin(25),parmin(26),parmin(27),parmin(28),parmin(29),parmin(30) !3
    read(15,*)parmin(31),parmin(32),parmin(33),parmin(34),parmin(35),parmin(36),parmin(37),parmin(38),parmin(39),parmin(40) !4
    read(15,*)parmin(41),parmin(42),parmin(43),parmin(44),parmin(45),parmin(46),parmin(47),parmin(48),parmin(49),parmin(50) !5
    read(15,*)parmin(51),parmin(52),parmin(53),parmin(54),parmin(55),parmin(56),parmin(57),parmin(58),parmin(59),parmin(60) !6
    read(15,*)parmin(61),parmin(62),parmin(63),parmin(64),parmin(65),parmin(66),parmin(67),parmin(68),parmin(69),parmin(70) !7
    read(15,*)parmin(71),parmin(72),parmin(73),parmin(74),parmin(75),parmin(76),parmin(77),parmin(78),parmin(79),parmin(80) !8
    read(15,*)parmin(81),parmin(82),parmin(83),parmin(84),parmin(85),parmin(86),parmin(87),parmin(88),parmin(89),parmin(90) !9
    read(15,*)parmin(91),parmin(92),parmin(93),parmin(94),parmin(95),parmin(96),parmin(97),parmin(98),parmin(99),parmin(100) !10
    read(15,*)parmin(101),parmin(102),parmin(103),parmin(104),parmin(105),parmin(106),        &
    &   parmin(107),parmin(108),parmin(109),parmin(110) !11
    
    
    read(15, *) ! skip the first line
    read(15,*)  ! this line is the header.
    
    !           thick,  wvc,    ac,         bc,     cc,     cond_th, cond_fr, cap_th, cap_fr
    read(15,*)parmax(1),parmax(2),parmax(3),parmax(4),parmax(5),parmax(6),parmax(7),parmax(8),parmax(9),parmax(10) ! layer 1
    read(15,*)parmax(11),parmax(12),parmax(13),parmax(14),parmax(15),parmax(16),parmax(17),parmax(18),parmax(19),parmax(20) !2
    read(15,*)parmax(21),parmax(22),parmax(23),parmax(24),parmax(25),parmax(26),parmax(27),parmax(28),parmax(29),parmax(30) !3
    read(15,*)parmax(31),parmax(32),parmax(33),parmax(34),parmax(35),parmax(36),parmax(37),parmax(38),parmax(39),parmax(40) !4
    read(15,*)parmax(41),parmax(42),parmax(43),parmax(44),parmax(45),parmax(46),parmax(47),parmax(48),parmax(49),parmax(50) !5
    read(15,*)parmax(51),parmax(52),parmax(53),parmax(54),parmax(55),parmax(56),parmax(57),parmax(58),parmax(59),parmax(60) !6
    read(15,*)parmax(61),parmax(62),parmax(63),parmax(64),parmax(65),parmax(66),parmax(67),parmax(68),parmax(69),parmax(70) !7
    read(15,*)parmax(71),parmax(72),parmax(73),parmax(74),parmax(75),parmax(76),parmax(77),parmax(78),parmax(79),parmax(80) !8
    read(15,*)parmax(81),parmax(82),parmax(83),parmax(84),parmax(85),parmax(86),parmax(87),parmax(88),parmax(89),parmax(90) !9
    read(15,*)parmax(91),parmax(92),parmax(93),parmax(94),parmax(95),parmax(96),parmax(97),parmax(98),parmax(99),parmax(100) !10
    read(15,*)parmax(101),parmax(102),parmax(103),parmax(104),parmax(105),parmax(106),        &
    &   parmax(107),parmax(108),parmax(109),parmax(110) !11

    
    close(15)
end subroutine GetDAcheckbox
! subroutine GetDAcheckbox opens and reads the text file './in/EML_da_pars.txt', which contains three matrices defining A) which GIPL parameters are actively being parameterized using the MCMCM method, B) the lower end of their uniform prior distributions and C) the upper end of their uniform prior distributions. Each value corresponds to a specific parameter within a specific layer and is assigned a spot in data matrix.

subroutine param_write(parval)
! this is to intake the generated parameters via M-H, and write them to a datafile that matches the SoilProperties.txt
    integer, parameter :: partotal = 110
    integer isimu
    real,dimension(partotal):: parmin,parmax,parval
    character(len=140) param_file_out, paraestfile
    param_file_out = "./in/SoilProperties.txt"
    open(unit =15,  file = param_file_out, status = "old") 
!    print*,parval
    write(15, *) "11 // number of soil layers//"
    write(15, *) "  thick     wvc    ac       bc       cc    cond_th  cond_fr         cap_th                 cap_fr      FIT "
    write(15, "(F10.2, 1x, F7.4, 1x, F7.3, 1x, F7.3, 1x, F7.3, 1x, F7.2, 1x, F7.4, 1x, F20.0, 1x, F20.0, 1x, F7.2, 1x)")parval
    close(15)
end subroutine param_write
! subroutine paramwrite is responsible for writing a new "./in/SoilProperties.txt" file containing the parameter estimates output by the M-H algorithm inside the matrix parval. The output format of this file must match the original file readable by the GIPL model.


subroutine Gipl_Wrap()
    !           *** variable definition for reading observed values and init parameters:
    integer, parameter :: ilines = 5000 ! over-estimate the number of lines in the results text file
    integer, parameter :: jlines = 11 ! determine the number of soil layers
    integer, parameter :: partotal = 110 ! determine the total number of parameters inside the GIPL (11 soil layers, 10 parameters per layer = 110)

!           *** variable definition for reading observed calibration data **
    character(len=120) outdir, outfile, my_fmt
    real O5_cm(ilines),O10_cm(ilines),O20_cm(ilines),O40_cm(ilines),O50_cm(ilines),O60_cm(ilines), O100_cm(ilines),      &
    &    O200_cm(ilines),O400_cm(ilines),O470_cm(ilines),O800_cm(ilines),O1500_cm(ilines),O2400_cm(ilines),OALT(ilines)
    
!           *** variable definition for reading model result **
    character(len=120) modeloutputfile,Day,covfile
    real R0_cm(ilines),R5_cm(ilines),R10_cm(ilines),R20_cm(ilines),R40_cm(ilines),R50_cm(ilines),R60_cm(ilines),      &
    &    R100_cm(ilines),R200_cm(ilines),R400_cm(ilines),R470_cm(ilines),R800_cm(ilines),R1500_cm(ilines),R2400_cm(ilines),        &
    &       sFrost(ilines), talik(ilines), pfTab(ilines),pfBase(ilines),pfBase2(ilines),talik2(ilines),sFrost3(ilines),     &
    &   pfTab3(ilines), Air_T(ilines), Sn_Dep(ilines), Sn_Den(ilines),pfTab2(ilines)
!            *** variable definition for cost function:  
    logical, parameter :: do_soilt_da = .True.
    real J_dsoilt_5,j_ds_5,J_dsoilt_10,j_ds_10,J_dsoilt_20,j_ds_20,J_dsoilt_40,j_ds_40,J_dsoilt_50,j_ds_50,J_dsoilt_60,j_ds_60
    real J_dsoilt_100,j_ds_100,J_dsoilt_200,j_ds_200,J_dsoilt_400,j_ds_400,J_dsoilt_470,j_ds_470,J_dsoilt_800,j_ds_800
    real J_dsoilt_1500,j_ds_1500,J_dsoilt_2400,j_ds_2400, J_dsoilt_ALT, J_ds_ALT
    real J_new, J_last, r_num, Jout, J_newA, J_newB
     
    !       *** variable definition for reading the Data-Assimilation file
    integer,dimension(partotal):: DApar
    real,dimension(partotal):: parmin,parmax,parval
    character(len=50) DAparfile, paraestfile
    
    !       *** variable definition for Parameter generation:
    real paraest(11,40000)
    integer npara, upgraded, tmp_up, new, isimu, j, i, gaussgencount
    
    real, allocatable :: coef(:), coefac(:), coefnorm(:)
    real, allocatable :: coefmax(:),coefmin(:)
    real, allocatable :: gamma(:,:),gamnew(:,:)
    real, allocatable :: coefhistory(:,:)
    integer,allocatable :: coefindex(:)

    integer k1,k2,rejet,paraflag,k3
    integer, parameter :: nc=1000
    integer, parameter :: ncov=1500
    
    integer nObs
!            *** END variable definition ***  
    
    modeloutputfile='./out/Result_daily.txt'
    outdir = "./out/"
    DAparfile='./in/EML_da_pars.txt'
    paraestfile = './out/Paraest.txt'
    print *, 'Inside the GIPL Wrap Function...'
    
    !       *** read-in the observed data
    call GetObsData(O5_cm, O10_cm, O20_cm, O40_cm, O50_cm ,O60_cm, O100_cm, O200_cm, O400_cm, O470_cm,      &
    &       O800_cm, O1500_cm, O2400_cm, OALT)
    
    !       *** read-in the initial parameter values
    call GetParams(parval)
!       *** read-in the DA checkbox
    call GetDAcheckbox(DApar,parmin,parmax,DAparfile)
    npara=sum(DApar) ! determine how many parameters will be included in MCMC

    allocate(coef(npara),coefac(npara),coefnorm(npara))
    allocate(coefindex(npara))
    allocate(coefmax(npara),coefmin(npara))
    allocate(gamma(npara,npara),gamnew(npara,npara))
    allocate(coefhistory(ncov,npara))

    J_last=9000000.0 ! assign an unreasonably high J_last to start

    upgraded=0
    new=0
    k3=0
    j=0
    gaussgencount=0
    ! create the paraest file
    paraestfile = './out/Paraest.txt'
    open(15, file = paraestfile, action = "write", status = "old")
    write(15, 101) !"iteration thick wvc ac bc cc cond_th cond_fr cap_th cap_fr FIT"
101 format("upgraded",12x, "cost", 8x, "thick",4x, "wvc",6x, "ac",8x, "bc",6x, "cc",6x, "cond_th",2x, "cond_fr", 12x,     &
    &   "cap_th", 16x, "cap_fr",8x,"FIT")
    close(15)
    
    do i=1, partotal ! this moves all across DApar, which is the full size of the parameter space.

       if (DApar(i).eq. 1) then ! this pulls out the parameters ONLY for those indicated in the DApar matrix
           j=j+1
           coef(j)=parval(i) ! define initial value for parameters, which are equal to those assigned in the ./in/SoilProperties file.
           coefindex(j)=i           
           coefmin(j)=parmin(i)
           coefmax(j)=parmax(i)
       endif
    enddo
    
    covexist=0 ! assign to 1 if using a pre-defined covariance matrix
    if(covexist.eq.1)then      ! If prior covariance exists, read from file
        write(covfile,"(A120,A15)") "./out/covariance.txt"
        covfile = trim(covfile)
        covfile = adjustl(covfile)
        call getCov(gamma,covfile,npara)

        call racine_mat(gamma,gamnew,npara)
        gamma=gamnew
        do k1=1,npara
           coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1))
            coefac(k1)=coefnorm(k1)
        enddo
    else
        coefac=coef
    endif

    fact_rejet=2.4/sqrt(real(npara))
    search_length=0.05
    rejet = 0

    ! enter the iterative section of the model
    do isimu = 1, 150000
        ! generate parameters ! 
61 continue
        gaussgencount=0
        if(covexist.eq.1)then
            paraflag=1
            gaussgencount=0
            do while(paraflag.gt.0)
                call gengaussvect(fact_rejet*gamma,coefac,coefnorm,npara)
                paraflag=0
                do k1=1,npara ! 
                    if(coefnorm(k1).lt.0. .or. coefnorm(k1).gt.1.)then
                    paraflag=paraflag+1 ! paraflag shows which parameter is out of bounds
                    
                    write(*,*)'out of range',k1,coefnorm(k1) ! print out to screen
                    endif
                enddo
            enddo
            do k1=1,npara
                coef(k1)=coefmin(k1)+coefnorm(k1)*(coefmax(k1)-coefmin(k1))
            enddo
        else

            call coefgenerate(coefac,coefmax,coefmin,coef,search_length,npara)
        endif
                ! update parameters with the newly generated values
        do k1=1,npara
            parval(coefindex(k1))=coef(k1)
        enddo

        call param_write(parval) ! write a new SoilProperties.txt file
        call SYSTEM("./gipl2cpp"); !call the GIPL model to run using the new SoilProperties
        call GetResData(R0_cm, R5_cm, R10_cm, R20_cm, R40_cm, R50_cm, R60_cm, R100_cm, R200_cm,R400_cm,R470_cm,R800_cm,         &
    &       R1500_cm,R2400_cm,sFrost,talik,pfTab,pfBase,pfBase2,talik2,pfTab2,pfTab3,Air_t,Sn_Dep,Sn_Den, nObs)      
        
!       Begin MCMC ********************
! Likelihood function definition
        tmp_up=upgraded
         if (do_soilt_da) then
            
            !initialization of model variables
            J_dsoilt_5=0.0
            J_dsoilt_10=0.0
            J_dsoilt_20=0.0
            J_dsoilt_40=0.0
            J_dsoilt_ALT=0.0
            J_dsoilt_100=0.0
            J_dsoilt_200=0.0
            J_dsoilt_400=0.0
            J_dsoilt_800=0.0
            J_dsoilt_1500=0.0

            j_ds_5=0
            j_ds_10=0
            j_ds_20=0
            j_ds_40=0
            j_ds_ALT=0
            j_ds_100=0
            j_ds_200=0
            j_ds_400=0
            j_ds_800=0
            j_ds_1500=0

            ! grab and isolate the calibration values
            do j=1,nObs ! j is daily
                if(O5_cm(j).gt.-999.)then
                    j_ds_5 = j_ds_5+1 ! this counts the number of observations gt -999
                    dObsSim_5=R5_cm(j) - O5_cm(j) ! find the difference between observed and modeled temperature at time j
                    J_dsoilt_5=J_dsoilt_5+(dObsSim_5*dObsSim_5)/(2*0.9*0.9)! running sum of difference, normalizing to 2*sigma^2
                endif
                if(O10_cm(j).gt.-999.)then
                   j_ds_10=j_ds_10+1 
                   dObsSim_10=R10_cm(j)-O10_cm(j)
                   J_dsoilt_10=J_dsoilt_10+(dObsSim_10*dObsSim_10)/(2*0.9*0.9)
                endif
                if(O20_cm(j).gt.-999.)then
                   j_ds_20=j_ds_20+1 
                   dObsSim_20=R20_cm(j)-O20_cm(j)
                   J_dsoilt_20=J_dsoilt_20+(dObsSim_20*dObsSim_20)/(2*0.9*0.9)
                endif
                if(O40_cm(j).gt.-999.)then
                   j_ds_40=j_ds_40+1 
                   dObsSim_40=R40_cm(j)-O40_cm(j)
                   J_dsoilt_40=J_dsoilt_40+(dObsSim_40*dObsSim_40)/(2*0.8*0.8)
                endif
                if(OALT(j).gt.-999.)then
                   j_ds_ALT=j_ds_ALT+1 
                   dObsSim_ALT=pfTab3(j)-OALT(j)
                   J_dsoilt_ALT=J_dsoilt_ALT+(dObsSim_ALT*dObsSim_ALT)/(2*0.1*0.1)
                endif
                if(O100_cm(j).gt.-999.)then
                   j_ds_100=j_ds_100+1 
                   dObsSim_100=R100_cm(j)-O100_cm(j)
                   J_dsoilt_100=J_dsoilt_100+(dObsSim_100*dObsSim_100)/(2*0.3*0.3) ! changed from 0.75
                endif
                if(O200_cm(j).gt.-999.)then
                   j_ds_200=j_ds_200+1 
                   dObsSim_200=R200_cm(j)-O200_cm(j)
                   J_dsoilt_200=J_dsoilt_200+(dObsSim_200*dObsSim_200)/(2*0.15*0.15)
                endif
                if(O400_cm(j).gt.-999.)then
                   j_ds_400=j_ds_400+1 
                   dObsSim_400=R400_cm(j)-O400_cm(j)
                   J_dsoilt_400=J_dsoilt_400+(dObsSim_400*dObsSim_400)/(2*0.1*0.1)
                endif
                if(O800_cm(j).gt.-999.)then
                   j_ds_800=j_ds_800+1 
                   dObsSim_800=R800_cm(j)-O800_cm(j)
                   J_dsoilt_800=J_dsoilt_800+(dObsSim_800*dObsSim_800)/(2*0.05*0.05)
                endif
                if(O1500_cm(j).gt.-999.)then
                   j_ds_1500=j_ds_1500+1 
                   dObsSim_1500=R1500_cm(j)-O1500_cm(j)
                   J_dsoilt_1500=J_dsoilt_1500+(dObsSim_1500*dObsSim_1500)/(2*0.08*0.08)
                endif

           enddo  
            J_dsoilt_5=J_dsoilt_5/real(j_ds_5) ! this creates an average difference to account for differences in available data
            write(*,*)'5 Cost: ', J_dsoilt_5 ! print to screen for quick-check
            J_dsoilt_10=J_dsoilt_10/real(j_ds_10)
            write(*,*)'10 Cost: ', J_dsoilt_10
            J_dsoilt_20=J_dsoilt_20/real(j_ds_20)
            write(*,*)'20 Cost: ', J_dsoilt_20
            J_dsoilt_40=J_dsoilt_40/real(j_ds_40)
            write(*,*) '40 Cost: ', J_dsoilt_40
            J_dsoilt_ALT=J_dsoilt_ALT/real(j_ds_ALT)
            write(*,*)'ALT Cost: ', J_dsoilt_ALT
            J_dsoilt_100=J_dsoilt_100/real(j_ds_100)
            write(*,*)'100 Cost: ', J_dsoilt_100
            J_dsoilt_200=J_dsoilt_200/real(j_ds_200)
            write(*,*)'200 Cost: ', J_dsoilt_200
            J_dsoilt_400=J_dsoilt_400/real(j_ds_400)
            write(*,*)'400 Cost: ', J_dsoilt_400
            J_dsoilt_800=J_dsoilt_800/real(j_ds_800)
            write(*,*)'800 Cost: ', J_dsoilt_800
            J_dsoilt_1500=J_dsoilt_1500/real(j_ds_1500)
            write(*,*)'1500 Cost: ', J_dsoilt_1500

            ! sum the total disaggreement, or cost, of the data series
            J_newA=J_dsoilt_5+J_dsoilt_10+J_dsoilt_20+J_dsoilt_40+J_dsoilt_ALT+J_dsoilt_100+J_dsoilt_200

           J_newB=J_dsoilt_400+J_dsoilt_800+J_dsoilt_1500
            
           J_new = J_newA+J_newB
          endif

          delta_J=J_new-J_last    ! if delta_j < 0, means new cost was lower than previous cost.
          write(*,*) 'J_new: ', J_new, '    delta_J: ', delta_J
          call random_number(r_num)

            if(ISNAN(J_new))then ! Most likely an input read problem (faulty GIPL results file)
                write(*,*)'NaN return, upgraded', upgraded
                goto 61 ! rinse and repeat
            endif

            if(J_new<-2e+09)then ! Most likely an input read problem (faulty GIPL results file)
                write(*,*)'unbelievably low J_new, returning to 61, upgraded', upgraded
                goto 61 ! rinse and repeat
            endif

            if(delta_J.eq.0)then ! Most likely an input read problem (faulty GIPL results file)
                write(*,*)'delta_J equals 0, returning to 61, upgraded', upgraded
                goto 61 ! rinse and repeat
            endif

            ! Metropolis-Hastings Test Criteria
            if(AMIN1(1.0,exp(-delta_J)).gt.r_num)then
                upgraded=upgraded+1 
                J_last=J_new
            endif
            
            ! *** end of Likelihood function definition ***
            if(upgraded.gt.tmp_up)then ! true if newly generated parameters produced a lower cost than the previously generated parameters
                new=new+1
                open(15, file = paraestfile, position = "append", action = "write", status = "old") ! Update the ./in/SoilProperties.txt file with the accepted parameters.
                
                write(15, 201) upgraded, J_new, parval(1:10) !"iteration thick wvc ac bc cc cond_th cond_fr cap_th cap_fr IT"
                write(15, 201) upgraded, J_new, parval(11:20)
                write(15, 201) upgraded, J_new, parval(21:30)
                write(15, 201) upgraded, J_new, parval(31:40)
                write(15, 201) upgraded, J_new, parval(41:50)
                write(15, 201) upgraded, J_new, parval(51:60)
                write(15, 201) upgraded, J_new, parval(61:70)
                write(15, 201) upgraded, J_new, parval(71:80)
                write(15, 201) upgraded, J_new, parval(81:90)
                write(15, 201) upgraded, J_new, parval(91:100)
                write(15, 201) upgraded, J_new, parval(101:110)

            201 format(I6,1x, F20.5,1x, F10.4,2x, 4(F7.4,2x), 2(F7.4,2x), 2(F20.1,2x), F7.4,2x)
                close(15)


                if(covexist.eq.1)then
                    coefac=coefnorm
                    coefhistory(new,:)=coefnorm
                else
                    coefac=coef
                    do k1=1,npara
                        coefnorm(k1)=(coef(k1)-coefmin(k1))/(coefmax(k1)-coefmin(k1)) ! normalized coefficients.
                    enddo
                endif
            coefhistory(new,:)=coefnorm ! store normalized coefficients in coefhistory.
            
            if(new.ge.ncov)new=0 !
        else
            rejet=rejet+1
        endif      
        
        write(*,*)'isimu',isimu,'upgraded',upgraded
        
        
        if(covexist.eq.1)then
            if(mod(isimu,nc).eq.0)then
                if ((1. - real(rejet)/real(nc)) < 0.23) then
                    fact_rejet = fact_rejet*0.9
                else
                    if ((1. - real(rejet)/real(nc)) > 0.44) then
                    fact_rejet = fact_rejet * 1.1
                    endif
                endif
            rejet=0
            write(*,*)'search length is', search_length
            endif
        else
            if(mod(isimu,nc).eq.0)then
                if(real(upgraded)/real(isimu) .lt. 0.23)then ! this is checking for the acceptance rate and changing the search length accordingly
                    search_length=search_length*0.9
                else
                    if(real(upgraded)/real(isimu) .gt. 0.44)then
                        search_length=search_length*1.1
                    endif
                endif
                rejet=0
                write(*,*)'search length is', search_length
            endif
        endif
    
        if(covexist.eq.0 .and. mod(upgraded,ncov).eq.0)then
            covexist=1
            coefac=coefnorm ! coefnorm has been storing the normalized accepted coefficients
            call varcov(coefhistory,gamnew,npara,ncov)
            if (.not.(all(gamnew==0.))) then
                gamma=gamnew

                call racine_mat(gamma,gamnew,npara)
                gamma=gamnew
            endif
        endif
        
	if (mod(upgraded,ncov).eq.0 .and. covexist.eq.1) then
            call varcov(coefhistory,gamnew,npara,ncov)
            if (.not.(all(gamnew==0.))) then
                gamma=gamnew

                call racine_mat(gamma,gamnew,npara)
                gamma=gamnew
            endif
	endif

    enddo
        write(outfile,"(A120,A21)")"./out/covariance_temp.txt"
        outfile = trim(outfile)
        outfile = adjustl(outfile)
        open(72,file=outfile)
        do i=1,npara
            write(72,*) (gamma(j,i),j=1,npara)

        enddo
        close(72) 

    
end subroutine Gipl_Wrap


program main
    implicit none
!       *** variable definition for reading model result:
    character(len=120) modeloutputfile,Day,Air_T,Sn_Dep,Sn_Den,sFrost,talik,pfTab,pfBase,   &
    &   initialParamfile, SoilProperties, paraestfile, outdir
    
    !       *** variable definition for Parameter generation:
    real paraest(11,40000)
    integer npara ! number of parameters to be estimated
    integer,parameter :: partotal= 110
    integer,dimension(partotal):: DApar
    real, allocatable :: coef(:), coefac(:), coefnorm(:)
    real, allocatable :: coefmax(:),coefmin(:)
    real, allocatable :: gamma(:,:),gamnew(:,:)
    integer,allocatable :: coefindex(:)
    
    integer k1,k2,rejet,paraflag,k3
    integer, parameter :: nc=500
    integer, parameter :: ncov=1000
    real,dimension(partotal) :: parval,parmin,parmax
    
    !       *** variable definition for Parameter generation:
    integer, parameter :: ilines = 5000
    real R0_cm(ilines),R5_cm(ilines),R10_cm(ilines),R20_cm(ilines),R40_cm(ilines),R50_cm(ilines),R60_cm(ilines),      &
    &    R100_cm(ilines),R200_cm(ilines),R400_cm(ilines),R470_cm(ilines),R800_cm(ilines),R1500_cm(ilines),R2400_cm(ilines)
    modeloutputfile='./out/Result_daily.txt'
    initialParamfile='./in/Initial_Properties.txt'
    SoilProperties = './in/SoilProperties.txt'
    paraestfile = './out/Paraest.txt'
    outdir = './out'
!       *** 
! Initialize the soilPropertiex.txt file using the Initial_Properties.txt file
     call SYSTEM('cp ./in/Initial_Properties.txt ./in/SoilProperties.txt')
    print *, "starting main function!"
!       *** Call the wrapper function
    call Gipl_Wrap()
     
end program main

! *********************************** functions:

!******************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Square root of a matrix							  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine racine_mat(M, Mrac,npara)

    integer npara,i
    real M(npara,npara),Mrac(npara,npara)
    real valpr(npara),vectpr(npara,npara)

    Mrac=0. ! this resets gamnew to all 0's

    call jacobi(M,npara,npara,valpr,vectpr,nrot)
    do i=1,npara
	if(valpr(i).ge.0.) then
            Mrac(i,i)=sqrt(valpr(i)) ! this is filling the diagonal, I.E. The variance vector of covar. mat.
	else
            print*, 'WARNING!!! Square root of the matrix is undefined.'
            print*, ' A negative eigenvalue has been set to zero - results may be wrong'
 
            Mrac=M
            return ! this ejects out of the subroutine, preventing fill of the rest of the diagonal
	endif

    enddo
    Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr)) ! this is a symmetric matrix.
end subroutine racine_mat      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extraction of the eigenvalues and the eigenvectors !!
!! of a matrix (Numerical Recipes)					  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE jacobi(a,n,np,d,v,nrot) ! a is M, which is gamma, d is valpr()), v is vectpr(). 
INTEGER :: n,np,nrot
REAL :: a(np,np),d(np),v(np,np)
INTEGER, PARAMETER :: NMAX=500
INTEGER :: i,ip,iq,j
REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

do ip=1,n
	do iq=1,n
		v(ip,iq)=0.
	end do
	v(ip,ip)=1.
end do
 
do ip=1,n
	b(ip)=a(ip,ip)
	d(ip)=b(ip)
	z(ip)=0.
end do

nrot=0
do i=1,50
	sm=0.
	do ip=1,n-1
		do iq=ip+1,n
			sm=sm+abs(a(ip,iq))
		end do
	end do

	if(sm.eq.0.)then

            return
        endif
	if(i.lt.4)then
		tresh=0.2*sm/n**2
	else
		tresh=0.
	endif
	do ip=1,n-1
		do iq=ip+1,n
			g=100.*abs(a(ip,iq))
			if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
                                a(ip,iq)=0.
                        else if(abs(a(ip,iq)).gt.tresh)then
				h=d(iq)-d(ip)
				if(abs(h)+g.eq.abs(h))then
                                        t=a(ip,iq)/h
				else

					theta=0.5*h/a(ip,iq)
					t=1./(abs(theta)+sqrt(1.+theta**2))
					if(theta.lt.0.) then
						t=-t
					endif
				endif
				c=1./sqrt(1+t**2)
				s=t*c
				tau=s/(1.+c)
				h=t*a(ip,iq)
				z(ip)=z(ip)-h
				z(iq)=z(iq)+h
				d(ip)=d(ip)-h
				d(iq)=d(iq)+h
				a(ip,iq)=0.
				do j=1,ip-1
					g=a(j,ip)
					h=a(j,iq)
					a(j,ip)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				end do
				do j=ip+1,iq-1
					g=a(ip,j)
					h=a(j,iq)
					a(ip,j)=g-s*(h+g*tau)
					a(j,iq)=h+s*(g-h*tau)
				end do
				do j=iq+1,n
					g=a(ip,j)
					h=a(iq,j)
					a(ip,j)=g-s*(h+g*tau)
					a(iq,j)=h+s*(g-h*tau)
				end do
				do j=1,n
					g=v(j,ip)
					h=v(j,iq)
					v(j,ip)=g-s*(h+g*tau)
					v(j,iq)=h+s*(g-h*tau)
				end do
				nrot=nrot+1
			endif
		end do
	end do
	do ip=1,n
		b(ip)=b(ip)+z(ip)
		d(ip)=b(ip)
		z(ip)=0.
	end do
end do
print*, 'too many iterations in jacobi' 
return
END subroutine jacobi


!===================================================
!       generate new coefficents
        subroutine coefgenerate(coefac,coefmax,coefmin,coef,search_length,npara)
        integer npara
        real coefac(npara),coefmax(npara),coefmin(npara),coef(npara)
        real r,coefmid,random_harvest, cond_th(7), cap_th(10), cond_fr(7), cap_fr(10)
        real cond_th_min(7), cap_th_min(10), cond_fr_min(7), cap_fr_min(10)
        real cond_th_max(7), cap_th_max(10), cond_fr_max(7), cap_fr_max(10)
        integer i
        real search_length

        do i=1,npara
999         continue
            CALL random_number(random_harvest)
            r=random_harvest-0.5

            coef(i)=coefac(i)+r*(coefmax(i)-coefmin(i))*search_length
            if(coef(i).gt.coefmax(i).or.coef(i).lt.coefmin(i))goto 999
        enddo

        ! When generating parameters, it is crucial to restrict the values of the frozen and thawed conductivities and capacities of soil layers. Realistically, the frozen thermal conductivity will always be higher than the same layer when thawed. And, the frozen thermal capacity will always be lower than the capacity for the same layer when thawed.
        ! see the DA file to ensure that the variables pointed at here always select the conductivites
        do i=2,32,5
            do while(coef(i).gt.coef(i+1))
998         continue
            CALL random_number(random_harvest)
            r=random_harvest-0.5

            coef(i)=coefac(i)+r*(coefmax(i)-coefmin(i))*search_length ! change the thawed conductivities
            if(coef(i).gt.coefmax(i).or.coef(i).lt.coefmin(i))goto 998
            CALL random_number(random_harvest)
            r=random_harvest-0.5
            coef(i+1)=coefac(i+1)+r*(coefmax(i+1)-coefmin(i+1))*search_length ! change the frozen conductivities
            
            if(coef(i+1).gt.coefmax(i+1).or.coef(i+1).lt.coefmin(i+1))goto 998
            enddo
        enddo
        
        do i=4,34,5
            do while(coef(i).lt.coef(i+1))
997         continue                
            CALL random_number(random_harvest)
            r=random_harvest-0.5

            coef(i)=coefac(i)+r*(coefmax(i)-coefmin(i))*search_length ! change the thawed capacities
            if(coef(i).gt.coefmax(i).or.coef(i).lt.coefmin(i))goto 997
            CALL random_number(random_harvest)
            r=random_harvest-0.5

            coef(i+1)=coefac(i+1)+r*(coefmax(i+1)-coefmin(i+1))*search_length ! change the frozen capacities
            if(coef(i+1).gt.coefmax(i+1).or.coef(i+1).lt.coefmin(i+1))goto 997
            enddo
        enddo

        
        return
        end
!============
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random vector from a multivariate  !!
!! normal distribution with mean zero and covariance  !!
!! matrix gamma.									  !!
!! Beware!!! In order to improve the speed of the	  !!
!! algorithms, the subroutine use the Square root	  !!
!! matrix of gamma									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! what does this subroutine do?
subroutine gengaussvect(gamma_racine,xold,xnew,npara)
integer npara
real gamma_racine(npara,npara)
real x(npara),xold(npara),xnew(npara)

do i=1,npara
    x(i)=rangauss(25)
enddo

x = matmul(gamma_racine, x)
xnew = xold + x
end subroutine gengaussvect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Generation of a random number from a standard      !!
!! normal distribution. (Numerical Recipes)           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! but why 25
function rangauss(idum)


integer idum
real v1, v2, r, fac, gset
real r_num

data iset/0/

if(iset==0) then
1	CALL random_number(r_num)
            v1=2.*r_num-1
            CALL random_number(r_num)
	v2=2.*r_num-1
	r=(v1)**2+(v2)**2
	if(r>=1) go to 1
	fac=sqrt(-2.*log(r)/r)
	gset=v1*fac
	rangauss=v2*fac
	iset=1
else
	rangauss=gset
	iset=0
end if

return
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! variance matrix of a matrix of data				  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine varcov(tab,varcovar,npara,ncov)

integer npara,ncov
real tab(ncov,npara),tab2(ncov,npara)
real varcovar(npara,npara)

call centre(tab,tab2,npara,ncov)

varcovar = matmul(transpose(tab2), tab2)*(1./real(ncov)) 

end subroutine varcov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute the centered matrix, ie. the matrix minus  !!
!! the column means									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine centre(mat,mat_out,npara,ncov)

    integer npara,ncov
    real mat(ncov,npara),mat_out(ncov,npara)
    real mean
           
do i=1,npara
    mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
enddo

end subroutine centre

! **********************************************************    
subroutine getCov(gamma,covfile,npara)
    implicit none
    integer npara,i,k
    real gamma(npara,npara)    
    character(len=80) covfile

    open(14,file=covfile,status='old')

    do i=1,npara
        read (14,*)(gamma(i,k),k=1,npara)
    enddo  
    return
end subroutine getCov        


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mean of a vector									  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function mean(tab,ncov)
    integer ncov
real tab(ncov)
real mean,mean_tt
mean_tt=0.
do i=1,ncov	
mean_tt=mean_tt+tab(i)/real(ncov)
enddo
mean=mean_tt
End Function

