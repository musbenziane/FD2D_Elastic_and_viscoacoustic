<h2>FD Solver for the Elastic Wave Equation in 2D using Staggered Grids Scheme, as per Virieux 1986</h2>
<br> Changes are still being made ...<br>
Fot input, "parameters.in" in build/. The matlab script generates the receivers positions ASCII file.
<br> Same matlab script can be used to generate shot positions.
<br> OpenMP flag needs to be used during compliation. Otherwise the program will run in serial.
     For CMAKE use the Cmake file as it is. 
<pre>
    !####################################################################################################################
    ! FD2D: Finite difference solver for the elastic wave equation in cartesian coordinates with flat free surface.
    ! NOTE: This program has not been validated.
    ! 
    ! Language: Fortran 90, with parralel impelementation using OpenMP API
    ! 
    ! Author: Mus Benziane
    ! Source: Virieux 1986
    !         other sources [2] are mentioned as comments.
    !
    ! Input file: [Example]
    ! testing              ! prefix of model file names, suffixes: _vp _vs _rho
    ! acqui_src            ! sources positions ascii file name
    ! acqui_rcv            ! receiver positions ascii file name
    ! 1996                 ! zmax
    ! 3196                 ! xmax
    ! 500                  ! nz
    ! 800                  ! nx
    ! 4                    ! dz
    ! 4                    ! dx
    ! 0.0004               ! dt
    ! 4000                 ! nt
    ! 6.                   ! Wavelet's peak frequency
    ! 50                   ! Snapshot interval
    ! 1                    ! [0/1] explosive, otherwise vertical
    ! 1                    ! [0/1] Free Surface - if set to 0, rigid BC are used, i.e no displacement at boundary
    ! 1                    ! [0/1] Absorbing
    ! .65                  ! Att constant for sponge layer
    ! 70                   ! Sponge layer width
    !
    ! -> Acquisition file in ASCII Format: 
    !  1) receivers ASCII file [N: Receivers number]
    !
    !  z1 x1
    !  z2 x2
    !  .  .
    !  zN xN
    !
    !  2) Shot points ASCII file [M: Shot points number]
    !
    !  z1 x1
    !  z2 x2
    !  .  .
    !  zM xM
    ! 
    ! -> Model files in C-Style binary floats [doubles]: Vp, Vs, Rho files are needed.
    !                                                  : For simple models, use create2Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    ! Output files in OUTPUT directory:
    !
    !####################################################################################################################
</pre>

