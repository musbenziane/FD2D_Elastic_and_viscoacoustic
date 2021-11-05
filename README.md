<h2>FD Solver for the Elastic Wave Equation in 2D using Staggered Grids Scheme, as per Virieux 1986<h2/>



<br>
testing              ! prefix of model file names, suffixes: vp_ vs_ rho_  <br>
acqui_src            ! sources positions ascii file name<br>
acqui_rcv            ! receiver positions ascoo file nam<br>e
1996                 ! zmax<br>
1996                 ! xmax<br>
500                  ! nz<br>
500                  ! nx<br>
4                    ! dz<br>
4                    ! dx<br>
0.0005               ! dt<br>
3000                 ! nt<br>
10.                  ! Wavelet's peak frequency<br>
50                   ! Snapshot interval<br>
1                    ! [0/1] explosive, otherwise v<br>ertical
1                    ! [0/1] Free Surface<br>
1                    ! [0/1] Absorbing<br>
.6                   ! Att constant for sponge layer<br>
30                   ! Sponge layer width<br>
