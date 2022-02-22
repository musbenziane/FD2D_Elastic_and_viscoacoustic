#!/bin/bash



mpirun -n 1 FD2DVAC  << OUT

testing              ! prefix of model file names, suffixes: _vp _vs _rho
acqui_src            ! sources positions ascii file name
acqui_rcv            ! receiver positions ascii file name
200                  ! nz
200                  ! nx
2.5                  ! dz, dx
0.00012              ! dt
5200                 ! nt
15.                  ! Wavelet's peak frequency
20                   ! Snapshot interval
1                    ! is_FS
1                    ! is_PML
60                   ! nPML

OUT
