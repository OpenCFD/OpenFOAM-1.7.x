#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
#
# Script
#     aliases.csh
#
# Description
#     Aliases for working with OpenFOAM
#     Sourced from OpenFOAM-??/etc/cshrc and/or ~/.cshrc
#
#------------------------------------------------------------------------------

# Change compiled version aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wm64 'setenv WM_ARCH_OPTION 64; source $WM_PROJECT_DIR/etc/cshrc'
alias wm32 'setenv WM_ARCH_OPTION 32; source $WM_PROJECT_DIR/etc/cshrc'
alias wmSP 'setenv WM_PRECISION_OPTION SP; source $WM_PROJECT_DIR/etc/cshrc'
alias wmDP 'setenv WM_PRECISION_OPTION DP; source $WM_PROJECT_DIR/etc/cshrc'

# Toggle wmakeScheduler on/off
#  - also need to set WM_HOSTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wmSchedON 'setenv WM_SCHEDULER $WM_PROJECT_DIR/wmake/wmakeScheduler'
alias wmSchedOFF 'unsetenv WM_SCHEDULER'

# Change directory aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~
alias src 'cd $FOAM_SRC'
alias lib 'cd $FOAM_LIB'
alias run 'cd $FOAM_RUN'
alias foam 'cd $WM_PROJECT_DIR'
alias foamsrc 'cd $FOAM_SRC/$WM_PROJECT'
alias foamfv 'cd $FOAM_SRC/finiteVolume'
alias app 'cd $FOAM_APP'
alias util 'cd $FOAM_UTILITIES'
alias sol 'cd $FOAM_SOLVERS'
alias tut 'cd $FOAM_TUTORIALS'
alias foam3rdParty 'cd $WM_THIRD_PARTY_DIR'

# -----------------------------------------------------------------------------
