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
#     aliases.sh
#
# Description
#     Aliases for working with OpenFOAM
#     Sourced from OpenFOAM-??/etc/bashrc and/or ~/.bashrc
#
#------------------------------------------------------------------------------

# Change compiled version aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wm64='export WM_ARCH_OPTION=64; . $WM_PROJECT_DIR/etc/bashrc'
alias wm32='export WM_ARCH_OPTION=32; . $WM_PROJECT_DIR/etc/bashrc'
alias wmSP='export WM_PRECISION_OPTION=SP; . $WM_PROJECT_DIR/etc/bashrc'
alias wmDP='export WM_PRECISION_OPTION=DP; . $WM_PROJECT_DIR/etc/bashrc'

# Toggle wmakeScheduler on/off
#  - also need to set WM_HOSTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alias wmSchedON='export WM_SCHEDULER=$WM_PROJECT_DIR/wmake/wmakeScheduler'
alias wmSchedOFF='unset WM_SCHEDULER'

# Change directory aliases
# ~~~~~~~~~~~~~~~~~~~~~~~~
alias src='cd $FOAM_SRC'
alias lib='cd $FOAM_LIB'
alias run='cd $FOAM_RUN'
alias foam='cd $WM_PROJECT_DIR'
alias foamsrc='cd $FOAM_SRC/$WM_PROJECT'
alias foamfv='cd $FOAM_SRC/finiteVolume'
alias app='cd $FOAM_APP'
alias util='cd $FOAM_UTILITIES'
alias sol='cd $FOAM_SOLVERS'
alias tut='cd $FOAM_TUTORIALS'
alias foam3rdParty='cd $WM_THIRD_PARTY_DIR'

# -----------------------------------------------------------------------------
