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
#     etc/settings.sh
#
# Description
#     Startup file for OpenFOAM
#     Sourced from OpenFOAM-??/etc/bashrc
#
#------------------------------------------------------------------------------

# prefix to PATH
_foamAddPath()
{
    while [ $# -ge 1 ]
    do
        export PATH=$1:$PATH
        shift
    done
}

# prefix to LD_LIBRARY_PATH
_foamAddLib()
{
    while [ $# -ge 1 ]
    do
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
        shift
    done
}


# location of the jobControl directory
export FOAM_JOB_DIR=$WM_PROJECT_INST_DIR/jobControl

# wmake configuration
export WM_DIR=$WM_PROJECT_DIR/wmake
export WM_LINK_LANGUAGE=c++
export WM_OPTIONS=$WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_COMPILE_OPTION

# base configuration
export FOAM_APP=$WM_PROJECT_DIR/applications
export FOAM_APPBIN=$WM_PROJECT_DIR/applications/bin/$WM_OPTIONS
export FOAM_LIB=$WM_PROJECT_DIR/lib
export FOAM_LIBBIN=$WM_PROJECT_DIR/lib/$WM_OPTIONS
export FOAM_SRC=$WM_PROJECT_DIR/src

# shared site configuration - similar naming convention as ~OpenFOAM expansion
export FOAM_SITE_APPBIN=$WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION/bin/$WM_OPTIONS
export FOAM_SITE_LIBBIN=$WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION/lib/$WM_OPTIONS

# user configuration
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/applications/bin/$WM_OPTIONS
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/lib/$WM_OPTIONS

# convenience
export FOAM_TUTORIALS=$WM_PROJECT_DIR/tutorials
export FOAM_UTILITIES=$FOAM_APP/utilities
export FOAM_SOLVERS=$FOAM_APP/solvers
export FOAM_RUN=$WM_PROJECT_USER_DIR/run

# add OpenFOAM scripts and wmake to the path
export PATH=$WM_DIR:$WM_PROJECT_DIR/bin:$PATH

_foamAddPath $FOAM_APPBIN $FOAM_SITE_APPBIN $FOAM_USER_APPBIN
 # Make sure to pick up dummy versions of external libraries last
_foamAddLib  $FOAM_LIBBIN/dummy
_foamAddLib  $FOAM_LIBBIN $FOAM_SITE_LIBBIN $FOAM_USER_LIBBIN


# Compiler settings
# ~~~~~~~~~~~~~~~~~
unset compilerBin compilerLib

# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compilerInstall = OpenFOAM | System
compilerInstall=OpenFOAM

case "${compilerInstall:-OpenFOAM}" in
OpenFOAM)
    case "$WM_COMPILER" in
    Gcc)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/gcc-4.3.3/platforms/$WM_ARCH$WM_COMPILER_ARCH
        _foamAddLib $WM_THIRD_PARTY_DIR/mpfr-2.4.1/platforms/$WM_ARCH$WM_COMPILER_ARCH/lib
        _foamAddLib $WM_THIRD_PARTY_DIR/gmp-4.2.4/platforms/$WM_ARCH$WM_COMPILER_ARCH/lib
        ;;
    Gcc43)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/gcc-4.3.3/platforms/$WM_ARCH$WM_COMPILER_ARCH
        _foamAddLib $WM_THIRD_PARTY_DIR/mpfr-2.4.1/platforms/$WM_ARCH$WM_COMPILER_ARCH/lib
        _foamAddLib $WM_THIRD_PARTY_DIR/gmp-4.2.4/platforms/$WM_ARCH$WM_COMPILER_ARCH/lib
        ;;
    Gcc42)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/gcc-4.2.4/platforms/$WM_ARCH$WM_COMPILER_ARCH
        ;;
    esac

    # Check that the compiler directory can be found
    if [ ! -d "$WM_COMPILER_DIR" ]
    then
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.sh:"
        echo "    Cannot find $WM_COMPILER_DIR installation."
        echo "    Please install this compiler version or if you wish to use the system compiler,"
        echo "    change the 'compilerInstall' setting to 'System' in this file"
        echo
    fi

    compilerBin=$WM_COMPILER_DIR/bin
    compilerLib=$WM_COMPILER_DIR/lib$WM_COMPILER_LIB_ARCH:$WM_COMPILER_DIR/lib
    ;;
esac

if [ -d "$compilerBin" ]
then
    _foamAddPath $compilerBin
    _foamAddLib  $compilerLib
fi

unset compilerBin compilerLib compilerInstall

# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH

case "$WM_MPLIB" in
OPENMPI)
    mpi_version=openmpi-1.3.3
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$mpi_version
    export MPI_ARCH_PATH=$MPI_HOME/platforms/$WM_OPTIONS

    # Tell OpenMPI where to find its install directory
    export OPAL_PREFIX=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

SYSTEMOPENMPI)
    mpi_version=openmpi-system

    # Set compilation flags here instead of in wmake/rules/../mplibSYSTEMOPENMPI
    export PINC=`mpicc --showme:compile` 
    export PLIBS=`mpicc --showme:link`
    libDir=`echo "$PLIBS" | sed -e 's/.*-L\([^ ]*\).*/\1/'`

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using system installed MPI:"
        echo "    compile flags : $PINC"
        echo "    link flags    : $PLIBS"
        echo "    libmpi dir    : $libDir"
    fi

    _foamAddLib $libDir


    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

MPICH)
    mpi_version=mpich-1.2.4
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$mpi_version
    export MPI_ARCH_PATH=$MPI_HOME/platforms/$WM_OPTIONS
    export MPICH_ROOT=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

MPICH-GM)
    export MPI_ARCH_PATH=/opt/mpi
    export MPICH_PATH=$MPI_ARCH_PATH
    export MPICH_ROOT=$MPI_ARCH_PATH
    export GM_LIB_PATH=/opt/gm/lib64

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib
    _foamAddLib  $GM_LIB_PATH

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/mpich-gm
    ;;

HPMPI)
    export MPI_HOME=/opt/hpmpi
    export MPI_ARCH_PATH=$MPI_HOME
    export MPICH_ROOT=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin

    case `uname -m` in
    i686)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia32
        ;;

    x86_64)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_amd64
        ;;
    ia64)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia64
        ;;
    *)
        echo Unknown processor type `uname -m` for Linux
        ;;
    esac

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/hpmpi
    ;;

GAMMA)
    export MPI_ARCH_PATH=/usr
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/gamma
    ;;

MPI)
    export MPI_ARCH_PATH=/opt/mpi
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/mpi
    ;;

FJMPI)
    export MPI_ARCH_PATH=/opt/FJSVmpi2
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/mpi

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib/sparcv9
    _foamAddLib  /opt/FSUNf90/lib/sparcv9
    _foamAddLib  /opt/FJSVpnidt/lib
    ;;

QSMPI)
    export MPI_ARCH_PATH=/usr/lib/mpi
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/qsmpi

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib $MPI_ARCH_PATH/lib

    ;;

*)
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/dummy
    ;;
esac

_foamAddLib $FOAM_MPI_LIBBIN


# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
minBufferSize=20000000

if [ "${MPI_BUFFER_SIZE:=$minBufferSize}" -lt $minBufferSize ]
then
    MPI_BUFFER_SIZE=$minBufferSize
fi
export MPI_BUFFER_SIZE


# CGAL library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~
[ -d "$CGAL_LIB_DIR" ] && _foamAddLib $CGAL_LIB_DIR


# Switch on the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if [ -f $FOAM_LIBBIN/libhoard.so ]
#then
#    export LD_PRELOAD=$FOAM_LIBBIN/libhoard.so:$LD_PRELOAD
#fi


# cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
unset _foamAddPath _foamAddLib minBufferSize

# -----------------------------------------------------------------------------
