#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
# File
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

# prefix to MANPATH
_foamAddMan()
{
    while [ $# -ge 1 ]
    do
        export MANPATH=$1:$MANPATH
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

_foamAddPath $FOAM_USER_APPBIN:$FOAM_SITE_APPBIN:$FOAM_APPBIN
 # Make sure to pick up dummy versions of external libraries last
_foamAddLib  $FOAM_USER_LIBBIN:$FOAM_SITE_LIBBIN:$FOAM_LIBBIN:$FOAM_LIBBIN/dummy

# Compiler settings
# ~~~~~~~~~~~~~~~~~
unset gcc_version gmp_version mpfr_version mpc_version
unset MPFR_ARCH_PATH

# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compilerInstall = OpenFOAM | system
: ${compilerInstall:=OpenFOAM}

case "${compilerInstall:-OpenFOAM}" in
OpenFOAM)
    case "$WM_COMPILER" in
    Gcc)
        gcc_version=gcc-4.4.3
        gmp_version=gmp-5.0.1
        mpfr_version=mpfr-2.4.2
        ;;
    Gcc45)
        gcc_version=gcc-4.5.0
        gmp_version=gmp-5.0.1
        mpfr_version=mpfr-2.4.2
        mpc_version=mpc-0.8.1
        ;;
    Gcc44)
        gcc_version=gcc-4.4.3
        gmp_version=gmp-5.0.1
        mpfr_version=mpfr-2.4.2
        ;;
    Gcc43)
        gcc_version=gcc-4.3.3
        gmp_version=gmp-4.2.4
        mpfr_version=mpfr-2.4.1
        ;;
    *)
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.sh:"
        echo "    Unknown OpenFOAM compiler type '$WM_COMPILER'"
        echo "    Please check your settings"
        echo
        ;;
    esac

    if [ -n "$gcc_version" ]
    then
        gccDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gcc_version
        gmpDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gmp_version
        mpfrDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpfr_version
        mpcDir=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpc_version

        # Check that the compiler directory can be found
        [ -d "$gccDir" ] || {
            echo
            echo "Warning in $WM_PROJECT_DIR/etc/settings.sh:"
            echo "    Cannot find $gccDir installation."
            echo "    Please install this compiler version or if you wish to use the system compiler,"
            echo "    change the 'compilerInstall' setting to 'system' in this file"
            echo
        }

        _foamAddMan     $gccDir/man
        _foamAddPath    $gccDir/bin

        # 64-bit needs lib64, but 32-bit needs lib (not lib32)
        if [ "$WM_ARCH_OPTION" = 64 ]
        then
            _foamAddLib     $gccDir/lib$WM_COMPILER_LIB_ARCH
        else
            _foamAddLib     $gccDir/lib
        fi

        # add in gmp/mpfr libraries
        _foamAddLib     $gmpDir/lib
        _foamAddLib     $mpfrDir/lib

        # add in mpc libraries (not need for older gcc)
        if [ -n "$mpc_version" ]
        then
            _foamAddLib     $mpcDir/lib
        fi

        # used by boost/CGAL:
        export MPFR_ARCH_PATH=$mpfrDir
    fi
    unset gcc_version gccDir
    unset gmp_version gmpDir  mpfr_version mpfrDir  mpc_version mpcDir
    ;;
esac


# boost and CGAL
# ~~~~~~~~~~~~~~

boost_version=boost_1_42_0
cgal_version=CGAL-3.5.1

export BOOST_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$boost_version
export CGAL_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$cgal_version

# enabled if CGAL is available
if [ -d "$CGAL_ARCH_PATH" ]
then
    if [ -d "$BOOST_ARCH_PATH" ]
    then
        _foamAddLib $BOOST_ARCH_PATH/lib
    else
        unset BOOST_ARCH_PATH
    fi
    _foamAddLib $CGAL_ARCH_PATH/lib
else
    unset BOOST_ARCH_PATH CGAL_ARCH_PATH MPFR_ARCH_PATH
fi

unset boost_version cgal_version


# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH MPI_HOME

case "$WM_MPLIB" in
OPENMPI)
    mpi_version=openmpi-1.4.1
    export MPI_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$mpi_version

    # Tell OpenMPI where to find its install directory
    export OPAL_PREFIX=$MPI_ARCH_PATH

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    _foamAddMan     $MPI_ARCH_PATH/man

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

SYSTEMOPENMPI)
    # use the system installed openmpi, get library directory via mpicc
    mpi_version=openmpi-system
    libDir=`mpicc --showme:link | sed -e 's/.*-L\([^ ]*\).*/\1/'`

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using system installed OpenMPI:"
        echo "    compile flags : `mpicc --showme:compile`"
        echo "    link flags    : `mpicc --showme:link`"
        echo "    libmpi dir    : $libDir"
    fi

    _foamAddLib     $libDir

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version libDir
    ;;

MPICH)
    mpi_version=mpich2-1.1.1p1
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$mpi_version
    export MPI_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$mpi_version

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    _foamAddMan     $MPI_ARCH_PATH/share/man

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/$mpi_version
    unset mpi_version
    ;;

MPICH-GM)
    export MPI_ARCH_PATH=/opt/mpi
    export MPICH_PATH=$MPI_ARCH_PATH
    export GM_LIB_PATH=/opt/gm/lib64

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    _foamAddLib     $GM_LIB_PATH

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/mpich-gm
    ;;

HPMPI)
    export MPI_HOME=/opt/hpmpi
    export MPI_ARCH_PATH=$MPI_HOME

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

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib/sparcv9
    _foamAddLib     /opt/FSUNf90/lib/sparcv9
    _foamAddLib     /opt/FJSVpnidt/lib
    ;;

QSMPI)
    export MPI_ARCH_PATH=/usr/lib/mpi
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/qsmpi

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib

    ;;

*)
    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/dummy
    ;;
esac

_foamAddLib $FOAM_MPI_LIBBIN


# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
: ${minBufferSize:=20000000}

if [ "${MPI_BUFFER_SIZE:=$minBufferSize}" -lt $minBufferSize ]
then
    MPI_BUFFER_SIZE=$minBufferSize
fi
export MPI_BUFFER_SIZE


# Enable the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if [ -f $FOAM_LIBBIN/libhoard.so ]
#then
#    export LD_PRELOAD=$FOAM_LIBBIN/libhoard.so:$LD_PRELOAD
#fi


# cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
unset _foamAddPath _foamAddLib _foamAddMan compilerInstall minBufferSize

# ----------------------------------------------------------------- end-of-file
