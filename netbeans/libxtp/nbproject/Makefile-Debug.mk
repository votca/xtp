#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1ce08c7a/ERIs.o \
	${OBJECTDIR}/_ext/1ce08c7a/aobasis.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aodipole_potential.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o \
	${OBJECTDIR}/_ext/2a3bfc3d/aoquadrupole_potential.o \
	${OBJECTDIR}/_ext/1ce08c7a/aoshell.o \
	${OBJECTDIR}/_ext/1ce08c7a/basisset.o \
	${OBJECTDIR}/_ext/1ce08c7a/bfgs-trm.o \
	${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o \
	${OBJECTDIR}/_ext/1ce08c7a/bulkesp.o \
	${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o \
	${OBJECTDIR}/_ext/50c6f79c/jobwriter.o \
	${OBJECTDIR}/_ext/50c6f79c/kmclifetime.o \
	${OBJECTDIR}/_ext/50c6f79c/kmcmultiple.o \
	${OBJECTDIR}/_ext/5e78919f/dftengine.o \
	${OBJECTDIR}/_ext/1ce08c7a/diis.o \
	${OBJECTDIR}/_ext/1ce08c7a/esp2multipole.o \
	${OBJECTDIR}/_ext/1ce08c7a/espfit.o \
	${OBJECTDIR}/_ext/1ce08c7a/exchange_correlation.o \
	${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o \
	${OBJECTDIR}/_ext/1ce08c7a/forces.o \
	${OBJECTDIR}/_ext/1ce08c7a/fourcenter_rep.o \
	${OBJECTDIR}/_ext/1ce08c7a/fourcenters_dft.o \
	${OBJECTDIR}/_ext/1ce08c7a/gdma.o \
	${OBJECTDIR}/_ext/1ce08c7a/geometry_optimization.o \
	${OBJECTDIR}/_ext/1ce08c7a/gnode.o \
	${OBJECTDIR}/_ext/1ce08c7a/grid.o \
	${OBJECTDIR}/_ext/1ce08c7a/gridbox.o \
	${OBJECTDIR}/_ext/6916996f/bse.o \
	${OBJECTDIR}/_ext/6916996f/bse_analysis.o \
	${OBJECTDIR}/_ext/6916996f/gwa.o \
	${OBJECTDIR}/_ext/6916996f/gwbse.o \
	${OBJECTDIR}/_ext/6916996f/gwbseengine.o \
	${OBJECTDIR}/_ext/6916996f/rpa.o \
	${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o \
	${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o \
	${OBJECTDIR}/_ext/f632c409/dma.o \
	${OBJECTDIR}/_ext/f632c409/egwbse.o \
	${OBJECTDIR}/_ext/f632c409/idft.o \
	${OBJECTDIR}/_ext/f632c409/iexcitoncl.o \
	${OBJECTDIR}/_ext/f632c409/igwbse.o \
	${OBJECTDIR}/_ext/1ce08c7a/kmccalculator.o \
	${OBJECTDIR}/_ext/1ce08c7a/lowdin.o \
	${OBJECTDIR}/_ext/1ce08c7a/mixing.o \
	${OBJECTDIR}/_ext/1ce08c7a/mulliken.o \
	${OBJECTDIR}/_ext/1ce08c7a/nbo.o \
	${OBJECTDIR}/_ext/4d261038/numerical_integrations.o \
	${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o \
	${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o \
	${OBJECTDIR}/_ext/1ce08c7a/orbitals.o \
	${OBJECTDIR}/_ext/1ce08c7a/overlap.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o \
	${OBJECTDIR}/_ext/1ce08c7a/qminterface.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmiter.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o \
	${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o \
	${OBJECTDIR}/_ext/60851bbe/cpmd.o \
	${OBJECTDIR}/_ext/60851bbe/gaussian.o \
	${OBJECTDIR}/_ext/60851bbe/nwchem.o \
	${OBJECTDIR}/_ext/60851bbe/orca.o \
	${OBJECTDIR}/_ext/60851bbe/turbomole.o \
	${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o \
	${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o \
	${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o \
	${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o \
	${OBJECTDIR}/_ext/69ca5806/exciton.o \
	${OBJECTDIR}/_ext/1ce08c7a/version.o \
	${OBJECTDIR}/_ext/1ce08c7a/version_nb.o \
	${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ../../src/libxtp/libxtp.a

../../src/libxtp/libxtp.a: ${OBJECTFILES}
	${MKDIR} -p ../../src/libxtp
	${RM} ../../src/libxtp/libxtp.a
	${AR} -rv ../../src/libxtp/libxtp.a ${OBJECTFILES} 
	$(RANLIB) ../../src/libxtp/libxtp.a

${OBJECTDIR}/_ext/1ce08c7a/ERIs.o: ../../src/libxtp/ERIs.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/ERIs.o ../../src/libxtp/ERIs.cc

${OBJECTDIR}/_ext/1ce08c7a/aobasis.o: ../../src/libxtp/aobasis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/aobasis.o ../../src/libxtp/aobasis.cc

${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o: ../../src/libxtp/aomatrices/aodipole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aodipole.o ../../src/libxtp/aomatrices/aodipole.cc

${OBJECTDIR}/_ext/2a3bfc3d/aodipole_potential.o: ../../src/libxtp/aomatrices/aodipole_potential.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aodipole_potential.o ../../src/libxtp/aomatrices/aodipole_potential.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o: ../../src/libxtp/aomatrices/aoecp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoecp.o ../../src/libxtp/aomatrices/aoecp.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o: ../../src/libxtp/aomatrices/aoesp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoesp.o ../../src/libxtp/aomatrices/aoesp.cc

${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o: ../../src/libxtp/aomatrices/aokinetic.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aokinetic.o ../../src/libxtp/aomatrices/aokinetic.cc

${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o: ../../src/libxtp/aomatrices/aomatrix.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aomatrix.o ../../src/libxtp/aomatrices/aomatrix.cc

${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o: ../../src/libxtp/aomatrices/aomomentum.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aomomentum.o ../../src/libxtp/aomatrices/aomomentum.cc

${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o: ../../src/libxtp/aomatrices/aooverlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aooverlap.o ../../src/libxtp/aomatrices/aooverlap.cc

${OBJECTDIR}/_ext/2a3bfc3d/aoquadrupole_potential.o: ../../src/libxtp/aomatrices/aoquadrupole_potential.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/2a3bfc3d
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/2a3bfc3d/aoquadrupole_potential.o ../../src/libxtp/aomatrices/aoquadrupole_potential.cc

${OBJECTDIR}/_ext/1ce08c7a/aoshell.o: ../../src/libxtp/aoshell.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/aoshell.o ../../src/libxtp/aoshell.cc

${OBJECTDIR}/_ext/1ce08c7a/basisset.o: ../../src/libxtp/basisset.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/basisset.o ../../src/libxtp/basisset.cc

${OBJECTDIR}/_ext/1ce08c7a/bfgs-trm.o: ../../src/libxtp/bfgs-trm.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/bfgs-trm.o ../../src/libxtp/bfgs-trm.cc

${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o: ../../src/libxtp/bsecoupling.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/bsecoupling.o ../../src/libxtp/bsecoupling.cc

${OBJECTDIR}/_ext/1ce08c7a/bulkesp.o: ../../src/libxtp/bulkesp.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/bulkesp.o ../../src/libxtp/bulkesp.cc

${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o: ../../src/libxtp/calculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/calculatorfactory.o ../../src/libxtp/calculatorfactory.cc

${OBJECTDIR}/_ext/50c6f79c/jobwriter.o: ../../src/libxtp/calculators/jobwriter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/50c6f79c
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/50c6f79c/jobwriter.o ../../src/libxtp/calculators/jobwriter.cc

${OBJECTDIR}/_ext/50c6f79c/kmclifetime.o: ../../src/libxtp/calculators/kmclifetime.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/50c6f79c
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/50c6f79c/kmclifetime.o ../../src/libxtp/calculators/kmclifetime.cc

${OBJECTDIR}/_ext/50c6f79c/kmcmultiple.o: ../../src/libxtp/calculators/kmcmultiple.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/50c6f79c
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/50c6f79c/kmcmultiple.o ../../src/libxtp/calculators/kmcmultiple.cc

${OBJECTDIR}/_ext/5e78919f/dftengine.o: ../../src/libxtp/dftengine/dftengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/5e78919f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/5e78919f/dftengine.o ../../src/libxtp/dftengine/dftengine.cc

${OBJECTDIR}/_ext/1ce08c7a/diis.o: ../../src/libxtp/diis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/diis.o ../../src/libxtp/diis.cc

${OBJECTDIR}/_ext/1ce08c7a/esp2multipole.o: ../../src/libxtp/esp2multipole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/esp2multipole.o ../../src/libxtp/esp2multipole.cc

${OBJECTDIR}/_ext/1ce08c7a/espfit.o: ../../src/libxtp/espfit.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/espfit.o ../../src/libxtp/espfit.cc

${OBJECTDIR}/_ext/1ce08c7a/exchange_correlation.o: ../../src/libxtp/exchange_correlation.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/exchange_correlation.o ../../src/libxtp/exchange_correlation.cc

${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o: ../../src/libxtp/extractorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/extractorfactory.o ../../src/libxtp/extractorfactory.cc

${OBJECTDIR}/_ext/1ce08c7a/forces.o: ../../src/libxtp/forces.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/forces.o ../../src/libxtp/forces.cc

${OBJECTDIR}/_ext/1ce08c7a/fourcenter_rep.o: ../../src/libxtp/fourcenter_rep.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/fourcenter_rep.o ../../src/libxtp/fourcenter_rep.cc

${OBJECTDIR}/_ext/1ce08c7a/fourcenters_dft.o: ../../src/libxtp/fourcenters_dft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/fourcenters_dft.o ../../src/libxtp/fourcenters_dft.cc

${OBJECTDIR}/_ext/1ce08c7a/gdma.o: ../../src/libxtp/gdma.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/gdma.o ../../src/libxtp/gdma.cc

${OBJECTDIR}/_ext/1ce08c7a/geometry_optimization.o: ../../src/libxtp/geometry_optimization.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/geometry_optimization.o ../../src/libxtp/geometry_optimization.cc

${OBJECTDIR}/_ext/1ce08c7a/gnode.o: ../../src/libxtp/gnode.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/gnode.o ../../src/libxtp/gnode.cc

${OBJECTDIR}/_ext/1ce08c7a/grid.o: ../../src/libxtp/grid.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/grid.o ../../src/libxtp/grid.cc

${OBJECTDIR}/_ext/1ce08c7a/gridbox.o: ../../src/libxtp/gridbox.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/gridbox.o ../../src/libxtp/gridbox.cc

${OBJECTDIR}/_ext/6916996f/bse.o: ../../src/libxtp/gwbse/bse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/bse.o ../../src/libxtp/gwbse/bse.cc

${OBJECTDIR}/_ext/6916996f/bse_analysis.o: ../../src/libxtp/gwbse/bse_analysis.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/bse_analysis.o ../../src/libxtp/gwbse/bse_analysis.cc

${OBJECTDIR}/_ext/6916996f/gwa.o: ../../src/libxtp/gwbse/gwa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwa.o ../../src/libxtp/gwbse/gwa.cc

${OBJECTDIR}/_ext/6916996f/gwbse.o: ../../src/libxtp/gwbse/gwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwbse.o ../../src/libxtp/gwbse/gwbse.cc

${OBJECTDIR}/_ext/6916996f/gwbseengine.o: ../../src/libxtp/gwbse/gwbseengine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/gwbseengine.o ../../src/libxtp/gwbse/gwbseengine.cc

${OBJECTDIR}/_ext/6916996f/rpa.o: ../../src/libxtp/gwbse/rpa.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/6916996f
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/6916996f/rpa.o ../../src/libxtp/gwbse/rpa.cc

${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o: ../../src/libxtp/jobapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/jobapplication.o ../../src/libxtp/jobapplication.cc

${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o: ../../src/libxtp/jobcalculatorfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/jobcalculatorfactory.o ../../src/libxtp/jobcalculatorfactory.cc

${OBJECTDIR}/_ext/f632c409/dma.o: ../../src/libxtp/jobcalculators/dma.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/dma.o ../../src/libxtp/jobcalculators/dma.cc

${OBJECTDIR}/_ext/f632c409/egwbse.o: ../../src/libxtp/jobcalculators/egwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/egwbse.o ../../src/libxtp/jobcalculators/egwbse.cc

${OBJECTDIR}/_ext/f632c409/idft.o: ../../src/libxtp/jobcalculators/idft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/idft.o ../../src/libxtp/jobcalculators/idft.cc

${OBJECTDIR}/_ext/f632c409/iexcitoncl.o: ../../src/libxtp/jobcalculators/iexcitoncl.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/iexcitoncl.o ../../src/libxtp/jobcalculators/iexcitoncl.cc

${OBJECTDIR}/_ext/f632c409/igwbse.o: ../../src/libxtp/jobcalculators/igwbse.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/f632c409
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/f632c409/igwbse.o ../../src/libxtp/jobcalculators/igwbse.cc

${OBJECTDIR}/_ext/1ce08c7a/kmccalculator.o: ../../src/libxtp/kmccalculator.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/kmccalculator.o ../../src/libxtp/kmccalculator.cc

${OBJECTDIR}/_ext/1ce08c7a/lowdin.o: ../../src/libxtp/lowdin.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/lowdin.o ../../src/libxtp/lowdin.cc

${OBJECTDIR}/_ext/1ce08c7a/mixing.o: ../../src/libxtp/mixing.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/mixing.o ../../src/libxtp/mixing.cc

${OBJECTDIR}/_ext/1ce08c7a/mulliken.o: ../../src/libxtp/mulliken.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/mulliken.o ../../src/libxtp/mulliken.cc

${OBJECTDIR}/_ext/1ce08c7a/nbo.o: ../../src/libxtp/nbo.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/nbo.o ../../src/libxtp/nbo.cc

${OBJECTDIR}/_ext/4d261038/numerical_integrations.o: ../../src/libxtp/numerical_integration/numerical_integrations.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/numerical_integrations.o ../../src/libxtp/numerical_integration/numerical_integrations.cc

${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o: ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/radial_euler_maclaurin_rule.o ../../src/libxtp/numerical_integration/radial_euler_maclaurin_rule.cc

${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o: ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/4d261038
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/4d261038/sphere_lebedev_rule.o ../../src/libxtp/numerical_integration/sphere_lebedev_rule.cc

${OBJECTDIR}/_ext/1ce08c7a/orbitals.o: ../../src/libxtp/orbitals.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/orbitals.o ../../src/libxtp/orbitals.cc

${OBJECTDIR}/_ext/1ce08c7a/overlap.o: ../../src/libxtp/overlap.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/overlap.o ../../src/libxtp/overlap.cc

${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o: ../../src/libxtp/qmapemachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmapemachine.o ../../src/libxtp/qmapemachine.cc

${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o: ../../src/libxtp/qmdatabase.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmdatabase.o ../../src/libxtp/qmdatabase.cc

${OBJECTDIR}/_ext/1ce08c7a/qminterface.o: ../../src/libxtp/qminterface.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qminterface.o ../../src/libxtp/qminterface.cc

${OBJECTDIR}/_ext/1ce08c7a/qmiter.o: ../../src/libxtp/qmiter.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmiter.o ../../src/libxtp/qmiter.cc

${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o: ../../src/libxtp/qmmachine.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmmachine.o ../../src/libxtp/qmmachine.cc

${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o: ../../src/libxtp/qmpackagefactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/qmpackagefactory.o ../../src/libxtp/qmpackagefactory.cc

${OBJECTDIR}/_ext/60851bbe/cpmd.o: ../../src/libxtp/qmpackages/cpmd.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/cpmd.o ../../src/libxtp/qmpackages/cpmd.cc

${OBJECTDIR}/_ext/60851bbe/gaussian.o: ../../src/libxtp/qmpackages/gaussian.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/gaussian.o ../../src/libxtp/qmpackages/gaussian.cc

${OBJECTDIR}/_ext/60851bbe/nwchem.o: ../../src/libxtp/qmpackages/nwchem.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/nwchem.o ../../src/libxtp/qmpackages/nwchem.cc

${OBJECTDIR}/_ext/60851bbe/orca.o: ../../src/libxtp/qmpackages/orca.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/orca.o ../../src/libxtp/qmpackages/orca.cc

${OBJECTDIR}/_ext/60851bbe/turbomole.o: ../../src/libxtp/qmpackages/turbomole.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/60851bbe
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/60851bbe/turbomole.o ../../src/libxtp/qmpackages/turbomole.cc

${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o: ../../src/libxtp/sqlapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/sqlapplication.o ../../src/libxtp/sqlapplication.cc

${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o: ../../src/libxtp/statesaversqlite.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/statesaversqlite.o ../../src/libxtp/statesaversqlite.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o: ../../src/libxtp/threecenter_rep.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenter_rep.o ../../src/libxtp/threecenter_rep.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters.o: ../../src/libxtp/threecenters.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters.o ../../src/libxtp/threecenters.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o: ../../src/libxtp/threecenters_dft.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_dft.o ../../src/libxtp/threecenters_dft.cc

${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o: ../../src/libxtp/threecenters_tools.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/threecenters_tools.o ../../src/libxtp/threecenters_tools.cc

${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o: ../../src/libxtp/toolfactory.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/toolfactory.o ../../src/libxtp/toolfactory.cc

${OBJECTDIR}/_ext/69ca5806/exciton.o: ../../src/libxtp/tools/exciton.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/69ca5806
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/69ca5806/exciton.o ../../src/libxtp/tools/exciton.cc

${OBJECTDIR}/_ext/1ce08c7a/version.o: ../../src/libxtp/version.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/version.o ../../src/libxtp/version.cc

${OBJECTDIR}/_ext/1ce08c7a/version_nb.o: ../../src/libxtp/version_nb.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/version_nb.o ../../src/libxtp/version_nb.cc

${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o: ../../src/libxtp/xtpapplication.cc 
	${MKDIR} -p ${OBJECTDIR}/_ext/1ce08c7a
	${RM} "$@.d"
	$(COMPILE.cc) -g -I../../include -I../../../tools/include -I../../../csg/include -I../../../xtp/include -I/usr/include/libxml2 -I/sw/linux/intel/XE13u2/mkl/include -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1ce08c7a/xtpapplication.o ../../src/libxtp/xtpapplication.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ../../src/libxtp/libxtp.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
