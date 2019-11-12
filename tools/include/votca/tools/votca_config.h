/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_TOOLS_VOTCA_CONFIG_H
#define __VOTCA_TOOLS_VOTCA_CONFIG_H

/* OpenMP */
#if defined(_OPENMP)
#include <omp.h>
#endif

/* Linear algebra packages */
/* #undef MKL_FOUND */

/* FFT library */
#define FFTW3_FOUND

/* Version number of package */
#define TOOLS_VERSION "1.6-dev"

#endif  // __VOTCA_TOOLS_VOTCA_CONFIG_H
