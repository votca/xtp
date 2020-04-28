/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_VOTCA_HDF5_H
#define VOTCA_XTP_VOTCA_HDF5_H

#include <votca/xtp/eigen.h>  //get the MACROS for gcc

#define DO_PRAGMA(X) _Pragma(#X)

#define DISABLE_WARNING(warningName) \
DO_PRAGMA(GCC diagnostic ignored #warningName)

#if (defined STRICT_GNUC) && GCC_VERSION > 90000
#define DISABLE_DEPRECATE_COPY_WARNING(X) \
DO_PRAGMA(GCC diagnostic push) \
DISABLE_WARNING(-Wdeprecated-copy) \
X \
DO_PRAGMA(GCC diagnostic pop)
#else
#define DISABLE_DEPRECATE_COPY_WARNING(X) X
#endif

#include <H5Cpp.h>

#endif  // VOTCA_XTP_VOTCA_HDF5_H
