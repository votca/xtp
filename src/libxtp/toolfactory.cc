/*
 *            Copyright 2009-2017 The VOTCA Development Team
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


#include <votca/xtp/toolfactory.h>

#include "tools/pdb2map.h"
#include "tools/log2mps.h"
#include "tools/pdb2top.h"
#include "tools/exciton.h"
#include "tools/qmanalyze.h"
#include "tools/qmsandbox.h"
#include "tools/spectrum.h"
#include "tools/excitoncoupling.h"
#include "tools/dft.h"
#include "tools/gencube.h"
#include "tools/partialcharges.h"
#include "tools/matrixproduct.h"
#include "tools/densityanalysis.h"
#include "tools/coupling.h"

namespace votca { namespace xtp {

void QMToolFactory::RegisterAll(void)
{

        QMTools().Register<PDB2Map>            ("pdb2map");
        QMTools().Register<Log2Mps>            ("log2mps");
        QMTools().Register<Exciton>            ("exciton");
        QMTools().Register<QMAnalyze>          ("qmanalyze");
        QMTools().Register<QMSandbox>          ("qmsandbox");
        QMTools().Register<Spectrum>           ("spectrum");
        QMTools().Register<ExcitonCoupling>    ("excitoncoupling");
        QMTools().Register<PDB2Top>            ("pdb2top");
        QMTools().Register<DFT>                ("dft");
        QMTools().Register<GenCube>            ("gencube");
        QMTools().Register<Partialcharges>     ("partialcharges");
        //QMTools().Register<MatProd>            ("matrixproduct");
        QMTools().Register<DensityAnalysis>    ("densityanalysis");
        QMTools().Register<Coupling>           ("coupling");

}

}}
