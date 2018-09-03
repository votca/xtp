#include<votca/xtp/internalcoords.h>
#include<boost/graph/connected_components.hpp>
#include<iostream>
#include<limits>
#include<algorithm>
#include<votca/tools/constants.h>

namespace votca { namespace xtp {

const std::vector<std::string> electroNegElements {"N", "O", "F", "P", "S", "Cl"};
tools::Elements elements;

inline bool IsElectronegative(const std::string& type){
    return (std::find(electroNegElements.begin(),
                      electroNegElements.end(),
                      type) != electroNegElements.end());
}

inline std::tuple<int, int, double> ClosestAtoms(const std::vector<int>& comp1,
                                                 const std::vector<int>& comp2,
                                                 const std::vector<QMAtom*>& qmm){

    double min = std::numeric_limits<double>::infinity();
    std::pair<int, int> closest = std::make_pair(-1, -1);

    for (const auto& i : comp1){
        auto posI = qmm[i]->getPos();
        for (const auto& j: comp2){
            auto posJ = qmm[j]->getPos();

            const double dist = abs(posI - posJ);

            if (dist < min){
                closest = std::make_pair(i, j);
                min = dist;
            }
        }
    }
    return std::make_tuple(closest.first, closest.second, min);
}


void  InternalCoords::ConnectBonds(){
    double threshFactor = 1.3;

    if (withAuxiliary){
        threshFactor = 2.5;
    }

    int numAtoms = qmMolecule.size();
    int numBonds = 0;

    for (int i = 0; i < numAtoms; ++i){
        auto atomI = qmMolecule[i];

        const double iCovRad = elements.getCovRad(atomI->getType(), "bohr");
        const tools::vec iPos = atomI->getPos();

        for (int j = i+1; j < numAtoms; ++j){
            auto atomJ = qmMolecule[j];

            const double jCovRad = elements.getCovRad(atomJ->getType(), "bohr");
            const tools::vec jPos = atomJ->getPos();

            double thresh = threshFactor*(iCovRad + jCovRad);

            double dist = abs(iPos - jPos);

            if (dist < thresh){
                bondMatrix(i,j) = dist;
                bondMatrix(j,i) = dist;

                boost::add_edge(i, j, bondGraph);
                numBonds += 1;
            }
        }
    }
}


void InternalCoords::ConnectMolecules(){
    // This part of the algorithm is a bit involved...
    // First we find the number of connected components
    // These are not necessarily real molecules in our case. But it is a
    // starting point for a topology...
    // Rather, they are atoms that are within covalent range of each other
    // times a factor.
    // Since this code is for calculating internal coordinates, we must only
    // have a single connected component in our graph.
    // This probably has some name in Graph theory, but I don't know it.

    // this will contain the component index that an atom belongs to.
    std::vector<int> idxInComponent(boost::num_vertices(bondGraph));

    int numComponents = boost::connected_components(bondGraph,
                                                    idxInComponent.data());
    possibleNumMols = numComponents;

    if (numComponents > 1){
        std::vector<std::vector<int>> components(possibleNumMols);

        for (int atomIdx = 0; atomIdx < numAtoms; ++atomIdx){
            int componentIdx = idxInComponent[atomIdx];
            components[componentIdx].push_back(atomIdx);
        }


        // Now connect the closest two atoms in each component with a bond
        for (int compI = 0; compI < numComponents; ++compI){
            std::vector<int> compIAtoms = components[compI];

            for (int compJ = compI + 1; compJ < numComponents; ++compJ){
                std::vector<int> compJAtoms = components[compJ];

                auto closest = ClosestAtoms(compIAtoms, compJAtoms,
                                            qmMolecule);


                const int i = std::get<0>(closest);
                const int j = std::get<1>(closest);
                const double dist = std::get<2>(closest);

                bondMatrix(i,j) = dist;
                bondMatrix(j,i) = dist;

                boost::add_edge(i, j, bondGraph);
                numInterMolBonds += 1;

                if (withAuxiliary){
                    // auxiliary interfragment bonds if needed

                    // the paper says the threshold for this can be
                    // 2.0 ang or less than 1.3 times the minimum
                    // interfragment distance. I'll go with the 1.3
                    // factor...

                    double thresholdDist = 1.3*dist;

                    // Now connect atoms within that threshold belonging
                    // to different components
                    for (const auto& iAtom : compIAtoms){
                        for (const auto& jAtom : compJAtoms){
                            double dist = abs(qmMolecule[iAtom]->getPos() -
                                              qmMolecule[jAtom]->getPos());
                            if (dist <= thresholdDist){
                                boost::add_edge(iAtom, jAtom, bondGraph);
                                numInterMolBonds += 1;
                            }
                        }
                    }
                }
            }
        }
    }

    numComponents = boost::connected_components(bondGraph,
                                                idxInComponent.data());

    if (numComponents != 1){
        throw std::runtime_error("Failed to create a single connected component");
    }
}


void InternalCoords::ConnectHBonds(){
    // we can proceed to calculate the hydrogen bonds
    // Assume we have this situation:
    //                H-----B
    //               /
    //              A
    // Where A,B are some elements, A is bonded to H, and we are checking
    // if B and H can form a hydrogen bond.
    // Then:
    // 1. A, B must be electronegative
    // 2. The distance between H, B must be greater than 0.9 times the sum of their
    // (H and B) covalent radii and less than the sum of their VdW radii
    // 3. The angle AHB must be > 90 deg

    double hVdWRad = tools::conv::ang2bohr*elements.getVdWChelpG("H");
    for (int i = 0; i<numAtoms; ++i){
        // first, for each H atom...
        if (qmMolecule[i]->getType() == "H"){
            const int HAtomInd = i;
            BglGraph::adjacency_iterator it, it_end;
            boost::tie(it, it_end) = boost::adjacent_vertices(HAtomInd, bondGraph);

            const tools::vec HAtomPos = qmMolecule[i]->getPos();

            for (;it != it_end; ++it){
                // if a BONDED neighbour is electronegative
                if ( IsElectronegative( qmMolecule[*it]->getType() ) ){
                    const int neighAInd = *it;
                    const tools::vec neighAPos = qmMolecule[neighAInd]->getPos();

                    // for each neighbour within range
                    for (int neighBInd = 0; neighBInd < numAtoms; ++neighBInd){

                        if (neighBInd == HAtomInd || neighBInd == neighAInd ||
                            // if the neighbour is electronegative
                            !IsElectronegative(qmMolecule[neighBInd]->getType()))
                            continue;

                        const tools::vec neighBPos = qmMolecule[neighBInd]->getPos();

                        const double dist = abs(HAtomPos - neighBPos);

                        double lowerBound = 0.9*(elements.getCovRad(qmMolecule[neighBInd]->getType(), "bohr") + elements.getCovRad("H", "bohr"));


                        double upperBound = tools::conv::ang2bohr*(elements.getVdWChelpG(qmMolecule[neighBInd]->getType()) + elements.getVdWChelpG("H"));
                        // And it is within range
                        if (lowerBound < dist && dist < upperBound){

                            const tools::vec HABond = qmMolecule[neighAInd]->getPos() - HAtomPos;
                            const tools::vec HBBond = qmMolecule[neighBInd]->getPos() - HAtomPos;
                            const double cosTheta = HABond*HBBond;
                            // And the angle is greater than 90
                            // ie cos(angle) < 0
                            if (cosTheta < 0){
                                boost::add_edge(HAtomInd, neighBInd, bondGraph);
                                bondMatrix(HAtomInd, neighBInd) = dist;
                                bondMatrix(neighBInd, HAtomInd) = dist;
                                numHBonds+=1;
                            }
                       }
                    }
                }
            }
        }
    }
}

InternalCoords::InternalCoords(const Orbitals& orb):
    InternalCoords(orb, false){};

InternalCoords::InternalCoords(const std::vector<QMAtom*>& _qmm):
    InternalCoords(_qmm, false){};

InternalCoords::InternalCoords(const Orbitals& orb, const bool _withAux):
    InternalCoords(orb.QMAtoms(), _withAux){};

InternalCoords::InternalCoords(const std::vector<QMAtom*>& _qmm, const bool _withAux):
    CoordBase(INTERNAL, _qmm), withAuxiliary(_withAux),
    numBonds(0), numInterMolBonds(0), numHBonds(0),
    numAngles(0), numDihedrals(0){

    // This code implements the algorithm described in
    // https://doi.org/10.1063/1.1515483

    bondMatrix = Eigen::MatrixXd::Zero(numAtoms, numAtoms);

    double threshFactor = 1.1;

    if (withAuxiliary){
        threshFactor = 2.5;
    }

    // covalent bonds
    ConnectAtomsWithin(threshFactor);

    // Intermolecule bonds
    ConnectMolecules();

    // Hydrogen bonds
    ConnectHBonds();

}


int InternalCoords::getPossibleNumMols(){
    return possibleNumMols;
}

int InternalCoords::getNumBonds(){
    return numBonds;
}

int InternalCoords::getNumHBonds(){
    return numHBonds;
}


} // xtp
} // votca
