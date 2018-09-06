#include<votca/xtp/internalcoords.h>
#include<boost/graph/connected_components.hpp>
#include<iostream>
#include<limits>
#include<algorithm>
#include<votca/tools/constants.h>
#include<cmath>
#include<random>
namespace votca { namespace xtp {

const std::vector<std::string> electroNegElements {"N", "O", "F", "P", "S", "Cl"};
tools::Elements elements;

template < typename T>
inline bool VectorContains (const T& item, const std::vector<T>& vec){
    return (std::find(vec.begin(), vec.end(), item) != vec.end());
}

inline bool IsElectronegative(const std::string& type){
    return VectorContains(type, electroNegElements);
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
    const double threshFactor = 1.3;
    const double auxThreshFactor = 2.5;

    int numAtoms = qmMolecule.size();

    for (int i = 0; i < numAtoms; ++i){
        auto atomI = qmMolecule[i];

        const double iCovRad = elements.getCovRad(atomI->getType(), "bohr");
        const tools::vec iPos = atomI->getPos();

        for (int j = i+1; j < numAtoms; ++j){
            auto atomJ = qmMolecule[j];

            const double jCovRad = elements.getCovRad(atomJ->getType(), "bohr");
            const tools::vec jPos = atomJ->getPos();

            double thresh = (iCovRad + jCovRad);

            double dist = abs(iPos - jPos);

            if (dist < threshFactor*thresh){
                _bondMatrix(i,j) = dist;
                _bondMatrix(j,i) = dist;

                boost::add_edge(i, j, _bondGraph);
                _numBonds += 1;
                vector.emplace_back(dist);

            } else if (_withAuxiliary && dist < auxThreshFactor*thresh){
                _bondMatrix(i,j) = dist;
                _bondMatrix(j,i) = dist;

                boost::add_edge(i, j, _bondGraph);
                _numAuxBonds += 1;
                _auxBonds.emplace_back(std::make_pair(i,j));
                vector.emplace_back(dist);
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
    std::vector<int> idxInComponent(boost::num_vertices(_bondGraph));

    int numComponents = boost::connected_components(_bondGraph,
                                                    idxInComponent.data());
    _possibleNumMols = numComponents;

    if (numComponents > 1){
        std::vector<std::vector<int>> components(_possibleNumMols);

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

                _bondMatrix(i,j) = dist;
                _bondMatrix(j,i) = dist;

                boost::add_edge(i, j, _bondGraph);
                numInterMolBonds += 1;
                vector.emplace_back(dist);

                if (_withAuxiliary){
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
                            const double dist = abs(qmMolecule[iAtom]->getPos() -
                                                    qmMolecule[jAtom]->getPos());
                            if (dist <= thresholdDist){
                                boost::add_edge(iAtom, jAtom, _bondGraph);
                                numInterMolBonds += 1;
                                _auxBonds.emplace_back(std::make_pair(iAtom, jAtom));
                                vector.emplace_back(dist);
                            }
                        }
                    }
                }
            }
        }
    }

    numComponents = boost::connected_components(_bondGraph,
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
            boost::tie(it, it_end) = boost::adjacent_vertices(HAtomInd, _bondGraph);

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
                                boost::add_edge(HAtomInd, neighBInd, _bondGraph);
                                _bondMatrix(HAtomInd, neighBInd) = dist;
                                _bondMatrix(neighBInd, HAtomInd) = dist;
                                _numHBonds+=1;
                                vector.emplace_back(dist);
                            }
                        }
                    }
                }
            }
        }
    }
}

void InternalCoords::CalculateAnglesDihedrals(){


    // There will be an angle between every triplet of bonded atoms
    // Covalent, interfragment, and hbonds will generate angles
    // auxiliary bonds of all types will not

    double tol = 1e-6;

    auto IsAnAuxBond = [&](const int& i, const int& j) -> bool {
        return (VectorContains(std::make_pair(i, j), _auxBonds) ||
                VectorContains(std::make_pair(j, i), _auxBonds));
    };

    for (int atomAIdx = 0; atomAIdx < numAtoms; ++atomAIdx){
        BglGraph::adjacency_iterator itA, itA_end;
        boost::tie(itA, itA_end) = boost::adjacent_vertices(atomAIdx, _bondGraph);

        const tools::vec atomAPos = qmMolecule[atomAIdx]->getPos();

        for (; itA!=itA_end; ++itA){
            const int atomBIdx = *itA;
            if (IsAnAuxBond(atomAIdx, atomBIdx)) continue;

            const tools::vec atomBPos = qmMolecule[atomBIdx]->getPos();
            const tools::vec BAVec = (atomAPos - atomBPos).normalize();

            BglGraph::adjacency_iterator itB, itB_end;
            boost::tie(itB, itB_end) = boost::adjacent_vertices(atomBIdx, _bondGraph);

            for (; itB != itB_end; ++itB){
                const int atomCIdx = *itB;
                if (IsAnAuxBond(atomBIdx, atomCIdx) ||
                    atomCIdx == atomAIdx) continue;

                const tools::vec atomCPos = qmMolecule[atomCIdx]->getPos();

                auto index = std::make_tuple(atomAIdx, atomBIdx, atomCIdx);
                auto index2 = std::make_tuple(atomCIdx, atomBIdx, atomAIdx);

                auto NotAnAngle = [&] (std::tuple<int,int,int>& index) -> bool {
                    return (_angles.find(index)==_angles.end());
                };

                const tools::vec BCVec = (atomCPos - atomBPos).normalize();

                double cosABC = BAVec*BCVec;

                if (NotAnAngle(index) && NotAnAngle(index2)){
                    // If the angle is > 175, then a special orthogonal angle
                    // should be added, but I don't know how to do that...
                    _angles[index] = std::acos(cosABC);
                    vector.emplace_back(std::acos(cosABC));
                }


                // Now add proper dihedrals
                BglGraph::adjacency_iterator itC, itC_end;
                boost::tie(itC, itC_end) = boost::adjacent_vertices(atomCIdx, _bondGraph);
                tools::vec normPlaneA = BAVec^BCVec;

                for (; itC != itC_end; ++itC){
                    const int atomDIdx = *itC;
                    if (IsAnAuxBond(atomCIdx, atomDIdx) ||
                        atomDIdx == atomBIdx ||
                        atomDIdx == atomAIdx) continue;


                    const tools::vec atomDPos = qmMolecule[atomDIdx]->getPos();
                    auto index = std::make_tuple(atomAIdx, atomBIdx, atomCIdx, atomDIdx);
                    auto index2 = std::make_tuple(atomDIdx, atomCIdx, atomBIdx, atomAIdx);

                    auto NotADihedral = [&](std::tuple<int,int,int,int>& index) -> bool {
                        return (_dihedralss.find(index)==_dihedralss.end());
                    };

                    if (NotADihedral(index) && NotADihedral(index2)){
                        const tools::vec CBVec = -BCVec;
                        const tools::vec CDVec = (atomDPos - atomCPos).normalize();

                        const double cosBCD = BCVec*CDVec;

                        // ABC and BCD must not be 180 degrees
                        if (std::abs(-1 - cosABC) > tol && std::abs(-1 - cosBCD) > tol){

                            tools::vec normPlaneB = (CBVec^CDVec);

                            const double cosPhi = normPlaneA*normPlaneB;

                            _dihedralss[index] = std::acos(cosPhi);
                            vector.emplace_back(std::acos(cosPhi));
                        }
                    }
                }
            }
        }
    }

    _numAngles = _angles.size();
    _numDihedrals = _dihedralss.size();
    if (numAtoms > 3 && _numDihedrals == 0){
        // Need to do something here...
        // arbitrarily choose 4 atoms and check if they form a
        // valid dihedral

        std::random_device rd;
        std::mt19937 g(rd());
        std::uniform_int_distribution<int> dist(0, numAtoms-1);

        auto IsABond = [&] (const int& i, const int& j) -> bool{
            return (_bondMatrix(i,j) > 0);
        };

        std::vector<std::vector<int>> alreadyChosen;
        auto RandomSelector = [&]() -> std::vector<int> {
            std::vector<int> inds;
            do {
                inds = std::vector<int>();
                do {
                    int ind = dist(g);
                    if (!VectorContains(ind, inds)){
                        inds.emplace_back(ind);
                    }
                } while (inds.size() < 4);
            } while (VectorContains(inds, alreadyChosen));

            std::sort(inds.begin(), inds.end());

            do {
                alreadyChosen.emplace_back(inds);
            } while (std::next_permutation(inds.begin(), inds.end()));

            std::reverse(inds.begin(), inds.end());
            return inds;
        };

        std::function<int(int)> factorial =  [&](int num) -> int {
            if (num == 1) return 1;
            return num*factorial(num-1);
        };

        int checkCount = 0;
        int checkCountMax = factorial(numAtoms)/(factorial(4));
        do {
            std::vector<int> inds = RandomSelector();
            do {
                const int atomAIdx = inds[0];
                const int atomBIdx = inds[1];
                const int atomCIdx = inds[2];
                const int atomDIdx = inds[3];

                if (!(IsABond(atomAIdx, atomBIdx) && IsABond(atomBIdx, atomCIdx) &&
                      IsABond(atomCIdx, atomDIdx))) continue;

                const tools::vec BAVec = (qmMolecule[atomAIdx]->getPos() - qmMolecule[atomBIdx]->getPos()).normalize();
                const tools::vec BCVec = (qmMolecule[atomCIdx]->getPos() - qmMolecule[atomBIdx]->getPos()).normalize();

                const tools::vec CBVec = -BCVec;
                const tools::vec CDVec = (qmMolecule[atomDIdx]->getPos() - qmMolecule[atomCIdx]->getPos()).normalize();

                const double cosABC = BAVec*BCVec;
                const double cosBCD = CBVec*CDVec;


                if(abs(-1 - cosABC) > tol && abs(-1 - cosBCD) > tol){
                    const tools::vec normPlaneA = BAVec^BCVec;
                    const tools::vec normPlaneB = CBVec^CDVec;
                    const double cosPhi = normPlaneA*normPlaneB;
                    auto index = std::make_tuple(atomAIdx, atomBIdx, atomCIdx, atomDIdx);
                    auto index2 = std::make_tuple(atomDIdx, atomCIdx, atomBIdx, atomAIdx);

                    auto NotADihedral = [&](std::tuple<int,int,int,int>& index) -> bool {
                        return (_dihedralss.find(index)==_dihedralss.end());
                    };
                    if (NotADihedral(index)){
                        _dihedralss[index] = std::acos(cosPhi);
                        vector.emplace_back(std::acos(cosPhi));
                        _numDihedrals += 1;
                    }
                }

            } while (std::next_permutation(inds.begin(), inds.end()));
            ++checkCount;
        } while (_numDihedrals < 1 && checkCount < checkCountMax);
    }
}

InternalCoords::InternalCoords(const Orbitals& orb):
    InternalCoords(orb, true){};

InternalCoords::InternalCoords(const std::vector<QMAtom*>& _qmm):
    InternalCoords(_qmm, true){};

InternalCoords::InternalCoords(const Orbitals& orb, const bool _withAux):
    InternalCoords(orb.QMAtoms(), _withAux){};

InternalCoords::InternalCoords(const std::vector<QMAtom*>& _qmm, const bool _withAux):
    CoordBase(INTERNAL, _qmm), _withAuxiliary(_withAux),
    _numBonds(0), numInterMolBonds(0), _numHBonds(0),
    _numAngles(0), _numDihedrals(0), _numAuxBonds(0), _bondGraph(numAtoms),
    _bondMatrix(Eigen::MatrixXd::Zero(numAtoms, numAtoms))
{
    // This code implements the algorithm described in
    // https://doi.org/10.1063/1.1515483

    // covalent bonds
    ConnectBonds();

    // Intermolecule bonds
    ConnectMolecules();

    // Hydrogen bonds
    ConnectHBonds();

    CalculateAnglesDihedrals();

}

int InternalCoords::getPossibleNumMols(){
    return _possibleNumMols;
}

int InternalCoords::getNumBonds(){
    return _numBonds;
}

int InternalCoords::getNumHBonds(){
    return _numHBonds;
}

int InternalCoords::getNumAngles(){
    return _numAngles;
}

int InternalCoords::getNumAuxBonds(){
    return _numAuxBonds;
}

int InternalCoords::getNumDihedrals(){
    return _numDihedrals;
}

} // xtp
} // votca
