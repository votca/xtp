#include<votca/xtp/internalcoords.h>
#include<boost/graph/adjacency_list.hpp>
#include<boost/graph/connected_components.hpp>
#include<iostream>
#include<limits>

namespace votca { namespace xtp {

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> BglGraph;


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

InternalCoords::InternalCoords(const Orbitals& orb):
    InternalCoords(orb, false){};

InternalCoords::InternalCoords(const std::vector<QMAtom*>& _qmm):
    InternalCoords(_qmm, false){};

InternalCoords::InternalCoords(const Orbitals& orb, const bool _withAux):
    InternalCoords(orb.QMAtoms(), _withAux){};

InternalCoords::InternalCoords(const std::vector<QMAtom*>& _qmm, const bool _withAux):
    CoordBase(INTERNAL, _qmm), withAuxiliary(_withAux),
    numBonds(0), numAngles(0), numDihedrals(0){

    // This code implements the algorithm described in
    // https://doi.org/10.1063/1.1515483

    tools::Elements elements;
    bondMatrix = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    BglGraph bondGraph;

    // calculate the fabulous internal coordinates here

    // covalent bonds
    for (int i = 0; i < numAtoms; ++i){
        auto atomI = qmMolecule[i];

        const double iCovRad = elements.getCovRad(atomI->getType(), "bohr");
        const tools::vec iPos = atomI->getPos();

        for (int j = i+1; j < numAtoms; ++j){
            auto atomJ = qmMolecule[j];

            const double jCovRad = elements.getCovRad(atomJ->getType(), "bohr");
            const tools::vec jPos = atomJ->getPos();

            double thresh = 1.3*(iCovRad + jCovRad);

            double dist = abs(iPos - jPos);

            if (dist < thresh){
                bondMatrix(i,j) = dist;
                bondMatrix(j,i) = dist;

                boost::add_edge(i, j, bondGraph);
                numBonds += 1;
            }
        }
    }

    // check the number of connected components

    // This part of the algorithm is a bit involved...

    // First we find the number of connected components
    // These are not necessarily real molecules in our case.
    // Rather, they are atoms that are within covalent range of each other
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


                int i = std::get<0>(closest);
                int j = std::get<1>(closest);
                double dist = std::get<2>(closest);

                bondMatrix(i,j) = dist;
                bondMatrix(i,j) = dist;

                boost::add_edge(i, j, bondGraph);
                numBonds += 1;
            }

        }

        numComponents = boost::connected_components(bondGraph,
                                                    idxInComponent.data());

        if (numComponents != 1){
            throw std::runtime_error("Failed to create a single connected component");
        }

        // If everything was okay so far, we can proceed to calculate the
        // hydrogen bonds
    }


}

int InternalCoords::getPossibleNumMols(){
    return possibleNumMols;
}

} // xtp
} // votca
