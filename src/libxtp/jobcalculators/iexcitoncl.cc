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


#include "iexcitoncl.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/constants.h>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <votca/tools/propertyiomanipulator.h>

#include <votca/ctp/logger.h>
#include <votca/ctp/xinteractor.h>

using namespace boost::filesystem;
using namespace votca::tools;

namespace ub = boost::numeric::ublas;
    
namespace votca { namespace xtp {
    
// +++++++++++++++++++++++++++++ //
// IEXCITON MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IEXCITON::Initialize(tools::Property* options ) {
    

    
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;


    _induce= false;
    _statenumber=1;
     
    string key = "options."+Identify();

        if ( options->exists(key+".job_file")) {
            _jobfile = options->get(key+".job_file").as<string>();
        }
        else {
            throw std::runtime_error("Job-file not set. Abort.");
        }
    key = "options." + Identify();
    if ( options->exists(key+".mapping")) {
        _xml_file = options->get(key+".mapping").as<string>();
        }
    else {
            throw std::runtime_error("Mapping-file not set. Abort.");
        }
    if ( options->exists(key+".emp_file")) {
            _emp_file   = options->get(key+".emp_file").as<string>();
        }
    else {
            throw std::runtime_error("Emp-file not set. Abort.");
        }
    if ( options->exists(key+".statenumber")) {
            _statenumber=options->get(key+".statenumber").as<int>();
           
        }
    else {
        cout << endl << "Statenumber not specified, assume singlet s1 " << flush;
          _statenumber=1;
        }
    if ( options->exists(key+".induce")) {
            _induce   = options->get(key+".induce").as<bool>();
        }     
    
    cout << "done"<< endl;
}



void IEXCITON::PreProcess(ctp::Topology *top) {

    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
     cout << endl << "... ... Initialize MPS-mapper: "  << flush; 
    
    _mps_mapper.GenerateMap(_xml_file, _emp_file, top);
}


void IEXCITON::CustomizeLogger(ctp::QMThread *thread) {
    
    // CONFIGURE LOGGER
    ctp::Logger* log = thread->getLogger();
    log->setReportLevel(ctp::logDEBUG);
    log->setMultithreading(_maverick);

    log->setPreface(ctp::logINFO,    (boost::format("\nT%1$02d INF ...") % thread->getId()).str());
    log->setPreface(ctp::logERROR,   (boost::format("\nT%1$02d ERR ...") % thread->getId()).str());
    log->setPreface(ctp::logWARNING, (boost::format("\nT%1$02d WAR ...") % thread->getId()).str());
    log->setPreface(ctp::logDEBUG,   (boost::format("\nT%1$02d DBG ...") % thread->getId()).str());        
}



ctp::Job::JobResult IEXCITON::EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *opThread) {
    
     // report back to the progress observer
    ctp::Job::JobResult jres = ctp::Job::JobResult();
    
    // get the logger from the thread
    ctp::Logger* pLog = opThread->getLogger();   
    
    // get the information about the job executed by the thread
    int _job_ID = job->getId();
    Property _job_input = job->getInput();  
    list<Property*> segment_list = _job_input.Select( "segment" );    
    int ID_A   = segment_list.front()->getAttribute<int>( "id" );
    string type_A = segment_list.front()->getAttribute<string>( "type" );
    string mps_fileA = segment_list.front()->getAttribute<string>( "mps_file" );
    int ID_B   = segment_list.back()->getAttribute<int>( "id" );
    string type_B = segment_list.back()->getAttribute<string>( "type" );
    string mps_fileB = segment_list.back()->getAttribute<string>( "mps_file" );

    
  
    ctp::Segment *seg_A = top->getSegment( ID_A );   
    assert( seg_A->getName() == type_A );
    
    ctp::Segment *seg_B = top->getSegment( ID_B );
    assert( seg_B->getName() == type_B );
    
    CTP_LOG(ctp::logINFO,*pLog) << ctp::TimeStamp() << " Evaluating pair "  
            << _job_ID << " ["  << ID_A << ":" << ID_B << "]" << flush; 
    
   vector<ctp::APolarSite*> seg_A_raw=ctp::APS_FROM_MPS(mps_fileA,0,opThread);
   vector<ctp::APolarSite*> seg_B_raw=ctp::APS_FROM_MPS(mps_fileB,0,opThread);
   
   ctp::PolarSeg* seg_A_polar=_mps_mapper.MapPolSitesToSeg(seg_A_raw,seg_A);
   ctp::PolarSeg* seg_B_polar=_mps_mapper.MapPolSitesToSeg(seg_B_raw,seg_B);
   

   
   double JAB=EvaluatePair(top,seg_A_polar,seg_B_polar, pLog);
   
   std::vector< ctp::APolarSite* >::iterator it;
    
   for (it = seg_A_raw.begin() ; it !=seg_A_raw.end(); ++it){
         delete *it;
     }
     seg_A_raw.clear();
     
      for (it = seg_B_raw.begin() ; it !=seg_B_raw.end(); ++it){
         delete *it;
     }
     seg_B_raw.clear();
    
   
   delete seg_A_polar;
   delete seg_B_polar;
  
    
    
   
   Property _job_summary;
 

    Property *_job_output = &_job_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    string nameA = seg_A->getName();
    string nameB = seg_B->getName();
    _pair_summary->setAttribute("idA", ID_A);
    _pair_summary->setAttribute("idB", ID_B);
    _pair_summary->setAttribute("typeA", nameA);
    _pair_summary->setAttribute("typeB", nameB);
    Property *_coupling_summary = &_pair_summary->add("Coupling",""); 
    _coupling_summary->setAttribute("jABstatic", JAB);
     



 
   
    jres.setOutput( _job_summary );   
    jres.setStatus(ctp::Job::COMPLETE);
    
    return jres;
}

double IEXCITON::EvaluatePair(ctp::Topology *top,ctp::PolarSeg* Seg1,ctp::PolarSeg* Seg2, ctp::Logger* pLog ){
    
    ctp::XInteractor actor;
    actor.ResetEnergy();
    Seg1->CalcPos();
    Seg2->CalcPos();
    vec s=top->PbShortestConnect(Seg1->getPos(),Seg2->getPos())+Seg1->getPos()-Seg2->getPos();
    //CTP_LOG(logINFO, *pLog) << "Evaluate pair for debugging " << Seg1->getId() << ":" <<Seg2->getId() << " Distance "<< abs(s) << flush; 
    ctp::PolarSeg::iterator pit1;
    ctp::PolarSeg::iterator pit2;
    double E=0.0;
    for (pit1=Seg1->begin();pit1<Seg1->end();++pit1){
        for (pit2=Seg2->begin();pit2<Seg2->end();++pit2){
            actor.BiasIndu(*(*pit1), *(*pit2),s);
            (*pit1)->Depolarize();
            (*pit2)->Depolarize();

  
            E += actor.E_f(*(*pit1), *(*pit2));                           
           
    }
    }

    
  return E*conv::int2eV;  
}


void IEXCITON::WriteJobFile(ctp::Topology *top) {

    cout << endl << "... ... Writing job file " << flush;
    std::ofstream ofs;
    ofs.open(_jobfile.c_str(), std::ofstream::out);
    if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);

 
    ctp::QMNBList::iterator pit;
    ctp::QMNBList &nblist = top->NBList();    

    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    }    

    ofs << "<jobs>" << endl;    
    string tag = "";
    
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        if ((*pit)->getType()==3){
            int id1 = (*pit)->Seg1()->getId();
            string name1 = (*pit)->Seg1()->getName();
            int id2 = (*pit)->Seg2()->getId();
            string name2 = (*pit)->Seg2()->getName();   

            int id = ++jobCount;

            
            string mps_file1=(boost::format("MP_FILES/%s_n2s%d.mps") % name1  % _statenumber).str();
            string mps_file2=(boost::format("MP_FILES/%s_n2s%d.mps") % name1  % _statenumber).str();
            
            Property Input;
            Property *pInput = &Input.add("input","");
            Property *pSegment =  &pInput->add("segment" , boost::lexical_cast<string>(id1) );
            pSegment->setAttribute<string>("type", name1 );
            pSegment->setAttribute<int>("id", id1 );
            pSegment->setAttribute<string>("mps_file",mps_file1);
            pSegment =  &pInput->add("segment" , boost::lexical_cast<string>(id2) );
            pSegment->setAttribute<string>("type", name2 );
            pSegment->setAttribute<int>("id", id2 );
            pSegment->setAttribute<string>("mps_file",mps_file2);

            ctp::Job job(id, tag, Input, ctp::Job::AVAILABLE );
            job.ToStream(ofs,"xml");
        }
    }

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
    cout << endl << "... ... In total " << jobCount << " jobs" << flush;
    
}

void IEXCITON::ReadJobFile(ctp::Topology *top) {

    Property xml;

    vector<Property*> records;
    
    // gets the neighborlist from the topology
    ctp::QMNBList &nblist = top->NBList();
    int _number_of_pairs = nblist.size();
    int _current_pairs=0;
    
    // output using logger
    ctp::Logger _log;
    _log.setReportLevel(ctp::logINFO);
    

    // load the QC results in a vector indexed by the pair ID
    load_property_from_xml(xml, _jobfile);
    list<Property*> jobProps = xml.Select("jobs.job");
    records.resize( _number_of_pairs + 1  );
    
    //to skip pairs which are not in the jobfile
    for (unsigned i=0;i<records.size();i++){
        records[i]=NULL;
    }
    // loop over all jobs = pair records in the job file
    for (list<Property*> ::iterator  it = jobProps.begin(); it != jobProps.end(); ++it) {
        // if job produced an output, then continue with analysis
        if ( (*it)->exists("output") && (*it)->exists("output.pair") ) {
            _current_pairs++;
            // get the output records
            Property& poutput = (*it)->get("output.pair");
            // id's of two segments of a pair
            int idA = poutput.getAttribute<int>("idA");
            int idB = poutput.getAttribute<int>("idB");
            // segments which correspond to these ids           
            ctp::Segment *segA = top->getSegment(idA);
            ctp::Segment *segB = top->getSegment(idB);
            // pair that corresponds to the two segments
            ctp::QMPair *qmp = nblist.FindPair(segA,segB);
            
            if (qmp == NULL) { // there is no pair in the neighbor list with this name
                CTP_LOG_SAVE(ctp::logINFO, _log) << "No pair " <<  idA << ":" << idB << " found in the neighbor list. Ignoring" << flush; 
            }   else {
                //CTP_LOG(logINFO, _log) << "Store in record: " <<  idA << ":" << idB << flush; 
                records[qmp->getId()] = & ((*it)->get("output.pair"));
            }
        } else {
            Property thebadone = (*it)->get("id");
            throw runtime_error("\nERROR: Job file incomplete.\n Job with id "+thebadone.as<string>()+" is not finished. Check your job file for FAIL, AVAILABLE, or ASSIGNED. Exiting\n");
        }
    } // finished loading from the file


    // loop over all pairs in the neighbor list
    CTP_LOG_SAVE(ctp::logINFO, _log) << "Neighborlist size " << top->NBList().size() << flush; 
    for (ctp::QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {
        
        ctp::QMPair *pair = *ipair;
        
        if (records[ pair->getId() ]==NULL) continue; //skip pairs which are not in the jobfile
        
        //Segment* segmentA = pair->Seg1();
        //Segment* segmentB = pair->Seg2();
        
        double Jeff2 = 0.0;
        double jAB=0.0;
        
        //cout << "\nProcessing pair " << segmentA->getId() << ":" << segmentB->getId() << flush;
        
        
        if ( pair->getType() == ctp::QMPair::Excitoncl){
        Property* pair_property = records[ pair->getId() ];
 
        list<Property*> pCoupling = pair_property->Select("Coupling");
 
            for (list<Property*> ::iterator itCoupling = pCoupling.begin(); itCoupling != pCoupling.end(); ++itCoupling) {
                jAB = (*itCoupling)->getAttribute<double>("jABstatic");
            }
        Jeff2 = jAB*jAB;
        
    
        
        pair->setJeff2(Jeff2, 2);
        pair->setIsPathCarrier(true, 2);
        
        
      
        }
    }
                    
    CTP_LOG_SAVE(ctp::logINFO, _log) << "Pairs [total:updated] " <<  _number_of_pairs << ":" << _current_pairs  << flush; 
    cout << _log;
}



}};
