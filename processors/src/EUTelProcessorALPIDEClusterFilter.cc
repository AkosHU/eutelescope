// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// ROOT includes:
#include "TVector3.h"

// eutelescope includes ".h"

#include "EUTelGeometryTelescopeGeoDescription.h"

#include "EUTelProcessorALPIDEClusterFilter.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"

#include "EUTelSimpleVirtualCluster.h"
#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricClusterImpl.h"
#include "EUTelFFClusterImpl.h"
#include "EUTelDFFClusterImpl.h"
#include "EUTelBrickedClusterImpl.h"
#include "EUTelSparseClusterImpl.h"

#include "EUTelExceptions.h"
#include "EUTelAlignmentConstant.h"
#include "EUTelReferenceHit.h"


#include "EUTelVirtualCluster.h"
#include "EUTelMatrixDecoder.h"
#include "EUTelTrackerDataInterfacerImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Exceptions.h"


// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerPulseImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTime.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <list>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;
using namespace marlin;
using namespace eutelescope;
using namespace lcio;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#endif

EUTelProcessorALPIDEClusterFilter::EUTelProcessorALPIDEClusterFilter () : Processor("EUTelProcessorALPIDEClusterFilter"),
_nzsDataCollectionName(""),
_zsDataCollectionName(""),
_noiseCollectionName(""),
_nDeep(5),
_nShift(10),
_Range(0.1),
_pulseCollectionName(""),
_initialPulseCollectionSize(0),
allCluster(0),
notDouvbleCluster(0),
noiseCollectionVec(NULL),
_sparseClusterCollectionName(""),
_noOfDetector(0),
_isGeometryReady(),
_totClusterMap(),
_NoEvent()

{

    registerInputCollection (LCIO::TRACKERDATA, "NZSDataCollectionName",
                             "Input calibrated data not zero suppressed collection name",
                             _nzsDataCollectionName, string ("data"));

    registerInputCollection (LCIO::TRACKERDATA, "ZSDataCollectionName",
                             "Input of Zero Suppressed data",
                             _zsDataCollectionName, string ("zsdata") );

    registerInputCollection (LCIO::TRACKERDATA, "NoiseCollectionName",
                             "Noise (input) collection name",
                             _noiseCollectionName, string("noise"));

	registerOutputCollection(LCIO::TRACKERPULSE, "PulseCollectionName",
                             "Cluster (output) collection name",
                             _pulseCollectionName, string("cluster"));

	registerOutputCollection(LCIO::TRACKERPULSE, "sparseClusterCollectionName",
                             "Cluster (output) collection name to _sparseClusterCollectionName",
                             _sparseClusterCollectionName, string("filtered_zsdata"));

}


void EUTelProcessorALPIDEClusterFilter::init(){
	cout<<"IN INIT"<<endl;
	_isGeometryReady=false;
	_NoEvent=0;
}


bool EUTelProcessorALPIDEClusterFilter::SameCluster(int iEvent, int iCluster, int jEvent,int jCluster)
{
	int nSame=0;
	for(int iPixel=0; iPixel<PixelsOfEvents[iEvent][iCluster].size(); iPixel++)
	{
		for(int jPixel=0; jPixel<PixelsOfEvents[jEvent][jCluster].size(); jPixel++)
		{
			if(PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][0] && 
			PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][1] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][2] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][3])
			{
				nSame++;
				break;
			}
		}
	}
	if(nSame>PixelsOfEvents[iEvent][iCluster].size() * _Range || nSame>PixelsOfEvents[jEvent][jCluster].size() * _Range) return true;
	return false;
}

void EUTelProcessorALPIDEClusterFilter::AddCluster(int iEvent, int iCluster, int jEvent,int jCluster)
{
	for(int jPixel=0; jPixel<PixelsOfEvents[jEvent][jCluster].size(); jPixel++)
	{
		bool samePixel=false;
		for(int iPixel=0; iPixel<PixelsOfEvents[iEvent][iCluster].size(); iPixel++)
		{
			if(PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][0] && 
			PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][1] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][2] &&
			PixelsOfEvents[iEvent][iCluster][iPixel][0]==PixelsOfEvents[jEvent][jCluster][jPixel][3])
			{
				samePixel=true;
				break;
			}
		}
		if(!samePixel) PixelsOfEvents[iEvent][iCluster].push_back(PixelsOfEvents[jEvent][jCluster][jPixel]);
	}
}

void EUTelProcessorALPIDEClusterFilter::DeletCluster(int jEvent, int jCluster){
	PixelsOfEvents[jEvent].erase(PixelsOfEvents[jEvent].begin()+jCluster);
}

void EUTelProcessorALPIDEClusterFilter::readCollections (LCCollectionVec * zsInputDataCollectionVec) {
    vector<vector<vector<int>>>EventPixels;
    //cerr<<"EVENT"<<endl;
	//cerr<<"zsData size: "<<zsInputDataCollectionVec->size()<<endl;
	for ( size_t actualCluster=0 ; actualCluster<zsInputDataCollectionVec->size(); actualCluster++) {
		//cerr<<"actualCluster: "<<actualCluster<<endl;
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(actualCluster) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
		//cout<<"Sensor id: "<<zsData->getTime()<<endl;

		int clusterSize = zsData->getChargeValues().size()/4;
		vector<int> X(clusterSize);
		vector<int> Y(clusterSize);
		Cluster cluster;
		if ( type == kEUTelGenericSparsePixel )
		{
			vector<vector<int> > pixVector;
			auto sparseData = EUTelTrackerDataInterfacerImpl<EUTelGenericSparsePixel>(zsData);	
			for(size_t iPixel = 0; iPixel < sparseData.size(); iPixel++ )
			{
				auto& pixel = sparseData.at( iPixel );
				X[iPixel] = pixel.getXCoord();
				Y[iPixel] = pixel.getYCoord();
				vector<int> pix;
				pix.push_back(X[iPixel]);
				pix.push_back(Y[iPixel]);
				//pix[2] will be the sensor id.
				pix.push_back((int)cellDecoder(zsData)["sensorID"]);
				//pix[3] will be the pixel type.
				pix.push_back((int)cellDecoder( zsData )["sparsePixelType"]);
				//pix[4] will be the time (it is like an id) of the cluster
				pix.push_back(zsData->getTime());
				//pix[5] will bw Signal
				pix.push_back(pixel.getSignal());
				//pix[6] will be the time (from the pixel)
				pix.push_back(pixel.getTime());
				pixVector.push_back(pix);
			}
			EventPixels.push_back(pixVector);
		}
	}
	if(EventPixels.size()!=0) {
 	   	PixelsOfEvents.push_back(EventPixels);
    }

}

void EUTelProcessorALPIDEClusterFilter::writeCollection (LCCollectionVec * sparseClusterCollectionVec, LCCollectionVec * pulseCollection) {
	CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );
	CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
	if(PixelsOfEvents.size()>_nDeep) {
		for(int iCluster=0; iCluster<PixelsOfEvents[0].size();iCluster++) {
			// prepare a TrackerData to store the cluster candidate
    	    auto zsCluster = std::make_unique<TrackerDataImpl>();
        	// prepare a reimplementation of sparsified cluster
        	auto sparseCluster = std::make_unique<EUTelSparseClusterImpl<EUTelGenericSparsePixel>>(zsCluster.get());
			int sensorID;
			int TYPE;
			int TIME;

			while(PixelsOfEvents[0][iCluster].size()>0)
			{
				EUTelGenericSparsePixel Pixel;
				Pixel.setXCoord(PixelsOfEvents[0][iCluster][0][0]);
				Pixel.setYCoord(PixelsOfEvents[0][iCluster][0][1]);
				Pixel.setTime(PixelsOfEvents[0][iCluster][0][6]);
				Pixel.setSignal(PixelsOfEvents[0][iCluster][0][5]);
				sensorID=PixelsOfEvents[0][iCluster][0][2];
				TYPE=PixelsOfEvents[0][iCluster][0][3];
				TIME=PixelsOfEvents[0][iCluster][0][4];
				PixelsOfEvents[0][iCluster].erase(PixelsOfEvents[0][iCluster].begin());
				sparseCluster->push_back( Pixel );
			}
			//cerr<<"CHECK POINT 4"<<endl;
			if ( sparseCluster->size() > 0)
			{
				//cerr<<"CHECK POINT  sparseCluster->size() > 0"<<endl;
				// set the ID for this zsCluster
    	        idZSClusterEncoder["sensorID"] = static_cast<int >(sensorID);
    	        //idZSClusterEncoder["sparsePixelType"] = static_cast<int> (type);		//cout<<"TIPE: "<<type<<endl;
				idZSClusterEncoder["sparsePixelType"]= static_cast<int >(TYPE);
    	        idZSClusterEncoder["quality"] = 0;
    	        idZSClusterEncoder.setCellID( zsCluster.get() );
    	        zsCluster->setTime(TIME);

				// add it to the cluster collection
    	        sparseClusterCollectionVec->push_back( zsCluster.get() );

				// prepare a pulse for this cluster
    	        auto zsPulse = std::make_unique<TrackerPulseImpl>();
    	        idZSPulseEncoder["sensorID"] = static_cast<int >(sensorID);
    	        idZSPulseEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
    	        idZSPulseEncoder.setCellID( zsPulse.get() );
    	        zsPulse->setTime(TIME);
    	        //zsPulse->setCharge( sparseCluster->getTotalCharge() );
    	        zsPulse->setTrackerData( zsCluster.release() );
    	        pulseCollection->push_back( zsPulse.release() );
				_totClusterMap [static_cast<int >(sensorID)] +=1;
			}
		}
		PixelsOfEvents.erase(PixelsOfEvents.begin());
	}
}

void EUTelProcessorALPIDEClusterFilter::filter () {
	if(PixelsOfEvents.size()>_nDeep)
	{
		//cerr<<"	IN NDEEP"<<endl;
		for(int iCluster=0; iCluster<PixelsOfEvents[0].size();iCluster++)
		{
			//notDouvbleCluster++;
			for(int jEvent=1; jEvent<_nDeep; jEvent++)
			{
				bool wasSameCluster=false;
				for(int jCluster=0; jCluster<PixelsOfEvents[jEvent].size(); jCluster++)
				{
					if(SameCluster(0,iCluster,jEvent,jCluster))
					{
						AddCluster(0,iCluster,jEvent,jCluster);
						DeletCluster(jEvent,jCluster);
						wasSameCluster=true;
					}
				}
				if(!wasSameCluster) break;
			}
		}
	}
}

void EUTelProcessorALPIDEClusterFilter::processEvent (LCEvent * evt) {
	EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt);
	if ( event->getEventType() == kEORE )
    {
        streamlog_out ( DEBUG4 ) <<  "EORE found: nothing else to do." <<  endl;
        return;
    }
    else if ( event->getEventType() == kUNKNOWN )
    {
        streamlog_out ( WARNING2 ) << "Event number " << event->getEventNumber()
                                   << " is of unknown type. Continue considering it as a normal Data Event." << endl;
    }
	//cout<<"EVT NUMB: "<<evt->getEventNumber()<<endl;
	_noOfDetector =0;
  	_clusterAvailable = true;
  	try {
	    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( "original_zsdata" ) ) ;
	    streamlog_out ( DEBUG5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << endl;
		//cerr<<"zsInputDataCollectionVec AVAILABLE!"<<endl;
		_noOfDetector += zsInputDataCollectionVec->getNumberOfElements();
		allCluster++;
		CellIDDecoder<TrackerDataImpl > cellDecoder( zsInputDataCollectionVec );
		for ( size_t i = 0; i < zsInputDataCollectionVec->size(); ++i )
        {
            TrackerDataImpl * data = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt( i ) ) ;
            _totClusterMap.insert( make_pair( cellDecoder( data )[ "sensorID" ] , 0 ));
        }
  	} catch ( lcio::DataNotAvailableException ) {
	    //streamlog_out ( DEBUG5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << endl;
	    _clusterAvailable = false;
		//throw SkipEventException( this );
		//cerr<<"zsInputDataCollectionVec NOT AVAILABLE!"<<endl;
  	}

	if ( _noOfDetector == 0 && _isGeometryReady==false) {
        //streamlog_out( WARNING2 ) << "Unable to initialize the geometry. Trying with the following event" << endl;
		_isGeometryReady=false;
        //throw SkipEventException( this );
    } else {
		_isGeometryReady=true;
    }
	bool isDummyAlreadyExisting = false;
  	LCCollectionVec * sparseClusterCollectionVec = NULL;
  	ID = 0;
  	int TYPE=0;
    try
    {
        sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( _sparseClusterCollectionName ) );
        isDummyAlreadyExisting = true ;
    }
    catch (lcio::DataNotAvailableException& e)
    {
        sparseClusterCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
        isDummyAlreadyExisting = false;
    }
	LCCollectionVec * pulseCollection;
    bool pulseCollectionExists = false;
    _initialPulseCollectionSize = 0;
    try
    {
        pulseCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _pulseCollectionName ) );
        pulseCollectionExists = true;
        _initialPulseCollectionSize = pulseCollection->size();
		//cerr<<"pulseCollection AVAILABLE!"<<endl;
    }
    catch ( lcio::DataNotAvailableException& e )
    {
        pulseCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
		//cerr<<"pulseCollection DONE!"<<endl;
    }
	    // prepare an encoder also for the pulse collection
	noiseCollectionVec = 0;
    try
    {
        noiseCollectionVec  = dynamic_cast < LCCollectionVec * > (evt->getCollection( _noiseCollectionName ));
        streamlog_out ( DEBUG4 ) << "noiseCollectionName: " << _noiseCollectionName.c_str() << " found " << endl;
    }
    catch (lcio::DataNotAvailableException& e )
    {
        streamlog_out ( DEBUG4 ) << "No noise pixel DB collection found in the event" << endl;
    }
	if(_clusterAvailable) {
		readCollections(zsInputDataCollectionVec);
		filter();
		writeCollection(sparseClusterCollectionVec, pulseCollection);
	}
	if ( ! isDummyAlreadyExisting )
    {
		//cerr<<"CHECK POINT ! isDummyAlreadyExisting"<<endl;
        if ( sparseClusterCollectionVec->size() != 0 )
        {
			//notDouvbleCluster+=1;
			//cerr<<"CHECK POINT sparseClusterCollectionVec"<<endl;
            evt->addCollection( sparseClusterCollectionVec, _sparseClusterCollectionName );
        }
        else
        {
            delete sparseClusterCollectionVec;
        }
    }
	//cerr<<"CHECK POINT 7"<<endl;
	// if the pulseCollection is not empty add it to the event
    if ( ! pulseCollectionExists && ( pulseCollection->size() != _initialPulseCollectionSize ))
    {
		notDouvbleCluster+=1;
        evt->addCollection( pulseCollection, _pulseCollectionName );
		//cout<<_pulseCollectionName<<endl;
    }

    if ( ! pulseCollectionExists && ( pulseCollection->size() == _initialPulseCollectionSize ) )
    {
        delete pulseCollection;
    }
	_NoEvent++;
}

void EUTelProcessorALPIDEClusterFilter::end() 
{
	cerr<<"IN END"<<endl;
	cerr<<allCluster<<"; "<<notDouvbleCluster<<endl;
	map< int, int >::iterator iter = _totClusterMap.begin();
	while ( iter != _totClusterMap.end() ) {
    	streamlog_out ( MESSAGE2 ) << "Found " << iter->second << " clusters on detector " << iter->first << endl;
    	++iter;
    }
}
