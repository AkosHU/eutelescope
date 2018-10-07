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


//Akos
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
_Range(0.1),
_pulseCollectionName(""),
_initialPulseCollectionSize(0),
allCluster(0),
notDouvbleCluster(0),
noiseCollectionVec(NULL)

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
}


void EUTelProcessorALPIDEClusterFilter::init(){
	cout<<"IN INIT"<<endl;
}


bool EUTelProcessorALPIDEClusterFilter::SameCluster(int iEvent, int iCluster, int jEvent,int jCluster)
{
	int nSame=0;
	for(int iPixel=0; iPixel<PixelsOfEvents[iEvent][iCluster].size(); iPixel++)
	{
		for(int jPixel=0; jPixel<PixelsOfEvents[jEvent][jCluster].size(); jPixel++)
		{
			if(PixelsOfEvents[iEvent][iCluster][iPixel]==PixelsOfEvents[jEvent][jCluster][jPixel])
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
			if(PixelsOfEvents[iEvent][iCluster][iPixel]==PixelsOfEvents[jEvent][jCluster][jPixel])
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

void EUTelProcessorALPIDEClusterFilter::readCollections (LCEvent * evt) {

	//In this step we can not handbe, if pulseCollection is alredy exist, so if it is exist, we skip all of the event.
    try
    {
        LCCollectionVec * pulseCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _pulseCollectionName ) );
		streamlog_out ( WARNING5 ) << "pulseCollection is exist. We can not handle this situation, so we skip this event." << endl;
		//throw SkipEventException( this );
    }
    catch ( lcio::DataNotAvailableException& e )
    {
		//Don't do anything.
    }

	//In this step we can not handbe, if sparseClusterCollectionVec is alredy exist, so if it is exist, we skip all of the event.
    try
    {
        LCCollectionVec * sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "oOriginal_zsdata") );
		streamlog_out ( WARNING5 ) << "sparseClusterCollectionVec is exist. We can not handle this situation, so we skip this event." << endl;
		throw SkipEventException( this );
    }
    catch (lcio::DataNotAvailableException& e)
    {
        //Don't do anything.
    }

	//In this step we can not handbe, if zsInputDataCollectionVec is not available, so we skip all of the event, when it is not available.
	try {
	    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( "original_zsdata" ) ) ;
	    streamlog_out ( DEBUG5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << endl;
		cerr<<"zsInputDataCollectionVec AVAILABLE!"<<endl;
  	}
	catch ( lcio::DataNotAvailableException ) {
	    streamlog_out ( WARNING5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << endl;
		cerr<<"zsInputDataCollectionVec NOT AVAILABLE!"<<endl;
		throw SkipEventException( this );
  	}
	notDouvbleCluster++;

    vector<vector<vector<int>>>EventPixels;
    cerr<<"EVENT"<<endl;
	cerr<<"zsData size: "<<zsInputDataCollectionVec->size()<<endl;
	for ( size_t actualCluster=0 ; actualCluster<zsInputDataCollectionVec->size(); actualCluster++) {
		cerr<<"actualCluster: "<<actualCluster<<endl;
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(actualCluster) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );

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
				pixVector.push_back(pix);
			}
			EventPixels.push_back(pixVector);
		}
	}
	PixelsOfEvents.push_back(EventPixels);
	evtVec.push_back(evt);
	if(PixelsOfEvents.size()!=evtVec.size()) {
		streamlog_out ( WARNING6 ) << "PixelsOfEvents.size()!=evtVec.size() !!!! We will enpty both of them!!!!" << endl;
		while(PixelsOfEvents.size()!=0) PixelsOfEvents.erase(PixelsOfEvents.begin());
		while(evtVec.size()!=0) evtVec.erase(evtVec.begin());
	}
}

void EUTelProcessorALPIDEClusterFilter::writeCollection () {
	
}
		

void EUTelProcessorALPIDEClusterFilter::processEvent (LCEvent * evt) {
	cerr<<"IN PROCESSEVENT"<<endl;
	allCluster++;
	readCollections(evt);

	//Search double firing pixels in PixelsOfEvents	
	if(PixelsOfEvents.size()>_nDeep)
	{
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
		writeCollection();
		PixelsOfEvents.erase(PixelsOfEvents.begin());
	}
}

void EUTelProcessorALPIDEClusterFilter::end() 
{
cerr<<"IN END"<<endl;
cerr<<allCluster<<"; "<<notDouvbleCluster<<endl;
}
