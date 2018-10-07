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

void EUTelProcessorALPIDEClusterFilter::DeletCluster(int jEvent, int jCluster)
{
	PixelsOfEvents[jEvent].erase(PixelsOfEvents[jEvent].begin()+jCluster);
}

void EUTelProcessorALPIDEClusterFilter::processEvent (LCEvent * evt) {
  cerr<<"IN PROCESSEVENT"<<endl;
  //only for noise
  //CellIDDecoder<TrackerDataImpl> noiseDecoder( noiseCollectionVec );
  bool isDummyAlreadyExisting = false;
  LCCollectionVec * sparseClusterCollectionVec = NULL;
  ID = 0;
  int TYPE=0;
    try
    {
        sparseClusterCollectionVec = dynamic_cast< LCCollectionVec* > ( evt->getCollection( "oOriginal_zsdata") );
        isDummyAlreadyExisting = true ;
		cerr<<"sparseClusterCollection Vec AVAILABLE!"<<endl;
    }
    catch (lcio::DataNotAvailableException& e)
    {
        sparseClusterCollectionVec = new LCCollectionVec(LCIO::TRACKERDATA);
        isDummyAlreadyExisting = false;
		cerr<<"sparseClusterCollection Vec DONE!"<<endl;
    }
	CellIDEncoder<TrackerDataImpl> idZSClusterEncoder( EUTELESCOPE::ZSCLUSTERDEFAULTENCODING, sparseClusterCollectionVec );
	LCCollectionVec * pulseCollection;
    bool pulseCollectionExists = false;
    _initialPulseCollectionSize = 0;
    try
    {
        pulseCollection = dynamic_cast< LCCollectionVec * > ( evt->getCollection( _pulseCollectionName ) );
        pulseCollectionExists = true;
        _initialPulseCollectionSize = pulseCollection->size();
		cerr<<"pulseCollection AVAILABLE!"<<endl;
    }
    catch ( lcio::DataNotAvailableException& e )
    {
        pulseCollection = new LCCollectionVec(LCIO::TRACKERPULSE);
		cerr<<"pulseCollection DONE!"<<endl;
    }
	    // prepare an encoder also for the pulse collection
    CellIDEncoder<TrackerPulseImpl> idZSPulseEncoder(EUTELESCOPE::PULSEDEFAULTENCODING, pulseCollection);
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

  _clusterAvailable = true;
  try {
	    zsInputDataCollectionVec = dynamic_cast< LCCollectionVec * > ( evt->getCollection( "original_zsdata" ) ) ;
	    streamlog_out ( DEBUG5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " found " << endl;
		cerr<<"zsInputDataCollectionVec AVAILABLE!"<<endl;
  } catch ( lcio::DataNotAvailableException ) {
	    streamlog_out ( DEBUG5 ) << "zsInputDataCollectionVec: " << _zsDataCollectionName.c_str() << " not found " << endl;
	    _clusterAvailable = false;
		cerr<<"zsInputDataCollectionVec NOT AVAILABLE!"<<endl;
  }
  //cerr<<"Cluster available? : "<<_clusterAvailable<<endl;
  if (_clusterAvailable)
  {
    vector<vector<vector<int>>>EventPixels;
    cerr<<"EVENT"<<endl;
	cerr<<"zsData size: "<<zsInputDataCollectionVec->size()<<endl;
	int sensorID;
	for ( size_t actualCluster=0 ; actualCluster<zsInputDataCollectionVec->size(); actualCluster++)
	{
		//cerr<<"actualCluster: "<<actualCluster<<endl;
		CellIDDecoder<TrackerDataImpl> cellDecoder( zsInputDataCollectionVec );
		TrackerDataImpl * zsData = dynamic_cast< TrackerDataImpl * > ( zsInputDataCollectionVec->getElementAt(actualCluster) );
		SparsePixelType   type   = static_cast<SparsePixelType> ( static_cast<int> (cellDecoder( zsData )["sparsePixelType"]) );
		TYPE=static_cast<int> (cellDecoder( zsData )["sparsePixelType"]);
		sensorID = static_cast<int > ( cellDecoder( zsData )["sensorID"] );

		int clusterSize = zsData->getChargeValues().size()/4;
		vector<int> X(clusterSize);
		vector<int> Y(clusterSize);
		Cluster cluster;
		if ( type == kEUTelGenericSparsePixel )
		{
			//starting actual cluster analysis
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
			//cerr<<"dutID: "<<(int)cellDecoder(zsData)["sensorID"]<<endl;
			//cerr<<"sensorID: "<<pixVector[0][2]<<endl;
			//cerr<<"clustersize: "<<pixVector.size()<<endl;
			EventPixels.push_back(pixVector);
			allCluster++;
		}

	}
	//cerr<<"CHECK POINT 1"<<endl;
	//Add the new event to PixelsOfEvents: PixelsOfEvent[event][cluster][pixel][x/y/layer]
	PixelsOfEvents.push_back(EventPixels);

	//Search double firing pixels in PixelsOfEvents
	if(PixelsOfEvents.size()>_nDeep)
	{
		//cerr<<"CHECK POINT 2"<<endl;
		for(int iCluster=0; iCluster<PixelsOfEvents[0].size();iCluster++)
		{
			notDouvbleCluster++;
			// prepare a TrackerData to store the cluster candidate
        	auto zsCluster = std::make_unique<TrackerDataImpl>();
        	// prepare a reimplementation of sparsified cluster
        	auto sparseCluster = std::make_unique<EUTelSparseClusterImpl<EUTelGenericSparsePixel>>(zsCluster.get());
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
			//cerr<<"CHECK POINT 3"<<endl;
			cerr<<PixelsOfEvents.size()<<"; "<<PixelsOfEvents[0].size()<<"; "<<PixelsOfEvents[0][iCluster].size()<<"; "<<PixelsOfEvents[0][iCluster][0].size()<<endl;
			while(PixelsOfEvents[0][iCluster].size()>0)
			{
				cerr<<"CHECK POINT PixelsOfEvents[0][iCluster].size()>0"<<endl;
				// get the noise matrix with the right detectorID
        		//TrackerDataImpl* noise  = dynamic_cast<TrackerDataImpl*>   (noiseCollectionVec->getElementAt( _ancillaryIndexMap[ sensorID ] ));
        		// prepare the matrix decoder
        		//EUTelMatrixDecoder matrixDecoder( noiseDecoder , noise );
        		// prepare a vector to store the noise values
        		//vector<float> noiseValueVec;


				/*while(!cluCandidate.empty())
        		{
        	    EUTelGenericSparsePixel pixel = cluCandidate.front();
            	cluCandidate.erase( cluCandidate.begin() );

            	int index = matrixDecoder.getIndexFromXY( pixel.getXCoord(), pixel.getYCoord() );
            	if( _hitIndexMapVec[idetector].find( index ) != _hitIndexMapVec[idetector].end() )
            	{
            	    // do nothing
            	}
            	else
            	{
            	    sparseCluster->push_back( pixel );
            	    //noiseValueVec.push_back(noise->getChargeValues()[ index ]);
            	}
        		}*/

				EUTelGenericSparsePixel Pixel;
				Pixel.setXCoord(PixelsOfEvents[0][iCluster][0][0]);
				Pixel.setYCoord(PixelsOfEvents[0][iCluster][0][1]);
				PixelsOfEvents[0][iCluster].erase(PixelsOfEvents[0][iCluster].begin());
				sparseCluster->push_back( Pixel );
			}
			cerr<<"CHECK POINT 4"<<endl;
			if ( sparseCluster->size() > 0)
			{
				cerr<<"CHECK POINT  sparseCluster->size() > 0"<<endl;
				// set the ID for this zsCluster
                idZSClusterEncoder["sensorID"] = sensorID;
                //idZSClusterEncoder["sparsePixelType"] = static_cast<int> (type);
				//cout<<"TIPE: "<<type<<endl;
				idZSClusterEncoder["sparsePixelType"]=TYPE;
                idZSClusterEncoder["quality"] = 0;
                idZSClusterEncoder.setCellID( zsCluster.get() );
                zsCluster->setTime(ID);

				// add it to the cluster collection
                sparseClusterCollectionVec->push_back( zsCluster.get() );

				// prepare a pulse for this cluster
                auto zsPulse = std::make_unique<TrackerPulseImpl>();
                idZSPulseEncoder["sensorID"] = sensorID;
                idZSPulseEncoder["type"] = static_cast<int>(kEUTelSparseClusterImpl);
                idZSPulseEncoder.setCellID( zsPulse.get() );

                zsPulse->setTime(ID);
                ID++;
                //zsPulse->setCharge( sparseCluster->getTotalCharge() );
                zsPulse->setTrackerData( zsCluster.release() );
                pulseCollection->push_back( zsPulse.release() );
			}

			
		}
		PixelsOfEvents.erase(PixelsOfEvents.begin());
	}
	cerr<<"CHECK POINT 5"<<endl;
  }
  cerr<<"CHECK POINT 6"<<endl;
  
  if ( ! isDummyAlreadyExisting )
    {
		cerr<<"CHECK POINT ! isDummyAlreadyExisting"<<endl;
        if ( sparseClusterCollectionVec->size() != 0 )
        {
			cerr<<"CHECK POINT sparseClusterCollectionVec"<<endl;
            evt->addCollection( sparseClusterCollectionVec, "oOriginal_zsdata" );
        }
        else
        {
            delete sparseClusterCollectionVec;
        }
    }
	cerr<<"CHECK POINT 7"<<endl;
	// if the pulseCollection is not empty add it to the event
    if ( ! pulseCollectionExists && ( pulseCollection->size() != _initialPulseCollectionSize ))
    {
        evt->addCollection( pulseCollection, _pulseCollectionName );
    }

    /*if ( pulseCollection->size() != _initialPulseCollectionSize )
    {
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
        if ( _fillHistos ) fillHistos(event);
#endif
    }*/
    if ( ! pulseCollectionExists && ( pulseCollection->size() == _initialPulseCollectionSize ) )
    {
        delete pulseCollection;
    }

    _isFirstEvent = false;
	cerr<<"CHECK POINT 8"<<endl;
}

void EUTelProcessorALPIDEClusterFilter::end() 
{
cerr<<"IN END"<<endl;
cerr<<allCluster<<"; "<<notDouvbleCluster<<endl;
}
