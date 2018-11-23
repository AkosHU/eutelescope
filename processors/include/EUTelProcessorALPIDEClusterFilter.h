/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTelProcessorALPIDEClusterFilter_H
#define EUTelProcessorALPIDEClusterFilter_H
#include "cluster.h"

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"
#include "EUTelUtility.h"
#include "EUTelEventImpl.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

#include <IMPL/LCCollectionVec.h>


// system includes <>
#include <string>
#include <vector>
#include <map>
#include <set>

namespace eutelescope {

  

  class EUTelProcessorALPIDEClusterFilter : public marlin::Processor {

  private:
      DISALLOW_COPY_AND_ASSIGN(EUTelProcessorALPIDEClusterFilter)
    bool _clusterAvailable;
    LCCollectionVec *zsInputDataCollectionVec;
    std::string _zsDataCollectionName;
    std::vector<std::vector<std::vector<std::vector<int>>>>PixelsOfEvents;
    std::vector<std::vector<std::vector<std::vector<int>>>>PixelsOfEventsAfterShift;
    int _nDeep;
    float _Range;


  public:


   
    virtual Processor * newProcessor() {
      return new EUTelProcessorALPIDEClusterFilter;
    }

    EUTelProcessorALPIDEClusterFilter ();

  
    virtual void init ();

    virtual bool SameCluster(int iEvent, int iCluster, int jEvent,int jCluster);

    virtual void AddCluster(int iEvent, int iCluster, int jEvent,int jCluster);

    virtual void DeletCluster(int jEvents, int jCluster);

    
    virtual void processEvent (LCEvent * evt);


    virtual void end();

    virtual void readCollections(LCCollectionVec * zsInputDataCollectionVec);

    virtual void writeCollection(LCCollectionVec * sparseClusterCollectionVec, LCCollectionVec * pulseCollection);

    virtual void filter();

    virtual void effOfALPIDE();


    


   
    
  protected:


    
    //! Pulse collection size
    size_t _initialPulseCollectionSize;

    //! Pulse collection name.
    /*! This is the name used to store the output cluster
     *  collection.
     */
    std::string _pulseCollectionName;

        //! Input collection name for NZS data
    /*! The input collection is the calibrated data one coming from
     *  the EUTelCalibrateEventProcessor. It is, usually, called
     *  "nzsdata" and it is a collection of TrackerData
     */
    std::string _nzsDataCollectionName;

    //! Noise collection name.
    /*! See _pedestalCollectionName for the detailed description
     */
    std::string _noiseCollectionName;



    //
    //! noise Collection
    LCCollectionVec *noiseCollectionVec;

    std::string _sparseClusterCollectionName;

    int allCluster;
    int notDouvbleCluster;
    int ID;
    int _nShift;
    bool _shiftedRun;
    int _nLayers;
    int _noOfDetector;
    bool _isGeometryReady;
    std::map< int, int > _totClusterMap;
    int _NoEvent;
    int _nOfAll;
    int _nOfNoise;
    int _nOfGood;
    int _nOfFalse;

  private:

    
 
  };


  EUTelProcessorALPIDEClusterFilter gEUTelProcessorALPIDEClusterFilter;

}
#endif
#endif