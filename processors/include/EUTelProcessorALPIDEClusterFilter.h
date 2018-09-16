#ifndef EUTelProcessorALPIDEClusterFilter_h
#define EUTelProcessorALPIDEClusterFilter_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <algorithm>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerDataImpl.h>
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include <UTIL/LCTime.h>
#include <UTIL/CellIDEncoder.h>

#include "TH1.h"
#include "TH2.h"

class EUTelProcessorALPIDEClusterFilter : public marlin::Processor {
public:
  virtual Processor* newProcessor() {return new EUTelProcessorALPIDEClusterFilter;}
  EUTelProcessorALPIDEClusterFilter();
  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ;
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
  void bookHistos();
#endif

  virtual void end();

protected:

private:
bool isDead;

};
#endif