#include "EUTelProcessorALPIDEClusterFilter.h"
#include "EUTELESCOPE.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTelGenericSparsePixel.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "marlin/Global.h"

using namespace lcio;
using namespace marlin;
using namespace std;
using namespace eutelescope;
using namespace gear;

EUTelProcessorALPIDEClusterFilter aEUTelProcessorALPIDEClusterFilter;

EUTelProcessorALPIDEClusterFilter::EUTelProcessorALPIDEClusterFilter()
: Processor("EUTelProcessorALPIDEClusterFilter"),
isDead(0)
//ismeretlenek
  {
//ismeretlenek beolvasása
  }

void EUTelProcessorALPIDEClusterFilter::init() {
//első esemény előtt
}

void EUTelProcessorALPIDEClusterFilter::processEvent(LCEvent *evt)
{
//minden eseménynél
}

void EUTelProcessorALPIDEClusterFilter::bookHistos()
{
//histogrammok foglalása
}

void EUTelProcessorALPIDEClusterFilter::end()
{
//utolsó esemény után
}
