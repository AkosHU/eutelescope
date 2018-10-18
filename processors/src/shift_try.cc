void EUTelProcessorALPIDEClusterFilter::shift() {
	if(_shiftedRun) {
		if(PixelsOfEvents.size()>_nLayers*_nShift) {
			vector<vector<vector<int>>>ActualEvent;
			for(int iLayer=0; iLayer<_nLayers; iLayer++) {
				for(int iCluster=0; iCluster<PixelsOfEvents[_nShift*iLayer].size(); iCluster++) {
					if(PixelsOfEvents[_nShift*iLayer][iCluster][0][2]==iLayer) {
						ActualEvent.push_back(PixelsOfEvents[_nShift*iLayer][iCluster]);
					}
				}
			}
			PixelsOfEventsAfterShift.push_back(ActualEvent);
		}
	} 
	else {
		while(PixelsOfEvents.size()>0) {
			PixelsOfEventsAfterShift.push_back(PixelsOfEvents[0]);
			PixelsOfEvents.erase(PixelsOfEvents.begin());
		}
	}
	cerr<<"SIZE!!!!! "<<PixelsOfEventsAfterShift.size()<<"!!!!!!"<<endl;
}