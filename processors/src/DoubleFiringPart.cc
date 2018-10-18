<vector<vector<vector<int>>>pixVectorAll;

pixVectorAll.push_back(pixVector);
if(pixVectorAll.size()>n)
{
    int Nsame=0;
    for(int i=1; i<n && SameCluster(pixVectorAll[0];pixVectorAll[i]); i++) Nsame=i;
    for(int i=1; i<Nsame && i<n; i++) AddCluster(&pixVectorAll[0]; pixVectorAll[i])
    pixVectorToSave=pixVectorAll[0];
    pixVectorAll.erase(pixVectorAll.begin());
}






if(_iEvt<1000) cerr<<"CLUSTER SIZE(3): "<<cluster->size()<<endl;
					vector<int>X;
					vector<int>Y;
					for(int iClusterrrr=0; iClusterrrr<cluster->size();iClusterrrr++)
					{
						auto& pixel = cluster->at( iClusterrrr );
						X.push_back(pixel.getXCoord());
						Y.push_back(pixel.getYCoord());
						//cout<<"X: "<<X[iClusterrrr]<<", Y: "<<Y[iClusterrrr]<<endl;
					}
					vector<vector<int>>pixVector;
					vector<int>pix;
					for(iPix=0; iPix<X.size();iPix++)
					{
						pix.push_back(X[iPixel]);
						pix.push_back(Y[iPixel]);
						pixVector.push_back(pix);
					}
					
					int same=0;
					for(iPix=0;iPix<pixVector.size();iPix++)
					{
						for(jPix=0;jPix<pixVectorBefore.size();jPix++)
						{
							if(pixVector[iPix]==pixVectorBefore[jPix]) { same++; break; }
						}
					}

					if(same>pixVector.size()*0.1)
					{
						cerr<<"Same"<<endl;
						for(iPix=0;iPix<pixVector.size();iPix++)
						{
							bool samePix=false;
							for(jPix=0;jPix<pixVectorBefore.size();jPix++)
							{
								if(pixVector[iPix]==pixVectorBefore[jPix]) { samePix=true; break; }
							}
							if(!samePix)
							{
								pixVectorBefore.push_back(pixVector[iPix]);
							}
						}
					}
					else
					{
						cerr<<"Not same"<<endl;
					}
