#include<iostream>
#include<vector>

using namespace std;

vector<int> MyVectorFunction()
{
    vector<int> vvv;
    vvv.push_back(1);
    vvv.push_back(2);
    return vvv;
}

int main()
{
    vector<int>v =MyVectorFunction();
    cout<<v[0]<<", "<<v[1]<<endl;
    return 0;
}