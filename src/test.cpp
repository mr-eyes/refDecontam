#include "kDataFrame.hpp"

int main(){
    auto * kf = new kDataFramePHMAP(9);
    kf->insert("ACTGACAAG");
    cout << kf->size() << endl;
    
}