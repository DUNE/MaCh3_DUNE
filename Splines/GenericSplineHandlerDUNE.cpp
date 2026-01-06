#include "GenericSplineHandlerDUNE.h"

GenericSplineHandlerDUNE::GenericSplineHandlerDUNE(ParameterHandlerGeneric* xsec, MaCh3Modes* Modes_) 
{
    YAML::Node ConfigCurrent = xsec->GetConfig();
    //Dump the whole config
    std::cout << ConfigCurrent << std::endl;
}

GenericSplineHandlerDUNE::~GenericSplineHandlerDUNE() {
}



