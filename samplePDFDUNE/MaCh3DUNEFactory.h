#pragma once

// Include the input manager
#include "manager/manager.h"

// Include the samplePDFs
#include "samplePDF/samplePDFFDBase.h"

void MakeMaCh3DuneInstance(manager *fitMan, std::vector<samplePDFFDBase*> &sample_vec,  covarianceXsec *&xsec, covarianceOsc *&osc);
samplePDFFDBase* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, covarianceXsec* &xsec);
