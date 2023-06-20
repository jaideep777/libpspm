// Extend tensor to add boundary cohort

#include "tensor.h"

template <class T>
class TensorCohort : public Tensor{

    public: 
    
    TensorCohort(std::vector<int> _dim) : Cohort (std::vector<int> _dim){
        nelem = nelem + 1;
		vec.resize(nelem);
    }

};