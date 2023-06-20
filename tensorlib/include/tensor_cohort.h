// Extend tensor to add boundary cohort

#include "tensor.h"

template <class T>
class TensorCohort : public Tensor<T> {

    public: 
    
    TensorCohort(std::vector<int> _dim) : Tensor<T> (_dim){
        this->nelem = this->nelem + 1;
		this->vec.resize(this->nelem);
    }

};