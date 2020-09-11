#include <iterator_set.h>
#include <string>

using namespace std;


int main(){
    std::vector<double> state = {0,0.1,0.2,0.3, 1,1.1,1.2,1.3, 2,2.1,2.2,2.3, 3,3.1,3.2,3.3};
    std::vector<std::string> names = {"x0", "x1", "x2", "x3"};

    IteratorSet<std::vector<double>::iterator> iset(state.begin(), names, 4);
    for (int i=0; i<4; ++i){
		auto it = iset.get(); // get invalidates after ++iset
        std::cout << *it[0] << "\t" << *it[1] << "\t" << *it[2] << "\t" << *it[3] << "\n";
        ++iset;
    }
    cout << "----------------------------------------------------\n";

	auto& it0 = iset.get("x0");
	auto& it1 = iset.get("x1");
	auto& it2 = iset.get("x2");
	auto& it3 = iset.get("x3");
    for (iset.begin(); !iset.end(); ++iset){
        std::cout << *it0 << "\t" << *it1 << "\t" << *it2 << "\t" << *it3 << "\n";
    }
    cout << "----------------------------------------------------\n";

    IteratorSet<std::vector<double>::iterator> iset1(state.begin(), names, 4);
    iset1.printHeader();
    for (int i=0; i<4; ++i){
    	iset1.print();
    	++iset1;
    }
	return 0;    
}

