# MAKEFILE FOR CUDA

#-------------------------------------------------------------------------------
# executable name
TARGET := 

# files
CCFILES  :=  $(filter-out ../src/main.cpp, $(wildcard ../src/*.cpp))

# ------------------------------------------------------------------------------

# paths
#CUDA_INSTALL_PATH ?= /usr/local/cuda#-5.0

# compilers

# include and lib dirs (esp for cuda)
INC_PATH := -I../include #-I/home/jaideep/codes/quickGL/include
LIB_PATH := #-L/home/jaideep/codes/quickGL/lib
GLLIB_PATH := 

# flags
COMMONFLAGS = 
CPPFLAGS = -O3 -std=c++11 -g -pg -Wno-c++11narrowing # -Wall -Wno-unused-variable  #-fno-inline-functions -fno-optimize-sibling-calls  
LINKFLAGS += -pg $(COMMONFLAGS) 

# libs
#LIBS = -lcudart 					# cuda libs 		-lcutil_x86_64 -lshrutil_x86_64
GLLIBS = #-lquickgl -lglut 				# openGL libs       -lGL -lGLEW  #-lX11 -lXi -lXmu 		
LIBS = 	 -lgsl -lgslcblas 	# additional libs

# files
OBJECTS = $(patsubst ../src/%.cpp, ../build/%.o, $(CCFILES))
TARGET_OBJ = build/$(TARGET).o


all: dir $(TARGET)

dir:	
	mkdir -p build
	rm -f $(TARGET_OBJ)

$(TARGET): $(TARGET_OBJ) $(OBJECTS) 
	g++ $(LINKFLAGS) -o exec_$(TARGET) $(LIB_PATH) $(GLLIB_PATH) $(OBJECTS) $(TARGET_OBJ) $(LIBS) $(GLLIBS)

$(OBJECTS): ../build/%.o : ../src/%.cpp 
	g++ -c $(CPPFLAGS) $(INC_PATH) $< -o $@ 

$(TARGET_OBJ):
	g++ -c $(CPPFLAGS) $(INC_PATH) $(TARGET).cpp -o $(TARGET_OBJ)

superclean:
	rm -f exec_$(TARGET) ../build/*.o build/*.o 

clean: 
	rm -f exec_$(TARGET) build/*.o

re: clean all

# ------------------------------------------------------------------------------





#-gencode=arch=compute_10,code=\"sm_10,compute_10\"  -gencode=arch=compute_20,code=\"sm_20,compute_20\"  -gencode=arch=compute_30,code=\"sm_30,compute_30\" 

#-W -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function 
#-m64 -fno-strict-aliasing 
#-I. -I/usr/local/cuda/include -I../../common/inc -I../../../shared//inc 
#-DUNIX -O2


#g++ -fPIC   -m64 -o ../../bin/linux/release/swarming_chasing_predator obj/x86_64/release/genmtrand.cpp.o  obj/x86_64/release/simpleGL.cu.o  -L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart   
#-lGL -lGLU -lX11 -lXi -lXmu -lGLEW_x86_64 -L/usr/X11R6/lib64 -lGLEW_x86_64 -L/usr/X11R6/lib64 -lglut 
#-L/usr/local/cuda/lib64 -L../../lib -L../../common/lib/linux -L../../../shared//lib -lcudart 
#-L/usr/lib -lgsl -lgslcblas 
#-lcutil_x86_64  -lshrutil_x86_64 




#CXXWARN_FLAGS := \
#	-W -Wall \
#	-Wimplicit \
#	-Wswitch \
#	-Wformat \
#	-Wchar-subscripts \
#	-Wparentheses \
#	-Wmultichar \
#	-Wtrigraphs \
#	-Wpointer-arith \
#	-Wcast-align \
#	-Wreturn-type \
#	-Wno-unused-function \
#	$(SPACE)

#CWARN_FLAGS := $(CXXWARN_FLAGS) \
#	-Wstrict-prototypes \
#	-Wmissing-prototypes \
#	-Wmissing-declarations \
#	-Wnested-externs \
#	-Wmain \
#	
	
#HEADERS  := $(wildcard *.h)
	
