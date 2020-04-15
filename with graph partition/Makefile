exe=trianglecounting.bin
N=1
cucc= "$(shell which nvcc)"
cc= "$(shell which mpicxx)"
commflags=-lcudart -L"$(shell dirname $(cucc))"/../lib64  -O3 -W -Wall -Wno-unused-function -Wno-unused-parameter
cuflags= --compiler-options -Wall --gpu-architecture=compute_70 --gpu-code=sm_70  -O3   # --resource-usage 
# cuflags= --compiler-options -v -Xcudafe  -\# --resource-usage 

.SILENT: cucc
.SILENT: cc
.SILENT: cuflags
.SILENT: %.o


objs	= $(patsubst %.cu,%.o,$(wildcard *.cu)) \
	$(patsubst %.cpp,%.o,$(wildcard *.cpp))

deps	= $(wildcard ./*.cuh) \
	$(wildcard ./*.hpp) \
	$(wildcard ./*.h) \


# %.o:%.cu $(deps)
# 	$(cucc) -c  $(cuflags)  $< -o $@ -g -G 

# %.o:%.cpp $(deps)
# 	$(cc) -c -g  $< -o $@ 

# $(exe):$(objs)
# 	$(cc) $(objs) $(commflags) -O3 -o $(exe)

%.o:%.cu $(deps)
	$(cucc) -c $(cuflags) $<  -o $@ 

%.o:%.cpp $(deps)
	$(cc) -c  $(commflags) $< -o $@ 

$(exe):$(objs)
	$(cc) $(objs) $(commflags) -o $(exe)

test:$(exe)
	@#mpirun -n $(N) $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/roadNet-CA_adj/ $(N) 128 1280 32 0 0 60

	@#mpirun -n $(N) $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/p2p-Gnutella08_adj/ $(N) 32 1 32 0 0 60
	mpirun -n $(N) $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/201512012345.v18571154_e38040320/ $(N) 128 1280 32 0 0 60
	@#mpirun -n $(N) $(exe) cuda-memcheck /data/p2p/input/ $(N) 128 5 32 0 0 
	@#mpirun -n $(N) $(exe) /gpfs/alpine/proj-shared/csc289/trianglecounting/snap/email-EuAll_adj/ $(N) 32 10 32 0 0 


clean:
	rm -rf *.o ${exe}
