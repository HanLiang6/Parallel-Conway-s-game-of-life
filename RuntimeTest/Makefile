CXX=mpic++
# CCFLAGS=-Wall -g
# activate for compiler optimizations:
CCFLAGS=-Wall --std=c++17 -O3 -g -Wall
LDFLAGS=

all: life

.PHONY : tiny_test

life: life.cpp
	$(CXX) $(CCFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o life

tiny_test: life 5_example.txt 3_example.txt
	mpirun -np 1 ./life 5 5 1 5_example.txt output.txt
	@cmp --silent output.txt 5_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS TINY TEST"; else echo "Fail tiny test"; fi

	mpirun -np 1 ./life 3 3 1 3_example.txt output.txt
	@cmp --silent output.txt 3_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS TINY TEST"; else echo "Fail tiny test"; fi


non_square_test: life 3x8_example.txt 8x3_example.txt
	mpirun -np 1 ./life 3 8 1 3x8_example.txt output.txt
	@cmp --silent output.txt 3x8_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS NON-SQUARE TEST"; else echo "Fail non-square test"; fi

	mpirun -np 1 ./life 8 3 1 8x3_example.txt output.txt
	@cmp --silent output.txt 8x3_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS NON-SQUARE TEST"; else echo "Fail non-square test"; fi

test:
	mpirun -np 1 ./life 5000 5000 1 new_grid_5K.txt new_output.text
	mpirun -np 4 ./life 5000 5000 1 new_grid_5K.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 5k TEST"; else echo "Fail 5k test"; fi
	
	mpirun -np 1 ./life 2000 2000 1 new_grid_2K.txt new_output.text
	mpirun -np 4 ./life 2000 2000 1 new_grid_2K.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 2k TEST"; else echo "Fail 2k test"; fi
	
	mpirun -np 1 ./life 1000 1000 1 new_grid_1K.txt new_output.text
	mpirun -np 4 ./life 1000 1000 1 new_grid_1K.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 1k TEST"; else echo "Fail 1k test"; fi
	
	mpirun -np 1 ./life 500 500 1 new_grid_500.txt new_output.text
	mpirun -np 4 ./life 500 500 1 new_grid_500.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 500 TEST"; else echo "Fail 500 test"; fi
	
	mpirun -np 1 ./life 250 250 1 new_grid_250.txt new_output.text
	mpirun -np 4 ./life 250 250 1 new_grid_250.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 250 TEST"; else echo "Fail 250 test"; fi
	
	mpirun -np 1 ./life 100 100 1 new_grid_100.txt new_output.text
	mpirun -np 4 ./life 100 100 1 new_grid_100.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 100 TEST"; else echo "Fail 100 test"; fi

non_square_test_np: life 3x8_example.txt 8x3_example.txt
	mpirun --use-hwthread-cpus -np 4 ./life 3 8 1 3x8_example.txt output.txt
	@cmp --silent output.txt 3x8_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS NON-SQUARE TEST"; else echo "Fail non-square test"; fi
	
	mpirun --use-hwthread-cpus -np 4 ./life 8 3 1 8x3_example.txt output.txt
	@cmp --silent output.txt 8x3_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS NON-SQUARE TEST"; else echo "Fail non-square test"; fi

tiny_test_np: life 5_example.txt 3_example.txt
	mpirun --use-hwthread-cpus -np 4 ./life 5 5 1 5_example.txt output.txt
	@cmp --silent output.txt 5_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS TINY TEST"; else echo "Fail tiny test"; fi
	
	mpirun --use-hwthread-cpus -np 4 ./life 3 3 1 3_example.txt output.txt
	@cmp --silent output.txt 3_ref.txt; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS TINY TEST"; else echo "Fail tiny test"; fi

small_test:
	mpirun --use-hwthread-cpus -np 1 ./life 5000 5000 1 new_grid_5K.txt new_output.text
	mpirun --use-hwthread-cpus -np 4 ./life 5000 5000 1 new_grid_5K.txt new_output_np.text
	@cmp --silent new_output.text new_output_np.text; RETVAL=$$?; if [ $$RETVAL -eq 0 ]; then echo "PASS 5k TEST"; else echo "Fail 5k test"; fi
