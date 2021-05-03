# sequence-alignment-gpu


### Run on Barkla (no prerequisites required)
       1. Get code base:
              cd
              tar -xf ~sgrszafa/sequence-alignment-gpu.tgz
              cd sequence-alignment-gpu
       2. Load Cuda and a C compiler:
              module load libs/nvidia-cuda/10.1.168/bin && module load compilers/gcc/8.3.0
       3. Run:
              // Will run on Nvidia Quadro P4000 (intended for lighter use).
              ./test
              ./benchmark
              ./alignSequence

              // Will run on Nvidia V100 (or P100) using sbatch.
              // Fill out barkla_alignSequence.sh with your desired input sequences.
              sbatch barkla_runBenchmark.sh
              sbatch barkla_runTest.sh
              sbatch barkla_alignSequence.sh


### Run locally
       1. Check prerequisites.
       2. Unzip code
       3. Compile with:
              make -j
       4. Run:
              // Main program
              Usage: ./alignSequence [-d|-p] [-c|-g] [--global|--local] [-s <file>] [--gap-penalty <int>] <file> <file>
                     -d, --dna             - align dna sequences (default)
                     -p, --protein         - align protein sequence
                     -c, --cpu             - use cpu device (default)
                     -g, --gpu             - use gpu device
                     --global              - use global alignment (default)
                     --local               - use local alignment
                     -s, --score-matrix    - next argument is a score matrix file
                     --gap-penalty         - next argument is a gap open penalty (default 5)

              // test suite
              ./test

              // performance tests
              ./benchmark


### Prerequisites
       - Nvidia GPU with compute capability >3.0 (Kepler, and later)
       - CUDA (tested with 10.1 and 11.1)
       - A C++14 compliant host compiler (tested with clang and gcc)
       - Make
