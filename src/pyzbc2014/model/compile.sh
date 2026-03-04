# This is the compilation script that I use when developing the wrapper Each .c file is
# compiled separately and then they are linked into a shared library file called
# `libzbc2014.so`. Flags are used as follows:
#
# `-c`:      Compile into object files without linking into executable
# `-Wall`:   Enable all warnings
# `-fPIC`:   Generate position-independent code (suitable for shared library)
# `-03`:     Aggressively optimize
# `-shared`: Generate a shared library instead of an executable

gcc -c -Wall -fPIC -O3 complex.c 
gcc -c -Wall -fPIC -O3 model_IHC.c
gcc -c -Wall -fPIC -O3 model_Synapse.c
gcc -shared -o libzbc2014.so complex.o model_IHC.o model_Synapse.o