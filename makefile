
OUTPUT = compchem

CC = g++
FLAGS = -Wall -Werror -Wextra -O3 -pedantic -ffast-math -std=c++20

SOURCE = TN_Atom.cpp TN_Experiment.cpp compchem.cpp 

all:

	@$(CC) $(FLAGS) $(SOURCE) -o $(OUTPUT)
	@echo "project compiled!"

clean:
	@rm -f *.o
	@echo "object files cleaned!"

fclean: clean
	@rm -f $(OUTPUT)
	@echo "project cleaned!"

re: fclean clean all