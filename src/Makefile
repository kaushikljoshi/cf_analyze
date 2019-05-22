# Project: connection_table

CPP  = g++
CC   = gcc
OBJ  = atom.o emptyBin.o frame_read_xmolout.o get_settings.o string_conversion_functions.o make_connection_table.o identify_header_lines.o connection_table.o read_table.o read_reax_table.o find_no_connt_lines.o bin_atoms.o ring_analysis.o identify_voids.o $(RES)
LINKOBJ  = atom.o emptyBin.o frame_read_xmolout.o get_settings.o string_conversion_functions.o make_connection_table.o identify_header_lines.o connection_table.o read_table.o read_reax_table.o find_no_connt_lines.o bin_atoms.o ring_analysis.o identify_voids.o $(RES)
LIBS =  -L /usr/lib  
INCS =  -I /usr/include/ 
#CXXINCS = -I /usr/include/c++/3.4.6/backward -I /usr/include/c++/3.4.2 -g 
CXXINCS = -I /usr/include/c++/3.4.6/backward -I /usr/include/c++/3.4.2 -std=c++11
BIN  = cf_analyze
CXXFLAGS = $(CXXINCS)  
CFLAGS = $(INCS)  
RM = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before cf_analyze all-after


clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o "cf_analyze" $(LIBS)

atom.o: atom.cpp
	$(CPP) -c atom.cpp -o atom.o $(CXXFLAGS)

emptyBin.o: emptyBin.cpp
	$(CPP) -c emptyBin.cpp -o emptyBin.o $(CXXFLAGS)

frame_read_xmolout.o: frame_read_xmolout.cpp
	$(CPP) -c frame_read_xmolout.cpp -o frame_read_xmolout.o $(CXXFLAGS)

get_settings.o: get_settings.cpp
	$(CPP) -c get_settings.cpp -o get_settings.o $(CXXFLAGS)

string_conversion_functions.o: string_conversion_functions.cpp
	$(CPP) -c string_conversion_functions.cpp -o string_conversion_functions.o $(CXXFLAGS)

make_connection_table.o: make_connection_table.cpp
	$(CPP) -c make_connection_table.cpp -o make_connection_table.o $(CXXFLAGS)

identify_header_lines.o: identify_header_lines.cpp
	$(CPP) -c identify_header_lines.cpp -o identify_header_lines.o $(CXXFLAGS)

connection_table.o: connection_table.cpp
	$(CPP) -c connection_table.cpp -o connection_table.o $(CXXFLAGS)

read_table.o: read_table.cpp
	$(CPP) -c read_table.cpp -o read_table.o $(CXXFLAGS)

read_reax_table.o: read_reax_table.cpp
	$(CPP) -c read_reax_table.cpp -o read_reax_table.o $(CXXFLAGS)

find_no_connt_lines.o: find_no_connt_lines.cpp
	$(CPP) -c find_no_connt_lines.cpp -o find_no_connt_lines.o $(CXXFLAGS)

bin_atoms.o: bin_atoms.cpp
	$(CPP) -c bin_atoms.cpp -o bin_atoms.o $(CXXFLAGS)

ring_analysis.o: ring_analysis.cpp
	$(CPP) -c ring_analysis.cpp -o ring_analysis.o $(CXXFLAGS)

identify_voids.o: identify_voids.cpp
	$(CPP) -c identify_voids.cpp -o identify_voids.o $(CXXFLAGS)
