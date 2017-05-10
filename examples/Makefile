-include /home1/03838/rsatija/pvfmm/MakeVariables

ifndef CXX_PVFMM
-include MakeVariables
endif

ifndef CXXFLAGS_PVFMM
$(error Cannot find file: MakeVariables)
endif

CXX=$(CXX_PVFMM)
CXXFLAGS=$(CXXFLAGS_PVFMM) -std=c++11
LDLIBS=$(LDLIBS_PVFMM)

RM = rm -f
MKDIRS = mkdir -p

BINDIR = ./bin
SRCDIR = ./src
OBJDIR = ./obj
INCDIR = ./include

TARGET_BIN = \
       $(BINDIR)/runMD
#       $(BINDIR)/example1 \
       
SOURCE_FILES =  energy.cpp forces.cpp killer.cpp parser.cpp forces_pvfmm.cpp init.cpp old_verlet.cpp verlet.cpp init_par.cpp potentials.cpp #verlet_par.cpp

OBJ_FILES =  obj/energy.o obj/forces.o obj/killer.o obj/parser.o obj/forces_pvfmm.o obj/init.o obj/verlet.o obj/init_par.o obj/potentials.o #obj/verlet_par.o

# $(addprefix obj/, $($(SOURCE_FILES):.cpp=.o))


all : $(TARGET_BIN)

$(BINDIR)/%: $(OBJ_FILES) $(OBJDIR)/%.o
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) $^ $(LDLIBS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-@$(MKDIRS) $(dir $@)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $^ -o $@

clean:
	$(RM) -r $(BINDIR)/* $(OBJDIR)/*
	$(RM) *~ */*~

