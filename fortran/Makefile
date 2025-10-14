FC       = gfortran
FFLAGS   = -O3 -march=native -mtune=native -fopenmp -flto -std=f2018

BAOAB_SRC    = NVT_Langevin_BAOAB/md_baoab.f90
BAOAB_BIN    = NVT_Langevin_BAOAB/bin/md_baoab

MC_BUILD_DIR = MuVT_MonteCarlo/build
MC_BIN       = MuVT_MonteCarlo/bin/mc_muvt

MC_MODULES = rng_mod system_mod potential_mod mc_mod analysis_mod io_mod
MC_MODULE_SRC = $(addprefix MuVT_MonteCarlo/modules/,$(addsuffix .f90,$(MC_MODULES)))
MC_MAIN = MuVT_MonteCarlo/mc_muvt.f90

.PHONY: all clean help

all: $(BAOAB_BIN) $(MC_BIN)

# ----------------------
# NVT_Langevin_BAOAB
# ----------------------
$(BAOAB_BIN): $(BAOAB_SRC)
	@mkdir -p $(dir $@)
	$(FC) $(FFLAGS) -o $@ $<

# ----------------------
# MuVT_MonteCarlo
# ----------------------
$(MC_BUILD_DIR):
	mkdir -p $(MC_BUILD_DIR)

$(MC_BUILD_DIR)/%.o: MuVT_MonteCarlo/modules/%.f90 | $(MC_BUILD_DIR)
	$(FC) $(FFLAGS) -c -J$(MC_BUILD_DIR) $< -o $@

$(MC_BIN): $(MC_MODULE_SRC:MuVT_MonteCarlo/modules/%.f90=$(MC_BUILD_DIR)/%.o) $(MC_MAIN) | $(MC_BUILD_DIR)
	$(FC) $(FFLAGS) -c -J$(MC_BUILD_DIR) $(MC_MAIN) -o $(MC_BUILD_DIR)/mc_muvt.o
	$(FC) $(FFLAGS) -o $(MC_BIN) $(MC_BUILD_DIR)/*.o

clean:
	rm -rf $(MC_BIN) $(MC_BUILD_DIR)/*.o $(MC_BUILD_DIR)/*.mod $(BAOAB_BIN)

help:
	@echo "Usage:"
	@echo "  make        -> compila ambos proyectos (NVT_Langevin_BAOAB y MuVT_MonteCarlo)"
	@echo "  make clean  -> elimina binarios y objetos"


# ../fortran/simulate 10000 2000 0.005 50.0 2.5 0.3 12345
# gfortran -O3 -march=native -mtune=native -fopenmp -flto -std=f2018 -o simulate simulate.f90
# docker compose up --build
# ../../fortran/MuVT_MonteCarlo/bin/mc_muvt 20 300 5.0 0.1 5000 test.bin 0.5 5.0 0.1 1.0
