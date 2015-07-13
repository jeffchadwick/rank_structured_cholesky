
PYTHON=python
UFBASE=http://www.cise.ufl.edu/research/sparse/mat
RUNNER=$(PYTHON) runner.py
CSUB=csub -f pac.txt $(RUNNER) run

.PHONY: download setup

setup:
	$(RUNNER) setup HB/1138_bus
	$(RUNNER) setup Janna/Flan_1565

run:
	$(RUNNER) run HB/1138_bus
	$(RUNNER) run Janna/Flan_1565

sub_1138_bus:
	$(CSUB) HB/1138_bus

sub_Flan_1565:
	$(CSUB) Janna/Flan_1565

clean:
	rm csub*

distclean: clean
	rm -rf Flan_1565*
	rm -rf 1138_bus*
